"""Offline tests: no ECDS API, no Cloudflare network, no Cartopy downloads."""
from datetime import datetime, timezone
import tempfile
from pathlib import Path
import os
import sys
from unittest.mock import patch
import types
import numpy as np
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import pipeline


def fake_arrays():
    rng = np.random.default_rng(42)
    n_members = 12
    hours = np.arange(24, 1105, 24)
    lat = np.arange(55., 27., -1.5)
    lon = np.arange(-15., 45., 1.5)
    temp = np.empty((n_members, len(hours), len(lat), len(lon)), dtype=np.float32)
    prec = np.empty_like(temp)
    for i, hour in enumerate(hours):
        temp[:, i] = 273.15 + 10 + .25 * (lat[None, :, None] - 40) + \
                     rng.normal(0, 2, (n_members, len(lat), len(lon)))
        prec[:, i] = (hour / 24) * (2 + .1 * np.arange(n_members)[:, None, None]) \
                     + .01 * np.arange(len(lon))[None, None, :]
    coords = dict(member=np.arange(1, n_members + 1), step=hours.astype('timedelta64[h]'),
                  latitude=lat, longitude=lon)
    return (xr.DataArray(temp, dims=tuple(coords), coords=coords, attrs={'units': 'K'}),
            xr.DataArray(prec, dims=tuple(coords), coords=coords,
                         attrs={'units': 'kg m**-2'}))


def test_dates_and_ecd_api():
    for value, expected in [
        (datetime(2026, 10, 24, 12, 15, tzinfo=timezone.utc), '2026-10-22'),
        (datetime(2026, 10, 25, 12, 15, tzinfo=timezone.utc), '2026-10-23'),
        (datetime(2026, 11, 1, 12, 15, tzinfo=timezone.utc), '2026-10-30'),
    ]:
        assert pipeline.date_from_args(now=value).strftime('%Y-%m-%d') == expected
    req = pipeline.request_for(datetime(2026, 10, 22, tzinfo=timezone.utc))
    assert req['variable'] == ['2_m_temperature', 'total_precipitation']
    assert req['origin'] == 'ecmwf'
    assert req['forecast_type'] == 'perturbed_forecast'
    assert req['data_format'] == 'grib'
    assert req['area'] == [54, -14, 28, 44]
    assert req['year'] == '2026' and req['month'] == '10' and req['day'] == '22'
    assert req['time'] == '00:00' and req['level_type'] == 'single_level'
    assert len(req['leadtime_hour']) == 92
    assert req['leadtime_hour'][:2] == ['0_24', '24_48']
    assert req['leadtime_hour'][46:48] == ['24', '48']
    assert req['leadtime_hour'][-1] == '1104'
    assert 'class' not in req and 'param' not in req and 'step' not in req


def test_cds_api_download_mock():
    """Download MUST make one ECDS request containing both parameters."""
    requests = []
    class FakeResult:
        def download(self, filename):
            Path(filename).write_bytes(b'GRIB' + b'\0' * 128)
    class FakeClient:
        def __init__(self, url, key):
            assert url == 'https://ecds.ecmwf.int/api' and key == 'dummy-token'
        def retrieve(self, dataset, request):
            assert dataset == 's2s-forecasts'
            requests.append(request)
            return FakeResult()
    with tempfile.TemporaryDirectory() as td, \
         patch.dict(os.environ, {'ECDS_API_KEY': 'dummy-token'}), \
         patch.dict(sys.modules, {'cdsapi': types.SimpleNamespace(Client=FakeClient)}):
        file = pipeline.download_grib(datetime(2026,10,22,tzinfo=timezone.utc), Path(td))
        assert file.is_file() and file.read_bytes().startswith(b'GRIB')
        assert file.name == 'ecmwf_s2s_t2m_tp_20261022.grib'
        assert len(requests) == 1
        assert requests[0]['variable'] == ['2_m_temperature', 'total_precipitation']
        # Existing cached file -> no second API request.
        pipeline.download_grib(datetime(2026,10,22,tzinfo=timezone.utc), Path(td))
        assert len(requests) == 1


def test_split_combined_grib_mock():
    """Check GRIB metadata routing, order independence, and ZIP input.

    The ecCodes API is mocked, so no remote service or system GRIB libraries
    are needed. Each output stream must contain ONLY its parameter.
    """
    from zipfile import ZipFile
    init = datetime(2026,10,22,tzinfo=timezone.utc)
    messages = [
        {'shortName': 'tp', 'stepType': 'accum', 'blob': b'GRIB-P1'},
        {'shortName': 'mean2t24', 'stepType': 'avg', 'blob': b'GRIB-T1'},
        {'shortName': 'tp', 'stepType': 'accum', 'blob': b'GRIB-P2'},
        {'shortName': 'mean2t24', 'stepType': 'avg', 'blob': b'GRIB-T2'},
    ]
    def fake_new(f):
        # Simulate a low-memory GRIB stream parser; a real ecCodes parser
        # extracts messages from the file rather than using this list.
        nextmsg = getattr(f, '_mock_idx', 0)
        if nextmsg >= len(messages):
            return None
        f._mock_idx = nextmsg + 1
        return messages[nextmsg]
    def fake_get(gid, key):
        return gid[key]
    fake_codes = types.SimpleNamespace(
        codes_grib_new_from_file=fake_new,
        codes_get=fake_get,
        codes_get_message=lambda gid: gid['blob'],
        codes_release=lambda gid: None,
    )
    with tempfile.TemporaryDirectory() as td, \
         patch.dict(sys.modules, {'eccodes': fake_codes}):
        directory = Path(td)
        source = directory / 'joint.grib'
        source.write_bytes(b'GRIB' + b'\x00'*128)
        paths = pipeline.split_combined_grib(source, directory, init)
        assert paths['t2m'].read_bytes() == b'GRIB-T1GRIB-T2'
        assert paths['tp'].read_bytes() == b'GRIB-P1GRIB-P2'
        assert source.is_file()  # combined file preserved
        # ZIP response: GRIB extraction does not modify the source archive.
        archive = directory / 'joint.zip'
        with ZipFile(archive, 'w') as z:
            z.writestr('part.grib', b'GRIB' + b'\x00'*128)
        paths_zip = pipeline.split_combined_grib(archive, directory, init)
        assert paths_zip['t2m'].read_bytes() == b'GRIB-T1GRIB-T2'
        assert paths_zip['tp'].read_bytes() == b'GRIB-P1GRIB-P2'
        assert not list(directory.glob('.ecds_archive_part_*'))

    # Strictly forbid maximum temperatures from being misused with mean bias.
    bad_temperature = {'shortName': 'mx2t24', 'stepType': 'max'}
    try:
        pipeline._grib_product(bad_temperature, lambda gid, key: gid[key])
    except ValueError as ex:
        assert 'non è media' in str(ex)
    else:
        raise AssertionError('La massima giornaliera non è una media')


class FakeR2:
    def __init__(self):
        self.calls = []
    def upload_file(self, path, bucket, key, ExtraArgs):
        self.calls.append((path, bucket, key, ExtraArgs))
    def head_object(self, Bucket, Key):
        matches = [c for c in self.calls if c[1] == Bucket and c[2] == Key]
        if not matches:
            raise AssertionError(f'HEAD: missing key {Key}')
        return {'ContentLength': Path(matches[-1][0]).stat().st_size}


def test_plot_webp_and_upload(city_data, thresholds):
    import plotting
    with tempfile.TemporaryDirectory() as td:
        out = Path(td)
        city = 'Bologna'
        weeks = sorted(city_data[city])
        first = out / f'{city}_{weeks[0]}.webp'
        plotting.plot_city_combined(city, weeks, city_data[city],
                                    thresholds[city], str(first), '22/10/2026')
        assert first.read_bytes()[:4] == b'RIFF'
        assert first.read_bytes()[8:12] == b'WEBP'
        from PIL import Image
        with Image.open(first) as img:
            assert img.format == 'WEBP'
            assert img.width > 0 and img.height > 0
        second = out / 'subseasonal_20261026_2tp.webp'
        second.write_bytes(first.read_bytes())
        client = FakeR2()
        with patch.dict(os.environ, {'R2_ACCOUNT_ID': 'example',
                                     'R2_ACCESS_KEY_ID': 'test',
                                     'R2_SECRET_ACCESS_KEY': 'test',
                                     'R2_BUCKET': 'settimanali', 'R2_PREFIX': '',
                                     'R2_VERIFY_OBJECTS': '1'}):
            assert pipeline.upload_webps_to_r2(out, 2, client=client) == 2
            assert {c[2] for c in client.calls} == {first.name, second.name}
            assert all(c[1] == 'settimanali' for c in client.calls)
            assert all(c[3]['ContentType'] == 'image/webp' for c in client.calls)
        assert not list(out.glob('*.nc')) and not list(out.glob('*.json'))
        try:
            with patch.dict(os.environ, {'R2_ACCOUNT_ID': 'x',
                                         'R2_ACCESS_KEY_ID': 'x',
                                         'R2_SECRET_ACCESS_KEY': 'x'}):
                pipeline.upload_webps_to_r2(out, 3, client=client)
        except RuntimeError as ex:
            assert 'Upload annullato' in str(ex)
        else:
            raise AssertionError('Missing WebP count was not rejected')



def test_temperature_bias_correction():
    """True NetCDF asset, leap calendar, forecast-week weighting and members."""
    from datetime import timedelta
    with __import__('h5py').File(pipeline.TEMPERATURE_BIAS_FILE, 'r') as f:
        assert list(f.keys()) == ['bias', 'latitude', 'longitude', 'week'] or set(f.keys()) == {'bias', 'latitude', 'longitude', 'week'}
        assert f['bias'].shape == (52, 58, 135)
        assert f['week'][:].tolist() == list(range(1, 53))
    bias = pipeline.load_temperature_bias()
    assert bias.shape == (52, 58, 136)
    assert np.all(np.isfinite(bias))

    # Jan 1 fixed weekly blocks, with leap day removed from the calendar.
    assert pipeline.climatological_week_for_day(datetime(2026,1,1)) == 0
    assert pipeline.climatological_week_for_day(datetime(2026,1,8)) == 1
    assert pipeline.climatological_week_for_day(datetime(2028,2,28)) == 8
    assert pipeline.climatological_week_for_day(datetime(2028,2,29)) == 8
    assert pipeline.climatological_week_for_day(datetime(2028,3,1)) == 8
    assert pipeline.climatological_week_for_day(datetime(2026,12,31)) == 51
    assert pipeline.climatological_week_for_day(datetime(2027,1,1)) == 0

    # Several forecast weeks, including crossings of two climatological bins
    # and calendar years. Unique member offsets must remain intact.
    starts = [datetime(2026,10,26,tzinfo=timezone.utc),
              datetime(2026,12,28,tzinfo=timezone.utc),
              datetime(2028,2,28,tzinfo=timezone.utc)]
    offset = np.array([[-2.,0.,3.],[-4.,2.,6.],[-1.,3.,9.]], dtype=np.float32)
    dummy = np.broadcast_to(offset[:, :, None, None],
                            (len(starts), 3, 58, 136)).copy() + 15
    coords = {'week': np.arange(1, len(starts)+1),
              'member': np.array([1,2,3]),
              'latitude': pipeline.LATS, 'longitude': pipeline.LONS}
    da = xr.DataArray(dummy, dims=('week','member','latitude','longitude'), coords=coords)
    da_original = da.copy(deep=True)
    corrected = pipeline.correct_ensemble_temperature(da, starts, bias)

    for n, start in enumerate(starts):
        week_indices = [pipeline.climatological_week_for_day(start + timedelta(days=j)) for j in range(7)]
        expected = np.mean(bias[week_indices], axis=0)
        for member in range(3):
            assert np.allclose(corrected.isel(week=n, member=member).values,
                               dummy[n, member] - expected, atol=3e-6)
        # The additive correction MUST NOT alter the ensemble spread.
        assert np.allclose((corrected.isel(week=n,member=2) - corrected.isel(week=n,member=0)).values,
                           dummy[n,2] - dummy[n,0], atol=3e-6)
    # The original forecast is unchanged (temperature only in new output).
    assert np.array_equal(da.values, da_original.values)
    assert np.isfinite(corrected.values).all()

    # A forecast exactly 2 deg warmer than ERA5 must become unbiased when
    # S2S-ERA5 climatological bias is +2 deg, on ALL members and pixels.
    synthetic_bias = np.full((52,58,136), 2., dtype=np.float32)
    uniform = pipeline.correct_ensemble_temperature(da, starts, synthetic_bias)
    assert np.allclose(uniform.values, da.values - 2, atol=1e-6)
    print('TEST BIAS PASS: file, calendario, remapping, tutti i membri, segno e spread')



def test_gridded_percentile_classification(tw, pw, indices, thresholds, starts,
                                           era5, era5_bounds, city_data, first_map):
    """Probabilities use P20/P40/P60/P80 at the pixel; city JSON for labels only."""
    import h5py
    from datetime import timedelta
    assert set(era5) == {'T', 'P'}
    for kind in ('T', 'P'):
        assert era5[kind].shape == (52, 4, 58, 136)
        assert era5_bounds[kind].shape == (len(starts), 4, 58, 136)
        assert np.isfinite(era5_bounds[kind]).all()
        assert np.all(np.diff(era5_bounds[kind], axis=1) >= -1e-4)

    # Real data check: exact source grid centre at 51.8N/-11.8E matches the
    # NetCDF field, even after latitude flip and interpolation.
    with h5py.File(pipeline.ERA5_PERCENTILES_FILE) as h:
        nlat, nlon = 51.8, -11.8
        # Target grid is 52, 51.6,... and -12, -11.6,... .
        # Compare interpolation at the internal midpoint to nearby source cells.
        r, c = 1, 1  # target 51.6, -11.6
        assert np.isclose(pipeline.LATS[r], 51.6)
        assert np.isclose(pipeline.LONS[c], -11.6)
        # 51.6 is midpoint of 51.8 and 51.4; -11.6 between -11.8/-11.4
        expected = np.mean(h['t2m_p20'][0, 0:2, 0:2])
        assert np.isclose(era5['T'][0,0,r,c], expected, atol=1e-4)

    # Matches the same daily-weighted climatology used by the bias correction.
    idx = [pipeline.climatological_week_for_day(starts[0] + timedelta(days=d)) for d in range(7)]
    assert np.allclose(era5_bounds['T'][0], np.mean(era5['T'][idx],axis=0), atol=1e-5)

    city = 'Bologna'
    i, j = indices[city]
    expected_t = pipeline.quintile_probs(tw.isel(week=0,latitude=i,longitude=j).values,
                                         era5_bounds['T'][0,:,i,j])
    expected_p = pipeline.quintile_probs(pw.isel(week=0,latitude=i,longitude=j).values,
                                         era5_bounds['P'][0,:,i,j])
    k = starts[0].strftime('%Y%m%d')
    actual = city_data[city][k]
    for kq in range(5):
        assert np.isclose(actual[f'T_{kq+1}'],pipeline.rounded_proba(expected_t)[kq],atol=1e-3)
        assert np.isclose(actual[f'P_{kq+1}'],pipeline.rounded_proba(expected_p)[kq],atol=1e-3)
        assert np.isclose(first_map['2tp'][kq,i,j],expected_t[kq],atol=1e-3)
        assert np.isclose(first_map['tpp'][kq,i,j],expected_p[kq],atol=1e-3)

    # When the historical city JSON is changed, classification is unchanged.
    # No old city label table is passed into the probability functions:
    # ERA5 gridded thresholds are the only classification inputs.
    again = pipeline.city_probabilities(tw, pw, indices, era5_bounds, starts)
    assert again[city][k] == actual
    print('TEST ERA5 PERCENTILES PASS: both variables, same pixel for city/map, daily mix, city labels only')


def test_all():
    test_dates_and_ecd_api()
    test_cds_api_download_mock()
    test_split_combined_grib_mock()
    test_temperature_bias_correction()
    init = datetime(2026, 10, 22, tzinfo=timezone.utc)
    starts = pipeline.week_starts(init)
    indices, thresholds = pipeline.load_assets()
    pipeline.check_threshold_weeks(indices, thresholds, starts)
    t, p = fake_arrays()
    tn = pipeline.normalize_da(t, init, 't2m')
    pn = pipeline.normalize_da(p, init, 'tp')
    assert np.isclose(float(pn.isel(member=0, day=5, latitude=0, longitude=0)), 2, atol=.01)
    tw = pipeline.weekly_interpolated(tn, starts, 't2m')
    pw = pipeline.weekly_interpolated(pn, starts, 'tp')
    # Longitude may be encoded as 0..360 by ecCodes, including western longitudes.
    tn_0360 = tn.assign_coords(longitude=np.mod(tn.longitude.values, 360))
    tw_0360 = pipeline.weekly_interpolated(tn_0360, starts, 't2m')
    assert np.allclose(tw.values, tw_0360.values, equal_nan=True)
    assert tw.shape == (5, 12, 58, 136) and pw.shape == tw.shape
    assert np.isclose(float(pw.isel(week=0, member=0, latitude=0, longitude=0)), 14, atol=.1)
    era5 = pipeline.load_era5_weekly_percentiles()
    era5_bounds = pipeline.thresholds_for_forecast_weeks(era5, starts)
    city_data = pipeline.city_probabilities(tw, pw, indices, era5_bounds, starts)
    assert len(city_data) == 151
    for rec in city_data['Bologna'].values():
        assert abs(sum(rec[f'T_{i}'] for i in range(1, 6)) - 100) < .001
        assert abs(sum(rec[f'P_{i}'] for i in range(1, 6)) - 100) < .001
    _, _, arrays, _, _ = next(pipeline.build_map_probabilities(tw, pw, era5_bounds, starts))
    assert arrays['2tp'].shape == (13, 58, 136)
    assert np.allclose(arrays['2tp'][:5].sum(axis=0), 100, atol=.01)
    test_gridded_percentile_classification(tw, pw, indices, thresholds, starts, era5, era5_bounds, city_data, arrays)
    test_plot_webp_and_upload(city_data, thresholds)
    print('TEST OFFLINE PASS: date, ECDS API request, GRIB maths, percentili ERA5 grigliati, WebP, R2 mock')


if __name__ == '__main__':
    test_all()
