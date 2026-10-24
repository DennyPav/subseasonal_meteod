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
    for variable, api_name in [('t2m', '2_m_temperature'), ('tp', 'total_precipitation')]:
        req = pipeline.request_for(datetime(2026, 10, 22, tzinfo=timezone.utc), variable)
        assert req['variable'] == [api_name]
        assert req['origin'] == 'ecmwf'
        assert req['forecast_type'] == 'perturbed_forecast'
        assert req['data_format'] == 'grib'
        assert req['area'] == [54, -14, 28, 44]
        assert req['year'] == '2026' and req['month'] == '10' and req['day'] == '22'
        assert req['time'] == '00:00' and req['level_type'] == 'single_level'
        assert len(req['leadtime_hour']) == 46
        assert 'class' not in req and 'param' not in req and 'step' not in req
    assert pipeline.request_for(datetime(2026, 10, 22, tzinfo=timezone.utc), 't2m')['leadtime_hour'][:2] == ['0_24','24_48']
    assert pipeline.request_for(datetime(2026, 10, 22, tzinfo=timezone.utc), 'tp')['leadtime_hour'][:2] == ['24','48']


def test_cds_api_download_mock():
    """The ECDS client must use the official cdsapi retrieve/result/download flow."""
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
        files = pipeline.download_grib(datetime(2026,10,22,tzinfo=timezone.utc), Path(td))
        assert set(files) == {'t2m', 'tp'}
        assert len(requests) == 2
        assert all(f.exists() and f.read_bytes()[:4] == b'GRIB' for f in files.values())
        assert requests[0]['variable'] == ['2_m_temperature']
        assert requests[1]['variable'] == ['total_precipitation']


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


def test_all():
    test_dates_and_ecd_api()
    test_cds_api_download_mock()
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
    city_data = pipeline.city_probabilities(tw, pw, indices, thresholds, starts)
    assert len(city_data) == 151
    for rec in city_data['Bologna'].values():
        assert abs(sum(rec[f'T_{i}'] for i in range(1, 6)) - 100) < .001
        assert abs(sum(rec[f'P_{i}'] for i in range(1, 6)) - 100) < .001
    _, _, arrays, _, _ = next(pipeline.build_map_probabilities(tw, pw, thresholds, indices, starts))
    assert arrays['2tp'].shape == (13, 58, 136)
    assert np.allclose(arrays['2tp'][:5].sum(axis=0), 100, atol=.01)
    test_plot_webp_and_upload(city_data, thresholds)
    print('TEST OFFLINE PASS: date, ECDS API request, GRIB maths, quintiles, WebP, R2 mock')


if __name__ == '__main__':
    test_all()
