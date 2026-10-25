#!/usr/bin/env python3
"""MeteoD S2S 0.4°: ECDS ECMWF -> ensemble weekly -> WebP -> Cloudflare R2.

ERA5 gridded weekly percentiles (P20/P40/P60/P80) are the sole quintile
classification thresholds for both maps and city histograms. Historical
city thresholds are retained ONLY as displayed boundary labels in charts.
"""
import argparse
import calendar
import gzip
import json
import os
import tempfile
import zipfile
from pathlib import Path
from datetime import datetime, timedelta, timezone

import h5py
import numpy as np
import pandas as pd
import xarray as xr
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import distance_transform_edt

BASE = Path(__file__).resolve().parent
DATA_DIR = BASE / 'assets'
TEMPERATURE_BIAS_FILE = DATA_DIR / 'bias_s2s_minus_era5_2015_2020_04_smooth.nc'
ERA5_PERCENTILES_FILE = DATA_DIR / 'era5_week_percentiles_2006_2025_conservative_04_smooth.nc'
LATS = np.round(52.0 - np.arange(58, dtype=np.float64) * 0.4, 6)
LONS = np.round(-12.0 + np.arange(136, dtype=np.float64) * 0.4, 6)
GRID_SHAPE = (len(LATS), len(LONS))
QUANTILES = (20, 40, 60, 80)


def date_from_args(value=None, now=None):
    """Default to the latest 00 UTC ECMWF S2S date permitted by the 48h ECDS delay.

    ECMWF S2S runs DAILY since the 2024 model upgrade: do NOT force Mon/Thu.
    The schedule is at 12 UTC, so the cycle at (today - 2 days) 00 UTC
    has already aged at least 60 hours when requested.
    """
    if value:
        return datetime.strptime(value, '%Y-%m-%d').replace(tzinfo=timezone.utc)
    now = now or datetime.now(timezone.utc)
    if now.tzinfo is None:
        raise ValueError('La data/ora di riferimento deve includere il fuso orario')
    return (now.astimezone(timezone.utc) - timedelta(days=2)).replace(
        hour=0, minute=0, second=0, microsecond=0
    )


def week_starts(init):
    monday = init + timedelta(days=(7-init.weekday()) % 7)
    return [monday + timedelta(days=7*i) for i in range(5)]


def load_assets():
    with open(DATA_DIR / 'province_index_mediterraneo_0.4.json', encoding='utf-8') as f:
        indices = json.load(f)
    with gzip.open(DATA_DIR / 'soglie_settimanali_med.json.gz','rt',encoding='utf-8') as f:
        thresholds = json.load(f)
    for city, (i, j) in indices.items():
        if not (0 <= i < GRID_SHAPE[0] and 0 <= j < GRID_SHAPE[1]):
            raise ValueError(f'Indice non compatibile con la griglia 0.4°: {city} = {[i,j]}')
        if city not in thresholds:
            raise ValueError(f'Soglie mancanti per {city}')
    return indices, thresholds


def check_threshold_weeks(indices, thresholds, starts):
    # Avoid silent skipped cities or truncated forecast if a date is outside the 2026-2030 table.
    wanted = [d.strftime('%Y%m%d') for d in starts]
    missing = [(city, week) for city in indices for week in wanted if week not in thresholds[city]]
    if missing:
        raise ValueError(f'Soglie non disponibili per {len(missing)} località/settimane. Esempio: {missing[:5]}')
    for city in indices:
        for datecode in wanted:
            rec = thresholds[city][datecode]
            for v in ('T','P'):
                a = np.array([rec[f'{v}_p{p}'] for p in (0,20,40,60,80,100)], dtype=float)
                if not np.all(np.isfinite(a)) or not np.all(np.diff(a)>=0):
                    raise ValueError(f'Soglie non valide {city} {datecode} {v}: {a}')


# A single CDS-API selection can contain the two parameters and both kinds
# of lead time: the daily mean-temperature intervals, and the accumulated
# precipitation forecast endpoints. They are then split LOCALLY by GRIB
# metadata so cfgrib can decode the different step types independently.
TEMP_STEPS = [f'{h}_{h+24}' for h in range(0, 1104, 24)]
PRECIP_STEPS = [str(h) for h in range(24, 1105, 24)]


def request_for(init, param=None):
    """One ECDS request for temperature AND precipitation, all 46 daily steps.

    A `param` argument is not supported: callers must not accidentally
    issue a second request or lose one variable from the retrieval.
    """
    if param is not None:
        raise ValueError('Il download è congiunto: chiamare request_for(init) senza parametro')
    return {
        'origin': 'ecmwf',
        'year': init.strftime('%Y'),
        'month': init.strftime('%m'),
        'day': init.strftime('%d'),
        'time': '00:00',
        'level_type': 'single_level',
        'variable': ['2_m_temperature', 'total_precipitation'],
        'forecast_type': 'perturbed_forecast',
        # Keep each source product's native step representation.
        # Numeric cumulative endpoints are not interchangeable with
        # the temperature intervals (0_24, 24_48, ...).
        'leadtime_hour': TEMP_STEPS + PRECIP_STEPS,
        'data_format': 'grib',
        'area': [54, -14, 28, 44],
    }


def download_grib(init, raw_dir):
    """Submit exactly ONE ECDS job; cache its combined GRIB or ZIP result."""
    token = os.getenv('ECDS_API_KEY')
    if not token:
        raise RuntimeError('ECDS_API_KEY non configurata. Usare GitHub Secrets (o variabile ambiente).')
    import cdsapi
    client = cdsapi.Client(url='https://ecds.ecmwf.int/api', key=token)
    raw_dir.mkdir(parents=True, exist_ok=True)
    path = raw_dir / f'ecmwf_s2s_t2m_tp_{init:%Y%m%d}.grib'
    if not path.exists() or path.stat().st_size == 0:
        request = request_for(init)
        print('ECDS download UNICO', init.strftime('%Y%m%d'), flush=True)
        print('  variable =', request['variable'], flush=True)
        print('  leadtime_hour =', len(request['leadtime_hour']),
              f'({len(TEMP_STEPS)} temperatura + {len(PRECIP_STEPS)} precipitazione)', flush=True)
        tmp = path.with_suffix('.grib.partial')
        if tmp.exists():
            tmp.unlink()
        # One request and one result, even though it contains two variables.
        client.retrieve('s2s-forecasts', request).download(str(tmp))
        if not tmp.exists() or tmp.stat().st_size < 100:
            raise RuntimeError(f'Download ECDS vuoto: {tmp}')
        with tmp.open('rb') as f:
            magic = f.read(4)
        if magic != b'GRIB' and not zipfile.is_zipfile(tmp):
            raise RuntimeError(
                f'La risposta ECDS non contiene GRIB/ZIP (magic={magic!r}): {tmp}. '
                'Controlla log ECDS e data_format.'
            )
        tmp.replace(path)
    print('DOWNLOAD', path, 'MB=', round(path.stat().st_size / 1e6, 1), flush=True)
    return path


def _grib_product(gid, codes_get):
    """Classify GRIB messages by physical parameter, not file order.

    Temperature must be a 2-m MEAN (not tmax/tmin): otherwise subtracting
    an ERA5 t2m climatological bias would be physically incorrect.
    """
    def field(name, default=None):
        try:
            return codes_get(gid, name)
        except Exception:
            return default

    short = str(field('shortName', '')).lower()
    name = str(field('name', '')).lower()
    stat = str(field('stepType', '')).lower()
    if short in ('mx2t24', 'mn2t24', 'mx2t6', 'mn2t6') or stat in ('max', 'min'):
        if '2t' in short or 'temperature' in name:
            raise ValueError(f'La temperatura scaricata non è media: {short}, stepType={stat}')

    discipline = field('discipline')
    category = field('parameterCategory')
    number = field('parameterNumber')
    param_id = field('paramId')
    type_of_level = str(field('typeOfLevel', ''))
    level = field('level')
    at_2m = (type_of_level == 'heightAboveGround' and str(level) in ('2', '2.0'))

    # ECMWF's mean2t24 is a daily averaged 2m temperature; depending on
    # encoding it may retain generic shortName '2t'.
    is_temp = short in ('mean2t24', '2t', 't2m', '2_m_temperature') or (
        discipline == 0 and category == 0 and number == 0 and at_2m
    )
    if is_temp:
        if short in ('2t', 't2m', '2_m_temperature') and not at_2m:
            raise ValueError(f'Temperatura non a 2 m: {short}, level={type_of_level}:{level}')
        # For generic instantaneous 2t the 24-hour mean is NOT guaranteed.
        if stat not in ('avg', 'average', 'mean'):
            raise ValueError(f'S2S: attesa temperatura MEDIA giornaliera, trovato {short}, stepType={stat}')
        return 't2m'

    is_tp = short in ('tp', 'total_precipitation') or (
        discipline == 0 and category == 1 and number == 8
    )
    if is_tp:
        if stat not in ('accum', 'accumulation', 'sum'):
            raise ValueError(f'S2S: precipitazione non accumulata: {short}, stepType={stat}')
        return 'tp'

    # Do not silently drop unexpected GRIB messages from a mixed request.
    raise ValueError(f'Parametro GRIB inatteso: shortName={short}, paramId={param_id}, '
                     f'parameterCategory={category}, parameterNumber={number}, '
                     f'stepType={stat}')


def split_combined_grib(source, raw_dir, init):
    """Split mixed GRIB into two parameter-specific streams locally.

    This avoids the multiple-hypercube/open_dataset errors that cfgrib may
    produce when averaged temperature and accumulated precipitation have
    different forecast-step keys. GRIB messages are copied without changes.
    Supports a single GRIB stream and a ZIP containing one/more GRIB files.
    """
    from eccodes import (codes_grib_new_from_file, codes_get, codes_get_message,
                         codes_release)

    source = Path(source)
    raw_dir = Path(raw_dir)
    raw_dir.mkdir(parents=True, exist_ok=True)
    paths = {key: raw_dir / f'ecmwf_s2s_{key}_{init:%Y%m%d}.grib'
             for key in ('t2m', 'tp')}
    tmps = {key: path.with_suffix('.grib.partial') for key, path in paths.items()}
    for path in tmps.values():
        path.unlink(missing_ok=True)
    counts = {'t2m': 0, 'tp': 0}

    def process_grib(fileobj, outfiles):
        while True:
            gid = codes_grib_new_from_file(fileobj)
            if gid is None:
                break
            try:
                var = _grib_product(gid, codes_get)
                outfiles[var].write(codes_get_message(gid))
                counts[var] += 1
            finally:
                codes_release(gid)

    try:
        with tmps['t2m'].open('wb') as temp_out, tmps['tp'].open('wb') as precip_out:
            streams = {'t2m': temp_out, 'tp': precip_out}
            if zipfile.is_zipfile(source):
                # Some CDS-API backends package GRIB outputs as a ZIP.
                with zipfile.ZipFile(source) as archive:
                    members = [item for item in archive.infolist() if not item.is_dir()]
                    if not members:
                        raise ValueError(f'ZIP ECDS vuoto: {source}')
                    # A ZipExtFile is not guaranteed to support the seek/fileno
                    # operations expected by ecCodes; unpack each member to
                    # a temporary file, then delete it immediately.
                    import shutil
                    grib_members = 0
                    for i, item in enumerate(members):
                        scratch = raw_dir / f'.ecds_archive_part_{i}.grib'
                        try:
                            with archive.open(item) as src, scratch.open('wb') as dst:
                                shutil.copyfileobj(src, dst, length=1024*1024)
                            with scratch.open('rb') as f:
                                if f.read(4) != b'GRIB':
                                    # ZIP responses may contain a manifest/text file.
                                    continue
                                grib_members += 1
                                f.seek(0)
                                process_grib(f, streams)
                        finally:
                            scratch.unlink(missing_ok=True)
                    if grib_members == 0:
                        raise ValueError(f'ZIP ECDS senza contenuto GRIB: {source}')
            else:
                with source.open('rb') as f:
                    if f.read(4) != b'GRIB':
                        raise ValueError(f'File combinato non GRIB: {source}')
                    f.seek(0)
                    process_grib(f, streams)
        if not all(counts.values()):
            raise ValueError(f'Download congiunto incompleto: messaggi {counts}. '
                             'Verificare la selezione dei leadtime_hour nel Data Store.')
        for var in ('t2m', 'tp'):
            tmps[var].replace(paths[var])
    except Exception:
        for path in tmps.values():
            path.unlink(missing_ok=True)
        for path in paths.values():
            path.unlink(missing_ok=True)
        raise

    print(f'SPLIT LOCALE GRIB: temperatura={counts["t2m"]} messaggi; '
          f'precipitazione={counts["tp"]} messaggi', flush=True)
    return paths


def read_grib_da(path):
    # Each downloaded file contains exactly one GRIB parameter.
    ds=xr.open_dataset(path, engine='cfgrib',backend_kwargs={'indexpath':''})
    names=list(ds.data_vars)
    if len(names)!=1:
        raise RuntimeError(f'{path}: attesa una variabile, ricevute {names}')
    da=ds[names[0]]
    ren={}
    for src,dst in [('number','member'),('ens','member'),('latitude','latitude'),('lat','latitude'),('lon','longitude')]:
        if src in da.dims and src!=dst:ren[src]=dst
    if ren:da=da.rename(ren)
    if 'member' not in da.dims or da.sizes['member']<2:
        raise RuntimeError(f'{path}: ensemble perturbato non trovato! dims={da.dims}')
    required=('step','latitude','longitude')
    if not all(d in da.dims for d in required):
        raise RuntimeError(f'{path}: dimensioni GRIB non attese {da.dims}; richieste {required}')
    da=da.transpose('member','step','latitude','longitude')
    if da.sizes['member'] < 10:
        raise RuntimeError(f'{path}: soltanto {da.sizes["member"]} membri, download apparentemente incompleto')
    # Loaded into memory so the GRIB file is not kept open while xarray operates.
    da=da.load()
    ds.close()
    return da


def step_end_hours(da):
    a=da['step'].values
    if np.issubdtype(a.dtype,np.timedelta64):
        hours=(a/np.timedelta64(1,'h')).astype(float)
    elif np.issubdtype(a.dtype,np.number):
        hours=a.astype(float)
    else:
        raise ValueError(f'Tipo step inaspettato: {a.dtype}')
    if not np.all(np.isclose(hours, np.round(hours))):
        raise ValueError(f'Step non interi in ore: {hours}')
    return np.round(hours).astype(int)


def normalize_da(da, init, var):
    hours=step_end_hours(da)
    expected=np.arange(24,1105,24,dtype=int)
    if len(hours)!=len(expected) or not np.array_equal(hours,expected):
        raise RuntimeError(f'{var}: step incompleti/errati, trovati {len(hours)}: {hours.tolist()}')
    data=da.astype('float32')
    units=str(da.attrs.get('units','')).strip().lower().replace(' ','')
    if var=='t2m':
        if units not in ('k','kelvin'):
            raise ValueError(f'Temperatura con unità non riconosciute: {da.attrs.get("units")}')
        data=data-273.15
    else:
        # These values are cumulative from initialization. Each 24-hour amount is a difference.
        # S2S uses kg m-2, numerically equivalent to millimetres; some ecCodes releases may expose m.
        if units in ('kgm**-2','kgm-2','kg/m2','kg/m^2','mm'):
            factor=1.0
        elif units in ('m','metre','meter'):
            factor=1000.0
        else:
            raise ValueError(f'Unità precipitazione GRIB sconosciute: {da.attrs.get("units")}')
        data=data*factor
        # keep the true first accumulation [0,24] and then consecutive increments
        diff=xr.concat([data.isel(step=slice(0,1)), data.diff('step')],dim='step')
        diff=diff.assign_coords(step=data['step'])
        minimum=float(diff.min())
        if minimum < -0.1:
            raise ValueError(f'Precipitazione cumulata in diminuzione: minimo incremento {minimum:.3f} mm')
        data=diff.clip(min=0.0)
    # Each field for step ending 24,48,... is assigned to its interval START at 00 UTC.
    data=data.assign_coords(day=('step', [np.datetime64((init+timedelta(hours=int(h)-24)).replace(tzinfo=None)) for h in hours]))
    data=data.swap_dims({'step':'day'}).drop_vars('step',errors='ignore')
    return data.rename(var)


def weekly_interpolated(da, starts, var):
    weeks=[]
    for start in starts:
        days=np.array([(start+timedelta(days=i)).replace(tzinfo=None) for i in range(7)],dtype='datetime64[ns]')
        chunk=da.sel(day=days)
        if chunk.sizes['day'] != 7:
            raise RuntimeError(f'{var}: settimana {start:%Y%m%d} incompleta')
        agg=chunk.mean('day') if var=='t2m' else chunk.sum('day')
        weeks.append(agg)
    result=xr.concat(weeks,dim='week').assign_coords(week=np.arange(1,6))
    # GRIB may use 0..360 longitude coordinates. Convert to -180..180 before
    # interpolating across the Greenwich meridian (our target spans -12..42).
    lon = ((result['longitude'].values + 180.0) % 360.0) - 180.0
    result = result.assign_coords(longitude=lon).sortby('longitude')
    if np.unique(lon).size != lon.size:
        raise RuntimeError(f'{var}: longitude duplicate after 0..360 conversion')
    result=result.interp(latitude=LATS,longitude=LONS,method='linear',assume_sorted=False).astype('float32')
    if not np.all(np.isfinite(result.values)):
        raise ValueError(f'{var}: valori NaN dopo interpolazione 0.4°. Verificare area ECDS.')
    return result.transpose('week','member','latitude','longitude')



def climatological_week_for_day(day):
    """Index 0..51: Jan 1-7, Jan 8-14, ... on a no-leap calendar.

    Feb 29 takes Feb 28's bin; Dec 31 takes week 52 because the
    52 seven-day bins cover only 364 calendar days. This matches the
    calendar used to build the attached S2S/ERA5 bias climatology.
    """
    if isinstance(day, datetime):
        day = day.date()
    day_of_year = (day - day.replace(month=1, day=1)).days
    if calendar.isleap(day.year):
        if day.month == 2 and day.day == 29:
            day_of_year -= 1  # same 0-based index as Feb 28
        elif (day.month, day.day) > (2, 29):
            day_of_year -= 1
    return min(day_of_year // 7, 51)


def load_temperature_bias():
    """Read 52 weekly S2S-minus-ERA5 maps and resample to pipeline grid.

    The NetCDF is cell-centred (58x135), whereas our forecast maps are
    point-centred (58x136). Interpolation is bilinear. In the supplied
    reference, the northern/southern three rows have missing bias: we
    extend the closest valid bias there rather than leaving ensemble
    temperature uncorrected or allowing missing values into probabilities.
    """
    if not TEMPERATURE_BIAS_FILE.is_file():
        raise FileNotFoundError(f'File bias temperatura mancante: {TEMPERATURE_BIAS_FILE}')

    with h5py.File(TEMPERATURE_BIAS_FILE, 'r') as nc:
        expected_vars = ('bias', 'week', 'latitude', 'longitude')
        missing = [key for key in expected_vars if key not in nc]
        if missing:
            raise ValueError(f'NetCDF bias incompleto; variabili mancanti: {missing}')
        weeks = np.asarray(nc['week'][:])
        source_lat = np.asarray(nc['latitude'][:], dtype=np.float64)
        source_lon = np.asarray(nc['longitude'][:], dtype=np.float64)
        values = np.asarray(nc['bias'][:], dtype=np.float64)
        if values.shape != (52, len(source_lat), len(source_lon)):
            raise ValueError(f'Dimensioni bias inattese: {values.shape}')
        if not np.array_equal(weeks, np.arange(1, 53)):
            raise ValueError('La coordinata week del bias deve essere 1..52')
        if not (np.all(np.diff(source_lat) < 0) or np.all(np.diff(source_lat) > 0)):
            raise ValueError('Latitudini bias non monotone')
        if not (np.all(np.diff(source_lon) < 0) or np.all(np.diff(source_lon) > 0)):
            raise ValueError('Longitudini bias non monotone')
        units = nc['bias'].attrs.get('units', '')
        if isinstance(units, bytes):
            units = units.decode('utf-8')
        if str(units).strip().lower() not in ('degc', 'c', 'celsius', 'degree_celsius'):
            raise ValueError(f'Unita bias non attese: {units!r} (servono gradi Celsius)')
        definition = nc.attrs.get('bias_definition', '')
        if isinstance(definition, bytes):
            definition = definition.decode('utf-8')
        if str(definition).strip() != 'S2S t2m - ERA5 t2m':
            raise ValueError(f'Definizione bias non compatibile con sottrazione: {definition!r}')

    # Sort source axes, ensuring strictly ascending coordinates for scipy.
    lat_order = np.argsort(source_lat)
    lon_order = np.argsort(source_lon)
    source_lat = source_lat[lat_order]
    source_lon = source_lon[lon_order]
    if np.any(np.diff(source_lat) <= 0) or np.any(np.diff(source_lon) <= 0):
        raise ValueError('Coordinate bias duplicate')
    values = values[:, lat_order, :][:, :, lon_order]

    # Allow just the half-cell overhang between target cell-centres and
    # forecast points; never project a bias across an unrelated domain.
    step = 0.4
    if (LATS.max() > source_lat.max() + 2 * step
        or LATS.min() < source_lat.min() - 2 * step
        or LONS.max() > source_lon.max() + 2 * step
        or LONS.min() < source_lon.min() - 2 * step):
        raise ValueError('Griglia bias non compatibile con il dominio delle previsioni S2S')

    lat_mesh, lon_mesh = np.meshgrid(
        np.clip(LATS, source_lat[0], source_lat[-1]),
        np.clip(LONS, source_lon[0], source_lon[-1]), indexing='ij'
    )
    targets = np.stack((lat_mesh.ravel(), lon_mesh.ravel()), axis=-1)

    correction = np.empty((52, *GRID_SHAPE), dtype=np.float32)
    missing_count = 0
    for i in range(52):
        field = values[i]
        finite = np.isfinite(field)
        if not np.any(finite):
            raise ValueError(f'Settimana bias {i + 1}: nessun pixel valido')
        missing_count += int((~finite).sum())
        if not np.all(finite):
            # Only fill spatial holes in the bias field. In particular the
            # provided dataset has exactly three empty northern and three
            # empty southern rows, common to every week.
            nearest_indices = distance_transform_edt(~finite, return_distances=False,
                                                     return_indices=True)
            field = field[tuple(nearest_indices)]
        interpolator = RegularGridInterpolator(
            (source_lat, source_lon), field,
            method='linear', bounds_error=True
        )
        correction[i] = interpolator(targets).reshape(GRID_SHAPE).astype(np.float32)

    if not np.isfinite(correction).all():
        raise ValueError('Remapping bias ha prodotto valori non finiti')
    print(f'BIAS TEMPERATURA: {TEMPERATURE_BIAS_FILE.name}; '
          f'{len(correction)} settimane, griglia {GRID_SHAPE}, '
          f'pixel mancanti sorgente {missing_count}/{values.size} '
          '(riempiti dal pixel valido più vicino prima del remapping)', flush=True)
    return correction


def correct_ensemble_temperature(tweek, starts, weekly_bias):
    """Subtract S2S-ERA5 bias at every pixel of EVERY ensemble member.

    Forecast weeks are Monday-Sunday, unlike the 7-day fixed climatological
    bins from January 1. Average the seven daily corresponding bias maps
    to give the correct bias for the whole forecast week, including weeks
    crossing a climatological bin or year boundary.
    """
    if tweek.sizes.get('week') != len(starts):
        raise ValueError('Numero settimane previsionali incompatibile con le date')
    if not np.allclose(tweek.latitude.values, LATS, atol=1e-6, rtol=0):
        raise ValueError('La griglia latitudine della previsione non corrisponde al bias')
    if not np.allclose(tweek.longitude.values, LONS, atol=1e-6, rtol=0):
        raise ValueError('La griglia longitudine della previsione non corrisponde al bias')
    if weekly_bias.shape != (52, *GRID_SHAPE) or not np.isfinite(weekly_bias).all():
        raise ValueError('Campo bias 52 settimane non valido')

    corrections = np.stack([
        np.mean([weekly_bias[climatological_week_for_day(start + timedelta(days=d))]
                 for d in range(7)], axis=0)
        for start in starts
    ], axis=0).astype(np.float32)  # (week, latitude, longitude)

    # Subtraction broadcasts to all members, never changes their spread.
    corrected = (tweek.values - corrections[:, None, :, :]).astype(np.float32)
    output = tweek.copy(data=corrected)
    output.attrs.update({
        'temperature_bias_correction': 'T_corrected = T_S2S - (S2S t2m - ERA5 t2m)',
        'temperature_bias_file': TEMPERATURE_BIAS_FILE.name,
        'bias_time_mapping': '7 daily no-leap climatological week indices per forecast week',
    })
    for n, start in enumerate(starts):
        bins = [climatological_week_for_day(start + timedelta(days=d)) + 1
                for d in range(7)]
        print(f'  Bias W{n+1} ({start:%Y-%m-%d}): '
              f'settimane climatologiche {bins}, '
              f'media bias {np.mean(corrections[n]):+.2f} °C', flush=True)
    return output


def load_era5_weekly_percentiles():
    """Read ERA5 per-pixel P20/P40/P60/P80 and align to the forecast's nodes.

    Input grid: 58 x 135 cell centres (longitude -11.8 ... 41.8).
    Forecast:   58 x 136 sample nodes (longitude -12.0 ... 42.0).
    Edge nodes are clamped by <= one half-cell to the nearest valid centre;
    internal nodes are bilinearly interpolated. No city values are used.

    Returns {'T': (52,4,58,136), 'P': (52,4,58,136)} float32.
    """
    if not ERA5_PERCENTILES_FILE.is_file():
        raise FileNotFoundError(f'NetCDF percentili ERA5 mancante: {ERA5_PERCENTILES_FILE}')

    with h5py.File(ERA5_PERCENTILES_FILE, 'r') as nc:
        variable_names = {
            'T': [f't2m_p{q}' for q in QUANTILES],
            'P': [f'tp_p{q}' for q in QUANTILES],
        }
        required = {'latitude', 'longitude', 'week',
                    *(name for names in variable_names.values() for name in names)}
        missing = sorted(required - set(nc.keys()))
        if missing:
            raise ValueError(f'File percentili incompleto: {missing}')
        weeks = np.asarray(nc['week'][:])
        if not np.array_equal(weeks, np.arange(1, 53)):
            raise ValueError('Percentili ERA5: week deve essere 1..52')
        src_lat = np.asarray(nc['latitude'][:], dtype=np.float64)
        src_lon = np.asarray(nc['longitude'][:], dtype=np.float64)
        if src_lat.ndim != 1 or src_lon.ndim != 1:
            raise ValueError('Percentili ERA5: coordinate latitude/longitude devono essere 1D')
        if len(src_lat) < 2 or len(src_lon) < 2:
            raise ValueError('Percentili ERA5: griglia troppo piccola')
        order_lat = np.argsort(src_lat)
        order_lon = np.argsort(src_lon)
        src_lat = src_lat[order_lat]
        src_lon = src_lon[order_lon]
        if np.any(np.diff(src_lat) <= 0) or np.any(np.diff(src_lon) <= 0):
            raise ValueError('Percentili ERA5: coordinate duplicate/non monotone')
        if (LATS.max() > src_lat.max() + 0.4 or LATS.min() < src_lat.min() - 0.4
            or LONS.max() > src_lon.max() + 0.4 or LONS.min() < src_lon.min() - 0.4):
            raise ValueError('Percentili ERA5: griglia geografica non compatibile con le mappe')
        target_lat, target_lon = np.meshgrid(
            np.clip(LATS, src_lat[0], src_lat[-1]),
            np.clip(LONS, src_lon[0], src_lon[-1]), indexing='ij'
        )
        targets = np.stack((target_lat.ravel(), target_lon.ravel()), axis=-1)
        output = {}
        for kind, names in variable_names.items():
            cubes = []
            for name in names:
                var = nc[name]
                units = var.attrs.get('units', b'')
                if isinstance(units, bytes):
                    units = units.decode('utf-8')
                expected = ('degc', 'c', 'celsius', 'degree_celsius') if kind == 'T' else ('mm',)
                if str(units).strip().lower() not in expected:
                    raise ValueError(f'Percentili ERA5: unità inattesa {name}={units!r}')
                values = np.asarray(var[:], dtype=np.float32)
                if values.shape != (52, len(order_lat), len(order_lon)):
                    raise ValueError(f'{name}: shape inattesa {values.shape}')
                values = values[:, order_lat, :][:, :, order_lon]
                if not np.isfinite(values).all():
                    raise ValueError(f'{name}: soglie ERA5 non finite (NaN): non uso fallback cittadino')
                dest = np.empty((52, *GRID_SHAPE), dtype=np.float32)
                for wi, field in enumerate(values):
                    interp = RegularGridInterpolator(
                        (src_lat, src_lon), field, method='linear', bounds_error=True
                    )
                    dest[wi] = interp(targets).reshape(GRID_SHAPE).astype(np.float32)
                cubes.append(dest)
            cube = np.stack(cubes, axis=1)  # 52 weeks, 4 thresholds, lat, lon
            if (np.diff(cube, axis=1) < -1e-4).any():
                raise ValueError(f'Percentili ERA5 {kind}: P20/P40/P60/P80 non ordinati')
            if not np.isfinite(cube).all():
                raise ValueError(f'Percentili ERA5 {kind}: problemi nel remapping')
            output[kind] = cube

    print(f'SOGLIE QUINTILI ERA5: {ERA5_PERCENTILES_FILE.name}; '
          f'52 settimane, T e P, 4 percentili, griglia {GRID_SHAPE}; '
          'estremi spaziali limitati al centro-cella più vicino', flush=True)
    return output


def thresholds_for_forecast_weeks(weekly_percentiles, starts):
    """Average daily-matched climatological boundaries across each forecast week.

    For a Mon-Sun forecast week straddling fixed Jan-1 ERA5 week boundaries,
    a 7-day weighted mixture of weekly percentile MAPS is used, just as for
    the temperature bias. This is an approximation to the percentile of the
    shifted 7-day window, which cannot be recovered from weekly quantiles alone.

    Returns (nforecastweeks, 4 thresholds, lat, lon) for T and P.
    """
    indices = np.array([
        [climatological_week_for_day(start + timedelta(days=d)) for d in range(7)]
        for start in starts
    ], dtype=np.int16)
    out = {}
    for kind in ('T', 'P'):
        source = np.asarray(weekly_percentiles[kind])
        if source.shape != (52, 4, *GRID_SHAPE) or not np.isfinite(source).all():
            raise ValueError(f'Percentili ERA5 {kind}: shape o valori non validi')
        weekly = source[indices].mean(axis=1).astype(np.float32)
        if (np.diff(weekly, axis=1) < -1e-4).any():
            raise ValueError(f'Percentili ERA5 {kind}: soglie risultanti non ordinate')
        out[kind] = weekly
    return out


def quintile_probs(fc, bounds):
    """fc: (members, lat, lon) or (members,), bounds: (4, lat,lon) or (4,).
       Outputs: (5,lat,lon) or (5,); equality goes to the lower quintile.
    """
    fc=np.asarray(fc)
    bounds=np.asarray(bounds)
    if fc.ndim==1:
        classes=np.sum(fc[:,None]>bounds[None,:],axis=1)
    else:
        classes=np.sum(fc[:,None,:,:]>bounds[None,:,:,:],axis=1)
    return np.stack([100.0*np.mean(classes==k,axis=0) for k in range(5)],axis=0).astype('float32')


def rounded_proba(vals):
    a=np.asarray(vals,dtype=float)
    r=np.round(a,2)
    r[-1]=np.round(100-r[:-1].sum(),2)
    return r.tolist()


def city_probabilities(tweek, pweek, indices, era5_bounds, starts):
    """Classify ensembles with the exact ERA5 pixel thresholds of each city.

    The legacy JSON city thresholds do NOT enter probability calculations:
    they are passed only to plot_city_combined() for boundary labels.
    """
    out = {}
    for city, (row, col) in indices.items():
        perweek = {}
        for n, start in enumerate(starts[:4]):
            key = start.strftime('%Y%m%d')
            arrs = {}
            for name, weather in (('T', tweek), ('P', pweek)):
                series = weather.isel(week=n, latitude=row, longitude=col).values
                bounds = era5_bounds[name][n, :, row, col]
                probs = quintile_probs(series, bounds)
                arrs.update({f'{name}_{k+1}': v for k, v in enumerate(rounded_proba(probs))})
            perweek[key] = arrs
        out[city] = perweek
    return out


def build_map_probabilities(tweek, pweek, era5_bounds, starts):
    """Map probabilities with gridded ERA5 thresholds, NEVER city interpolation."""
    lon2d, lat2d = np.meshgrid(LONS, LATS)
    for n, start in enumerate(starts):
        t_bounds = era5_bounds['T'][n]
        p_bounds = era5_bounds['P'][n]
        tprob = quintile_probs(tweek.isel(week=n).values, t_bounds)
        pprob = quintile_probs(pweek.isel(week=n).values, p_bounds)
        arrays = {}
        for field, probs in (('2tp', tprob), ('tpp', pprob)):
            arr = np.full((13, *GRID_SHAPE), np.nan, dtype=np.float32)
            arr[:5] = probs
            arr[10] = probs[0] + probs[1]
            arr[11] = probs[2]
            arr[12] = probs[3] + probs[4]
            arrays[field] = arr
        yield n, start, arrays, lat2d, lon2d


def render_histograms(city_json,thresholds,init,folder):
    from plotting import plot_city_combined
    folder.mkdir(exist_ok=True,parents=True)
    for city, weeks in city_json.items():
        starts=sorted(weeks.keys())
        out=folder / f'{city}_{starts[0]}.webp'
        plot_city_combined(city,starts,weeks,thresholds[city],str(out),init.strftime('%d/%m/%Y'))
    print(f'Grafici località prodotti: {len(city_json)}',flush=True)


def render_maps(tweek,pweek,era5_bounds,starts,init,folder):
    import plotting
    plotting.TARGET_DATE=init.replace(tzinfo=None)
    plotting.target_date=init.replace(tzinfo=None)
    folder.mkdir(exist_ok=True,parents=True)
    for n,start,arrays,lat2d,lon2d in build_map_probabilities(tweek,pweek,era5_bounds,starts):
        # Legacy plot infers start by reading final four chars prior to ____, end-exclusive.
        end_exclusive=start+timedelta(days=7)
        key=f'S2F{init:%m%d}{end_exclusive:%m%d}____ECDS'
        for v in ('2tp','tpp'):
            print('Plot',start.strftime('%Y%m%d'),v,flush=True)
            plotting.plot_terciles_week(key,arrays[v],v,str(folder),lat2d,lon2d)


def validate_r2_config():
    """Check credentials before initiating an expensive ECDS download."""
    names = ('R2_ACCOUNT_ID', 'R2_ACCESS_KEY_ID', 'R2_SECRET_ACCESS_KEY')
    absent = [n for n in names if not os.environ.get(n)]
    if absent:
        raise RuntimeError('Credenziali Cloudflare R2 mancanti: ' + ', '.join(absent) +
                           '. Configurarle nei GitHub Actions Secrets.')


def upload_webps_to_r2(image_dir, expected_count, client=None):
    """Publish *only* WebP images to the original `settimanali` bucket root.

    The previous Leonardo workflow used `upload_file(..., bucket, filename)`;
    the object key was the bare filename, NOT a prefix such as `PLOT/`.
    An explicit R2_PREFIX can override this for a differently configured site.
    """
    from pathlib import Path
    image_dir = Path(image_dir)
    images = sorted(image_dir.rglob('*.webp'))
    if len(images) != expected_count:
        raise RuntimeError(f'WebP prodotti: {len(images)}, attesi: {expected_count}. Upload annullato.')
    names = [p.name for p in images]
    if len(set(names)) != len(names):
        raise RuntimeError('Nomi WebP duplicati fra mappe e località: upload annullato.')
    if any(p.stat().st_size < 100 for p in images):
        raise RuntimeError('WebP vuoti o incompleti: upload annullato.')
    validate_r2_config()

    bucket = os.getenv('R2_BUCKET', 'settimanali').strip()
    if not bucket:
        raise RuntimeError('R2_BUCKET vuoto')
    prefix = os.getenv('R2_PREFIX', '').strip().strip('/')
    if client is None:
        import boto3
        from botocore.config import Config
        client = boto3.client(
            's3', endpoint_url=f"https://{os.environ['R2_ACCOUNT_ID']}.r2.cloudflarestorage.com",
            aws_access_key_id=os.environ['R2_ACCESS_KEY_ID'],
            aws_secret_access_key=os.environ['R2_SECRET_ACCESS_KEY'],
            region_name='auto',
            config=Config(retries={'max_attempts': 5, 'mode': 'standard'}),
        )
    extra = {'ContentType': 'image/webp', 'CacheControl': 'public, max-age=300'}
    print(f'UPLOAD R2: {len(images)} WebP nel bucket {bucket}, prefix={prefix or "(radice)"}', flush=True)
    objects = []
    for i, path in enumerate(images, 1):
        key = f'{prefix}/{path.name}' if prefix else path.name
        # Exceptions must propagate: GitHub must report a failed upload.
        client.upload_file(str(path), bucket, key, ExtraArgs=extra)
        objects.append((key, path.stat().st_size))
        if i % 20 == 0 or i == len(images):
            print(f'Caricati {i}/{len(images)} WebP', flush=True)

    # Optional remote verification: a successful boto3 upload_file() normally
    # confirms completion, but HEAD provides an independent size check before
    # the GitHub action exits successfully.
    if os.environ.get('R2_VERIFY_OBJECTS', '').lower() in ('1', 'true', 'yes'):
        from concurrent.futures import ThreadPoolExecutor

        def verify_object(entry):
            key, size = entry
            response = client.head_object(Bucket=bucket, Key=key)
            remote_size = int(response['ContentLength'])
            if remote_size != size:
                raise RuntimeError(
                    f'Upload R2 non verificato per {key}: locale={size} B, remoto={remote_size} B'
                )
            return key

        with ThreadPoolExecutor(max_workers=8) as executor:
            confirmed = list(executor.map(verify_object, objects))
        print(f'Verifica R2 completata: {len(confirmed)}/{len(objects)} oggetti presenti con dimensione corretta', flush=True)
    return len(images)


def main(argv=None):
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--date', help='Emissione ECMWF YYYY-MM-DD; omesso = oggi UTC meno 2 giorni')
    cli.add_argument('--output', default='output', help='Directory locale per i SOLI WebP')
    cli.add_argument('--grib-dir', help='GRIB precedentemente scaricati; evita la richiesta ECDS')
    cli.add_argument('--no-maps', action='store_true', help='Genera soltanto gli istogrammi provinciali')
    cli.add_argument('--no-upload', action='store_true', help='Genera WebP senza inviarli a Cloudflare (prove locali)')
    cli.add_argument('--dry-run', action='store_true', help='Valida asset e stampa richieste senza scaricare o caricare')
    args = cli.parse_args(argv)
    init = date_from_args(args.date)
    starts = week_starts(init)
    indices, thresholds = load_assets()
    check_threshold_weeks(indices, thresholds, starts)
    bias_weekly = load_temperature_bias()  # fail before downloading if asset is invalid
    era5_percentiles = load_era5_weekly_percentiles()  # also validates all gridded thresholds
    era5_bounds = thresholds_for_forecast_weeks(era5_percentiles, starts)
    print('RUN ECMWF S2S UTC', init.strftime('%Y-%m-%d'),
          'Grid', GRID_SHAPE, 'Locations', len(indices), flush=True)
    print('Weeks:', ', '.join(s.strftime('%Y-%m-%d') for s in starts), flush=True)
    if args.dry_run:
        print(json.dumps(request_for(init), indent=2))
        return
    if not args.no_upload:
        validate_r2_config()  # fail early if secrets are missing

    # The only persistent products are WebP. Raw GRIB are downloaded to a
    # TemporaryDirectory and deleted automatically even when processing fails.
    image_dir = Path(args.output).resolve() / init.strftime('%Y%m%d')
    image_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='meteod_ecdsv4_') as temp:
        if args.grib_dir:
            raw = Path(args.grib_dir).resolve()
            combined = raw / f'ecmwf_s2s_t2m_tp_{init:%Y%m%d}.grib'
            if combined.is_file():
                # Never alter source files supplied by --grib-dir.
                paths = split_combined_grib(combined, Path(temp), init)
            else:
                # Backward compatibility with old saved GRIB pairs.
                paths = {v: raw / f'ecmwf_s2s_{v}_{init:%Y%m%d}.grib'
                         for v in ('t2m', 'tp')}
                for path in paths.values():
                    if not path.is_file():
                        raise FileNotFoundError(path)
        else:
            combined = download_grib(init, Path(temp))
            paths = split_combined_grib(combined, Path(temp), init)
        t = normalize_da(read_grib_da(paths['t2m']), init, 't2m')
        p = normalize_da(read_grib_da(paths['tp']), init, 'tp')
        if not np.array_equal(t['member'].values, p['member'].values):
            raise RuntimeError('Membri temperatura e precipitazione non allineati')
        tw = weekly_interpolated(t, starts, 't2m')
        pw = weekly_interpolated(p, starts, 'tp')
        del t, p
        # Correct each ensemble member BEFORE both city and map probabilities.
        # Precipitation remains untouched.
        tw = correct_ensemble_temperature(tw, starts, bias_weekly)
        city_data = city_probabilities(tw, pw, indices, era5_bounds, starts)
        # Only the printed values at the quintile boundaries use old city JSON.
        render_histograms(city_data, thresholds, init, image_dir / 'citta')
        if not args.no_maps:
            render_maps(tw, pw, era5_bounds, starts, init, image_dir / 'mappe')
    expected = len(indices) + (0 if args.no_maps else 2 * len(starts))
    if args.no_upload:
        webps = list(image_dir.rglob('*.webp'))
        if len(webps) != expected:
            raise RuntimeError(f'WebP generati: {len(webps)}; attesi: {expected}')
        print(f'DONE (senza upload): {len(webps)} WebP in {image_dir}', flush=True)
    else:
        upload_webps_to_r2(image_dir, expected)
        print(f'DONE: {expected} WebP caricati su R2', flush=True)


if __name__ == '__main__':
    main()
