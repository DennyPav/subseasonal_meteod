#!/usr/bin/env python3
"""MeteoD S2S 0.4°: ECDS ECMWF -> ensemble weekly -> PNG -> Cloudflare R2.

IMPORTANT: spatially interpolated city thresholds used for maps are an APPROXIMATION,
not a gridded ERA5 climatology. City histograms use the ORIGINAL per-city thresholds.
"""
import argparse
import gzip
import json
import os
import tempfile
from pathlib import Path
from datetime import datetime, timedelta, timezone

import numpy as np
import pandas as pd
import xarray as xr
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator
from scipy.spatial import Delaunay

BASE = Path(__file__).resolve().parent
DATA_DIR = BASE / 'assets'
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


def request_for(init, param):
    """Select S2S data using the *official ECDS form fields*, not raw MARS syntax.

    Mirrors the user-supplied `s2s-forecasts` example:
        origin / year / month / day / time / level_type / variable /
        forecast_type / leadtime_hour / data_format / area.

    The S2S daily 2m temperature is archived as intervals `0_24`, etc.
    Total precipitation is cumulative since initialization. Request 24-hour
    endpoints only: subtraction of adjacent endpoints yields daily mm while
    avoiding downloading 6-hourly steps we do not need for weekly products.
    """
    if param == 't2m':
        variable = '2_m_temperature'
        leadtimes = [f'{h}_{h+24}' for h in range(0, 1104, 24)]
    elif param == 'tp':
        variable = 'total_precipitation'
        leadtimes = [str(h) for h in range(24, 1105, 24)]
    else:
        raise ValueError(f'Parametro non previsto: {param}')

    return {
        'origin': 'ecmwf',
        'year': init.strftime('%Y'),
        'month': init.strftime('%m'),
        'day': init.strftime('%d'),
        'time': '00:00',
        'level_type': 'single_level',
        'variable': [variable],
        'forecast_type': 'perturbed_forecast',
        'leadtime_hour': leadtimes,
        'data_format': 'grib',
        # Margin for bilinear interpolation at all target-grid boundaries.
        'area': [54, -14, 28, 44],
    }


def download_grib(init, raw_dir):
    token = os.getenv('ECDS_API_KEY')
    if not token:
        raise RuntimeError('ECDS_API_KEY non configurata. Usare GitHub Secrets (o variabile ambiente).')
    import cdsapi
    client = cdsapi.Client(url='https://ecds.ecmwf.int/api', key=token)
    paths = {}
    raw_dir.mkdir(parents=True, exist_ok=True)
    for var in ('t2m', 'tp'):
        path = raw_dir / f'ecmwf_s2s_{var}_{init:%Y%m%d}.grib'
        if not path.exists() or path.stat().st_size == 0:
            request = request_for(init, var)
            print('ECDS download', var, init.strftime('%Y%m%d'), flush=True)
            print('  variable =', request['variable'], flush=True)
            print('  leadtime_hour =', len(request['leadtime_hour']), 'step', flush=True)
            tmp = path.with_suffix('.grib.partial')
            if tmp.exists():
                tmp.unlink()
            # ECDS CDS-API returns a job/result whose download() retrieves the GRIB.
            client.retrieve('s2s-forecasts', request).download(str(tmp))
            if not tmp.exists() or tmp.stat().st_size < 100:
                raise RuntimeError(f'Download ECDS vuoto: {tmp}')
            with tmp.open('rb') as f:
                magic = f.read(4)
            if magic != b'GRIB':
                raise RuntimeError(
                    f'La risposta ECDS non è un GRIB (magic={magic!r}): {tmp}. '
                    'Controlla log ECDS e data_format.'
                )
            tmp.replace(path)
        print('GRIB', path, 'MB=', round(path.stat().st_size / 1e6, 1), flush=True)
        paths[var] = path
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


def city_probabilities(tweek,pweek, indices, thresholds, starts):
    """Uses EACH city's OWN climatological percentiles without interpolation."""
    out={}
    for city,(row,col) in indices.items():
        perweek={}
        for n,start in enumerate(starts[:4]):
            key=start.strftime('%Y%m%d')
            record=thresholds[city][key]
            arrs={}
            for name,weather,kind in [('T',tweek,'T'),('P',pweek,'P')]:
                series=weather.isel(week=n,latitude=row,longitude=col).values
                bounds=np.array([record[f'{kind}_p{k}'] for k in QUANTILES])
                probs=quintile_probs(series,bounds)
                arrs.update({f'{name}_{k+1}':v for k,v in enumerate(rounded_proba(probs))})
            perweek[key]=arrs
        out[city]=perweek
    return out


class ApproximateMapThresholds:
    """Interpolate climatology *between city points* to build map-only approximations.
       Not a substitute for gridded ERA5 climatology. Duplicate locations are averaged.
    """
    def __init__(self,indices):
        groups={}
        for city,(row,col) in indices.items():
            groups.setdefault((row,col),[]).append(city)
        self.groups=groups
        self.points=np.array([[LONS[j],LATS[i]] for i,j in groups],dtype=float)
        self.lon2d,self.lat2d=np.meshgrid(LONS,LATS)
        self.tri=Delaunay(self.points)
    def bounds(self, city_thresholds, start, kind):
        key=start.strftime('%Y%m%d')
        source=np.array([[np.mean([city_thresholds[city][key][f'{kind}_p{q}'] for city in cities]) for q in QUANTILES] for cities in self.groups.values()],dtype=float)
        linear=LinearNDInterpolator(self.tri,source,fill_value=np.nan)(self.lon2d,self.lat2d)
        if np.isnan(linear).any():
            nearest=NearestNDInterpolator(self.points,source)(self.lon2d,self.lat2d)
            linear=np.where(np.isnan(linear),nearest,linear)
        bounds=np.maximum.accumulate(linear,axis=-1)
        return np.moveaxis(bounds,-1,0) # (4,lat,lon)


def build_map_probabilities(tweek,pweek, thresholds, indices, starts):
    interpolator=ApproximateMapThresholds(indices)
    for n,start in enumerate(starts):
        t_bounds=interpolator.bounds(thresholds,start,'T')
        p_bounds=interpolator.bounds(thresholds,start,'P')
        tprob=quintile_probs(tweek.isel(week=n).values,t_bounds)
        pprob=quintile_probs(pweek.isel(week=n).values,p_bounds)
        # layout exactly as plot_terciles_week() in original Leonardo script:
        # first 5 quintiles, positions 10..12 sum of Q1+Q2, Q3, Q4+Q5.
        arrays={}
        for field, probs in [('2tp',tprob),('tpp',pprob)]:
            arr=np.full((13,*GRID_SHAPE),np.nan,dtype=np.float32)
            arr[:5]=probs
            arr[10]=probs[0]+probs[1]
            arr[11]=probs[2]
            arr[12]=probs[3]+probs[4]
            arrays[field]=arr
        yield n,start,arrays,interpolator.lat2d,interpolator.lon2d


def render_histograms(city_json,thresholds,init,folder):
    from plotting import plot_city_combined
    folder.mkdir(exist_ok=True,parents=True)
    for city, weeks in city_json.items():
        starts=sorted(weeks.keys())
        out=folder / f'{city}_{starts[0]}.png'
        plot_city_combined(city,starts,weeks,thresholds[city],str(out),init.strftime('%d/%m/%Y'))
    print(f'Grafici località prodotti: {len(city_json)}',flush=True)


def render_maps(tweek,pweek,indices,thresholds,starts,init,folder):
    import plotting
    plotting.TARGET_DATE=init.replace(tzinfo=None)
    folder.mkdir(exist_ok=True,parents=True)
    for n,start,arrays,lat2d,lon2d in build_map_probabilities(tweek,pweek,thresholds,indices,starts):
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


def upload_pngs_to_r2(image_dir, expected_count, client=None):
    """Publish *only* PNG images to the original `settimanali` bucket root.

    The previous Leonardo workflow used `upload_file(..., bucket, filename)`;
    the object key was the bare filename, NOT a prefix such as `PLOT/`.
    An explicit R2_PREFIX can override this for a differently configured site.
    """
    from pathlib import Path
    image_dir = Path(image_dir)
    images = sorted(image_dir.rglob('*.png'))
    if len(images) != expected_count:
        raise RuntimeError(f'PNG prodotti: {len(images)}, attesi: {expected_count}. Upload annullato.')
    names = [p.name for p in images]
    if len(set(names)) != len(names):
        raise RuntimeError('Nomi PNG duplicati fra mappe e località: upload annullato.')
    if any(p.stat().st_size < 100 for p in images):
        raise RuntimeError('PNG vuoti o incompleti: upload annullato.')
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
    extra = {'ContentType': 'image/png', 'CacheControl': 'public, max-age=300'}
    print(f'UPLOAD R2: {len(images)} PNG nel bucket {bucket}, prefix={prefix or "(radice)"}', flush=True)
    for i, path in enumerate(images, 1):
        key = f'{prefix}/{path.name}' if prefix else path.name
        # Exceptions must propagate: GitHub must report a failed upload.
        client.upload_file(str(path), bucket, key, ExtraArgs=extra)
        if i % 20 == 0 or i == len(images):
            print(f'Caricati {i}/{len(images)} PNG', flush=True)
    return len(images)


def main(argv=None):
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--date', help='Emissione ECMWF YYYY-MM-DD; omesso = oggi UTC meno 2 giorni')
    cli.add_argument('--output', default='output', help='Directory locale per i SOLI PNG')
    cli.add_argument('--grib-dir', help='GRIB precedentemente scaricati; evita la richiesta ECDS')
    cli.add_argument('--no-maps', action='store_true', help='Genera soltanto gli istogrammi provinciali')
    cli.add_argument('--no-upload', action='store_true', help='Genera PNG senza inviarli a Cloudflare (prove locali)')
    cli.add_argument('--dry-run', action='store_true', help='Valida asset e stampa richieste senza scaricare o caricare')
    args = cli.parse_args(argv)
    init = date_from_args(args.date)
    starts = week_starts(init)
    indices, thresholds = load_assets()
    check_threshold_weeks(indices, thresholds, starts)
    print('RUN ECMWF S2S UTC', init.strftime('%Y-%m-%d'),
          'Grid', GRID_SHAPE, 'Locations', len(indices), flush=True)
    print('Weeks:', ', '.join(s.strftime('%Y-%m-%d') for s in starts), flush=True)
    if args.dry_run:
        print(json.dumps({p: request_for(init, p) for p in ('t2m', 'tp')}, indent=2))
        return
    if not args.no_upload:
        validate_r2_config()  # fail early if secrets are missing

    # The only persistent products are PNG. Raw GRIB are downloaded to a
    # TemporaryDirectory and deleted automatically even when processing fails.
    image_dir = Path(args.output).resolve() / init.strftime('%Y%m%d')
    image_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='meteod_ecdsv4_') as temp:
        if args.grib_dir:
            raw = Path(args.grib_dir).resolve()
            paths = {v: raw / f'ecmwf_s2s_{v}_{init:%Y%m%d}.grib' for v in ('t2m', 'tp')}
            for path in paths.values():
                if not path.is_file():
                    raise FileNotFoundError(path)
        else:
            paths = download_grib(init, Path(temp))
        t = normalize_da(read_grib_da(paths['t2m']), init, 't2m')
        p = normalize_da(read_grib_da(paths['tp']), init, 'tp')
        if not np.array_equal(t['member'].values, p['member'].values):
            raise RuntimeError('Membri temperatura e precipitazione non allineati')
        tw = weekly_interpolated(t, starts, 't2m')
        pw = weekly_interpolated(p, starts, 'tp')
        del t, p
        city_data = city_probabilities(tw, pw, indices, thresholds, starts)
        render_histograms(city_data, thresholds, init, image_dir / 'citta')
        if not args.no_maps:
            render_maps(tw, pw, indices, thresholds, starts, init, image_dir / 'mappe')
    expected = len(indices) + (0 if args.no_maps else 2 * len(starts))
    if args.no_upload:
        pngs = list(image_dir.rglob('*.png'))
        if len(pngs) != expected:
            raise RuntimeError(f'PNG generati: {len(pngs)}; attesi: {expected}')
        print(f'DONE (senza upload): {len(pngs)} PNG in {image_dir}', flush=True)
    else:
        upload_pngs_to_r2(image_dir, expected)
        print(f'DONE: {expected} PNG caricati su R2', flush=True)


if __name__ == '__main__':
    main()
