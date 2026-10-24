"""Deterministic offline tests: no CDS credentials, no web access, no map rendering."""
from datetime import datetime, timezone
import tempfile
from pathlib import Path
import numpy as np
import xarray as xr
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import pipeline


def fake_arrays(init):
    rng = np.random.default_rng(42)
    nm = 12
    hr=np.arange(24,1105,24)
    lat=np.arange(55.,27.,-1.5)
    lon=np.arange(-15.,45.,1.5)
    # data shape (member, step, lat,lon)
    shape=(nm,len(hr),len(lat),len(lon))
    temp=np.empty(shape,dtype=np.float32)
    prec=np.empty(shape,dtype=np.float32)
    latval=lat[None,None,:,None]
    for i,lead in enumerate(hr):
        temp[:,i,:,:]=273.15+10+0.25*(latval-40)+rng.normal(0,2,(nm,1,len(lon)))
        prec[:,i,:,:]=lead/24*(2+0.1*np.arange(nm)[:,None,None])+0.01*np.arange(len(lon))[None,None,:]
    coords={'member':np.arange(1,nm+1), 'step':hr.astype('timedelta64[h]'),'latitude':lat,'longitude':lon}
    t=xr.DataArray(temp,dims=('member','step','latitude','longitude'),coords=coords,attrs={'units':'K'})
    p=xr.DataArray(prec,dims=('member','step','latitude','longitude'),coords=coords,attrs={'units':'kg m**-2'})
    return t,p


def test_all():
    init=datetime(2026,10,22,tzinfo=timezone.utc)
    starts=pipeline.week_starts(init)
    indices,thresholds=pipeline.load_assets()
    pipeline.check_threshold_weeks(indices,thresholds,starts)
    t,p=fake_arrays(init)
    tn=pipeline.normalize_da(t,init,'t2m')
    pn=pipeline.normalize_da(p,init,'tp')
    # Correct precip totals, not accumulations of cumulative fields.
    assert np.isclose(float(pn.isel(member=0,day=5,latitude=0,longitude=0)),2.0,atol=.01)
    tw=pipeline.weekly_interpolated(tn,starts,'t2m')
    pw=pipeline.weekly_interpolated(pn,starts,'tp')
    assert tw.shape==(5,12,58,136),(tw.shape,pw.shape)
    assert np.isclose(float(pw.isel(week=0,member=0,latitude=0,longitude=0)),14.0,atol=.1)
    city=pipeline.city_probabilities(tw,pw,indices,thresholds,starts)
    assert len(city)==151
    for name in ('Bologna','Milano','Londra'):
        assert len(city[name])==4
        for rec in city[name].values():
            assert abs(sum(rec[f'T_{i}'] for i in range(1,6))-100)<0.001
            assert abs(sum(rec[f'P_{i}'] for i in range(1,6))-100)<0.001
    n,week,arrays,lat2d,lon2d=next(pipeline.build_map_probabilities(tw,pw,thresholds,indices,starts))
    assert arrays['2tp'].shape==(13,58,136)
    assert np.allclose(arrays['2tp'][:5].sum(axis=0),100,atol=.01)
    with tempfile.TemporaryDirectory() as td:
        out=Path(td)
        pipeline.write_probability_netcdf(tw,pw,indices,thresholds,starts,out/'weekly.nc')
        assert (out/'weekly.nc').stat().st_size>100
        import plotting
        src=city['Bologna']
        plotting.plot_city_combined('Bologna',sorted(src),src,thresholds['Bologna'],str(out/'bologna.webp'),'22/10/2026')
        assert (out/'bologna.webp').stat().st_size>100
    print('TEST OFFLINE PASS: griglia, precipitazione, 151 città, quintili, map array, NetCDF, istogramma')

if __name__=='__main__':test_all()
