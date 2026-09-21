#!/usr/bin/env python3
"""Extract one TCO319 LUH3 state year, crop partition, and corrected peat map."""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import numpy as np
from netCDF4 import Dataset, num2date

NATURAL=("primf","primn","secdf","secdn")
CROP=("c3ann","c3per","c3nfx","c4ann","c4per")
PASTURE=("pastr","range")
URBAN=("urban",)
REQUIRED=NATURAL+CROP+PASTURE+URBAN
ST_NAMES=("CC3ann","CC3per","CC3nfx","CC4ann","CC4per","CC3anni","CC3peri","CC3nfxi","CC4anni","CC4peri")

def sha256(path):
 d=hashlib.sha256()
 with path.open('rb') as f:
  for block in iter(lambda:f.read(8*1024*1024),b''): d.update(block)
 return d.hexdigest()

def time_index(dataset,year):
 t=dataset['time']; dates=num2date(t[:],units=t.units,calendar=getattr(t,'calendar','standard'),only_use_cftime_datetimes=True)
 found=[i for i,d in enumerate(dates) if d.year==year]
 if len(found)!=1: raise ValueError(f'year {year} occurs {len(found)} times')
 return found[0]

def read_peat(path,lon,lat):
 rows=np.loadtxt(path,skiprows=1,dtype=np.float64)
 if rows.shape!=(lon.size,3): raise ValueError(f'peat map shape {rows.shape}; expected {(lon.size,3)}')
 dl=np.abs(((rows[:,0]-lon+180)%360)-180); da=np.abs(rows[:,1]-lat)
 if dl.max()>1e-3 or da.max()>1e-3: raise ValueError('peat-map coordinates/order do not match LUH3')
 peat=rows[:,2]
 if np.any(~np.isfinite(peat)) or np.any((peat<0)|(peat>1)): raise ValueError('invalid peat')
 return peat

def main():
 p=argparse.ArgumentParser(); p.add_argument('--luh3',required=True,type=Path); p.add_argument('--management',required=True,type=Path); p.add_argument('--peat',required=True,type=Path); p.add_argument('--year',required=True,type=int); p.add_argument('--output',required=True,type=Path); p.add_argument('--summary',required=True,type=Path); a=p.parse_args()
 with Dataset(a.luh3) as ds, Dataset(a.management) as mg:
  si=time_index(ds,a.year); mi=time_index(mg,a.year)
  lon=np.asarray(ds['lon'][:],dtype=np.float64); lat=np.asarray(ds['lat'][:],dtype=np.float64)
  mlon=np.asarray(mg['lon'][:],dtype=np.float64); mlat=np.asarray(mg['lat'][:],dtype=np.float64)
  if not (np.array_equal(lon,mlon) and np.array_equal(lat,mlat)): raise ValueError('state/management grids differ')
  peat=read_peat(a.peat,lon,lat); values={}; valid_count=np.zeros(lon.size,dtype=np.int16)
  for name in REQUIRED:
   field=ds[name][si,:]; mask=np.ma.getmaskarray(field); data=np.asarray(np.ma.filled(field,np.nan),dtype=np.float64); ok=(~mask)&np.isfinite(data)&(data>=0)&(data<=1); valid_count+=ok; values[name]=np.where(ok,data,0.0)
  crop_st=[]
  for name in CROP:
   field=mg['irrig_'+name][mi,:]; irr=np.asarray(np.ma.filled(field,0.0),dtype=np.float64); irr=np.where(np.isfinite(irr),np.clip(irr,0,1),0.0); area=values[name]; crop_st.append(area*(1-irr))
  for name in CROP:
   field=mg['irrig_'+name][mi,:]; irr=np.asarray(np.ma.filled(field,0.0),dtype=np.float64); irr=np.where(np.isfinite(irr),np.clip(irr,0,1),0.0); crop_st.append(values[name]*irr)
 partial=(valid_count!=0)&(valid_count!=len(REQUIRED))
 if partial.any(): raise ValueError(f'partial LUH3 state at row {int(np.flatnonzero(partial)[0])}')
 valid=valid_count==len(REQUIRED); missing=~valid
 natural=sum(values[n] for n in NATURAL); crop=sum(values[n] for n in CROP); pasture=sum(values[n] for n in PASTURE); urban=values['urban']
 natural[missing]=1; crop[missing]=0; pasture[missing]=0; urban[missing]=0
 for v in crop_st: v[missing]=0
 if np.any(peat[missing]!=0): raise ValueError('missing LUH3 row has nonzero peat')
 matrix=np.column_stack((lon,lat,valid.astype(np.int8),natural,crop,pasture,urban,peat,*crop_st))
 header='Lon Lat LUH3Valid RawNatural RawCropland RawPasture RawUrban CorrectedPeat '+' '.join('Crop_'+n for n in ST_NAMES)
 a.output.parent.mkdir(parents=True,exist_ok=True); fmt=("%.17g","%.17g","%d")+("%.17g",)*15
 np.savetxt(a.output,matrix,fmt=fmt,delimiter='\t',header=header,comments='')
 raw=natural+crop+pasture+urban
 summary={'year':a.year,'rows':int(lon.size),'valid_luh3_rows':int(valid.sum()),'all_fields_missing_rows_using_natural_fallback':int(missing.sum()),'partial_rows':int(partial.sum()),'corrected_peat_nonzero_rows':int((peat>0).sum()),'maximum_corrected_peat':float(peat.max()),'raw_luh3_total_min_valid':float(raw[valid].min()),'raw_luh3_total_max_valid':float(raw[valid].max()),'luh3_source':str(a.luh3.resolve()),'luh3_sha256':sha256(a.luh3),'management_source':str(a.management.resolve()),'management_sha256':sha256(a.management),'peat_source':str(a.peat.resolve()),'peat_sha256':sha256(a.peat),'target_table':str(a.output.resolve()),'target_sha256':sha256(a.output),'crop_columns':list(ST_NAMES)}
 a.summary.write_text(json.dumps(summary,indent=2)+'\n'); print(json.dumps(summary,indent=2))
if __name__=='__main__': main()
