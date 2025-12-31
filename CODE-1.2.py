import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import xarray as xr
import matplotlib as mpl
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from palettable.colorbrewer.diverging import RdBu_11_r
import regionmask
import geopandas as gp
from math import pi, sin, cos, tan, log

# 定义目标网格的经度和纬度范围
lon_range = np.arange(-180, 180, 2)
lat_range = np.arange(-60, 92, 2)

# 插值
def bilinear_interpolation(data):
    resampled_data = data.interp(lon=lon_range, lat=lat_range, method='linear')
    return resampled_data


# def adjust_longitude(dataset):
#     lon_name = 'lon'
#     dataset = dataset.assign_coords({lon_name: (((dataset[lon_name] + 180) % 360) - 180)})
#     dataset = dataset.sortby(lon_name)
    
#     if 180 not in dataset[lon_name].values:
#         minus_180_data = dataset.sel({lon_name: -180}, method="nearest")
#         new_lon_values = np.append(dataset[lon_name].values, 180)
#         new_lon_values_sorted = np.sort(new_lon_values)
#         dataset = dataset.reindex(
#             {lon_name: new_lon_values_sorted},
#             method="nearest", 
#             tolerance=10
#         )
#         dataset.loc[{lon_name: 180}] = minus_180_data
    
#     return dataset


# 用于UKESM模式（-179.1, 179.1）
def adjust_longitude(dataset):
    lon_name = 'lon'
    dataset = dataset.assign_coords({lon_name: (((dataset[lon_name] + 180) % 360) - 180)})
    dataset = dataset.sortby(lon_name)
    lon_min = dataset[lon_name].min().item()
    lon_max = dataset[lon_name].max().item()
    
    if lon_min > -180 or lon_max < 180:
        left_data = dataset.sel({lon_name: lon_min}, method="nearest")
        right_data = dataset.sel({lon_name: lon_max}, method="nearest")
        new_lon_values = np.unique(
            np.concatenate([
                [-180], 
                dataset[lon_name].values,
                [180]  
            ])
        )
        dataset = dataset.reindex(
            {lon_name: new_lon_values},
            method="nearest",  
            tolerance=10 
        )
        dataset.loc[{lon_name: -180}] = right_data
        dataset.loc[{lon_name: 180}] = left_data
    
    return dataset

def save_as_netcdf(data, filename):
    data_ad = adjust_longitude(data)
    rdata = bilinear_interpolation(data_ad)
    rdata.to_netcdf(filename)



##### 降水
f_PCTL = r"E:\CMIP6\UKESM\pr_1pctCO2_1850-1999.nc"
PCTLq30 = xr.open_dataset(f_PCTL).pr.sel(time=slice('1851-01', '1880-12')) *60*60*24
PCTLh30 = xr.open_dataset(f_PCTL).pr.sel(time=slice('1955-01', '1984-12')) *60*60*24
f_Prad = r"E:\CMIP6\UKESM\pr_1pctCO2-rad_1850-1999.nc"
Pradq30 = xr.open_dataset(f_Prad).pr.sel(time=slice('1851-01', '1880-12')) *60*60*24
Pradh30 = xr.open_dataset(f_Prad).pr.sel(time=slice('1955-01', '1984-12')) *60*60*24
f_Pbgc = r"E:\CMIP6\UKESM\pr_1pctCO2-bgc_1850-1999.nc"
Pbgcq30 = xr.open_dataset(f_Pbgc).pr.sel(time=slice('1851-01', '1880-12')) *60*60*24
Pbgch30 = xr.open_dataset(f_Pbgc).pr.sel(time=slice('1955-01', '1984-12')) *60*60*24

month_q30 = PCTLq30.time.dt.days_in_month
month_h30 = PCTLh30.time.dt.days_in_month

PrCTLq30 = (PCTLq30 * month_q30).where(lambda x: x >= 0, 0)
PrCTLh30 = (PCTLh30 * month_h30).where(lambda x: x >= 0, 0)
Prradq30 = (Pradq30 * month_q30).where(lambda x: x >= 0, 0)
Prradh30 = (Pradh30 * month_h30).where(lambda x: x >= 0, 0)
Prbgcq30 = (Pbgcq30 * month_q30).where(lambda x: x >= 0, 0)
Prbgch30 = (Pbgch30 * month_h30).where(lambda x: x >= 0, 0)

months = PrCTLh30['time.month'].values

PCTL = PrCTLh30 - PrCTLq30.groupby('time.month').mean('time').sel(month=months).values
Prad = Prradh30 - Prradq30.groupby('time.month').mean('time').sel(month=months).values
Pbgc = Prbgch30 - Prbgcq30.groupby('time.month').mean('time').sel(month=months).values
Pint = PCTL - Prad - Pbgc

Prad = Prad.groupby('time.year').sum(dim='time').mean(dim='year')
Pbgc = Pbgc.groupby('time.year').sum(dim='time').mean(dim='year')
Pint = Pint.groupby('time.year').sum(dim='time').mean(dim='year')

save_as_netcdf(Prad, r"E:\CMIP6\UKESM\P_rad.nc")
save_as_netcdf(Pbgc, r"E:\CMIP6\UKESM\P_bgc.nc")
save_as_netcdf(Pint, r"E:\CMIP6\UKESM\P_int.nc")



# ##### 温度
# f_TmaxCTL = r"E:\CMIP6\UKESM\tasmax_1pctCO2_1850-1999.nc"
# TmaxCTLq30 = xr.open_dataset(f_TmaxCTL).tasmax.sel(time=slice('1851-01', '1880-12')) - 273.16
# TmaxCTLh30 = xr.open_dataset(f_TmaxCTL).tasmax.sel(time=slice('1955-01', '1984-12')) - 273.16
# f_Tmaxrad = r"E:\CMIP6\UKESM\tasmax_1pctCO2-rad_1850-1999.nc"
# Tmaxradq30 = xr.open_dataset(f_Tmaxrad).tasmax.sel(time=slice('1851-01', '1880-12')) - 273.16
# Tmaxradh30 = xr.open_dataset(f_Tmaxrad).tasmax.sel(time=slice('1955-01', '1984-12')) - 273.16
# f_Tmaxbgc = r"E:\CMIP6\UKESM\tasmax_1pctCO2-bgc_1850-1999.nc"
# Tmaxbgcq30 = xr.open_dataset(f_Tmaxbgc).tasmax.sel(time=slice('1851-01', '1880-12')) - 273.16
# Tmaxbgch30 = xr.open_dataset(f_Tmaxbgc).tasmax.sel(time=slice('1955-01', '1984-12')) - 273.16

# f_TminCTL = r"E:\CMIP6\UKESM\tasmin_1pctCO2_1850-1999.nc"
# TminCTLq30 = xr.open_dataset(f_TminCTL).tasmin.sel(time=slice('1851-01', '1880-12')) - 273.16
# TminCTLh30 = xr.open_dataset(f_TminCTL).tasmin.sel(time=slice('1955-01', '1984-12')) - 273.16
# f_Tminrad = r"E:\CMIP6\UKESM\tasmin_1pctCO2-rad_1850-1999.nc"
# Tminradq30 = xr.open_dataset(f_Tminrad).tasmin.sel(time=slice('1851-01', '1880-12')) - 273.16
# Tminradh30 = xr.open_dataset(f_Tminrad).tasmin.sel(time=slice('1955-01', '1984-12')) - 273.16
# f_Tminbgc = r"E:\CMIP6\UKESM\tasmin_1pctCO2-bgc_1850-1999.nc"
# Tminbgcq30 = xr.open_dataset(f_Tminbgc).tasmin.sel(time=slice('1851-01', '1880-12')) - 273.16
# Tminbgch30 = xr.open_dataset(f_Tminbgc).tasmin.sel(time=slice('1955-01', '1984-12')) - 273.16

# TCTLq30 = 0.5 * (TmaxCTLq30 + TminCTLq30)
# TCTLh30 = 0.5 * (TmaxCTLh30 + TminCTLh30)
# Tradq30 = 0.5 * (Tmaxradq30 + Tminradq30)
# Tradh30 = 0.5 * (Tmaxradh30 + Tminradh30)
# Tbgcq30 = 0.5 * (Tmaxbgcq30 + Tminbgcq30)
# Tbgch30 = 0.5 * (Tmaxbgch30 + Tminbgch30)

# months = TCTLh30['time.month'].values

# TCTL = TCTLh30 - TCTLq30.groupby('time.month').mean('time').sel(month=months).values
# Trad = Tradh30 - Tradq30.groupby('time.month').mean('time').sel(month=months).values
# Tbgc = Tbgch30 - Tbgcq30.groupby('time.month').mean('time').sel(month=months).values
# Tint = TCTL - Trad - Tbgc

# Trad = Trad.mean(dim='time')
# Tbgc = Tbgc.mean(dim='time')
# Tint = Tint.mean(dim='time')

# save_as_netcdf(Trad, r"E:\CMIP6\UKESM\T_rad.nc")
# save_as_netcdf(Tbgc, r"E:\CMIP6\UKESM\T_bgc.nc")
# save_as_netcdf(Tint, r"E:\CMIP6\UKESM\T_int.nc")



# ##### 净辐射
# f_rssCTL = r"E:\CMIP6\UKESM\rss_1pctCO2_1850-1999.nc"
# RnsCTLq30 = xr.open_dataset(f_rssCTL).rss.sel(time=slice('1851-01', '1880-12'))
# RnsCTLh30 = xr.open_dataset(f_rssCTL).rss.sel(time=slice('1955-01', '1984-12'))
# f_rlsCTL = r"E:\CMIP6\UKESM\rls_1pctCO2_1850-1999.nc"
# RnlCTLq30 = xr.open_dataset(f_rlsCTL).rls.sel(time=slice('1851-01', '1880-12'))
# RnlCTLh30 = xr.open_dataset(f_rlsCTL).rls.sel(time=slice('1955-01', '1984-12'))

# month_q30 = RnsCTLq30.time.dt.days_in_month
# month_h30 = RnsCTLh30.time.dt.days_in_month

# RNCTLq30 = ((RnsCTLq30 + RnlCTLq30) * 0.0864 * month_q30).where(lambda x: x >= 0, 0)
# RNCTLh30 = ((RnsCTLh30 + RnlCTLh30) * 0.0864 * month_h30).where(lambda x: x >= 0, 0)

# f_rssrad = r"E:\CMIP6\UKESM\rss_1pctCO2-rad_1850-1999.nc"
# Rnsradq30 = xr.open_dataset(f_rssrad).rss.sel(time=slice('1851-01', '1880-12'))
# Rnsradh30 = xr.open_dataset(f_rssrad).rss.sel(time=slice('1955-01', '1984-12'))
# f_rlsrad = r"E:\CMIP6\UKESM\rls_1pctCO2-rad_1850-1999.nc"
# Rnlradq30 = xr.open_dataset(f_rlsrad).rls.sel(time=slice('1851-01', '1880-12'))
# Rnlradh30 = xr.open_dataset(f_rlsrad).rls.sel(time=slice('1955-01', '1984-12'))

# RNradq30 = ((Rnsradq30 + Rnlradq30) * 0.0864 * month_q30).where(lambda x: x >= 0, 0)
# RNradh30 = ((Rnsradh30 + Rnlradh30) * 0.0864 * month_h30).where(lambda x: x >= 0, 0)

# f_rssbgc = r"E:\CMIP6\UKESM\rss_1pctCO2-bgc_1850-1999.nc"
# Rnsbgcq30 = xr.open_dataset(f_rssbgc).rss.sel(time=slice('1851-01', '1880-12'))
# Rnsbgch30 = xr.open_dataset(f_rssbgc).rss.sel(time=slice('1955-01', '1984-12'))
# f_rlsbgc = r"E:\CMIP6\UKESM\rls_1pctCO2-bgc_1850-1999.nc"
# Rnlbgcq30 = xr.open_dataset(f_rlsbgc).rls.sel(time=slice('1851-01', '1880-12'))
# Rnlbgch30 = xr.open_dataset(f_rlsbgc).rls.sel(time=slice('1955-01', '1984-12'))

# RNbgcq30 = ((Rnsbgcq30 + Rnlbgcq30) * 0.0864 * month_q30).where(lambda x: x >= 0, 0)
# RNbgch30 = ((Rnsbgch30 + Rnlbgch30) * 0.0864 * month_h30).where(lambda x: x >= 0, 0)

# months = RNCTLh30['time.month'].values

# RNCTL = RNCTLh30 - RNCTLq30.groupby('time.month').mean('time').sel(month=months).values
# RNrad = RNradh30 - RNradq30.groupby('time.month').mean('time').sel(month=months).values
# RNbgc = RNbgch30 - RNbgcq30.groupby('time.month').mean('time').sel(month=months).values
# RNint = RNCTL - RNrad - RNbgc

# RNrad = RNrad.groupby('time.year').sum(dim='time').mean(dim='year')
# RNbgc = RNbgc.groupby('time.year').sum(dim='time').mean(dim='year')
# RNint = RNint.groupby('time.year').sum(dim='time').mean(dim='year')

# save_as_netcdf(RNrad, r"E:\CMIP6\UKESM\RN_rad.nc")
# save_as_netcdf(RNbgc, r"E:\CMIP6\UKESM\RN_bgc.nc")
# save_as_netcdf(RNint, r"E:\CMIP6\UKESM\RN_int.nc")



# ##### 相对湿度
# f_RHCTL = r"E:\CMIP6\UKESM\hurs_1pctCO2_1850-1999.nc"
# RHCTLq30 = xr.open_dataset(f_RHCTL).hurs.sel(time=slice('1851-01', '1880-12'))
# RHCTLh30 = xr.open_dataset(f_RHCTL).hurs.sel(time=slice('1955-01', '1984-12'))
# f_RHrad = r"E:\CMIP6\UKESM\hurs_1pctCO2-rad_1850-1999.nc"
# RHradq30 = xr.open_dataset(f_RHrad).hurs.sel(time=slice('1851-01', '1880-12'))
# RHradh30 = xr.open_dataset(f_RHrad).hurs.sel(time=slice('1955-01', '1984-12'))
# f_RHbgc = r"E:\CMIP6\UKESM\hurs_1pctCO2-bgc_1850-1999.nc"
# RHbgcq30 = xr.open_dataset(f_RHbgc).hurs.sel(time=slice('1851-01', '1880-12'))
# RHbgch30 = xr.open_dataset(f_RHbgc).hurs.sel(time=slice('1955-01', '1984-12'))

# months = RHCTLh30['time.month'].values

# RHCTL = RHCTLh30 - RHCTLq30.groupby('time.month').mean('time').sel(month=months).values
# RHrad = RHradh30 - RHradq30.groupby('time.month').mean('time').sel(month=months).values
# RHbgc = RHbgch30 - RHbgcq30.groupby('time.month').mean('time').sel(month=months).values
# RHint = RHCTL - RHrad - RHbgc

# RHrad = RHrad.mean(dim='time')
# RHbgc = RHbgc.mean(dim='time')
# RHint = RHint.mean(dim='time')

# save_as_netcdf(RHrad, r"E:\CMIP6\UKESM\RH_rad.nc")
# save_as_netcdf(RHbgc, r"E:\CMIP6\UKESM\RH_bgc.nc")
# save_as_netcdf(RHint, r"E:\CMIP6\UKESM\RH_int.nc")



# ##### 风速
# f_UCTL = r"E:\CMIP6\UKESM\sfcWind_1pctCO2_1850-1999.nc"
# UCTLq30 = xr.open_dataset(f_UCTL).sfcWind.sel(time=slice('1851-01', '1880-12')) * (4.87 / log(67.8 * 10 - 5.42))
# UCTLh30 = xr.open_dataset(f_UCTL).sfcWind.sel(time=slice('1955-01', '1984-12')) * (4.87 / log(67.8 * 10 - 5.42))
# f_Urad = r"E:\CMIP6\UKESM\sfcWind_1pctCO2-rad_1850-1999.nc"
# Uradq30 = xr.open_dataset(f_Urad).sfcWind.sel(time=slice('1851-01', '1880-12')) * (4.87 / log(67.8 * 10 - 5.42))
# Uradh30 = xr.open_dataset(f_Urad).sfcWind.sel(time=slice('1955-01', '1984-12')) * (4.87 / log(67.8 * 10 - 5.42))
# f_Ubgc = r"E:\CMIP6\UKESM\sfcWind_1pctCO2-bgc_1850-1999.nc"
# Ubgcq30 = xr.open_dataset(f_Ubgc).sfcWind.sel(time=slice('1851-01', '1880-12')) * (4.87 / log(67.8 * 10 - 5.42))
# Ubgch30 = xr.open_dataset(f_Ubgc).sfcWind.sel(time=slice('1955-01', '1984-12')) * (4.87 / log(67.8 * 10 - 5.42))

# months = UCTLh30['time.month'].values

# UCTL = UCTLh30 - UCTLq30.groupby('time.month').mean('time').sel(month=months).values
# Urad = Uradh30 - Uradq30.groupby('time.month').mean('time').sel(month=months).values
# Ubgc = Ubgch30 - Ubgcq30.groupby('time.month').mean('time').sel(month=months).values
# Uint = UCTL - Urad - Ubgc

# Urad = Urad.mean(dim='time')
# Ubgc = Ubgc.mean(dim='time')
# Uint = Uint.mean(dim='time')

# save_as_netcdf(Urad, r"E:\CMIP6\UKESM\U_rad.nc")
# save_as_netcdf(Ubgc, r"E:\CMIP6\UKESM\U_bgc.nc")
# save_as_netcdf(Uint, r"E:\CMIP6\UKESM\U_int.nc")


