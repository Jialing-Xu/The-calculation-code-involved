import numpy as np
import xarray as xr
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


# #### 温度变化
# f_TmaxCTL = r"E:\CMIP6\ACCESS\tasmax_Amon_ACCESS-ESM1-5_1pctCO2-bgc_r1i1p1f1_gn_010101-025012.nc"
# TmaxCTLq30 = xr.open_dataset(f_TmaxCTL).tasmax.sel(time=slice('0102-01', '0131-12')) - 273.16
# TmaxCTLh30 = xr.open_dataset(f_TmaxCTL).tasmax.sel(time=slice('0206-01', '0235-12')) - 273.16

# # f_TminCTL = r"E:\CMIP6\UKESM\tasmin_1pctCO2-bgc_1850-1999.nc"
# # TminCTLq30 = xr.open_dataset(f_TminCTL).tasmin.sel(time=slice('1851-01', '1880-12')) - 273.16
# # TminCTLh30 = xr.open_dataset(f_TminCTL).tasmin.sel(time=slice('1955-01', '1984-12')) - 273.16

# # hou30 = 0.5 * (TmaxCTLh30 + TminCTLh30)
# # spei_30 = 0.5 * (TmaxCTLq30 + TminCTLq30)

# hou30 = TmaxCTLh30
# spei_30 = TmaxCTLq30


# ##### 降水变化
# f_CTL = r"E:\CMIP6\ACCESS\pr_Amon_ACCESS-ESM1-5_1pctCO2-bgc_r1i1p1f1_gn_010101-025012.nc"
# prq30 = xr.open_dataset(f_CTL).pr.sel(time=slice('0102-01', '0131-12')) *60*60*24
# prh30 = xr.open_dataset(f_CTL).pr.sel(time=slice('0206-01', '0235-12')) *60*60*24

# month_q30 = prq30.time.dt.days_in_month
# month_h30 = prh30.time.dt.days_in_month

# spei_30 = prq30 * month_q30
# hou30 = prh30 * month_h30

# spei_30 = spei_30.where(spei_30 >= 0, 0)
# hou30 = hou30.where(hou30 >= 0, 0)


# ##### 相对湿度变化
# f_CTL = r"E:\CMIP6\MPI\hurs_1pctCO2-bgc_1850-1989.nc"
# spei_30 = xr.open_dataset(f_CTL).hurs.sel(time=slice('1851-01', '1880-12'))
# hou30 = xr.open_dataset(f_CTL).hurs.sel(time=slice('1955-01', '1984-12'))


##### 净辐射变化
f_rlsCTL = r"E:\CMIP6\UKESM\rls_1pctCO2-bgc_1850-1999.nc"
rlsq30 = xr.open_dataset(f_rlsCTL).rls.sel(time=slice('1851-01', '1880-12'))
rlsh30 = xr.open_dataset(f_rlsCTL).rls.sel(time=slice('1955-01', '1984-12'))  

f_rssCTL = r"E:\CMIP6\UKESM\rss_1pctCO2-bgc_1850-1999.nc"
rssq30 = xr.open_dataset(f_rssCTL).rss.sel(time=slice('1851-01', '1880-12'))
rssh30 = xr.open_dataset(f_rssCTL).rss.sel(time=slice('1955-01', '1984-12'))

month_q30 = rlsq30.time.dt.days_in_month
month_h30 = rlsh30.time.dt.days_in_month

rnq30 = (rssq30 + rlsq30) * 0.0864 * month_q30
rnh30 = (rssh30 + rlsh30) * 0.0864 * month_h30

spei_30 = xr.where(rnq30 < 0, 0, rnq30) 
hou30 = xr.where(rnh30 < 0, 0, rnh30) 


# ##### 风速变化
# f_CTL = r"E:\CMIP6\MPI\sfcWind_1pctCO2-rad_1850-1989.nc"
# spei_30 = xr.open_dataset(f_CTL).sfcWind.sel(time=slice('1851-01', '1880-12')) * (4.87 / log(67.8 * 10 - 5.42))
# hou30 = xr.open_dataset(f_CTL).sfcWind.sel(time=slice('1955-01', '1984-12')) * (4.87 / log(67.8 * 10 - 5.42))


# ##### 二氧化碳变化
# f_CTL = r"E:\CMIP6\ACCESS\co2s_1%.nc"
# spei_30 = xr.open_dataset(f_CTL).co2s.sel(time=slice('1851-01', '1880-12'))
# hou30 = xr.open_dataset(f_CTL).co2s.sel(time=slice('1955-01', '1984-12'))




# # 计算 spei_30 和 hou30 的时间平均
# SPEI3 = spei_30.mean(dim='time')
# Hou30 = hou30.mean(dim='time')
# co2 = Hou30 - SPEI3

# data = adjust_longitude(co2)
# rdata = bilinear_interpolation(data)
# rdata.to_netcdf(r"E:\CMIP6\ACCESS\CO2_1pct.nc")



# 按年分组并计算年总和
SPEI3 = spei_30.groupby('time.year').sum(dim='time').mean(dim='year')
Hou30 = hou30.groupby('time.year').sum(dim='time').mean(dim='year')
co2 = Hou30 - SPEI3

data = adjust_longitude(co2)
rdata = bilinear_interpolation(data)
rdata.to_netcdf(r"E:\CMIP6\UKESM\RN_bgc.nc")


