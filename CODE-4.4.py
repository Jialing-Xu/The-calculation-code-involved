import numpy as np
import xarray as xr
from math import pi, sin, cos, tan, log

f_rn = r"E:\mingan_results\UKESM\1RN_no_int.nc"
Rn = xr.open_dataset(f_rn).__xarray_dataarray_variable__

f_uz = r"E:\CMIP6\UKESM\sfcWind_1pctCO2_1850-1999.nc"
u2 = xr.open_dataset(f_uz).sfcWind.sel(time=slice('1955-01', '1984-12')) * (4.87 / log(67.8 * 10 - 5.42))

f_P = r"E:\CMIP6\UKESM\ps_1pctCO2_1850-1999.nc"
P = xr.open_dataset(f_P).ps.sel(time=slice('1955-01', '1984-12')) / 1000

f_Tasmax = r"E:\CMIP6\UKESM\tasmax_1pctCO2_1850-1999.nc"
Tmax = xr.open_dataset(f_Tasmax).tasmax.sel(time=slice('1955-01', '1984-12')) - 273.16
f_Tasmin = r"E:\CMIP6\UKESM\tasmin_1pctCO2_1850-1999.nc"
Tmin = xr.open_dataset(f_Tasmin).tasmin.sel(time=slice('1955-01', '1984-12')) - 273.16

f_RH = r"E:\CMIP6\UKESM\hurs_1pctCO2_1850-1999.nc"
RH = xr.open_dataset(f_RH).hurs.sel(time=slice('1955-01', '1984-12'))

f_CO2 = r"E:\CMIP6\UKESM\co2s_1%.nc"
CO2 = xr.open_dataset(f_CO2).co2s.sel(time=slice('1955-01', '1984-12'))


def PM_ET0(Tmax, Tmin, RH, P, u2, Rn, CO2):
    Tmean = 0.5 * (Tmax + Tmin)
    G = 0
    ETMAX = 0.6108 * (np.exp((17.27 * Tmax) / (Tmax + 237.3)))
    ETMIN = 0.6108 * (np.exp((17.27 * Tmin) / (Tmin + 237.3)))
    es = 0.5 * (ETMAX + ETMIN)
    ea = RH * es / 100
    ETmean = 0.6108 * (np.exp((17.27 * Tmean) / (Tmean + 237.3)))
    delta = 4098 * ETmean / ((Tmean + 237.3) * (Tmean + 237.3))
    r = 0.665 * 0.001 * P
    # ET0 = ( 0.408*delta*(Rn-G) + r*900*u2*(es-ea)/(Tmean+273) )/( delta + r*(1 + (0.34+0.00024*(CO2-300))*u2) )
    ET0_1 = 0.408 * delta * (Rn - G)
    ET0_2 = r * 900 * u2 * (es - ea) / (Tmean + 273)
    ET0_3 = delta + r * (1 + (0.34 + 0.00024 * (CO2 - 300)) * u2)
    ET0 = (ET0_1 + ET0_2) / ET0_3
    return ET0


# 计算月平均PET:[mm/day]
PET = np.zeros((360, 144, 192), dtype=float, order='C') * np.nan
for i in range(360):
    PET[i, :, :] = PM_ET0(Tmax[i, :, :], Tmin[i, :, :], RH[i, :, :], P[i, :, :], u2[i, :, :], Rn[i, :, :], CO2[i, :, :])

# 计算每个月的天数
month_lengths = P.time.dt.days_in_month

# 将month_lengths转换为和PET相同的形状
month_lengths_expanded = np.expand_dims(np.expand_dims(month_lengths, axis=1), axis=2)
PET_total = PET * month_lengths_expanded

# 创建数据集
data = xr.Dataset(
    {
        "PET": (["time", "lat", "lon"], PET_total),
    },
    coords={
        "time": P.coords["time"],
        "lat": P.coords["lat"],
        "lon": P.coords["lon"],
    },
)

# 保存为nc文件
output_file = "E:\\PET_results\\UKESM\\1PET_rn_3.nc"
data.to_netcdf(output_file)

print(PET_total)
