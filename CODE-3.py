import numpy as np
import xarray as xr
from math import pi, sin, cos, tan, log

f_uz = r"E:\CMIP6\ACCESS\sfcWind_Amon_ACCESS-ESM1-5_1pctCO2_r1i1p1f1_gn_010101-025012.nc"
u2 = xr.open_dataset(f_uz).sfcWind.sel(time=slice('1850-01', '1989-12')) * (4.87 / log(67.8 * 10 - 5.42))

f_P = r"E:\CMIP6\ACCESS\ps_Amon_ACCESS-ESM1-5_1pctCO2_r1i1p1f1_gn_010101-025012.nc"
P = xr.open_dataset(f_P).ps.sel(time=slice('1850-01', '1989-12')) / 1000

f_rss = r"E:\CMIP6\ACCESS\rss_Emon_ACCESS-ESM1-5_1pctCO2_r1i1p1f1_gn_010101-025012.nc"
Rns = xr.open_dataset(f_rss).rss.sel(time=slice('1850-01', '1989-12')) * 0.0864
f_rls = r"E:\CMIP6\ACCESS\rls_Emon_ACCESS-ESM1-5_1pctCO2_r1i1p1f1_gn_010101-025012.nc"
Rnl = xr.open_dataset(f_rls).rls.sel(time=slice('1850-01', '1989-12')) * 0.0864

f_Tasmax = r"E:\CMIP6\ACCESS\tasmax_Amon_ACCESS-ESM1-5_1pctCO2_r1i1p1f1_gn_010101-025012.nc"
Tmax = xr.open_dataset(f_Tasmax).tasmax.sel(time=slice('1850-01', '1989-12')) - 273.16
f_Tasmin = r"E:\CMIP6\ACCESS\tasmin_Amon_ACCESS-ESM1-5_1pctCO2_r1i1p1f1_gn_010101-025012.nc"
Tmin = xr.open_dataset(f_Tasmin).tasmin.sel(time=slice('1850-01', '1989-12')) - 273.16

f_RH = r"E:\CMIP6\ACCESS\hurs_Amon_ACCESS-ESM1-5_1pctCO2_r1i1p1f1_gn_010101-025012.nc"
RH = xr.open_dataset(f_RH).hurs.sel(time=slice('1850-01', '1989-12'))

f_CO2 = r"E:\CMIP6\ACCESS\co2s_1%_0101-0240.nc"
CO2 = xr.open_dataset(f_CO2).co2s

# 计算净辐射Rn
rn = Rns + Rnl
Rn = xr.where(rn < 0, 0, rn)

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
PET = np.zeros((1680, 145, 192), dtype=float, order='C') * np.nan
for i in range(1680):
    PET[i, :, :] = PM_ET0(Tmax[i, :, :], Tmin[i, :, :], RH[i, :, :], P[i, :, :], u2[i, :, :], Rn[i, :, :], CO2[i, :, :])

# 计算每个月的天数
month_lengths = RH.time.dt.days_in_month

# 将month_lengths转换为和PET相同的形状
month_lengths_expanded = np.expand_dims(np.expand_dims(month_lengths, axis=1), axis=2)
PET_total = PET * month_lengths_expanded

# 创建数据集
data = xr.Dataset(
    {
        "PET": (["time", "lat", "lon"], PET_total),
    },
    coords={
        "time": RH.coords["time"],
        "lat": RH.coords["lat"],
        "lon": RH.coords["lon"],
    },
)

# 保存为nc文件
output_file = "E:\\PET_results\\ACCESS\\CTL_PET.nc"
data.to_netcdf(output_file)

print(PET_total)
