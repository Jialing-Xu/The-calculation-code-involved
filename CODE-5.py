import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

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


# 导入数据
f_pet = r"E:\PET_results\UKESM\CTL_PET.nc"
pet = xr.open_dataset(f_pet).PET.sel(time=slice('1851-01', '1880-12'))
hou30 = xr.open_dataset(f_pet).PET.sel(time=slice('1955-01', '1984-12'))
f_tas1 = "E:\\PET_results\\UKESM\\1PET_tas_1.nc"
pet_tas1 = xr.open_dataset(f_tas1).PET
f_tas2 = "E:\\PET_results\\UKESM\\1PET_tas_2.nc"
pet_tas2 = xr.open_dataset(f_tas2).PET
f_tas3 = "E:\\PET_results\\UKESM\\1PET_tas_3.nc"
pet_tas3 = xr.open_dataset(f_tas3).PET
f_rn1 = "E:\\PET_results\\UKESM\\1PET_rn_1.nc"
pet_rn1 = xr.open_dataset(f_rn1).PET
f_rn2 = "E:\\PET_results\\UKESM\\1PET_rn_2.nc"
pet_rn2 = xr.open_dataset(f_rn2).PET
f_rn3 = "E:\\PET_results\\UKESM\\1PET_rn_3.nc"
pet_rn3 = xr.open_dataset(f_rn3).PET
f_rh1 = "E:\\PET_results\\UKESM\\1PET_rh_1.nc"
pet_rh1 = xr.open_dataset(f_rh1).PET
f_rh2 = "E:\\PET_results\\UKESM\\1PET_rh_2.nc"
pet_rh2 = xr.open_dataset(f_rh2).PET
f_rh3 = "E:\\PET_results\\UKESM\\1PET_rh_3.nc"
pet_rh3 = xr.open_dataset(f_rh3).PET
f_u21 = "E:\\PET_results\\UKESM\\1PET_u2_1.nc"
pet_u21 = xr.open_dataset(f_u21).PET
f_u22 = "E:\\PET_results\\UKESM\\1PET_u2_2.nc"
pet_u22 = xr.open_dataset(f_u22).PET
f_u23 = "E:\\PET_results\\UKESM\\1PET_u2_3.nc"
pet_u23 = xr.open_dataset(f_u23).PET
f_co2 = r"E:\PET_results\UKESM\PET_co2.nc"
pet_co2 = xr.open_dataset(f_co2).PET

# 计算时间平均
PET = pet.groupby('time.year').sum(dim='time').mean(dim='year')
Hou30 = hou30.groupby('time.year').sum(dim='time').mean(dim='year')
Stas1 = pet_tas1.groupby('time.year').sum(dim='time').mean(dim='year')
Stas2 = pet_tas2.groupby('time.year').sum(dim='time').mean(dim='year')
Stas3 = pet_tas3.groupby('time.year').sum(dim='time').mean(dim='year')
Srn1 = pet_rn1.groupby('time.year').sum(dim='time').mean(dim='year')
Srn2 = pet_rn2.groupby('time.year').sum(dim='time').mean(dim='year')
Srn3 = pet_rn3.groupby('time.year').sum(dim='time').mean(dim='year')
Srh1 = pet_rh1.groupby('time.year').sum(dim='time').mean(dim='year')
Srh2 = pet_rh2.groupby('time.year').sum(dim='time').mean(dim='year')
Srh3 = pet_rh3.groupby('time.year').sum(dim='time').mean(dim='year')
Su21 = pet_u21.groupby('time.year').sum(dim='time').mean(dim='year')
Su22 = pet_u22.groupby('time.year').sum(dim='time').mean(dim='year')
Su23 = pet_u23.groupby('time.year').sum(dim='time').mean(dim='year')
Sco2 = pet_co2.groupby('time.year').sum(dim='time').mean(dim='year')

# 将每个时间平均后的结果转换为一维数组
PET_flat = PET.values.flatten()
Hou30_flat = Hou30.values.flatten()
Stas1_flat = Stas1.values.flatten()
Stas2_flat = Stas2.values.flatten()
Stas3_flat = Stas3.values.flatten()
Srn1_flat = Srn1.values.flatten()
Srn2_flat = Srn2.values.flatten()
Srn3_flat = Srn3.values.flatten()
Srh1_flat = Srh1.values.flatten()
Srh2_flat = Srh2.values.flatten()
Srh3_flat = Srh3.values.flatten()
Su21_flat = Su21.values.flatten()
Su22_flat = Su22.values.flatten()
Su23_flat = Su23.values.flatten()
Sco2_flat = Sco2.values.flatten()

# 定义X1-11
X1 = Stas1_flat - PET_flat
X2 = Stas2_flat - PET_flat
X3 = Stas3_flat - PET_flat
X4 = Srn1_flat - PET_flat
X5 = Srn2_flat - PET_flat
X6 = Srn3_flat - PET_flat
X7 = Srh1_flat - PET_flat
X8 = Srh2_flat - PET_flat
X9 = Srh3_flat - PET_flat
X10 = Su21_flat - PET_flat
X11 = Su22_flat - PET_flat
X12 = Su23_flat - PET_flat
X13 = Sco2_flat - PET_flat

# 构造方程组的系数矩阵
A = np.array([[0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
              [1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
              [1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
              [1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1],
              [1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1],
              [1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1],
              [1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1],
              [1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1],
              [1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1],
              [1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1],
              [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1],
              [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1],
              [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0]])

# 构造方程右侧的常数向量
X = np.array([X1, X2, X3, X4, X5, X6, X7, X8, X9, X10, X11, X12, X13])

# 求解方程组
result = np.linalg.solve(A, X)
C_T1 = result[0]
C_T2 = result[1]
C_T3 = result[2]
C_RN1 = result[3]
C_RN2 = result[4]
C_RN3 = result[5]
C_RH1 = result[6]
C_RH2 = result[7]
C_RH3 = result[8]
C_U1 = result[9]
C_U2 = result[10]
C_U3 = result[11]
C_CO2 = result[12]
C_all = Hou30_flat - PET_flat


# 创建结果数据集并保存
def save_as_netcdf(data, var_name, filename):
    lon = PET.lon
    lat = PET.lat
    
    ds = xr.Dataset(
        {var_name: (('lat', 'lon'), data.reshape(len(lat), len(lon)))},
        coords={'lat': lat, 'lon': lon}
    )
    
    data_ad = adjust_longitude(ds)
    rdata = bilinear_interpolation(data_ad)
    rdata.to_netcdf(filename)

# 保存各个分量结果
save_as_netcdf(C_T1, 'T1', 'E:\\PET_results\\UKESM\\1C_T1.nc')
save_as_netcdf(C_T2, 'T2', 'E:\\PET_results\\UKESM\\1C_T2.nc')
save_as_netcdf(C_T3, 'T3', 'E:\\PET_results\\UKESM\\1C_T3.nc')
save_as_netcdf(C_RN1, 'RN1', 'E:\\PET_results\\UKESM\\1C_RN1.nc')
save_as_netcdf(C_RN2, 'RN2', 'E:\\PET_results\\UKESM\\1C_RN2.nc')
save_as_netcdf(C_RN3, 'RN3', 'E:\\PET_results\\UKESM\\1C_RN3.nc')
save_as_netcdf(C_RH1, 'RH1', 'E:\\PET_results\\UKESM\\1C_RH1.nc')
save_as_netcdf(C_RH2, 'RH2', 'E:\\PET_results\\UKESM\\1C_RH2.nc')
save_as_netcdf(C_RH3, 'RH3', 'E:\\PET_results\\UKESM\\1C_RH3.nc')
save_as_netcdf(C_U1, 'U1', 'E:\\PET_results\\UKESM\\1C_U1.nc')
save_as_netcdf(C_U2, 'U2', 'E:\\PET_results\\UKESM\\1C_U2.nc')
save_as_netcdf(C_U3, 'U3', 'E:\\PET_results\\UKESM\\1C_U3.nc')
save_as_netcdf(C_CO2, 'CO2', 'E:\\PET_results\\UKESM\\1C_CO2.nc')
save_as_netcdf(C_all, 'CTL', 'E:\\PET_results\\UKESM\\1C_all.nc')


