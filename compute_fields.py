import numpy as np
import xarray as xr

# Define Earth's radius (meters)
R = 6.371e6  

# Load CMIP6 dataset (Modify file path accordingly)
file_path = "path_to_cmip6_data.nc"
ds = xr.open_dataset(file_path)

# Extract required variables (ensure correct variable names in CMIP6 dataset)
ua_850 = ds['ua'].sel(plev=85000)  # 850 hPa eastward wind (m/s)
va_850 = ds['va'].sel(plev=85000)  # 850 hPa northward wind (m/s)
ua_200 = ds['ua'].sel(plev=20000)  # 200 hPa eastward wind (m/s)
va_200 = ds['va'].sel(plev=20000)  # 200 hPa northward wind (m/s)
humidity_600 = ds['hus'].sel(plev=60000)  # 600 hPa specific humidity (kg/kg)
mpi = ds['mpi']  # Maximum Potential Intensity (if available)

# Convert lat/lon to radians
lat_radians = np.deg2rad(ds.lat)
lon_radians = np.deg2rad(ds.lon)

# Compute spacing in radians
dlat = np.gradient(lat_radians)
dlon = np.gradient(lon_radians)

# Compute spacing in meters
dx = R * np.cos(lat_radians[:, np.newaxis]) * dlon[np.newaxis, :]
dy = R * dlat[:, np.newaxis]

# Compute Relative Vorticity (ζ) at 850 hPa
vorticity_850 = (np.gradient(va_850, dx, axis=-1) - np.gradient(ua_850, dy, axis=-2))

# Compute Vertical Wind Shear (ΔU, ΔV)
shear_u = ua_200 - ua_850  # Shear in eastward wind
shear_v = va_200 - va_850  # Shear in northward wind
vertical_wind_shear = np.sqrt(shear_u**2 + shear_v**2)  # Magnitude of wind shear (m/s)

# Save output dataset
output_ds = xr.Dataset({
    'relative_vorticity_850': (["time", "lat", "lon"], vorticity_850),
    'vertical_wind_shear': (["time", "lat", "lon"], vertical_wind_shear),
    'humidity_600': (["time", "lat", "lon"], humidity_600),
    'maximum_potential_intensity': (["time", "lat", "lon"], mpi)
}, coords={"time": ds.time, "lat": ds.lat, "lon": ds.lon})

# Save to NetCDF
output_file = "processed_cmip6.nc"
output_ds.to_netcdf(output_file)

print(f"Processed data saved to {output_file}")
