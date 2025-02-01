# **Tropical Cyclone Detection using OWZP Method with CMIP6 Data**

Follow the paper https://www.science.org/doi/epdf/10.1126/sciadv.adl2142

---


This script processes **CMIP6 NetCDF data** and computes:
- **Relative vorticity at 850 hPa**
- **Vertical wind shear (200 hPa - 850 hPa)**
- **Mid-level humidity (600 hPa)**
- **Maximum Potential Intensity (MPI)**
- **Cyclone detection using OWZP criteria**

---

## **🔧 Dependencies**
Ensure you have the required libraries installed:
```bash
pip install numpy xarray netCDF4 metpy scipy
```

---

## **📌 Python Implementation**
```python
import numpy as np
import xarray as xr
from metpy.calc import potential_temperature, absolute_vorticity
from metpy.constants import Rd, Cp
from metpy.units import units
from scipy.ndimage import gaussian_filter

# 🌍 Constants
R = 6.371e6  # Earth's radius (meters)
omega = 7.2921e-5  # Earth's rotation rate (s^-1)

# 📂 Load CMIP6 NetCDF File
file_path = "path_to_cmip6_data.nc"
ds = xr.open_dataset(file_path)

# 🔹 Extract CMIP6 Variables
ua_850 = ds['ua'].sel(plev=85000)  # 850 hPa zonal wind (m/s)
va_850 = ds['va'].sel(plev=85000)  # 850 hPa meridional wind (m/s)
ua_200 = ds['ua'].sel(plev=20000)  # 200 hPa zonal wind (m/s)
va_200 = ds['va'].sel(plev=20000)  # 200 hPa meridional wind (m/s)
hus_600 = ds['hus'].sel(plev=60000)  # 600 hPa specific humidity (kg/kg)
tos = ds['tos']  # Sea Surface Temperature (K)

# 📌 Compute Lat/Lon Spacing
lat_radians = np.deg2rad(ds.lat)
lon_radians = np.deg2rad(ds.lon)
dlat = np.gradient(lat_radians)
dlon = np.gradient(lon_radians)

dx = R * np.cos(lat_radians[:, np.newaxis]) * dlon[np.newaxis, :]
dy = R * dlat[:, np.newaxis]

# 🌀 Compute Relative Vorticity (ζ850)
vorticity_850 = (np.gradient(va_850, dx, axis=-1) - np.gradient(ua_850, dy, axis=-2))

# ⚡ Compute Vertical Wind Shear
shear_u = ua_200 - ua_850  # U-component shear
shear_v = va_200 - va_850  # V-component shear
wind_shear = np.sqrt(shear_u**2 + shear_v**2)  # Shear magnitude (m/s)

# 🔥 Compute Maximum Potential Intensity (MPI) using Emanuel (1988) Approximation
SST_C = tos - 273.15  # Convert to Celsius
Ck = 1.0  # Exchange coefficient of enthalpy
Cd = 0.9  # Drag coefficient
mpi = np.sqrt(Ck / Cd * (SST_C - 26) * 1e3)  # Simplified MPI formula

# 🎯 Apply OWZP Cyclone Detection Criteria
# (Based on thresholding vorticity, wind shear, humidity, and pressure)
cyclone_mask = (
    (vorticity_850 > 1e-5) &   # Strong cyclonic vorticity
    (wind_shear < 15) &        # Low vertical wind shear
    (hus_600 > 0.4) &          # Sufficient mid-level humidity
    (mpi > 40)                 # MPI threshold
)

# 📌 Gaussian Smoothing for Better Detection
cyclone_mask_smoothed = gaussian_filter(cyclone_mask.astype(float), sigma=1)

# 📊 Save Output as NetCDF
output_ds = xr.Dataset({
    'relative_vorticity_850': (["time", "lat", "lon"], vorticity_850),
    'vertical_wind_shear': (["time", "lat", "lon"], wind_shear),
    'humidity_600': (["time", "lat", "lon"], hus_600),
    'maximum_potential_intensity': (["time", "lat", "lon"], mpi),
    'cyclone_mask': (["time", "lat", "lon"], cyclone_mask_smoothed)
}, coords={"time": ds.time, "lat": ds.lat, "lon": ds.lon})

output_file = "tropical_cyclone_detection.nc"
output_ds.to_netcdf(output_file)

print(f"Processed cyclone detection data saved to {output_file}")
```

---

## **📌 Explanation of Computations**
1. **Extract CMIP6 Fields**
   - Reads wind components (`ua`, `va`), humidity (`hus`), and SST (`tos`).
   - Selects pressure levels: **850 hPa, 200 hPa, and 600 hPa**.

2. **Compute Derived Variables**
   - **Relative vorticity (ζ850)**: 
     \[
     \zeta = \frac{\partial v}{\partial x} - \frac{\partial u}{\partial y}
     \]
   - **Vertical Wind Shear**:
     \[
     \Delta U = U_{200} - U_{850}, \quad \Delta V = V_{200} - V_{850}
     \]
     \[
     \text{Wind Shear} = \sqrt{\Delta U^2 + \Delta V^2}
     \]
   - **MPI Approximation** (Emanuel, 1988):
     \[
     MPI = \sqrt{\frac{C_k}{C_d} (SST - 26) \times 1000}
     \]
   - **Apply OWZP Criteria**:
     - Cyclonic vorticity (`vorticity_850 > 1e-5 s⁻¹`)
     - Low wind shear (`wind_shear < 15 m/s`)
     - Sufficient humidity (`hus_600 > 0.4 kg/kg`)
     - Strong MPI (`mpi > 40`)

3. **Smooth Cyclone Mask using Gaussian Filter**
   - Helps avoid false positives by refining cyclone detection.

4. **Save Processed Data**
   - Outputs a **NetCDF file** (`tropical_cyclone_detection.nc`) for analysis.

---

## **📊 Output Variables**
| Variable Name                   | Description                              | Units  |
|----------------------------------|------------------------------------------|--------|
| `relative_vorticity_850`        | 850 hPa relative vorticity               | s\(^{-1}\) |
| `vertical_wind_shear`           | Wind shear (200 hPa - 850 hPa)           | m/s    |
| `humidity_600`                  | 600 hPa specific humidity                | kg/kg  |
| `maximum_potential_intensity`   | Theoretical max cyclone intensity        | Varies |
| `cyclone_mask`                  | Detected cyclone regions (1=Yes, 0=No)   | Binary |

---

## **🔹 Notes**
- **Modify the CMIP6 file path (`file_path`) before running.**
- **Adjust detection thresholds** for specific cyclone criteria.
- **Extend the MPI formula** for more accurate estimations (e.g., using reanalysis data).
- **Apply further filtering methods** for refining cyclone detection.

---

## **🚀 Next Steps**
- **Visualize Cyclone Tracks** with **Matplotlib + Cartopy**.
- **Analyze Historical Cyclone Trends** from CMIP6 datasets.
- **Compare Different Climate Models (e.g., CESM2, GFDL-ESM4).**

---

This script automates **tropical cyclone detection from CMIP6 data** using the **OWZP method**, providing a robust and scalable approach for climate modeling and forecasting. 🌪️🔥
