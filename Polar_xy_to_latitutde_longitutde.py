# %%
# Use conda/analysis3-24.07 env
from polar_convert.constants import SOUTH
import netCDF4 as nc
import numpy as np

# %%
from polar_convert import polar_xy_to_lonlat

# %%
true_scale_lat = -70  # true-scale latitude in degrees
re = 6378.137  # earth radius in km
e = 0.01671 # earth eccentricity
hemisphere = 'south'

# %%
# Open the NetCDF file
dataset = nc.Dataset('/g/data/tm70/ek4684/Bedmachine_data/BedMachineAntarctica_2020-07-15_v02.nc', 'r')
# Extract the latitude and longitude arrays
y = dataset.variables['y'][:]/100.0  # X array (Convert meters to KMs)
x = dataset.variables['x'][:]/100.0  # Y array (Convert meters to KMs)
# Create 2D coordinate grids
X, Y = np.meshgrid(x, y)

# %%
# Initialize lists to store new latitudes and longitudes
# Initialize 2D arrays for new latitudes and longitudes
new_longitudes = np.zeros_like(X)
new_latitudes = np.zeros_like(Y)

# Loop through all (i, j) combinations
for i in range(x.shape[0]):
    for j in range(x.shape[0]):
        lon_result, lat_result = polar_xy_to_lonlat(X[i, j], Y[i, j], true_scale_lat, re, e, hemisphere)
        lon_result = (lon_result + 180) % 360 - 180  # Ensure longitude is in [-180, 180]
        
        new_longitudes[i, j] = lon_result
        new_latitudes[i, j] = lat_result
        print(lon_result,lat_result)
# %%
# # Initialize lists to store new latitudes and longitudes
# # Initialize 2D arrays for new latitudes and longitudes
# new_longitudes = np.zeros_like(X)
# new_latitudes = np.zeros_like(Y)

# # Vectorize the function to apply on entire arrays at once
# vectorized_polar_xy_to_lonlat = np.vectorize(polar_xy_to_lonlat)

# # Apply the function to all elements of X and Y at once
# lon_results, lat_results = vectorized_polar_xy_to_lonlat(X, Y, true_scale_lat, re, e, hemisphere)

# # Ensure longitude is in [-180, 180]
# lon_results = (lon_results + 180) % 360 - 180  

# # Store the results
# new_longitudes = lon_results
# new_latitudes = lat_results

# %%
lon_result

# %%
# Define new file path
new_file_path = '/g/data/tm70/ek4684/Bedmachine_data/Updated_BedMachineAntarctica_2020-07-15_v02.nc'
new_dataset = nc.Dataset(new_file_path, 'w', format='NETCDF4')

# Copy dimensions from the original dataset to the new file
for name, dimension in dataset.dimensions.items():
    new_dataset.createDimension(name, len(dimension) if not dimension.isunlimited() else None)

# Copy variables (excluding x and y) from the original dataset to the new file
for name, variable in dataset.variables.items():
    if name not in ['x', 'y']:  # Exclude 'x' and 'y', as we are replacing them with lat/lon
        new_var = new_dataset.createVariable(name, variable.datatype, variable.dimensions)
        new_var[:] = variable[:]
        # Copy all attributes except the fill value
        for attr_name in variable.ncattrs():
            if attr_name != '_FillValue':
                new_var.setncattr(attr_name, variable.getncattr(attr_name))

# Create new lat and lon variables with the correct dimensions
lat_var = new_dataset.createVariable('lat', 'f4', ('y', 'x'))  # Latitude (2D)
lon_var = new_dataset.createVariable('lon', 'f4', ('y', 'x'))  # Longitude (2D)

# Assign the computed latitudes and longitudes to the new variables
lat_var[:, :] = new_latitudes  # Ensure latitudes are assigned correctly as 2D
lon_var[:, :] = new_longitudes  # Ensure longitudes are assigned correctly as 2D

# Set attributes for the new variables
lat_var.units = "degrees_north"
lon_var.units = "degrees_east"

# Close the new dataset
new_dataset.close()

# Close the original dataset
dataset.close()

print(f"New file saved as: {new_file_path}")

# %%
print(new_longitudes)

# %%
dataset = nc.Dataset('/g/data1a/tm70/ek4684/Bedmachine_data/NSIDC0771_LatLon_PS_N6.25km_v1.0.nc', 'r')

# %%
latitude = dataset.variables['latitude'][:,:]
longitude = dataset.variables['longitude'][:,:]

# %%
latitude[-1,:]

# %%
