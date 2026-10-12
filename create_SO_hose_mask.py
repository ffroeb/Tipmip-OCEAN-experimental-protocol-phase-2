## IMPORT MODULES
from scipy.ndimage import distance_transform_edt
import cartopy.feature as cfeature
from shapely.geometry import Point
from shapely.ops import unary_union
import warnings
import sys
import numpy as np
import pandas as pd
import xarray as xr
import dask
import matplotlib.pyplot as plt
import netCDF4 as nc
import matplotlib.path as mpath

from scipy.ndimage import binary_dilation
from scipy.spatial import cKDTree

from cmcrameri import cm
from matplotlib.colors import LinearSegmentedColormap, ListedColormap, BoundaryNorm
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.util import add_cyclic_point
from mpl_toolkits.basemap import Basemap
import iris

## LOAD GRID INFORMATION NorESM2
grid_file = 'grid_tnx1v4.nc'
distance = 300.0
radius   = 6371.0
ant_lat    = -56.0  

ds      = nc.Dataset(grid_file)
cplmask = np.array(ds.variables["cplmask"][:])   # (ny, nx) — 1=ocean, 0=land
plat    = np.array(ds.variables["plat"][:])
plon    = np.array(ds.variables["plon"][:])
parea    = np.array(ds.variables["parea"][:])
ds.close()

## FIND ANTARCTIC COASTAL BOUNDARY CELLS
ant_region = plat < ant_lat

# Land/ice mask in the Antarctic region
land_ant = (cplmask == 0) & ant_region 

# Neighboring cells land / ocean
ocean_neighbor_kernel = np.ones((3, 3), dtype=bool)
land_dilated = binary_dilation(land_ant, structure=ocean_neighbor_kernel)
ocean_ant    = (cplmask == 1) & ant_region

# Coastal land cells = land cells adjacent to at least one ocean cell
coast_land = land_ant & land_dilated & (binary_dilation(ocean_ant, structure=ocean_neighbor_kernel))

# KD-tree in 2-D Cartesian space
def latlon_to_xyz(lat_deg, lon_deg):
    """Unit-sphere 3-D coords."""
    lat = np.radians(lat_deg)
    lon = np.radians(lon_deg)
    x = np.cos(lat) * np.cos(lon)
    y = np.cos(lat) * np.sin(lon)
    z = np.sin(lat)
    return np.column_stack([x, y, z])

coast_xyz = latlon_to_xyz(plat[coast_land], plon[coast_land])
tree = cKDTree(coast_xyz)

# Ocean cells south of the cutoff
ocean_xyz  = latlon_to_xyz(plat[ocean_ant], plon[ocean_ant])

# Nearest-neighbour
chord_dist, _ = tree.query(ocean_xyz, workers=-1)

# Chord length → great-circle distance
gc_dist_km = 2.0 * radius * np.arcsin(np.clip(chord_dist / 2.0, 0, 1))

# create mask file
dist_field = np.full(cplmask.shape, np.nan)
dist_field[ocean_ant] = gc_dist_km

# exponetial decay
hosing_rate = 0.3 # Sv 
hosing_rate_area = hosing_rate *10**6 / parea # m/s

weights = np.where(dist_field <= distance, np.exp(-dist_field / distance), 0.0)
mask_SO_hose = hosing_rate_area * weights / weights.sum()

mask_SO_hose = np.where(np.isnan(mask_SO_hose), 0.0, mask_SO_hose)

# change unit to kg/m2/s
dens = 1025 # kg / m3
mask_SO_hose_flux = mask_SO_hose * dens 

## SAVE NETCDF
out_file = 'SO_hosing_mask.nc'
with nc.Dataset(grid_file, 'r') as ds_in, nc.Dataset(out_file, 'w') as ds_out:
    
    # Copy dimensions
    for name, dim in ds_in.dimensions.items():
        ds_out.createDimension(name, len(dim) if not dim.isunlimited() else None)
    
    # Copy coordinate variables (lat, lon, etc.) 
    for name in ['plat', 'plon', 'y', 'x']:  
        if name in ds_in.variables:
            var_in  = ds_in.variables[name]
            var_out = ds_out.createVariable(name, var_in.datatype, var_in.dimensions)
            var_out[:] = var_in[:]
            var_out.setncatts({a: var_in.getncattr(a) for a in var_in.ncattrs()})
    
    # Copy global attributes
    ds_out.setncatts({a: ds_in.getncattr(a) for a in ds_in.ncattrs()})

    # Write hosing field
    hose_var = ds_out.createVariable('SO_hosing', 'f4', ('y', 'x'), fill_value=np.nan)
    hose_var[:] = mask_SO_hose
    hose_var.units        = 'm s^-1'
    hose_var.long_name    = 'Southern Ocean hosing mask (exponential decay)'

    
out_file = 'SO_hosing_mask_flux.nc'
with nc.Dataset(grid_file, 'r') as ds_in, nc.Dataset(out_file, 'w') as ds_out:
    
    # Copy dimensions
    for name, dim in ds_in.dimensions.items():
        ds_out.createDimension(name, len(dim) if not dim.isunlimited() else None)
    
    # Copy coordinate variables (lat, lon, etc.) 
    for name in ['plat', 'plon', 'y', 'x']:  
        if name in ds_in.variables:
            var_in  = ds_in.variables[name]
            var_out = ds_out.createVariable(name, var_in.datatype, var_in.dimensions)
            var_out[:] = var_in[:]
            var_out.setncatts({a: var_in.getncattr(a) for a in var_in.ncattrs()})
    
    # Copy global attributes
    ds_out.setncatts({a: ds_in.getncattr(a) for a in ds_in.ncattrs()})

    # Write hosing field
    hose_var = ds_out.createVariable('SO_hosing flux', 'f4', ('y', 'x'), fill_value=np.nan)
    hose_var[:] = mask_SO_hose_flux
    hose_var.units        = 'kg m^-2 s^-1'
    hose_var.long_name    = 'Southern Ocean hosing mask (exponential decay)'

## PLOT HOSING MAP
from matplotlib.ticker import ScalarFormatter
theta = np.linspace(0, 2*np.pi, 100)
map_circle = mpath.Path(np.vstack([np.sin(theta), np.cos(theta)]).T * 0.5 + [0.5, 0.5])

# --- data -----------------------------------------------------------------
data_so  = np.where(mask_SO_hose, mask_SO_hose, np.nan)

vmin, vmax = 0,  2.4*10**-7
pcmkwargs = {'transform': ccrs.PlateCarree(), 'cmap': cm.devon_r,
             'vmin': vmin, 'vmax': vmax}


fig = plt.figure(figsize=(14, 14))  # square figure = round circle

ax = fig.add_subplot(1, 1, 1, projection=ccrs.SouthPolarStereo(central_longitude=0))
ax.set_extent([-180, 180, -90, -55], ccrs.PlateCarree())
ax.set_boundary(map_circle, transform=ax.transAxes)
pcm = ax.pcolormesh(plon, plat, data_so, **pcmkwargs)
ax.add_feature(cfeature.LAND, facecolor='grey', zorder=100)
ax.gridlines(draw_labels=True, x_inline=False, y_inline=True,
              color='k', linestyle='dashed', linewidth=0.5);
ax.add_feature(cfeature.COASTLINE)

ax.set_title('Antarctica hosing', fontsize=15, pad=17)

cbar = plt.colorbar(pcm, ax=[ax], shrink=0.5, orientation='horizontal', extend='max')
cbar.ax.get_yaxis().labelpad = 0.01
cbar.ax.set_xlabel('Freshwater input [ms$^{-1}$]',  fontsize=15)
cbar.ax.tick_params(labelsize=15)

fmt = ScalarFormatter(useMathText=True)
fmt.set_powerlimits((-7, -7))  # forces scientific notation at exactly 10^-7
cbar.ax.xaxis.set_major_formatter(fmt)
cbar.ax.xaxis.get_offset_text().set_fontsize(15)


