
# coding: utf-8

# In[1]:

import matplotlib
matplotlib.use('agg')
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import Normalize
import matplotlib.colors as colors
import numpy.ma as ma
from metpy.plots import ctables
import scipy.ndimage as ndimage
from scipy.interpolate import interp1d
from datetime import datetime
import cmocean
import operator
import gc
import scipy
from scipy.ndimage import gaussian_filter
import matplotlib.patheffects as PathEffects
from scipy.interpolate import spline
from scipy.interpolate import RegularGridInterpolator
get_ipython().magic('matplotlib inline')

print('Running script for PPI only')


# In[2]:

#Creating colormaps for use in plots
def sftemp():
    sfc_cdict ={'red':      ((0.00, 0.20, 0.20),
                             (0.08, 0.40, 0.40),
                             (0.17, 0.27, 0.27),
                             (0.25, 0.80, 0.80),
                             (0.33, 0.20, 0.20),
                             (0.42, 0.20, 0.20),
                             (0.50, 0.00, 0.00),
                             (0.58, 0.99, 0.99),
                             (0.67, 1.00, 1.00),
                             (0.75, 0.82, 0.82),
                             (0.83, 0.53, 0.53),
                             (0.92, 0.95, 0.95),
                             (1.00, 1.00, 1.00)),
        
        'green':        ((0.00, 0.20, 0.20),
                         (0.08, 0.40, 0.40),
                         (0.17, 0.00, 0.00),
                         (0.25, 0.60, 0.60),
                         (0.33, 0.40, 0.40),
                         (0.42, 0.60, 0.60),
                         (0.50, 0.39, 0.39),
                         (0.58, 0.76, 0.76),
                         (0.67, 0.36, 0.36),
                         (0.75, 0.02, 0.02),
                         (0.83, 0.00, 0.00),
                         (0.92, 0.03, 0.03),
                         (1.00, 0.60, 0.60)),
            
            'blue':         ((0.00, 0.60, 0.60),
                             (0.08, 0.60, 0.60),
                             (0.17, 0.65, 0.65),
                             (0.25, 1.00, 1.00),
                             (0.33, 1.00, 1.00),
                             (0.42, 0.40, 0.40),
                             (0.50, 0.07, 0.07),
                             (0.58, 0.02, 0.02),
                             (0.67, 0.00, 0.00),
                             (0.75, 0.01, 0.01),
                             (0.83, 0.00, 0.00),
                             (0.92, 0.52, 0.52),
                             (1.00, 0.80, 0.80))}
                

    sfc_coltbl = LinearSegmentedColormap('SFC_COLTBL',sfc_cdict)
    return sfc_coltbl
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap
def discrete_cmap(N, base_cmap=None):
    """Create an N-bin discrete colormap from the specified input map"""

    # Note that if base_cmap is a string or None, you can simply do
    #    return plt.cm.get_cmap(base_cmap, N)
    # The following works for string, None, or a colormap instance:

    base = plt.cm.get_cmap(base_cmap)
    color_list = base(np.linspace(0, 1, N))
    cmap_name = base.name + str(N)
    return base.from_list(cmap_name, color_list, N)

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx


# In[3]:

vel_cmap = truncate_colormap(ctables.registry.get_colortable('Carbone42'), 0,1)
def reverse_colourmap(cmap, name = 'my_cmap_r'): #reverse vel_cmap if desired
    """
    In: 
    cmap, name 
    Out:
    my_cmap_r

    Explanation:
    t[0] goes from 0 to 1
    row i:   x  y0  y1 -> t[0] t[1] t[2]
                   /
                  /
    row i+1: x  y0  y1 -> t[n] t[1] t[2]

    so the inverse should do the same:
    row i+1: x  y1  y0 -> 1-t[0] t[2] t[1]
                   /
                  /
    row i:   x  y1  y0 -> 1-t[n] t[2] t[1]
    """        
    reverse = []
    k = []   

    for key in cmap._segmentdata:    
        k.append(key)
        channel = cmap._segmentdata[key]
        data = []

        for t in channel:                    
            data.append((1-t[0],t[2],t[1]))            
        reverse.append(sorted(data))    

    LinearL = dict(zip(k,reverse))
    my_cmap_r = matplotlib.colors.LinearSegmentedColormap(name, LinearL) 
    return my_cmap_r


vel_cmap_r = reverse_colourmap(vel_cmap, name = 'vel_cmap_r')


# In[4]:

startTime = datetime.now() #to time how long it takes the script to run.

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx

for enum,f in enumerate(range(450,451)): #looping through files
    #-----------------------------------------------------------------------------------

    radar_levels = [x/10. for x in range(-50,651)]
    tvlevels = [x / 10. for x in range(0,51,1)]

    filename = 'nc-ff.05550.000000.nc'
    print(filename)
    a = Dataset(filename, 'r')

print(a.variables.keys())    


# In[9]:

#setting variables from netcdf file

time = a.variables["time"][0]

u = a.variables["uinterp"][0,:,:,:] + 15.2 #use actual winds not storm-relative winds 
v = a.variables["vinterp"][0,:,:,:] + 10.5
w = a.variables["winterp"][0,:,:,:]

u_sr = a.variables["uinterp"][0,:,:,:]
v_sr = a.variables["vinterp"][0,:,:,:]

zvort = a.variables["zvort"][0,:,:,:]
yvort = a.variables["yvort"][0,:,:,:]
xvort = a.variables["xvort"][0,:,:,:]

x=a.variables["xh"][:]
y=a.variables["yh"][:]
z=a.variables["zh"][:] 

wndspd=np.sqrt(u**2+v**2+w**2)
svort=a.variables['streamvort'][0,:,:,:]

ref=a.variables["dbz"][0,:,:]
theta = a.variables["thrhopert"][0,:,:,:]
ppert = a.variables["prespert"][0,:,:,:]


# In[10]:

# New and improved PPI routine. Light years faster than previous iteration. The following class takes care of
# all calculations
ref_norm, ref_cmap = ctables.registry.get_with_range('NWSStormClearReflectivity', -5, 65)
ref_new_cmap = truncate_colormap(ref_cmap, 0.25,1.02)

x0,y0,z0 = 2, -4, 0. #kilometers from the origin, aka where the radar is located

class Radar(object):
    
    def __init__(self, **kwargs):
        """
        x0: radar x origin in km
        y0: radar y origin in km
        z0: radar z origin in km
        
        """
        self.x0 = kwargs.get("x0", 2) #default radar location relative to origin
        self.y0 = kwargs.get("y0", -4)
        self.z0 = kwargs.get("z0", 0)
        
        self.range_min = kwargs.get("rmin", 0.5) #default minimum radar range (cone of silence beginning)
        self.range_max = kwargs.get("rmax", 10.) #default max radar range
        self.range_step = kwargs.get("rstep", 0.175) #default range step
        
        self.az_min = kwargs.get("azmin", 0) #default starting azimuth
        self.az_max = kwargs.get("azmax", 360) #default end azimuth
        self.az_step = kwargs.get("azstep", 1) #default azimuth step
        
        self.elev_min = kwargs.get("elevmin", 0.5) #defualt minimum elevation
        self.elev_max = kwargs.get("elevmax", 60.5) #default max elevation
        self.elev_step = kwargs.get("elevstep", 1) #default elevation step
        
        self.range = np.arange(self.range_min, self.range_max+self.range_step, self.range_step)
        self.azimuth = np.arange(self.az_min, self.az_max+self.az_step, self.az_step)
        self.elevation = np.arange(self.elev_min, self.elev_max+self.elev_step, self.elev_step)
    
    def get_radar_ppi_grid(self, elevation=0):
        ## create a 2D polar coordinate grid
        ran, az = np.meshgrid(self.range, self.azimuth)

        ## get the cartesian points of the 
        ## polar grid for interpolation purposes
        xrad = x0+ran*np.cos(np.deg2rad(az))
        yrad = y0+ran*np.sin(np.deg2rad(az))
        ## fancy handling for the elevation - 
        ## compute the z coordinate for each elevation angle
        ## and stack it into a 3D array. However, for plotting/interpolating, we just 
        ## need a single index for the PPI.
        zrad = []
        for elev_idx in range(len(self.elevation)):
            zrad.append(ran*np.tan(np.deg2rad(self.elevation[elev_idx])))
        ## the zero index is for the first elevation. Change for higher elevation
        ## angles
        zrad = np.array(zrad)[elevation, :, :]
        
        return ran, az, xrad, yrad, zrad
    
    def interp_reflectivity_to_ppi(self, x, y, z, ref, elevation=0):
        ## pass the points to the interpolator in a 
        ## way it understands.
        interp3d = RegularGridInterpolator((z, y, x), ref, method="linear", bounds_error=False, fill_value=None)
        ran, az, xrad, yrad, zrad = self.get_radar_ppi_grid(elevation=elevation)
        points = [zrad, yrad, xrad]
        flat = np.array([m.flatten() for m in points])
        out_array_ref = interp3d(flat.T)
        refl_ppi = out_array_ref.reshape(*points[0].shape)
        
        return ran, az, refl_ppi
    
    def get_radial_velocity(self, interp_u, interp_v, interp_w, az, elev):
        """Get the radial velocity from the U and V components pre-interpolated
        onto our polar coordinate grid. Must have already been passed to the RegularGridInterpolator
        so that it's on the radar grid."""
        uproj = np.zeros(interp_u.shape)
        vproj = np.zeros(interp_v.shape)

        ## get the quadrant indices for maximum
        ## interpolation efficiency.
        quad1_idxs = np.where((az >= 0) & (az < 90))
        quad2_idxs = np.where((az >= 90.) & (az < 180))
        quad3_idxs = np.where((az >= 180.) & (az < 270))
        quad4_idxs = np.where((az >= 270.) & (az < 360))

        ## quadrant 1
        uproj[quad1_idxs] = interp_u[quad1_idxs]*(1.-np.sin(np.deg2rad(az[quad1_idxs])))
        vproj[quad1_idxs] = interp_v[quad1_idxs]*(1.-np.cos(np.deg2rad(az[quad1_idxs])))

        ## quadrant 2
        uproj[quad2_idxs] = -interp_u[quad2_idxs]*(1.-np.sin(np.deg2rad(az[quad2_idxs])))
        vproj[quad2_idxs] = interp_v[quad2_idxs]*(1.-np.cos(np.deg2rad(180-az[quad2_idxs])))

        ## quadrant 3
        uproj[quad3_idxs] = -interp_u[quad3_idxs]*(1.-np.sin(np.deg2rad(360-az[quad3_idxs])))
        vproj[quad3_idxs] = -interp_v[quad3_idxs]*(1.-np.cos(np.deg2rad(180-az[quad3_idxs])))

        ## quadrant 4
        uproj[quad4_idxs] = interp_u[quad4_idxs]*(1.-np.sin(np.deg2rad(360-az[quad4_idxs])))
        vproj[quad4_idxs] = -interp_v[quad4_idxs]*(1.-np.cos(np.deg2rad(az[quad4_idxs])))
        
        wproj = interp_w*(1.-np.cos(np.deg2rad(90-elev)))

        Vr = uproj + vproj + wproj
        return Vr
    
    def interp_radial_velocity_to_ppi(self, x, y, z, u, v, w, elevation=0):
        interp3d_u = RegularGridInterpolator((z, y, x), u, method="linear", bounds_error=False, fill_value=None)
        interp3d_v = RegularGridInterpolator((z, y, x), v, method="linear", bounds_error=False, fill_value=None)
        interp3d_w = RegularGridInterpolator((z, y, x), w, method="linear", bounds_error=False, fill_value=None)
        ran, az, xrad, yrad, zrad = self.get_radar_ppi_grid(elevation=elevation)
        
        points = [zrad, yrad, xrad]
        flat = np.array([m.flatten() for m in points])
        out_array_u = interp3d_u(flat.T)
        out_array_v = interp3d_v(flat.T)
        out_array_w = interp3d_w(flat.T)
        u_ppi = out_array_u.reshape(*points[0].shape)
        v_ppi = out_array_v.reshape(*points[0].shape)
        w_ppi = out_array_w.reshape(*points[0].shape)

        Vr = self.get_radial_velocity(u_ppi, v_ppi, w_ppi, az, elevation)
        return ran, az, Vr
        


# In[11]:

#plotter for reflectivity PPI

## construct the figure

print('Scanning for reflectivity')
fig = plt.figure(figsize=(19.5, 15))
ax = plt.gca()
plt.subplot(projection='polar')

## create and initialize our radar
my_radar = Radar(x0=4, y0=-5, azmin=0, azmax=360, rmin=2.1, rstep=0.0375, azstep=1)

# get the reflectivity plot
ran, az, refl_ppi = my_radar.interp_reflectivity_to_ppi(x, y, z, ref, elevation=0)

## plot the radar PPI
plt.pcolormesh(np.deg2rad(az), ran, refl_ppi, vmin=0, vmax=75, cmap=ref_new_cmap)

plt.colorbar()
plt.savefig('ppirefgrel.png')
plt.grid()
plt.show()

print('Reflectivity scan complete')


# In[16]:

#plotter for radial velocity PPI

print('Scanning for radial velocity')

## construct the figure
fig = plt.figure(figsize=(19.5, 15))
ax = plt.gca()
plt.subplot(projection='polar')

## create and initialize our radar
my_radar = Radar(x0=4, y0=-5, azmin=0, azmax=360, rmin=2.1, rstep=0.0375, azstep=1)

## get the radial velocity plot
ran, az, Vr = my_radar.interp_radial_velocity_to_ppi(x, y, z, u, v, w, elevation=0)

## plot the radar PPI
plt.pcolormesh(np.deg2rad(az), ran, Vr, vmin=-60, vmax=60, cmap=vel_cmap)

plt.colorbar()
plt.title('Ground-Relative Radial Velocity')
plt.savefig('ppiGRV.png')
plt.grid()
plt.show()

print('Radial velocity scan complete')
print('Done')


# In[17]:

#plotter for radial velocity PPI

print('Scanning for storm-relative radial velocity')

## construct the figure
fig = plt.figure(figsize=(19.5, 15))
ax = plt.gca()
plt.subplot(projection='polar')

## create and initialize our radar
my_radar = Radar(x0=4, y0=-5, azmin=0, azmax=360, rmin=2.1, rstep=0.0375, azstep=1)

## get the radial velocity plot
ran, az, Vr = my_radar.interp_radial_velocity_to_ppi(x, y, z, u_sr, v_sr, w, elevation=0)

## plot the radar PPI
plt.pcolormesh(np.deg2rad(az), ran, Vr, vmin=-60, vmax=60, cmap=vel_cmap)

plt.colorbar()
plt.title('Storm-realtive Radial Velocity')
plt.savefig('ppiSRV.png')
plt.grid()
plt.show()

print('Radial velocity scan complete')
print('Done')


# In[193]:

elev_min = 0.5 #defualt minimum elevation
elev_max = 60.5 #default max elevation
elev_step = 1 #default elevation step       
elevation = np.arange(elev_min, elev_max+elev_step, elev_step, dtype=float)  
print(elevation[1])
    
xrad = x0+ran*np.cos(np.deg2rad(az))
yrad = y0+ran*np.sin(np.deg2rad(az))
zrad = []
for elev_idx in range(len(elevation)):
    zrad.append(ran*np.tan(np.deg2rad(elevation[elev_idx])))
zrad = np.array(zrad)


# In[208]:

ppinetcdf = Dataset('ppiscan.nc', 'w', format='NETCDF4_CLASSIC')
rang = ppinetcdf.createDimension('range', ran.shape[1])
azimuth = ppinetcdf.createDimension('azimuth', ran.shape[0])

ranges = ppinetcdf.createVariable('ran', np.float32, ('azimuth','range'))
azimuths = ppinetcdf.createVariable('az', np.float32, ('azimuth', 'range'))
reflectivity = ppinetcdf.createVariable('reflectivity', np.float32, ('azimuth', 'range'))
velocities = ppinetcdf.createVariable('radialvels', np.float32, ('azimuth','range'))

cartx = ppinetcdf.createVariable('x', np.float32, ('azimuth', 'range'))
carty = ppinetcdf.createVariable('y', np.float32, ('azimuth', 'range'))
cartz = ppinetcdf.createVariable('z', np.float32, ('azimuth', 'range'))
xnot = ppinetcdf.createVariable('x0', np.float32)
ynot = ppinetcdf.createVariable('y0', np.float32)
znot = ppinetcdf.createVariable('z0', np.float32)
rangemin = ppinetcdf.createVariable('ranmin', np.float32)
rangemax = ppinetcdf.createVariable('ranmax', np.float32)
rangestep = ppinetcdf.createVariable('ranstep', np.float32)
azmin = ppinetcdf.createVariable('azmin', np.float32)
azmax = ppinetcdf.createVariable('azmax', np.float32)
azstep = ppinetcdf.createVariable('azstep', np.float32)
zenith = ppinetcdf.createVariable('zenith', np.float32)

ppinetcdf.set_auto_mask(False)
reflectivity[:] = refl_ppi
ranges[:] = ran
azimuths[:] = az
velocities[:] = Vr
cartx[:] = xrad
carty[:] = yrad
cartz[:] = zrad[0,:,:]
xnot[:] = 4
ynot[:] = -5
znot[:] = 0

rangemin[:] = 2.1
rangemax[:] = 10
rangestep[:] = 0.0375
azmin[:] = 0
azmax[:] = 360
azstep[:] = 1
zenith[:] = 0.5




ppinetcdf.close()

