
# coding: utf-8

# In[1]:

## This program computes and plots simulated doppler radar imagery in both RHI and PPI scans. It has been written in its 
## final form by Austin Dixon, Kelton Halbert, and Leigh Orf at the University of Wisconsin-Madison.
## For questions please email awdixon2@wisc.edu.
 
## Many of the RHI plots and computations were completed with the 
## generous assistance of Alex Schueth at Texas Tech University (Alex.Schueth@ttu.edu).

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

print('Running full script that plots both PPI and RHI')


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


# In[9]:

#plotting rectangular RHI scans (takes about 45 seconds)

print('Creating rectangular RHI scans')

x0,y0,z0 = 2, -4, 0. #kilometers from the origin, aka where the radar is located

RHIlength=7. #length of the Radar beam from center for RHI scans
dx=25. #horizontal grid spacing of the model used to calculate how many "bins" used to interpolate to
azimuth=155. #azimuth angle used for RHI scans
npts=RHIlength/(dx/1000)

print('Running computations')
zgradtv, ygradtv, xgradtv=np.gradient(theta[:,:,:])
maggradtv=np.sqrt((xgradtv)**2.+(ygradtv)**2+(zgradtv)**2)

a.close()

dudy,dudx=np.array(np.gradient(u,dx,axis=(1,2)))
dvdy,dvdx=np.array(np.gradient(v,dx,axis=(1,2)))
dwdy,dwdx=np.array(np.gradient(w,dx,axis=(1,2)))
dtdy,dtdx=np.array(np.gradient(theta,dx,axis=(1,2)))

Delt=np.sqrt(dtdx**2+dtdy**2)#+dtdz**2)

adia2dfronto=(1./Delt)*(dtdx*(-dudx*dtdx-dvdx*dtdy)+dtdy*(-dudy*dtdx-dvdy*dtdy))*1000.
adia2dfronto=np.nan_to_num(adia2dfronto)#replaces nans with zero's so interpolation doesn't fail   

#takes care of the special cases of geometry in the radar quadrants for RHI scans.
if azimuth >=0. and azimuth < 90.:
    uproj=u*(1.-np.sin(np.deg2rad(azimuth)))
    vproj=v*(1.-np.cos(np.deg2rad(azimuth)))
    xvortproj=xvort*(1.-np.cos(np.deg2rad(azimuth)))
    yvortproj=yvort*(1.-np.sin(np.deg2rad(azimuth)))
if azimuth >= 90. and azimuth < 180.:
    uproj=-u*(1.-np.sin(np.deg2rad(azimuth)))
    vproj=v*(1.-np.cos(np.deg2rad(180-azimuth)))
    xvortproj=-xvort*(1.-np.cos(np.deg2rad(azimuth)))
    yvortproj=yvort*(1.-np.sin(np.deg2rad(180-azimuth)))
if azimuth >= 180. and azimuth < 270.:
    uproj=-u*(1.-np.sin(np.deg2rad(360-azimuth)))
    vproj=-v*(1.-np.cos(np.deg2rad(180-azimuth)))
    xvortproj=-xvort*(1.-np.cos(np.deg2rad(360-azimuth)))
    yvortproj=-yvort*(1.-np.sin(np.deg2rad(180-azimuth)))
if azimuth >= 270. and azimuth <360.:
    uproj=u*(1.-np.sin(np.deg2rad(360-azimuth)))
    vproj=-v*(1.-np.cos(np.deg2rad(azimuth)))
    xvortproj=xvort*(1.-np.cos(np.deg2rad(360-azimuth)))
    yvortproj=-yvort*(1.-np.sin(np.deg2rad(azimuth)))


# the following block of code computes the necessary variables for RHI scans as a rectangular cross section,
# meaning it is not in the conical shape of traditional radar scans.

xplot=np.linspace(0.,RHIlength,npts)

crosswind=np.zeros((len(z),int(npts)))
crossref=np.zeros((len(z),int(npts)))
crosstvp=np.zeros((len(z),int(npts)))
crosstvgrad=np.zeros((len(z),int(npts)))
crossfronto=np.zeros((len(z),int(npts)))
crosssvortmag=np.zeros((len(z),int(npts)))
crossppert=np.zeros((len(z),int(npts)))

#finds the final position for the cross section for both RHI and PPI scans
x1cross=x0+RHIlength*np.cos(np.deg2rad(azimuth))
y1cross=y0+RHIlength*np.sin(np.deg2rad(azimuth))

crossx = np.linspace(x0,x1cross,npts)
crossy = np.linspace(y0,y1cross,npts)
crossxinterp=interp1d(x,np.arange(len(x)))
crossyinterp=interp1d(y,np.arange(len(y)))
crossxk = crossxinterp(crossx)
crossyk = crossyinterp(crossy)

for k in range(0, len(z)): 
    crosswind[k,:]= scipy.ndimage.map_coordinates(uproj[k,:,:]+vproj[k,:,:], np.vstack((crossyk, crossxk))) #adds the projection of the u component on the desired plane and the projection of the v component on the desired plane together
    crosssvortmag[k,:]= scipy.ndimage.map_coordinates(svort[k,:,:], np.vstack((crossyk, crossxk))) 
    crossref[k,:]= scipy.ndimage.map_coordinates(ref[k,:,:], np.vstack((crossyk, crossxk)))
    crosstvp[k,:]= scipy.ndimage.map_coordinates(theta[k,:,:], np.vstack((crossyk, crossxk)))
    crosstvgrad[k,:]= scipy.ndimage.map_coordinates(maggradtv[k,:,:], np.vstack((crossyk, crossxk)))
    crossfronto[k,:]= scipy.ndimage.map_coordinates(adia2dfronto[k,:,:], np.vstack((crossyk, crossxk)))
    crossppert[k,:]= scipy.ndimage.map_coordinates(ppert[k,:,:], np.vstack((crossyk, crossxk)))

fig = plt.figure(figsize=(30, 15))

print('Creating rectangular plots')
# plots basic reflectivity from netcdf file for viewing radar location and cross section
ax1 = plt.subplot2grid((3,3),(0,0),rowspan=2)
radar=plt.contourf(x,y,ref[0,:,:],radar_levels,cmap=ref_new_cmap, norm=ref_norm)
plt.plot([x0,x1cross],[y0,y1cross],color="black")
plt.plot(x0,y0,'ko', markersize=8)
plt.plot(x0,y0,'wo', markersize=5)
plt.colorbar(radar, label="Radar Reflectivity [dBZ]", ticks=[-5,0,5,10,15,20,25,30,35,40,45,50,55,60,65])
plt.xlabel('Kilometers')
plt.ylabel('Kilometers')
plt.title('Model Time: '+str(int(time))+' seconds')

# plots streamwise vorticity 
ax2 = plt.subplot2grid((3,3),(2,0))
pvort=ax2.pcolormesh(xplot,z,crosssvortmag,vmin=0,vmax=0.1)
plt.colorbar(pvort, label="Streamwise Vorticity Magnitude [1/s]")#, ticks=[-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30])
plt.title('Streamwise Vorticity Magnitude')
plt.ylim(0, 4.5)
plt.xlim(0, int(RHIlength))
plt.grid(True)

# Reflectivity RHI scan
ax3 = plt.subplot2grid((3,3),(0,1))
pref = ax3.pcolormesh(xplot,z,ma.array(crossref,mask=[crossref<-5]),cmap=ref_new_cmap, norm=ref_norm)
plt.colorbar(pref, label="Radar Reflectivity [dBZ]", ticks=[-5,0,5,10,15,20,25,30,35,40,45,50,55,60,65])
plt.title('Reflectivity')
plt.ylim(0, 4.5)
plt.xlim(0, int(RHIlength))
plt.grid('on')

# Radial velocity RHI
ax4 = plt.subplot2grid((3,3),(1,1))
pvel=ax4.pcolormesh(xplot,z,crosswind,cmap=vel_cmap, vmin=-15, vmax=15)
plt.colorbar(pvel, label="Velocity [m/s]", ticks=[-15,-12.5,-10,-7.5,-5,-2.5,0,2.5,5,7.5,10,12.5,15])
plt.title('Radial Velocity')
plt.ylim(0, 4.5)
plt.xlim(0, int(RHIlength))
plt.grid(True)

# Pressure perturbation RHI
ax5 = plt.subplot2grid((3,3),(2,1))
pppert=ax5.pcolormesh(xplot,z,crossppert,cmap='hot',vmin=-20.,vmax=0.)
plt.colorbar(pppert, label="Pressure Perturbation [hPa]")#, ticks=[-0.05,-0.04,-0.03,-0.02,-0.01,0,0.01,0.02,0.03,0.04,0.05])
plt.title('Pressure Perturbation')
plt.ylim(0, 4.5)
plt.xlim(0, int(RHIlength))
plt.grid(True)

# Perturbation Potential Temperature
ax6 = plt.subplot2grid((3,3),(0,2))
ptv=ax6.pcolormesh(xplot,z,crosstvp,cmap=sftemp(),vmin=-5,vmax=1)
plt.colorbar(ptv, label="Perturbation Potential Temperature [K]")#, ticks=[-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30])
plt.title('Perturbation Potential Temperature')
plt.ylim(0, 4.5)
plt.xlim(0, int(RHIlength))
plt.grid(True)

# Theta' gradient
ax7 = plt.subplot2grid((3,3),(1,2))
ptvgrad=ax7.pcolormesh(xplot,z,crosstvgrad,cmap=cmocean.cm.matter,vmin=0,vmax=2)
plt.colorbar(ptvgrad, label="Perturbation Potential Temperature Gradient [K]")#, ticks=[-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30])
plt.title('Perturbation Potential Temperature Gradient')
plt.ylim(0, 4.5)
plt.xlim(0, int(RHIlength))
plt.grid(True)

# Adiabatic 2D Frontogenesis
ax8 = plt.subplot2grid((3,3),(2,2))
pfronto=ax8.pcolormesh(xplot,z,crossfronto,cmap='seismic',vmin=-0.25,vmax=0.25)
plt.colorbar(pfronto, label="Adiabatic 2D Frontogenesis"+r'$[\frac{K}{km \cdot s}]$', ticks=[-0.05,-0.04,-0.03,-0.02,-0.01,0,0.01,0.02,0.03,0.04,0.05])
plt.title('Adiabatic 2D Frontogenesis')
plt.ylim(0, 4.5)
plt.xlim(0, int(RHIlength))
plt.grid(True)


fig.suptitle('Simulated RHI     Azimuth: '+str(azimuth)+'     '+str(int(time))+' seconds', fontsize=20, y=1.00)
plt.subplots_adjust(top=0.8)
plt.tight_layout(pad=3.0)#fixes the aspect ratio some

ax1.set_rasterized(True)
ax2.set_rasterized(True)
ax3.set_rasterized(True)
ax4.set_rasterized(True)

ax5.set_rasterized(True)
ax6.set_rasterized(True)
ax7.set_rasterized(True)
ax8.set_rasterized(True)

plt.savefig('RHIstuff.png')
plt.close()
fig.clf()
gc.collect()

print(' Rectangular plots complete')
print("Plotting Time: "+str(datetime.now() - startTime))


# In[ ]:

# This block creates true RHI scans in the conical shape of real radar scans. 
# It is currently quite slow (~25 min).

print('Creating conical RHI scans')
class MidpointNormalize(colors.Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))

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

radar_levels = [x/10. for x in range(0,651)]
tvlevels = [x / 10. for x in range(0,51)]

startTime = datetime.now()

# x0,y0,z0 = 1., 4., 0. #kilometers from the origin, aka where the radar is located
z0=0.
# azimuth=120.
RHIlength=11. #length of the RHI
bins=1000. #This is the number of points along the line, "bins" if you will, bin length 10-15 m in the Ka
zenithres=0.3 #angle resolution in degrees, this is the correct resolution of the Ka
zenithmax=20.1 #max angle of the RHI
azimuthmin=150.
azimuthmax=160.
# startfile=0
#fil=250 #file number
# path='D:\\'
# path='G:\\'

def label_count(array_label):
    unique, counts = np.unique(array_label, return_counts=True)
    dicounts=(dict(zip(unique, counts)))
    dicountsorted = sorted(dicounts.items(), key=operator.itemgetter(1),reverse=True)
    idbiggest= [i[0] for i in dicountsorted[1:]]
    sizebiggest= [i[1] for i in dicountsorted[1:]]
    return np.array(idbiggest)

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx

centroidslice=np.empty([750,3])


filename = 'nc-ff.05550.000000.nc'
print(filename)
a = Dataset(filename, 'r')

time = a.variables["time"][0]

u = a.variables["uinterp"][0,:,:,:]+11.
v = a.variables["vinterp"][0,:,:,:]+4.
w = a.variables["winterp"][0,:,:,:]

ref=a.variables["dbz"][0,:,:,:]

zvort = a.variables["zvort"][0,:,:,:]
yvort = a.variables["yvort"][0,:,:,:]
xvort = a.variables["xvort"][0,:,:,:]
# vortmag = np.sqrt(xvort**2+yvort**2)

print('Begin computations')

wndspd=np.sqrt(u**2+v**2+w**2)
svortmag=(u*xvort+v*yvort+w*zvort)/wndspd

uh=w*zvort

x=a.variables["xh"][:]
y=a.variables["yh"][:]
z=a.variables["zh"][:]

theta = a.variables["thrhopert"][0,:,:,:]
ppert = a.variables["prespert"][0,:,:,:]/100.

zgradtvp, ygradtvp, xgradtvp=np.gradient(theta[:,:,:])
maggradtvp=np.sqrt((xgradtvp)**2.+(ygradtvp)**2+(zgradtvp)**2)

zgradtv, ygradtv, xgradtv=np.gradient(theta[:,:,:])
maggradtv=np.sqrt((xgradtv)**2.+(ygradtv)**2+(zgradtv)**2)

a.close()

#---------algo to find the location of the radar that is storm relative. should work for the most part, may need some tweaking
uh25=uh[find_nearest(z,2):find_nearest(z,5),:,:]
print(find_nearest(z,2))
labelsmaxw,numlabelsmaxw=ndimage.label(gaussian_filter(uh25,sigma=50)>(0.75*np.amax(gaussian_filter(uh25,sigma=50))))
biggest=label_count(labelsmaxw)
mask=np.zeros_like(labelsmaxw)
mask[labelsmaxw==biggest[0]]=1

centroids=ndimage.measurements.center_of_mass(mask)

x0=x[int(centroids[2]+52)]
y0=y[int(centroids[1])]

#-------------

for azimuth in range(int(azimuthmin),int(azimuthmax),10):

    #takes care of the special cases of geometry
    if azimuth >=0. and azimuth < 90.:
        uproj=u*(1.-np.sin(np.deg2rad(azimuth)))
        vproj=v*(1.-np.cos(np.deg2rad(azimuth)))
        xvortproj=xvort*(1.-np.cos(np.deg2rad(azimuth)))
        yvortproj=yvort*(1.-np.sin(np.deg2rad(azimuth)))
    if azimuth >= 90. and azimuth < 180.:
        uproj=-u*(1.-np.sin(np.deg2rad(azimuth)))
        vproj=v*(1.-np.cos(np.deg2rad(180-azimuth)))
        xvortproj=-xvort*(1.-np.cos(np.deg2rad(azimuth)))
        yvortproj=yvort*(1.-np.sin(np.deg2rad(180-azimuth)))
    if azimuth >= 180. and azimuth < 270.:
        uproj=-u*(1.-np.sin(np.deg2rad(360-azimuth)))
        vproj=-v*(1.-np.cos(np.deg2rad(180-azimuth)))
        xvortproj=-xvort*(1.-np.cos(np.deg2rad(360-azimuth)))
        yvortproj=-yvort*(1.-np.sin(np.deg2rad(180-azimuth)))
    if azimuth >= 270. and azimuth <360.:
        uproj=u*(1.-np.sin(np.deg2rad(360-azimuth)))
        vproj=-v*(1.-np.cos(np.deg2rad(azimuth)))
        xvortproj=xvort*(1.-np.cos(np.deg2rad(360-azimuth)))
        yvortproj=-yvort*(1.-np.sin(np.deg2rad(azimuth)))


    zenith = np.arange(zenithmax/zenithres+1)*zenithres
    xplot = np.zeros((len(zenith),int(bins)))
    zplot = np.zeros((len(zenith),int(bins)))
    rhiwind=np.zeros((len(zenith),int(bins)))
    rhiref=np.zeros((len(zenith),int(bins)))
    rhitvp=np.zeros((len(zenith),int(bins)))
    rhitvpgrad=np.zeros((len(zenith),int(bins)))
    rhifronto=np.zeros((len(zenith),int(bins)))
    rhiw=np.zeros((len(zenith),int(bins)))
    rhisvortmag=np.zeros((len(zenith),int(bins)))
    rhivortnorm=np.zeros((len(zenith),int(bins)))
    rhippert=np.zeros((len(zenith),int(bins)))

    for iteration, zen in enumerate(zenith):

        #the distance value of the end of the RHI
        x1=x0+(RHIlength*np.cos(np.deg2rad(azimuth))*np.cos(np.deg2rad(zen)))
        y1=y0+(RHIlength*np.sin(np.deg2rad(azimuth))*np.cos(np.deg2rad(zen)))
        z1=z0+(RHIlength*np.sin(np.deg2rad(zen)))

        #similar to the distance value arrays above, but these are to plot the data
        xplot[iteration,:]=np.linspace(0.,np.sqrt(np.abs(x1-x0)**2+np.abs(y1-y0)**2),bins)
        zplot[iteration,:]=np.linspace(z0,z1,bins)

        #interpolates the distances to the indices
        xinterp=interp1d(x,np.arange(len(x)), fill_value='extrapolate')
        yinterp=interp1d(y,np.arange(len(y)), fill_value='extrapolate')
        zinterp=interp1d(z,np.arange(len(z)),bounds_error=False,fill_value=0.0)#takes care of the situation where the lowest level is below the lowest grid point so it can't interpolate

        #find the decimal indices where the rhix values exist, these indices is where the image mapping is going to take place
        rhixi = xinterp(np.linspace(x0,x1,bins))
        rhiyi = yinterp(np.linspace(y0,y1,bins))
        rhizi = zinterp(np.linspace(z0,z1,bins))

        #first input is the full three dimensional data array
        #second input is the indices of the data array where the observations should be
        rhiwind[iteration,:]= ndimage.map_coordinates(uproj*np.cos(np.deg2rad(zen))+vproj*np.cos(np.deg2rad(zen))-w*np.sin(np.deg2rad(zen)), [rhizi,rhiyi,rhixi]) #adds the projection of the u component on the desired plane and the projection of the v component on the desired plane together
        rhiref[iteration,:]= ndimage.map_coordinates(ref[:,:,:], [rhizi,rhiyi,rhixi])
        rhitvp[iteration,:]= ndimage.map_coordinates(theta[:,:,:], [rhizi,rhiyi,rhixi])
        rhitvpgrad[iteration,:]= ndimage.map_coordinates(maggradtvp[:,:,:], [rhizi,rhiyi,rhixi])
        rhiw[iteration,:]= ndimage.map_coordinates(w[:,:,:], [rhizi,rhiyi,rhixi])
        rhisvortmag[iteration,:]= ndimage.map_coordinates(svortmag[:,:,:], [rhizi,rhiyi,rhixi])
        rhivortnorm[iteration,:]= ndimage.map_coordinates(xvortproj[:,:,:]+yvortproj[:,:,:], [rhizi,rhiyi,rhixi])
        rhippert[iteration,:]= ndimage.map_coordinates(ppert[:,:,:], [rhizi, rhiyi, rhixi])
        print(zen)

print('Creating plots')
fig = plt.figure(figsize=(24, 12))
ref_norm, ref_cmap = ctables.registry.get_with_range('NWSStormClearReflectivity', 0, 65)
ref_new_cmap = truncate_colormap(ref_cmap, 0.25,1.02)
vel_cmap = truncate_colormap(ctables.registry.get_colortable('Carbone42'), 0,1)

plt.rc('font', size=12)          # controls default text sizes
plt.rc('axes', titlesize=12)     # fontsize of the axes title
plt.rc('axes', labelsize=12)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=12)    # fontsize of the tick labels
plt.rc('ytick', labelsize=12)    # fontsize of the tick labels
plt.rc('legend', fontsize=12)    # legend fontsize
plt.rc('figure', titlesize=28)  # fontsize of the figure title


ax1 = plt.subplot2grid((3,3),(0,0),rowspan=2)
radar=plt.contourf(x,y,ref[0,:,:],radar_levels,cmap=ref_new_cmap, norm=ref_norm)
plt.plot([x0,x0+(RHIlength*np.cos(np.deg2rad(azimuth))*np.cos(np.deg2rad(zen)))],[y0,y0+(RHIlength*np.sin(np.deg2rad(azimuth))*np.cos(np.deg2rad(zen)))],color="black")
plt.plot(x0,y0,'ko', markersize=8)
plt.plot(x0,y0,'wo', markersize=5)
plt.colorbar(radar, label="Radar Reflectivity [dBZ]", ticks=[0,5,10,15,20,25,30,35,40,45,50,55,60,65])
plt.xlabel('Kilometers')
plt.ylabel('Kilometers')
plt.title('Model Time: '+str(int(time))+' seconds')
A=plt.text(0.05, 0.95,'a)', ha='center', va='center', transform=ax1.transAxes, fontsize=30, fontweight='bold')
A.set_path_effects([PathEffects.withStroke(linewidth=4,foreground='w')])

ax2 = plt.subplot2grid((3,3),(0,1))
ptv=ax2.pcolormesh(xplot,zplot,rhitvp,cmap=sftemp(),vmin=-5,vmax=0)
plt.colorbar(ptv, label=r'$\theta_v \backprime$ '+'[K]')#, ticks=[-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30])
plt.title('Perturbation Potential Temperature')
plt.ylim(0, 4.5)
plt.xlim(0, 11)
plt.grid('on')
A=plt.text(0.05, 0.9,'b)', ha='center', va='center', transform=ax2.transAxes, fontsize=30, fontweight='bold')
A.set_path_effects([PathEffects.withStroke(linewidth=4,foreground='w')])
plt.ylabel("Kilometers")

ax3 = plt.subplot2grid((3,3),(0,2))
pref = ax3.pcolormesh(xplot,zplot,ma.array(rhiref,mask=[rhiref<0]),cmap=ref_new_cmap, norm=ref_norm)
plt.colorbar(pref, label="Radar Reflectivity [dBZ]", ticks=[0,5,10,15,20,25,30,35,40,45,50,55,60,65])
plt.title('Reflectivity')
plt.ylim(0, 4.5)
plt.xlim(0, 11)
plt.grid('on')
A=plt.text(0.05, 0.9,'c)', ha='center', va='center', transform=ax3.transAxes, fontsize=30, fontweight='bold')
A.set_path_effects([PathEffects.withStroke(linewidth=4,foreground='w')])

ax4 = plt.subplot2grid((3,3),(1,1))
ptvgrad=ax4.pcolormesh(xplot,zplot,rhitvpgrad,cmap=cmocean.cm.tempo,vmin=0,vmax=1.5)
plt.colorbar(ptvgrad, label=r'$\nabla \theta_v [\frac{K}{km}]$')#, ticks=[-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30])
plt.title('Potential Temperature Gradient')
plt.ylim(0, 4.5)
plt.xlim(0, 11)
plt.grid('on')
A=plt.text(0.05, 0.9,'d)', ha='center', va='center', transform=ax4.transAxes, fontsize=30, fontweight='bold')
A.set_path_effects([PathEffects.withStroke(linewidth=4,foreground='w')])
plt.ylabel("Kilometers")

ax5 = plt.subplot2grid((3,3),(1,2))
pw=ax5.pcolormesh(xplot,zplot,rhiw,vmin=-5,vmax=25,norm=MidpointNormalize(midpoint=0.), cmap='seismic')
plt.colorbar(pw, label="W "+r'$[\frac{m}{s}]$')
plt.title('Vertical Velocity')
plt.ylim(0, 4.5)
plt.xlim(0, 11)
plt.grid('on')
A=plt.text(0.05, 0.9,'e)', ha='center', va='center', transform=ax5.transAxes, fontsize=30, fontweight='bold')
A.set_path_effects([PathEffects.withStroke(linewidth=4,foreground='w')])

ax6 = plt.subplot2grid((3,3),(2,0))
pvort=ax6.pcolormesh(xplot,zplot,rhisvortmag,vmin=0,vmax=0.1)
plt.colorbar(pvort, label="Streamwise Vorticity "+r'$[\frac{1}{s}]$')#, ticks=[-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30])
plt.title('Streamwise Vorticity Magnitude')
plt.ylim(0, 4.5)
plt.xlim(0, 11)
plt.grid('on')
A=plt.text(0.05, 0.9,'f)', ha='center', va='center', transform=ax6.transAxes, fontsize=30, fontweight='bold')
A.set_path_effects([PathEffects.withStroke(linewidth=4,foreground='w')])
plt.ylabel("Kilometers")
plt.xlabel("Kilometers")

ax7 = plt.subplot2grid((3,3),(2,1))
pvel=ax7.pcolormesh(xplot,zplot,rhiwind,cmap=vel_cmap, vmin=-25, vmax=25)
plt.colorbar(pvel, label="Velocity "+r'$[\frac{m}{s}]$', ticks=[-25,-20,-15,-10,-5,0,5,10,15,20,25])
plt.title('Storm-Relative Radial Velocity')
plt.ylim(0, 4.5)
plt.xlim(0, 11)
plt.grid('on')
A=plt.text(0.05, 0.9,'g)', ha='center', va='center', transform=ax7.transAxes, fontsize=30, fontweight='bold')
A.set_path_effects([PathEffects.withStroke(linewidth=4,foreground='w')])
plt.xlabel("Kilometers")

ax8 = plt.subplot2grid((3,3),(2,2))
pppert=ax8.pcolormesh(xplot,zplot,rhippert,cmap='hot',vmin=-3.,vmax=0.)
plt.colorbar(pppert, label="Pressure Perturbation [hPa]")#, ticks=[-0.05,-0.04,-0.03,-0.02,-0.01,0,0.01,0.02,0.03,0.04,0.05])
plt.title('Pressure Perturbation')
plt.ylim(0, 4.5)
plt.xlim(0, 11)
plt.grid('on')
A=plt.text(0.05, 0.9,'h)', ha='center', va='center', transform=ax8.transAxes, fontsize=30, fontweight='bold')
A.set_path_effects([PathEffects.withStroke(linewidth=4,foreground='w')])
plt.xlabel("Kilometers")


fig.suptitle('Simulated RHI     Azimuth: '+str(azimuth)+'     '+str(int(time))+' seconds', fontsize=20, y=1.00)
plt.subplots_adjust(top=0.8)
plt.tight_layout(pad=3.0)#fixes the aspect ratio some

ax1.set_rasterized(True)
ax2.set_rasterized(True)
ax3.set_rasterized(True)
ax4.set_rasterized(True)
ax5.set_rasterized(True)
ax6.set_rasterized(True)
ax7.set_rasterized(True)
ax8.set_rasterized(True)

plt.savefig('RHINewStuff.png')
plt.close()
fig.clf()
gc.collect()

print('Conical plots complete')
print("Plotting Time: "+str(datetime.now() - startTime))

print('Done!')

