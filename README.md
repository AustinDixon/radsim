# radsim
Python scripts for simulating doppler radar data from netcdf files of simulated thunderstorms.

This likely will not work out of the box with every netcdf file. Variable names will need to be changed accordingly as well as the name of your working file. 

x0, y0, and z0 currently need to be directly modified in the code for RHI scans. 

In the PPI plotter blocks, the "my_radar" line calls the radar class that compiles the algorithms. This line will need to be edited for the desired radar specifications. Specifically, radar location relative to origin, range max/min and resolution (step), azimuth max/min and resolution (step), elevation (zenith) max/min and resolution (step).

Note that PPI plots currently only show winds accounting for U and V. The W component of the wind will be added shortly.
