#print __doc__

# Author: Scott Collis (scollis@anl.gov)
# License: BSD 3 clause

# History: Scott Collis created.

import numpy as np
import matplotlib.pyplot as plt
import pyart

# read in the file, create a RadarMapDisplay object
filename = '/home/jruiz/TMP/SCALE_TO_RADAR/radar001_0001.nc'
radar = pyart.io.read(filename)
display = pyart.graph.RadarMapDisplay(radar)

#Obtengo la maxima velocidad radial segun la estrategia de escaneo.
nyquistv=radar.get_nyquist_vel(0,check_uniform=True)

#Obtengo la cantidad de niveles verticales.
nlevels=np.size(radar.sweep_number['data'])

#Obtengo la fecha
fecha=radar.time['units']
fecha=fecha[14:18]+fecha[19:21]+fecha[22:24]+fecha[25:27]+fecha[28:30]+fecha[31:33]

lon_radar=radar.longitude['data'][0]
lat_radar=radar.latitude['data'][0]
mw=3.0

# plot the second  first tilt (0, second argument)

for ilev in range(0,nlevels):

	display.plot_ppi_map('V_model', ilev, vmin=-2.0*nyquistv, vmax=2.0*nyquistv ,
                     min_lon=lon_radar-mw, max_lon=lon_radar+mw, min_lat=lat_radar-mw, max_lat=lat_radar+mw,
                     lon_lines=np.arange(lon_radar-mw,lon_radar+mw, .5), projection='lcc',
                     lat_lines=np.arange(lat_radar-mw,lat_radar+mw, .5), resolution='h',
                     lat_0=radar.latitude['data'][0],
                     lon_0=radar.longitude['data'][0])

	# plot range rings at 10, 20, 30 and 40km
	display.plot_range_ring(60. , line_style='k-')
	display.plot_range_ring(120., line_style='k--')
	display.plot_range_ring(180., line_style='k-')
	display.plot_range_ring(240., line_style='k--')

	# plots cross hairs
	display.plot_line_xy(np.array([-40000.0, 40000.0]), np.array([0.0, 0.0]),
        line_style='k-')
	display.plot_line_xy(np.array([0.0, 0.0]), np.array([-20000.0, 200000.0]),
        line_style='k-')

	# Indicate the radar location with a point
	display.plot_point(radar.longitude['data'][0], radar.latitude['data'][0])

	plt.savefig('./test/V_model_' + fecha + 'elev' + str(ilev)  + '.png')

	plt.close()



#plt.show()



for ilev in range(0,nlevels):

        display.plot_ppi_map('dBZ_model', ilev, vmin=-20, vmax=70 ,
                     min_lon=lon_radar-mw, max_lon=lon_radar+mw, min_lat=lat_radar-mw, max_lat=lat_radar+mw,
                     lon_lines=np.arange(lon_radar-mw,lon_radar+mw, .5), projection='lcc',
                     lat_lines=np.arange(lat_radar-mw,lat_radar+mw, .5), resolution='h',
                     lat_0=radar.latitude['data'][0],
                     lon_0=radar.longitude['data'][0])

        # plot range rings at 10, 20, 30 and 40km
        display.plot_range_ring(60. , line_style='k-')
        display.plot_range_ring(120., line_style='k--')
        display.plot_range_ring(180., line_style='k-')
        display.plot_range_ring(240., line_style='k--')

        # plots cross hairs
        display.plot_line_xy(np.array([-40000.0, 40000.0]), np.array([0.0, 0.0]),
        line_style='k-')
        display.plot_line_xy(np.array([0.0, 0.0]), np.array([-20000.0, 200000.0]),
        line_style='k-')

        # Indicate the radar location with a point
        display.plot_point(radar.longitude['data'][0], radar.latitude['data'][0])

        plt.savefig('./test/Dbz_model_' + fecha + 'elev' + str(ilev)  + '.png')

        plt.close()


