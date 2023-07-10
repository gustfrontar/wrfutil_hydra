from wrf import getvar, ll_to_xy, CoordPair, GeoBounds
import numpy as np
from netCDF4 import Dataset
import glob
from copy import deepcopy
import sys

PATH_SIHVIGILA = str(sys.argv[1])

ll_lat = -38
ll_lon = -63
ur_lat = -32
ur_lon = -55

ll = CoordPair(lat = ll_lat, lon = ll_lon)
ur = CoordPair(lat = ur_lat, lon = ur_lon)
bounds = GeoBounds(ll, ur)


for archivo in glob.glob(PATH_SIHVIGILA + 'SIHpp*'):

    print(archivo)

    with Dataset(archivo) as infile, Dataset(archivo + '_R', mode = "w") as outfile:
        print ("reduce getting called!")

        lons = [bounds.bottom_left.lon, bounds.top_right.lon]
        lats = [bounds.bottom_left.lat, bounds.top_right.lat]
        orig_west_east = len(infile.dimensions["west_east"])
        orig_south_north = len(infile.dimensions["south_north"])

        x_y = ll_to_xy(infile, lats, lons, meta=False)
        start_x = 0 if x_y[0,0] == 0 else x_y[0,0] - 1
        end_x = orig_west_east - 1 if x_y[0,1] >= orig_west_east - 1 else x_y[0,1] + 1
        start_y = 0 if x_y[1,0] == 0 else x_y[1,0] - 1
        end_y = orig_south_north if x_y[1,1] >= orig_south_north - 1 else x_y[1,1] + 1

        west_east = end_x - start_x + 1
        west_east_stag = west_east + 1
        south_north = end_y - start_y + 1
        south_north_stag = south_north + 1


        # New dimension sizes
        dim_d = {"west_east" : west_east,
                 "west_east_stag" : west_east_stag,
                 "south_north" : south_north,
                 "south_north_stag" : south_north_stag
                }

        # Data slice sizes for the 2D dimensions
        slice_d = {"west_east" : slice(start_x, end_x + 1),
                   "west_east_stag" : slice(start_x, end_x + 2),
                   "south_north" : slice(start_y, end_y + 1),
                   "south_north_stag" : slice(start_y, end_y + 2)
                  }

        # Copy the global attributes
        atributos = deepcopy(infile.__dict__)
        atributos['WEST-EAST_GRID_DIMENSION'] = np.int32(west_east_stag)#self._west_east_stag
        atributos['SOUTH-NORTH_GRID_DIMENSION'] = np.int32(south_north_stag)#self._west_east_stag
        atributos['WEST-EAST_PATCH_END_UNSTAG'] = np.int32(west_east)#self._west_east_stag
        atributos['WEST-EAST_PATCH_END_STAG'] = np.int32(west_east_stag)#self._west_east_stag
        atributos['SOUTH-NORTH_PATCH_END_UNSTAG'] = np.int32(south_north)#self._west_east_stag
        atributos['SOUTH-NORTH_PATCH_END_STAG'] = np.int32(south_north_stag)#self._west_east_stag
        #atributos[''] = np.int32(dim_d['south_north_stag'])#self._west_east_stag
        outfile.setncatts(atributos)
        #outfile.setncatts(infile.__dict__)

        # Copy Dimensions, limiting south_north and west_east to desired domain
        for name, dimension in infile.dimensions.items():
            dimsize = dim_d.get(name, len(dimension))
            outfile.createDimension(name, dimsize)

        # Copy Variables  
        for name, variable in infile.variables.items():

            new_slices = tuple((slice_d.get(dimname, slice(None)) for dimname in variable.dimensions))

            outvar = outfile.createVariable(name, variable.datatype, variable.dimensions)

            outvar[:] = variable[new_slices]

            outvar.setncatts(variable.__dict__)


