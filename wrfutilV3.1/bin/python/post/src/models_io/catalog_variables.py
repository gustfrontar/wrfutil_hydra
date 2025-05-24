FILL_VALUE = 1e20

#####################
# DYNAMIC VARIABLES #
#####################
### 2D ###

CIN = {'name': 'CIN',
       'wrf': {'function': 'wrf.getvar', 'args': {'varname': 'cin2d_only'}},
       'gfs': {'name':None,
               'level':None},
       'dtype': 'float32',
       'attrs': {'units': 'J kg-1',
                 'standard_name': 'atmosphere_convective_inhibition',
                 'long_name': 'Maximum Convective Inhibition'}}

GRAUPELNC = {'name': 'GRAUPELNC',
             'wrf': {'function': 'wrf.getvar', 'args': {'varname': 'GRAUPELNC'}},
             'gfs': {'name':None,
                     'level':None},
             'dtype': 'float32',
             'attrs': {'units': 'mm',
                       'standard_name': 'graupel_fall_amount',
                       'long_name': 'Accumulated total graupel'}}


Gust10 = {'name': 'Gust10',
          'wrf': {'function': 'get_gust10', 'args': {}},
          'gfs': {'name':None,
                  'level':None},
          'dtype': 'float32',
          'attrs': {'units': 'm s-1',
                    'standard_name': 'wind_speed_of_gust',
                    'long_name': '10-m wind gust'}}

LCL = {'name': 'LCL',
       'wrf': {'function': 'wrf.getvar', 'args': {'varname': 'lcl'}},
       'gfs': {'name':None,
               'level':None},
       'dtype': 'float32',
       'attrs': {'units': 'm',
                 'standard_name': 'atmosphere_lifting_condensation_level',
                 'long_name': 'Lifted Condensation Level'}}

LFC = {'name': 'LFC',
       'wrf': {'function': 'wrf.getvar', 'args': {'varname': 'lfc'}},
       'gfs': {'name':None,
               'level':None},
       'dtype': 'float32',
       'attrs': {'units': 'm',
                 'standard_name': 'atmosphere_level_of_free_convection',
                 'long_name': 'Level of Free Convection'}}



MCAPE = {'name': 'MCAPE',
         'wrf': {'function': 'wrf.getvar', 'args': {'varname': 'cape2d_only'}},
         'gfs': {'name':None,
                 'level':None},
         'dtype': 'float32',
         'attrs': {'units': 'J kg-1',
                   'standard_name': 'atmosphere_convective_available_potential_energy',
                   'long_name': 'Maximum CAPE'}}


MDBZ = {'name': 'MDBZ',
        'wrf': {'function': 'wrf.getvar', 'args': {'varname': 'mdbz'}},
        'gfs': {'name': None, 
                'level': None},
        'dtype': 'float32',
        'attrs': {'units': 'dBZ',
                  'standard_name': 'equivalent_reflectivity_factor',
                  'long_name': 'Max. Reflectivity'}}


PP = {'name': 'PP', 
      'wrf': {'function': 'wrf.getvar', 'args': {'varname': 'RAINNC'}},
      'gfs': {'name': 'tp',
              'level': 'surface'},
      'dtype': 'float32',
      'attrs': {'units': 'mm', 
                'standard_name':'precipitation_amount',
                'long_name':'Accumulated Total Precipitation'}}


PP_calib = {'name': 'PP_calib', 
            'wrf': {'function': 'wrf.getvar', 'args': {'varname': None}},
            'gfs': {'name': None, 
                    'level': None},
            'dtype': 'float32',
            'attrs': {'units': 'mm', 
                      'standard_name':'precipitation_amount',
                      'long_name':'Calibrated Accumulated Total Precipitatioan'}}


PSFC = {'name': 'PSFC',
        'wrf': {'function': 'wrf.getvar', 'args': {'varname': 'PSFC'}},
        'gfs': {'name':'sp',
                'level': 'surface'},
        'dtype': 'float32',
        'attrs': {'units': 'hPa',
                  'standard_name': 'surface_air_pressure',
                  'long_name': 'Surface Pressure'}}

PW = {'name': 'PW',
      'wrf': {'function': 'wrf.getvar', 'args': {'varname': 'pw'}},
      'gfs': {'name':None,
              'level':None},
      'dtype': 'float32',
      'attrs': {'units': 'kg m-2',
                'standard_name': 'atmosphere_mass_content_of_water_vapor',
                'long_name': 'Precipitable Water'}}


Q2 = {'name': 'Q2',
      'wrf': {'function': 'wrf.getvar', 'args': {'varname': 'Q2'}},
      'gfs': {'name': None,
              'level': None},
      'dtype': 'float32',
      'attrs': {'units': 'g kg-1',
                'standard_name': 'specific_humidity',
                'long_name': '2-m Water Vapor Mixing Ratio'}}


REFL = {'name': 'REFL',
        'wrf': {'function': 'wrf.getvar', 'args': {'varname': 'dbz'}},
        'gfs': {'name': None,
                'level': None},
        'dtype': 'float32',
        'attrs': {'units': 'dBZ',
                  'standard_name': 'equivalent_reflectivity_factor',
                  'long_name': 'Reflectivity at X km agl'}}

RH2 = {'name': 'RH2',
      'wrf': {'function': 'wrf.getvar', 'args': {'varname': 'rh2'}},
      'gfs': {'name': None,
              'level': None},
      'dtype': 'float32',
      'attrs': {'units': '%',
                'standard_name': 'relative_humidity',
                'long_name': '2-m Relative Humidity'}}


SLP = {'name': 'SLP',
       'wrf': {'function': 'wrf.getvar', 'args': {'varname': 'slp'}},
       'gfs': {'name': 'prmsl',
               'level': 'meanSea'},
       'dtype': 'float32',
       'attrs': {'units': 'hPa',
                 'stdandard_name': 'air_pressure_at_mean_sea_level',
                 'long_name': 'Sea Level Pressure'}}

SRH_1000 = {'name': 'SRH_1000', 
            'wrf': {'function': 'wrf.getvar', 'args': {'varname': 'helicity', 'top':1000}},
            'gfs': {'name': None,
                    'level': None},
            'dtype': 'float32',
            'attrs': {'units': 'm2 s-2',
                      'stdandard_name': 'storm_relative_helicity', #este standard_name no esta, lo invento
                      'long_name': 'Storm relative helicity at 1 km agl'}}

SRH_3000 = {'name': 'SRH_3000', 
            'wrf': {'function': 'wrf.getvar', 'args': {'varname': 'helicity', 'top':3000}},
            'gfs': {'name': None,
                    'level': None},
            'dtype': 'float32',
            'attrs': {'units': 'm2 s-2',
                      'stdandard_name': 'storm_relative_helicity', #este standard_name no esta, lo invento
                      'long_name': 'Storm relative helicity at 3 km agl'}}

UH = {'name': 'UH',
            'wrf': {'function': 'wrf.getvar', 'args': {'varname': 'updraft_helicity','bottom':2000,'top':5000}},
            'gfs': {'name': None,
                    'level': None},
            'dtype': 'float32',
            'attrs': {'units': 'm2 s-2',
                      'stdandard_name': 'updraft helicity', #este standard_name no esta, lo invento
                      'long_name': 'Updraft Helicity in the layer 2-5 km'}}


T2 = {'name': 'T2',
      'wrf': {'function': 'wrf.getvar', 'args': {'varname': 'T2'}},
      'gfs': {'name':'t2m',
              'level': 'heightAboveGround'},
      'dtype': 'float32',
      'attrs': {'units': 'K',
                'standard_name': 'air_temperature',
                'long_name': '2-m Temperature'}}

TD2 = {'name': 'TD2',
       'wrf': {'function': 'wrf.getvar', 'args': {'varname': 'td2'}},
       'gfs': {'name':'td2m',
               'level': 'heightAboveGround'},
       'dtype': 'float32',
       'attrs': {'units': 'K',
                 'standard_name': 'dew_point_temperature',
                 'long_name': '2-m Dew Point Temperature'}}


Tmin = {'name': 'Tmin',
        'wrf': {None},
        'gfs': {'name':'tmin',
                'level': 'heightAboveGround'}, 
        'dtype': 'float32',
        'attrs': {'units': 'K', 
                  'standard_name':'air_temperature',
                  'long_name':'Minimum Temperature'}}


Tmax = {'name': 'Tmax',
        'wrf': {'name': None},
        'gfs': {'name':'tmax',
                'level': 'heightAboveGround'}, 
        'dtype': 'float32',
        'attrs': {'units': 'K', 
                  'standard_name':'air_temperature',
                  'long_name':'Maximum Temperature'}}


Tmax_cal = {'name': 'Tmax',
            'wrf': {'name': None},
            'gfs': {'name': None,
                    'level': None},
            'dtype': 'float32',
            'attrs': {'units': 'K', 
                      'standard_name':'air_temperature',
                      'long_name':'Calibrated Maximum Temperature'}}


Tmin_cal = {'name': 'Tmin',
            'wrf': {'name': None},
            'gfs': {'name': None,
                    'level': None},
            'dtype': 'float32',
            'attrs': {'units': 'K', 
                      'standard_name':'air_temperature',
                      'long_name':'Calibrated Minimum Temperature'}}

TSK = {'name': 'TSK',
       'wrf': {'function': 'wrf.getvar', 'args': {'varname': 'TSK'}},
       'gfs': {'name':None,
               'level':None},
       'dtype': 'float32',
       'attrs': {'units': 'K',
                 'standard_name': 'surface_temperature',
                 'long_name': 'Skin Temperature'}}

PBLH = {'name': 'PBLH',
        'wrf': {'function': 'wrf.getvar', 'args': {'varname': 'PBLH'}},
        'gfs': {'name':None,
               'level':None},
        'dtype': 'float32',
        'attrs': {'units': 'm',
                 'standard_name': 'atmosphere_boundary_layer_thickness',
                 'long_name': 'Planetary Boundary Layer Heigth'}}


Umet10 = {'name': 'Umet10',
          'wrf': {'function': 'wrf.getvar', 'args': {'varname': 'uvmet10'}},
          'gfs': {'name':'u10',
                  'level': 'heightAboveGround'},
          'dtype': 'float32',
          'attrs': {'units': 'm s-1',
                    'standard_name': 'eastward_wind',
                    'long_name': '10-m U Wind Component'}}


Vmet10 = {'name': 'Vmet10',
          'wrf': {'function': 'wrf.getvar', 'args': {'varname': 'uvmet10'}},
          'gfs': {'name':'v10',
                  'level': 'heightAboveGround'},
          'dtype': 'float32',
          'attrs': {'units': 'm s-1',
                    'standard_name': 'northward_wind',
                    'long_name': '10-m V Wind Component'}}

#####################
# DYNAMIC VARIABLES #
#####################
### 3D ###

DBZ = {'name': 'DBZ',
       'wrf': {'function': 'wrf.getvar', 'args': {'varname': 'dbz'}},
       'gfs': {'name': None,
               'level': None},
       'dtype': 'float32',
       'attrs': {'units': 'dBZ',
                 'standard_name': 'equivalent_reflectivity_factor',
                 'long_name': 'Reflectivity'}}


GEOPT = {'name': 'GEOPT',
         'wrf': {'function': 'wrf.getvar', 'args': {'varname': 'geopt'}},
         'gfs': {'name': 'gh',
                 'level': 'isobaricInhPa'},
         'dtype': 'float32',
         'attrs': {'units': 'm2 s-2',
                   'standard_name': 'geopotential_height',
                   'long_name': 'Geopotential Height'}}

PRES = {'name': 'PRES',
        'wrf': {'function': 'wrf.getvar', 'args': {'varname': 'pres'}},
        'gfs': {'name': None,
                'level': None},
        'dtype': 'float32',
        'attrs': {'units': 'Pa',
                  'standard_name': 'air_pressure',
                  'long_name': 'Air Pressure'}}

PRESSURE = {'name': 'PRESSURE',
            'wrf': {'function': 'wrf.getvar', 'args': {'varname': 'pressure'}},
            'gfs': {'name': None,
                    'level': None},
            'dtype': 'float32',
            'attrs': {'units': 'hPa',
                      'standard_name': 'air_pressure',
                      'long_name': 'Air Pressure'}}

Q = {'name': 'Q',
     'wrf': {'function': 'wrf.getvar', 'args': {'varname': 'QVAPOR'}},
     'gfs': {'name': 'q',
             'level': 'isobaricInhPa'},
     'dtype': 'float32',
     'attrs': {'units': 'g kg-1',
               'standard_name': 'specific_humidity',
               'long_name': 'Specific Humidity'}}


QCLOUD = {'name': 'QCLOUD',
          'wrf': {'function': 'wrf.getvar', 'args': {'varname': 'QCLOUD'}},
          'gfs': {'name': None,
                  'level': None},
          'dtype': 'float32',
          'attrs': {'units': 'g kg-1',
                    'standard_name': 'cloud_water_mixing_ratio',
                    'long_name': 'Cloud Water Mixing Ratio'}}


QGRAUP = {'name': 'QGRAUP',
          'wrf': {'function': 'wrf.getvar', 'args': {'varname': 'QGRAUP'}},
          'gfs': {'name': None,
                  'level': None},
          'dtype': 'float32',
          'attrs': {'units': 'g kg-1',
                    'standard_name': 'graupel_mixing_ratio',
                    'long_name': 'Graupel Mixing Ratio'}}


QICE = {'name': 'QICE',
        'wrf': {'function': 'wrf.getvar', 'args': {'varname': 'QICE'}},
        'gfs': {'name': None,
                'level': None},
        'dtype': 'float32',
        'attrs': {'units': 'g kg-1',
                  'standard_name': 'ice_mixing_ratio',
                  'long_name': 'Ice Mixing Ratio'}}


QRAIN = {'name': 'QRAIN',
         'wrf': {'function': 'wrf.getvar', 'args': {'varname': 'QRAIN'}},
         'gfs': {'name': None,
                 'level': None},
         'dtype': 'float32',
         'attrs': {'units': 'g kg-1',
                   'standard_name': 'rain_water_mixing_ratio',
                   'long_name': 'Rain Water Mixing Ratio'}}


QSNOW = {'name': 'QSNOW',
         'wrf': {'function': 'wrf.getvar', 'args': {'varname': 'QSNOW'}},
         'gfs': {'name': None,
                 'level': None},
         'dtype': 'float32',
         'attrs': {'units': 'g kg-1',
                   'standard_name': 'snow_mixing_ratio',
                   'long_name': 'Snow Mixing Ratio'}}

RH = {'name': 'RH',
      'wrf': {'function': 'wrf.getvar', 'args': {'varname': 'rh'}},
      'gfs': {'name': None,
              'level': None},
      'dtype': 'float32',
      'attrs': {'units': '%',
                'standard_name': 'relative_humidity',
                'long_name': 'Relative Humidity'}}

SMOIS = {'name': 'SMOIS',
         'wrf': {'function': 'wrf.getvar', 'args': {'varname': 'SMOIS'}},
         'gfs': {'name': None,
                 'level': None},
         'name_wrf': 'SMOIS',
         'dtype': 'float32',
         'attrs': {'units': 'm3 m-3',### como wrfout sale con m3/m3 lo pongo asi, pero segun CF deberia ir kg/m2
                   'standard_name': 'mass_content_of_water_in_soil_layer', 
                   'long_name': 'SOIL MOISTURE'}}

T = {'name': 'T',
     'wrf': {'function': 'wrf.getvar', 'args': {'varname': 'tk'}},
     'gfs': {'name': 't',
             'level': 'isobaricInhPa'},
     'dtype': 'float32',
     'attrs': {'units': 'K',
               'standard_name': 'air_temperature',
               'long_name': 'Temperature'}}

TD = {'name': 'TD',
      'wrf': {'function': 'wrf.getvar', 'args': {'varname': 'td'}},
      'gfs': {'name': 'td',
              'level': 'isobaricInhPa'},
      'dtype': 'float32',
      'attrs': {'units': 'K',
                'standard_name': 'dew_point_temperature',
                'long_name': 'Dew Point Temperature'}}

TPE = {'name': 'TPE',
       'wrf': {'function': 'wrf.getvar', 'args': {'varname': 'eth'}},
       'gfs': {'name': None,
             'level': None},
       'dtype': 'float32',
       'attrs': {'units': 'K',
               'standard_name': 'air_equivalent_potential_temperature',
               'long_name': 'Equivalent Potential Temperature'}}


TSLB = {'name': 'TSLB',
        'wrf': {'function': 'wrf.getvar', 'args': {'varname': 'TSLB'}},
        'gfs': {'name': None,
                'level': None},
        'dtype': 'float32',
        'attrs': {'units': 'K',
                  'standard_name': 'soil_temperature',
                  'long_name': 'Soil Temperature'}}

Umet = {'name': 'Umet',
        'wrf': {'function': 'wrf.getvar', 'args': {'varname': 'uvmet'}},
        'gfs': {'name': 'u',
                'level': 'isobaricInhPa'},
        'dtype': 'float32',
        'attrs': {'units': 'm s-1',
                  'standard_name': 'eastward_wind',
                  'long_name': 'Zonal Wind Component'}}

Vmet = {'name': 'Vmet',
        'wrf': {'function': 'wrf.getvar', 'args': {'varname': 'uvmet'}},
        'gfs': {'name': 'v',
                'level': 'isobaricInhPa'},
        'dtype': 'float32',
        'attrs': {'units': 'm s-1',
                  'standard_name': 'northward_wind',
                  'long_name': 'Meridional Wind Component'}}

W = {'name': 'W',
     'wrf': {'function': 'wrf.getvar', 'args': {'varname': 'wa'}},
     'gfs': {'name': None,
             'level': None},
     'dtype': 'float32',
     'attrs': {'units': 'm s-1', 
               'standard_name': 'upward_air_velocity', 
               'long_name': 'Z Wind Component'}}

Z_AGL = {'name': 'Z_AGL',
         'wrf': {'function': 'wrf.getvar', 'args': {'varname': 'height_agl'}},
         'gfs': {'name': None,
                 'level': None},
         'dtype': 'float32',
         'attrs': {'units': 'm',
                   'standard_name': 'height',
                   'long_name': 'Height above ground level'}}


####################
# STATIC VARIABLES #
####################

LANDMASK = {'name': 'LANDMASK',
            'wrf': {'function': 'wrf.getvar', 'args': {'varname': 'LANDMASK'}},
            'gfs': {'name':'lsm',
                    'level':'surface'},
            'dtype': 'int32',
            'attrs': {'units': '',
                      'standard_name': 'land_binary_mask',
                      'long_name': 'landmask'}}

HGT = {'name': 'HGT',
       'wrf': {'function': 'wrf.getvar', 'args': {'varname': 'HGT'}},
       'gfs': {'name': 'orog',
               'level': 'surface'},
       'dtype': 'float32',
       'attrs': {'units': 'm',
                 'standard_name': 'height_above_mean_sea_level',
                 'long_name': 'Terrain Height'}}


level_eta = {'name': 'level_eta',
             'wrf': {'function': 'wrf.getvar', 'args': {'varname': 'ZNU'}},
             'gfs': {'name': None,
                     'level': None},
             'dtype': 'float32',
             'attrs': {'units': '',
                       'standard_name': 'atmosphere_sigma_coordinate',
                       'long_name': 'Eta values on half (mass) levels',
                       'positive': 'down',
                       'axis': 'Z'}}


level_p = {'name': 'level_p',
           'wrf': {'name': None},
           'gfs': {'name': None,
                   'level': None},
           'dtype': 'float32',
           'attrs': {'units': 'hPa',
                     'standard_name': 'air_pressure',
                     'long_name': 'Pressure',
                     'positive': 'down',
                     'axis': 'Z'}}

level_t = {'name': 'level_t',
           'wrf': {'name': None},
           'gfs': {'name': None,
                   'level': None},
           'dtype': 'float32',
           'attrs': {'units': 'K',
                     'standard_name': 'air_temperature',
                     'long_name': 'Temperature',
                     'positive': 'down',
                     'axis': 'Z'}}


level_z = {'name': 'level_z',
           'wrf': {'name': None},
           'gfs': {'name': None,
                   'level': None},
           'dtype': 'float32',
           'attrs': {'units': 'm',
                     'standard_name': 'height',
                     'long_name': 'Height above ground level',
                     'positive': 'up',
                     'axis': 'Z'}}


MEMBER = {'name': 'MEMBER',
          'wrf': {'name': None},
          'gfs': {'name': 'number',
                  'level': None},
          'dtype': 'int32',
          'attrs': {'standard_name':'realization', 
                    'units':'', 
                    'long_name':'Ensemble member'}}


STEP = {'name': 'STEP',
        'wrf': {'name': None},
        'gfs': {'name':'step',
                'level': None},
        'dtype': 'int32',
        'attrs': {'standard_name': 'forecast_period',
                  'units': '', #Deberia completarse segun lo que se este haciendo
                  'long_name': 'Plazo de pronostico'}}


XLAT = {'name': 'XLAT',
        'wrf': {'function': 'wrf.getvar', 'args': {'varname': 'XLAT'}},
        'gfs': {'name': 'latitude',
                'level': None},
        'dtype': 'float32',
        'attrs': {'units': 'degrees_north',
                  'standard_name': 'latitude',
                  'long_name': 'Latitude'}}


XLONG = {'name': 'XLONG',
         'wrf': {'function': 'wrf.getvar', 'args': {'varname': 'XLONG'}},
         'gfs': {'name': 'longitude',
                 'level': None},
         'dtype': 'float32',
         'attrs': {'units': 'degrees_east',
                   'standard_name': 'longitude',
                   'long_name': 'Longitude'}}


x = {'name': 'x',
     'wrf': {'name': None},
     'gfs': {'name': None,
             'level': None},
     'dtype': 'float32',
     'attrs': {'units': 'm',
               'standard_name': 'projection_x_coordinate',
               'long_name': 'x-coordinate in projected coordinate system',
               'axis': 'X'}}


y = {'name': 'y',
     'wrf': {'name': None},
     'gfs': {'name': None,
             'level': None},
     'dtype': 'float32',
     'attrs': {'units': 'm',
               'standard_name': 'projection_y_coordinate',
               'long_name': 'y-coordinate in projected coordinate system',
               'axis': 'Y'}}


##################
# Time variables #
##################

#El atributo units va fuera del diccionario de atributos porque lo paso en el encoding

XTIME = {'name': 'XTIME',
         'wrf': {'function': 'wrf.getvar', 'args': {'varname': 'XTIME'}},
         'gfs': {'name': 'valid_time',
                 'level': None},
         'dtype': 'float64',
         'units': 'hours since',
         'attrs': {'standard_name': 'time',
                   'long_name': 'hours since',
                   'axis': 'T'}}


DATE = {'name': 'DATE',
        'wrf': {'name': None},
        'gfs': {'name': 'time',
                'level': None},
         'dtype': 'int32',
         'units': 'days since',
         'attrs': {'standard_name': 'forecast_reference_time', 
                   'long_name': 'Forecast initialization time',
                   'axis': 'T'}}

