#from in_data import InputRadar
import numpy as np
from collections import defaultdict
from utils.so_methods import so as cs

class SoOutput(object):
    """
    
    Parametros
    ----------
    radar : objeto
        Objeto radar original
    grid_dims : array
        Array con las dimensiones de la retícula donde se 
        hará el superobbing 
    variables : lista
        Lista de nombres de los campos que se van a procesar
    rays : lista
        Indice inicial y final del haz 
    var_opt : lista
        
    weight_vars : bool
        
    """
    
    def __init__(self, radar, grid_dims, variables, rays, var_opt, weight_vars):
        self.radar = SoRadar(radar, variables, rays)
        self.var_opt = var_opt
        self.grid = Grid(radar, *grid_dims)
        self.fields = defaultdict(dict)
        self.weight_vars = weight_vars
        self.compute_so(variables)
        
    def compute_so(self, variables):
        """
        Hacer el superobbing para las variables especificadas
        
        Parametros
        ----------
        variables : lista o str
            Nombre de variables del radar
        """
        vars_data = {}
        grid = self.grid
        radar = self.radar
        for var in variables:
            if var in list(radar.fields.keys()):
                # Asignarle las variables al atributo fields
                self.allocate_var(var, grid.nlon, grid.nlat, grid.nlev)
                # Obtener los datos del objeto radar
                vars_data[var] = radar.fields[var]
                # Convertir los dBZ en potencia
                if var == 'cref':
                    vars_data[var]['data'] = np.power(10., vars_data[var]['data']/10.)
                    vars_data[var]['data'].data[vars_data[var]['data'].mask] = vars_data[var]['_FillValue']
        self.compute_grid_boxmean(vars_data)
                
    def allocate_var(self, var, nx, ny, nz):
         """ Crea campos vacíos donde se van a guardar las variables luego del superobbing """
         self.fields['grid_' + var]['nobs'] = np.zeros((nz, ny, nx))
         self.fields['grid_' + var]['data'] = np.zeros((nz, ny, nx))
         self.fields['grid_' + var]['az'] = np.zeros((nz, ny, nx))
         self.fields['grid_' + var]['el'] = np.zeros((nz, ny, nx))
         self.fields['grid_' + var]['ra'] = np.zeros((nz, ny, nx))
         self.fields['grid_' + var]['id'] =  self.var_opt[0]
         self.fields['grid_' + var]['error'] =  self.var_opt[1]
            
         if len(self.var_opt) == 3:
             self.fields['grid_' + var]['min'] = self.var_opt[2]

                
    def compute_grid_boxmean(self, variables):
        """ Promedio en cajas para cada punto de reticula del superobbing """
        radar = self.radar
        nr = radar.ngates
        na = radar.nrays
        
        nvar = len(self.fields.keys())
        
        grid = self.grid
        lon_ini = grid.lon[0, 0]
        lat_ini = grid.lat[0, 0]
        z_ini = 0.0
        
        dlon = grid.dlon
        dlat = grid.dlat
        dz = grid.dz
        
        nlon = grid.nlon
        nlat = grid.nlat
        nz = grid.nlev
        
        local_undef = -8888.
        
        # Agrupa todas las variables que pasaran por el superobbing en un solo array
        # Esto ayuda a la parelizacion
        datain = np.zeros((na*nr, 4*nvar))
        latin = np.zeros((na*nr))
        lonin = np.zeros((na*nr))
        zin = np.zeros((na*nr))
        
        var_names = []
        weight = np.zeros(4*nvar).astype(bool)
        weight_ind = np.zeros(4*nvar)
        
        w_i = 0
        for iv, var in enumerate(variables):
            for ir in range(nr):
                datain[ir*na : (ir+1)*na, 4*iv+0] = radar.azimuth['data'][:]
                datain[ir*na : (ir+1)*na, 4*iv+1] = radar.elevation['data'][:]
                datain[ir*na : (ir+1)*na, 4*iv+2] = radar.range['data'][ir]
                datain[ir*na : (ir+1)*na, 4*iv+3] = variables[var]['data'].data[:, ir]

                latin[ir*na : (ir+1)*na] = radar.gate_latitude['data'][:, ir]
                lonin[ir*na : (ir+1)*na] = radar.gate_longitude['data'][:, ir]
                zin[ir*na : (ir+1)*na] = radar.gate_altitude['data'][:, ir]
                
            # Cambiar la variables undef a local_undef    
            tmp_mask = (datain[:, 4*iv+3] == variables[var]['_FillValue'])
            for ii in range(0, 4):
                datain[tmp_mask, 4*iv+ii] = local_undef
                
            if var == 'cref' and self.weight_vars:
                w_i = 4*iv + 3 # La reflectividad se puede usar como un peso
                weight = np.ones(4*nvar).astype(bool)
                weight[w_i-4 : w_i] = False # No pesar la reflectividad
                weight_ind = np.ones(4*nvar) * w_i
                
        is_angle = np.zeros(4*nvar) # este flag decide si la variable se trata como un angulo en el superobbing
        is_angle[range(0, 4*nvar, 4)] = True
        [data_ave, data_max, data_min, data_std, data_n, data_w] = cs.com_interp_boxavereg(
                             xini=lon_ini, dx=dlon, nx=nlon, yini=lat_ini, dy=dlat, ny=nlat,
                             zini=z_ini, dz=dz, nz=nz, nvar=4*nvar,
                             xin=lonin, yin=latin, zin=zin, datain=datain, undef=local_undef, nin=na*nr,   
                             weigth=weight, weigth_ind=weight_ind, is_angle=is_angle)
        #"Unpack" the superobbed data.
        for iv , var in enumerate(variables):
            self.fields['grid_' + var]['az'] = data_ave[:,:,:,0+iv*4]
            self.fields['grid_' + var]['el'] = data_ave[:,:,:,1+iv*4]
            self.fields['grid_' + var]['ra'] = data_ave[:,:,:,2+iv*4]
            self.fields['grid_' + var]['data'] = data_ave[:,:,:,3+iv*4]
            self.fields['grid_' + var]['nobs'] = data_n[:,:,:,3+iv*4]

class SoRadar(object):
    """
    A PyArt-Like radar object with necessary attributes and variables.

    Parameters
    ----------
    filename : str
        Radar file.
    variables: list
        Variables to extract from radar object.
    ray_interval : list
        Initial and final ray number to reduce radar fields data.

    Attributes
    ----------
    gate_longitude, gate_latitude : LazyLoadDict
        Geographic location of each gate.  The projection parameter(s) defined
        in the `projection` attribute are used to perform an inverse map
        projection from the Cartesian gate locations relative to the radar
        location to longitudes and latitudes.

    """
    def __init__(self, input, variables, ray_interval):
        self.__radar = input
        self.__variables = variables
        self.__ray_interval = ray_interval
        
        self.gate_longitude = {}
        self.gate_latitude = {}
        self.gate_altitude = {}
        self.azimuth = {}
        self.elevation = {}
        self.fields = defaultdict(dict)

        self.ngates = self.__radar.ngates
        self.nrays = ray_interval[-1] - ray_interval[0]
        self.longitude = self.__radar.longitude
        self.latitude = self.__radar.latitude
        self.altitude = self.__radar.altitude
        self.range = self.__radar.range

        self.get_data()

    def get_data(self):
        self.gate_longitude['data'] = self._extract_rays(self.__radar.gate_longitude['data']).copy()
        self.gate_latitude['data'] = self._extract_rays(self.__radar.gate_latitude['data']).copy()
        self.gate_altitude['data'] = self._extract_rays(self.__radar.gate_altitude['data']).copy()
        self.azimuth['data'] = self._extract_rays(self.__radar.azimuth['data']).copy()
        self.elevation['data'] = self._extract_rays(self.__radar.elevation['data']).copy()
        
        for var in self.__variables:
            self.fields[var] =  self.__radar.fields[var].copy()
            self.fields[var]['data'] = self._extract_rays(self.__radar.fields[var]['data']).copy()

    def _extract_rays(self, a):
        if a.ndim == 1:
            return a[self.__ray_interval[0]:self.__ray_interval[-1]]
        return a[self.__ray_interval[0]:self.__ray_interval[-1], :]

                
class Grid(object):
    """
    Reticula cartesiana.

    Parámetros
    ----------
    radar : Radar
        objeto radar de Pyart para usar como georeferencia.
    dx : int
        Dimensión horizontal de la reticula en metros.
    dz : int
        Dimensión vertical de la reticula en metros.
    maxz : int
        Altura máxima en metros.

    Atributos
    ---------
    dx : int
        Dimensión horizontal de la reticula en metros.
    dz : int
        Dimensión vertical de la reticula en metros.
    nlon, nlat, nlev : int
        Dimensiones del eje de longitud, latitud y altura.
    dlon, dlat : float
        Distnacia en grados entre los puntos de reticula en el eje de longitud y latitud, respectivamente.
    lon : numpy array
        Longitud de cada punto de reticula.
    lat : numpy array
        Latitud de cada punto de reticula.
    z : numpy array
        Altura de cada punto de reticula.
    """
    def __init__(self, radar, dx, dz, maxz , maxrange = None):

        self.dx = dx
        self.dz = dz

        # Compute possible value for `nlon` in order to cover the maximum radar range
        if maxrange == None  :
            maxrange = np.max(radar.range['data'])
            
        self.nlon = np.ceil(2.*maxrange/self.dx).astype('int')
        self.nlat = self.nlon
        self.nlev = np.ceil(maxz/self.dz).astype('int')

        # Force grid dimension to be odd
        if np.mod(self.nlon, 2) == 0:
            self.nlon += 1
        if np.mod(self.nlat, 2) == 0:
            self.nlat += 1

        self.dlon = None
        self.dlat = None
        self.lon = np.zeros((self.nlat, self.nlat))*np.nan
        self.lat = np.zeros((self.nlat, self.nlon))*np.nan
        self.z = np.zeros((self.nlev, self.nlat, self.nlon))*np.nan

        radar_lon = radar.longitude['data']
        radar_lat = radar.latitude['data']

        self.dx2ddeg(radar_lat)
        self.lat_lon(radar_lon, radar_lat)
        self.zlevs()

    def dx2ddeg(self, radar_lat):
        """
        Pasa un diferencial de distancia `dx` a un diferencial de latitud y longitud
        Hopefully, nobody will put a radar at the pole.
        """
        re = 6371e3
        self.dlon = float(np.rad2deg(self.dx/(re*np.cos(np.deg2rad(radar_lat)))))
        self.dlat = float(np.rad2deg(self.dx/re))

    def lat_lon(self, radar_lon, radar_lat):
        """ Define las latitudes y longitudes"""
        for i in range(self.nlat):
            for j in range(self.nlon):
                self.lon[i, j] = radar_lon + self.dlon*(-1.0-(self.nlon-1.0)/2.0 + (j+1))
                self.lat[i, j] = radar_lat + self.dlat*(-1.0-(self.nlat-1.0)/2.0 + (i+1))

    def zlevs(self):
        """ Definir la altura de los niveles """
        for k in range(self.nlev):
            self.z[k] = k*self.dz 
