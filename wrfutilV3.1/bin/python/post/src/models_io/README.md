# MODELS_IO

Herramienta de Python para postprocesar salidas de los pronósticos de WRF y GFS y generar un Dataset de [Xarray](https://docs.xarray.dev/en/stable/) listo para ser guardado como un NetCDF que cumple con las [Convenciones CF](https://cfconventions.org/) o seguir operando sobre los datos. 
Además del código se brinda un catálogo de variables que incluye, entre otras cosas, el nombre de la variable en los archivos de salida de los modelos y los atributos que va a tener la variable en el Dataset.

## Capacidades

* Leer una sola salida de pronóstico o varias de una misma inicialización y combinarlas en un solo dataset.
* Seleccionar variables específicas a leer del modelo.
* Seleccionar los datos en un nivel específico de presión en el caso de GFS e interpolar a niveles de presión, altura o temperatura o seleecionar un nivel sigma en particual para el caso de WRF.

## Limitaciones
* No se pueden seleccionar variables distintas con igual coordenada vertical pero en distintos niveles.
* No se puede obtener una misma variable en 2 coordenadas verticales distintas.

## Ejemplo de uso

### Lectura de un solo pronóstico del WRF

    # Se importa el módulo que contiene el código y el catálogo de variables
    import util_post as up
    import catalog_variables as catalog
    
    # Nombre del archivo a procesar
    filename = 'wrfout_d01_2022-01-01_00_00_00'

    # Diccionario con las variables y niveles que se quieren procesar. La "key" representa la variable que se quiere extraer (con el nombre que figura en el catálogo, no la kque figura en el modelo) y el valor la variable del catálogo que tiene los valores sobre los que se quiere interpolar. En caso de una variable 2D o que no se quiere extraer la coordenada vertical se pone None 
    variables = {'T': 'PRESSURE', 
                'Umet': 'Z_AGL', 
                'Vmet': 'Z_AGL',
                 'T2': None}

    # Diccionario con las variables que tienen los valores sobre los que se va a interpolar y los niveles que se quiere intepolar
    levels = {'PRESSURE': [1000, 850, 500, 250], 'Z_AGL': [500, 1000, 1500]}

    # En este caso se va a obtener los campos de temperatura interpolados a 1000, 850, 500, 250 hPa, los vientos zonal y meridional a 500, 1000 y 1500 m y la temperatura a 2m

    model_type = 'WRF'

    ds = util_post.smn_netcdf_from_file(filename, variables, catalog, levels, model_type)

    ds.to_netcdf('postprocessing.nc')

### Lectura de varios pronósticos del GEFS

    # Se importa el módulo que contiene el código y el catálogo de variables
    import util_post as up
    import catalog_variables as catalog
    
    import glob

    # Lista de archivos a postprocesar
    lista = glob.glob('gep*.t12z.pgrb2.0p50.f*')

    variables = {'T': 'PRESSURE', 
                 'Umet10': None, 
                 'Vmet10': None}

    levels = {'PRESSURE': [1000, 850, 500, 250]}

    model_type = 'GFS'

    ds = util_post.smn_netcdf_from_list(lista, variables, catalog, levels, model_type)


En este caso, mediante la librería [Dask](https://www.dask.org/), la lectura de los archivos se puede realizar en paralelo definiendo la variable de entorno ICORE con el número de procesos a utilizar (export ICORE=10 para utilizar 10 procesos por ejemplo). En caso de no definirla se utiliza un solo proceso.


## Agregado de variables al catálogo

Las variables en el catálogo están representadas de la siguiente manera

    T = {'name': 'T',
         'wrf': {'name': 'tk'},
         'gfs': {'name': 't',
                 'level': 'isobaricInhPa'},
         'dtype': 'float32',
         'attrs': {'units': 'K',
                   'standard_name': 'air_temperature',
                   'long_name': 'Temperature'}}

* name: representa el nombre asignado a la variable en el catálogo y debe ser igual al nombre de la variable que contiene el diccionario.
* wrf: diccionario en el que 'name' representa el nombre de la variable en la librería [WRF-python](https://wrf-python.readthedocs.io/en/latest/index.html) o el nombre de la variable en los wrfout<sup>1</sup>.
* gfs: diccionario en el que 'name' representa el atributo 'cfVarName' de la variable en los grib y level el nombre del nivel en que se encuentra la variable (atributo 'typeOfLevel' de la variable en el grib)<sup>1</sup>.
* dtype: tipo de dato al que se quiere convertir la variable. En algunos casos especiales (variables relacionadas a fechas por ejemplo) podría no ser necesariamente el tipo de dato que se quiere guardar.
* attrs: diccionario con los atributos obligatorios que deben tener las variable al guardarse en un NetCDF. Se pueden agregar más atributos en caso de necesitarse. 
    - 'units': unidades en que se quiere los datos.
    - 'standard_name': nombre estandar que tiene la variable segun las [Convenciones CF](https://cfconventions.org/). La lista de nombres se encuentra en https://cfconventions.org/Data/cf-standard-names/current/build/cf-standard-name-table.html
    - 'long_name': descripción de la variable.

<sup>1</sup> En caso de que la variable no se encuentre en algún modelo, se debe completar el diccionario con None.
 
El orden en que se encuentran las variables en el catálogo no es importante pero para mantener un orden se las ordenó según si son variables 2D como la precipitación, 3D como el geopotencial, estáticos en el tiempo como la topografía o variables que indican fechas y por orden alfabético. En lo posible, seguir este ordenamiento.
