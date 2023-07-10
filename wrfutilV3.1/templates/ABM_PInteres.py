import pandas as pd
import numpy as np
import os 
import datetime
import re

#### 
# Seccion de Configuracion
####
PRECISION=2
wrfutildir=os.environ["WRFUTILDIR"]
TABLAS={'tipos':{'file':"{}/templates/Tipo_{}.{}".format(wrfutildir,"{:%Y%m%d_%H%M%S}","{}"),
                 'link':"{}/templates/Tipo.dat".format(wrfutildir),
                 'columns':['tipo', 'descTipo', 'descID'],
                 'index':['tipo'],
                 'foreignkey':[]},
        'pinteres':{'file':"{}/templates/PInteres_{}.{}".format(wrfutildir,"{:%Y%m%d_%H%M%S}","{}"),
                    'link':"{}/templates/PInteres.dat".format(wrfutildir),
                    'columns':['tipo','ID','CIPI','nombre', 'desc'],
                    'index':['tipo','ID'],
                    'foreignkey':['tipo.tipos','CIPI.CIPI']},
        'CIPI':{'file':"{}/templates/CIPI_{}.{}".format(wrfutildir,"{:%Y%m%d_%H%M%S}","{}"),
                'link':"{}/templates/CIPI.dat".format(wrfutildir),
                'columns':['CIPI', 'lat', 'lon'],
                'index':['CIPI'],
                'foreignkey':[]}        
        }
####
# Funciones Auxiliares
####
def distancia(latA,lonA,latB,lonB):
    # distancia de dos corredenadas en metros
    d = np.sin(np.deg2rad(latA)) * np.sin(np.deg2rad(latB)) + (
      np.cos(np.deg2rad(latA)) * np.cos(np.deg2rad(latB)) * np.cos(
      np.deg2rad(lonA - lonB)))

    d[d>1] = 1 #a veces da mayor a 1 y el arccos de abajo rompe
    return np.arccos(d) * 6371000
    
def pandasRead(file,ext):
    if ext=="csv":
        res=pd.read_csv(file,encoding='latin-1')
        
    elif ext=="hf5":
        res=pd.read_hdf(file,encoding='latin-1')
    else :
        raise NameError('La extencion {} no esta soportada'.format(ext))
    return res

def pandasWrite(df,file,ext,link=None):
    now = datetime.datetime.now()
    actual=file.format(now,ext)
    if ext=="csv":
        df.to_csv(actual,index=True,encoding='latin-1')
    elif ext=="hf5":
        key="version_{:%Y%m%d_%H%M%S}".format(now)
        df.to_hdf(actual,key=key,append=True,format='table',index=True,encoding='latin-1')
    else :
        raise NameError('La extencion {} no esta soportada'.format(ext))
    if link is not None:
        os.symlink(actual,link+"tmp")
        os.rename(link+"tmp", link)

def inicializar(tablas,ext):
    for tabla in tablas.keys():
        if os.path.isfile(tablas[tabla]['link']):
            df=pandasRead(tablas[tabla]['link'],ext)
            if ext=="csv":
            	df.set_index(tablas[tabla]['index'],inplace=True,verify_integrity=True)
        else :
            df=pd.DataFrame(columns = tablas[tabla]['columns'])
            df.set_index(tablas[tabla]['index'],inplace=True,verify_integrity=True)
            pandasWrite(df,tablas[tabla]['file'],ext,tablas[tabla]['link'])
        tablas[tabla]['df']=df
    return tablas

def _alta_TIPO(tablas,tipo, desc, desc_id):
    clave=re.sub(r"[^A-Za-z0-9]+", '', tipo.upper())
    if not clave==tipo:
        print("WARNING: el tipo '{}' solicitado no cumple con el formato, usaremos '{}' ".format(tipo, clave))
    row=pd.DataFrame([[clave,desc,desc_id]],columns = tablas['tipos']['columns'])
    row.set_index(tablas['tipos']['index'],inplace=True)
    try:
        tablas['tipos']['df']=tablas['tipos']['df'].append(row,verify_integrity=True,ignore_index=False)
    except:
        return False
    return True

def _alta_CIPI(tablas, latitud, longitud):
    #Redondeamos las coordenadas para emparejar puntos cercanos
    lat=round(float(latitud),PRECISION)
    lon=round(float(longitud),PRECISION)
    #Verificamos si la coordenada ya tiene un CIPI asignado
    CIPIS=tablas['CIPI']['df'].loc[(tablas['CIPI']['df']['lat']==lat) & 
                                   (tablas['CIPI']['df']['lon']==lon)].index.values
    if len(CIPIS)==0:
        #LA coordenada es nueva, y hay que aregarla
        try :
            #calculamos el nuevo CIPI 
            CIPI=__nextCIPI(tablas['CIPI']['df'].index.values.max())
        except ValueError:
            # En el caso que la tabla esta vacia
            CIPI="AA00"
        #agregamos
        row=pd.DataFrame([[CIPI,lat,lon]],columns = tablas['CIPI']['columns'])
        row.set_index(tablas['CIPI']['index'],inplace=True)
        try:    
            tablas['CIPI']['df']=tablas['CIPI']['df'].append(row,verify_integrity=True,ignore_index=False)
        except :
            #No se pudo agregar por razones desconocidas
            raise 
    elif len(CIPIS)==1:
        # La coordenada ya existia y devolvemos el CIPI 
        CIPI=CIPIS[0]
    else :
        raise IndexError("La tabla CIPI tiene puntos duplicadas")
    return CIPI

def __nextCIPI(cipi):
    #verificamos la integridad del formato
    pattern = re.compile("^([A-Z]{2}[0-9]{2})$")
    try:
        codigo=re.fullmatch(pattern,cipi).groups()[0]
    except :
        raise ValueError("ERROR: Formato de CIPI invalido")
    res=''
    carry=1
    for char in codigo[::-1]:
        if char=="Z" or char=="9":
            res=("A" if char=="Z" else "0")+res
            carry=1
        else :  
            res=chr(ord(char) + carry)+res
            carry=0 
    return res

def __verificacionTIPO(tablas,tipo):
    clave=re.sub(r"[^A-Za-z0-9]+", '', tipo.upper())
    if not clave==tipo:
        print("WARNING: el tipo '{}' solicitado no cumple con el formato, usaremos '{}' ".format(tipo, clave))
    ## verificamos que el tipo exista
    if not (tablas['tipos']['df'].index.values==clave).any():
        print("ERROR: el tipo {} no existe, debe crearlo primero!".format(clave))
        return ''    
    return clave
def _alta_PINTERES(tablas,tipo,ID,nombre,desc,latitud,longitud):
    #Como tipo es clave, nos aseguramos que siempe este en el mismo formato
    clave=__verificacionTIPO(tablas,tipo)
    if clave=='':
        return False
    ## El tipo ya esta declarado, BIEN !
    ## Ahora vamos por el CIPI.
    CIPI=_alta_CIPI(tablas, latitud, longitud)       
    if (tablas['tipos']['df'].loc[[clave]]['descID'].values=='CIPI').any():
        if ID=='':
            ID=CIPI
        else : 
            raise ValueError("El ID es distinto de vacio y la descripcion del tipo indica 'CIPI'")
    if ID=='':
        raise ValueError("El ID es vacio pero la descripcion del tipo indicada no es 'CIPI'")
    
    #Ahora agegamos el pinteres
    row=pd.DataFrame([[clave,ID,CIPI,nombre,desc]],columns = tablas['pinteres']['columns'])
    row.set_index(tablas['pinteres']['index'],inplace=True)
    try:    
        tablas['pinteres']['df']=tablas['pinteres']['df'].append(row,verify_integrity=True,ignore_index=False)
    except ValueError:
        #No se pudo agregar, seguramente porque ya existia
        print("ERROR: insertando registro con indice ({},{}), probablemente ya existe".format(clave,ID))
        return False
    return True

def altaPInteresFile(tablas,file,fs=';',encabezado=False):
    res=[]
    if not os.path.isfile(file):
        print("ERROR: El file no existe ! intente nuevamente")
    with open(file,"r") as data:
        cont=0
        for pinteres in data:
            tipo, ID, nombre, desc, latitud , longitud,altitud =pinteres.split(fs)
            if encabezado and cont==0:
                cont+=1
                continue
            if not _alta_PINTERES(tablas, tipo, ID, nombre, desc, latitud, longitud):
                res.append(cont)
            cont+=1
    if not res==[]:
        print("los siguiente pinteres no pudieron agregarse:")
        print(res)


def altaPInteresManual(tablas):
    tipo=input("Ingrese el tipo del Pinteres que desea crear:")
    clave=__verificacionTIPO(tablas,tipo)
    if clave=='':
        return False
    ## El tipo ya esta declarado, BIEN !
    ID=''
    if not (tablas['tipos']['df'].loc[[clave]]['descID'].values=='CIPI').any():
        ID=input("Ingrese el ID correspondiente para el tipo={}\n".format(clave))
    latitud=input("Ingrese la latitud de PInteres a crear con 2 decimales de precision: ")
    longitud=input("Ingrese la longitud de PInteres a crear con 2 decimales de precision: ")
    nombre=input("Ingrese el nombre que se asignara al PInteres:")
    desc=input("De una descripcion de una linea del Pinteres:\n")
    return _alta_PINTERES(tablas,clave,ID,nombre,desc,latitud,longitud)
    
def altaTipoManual(tablas):
    ## Solicitamos el tipo
    tipo= input("Ingrese en indentificador de tipo (letras y numeros solo mayusculas):\n")
    clave=re.sub(r"[^A-Za-z0-9]+", '', tipo.upper())
    ## Nos ponemos de acuerdo en el formato del tipo
    if not clave==tipo:
        print("El tipo '{}' solicitado no cumple con el formato, usaremos '{}' ".format(tipo, clave))
        if not input("Esta deacuerdo (S/N):").upper()=="S":
            return False
    ## verificamos si el tipo existe
    if (tablas['tipos']['df'].index.values==clave).any():
        print("El tipo {} ya existe:".format(clave))
        print(tablas['tipos']['df'].loc[clave])
        return False   
    desc=input("Describa el tipo: ")
    descID=input("Describa el identificador de tipo, usar 'CIPI' si el tipo no tiene un identificador asociado: \n")
    if not _alta_TIPO(tablas, clave, desc, descID):
        if input("Hubo un error desconocido, desea intentarlo nuevamente (S/N)").upper()=="S":
            return altaTipoManual(tablas)
        else:
            return False
    return True

def bajaTIPO(tablas,tipo):
    clave=__verificacionTIPO(tablas,tipo)
    tablas['pinteres']['df'].drop(clave, level='tipo',inplace=True)
    tablas['tipos']['df'].drop(clave,inplace=True)

def bajaPInteres(tablas, tipo,ID):
    clave=__verificacionTIPO(tablas,tipo)
    tablas['pinteres']['df'].drop((clave,ID),inplace=True)
    
def saveAll(tablas,ext,link=True):
    for tabla in tablas.keys():
        pandasWrite(tablas[tabla]['df'], tablas[tabla]['file'], ext,link=tablas[tabla]['link'] if link else None)
    
if __name__ == "__main__" :
    try: 
        ext="csv"
        tablas=inicializar(TABLAS,ext)
    except:
        ext="hf5"
        tablas=inicializar(TABLAS,ext)
    menu="""
    1) Baja de Tipo       -> necesita el tipo
    2) Baja de PInteres   -> necesita el tipo y ID del tipo
    3) Alta de Tipo       -> necesita el tipo, descripcion y descripcion de tipo
    4) Alta de Pinteres   -> necesita el tipo, ID, nombre, descripcion, latitud y longitud
    5) Alta de Pinteres   -> necesita un archivo CSV con formato "tipo, ID, nombre, descripcion, latitud y longitud"
    6) Guardar y Salir    -> Guarda CSV y sale
    7) Salir SIN guardar
  
    """
    salir=False
    while not salir:
        print(menu)
        res=input("Elija Opcion: ")
        if res=="1" :
            tipo=input("Ingrese el tipo que desea eliminar: ")
            bajaTIPO(tablas, tipo)
        elif res=="2":
            tipo=input("Ingrese el tipo que desea eliminar: ")
            ID=input("Ingrese el ID correspondiente al tipo del PInteres: ")
            bajaPInteres(tablas, tipo, ID)
        elif res=="3":
            altaTipoManual(tablas)
        elif res=="4":
            altaPInteresManual(tablas)
        elif res=="5":
            file=input("ingrese el archivo: ")
            encabezado=input("El archivo tiene un encabezado? ")
            fs=input("Ingrese el separador de campo:")
            altaPInteresFile(tablas, file,fs,encabezado in ["s","S","SI","Si","si","y","yes","Y","YES"])
        elif res=="6":
            ext="csv"
            saveAll(tablas,ext)
            salir=True
        elif res=="7":
            salir=True
        else: 
            print("ERROR: Opcion no disponible!!! vuelva a intentarlo")
            
    
    
    
    
