###################################################################
# 19/12/2019
# graficar los analisis de las variables: varList_sup y varList_lev 
# para todos los miembros. Graficos de 40 panales
# Para correr: source activate wrfutil
# export ANAFILE=$(ls /data/yanina/wrfutilV3.1/RUNs/HIST/POST/20191216_190000/prueba/*ANA* )  

# FALTA: poner logo, ver estetica, el export de anafiles y hacer 
#        dinamica la cantidad de miembros.

# JPG pesa mas que PNG
###################################################################
from time import process_time
#tic = process_time()
import os
import numpy as np
from netCDF4 import Dataset
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import cartopy.feature as cfeature
import cartopy.crs as ccrs

#### CMAPS Y LEVELS  #### 
cmap_T2m=plt.get_cmap('jet')   ;levels_T2m=np.arange(-10,40,2)
cmap_PP=plt.get_cmap('jet')    ;levels_PP=[1,5,10,15,20,25,30,35,40,45,50]
cmap_U200=plt.get_cmap('jet')  ;levels_U200=np.arange(-20,75,5)
cmap_T500=plt.get_cmap('jet')  ;levels_T500=np.arange(-50,12,2)
cmap_Q850=plt.get_cmap('jet')  ;levels_Q850=np.arange(5,27,2)

cmap_T2m.set_under('pink')     ;cmap_T2m.set_over('indigo') 
cmap_PP.set_over('indigo')     ;
cmap_U200.set_under('pink')    ;cmap_U200.set_over('indigo')
cmap_T500.set_under('pink')    ;cmap_T500.set_over('indigo')
cmap_Q850.set_over('indigo')   ;

figsize = (20,10)
fontsize_title=20
fontsize_cbar=15

## USO EL CARTOPY ##
states_provinces = cfeature.NaturalEarthFeature(category='cultural',name='admin_1_states_provinces_lines',scale='10m',facecolor='none')
countries = cfeature.NaturalEarthFeature(category='cultural', name='admin_0_countries',scale='50m',facecolor='none')
# Definimos un sistema de referencia en una proyeccin determinada
crs_latlon = ccrs.PlateCarree()

#### VARIABLES A GRAFICAR ###  (alias en wrf, nivel-si es varList_lev, short name, titulo de figura, aux para cambiar de unidades, unidades, cmap, levels)
# Temperatura en Kelvin, la paso a Centigrados
varList_sup=[("T2"           ,"T2m"    ,"Análisis temperatura a 2m"                 ,273.1   ,"[°C]"       ,cmap_T2m    ,levels_T2m  ,'both'),
             ("PP"           ,"PP"     ,"Análisis Precipitación acumulada"          ,0       ,"[mm]"      ,cmap_PP     ,levels_PP   ,'max')]

varList_lev=[( "Umet"  ,2    ,"u200"   ,"Análisis Viento zonal en 200 hPa"          ,0       ,"[m s-1]"   ,cmap_U200   ,levels_U200 ,'both'),
             ("T"      ,1    ,"T500"   ,"Análisis Temperatura en 500 hPa"           ,273.1   ,"[°C]"       ,cmap_T500   ,levels_T500 ,'both' ),
             ("Q"      ,0    ,"Q850"   ,"Análisis humedad epecifica en 850 hPa"     ,0       ,"[g Kg-1]"  ,cmap_Q850   ,levels_Q850 ,'max')]

#### LISTA DE ARCHIVOS Y LOS ABRE, LOS GUARDA EN UN DICCIONARIO ####
pathout=os.environ['PATHOUT_PLOT']
anafiles=os.environ['ANAFILES']
anafiles=anafiles.split("\n")
netfiles={} 
for afile in anafiles: ### Este for genera un diccionario con clave:Miembro, valor:(archivo abrierto,nombre archivo, fecha,miembro)
      miem=int(afile.split('.')[-2].split("OPER")[-1])
      netfiles[miem]=(Dataset(afile,'r'), afile, afile.split('.')[-4], miem)

#### VARIABLES EN SUPERFICIE ####
for tupla in varList_sup:
    (var_wrf,sname,lname,aux,units,cmap,levels,extend)=tupla

    fig=plt.figure(figsize = figsize)
    fig.subplots_adjust(top=0.92,right=0.87,wspace=0.1, hspace=0.1)
    lats = netfiles[list(netfiles.keys())[0]][0].variables["XLAT"][:] # Latitudes de las variables de masa
    lons = netfiles[list(netfiles.keys())[0]][0].variables["XLONG"][:] # Longitudes de las variables de masa
    for miem in range(1,41):
         ## que pasa si falta un miembro
         try:
                dato = netfiles[miem][0].variables[var_wrf][0,:,:]
         except :
                dato = np.zeros_like(lats) # se completa con ceros
	 ## ploteo
         norm = mpl.colors.BoundaryNorm(levels, ncolors=cmap.N)
         ax1=plt.subplot(5,8,miem,projection=crs_latlon)
         cs=ax1.contourf(lons, lats, dato-aux, cmap=cmap, levels=levels, extend=extend, norm=norm)
         plt.title('Miembro '+str(miem))
         ax1.add_feature(states_provinces, linewidth=0.5, edgecolor="black",zorder=10)
         ax1.coastlines('50m', linewidth=0.8,zorder=10)
         ax1.add_feature(countries,edgecolor="black", linewidth=1.0,zorder=10)

    cbar_ax = fig.add_axes([0.9, 0.15, 0.02, 0.75])
    fig.colorbar(cs, cax=cbar_ax,ticks=levels)
    cbar_ax.tick_params(labelsize=fontsize_cbar)
    plt.suptitle(lname + ' ' + units , fontsize=fontsize_title)
    plt.savefig(pathout +'/ANA_'+sname+'_'+netfiles[miem][2]+'.png', quality = 95, bbox_inches = 'tight')

#### VARIABLES EN NIVELES ####
for tupla in varList_lev:
    (var_wrf,lev,sname,lname,aux,units,cmap,levels,extend)=tupla

    fig=plt.figure(figsize = figsize)
    fig.subplots_adjust(top=0.92,right=0.87,wspace=0.1, hspace=0.1)
    lats = netfiles[list(netfiles.keys())[0]][0].variables["XLAT"][:] # Latitudes de las variables de masa
    lons = netfiles[list(netfiles.keys())[0]][0].variables["XLONG"][:] # Longitudes de las variables de masa
    for miem in range(1,41):
         try:
                dato = netfiles[miem][0].variables[var_wrf][0,lev,:,:]
         except :
                dato = np.zeros_like(lats)
         ## ploteo
         norm = mpl.colors.BoundaryNorm(levels, ncolors=cmap.N)
         ax1=plt.subplot(5,8,miem,projection=crs_latlon)
         cs=ax1.contourf(lons, lats,dato-aux, cmap = cmap, levels=levels, extend=extend, norm=norm)
         plt.title('Miembro '+str(miem))
         ax1.add_feature(states_provinces, linewidth=0.5, edgecolor="black",zorder=10)
         ax1.coastlines('50m', linewidth=0.8,zorder=10)
         ax1.add_feature(countries,edgecolor="black", linewidth=1.0,zorder=10)

    cbar_ax = fig.add_axes([0.9, 0.15, 0.02, 0.75])
    fig.colorbar(cs, cax=cbar_ax,ticks=levels)
    cbar_ax.tick_params(labelsize=fontsize_cbar)
    plt.suptitle(lname + ' ' + units , fontsize=fontsize_title)
    plt.savefig(pathout+'/ANA_'+sname+'_'+netfiles[miem][2]+'.png', quality = 95, bbox_inches = 'tight')

#### CIERRO LOS ARCHIVOS ####
for i in range(1,len(anafiles)): ### Este for genera un diccionario con clave:Miembro, valor:(archivo abrierto,nombre archivo, fecha,miembro)
      netfiles[i][0].close() 

#toc = process_time()
#print((toc -tic)/60.)
