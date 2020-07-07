# Descarga de datos clima usando GEE 
<img src="https://ciat.cgiar.org/wp-content/uploads/Alliance_logo.png" alt="CIAT" id="logo" data-height-percentage="90" data-actual-width="140" data-actual-height="55">

## 1. Propósito

El propósito de este código es poder realizar la consulta y descarga de datos meteorológicos medidos via sensoramiento remoto y disponibles en Google Earth Engine (GEE). Las ventajas de utilizar este repositorio es que permite aumentar el número de datos a descargar, el cual esta limitado en 5000 en GEE. La extracción de la información, se lleva a cabo mediante coordenadas latitud/longitud en sistema geográfico WGS84.

## 2. Pasos iniciales

* **[Registrarse](https://earthengine.google.com/)**: La consulta y descarga de los datos se realiza a través de Google Earth Engine (GEE), para ello es necesario crear una cuenta con dicha aplicación.
* **[Configuración GEE](https://developers.google.com/earth-engine/python_install-conda)**: el código fue implementado en python, por lo cual se sugiere instalar conda o miniconda, para crear un ambiente dedicado a GEE. Adicional GEE requiere autenticar las credenciales del usuario.
* **[Descarga del repositorio](https://github.com/anaguilarar/gee_NOAA.git)**

## 3. Misiones Disponibles

* [CFS](https://developers.google.com/earth-engine/datasets/catalog/NOAA_CFSV2_FOR6H)
* [GLDAS](https://developers.google.com/earth-engine/datasets/catalog/NASA_GLDAS_V021_NOAH_G025_T3H)
* [CHIRPS](https://developers.google.com/earth-engine/datasets/catalog/UCSB-CHG_CHIRPS_DAILY)
* [DAYMET](https://developers.google.com/earth-engine/datasets/catalog/NASA_ORNL_DAYMET_V3)
* [ERA5](https://developers.google.com/earth-engine/datasets/catalog/ECMWF_ERA5_DAILY)

## 4. Requisitos

  1. python 3.6, 3.7
  2. Google earth Engine API
  3. Un archivo CSV el cual disponga de las coordenadas para cada punto de interés.
  4. En cuanto a librerias, se recomiendan las siguientes:
  ```txt
  earthengine_api==0.1.211
numpy==1.18.5
wget==3.2
folium==0.11.0
geehydro==0.2.0
pandas==1.0.5
ee==0.4
geopandas==0.8.0
Shapely==1.7.0
  ```

 ## 5. Como usarlo

El siguiente ejemplo muestra como visualizar y descargar datos diarios de precipitación de la misión CHIRPS. Para más ejemplos por favor dirigirse a el notebook [Descargar clima de GEE](https://github.com/anaguilarar/gee_NOAA/blob/master/Descargar%20clima%20de%20GEE.ipynb).

Para poder realizar la consulta inicial de la información, es necesario especificar el tiempo, la región de interés y la misión. Para el tiempo se señala la fecha inicial y final; en cuanto a la región, se indica la ruta del archivo csv; por último, se indica cual de los siguientes nombres corresponde a la misión de interés: "cfs", "gldas", "chirsp", "daymet" y "era5"

```python
### importando el script

import get_geedata

## realizando la consulta
datachirps = get_geedata.gee_weatherdata("2017-01-01", # Fecha Inicial
                                         "2020-04-20", # Fecha Final
                                         "data/truestationlocation.csv", # directorio
                                         "chirps") # misión
                                         
```

Con la intención de evitar la restricción de la descarga de máximo 5000 atributos, los puntos de intereés se subdividen y se realiza la descarga de cada subconjunto, al final se une la base de datos en un dataframe.

```python

dfchirps = datachirps.CHIRPSdata_asdf()

### gráficar resultados
import matplotlib.pyplot as plt

ref_long = datachirps.features.longitude.loc[0]
plotdata = dfchirps[np.round(dfchirps.longitude, 3) == np.round(ref_long, 3)]

plt.figure(figsize=[12,5])
plt.plot(plotdata.date, plotdata['precipitation'].values)

plt.show()

```

Finalmente el dataframe se guarda como csv, señalando el directorio de destino y el nombre del archivo.

```python
### exportar datos
import os
outputdir = "../results/
dfchirps.to_csv(os.path.join(outputdir,"chirpsdata.csv"))
```


 