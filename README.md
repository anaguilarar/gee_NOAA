# Descarga de datos NOAA 

* El propósito de este código es poder realizar la consulta y descarga de datos NOAA (https://developers.google.com/earth-engine/datasets/catalog/NOAA_CFSV2_FOR6H) y CHIRPS (https://developers.google.com/earth-engine/datasets/catalog/UCSB-CHG_CHIRPS_DAILY). La extracción de los datos se lleva a cabo mediante coordenadas latitud/longitud en proyección WGS84.  
* La consulta y descarga de los datos se realiza a través de Google Earth Engine (GEE), para ello es necesario crear una cuenta con dicha aplicación (https://earthengine.google.com/).
* El código fue implementado en python, por lo cual se sugiere instalar conda o miniconda, para crear un ambiente dedicado a GEE. Adicional GEE requiere autenticar las credenciales del usuario, para dicha configuración por favor consultar: https://developers.google.com/earth-engine/python_install-conda.

## Requisitos

  1. python 3.6, 3.7
  2. Google earth Engine API
 
