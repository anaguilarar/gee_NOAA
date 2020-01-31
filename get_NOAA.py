## todo: date is not compatible
## TODO: CALENDAR RANGE MAKE THE QUERY PER DAYS

import ee
import numpy as np
import pandas as pd

ee.Initialize()

NOAA_bands = ['Downward_Long-Wave_Radp_Flux_surface_6_Hour_Average', 'Downward_Short-Wave_Radiation_Flux_surface_6_Hour_Average',
              'Maximum_temperature_height_above_ground_6_Hour_Interval','Minimum_temperature_height_above_ground_6_Hour_Interval',
              'Precipitation_rate_surface_6_Hour_Average', 'Pressure_surface',
              'Potential_Evaporation_Rate_surface_6_Hour_Average','Specific_humidity_height_above_ground',
              'Upward_Long-Wave_Radp_Flux_surface_6_Hour_Average', 'Upward_Short-Wave_Radiation_Flux_surface_6_Hour_Average']
###
def organice_coordinates(dataframe, longcolname="longitude", latcolname="latitude"):
    '''organice the coordinates as long, lat per list'''

    return [[dataframe[longcolname][row_i], dataframe[latcolname][row_i]] for row_i in range(dataframe.shape[0])]
###
def _reduce_region(image,geometry):
    """Spatial aggregation function for a single image and a polygon feature"""
    stat_dict = image.reduceRegion(ee.Reducer.mean(), geometry,  30,
     crs= 'EPSG:4326');
    return ee.Feature(None, stat_dict)

def _reduce_regions(image,geometry):
    """Spatial aggregation function for a single image and a polygon feature"""
    image = ee.Image(image)
    stat_dict = image.reduceRegions(geometry,'mean',10);
    return stat_dict#.map(lambda y: y.set('date', image.get('system:index')))

def _get_dates(imagecol):
    def _iter_func (image, newlist):
        date = ee.Number.parse(image.date().format("YYYYMMdd"))
        newlist = ee.List(newlist)

        return ee.List(newlist.add(date).sort())

    return imagecol.iterate(_iter_func, ee.List([])).getInfo()

def get_dates(featurecol, colname = 'id'):
    """get dates from a feature collection"""
    return [featurecol['features'][feature][colname] for feature in range(len(featurecol['features']))]

def get_coordinates(featurecol, colname = 'geometry'):

    """get dates from a feature collection"""
    longitudes = [featurecol['features'][feature][colname]['coordinates'][0] for feature in range(len(featurecol['features']))]
    latitudes = [featurecol['features'][feature][colname]['coordinates'][1] for feature in range(len(featurecol['features']))]
    return longitudes, latitudes

def get_values(featurecol, missionbands, colname = 'properties'):
    """get band values from a feature collection"""
    datappend = []

    for i in missionbands:
        datappend.append(
            [featurecol['features'][feature][colname][i]
             for feature in range(len(featurecol['features']))]
        )

    return datappend

###
def extract_NOAA_data_pergeom(init_date, end_date, sp_features, bands = NOAA_bands):
    """Spatial extraction function for a image collection using spatial points"""

    satellite_mission = "NOAA/CFSV2/FOR6H"
    ## mission data query
    imgcoll = ee.ImageCollection(satellite_mission).filterDate(init_date, end_date).select(bands);
    ## transform spatial data ee feature collection
    
    ee_sp = ee.FeatureCollection([ee.Geometry.Point(sp_features[i]) for i in range(len(sp_features))])

    ##extract data
    dataextracted = imgcoll.map(lambda x: _reduce_regions(x, ee_sp))
    dataextracted = dataextracted.flatten().getInfo()

    #if 'date' not in bands:
    #    bands.append('date')

    #dates = get_dates(dataextracted)
    ## get values
    band_values = pd.DataFrame(np.transpose(get_values(dataextracted, bands)), columns = bands)
    #print(band_values.head())

    ## adding dates
    band_values['dates'] = _get_dates(imgcoll)*len(sp_features)
    ## add spatial coordinates
    band_values['longitude'], band_values['latitude'] = get_coordinates(dataextracted)

    return band_values

###
def _summarisedata(df, cumvars, avgvars):
    cumdata = pd.DataFrame(df[cumvars].sum())
    avdata = pd.DataFrame(df[avgvars].mean())

    dfsumm = pd.DataFrame(pd.concat([cumdata, avdata], ignore_index=True, axis=1))
    dfsumm.columns = cumvars + avgvars
    return dfsumm

def summarise_noaa(noaadf, by = "days", averagecolumns = None, cummulativecolumns = None, noaabands = NOAA_bands.copy()):
    ##
    noaabands = [x for x in noaabands if not x in ['dates']]
    ## add dates as a new column
    noaadf['dates']
    noaadf['year_month'] = [str(i)[:-4] + "_" + str(i)[-4:-2] for i in noaadf['dates']]
    noaadf['year'] = [str(i)[:-4] for i in noaadf['dates']]

    ## convert str to numeric
    noaadf[noaabands] = noaadf.apply(lambda x: pd.to_numeric(x[noaabands]), axis=1)
    ## select the agrupation type
    if by == "days":
        columntogroup = "dates"
        noaadf['dates'] = pd.to_datetime(noaadf['dates'], format='%Y%m%d', errors='ignore')
    if by == "month":
        columntogroup = "year_month"
        noaadf['year_month'] = pd.to_datetime(noaadf['year_month'], format='%Y%m', errors='ignore')
    if by == "year":
        columntogroup = "year"
    try:

        ## create group variable
        datagrouped = noaadf.groupby([columntogroup, 'longitude', 'latitude'])
    except:
        raise ValueError("check by options")

    if cummulativecolumns is None and averagecolumns is None:
        datasumm = pd.DataFrame(datagrouped[noaabands].mean())

    else:
        if cummulativecolumns is not None and averagecolumns is not None:
            datasumm = _summarisedata(datagrouped, cummulativecolumns,averagecolumns )
        if cummulativecolumns is not None:
            datasumm = _summarisedata(datagrouped, cummulativecolumns, [x for x in noaabands if not x in cummulativecolumns])
        if averagecolumns is not None:
            datasumm = _summarisedata(datagrouped, [x for x in noaabands if not x in averagecolumns], averagecolumns)

    return datasumm








