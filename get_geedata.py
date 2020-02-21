import ee
import numpy as np
import pandas as pd
import datetime as dt
import matplotlib.pyplot as plt

ee.Initialize()

NOAA_bands = ['Downward_Long-Wave_Radp_Flux_surface_6_Hour_Average',
              'Downward_Short-Wave_Radiation_Flux_surface_6_Hour_Average',
              'Maximum_temperature_height_above_ground_6_Hour_Interval',
              'Minimum_temperature_height_above_ground_6_Hour_Interval',
              'Precipitation_rate_surface_6_Hour_Average',
              'Pressure_surface',
              'Potential_Evaporation_Rate_surface_6_Hour_Average',
              'Specific_humidity_height_above_ground',
              'Upward_Long-Wave_Radp_Flux_surface_6_Hour_Average',
              'Upward_Short-Wave_Radiation_Flux_surface_6_Hour_Average']


class gee_weatherdata:
    """Download weathaer data from the Google Earth Engine platform.

           the final output will be table in which the user can quaery the weather data foor especifical spatial point throughout time.


           Parameters
           ----------

           start_date : str
                   The start of the time period used for data extraction, it must have the folowing format "YYYY-MM-DD"

           end_date : str
                  The end of the time period used for data extraction., it must have the following format "YYYY-MM-DD"

           roi_filename : str
                   string path to a csv file that must conatins the coordinates in longitude and latitude format

           output_path : str
                   string path to a destination folder

           mission : str
                   id reference to the stallite which will be processed:
                       - NOAA: "noaa"
                       - CHIRPS: "chirsp"

           Attributes
           ----------
           products : dict
               Filtered copy of `product_list` passed to the object containing only
               products generated between `start_date` and `end_date`.
           product_boundaries : dict
               Contains `shapely.geometry.Polygon` objects describing the boundaries
               of each of the products in `products`.
    """

    def __init__(self, start_date,
                 end_date,
                 roi_filename,
                 mission,
                 output_path=None,
                 bands=None):

        ### mission reference setting

        if (mission == "noaa"):
            self._mission = 'NOAA/CFSV2/FOR6H'
            self._bands = NOAA_bands
        else:
            self._mission = 'UCSB-CHG/CHIRPS/DAILY'
            self._bands = ["precipitation"]

        self._dates = [start_date, end_date]
        ## get spatial points
        self._ee_sp = read_pointsas_ee_sp(roi_filename)
        self.features = read_pointsas_df(roi_filename)
        ###
        self.image_collection = query_image_collection(ee.Date(start_date),
                                                       ee.Date(end_date).advance(1, 'day'),
                                                       self._mission,
                                                       self._ee_sp)

    def _extract_data(self, bands, function='mean'):
        if self._mission == 'NOAA/CFSV2/FOR6H':
            imagecoll = imagecollection_fromlistiteration(self.image_collection.select(bands), self._dates, function)
        else:
            imagecoll = self.image_collection.filterDate(self._dates[0], self._dates[1]).select(bands)

        dataextracted = imagecoll.map(lambda x: reduce_region(x, self._ee_sp))
        dataextracted = dataextracted.flatten().getInfo()
        return dataextracted

    def CHIRPSdata_asdf(self, by="days"):
        if self._mission == 'UCSB-CHG/CHIRPS/DAILY':
            dataperday = self._extract_data(self._bands)

            date = pd.Series(getfeature_fromeedict(dataperday, 'properties', 'date'))

            coords = pd.DataFrame(getfeature_fromeedict(dataperday, 'geometry', 'coordinates'),
                                  columns=['longitude', 'latitude'])

            df = fromeedict_todataframe(dataperday, self._bands)
            df['date'] = date.apply(lambda date_i:
                                    dt.datetime.strptime(date_i, '%Y%m%d'))
            df = pd.concat([df, coords], axis=1)
            return df

    def plot_CHIRPS(self, feature_index = 1, fig_size = [12,5]):
        if self._mission == 'UCSB-CHG/CHIRPS/DAILY':

            ref_long = self.features.longitude.loc[feature_index-1]
            plotdata = self.CHIRPSdata_asdf()
            plotdata = plotdata.loc[np.round(plotdata.longitude, 5) == np.round(ref_long, 5)]

            plt.figure(figsize=fig_size)
            plt.plot(plotdata.date, plotdata[self._bands].values)
            plt.ylabel(self._bands[0])
            plt.xlabel("Dates")
            plt.title("longitude: {}; latitude: {}".format(np.round(ref_long, 4),
                                                           np.round(self.features.latitude.loc[feature_index], 4)))
            plt.show()


    def summarise_noaa(self,
                       averagecols=None,
                       cummulativecols=None,
                       minimumcols=None,
                       maximumcols=None,
                       by="days"):

        '''resume noaa data per a specific time'''

        ### group data by days
        ##TODO: create monthly and yearly module
        if by == "days":
            if cummulativecols is None and averagecols is None and minimumcols is None and maximumcols is None:
                avgdata = self._extract_data(self._bands)


            else:
                if cummulativecols is not None and averagecols is not None and minimumcols is not None and maximumcols is not None:
                    avgdata = self._extract_data(averagecols)
                    avgdata = fromeedict_todataframe(avgdata, averagecols)
                    print('average features processed')
                    sumdata = self._extract_data(cummulativecols, 'sum')
                    sumdata = fromeedict_todataframe(sumdata, cummulativecols)
                    print('cummulative features processed')
                    mindata = self._extract_data(minimumcols, 'min')
                    mindata = fromeedict_todataframe(mindata, minimumcols)
                    print('minimum features processed')
                    maxata = self._extract_data(maximumcols, 'max')
                    maxata = fromeedict_todataframe(maxata, maximumcols)
                    print('maximum features processed')

                    summarised = pd.concat([sumdata, avgdata,
                                            maxata, mindata],
                                           axis=1)

        date = pd.Series(getfeature_fromeedict(self._extract_data(averagecols), 'properties', 'date'))

        summarised['date'] = date.apply(lambda date:
                                        dt.timedelta(days=int(date)) + dt.datetime.strptime(self._dates[0],
                                                                                            '%Y-%m-%d'))
        coords = pd.DataFrame(getfeature_fromeedict(self._extract_data(averagecols), 'geometry', 'coordinates'),
                              columns=['longitude', 'latitude'])

        summarised = pd.concat([summarised, coords], axis=1)

        return summarised


'''
            if cummulativecolumns is not None:
                averagecolumns = [x for x in self._bands if not x in cummulativecolumns]
                bydays_cumm = self._imagecollection_fromlistiteration(cummulativecolumns, 'sum')
                bydays_avg = self._imagecollection_fromlistiteration(averagecolumns, 'mean')

                if averagecolumns is not None:
                    cummulativecolumns = [x for x in self._bands if not x in averagecolumns]
                    bydays_cumm = self._imagecollection_fromlistiteration(cummulativecolumns, 'sum')
                    bydays_avg = self._imagecollection_fromlistiteration(averagecolumns, 'mean')
'''


### extract data using a gee feature collection


def get_coordinates(featurecol):
    """get dates from a feature collection"""
    longitudes = [featurecol['features'][feature]['geometry']['coordinates'][0] for feature in
                  range(len(featurecol['features']))]
    latitudes = [featurecol['features'][feature]['geometry']['coordinates'][1] for feature in
                 range(len(featurecol['features']))]
    return longitudes, latitudes


def query_image_collection(initdate, enddate, satellite_mission, ee_sp):
    '''mission data query'''

    ## mission data query
    return ee.ImageCollection(satellite_mission).filterDate(initdate, enddate).filterBounds(ee_sp)


def imagecollection_fromlistiteration(imagecollection, dates, function):
    ## calculate diff days

    diffdays = ee.Number(ee.Date(dates[1]).millis()).subtract(
        ee.Number(ee.Date(dates[0]).millis())).divide(1000).divide(86400)

    return ee.ImageCollection(
        ee.List.sequence(0, diffdays).map(
            lambda n: summarisebydates(imagecollection, n, dates,
                                       function=function)))


def getfeature_fromeedict(eecollection, attribute, featname):
    return [eecollection['features'][feature][attribute][featname]
            for feature in range(len(eecollection['features']))]


def check_features(featurecoll, featurenames):
    featurenotinlist = []
    for feature in range(len(featurecoll['features'])):
        featurelist = featurecoll['features'][feature]['properties']
        difflist = [x for x in featurenames if x not in list(featurelist.keys())[:-1]]
        if len(difflist) > 0:
            featurenotinlist.append(difflist)

    return np.unique(featurenotinlist)


def fromeedict_todataframe(featurecollection, bands, colname='properties'):
    """get band values from a feature collection"""

    listinfo = []
    if len(bands) > 1:
        notinall = check_features(featurecollection, bands)

        bands = [x for x in bands if x not in notinall]

    for i in bands:

        if len(bands) < 2:
            i = 'mean'

        listinfo.append(
            getfeature_fromeedict(featurecollection, colname, i)
        )

    dfinfo = pd.DataFrame(np.transpose(listinfo), columns=bands)

    return dfinfo


def summarisebydates(image_collection, step, dates, function='mean'):
    ini = ee.Date(dates[0]).advance(step, 'day')
    end = ini.advance(1, 'day')

    if (function == 'mean'):
        resumeddata = image_collection.filterDate(ini, end).mean().set('system:time_start', ini)

    if (function == 'sum'):
        resumeddata = image_collection.filterDate(ini, end).sum().set('system:time_start', ini)

    if (function == 'max'):
        resumeddata = image_collection.filterDate(ini, end).max().set('system:time_start', ini)

    if (function == 'min'):
        resumeddata = image_collection.filterDate(ini, end).min().set('system:time_start', ini)
    return resumeddata


def read_pointsas_ee_sp(filename):
    '''organice coordinates as gee spatial feature collection'''
    ### read csv file
    sp_features = pd.read_csv(filename)
    sp_features = organice_coordinates(sp_features)
    return ee.FeatureCollection([ee.Geometry.Point(sp_features[i]) for i in range(len(sp_features))])

def read_pointsas_df(filename, colnames = ["longitude", "latitude"]):
    '''organice coordinates as gee spatial feature collection'''
    ### read csv file
    sp_features = pd.read_csv(filename)
    sp_features = organice_coordinates(sp_features)
    sp_df = pd.DataFrame(sp_features, columns= colnames)
    sp_df['index'] = [i+1 for i in range(len(sp_features))]
    return sp_df


def organice_coordinates(dataframe, longcolname="longitude", latcolname="latitude"):
    '''organice the coordinates as long, lat per list'''

    return [[dataframe[longcolname][row_i], dataframe[latcolname][row_i]] for row_i in range(dataframe.shape[0])]


def reduce_region(image, geometry):
    """Spatial aggregation function for a single image and a polygon feature"""
    stat_dict = image.reduceRegions(geometry, 'mean', 10, crs='EPSG:4326');
    return stat_dict.map(lambda y: y.set('date', image.get('system:index')))
