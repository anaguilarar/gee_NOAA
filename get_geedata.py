import ee
import functools
import operator
import numpy as np
import pandas as pd
import datetime as dt
import matplotlib.pyplot as plt

## TODO: check if multiple points are located in the same pixel
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

GLDA_bands = ['AvgSurfT_inst',
              'Tair_f_inst',
              'SoilMoi0_10cm_inst',
              'Evap_tavg',
              'PotEvap_tavg',
              'Lwnet_tavg',
              'LWdown_f_tavg',
              'Qair_f_inst',
              'Psurf_f_inst',
              'Rainf_f_tavg']



class gee_weatherdata:
    """Download weathaer data from the Google Earth Engine platform.

           the final output will be a table in which the user can quaery the weather data foor especifical spatial point throughout time.


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
                 mission):

        ### mission reference setting

        if mission == "noaa":
            self._mission = 'NOAA/CFSV2/FOR6H'
            self._bands = NOAA_bands
        if mission == 'gldas':
            self._mission = 'NASA/GLDAS/V021/NOAH/G025/T3H'
            self._bands = GLDA_bands
        if mission == 'chirps':
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

    def _extract_multifunction(self, ee_sp, bands):

        if len(bands) == 4:
            avgdata = self._extract_data(bands[0], ee_sp)
            dfavgdata = fromeedict_todataframe(avgdata, bands[0])
            print('average features processed')
            sumdata = self._extract_data(bands[1], ee_sp, 'sum')
            dfsumdata = fromeedict_todataframe(sumdata, bands[1])
            print('cummulative features processed')
            mindata = self._extract_data(bands[2], ee_sp, 'min')
            dfmindata = fromeedict_todataframe(mindata, bands[2])
            print('minimum features processed')
            maxata = self._extract_data(bands[3], ee_sp, 'max')
            dfmaxata = fromeedict_todataframe(maxata, bands[3])
            print('maximum features processed')

            datasummarised = pd.concat([dfsumdata, dfavgdata,
                                        dfmaxata, dfmindata], axis=1)
        else:
            avgdata = self._extract_data(self._bands, ee_sp)
            datasummarised = fromeedict_todataframe(avgdata,self._bands)

        return [datasummarised, avgdata]

    def _extract_data(self, bands, ee_sp, function='mean'):

        if self._mission == 'NOAA/CFSV2/FOR6H':
            imagecoll = imagecollection_fromlistiteration(self.image_collection.select(bands), self._dates, function)
        if self._mission == 'NASA/GLDAS/V021/NOAH/G025/T3H':
            imagecoll = imagecollection_fromlistiteration(self.image_collection.select(bands), self._dates, function)
        if self._mission == 'UCSB-CHG/CHIRPS/DAILY':
            imagecoll = self.image_collection.filterDate(self._dates[0], self._dates[1]).select(bands)

        dataextracted = imagecoll.map(lambda x: reduce_region(x, ee_sp))
        dataextracted = dataextracted.flatten().getInfo()

        return dataextracted

    def _get_dates_coordinatesfromee(self, geedict):

        date = pd.Series(
            getfeature_fromeedict(geedict, 'properties', 'date'))

        if self._mission == 'NOAA/CFSV2/FOR6H' or self._mission == 'NASA/GLDAS/V021/NOAH/G025/T3H':
            dates = date.apply(lambda x:
                                            dt.timedelta(days=int(x)) + dt.datetime.strptime(self._dates[0],
                                                                                             '%Y-%m-%d'))
        if self._mission == 'UCSB-CHG/CHIRPS/DAILY':
            dates = date

        coords = pd.DataFrame(getfeature_fromeedict(geedict, 'geometry', 'coordinates'),
                              columns=['longitude', 'latitude'])

        return [dates, coords]

    def _extract_databypieces(self,  bands, features, steps=2):

        dataextracted = []

        for spoint in range(0, len(features), steps):
            sp_features = organice_coordinates(features[spoint: spoint + steps])
            ee_sp = ee.FeatureCollection([ee.Geometry.Point(sp_features[i]) for i in range(len(sp_features))])

            summarised, eedict = self._extract_multifunction(ee_sp, bands)

            date, coordinates = self._get_dates_coordinatesfromee(eedict)
            summarised['date'] = date

            dataextracted.append(pd.concat([summarised, coordinates], axis=1))

            print("Points from {0} to {1} were extracted".format(spoint, (spoint + steps) - 1))

        return pd.concat(dataextracted)




    def CHIRPSdata_asdf(self, by ="days"):
        if self._mission == 'UCSB-CHG/CHIRPS/DAILY':
            try:
                dataperday = self._extract_data(self._bands, self._ee_sp)

                date = pd.Series(getfeature_fromeedict(dataperday, 'properties', 'date'))

                coords = pd.DataFrame(getfeature_fromeedict(dataperday, 'geometry', 'coordinates'),
                                      columns=['longitude', 'latitude'])

                df = fromeedict_todataframe(dataperday, self._bands)
                df['date'] = date.apply(lambda date_i:
                                        dt.datetime.strptime(date_i, '%Y%m%d'))
                df = pd.concat([df, coords], axis=1)

            except:
                datedif = dt.datetime.strptime(self._dates[1], "%Y-%m-%d") - dt.datetime.strptime(self._dates[0],
                                                                                                  "%Y-%m-%d")
                step = int(np.ceil(4900 / datedif.days))
                print(
                    'an exception was genereted, query aborted after accumulating over 5000 elements, running by {} features'
                    .format(step))

                listindex, featuresreduced = self.check_duplicatedvalues()

                df = self._extract_databypieces(self._bands, featuresreduced, step)
                df = organice_duplicatedf(df, listindex, featuresreduced, self.features)

            return df

    def plot_CHIRPS(self, feature_index=1, fig_size=[12, 5]):
        if self._mission == 'UCSB-CHG/CHIRPS/DAILY':

            ref_long = self.features.longitude.loc[feature_index - 1]
            plotdata = self.CHIRPSdata_asdf()
            plotdata = plotdata.loc[np.round(plotdata.longitude, 5) == np.round(ref_long, 5)]

            plt.figure(figsize=fig_size)
            plt.plot(plotdata.date, plotdata[self._bands].values)
            plt.ylabel(self._bands[0])
            plt.xlabel("Dates")
            plt.title("longitude: {}; latitude: {}".format(np.round(ref_long, 4),
                                                           np.round(self.features.latitude.loc[feature_index], 4)))
            plt.show()

    def _summarisedatafromcheck(self, imagecoll, ee_sp):
        dataextracted = imagecoll.map(lambda x: reduce_region(x, ee_sp))
        dataextracted = dataextracted.flatten().getInfo()

        summarised = fromeedict_todataframe(dataextracted, self._bands)
        date, coordinates = self._get_dates_coordinatesfromee(dataextracted)
        summarised['date'] = date
        summarised = summarised.loc[summarised.date == summarised.date.values[0]]

        listduplicated = find_duplicatedvalues(summarised[self._bands].mean(axis=1))
        return [listduplicated, self.features.iloc[[x[0] for x in listduplicated]]]

    def check_duplicatedvalues(self):
        """
        query data in gee for two dates only, just to check which coordinades share similar data
        :return:
        """
        listdup =None
        featuresred = None
        ee_sp = self._ee_sp

        if self._mission == 'NOAA/CFSV2/FOR6H' or self._mission == 'NASA/GLDAS/V021/NOAH/G025/T3H':
            imagecoll = imagecollection_fromlistiteration(self.image_collection.select(self._bands), [ee.Date(self._dates[0]),ee.Date(self._dates[0]).advance(1, 'day')], 'mean')

            listdup, featuresred = self._summarisedatafromcheck(imagecoll, ee_sp)

        if self._mission == 'UCSB-CHG/CHIRPS/DAILY':
            imagecoll = self.image_collection.filterDate(ee.Date(self._dates[0]), ee.Date(self._dates[0]).advance(150, 'day')).select(self._bands)
            dataextracted = ee.Image(imagecoll.sum()).reduceRegions(ee_sp, 'mean', 10, crs='EPSG:4326')
            dataextracted = dataextracted.getInfo()
            cummulative_features = getfeature_fromeedict(dataextracted, 'properties', 'mean')
            listdup = find_duplicatedvalues(cummulative_features)
            featuresred = self.features.iloc[[x[0] for x in listdup]]

        return [listdup, featuresred]



    def summarise_hourlydata(self,
                       averagecols=None,
                       cummulativecols=None,
                       minimumcols=None,
                       maximumcols=None,
                       by="days"):
        if self._mission == 'UCSB-CHG/CHIRPS/DAILY' or self._mission =='NASA/GLDAS/V021/NOAH/G025/T3H' or self._mission == 'NOAA/CFSV2/FOR6H':
            '''resume noaa  and gldas data per a specific period'''

            ### group data by days
            ##TODO: create monthly and yearly module

            ### due to the pixel size some features will have same values, so the code will download only one value for multiple features

            if by == "days":
                if cummulativecols is None and averagecols is None and minimumcols is None and maximumcols is None:
                    ee_sp = self._ee_sp
                    summarised = self._extract_multifunction(ee_sp, self._bands)

                else:
                    if cummulativecols is not None and averagecols is not None and minimumcols is not None and maximumcols is not None:
                        try:
                            ee_sp = self._ee_sp

                            summarised, eedict = self._extract_multifunction(ee_sp, [averagecols, cummulativecols,
                                                                                    minimumcols, maximumcols])
                            date, coordinates = self._get_dates_coordinatesfromee(eedict)
                            summarised['date'] = date
                            summarised = pd.concat([summarised, coordinates], axis=1)


                        except:
                            ### Looking for those points that have repeated information
                            listindex, featuresreduced = self.check_duplicatedvalues()
                            datedif = dt.datetime.strptime(self._dates[1], "%Y-%m-%d") - dt.datetime.strptime(self._dates[0], "%Y-%m-%d")
                            step = int(np.floor(4900 / datedif.days))

                            print('generated an exception, query aborted after accumulating over 5000 elements, running by {} features'
                                  .format(step))
                            summarised = self._extract_databypieces([averagecols, cummulativecols,
                                                                     minimumcols, maximumcols], featuresreduced, step)

                            summarised = organice_duplicatedf(summarised, listindex, featuresreduced, self.features)

        return summarised



### extract data using a gee feature collection


def find_duplicatedvalues(listvalues):
    """
    function to identify which indexes are repeated in a list, it returns which values are repeated
    :param listvalues: list
    :return: list
    """
    listduplicated = []
    rangefor = range(0,len(listvalues))

    while(len(rangefor)>0):
        i = rangefor[0]
        similar = [i]

        for j in rangefor:
            if listvalues[i] == listvalues[j] and i != j:
                similar.append(j)

        listduplicated.append(similar)
        rangefor = [x for x in rangefor if x not in functools.reduce(operator.concat, listduplicated)]

    return listduplicated

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
    aux = []
    for feature in range(len(eecollection['features'])):
        ## get data from the dictionary
        datadict = eecollection['features'][feature][attribute][featname]
        ## check if it has info
        aux.append(datadict)
    return (aux)

def check_features(featurecoll, featurenames):
    featurenotinlist = []
    for feature in range(len(featurecoll['features'])):
        featurelist = featurecoll['features'][feature]['properties']
        ## check if all the features were captured by the mission
        difflist = [x for x in featurenames if x not in list(featurelist.keys())[:-1]]
        ## if there is a feature that was not measured by that date update the dictionary with NA
        if len(difflist) > 0:
            if len(list(featurelist.keys())[:-1]) > 0:
                for i in difflist:
                    featurecoll['features'][feature]['properties'].update({i: np.nan})
            else:
                for i in featurenames:
                    featurecoll['features'][feature]['properties'].update({i: np.nan})

def check_1feature(featurecoll):
    for feature in range(len(featurecoll['features'])):
        featurelist = featurecoll['features'][feature]['properties']
        ## if there is are no features update a new value with NA
        if len(list(featurelist.keys())) == 1:
            featurecoll['features'][feature]['properties'].update({'mean': np.nan})


def multiple_eedict_todf(eelistcollection, bands, colname='properties'):
    featuresvalues = []
    for i in range(len(eelistcollection)):
        featuresvalues.append(fromeedict_todataframe(eelistcollection[i], bands, colname))

    return pd.concat(featuresvalues, axis=1)


def fromeedict_todataframe(featurecollection, bands, colname='properties'):
    """get band values from a feature collection"""

    listinfo = []
    ## check if there are empty values in the dictionary
    if len(bands) > 1:
        check_features(featurecollection, bands)
    if len(bands) == 1:
        check_1feature(featurecollection)

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


def read_pointsas_ee_sp(filename, featurestoselect = None):
    '''organice coordinates as gee spatial feature collection'''
    ### read csv file
    try:
        sp_features = pd.read_csv(filename)
    except:
        sp_features = pd.read_csv(filename, encoding="ISO-8859-1")

    if featurestoselect is not None:
        sp_features = sp_features[featurestoselect[0]:featurestoselect[1]]

    sp_features = organice_coordinates(sp_features)
    return ee.FeatureCollection([ee.Geometry.Point(sp_features[i]) for i in range(len(sp_features))])


def read_pointsas_df(filename, colnames=["longitude", "latitude"]):
    '''organice coordinates as gee spatial feature collection'''
    ### read csv file
    sp_features = pd.read_csv(filename)
    sp_features = organice_coordinates(sp_features)
    sp_df = pd.DataFrame(sp_features, columns=colnames)
    sp_df['index'] = [i + 1 for i in range(len(sp_features))]
    return sp_df


def organice_coordinates(dataframe, longcolname="longitude", latcolname="latitude"):
    '''organice the coordinates as long, lat per list'''

    return [[dataframe[longcolname].iloc[row_i], dataframe[latcolname].iloc[row_i]] for row_i in
            range(dataframe.shape[0])]


def reduce_region(image, geometry):
    """Spatial aggregation function for a single image and a polygon feature"""
    stat_dict = image.reduceRegions(geometry, 'mean', 10, crs='EPSG:4326');
    return stat_dict.map(lambda y: y.set('date', image.get('system:index')))

def organice_duplicatedf(df, list_index, reduced_features, orig_features):
    idcoords = [str(df.longitude.values[i]) + str(df.latitude.values[i]) for i in
                range(len(df.latitude.values))]
    dataaux = []
    for j in range(len(list_index)):
        alldata = []
        refpos = str(reduced_features.longitude.values[j]) + str(reduced_features.latitude.values[j])
        pddata = df.loc[np.array(idcoords) == refpos]
        ### assign values to the repeated features
        for i in list_index[j]:
            pddata.longitude = orig_features.iloc[i].longitude
            pddata.latitude = orig_features.iloc[i].latitude
            alldata.append(pddata)

        dataaux.append(pd.concat(alldata))

    return pd.concat(dataaux)
