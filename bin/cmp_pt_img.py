#!/usr/bin/env python3
"""
v1.0 20220315 Yu Morishita

This script compares GeoTIFF image with pointwise data.

=====
Usage
=====
cmp_pt_img.py -i GeoTIFF -p points.txt [-o cmp.txt] [--deramp deg]

 -i  Input GeoTIFF image
 -p  Input point data txt file
     (Format: stcode lat lon data)
     (Lines start with # are skipped as comments)
 -o  Output txt (Default: cmp.txt)
     *.geojson is also created
     (Format: stcode lat lon pt_data img_data diff)
 --deramp  Deramp with specified degree (0, 1, bl, 2) [not supported yet]

"""
#%% Change log
'''
v1.0 20220315 Yu Morishita
 - Original implementation
'''
"""
To do:
    - deramp

"""


#%% Import
import getopt
import os
import sys
import numpy as np
import time
import json

import GNSSlib
from osgeo import gdal

class Usage(Exception):
    """Usage context manager"""
    def __init__(self, msg):
        self.msg = msg


#%% bl2xy
def bl2xy(lon, lat, width, length, lat1, postlat, lon1, postlon):
    """
    lat1 is north edge and postlat is negative value.
    lat lon values are in grid registration
    x/y index start from 0, end with width-1
    """
    x = int(np.round((lon - lon1)/postlon))
    y = int(np.round((lat - lat1)/postlat))

    return [x, y]


#%% Main
def main(argv=None):

    #%% Check argv
    if argv == None:
        argv = sys.argv

    start = time.time()
    ver='1.0'; date=20220315; author="Y. Morishita"
    prog_name = os.path.basename(argv[0])
    print(f"\n{prog_name} ver{ver} {date} {author}")
    print(f"{prog_name} {' '.join(argv[1:])}")


    #%% Set default
    geotiff_file = []
    pt_file = []
    out_file = 'cmp.txt'
    deg = []


    #%% Read options
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "hi:p:o:",
                                       ["help", "deramp="])
        except getopt.error as msg:
            raise Usage(msg)
        for o, a in opts:
            if o == '-h' or o == '--help':
                print(__doc__)
                return 0
            elif o == '-i':
                geotiff_file = a
            elif o == '-p':
                pt_file = a
            elif o == '-o':
                out_file = True
            elif o == '--deramp':
                deg = a

        if not geotiff_file:
            raise Usage('No GeoTIFF given, -i is not optional!')
        elif not os.path.exists(geotiff_file):
            raise Usage(f'{geotiff_file} does not exist!')
        if not pt_file:
            raise Usage('No point txt given, -p is not optional!')
        elif not os.path.exists(pt_file):
            raise Usage('{} does not exist!'.format(pt_file))
        if deg and not deg in ['0', '1', 'bl', '2']:
            raise Usage('deg must be either 0, 1, bl, or 2!')

    except Usage as err:
        print("\nERROR:", file=sys.stderr, end='')
        print("  "+str(err.msg), file=sys.stderr)
        print("\nFor help, use -h or --help.\n", file=sys.stderr)
        return 2


    #%%  Settings
    out_geojson = out_file.replace('.txt', '.geojson')


    #%% Read GeoTIFF file
    geotiff = gdal.Open(geotiff_file)
    img_data = geotiff.ReadAsArray()
    width = geotiff.RasterXSize
    length = geotiff.RasterYSize
    lon_w_p, dlon, _, lat_n_p, _, dlat = geotiff.GetGeoTransform()
    ## lat lon are in pixel registration. dlat is negative
    lon_w = lon_w_p + dlon/2 # grid registration
    lat_n = lat_n_p + dlat/2

    try:
        nodata = geotiff.GetRasterBand(1).GetNoDataValue()
        if ~np.isnan(nodata) :
            img_data[img_data==nodata] = np.nan
    except:
        pass

    df_pt = GNSSlib.read_stations(pt_file) # h col has data


    # %% Calc difference and output txt
    features_out_list = []
    diff_list = []
    f = open(out_file, 'w')

    for i in range(len(df_pt)):
        x, y = bl2xy(df_pt.lon[i], df_pt.lat[i], width, length,
                                   lat_n, dlat, lon_w, dlon)
        if x < 0 or y < 0 or x >= width or y >= length: # out of area
            continue
        img_data1 = float(img_data[y, x]) # need float64 instead of float32 for json
        diff = df_pt.h[i] - img_data1
        diff_list.append(diff)
        print(f'{df_pt.index[i]:6s} {df_pt.lat[i]:10.7f} {df_pt.lon[i]:11.7f}'
              f'{df_pt.h[i]:7.2f} {img_data1:7.2f} {diff:7.2f}', file=f)

        # JSON
        if ~np.isnan(diff):
            feature_prop = {}
            feature_geom = {}
            feature_prop["name"] = df_pt.index[i]
            feature_prop["pt"] = df_pt.h[i]
            feature_prop["img"] = img_data1
            feature_prop["diff"] = diff
            feature_geom["type"] = "Point"
            feature_geom["coordinates"] = [df_pt.lon[i], df_pt.lat[i]]

            feature_out = {"type": "Feature", "properties": feature_prop,
                           "geometry": feature_geom}

            features_out_list.append(feature_out)

    diff_ar = np.array(diff_list)

    print(f'# Ave: {np.nanmean(diff_ar):7.2f}, Std: {np.nanstd(diff_ar):7.2f}',
          file=f)

    f.close()

    # %% Create GeoJSON
    jsonout_dict = {'type': 'FeatureCollection',
                    'features': features_out_list}

    with open(out_geojson, 'w') as f:
        json.dump(jsonout_dict, f, indent=2)


    #%% Finish
    elapsed_time = time.time()-start
    hour = int(elapsed_time/3600)
    minite = int(np.mod((elapsed_time/60),60))
    sec = int(np.mod(elapsed_time,60))
    print("\nElapsed time: {0:02}h {1:02}m {2:02}s".format(hour,minite,sec))
    print('\n{} Successfully finished!!\n'.format(os.path.basename(argv[0])))
    print('Output: {}'.format(out_file))
    print('        {}'.format(out_geojson))
    print('', flush=True)


#%% main
if __name__ == "__main__":
    sys.exit(main())
