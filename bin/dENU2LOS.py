#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
v1.0.0  20220218  Yu Morishita

Convert ENU displacements (or velocities) to LOS.

=====
Usage
=====
dENU2LOS.py -i displ_file -e LOS_E_GeoTIFF -n LOS_N_GeoTIFF

-i  Input ENU displacement txt file
    Format: stcode lat lon dE dN dU
-e  GeoTIFF file of LOS EW component
-n  GeoTIFF file of LOS NS component

Output: ${displ_file%.txt}_LOS.txt
            Format: stcode lat lon dE dN dU dLOS LOS_E LOS_N LOS_U
        ${displ_file%.txt}_LOS.geojson

"""
# %% Change log
'''
v1.0.0  20220218  Yu Morishita
- Original implementation
'''


# %% Import
import getopt
import os
import sys
import time
import json
import numpy as np
from osgeo import gdal
import GNSSlib

class Usage(Exception):
    """Usage context manager"""

    def __init__(self, msg):
        self.msg = msg


# %% Main
def main(argv=None):

    # %% Check argv
    if argv == None:
        argv = sys.argv

    start = time.time()
    ver = "1.0.0"
    date = 20220218
    author = "Y. Morishita"
    print("\n{} ver{} {} {}".format(os.path.basename(
        argv[0]), ver, date, author), flush=True)
    print("{} {}".format(os.path.basename(
        argv[0]), ' '.join(argv[1:])), flush=True)

    # %% Set default
    displ_file = None
    LOS_E_file = None
    LOS_N_file = None

    # %% Read options
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "hi:e:n:", ["help"])
        except getopt.error as msg:
            raise Usage(msg)
        for o, a in opts:
            if o == '-h' or o == '--help':
                print(__doc__)
                return 0
            elif o == '-i':
                displ_file = a
            elif o == '-e':
                LOS_E_file = a
            elif o == '-n':
                LOS_N_file = a

        if displ_file is None:
            raise Usage('No input file given, -i is not optional!')
        elif not os.path.exists(displ_file):
            raise Usage('No {} exists!'.format(displ_file))
        if LOS_E_file is None:
            raise Usage('No LOS_E_file file given, -e is not optional!')
        elif not os.path.exists(LOS_E_file):
            raise Usage('No {} exists!'.format(LOS_E_file))
        if LOS_N_file is None:
            raise Usage('No LOS_N_file file given, -n is not optional!')
        elif not os.path.exists(LOS_N_file):
            raise Usage('No {} exists!'.format(LOS_N_file))

    except Usage as err:
        print("\nERROR:", file=sys.stderr, end='')
        print("  "+str(err.msg), file=sys.stderr)
        print("\nFor help, use -h or --help.\n", file=sys.stderr)
        return 2

    # %% Directory and file setting
    LOS_txt = displ_file.replace('.txt', '')+'_LOS.txt'
    LOS_geojson = displ_file.replace('.txt', '')+'_LOS.geojson'

    # %% Read LOS
    geotiff = gdal.Open(LOS_E_file)
    LOS_E = geotiff.ReadAsArray()
    width = geotiff.RasterXSize
    length = geotiff.RasterYSize
    lon_w_p, dlon, _, lat_n_p, _, dlat = geotiff.GetGeoTransform()
    ## lat lon are in pixel registration. dlat is negative
    lon_w = lon_w_p + dlon/2 # grid registration
    lat_n = lat_n_p + dlat/2

    LOS_N = gdal.Open(LOS_N_file).ReadAsArray()
    LOS_E[LOS_E==0] = np.nan
    LOS_N[LOS_N==0] = np.nan

    if LOS_E.shape != LOS_E.shape:
        raise Exception('Size of {} and {} are not identical!'.format(
            LOS_E_file, LOS_N_file))

    # %% Read displ_file
    df_displ = GNSSlib.read_dENU(displ_file)

    # %% Each pos
    features_out_list = []
    f1 = open(LOS_txt, 'w')
    print('{:6s} {:>10s} {:>11s} {:>5s} {:>5s} {:>5s} {:>7s} '
          '{:>5s} {:>5s} {:>5s}'.format('#code', 'lat', 'lon',
                                        'vel_E', 'vel_N', 'vel_U', 'vel_LOS',
                                        'LOS_E', 'LOS_N', 'LOS_U'), file=f1)

    for st in df_displ.itertuples():
        stcode = st.Index

        ### Identify x/y from lat/lon
        x = int(np.round((st.lon-lon_w)/dlon))
        y = int(np.round((st.lat-lat_n)/dlat))

        if x >= width or x < 0 or y >= length or y < 0: ### If outside of area
            _LOS_E = _LOS_N = _LOS_U = np.nan
        else: ### Inside
            _LOS_E = LOS_E[y, x]
            _LOS_N = LOS_N[y, x]
            _LOS_U = np.sqrt(1 - (_LOS_E**2+_LOS_N**2))
            if np.iscomplex(_LOS_U):
                _LOS_U = 0

        dLOS = st.dE*_LOS_E + st.dN*_LOS_N + st.dU*_LOS_U

        # txt
        print('{:6s} {:10.7f} {:11.7f} {:5.1f} {:5.1f} {:5.1f} {:7.1f} '
              '{:5.3f} {:5.3f} {:5.3f}'.format(stcode, st.lat, st.lon,
                                     st.dE, st.dN, st.dU, dLOS,
                                     _LOS_E, _LOS_N, _LOS_U), file=f1)

        if not np.isnan(dLOS):
            # JSON
            feature_prop = {}
            feature_geom = {}
            feature_prop["name"] = stcode
            feature_prop["dE"] = st.dE
            feature_prop["dN"] = st.dN
            feature_prop["dU"] = st.dU
            feature_prop["dLOS"] = dLOS
            feature_geom["type"] = "Point"
            feature_geom["coordinates"] = [st.lon, st.lat]

            feature_out = {"type": "Feature", "properties": feature_prop,
                           "geometry": feature_geom}

            features_out_list.append(feature_out)

    f1.close()

    # %% Create GeoJSON
    jsonout_dict = {'type': 'FeatureCollection',
                    'features': features_out_list}

    with open(LOS_geojson, 'w') as f:
        json.dump(jsonout_dict, f, indent=2)

    # %% Finish
    elapsed_time = time.time()-start
    hour = int(elapsed_time/3600)
    minite = int(np.mod((elapsed_time/60), 60))
    sec = int(np.mod(elapsed_time, 60))
    print("\nElapsed time: {0:02}h {1:02}m {2:02}s".format(hour, minite, sec))

    print('\n{} Successfully finished!!\n'.format(os.path.basename(argv[0])))
    print('Output: {}'.format(LOS_txt))
    print('        {}\n'.format(LOS_geojson))


# %% main
if __name__ == "__main__":
    sys.exit(main())
