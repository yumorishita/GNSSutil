#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
v0.0.1  20220223  Yu Morishita

Calculate displacements from *.pos files.
Not completed.

====================
Input & output files
====================
Inputs in POSDIR:
- POSCOR/*.pos
  Format: yyyy mm dd time lat lon h X Y Z
- stations.txt
  Format: stcode  lat  lon  h

Outputs:
- displacement.txt
  Format: stcode dE dN dU s_date1 s_date2 s_n_day e_date1 e_date2 e_n_day
- displacement.geojson


=====
Usage
=====
calc_displacement.py -s s_date1 [--s2 s_date2] -e e_date1 [--e2 e_date2]
 [-d n_day_thre] [-a lon_w/lon_e/lat_s/lat_n] [-r stcode_ref] [--POSDIR POSDIR]

-s  Start date of reference period (yyyymmdd) for displacement calculation
--s2  End date of reference period (yyyymmdd) (Default: s_date1+13days)
-e  Start date of secondary period (yyyymmdd)
--e2  End date of secondary period (yyyymmdd) (Default: e_date1+13days)
-d  Threshold of the number of available days (Default: 0)
-a  Area (Default: all stations)
-r  Reference station code (Default: no reference)
--POSDIR  Path to directory containing POSCOR and station.txt
          (Default: $POSDIR)

"""
# %% Change log
'''
v1.0.0  20220223  Yu Morishita
- Original implementation (Not completed)
'''


# %% Import
import getopt
import os
import sys
import time
import json
import numpy as np
import datetime as dt
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
    date = 20220223
    author = "Y. Morishita"
    print("\n{} ver{} {} {}".format(os.path.basename(
        argv[0]), ver, date, author), flush=True)
    print("{} {}".format(os.path.basename(
        argv[0]), ' '.join(argv[1:])), flush=True)

    # %% Set default
    s_date1 = None
    s_date2 = None
    e_date1 = None
    e_date2 = None
    n_day_thre = 0
    area = None
    ref = None
    POSDIR = os.environ['POSDIR']

    # %% Read options
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "hs:e:d:a:r:",
                                       ["help", "s2=", "e2=", "POSDIR="])
        except getopt.error as msg:
            raise Usage(msg)
        for o, a in opts:
            if o == '-h' or o == '--help':
                print(__doc__)
                return 0
            elif o == '-s':
                s_date = a
            elif o == '-e':
                e_date = a
            elif o == '-d':
                n_day_thre = int(a)
            elif o == '-a':
                area = a
            elif o == '-r':
                ref = a
            elif o == '--s2':
                s_date2 = a
            elif o == '--e2':
                e_date2 = a
            elif o == '--POSDIR':
                POSDIR = a

        if s_date1 is None:
            raise Usage('No s_date1 given, -s is not optional!')
        if e_date1 is None:
            raise Usage('No e_date1 given, -e is not optional!')
        if not os.path.exists(POSDIR):
            raise Usage('No {} exists, specify by --POSDIR option!'.format(
                POSDIR))

    except Usage as err:
        print("\nERROR:", file=sys.stderr, end='')
        print("  "+str(err.msg), file=sys.stderr)
        print("\nFor help, use -h or --help.\n", file=sys.stderr)
        return 2

    # %% Directory and file setting
    POSDIR = os.path.abspath(POSDIR)
    poscor_dir = os.path.join(POSDIR, 'POSCOR')
    station_file = os.path.join(POSDIR, 'stations.txt')

    d_day_default = 13
    if s_date2 is None:
        s_date2 = dt.datetime.strptime(s_date1, "%Y%m%d") + \
                    dt.timedelta(days=d_day_default)
    if e_date2 is None:
        e_date2 = dt.datetime.strptime(e_date1, "%Y%m%d") + \
                    dt.timedelta(days=d_day_default)

    displ_txt = 'displacement.txt'
    displ_geojson = 'displacement.geojson'

    if not os.path.exists(poscor_dir):
        raise FileNotFoundError('No {} exists!'.format(poscor_dir))
    if not os.path.exists(station_file):
        raise FileNotFoundError('No {} exists!'.format(station_file))

    if area is not None:
        lon_w, lon_e, lat_s, lat_n = [float(s) for s in area.split('/')]
    else:
        lat_s = -90
        lat_n = 90
        lon_w = -180
        lon_e = 180

    # %% Read stations
    df_station = GNSSlib.read_stations(station_file, lat_s=lat_s, lat_n=lat_n,
                                       lon_w=lon_w, lon_e=lon_e)

    # %% Read reference
    if ref is not None:
        refpos_file = os.path.join(poscor_dir, '{}.pos'.format(ref))
        if not os.path.exists(refpos_file):
            raise FileNotFoundError('No {} exists!'.format(refpos_file))

        df_refpos = GNSSlib.read_pos(refpos_file, s_date=s_date, e_date=e_date2)
        lat_ref = df_refpos.lat[0] # Not read from df_station because ref could be outside of area
        lon_ref = df_refpos.lon[0]

        s_n_day = len(df_refpos.loc[s_date1:s_date2])
        if s_n_day < n_day_thre:
            raise ValueError("Number of available days at {} between {} and {}"
                             " is {} < n_day_thre {}".format(
                         ref, s_date1, s_date2, s_n_day, n_day_thre))
        e_n_day = len(df_refpos.loc[e_date1:e_date2])
        if e_n_day < n_day_thre:
            raise ValueError("Number of available days at {} between {} and {}"
                             " is {} < n_day_thre {}".format(
                         ref, e_date1, e_date2, e_n_day, n_day_thre))

        dENU_ref, dENUstd_ref, s_n_day, e_n_day = GNSSlib.calc_displacement(
            df_refpos, s_date1, s_date2, e_date1, e_date2)


    # %% Each pos
    features_out_list = []
    f1 = open(displ_txt, 'w')
    print('{:6s} {:>10s} {:>11s} {:>7s} {:>5s} {:>5s} {:>5s} '
          '{:>8s} {:>8s} {:>8s} {:>8s} {:>4s}'.format('#code', 'lat', 'lon', 'height',
                                        'dE', 'dN', 'dU', 's_date1', 's_date2',
                                        'e_date1', 'e_date2', 'nday'), file=f1)

    for pos in df_station.itertuples():
        stcode = pos.Index
        pos_file = os.path.join(poscor_dir, '{}.pos'.format(stcode))
        if not os.path.exists(pos_file):
            print("No {} exist! Skip".format(pos_file))
            continue
        else:
            print("Calculate displacement at {}...".format(stcode))

        df_pos = GNSSlib.read_pos(pos_file, s_date=s_date1, e_date=e_date2)

        if len(df_pos) < 2:
            print("  Available n_day is {}. Skip.".format(len(df_pos)))
            continue

        dENU, dENUstd, s_n_day, e_n_day = GNSSlib.calc_displacement(
            df_pos, s_date1, s_date2, e_date1, e_date2)

        if s_n_day < n_day_thre:
            print("s_n_day ({}) < {}. Skip".format(s_n_day, n_day_thre))
            continue
        elif e_n_day < n_day_thre:
            print("e_n_day ({}) < {}. Skip".format(e_n_day, n_day_thre))
            continue

        # Relative to ref
        ### 20220223 Ref must be subtracted for each day before computing std!
        if ref is not None:
            dENU = np.array(GNSSlib.dENU_sub_ref(pos.lat, pos.lon,
                           dENU[0], dENU[1], dENU[2], lat_ref, lon_ref,
                           dENU_ref[0], dENU_ref[1], dENU_ref[2]))

        dENUmm = dENU*1000 # m/yr -> mm/yr

        # txt
        print('{:6s} {:10.7f} {:11.7f} {:7.2f} {:5.1f} {:5.1f} {:5.1f} '
              '{:8s} {:8s} {:4d}'.format(stcode, pos.lat, pos.lon, pos.h,
                                     dENUmm[0], dENUmm[1], dENUmm[2],
                                     s_date1, s_date2, e_date1, e_date2, n_day), file=f1)

        # JSON
        feature_prop = {}
        feature_geom = {}
        feature_prop["name"] = stcode
        feature_prop["vel_E"] = vel_ENUmm[0]
        feature_prop["vel_N"] = vel_ENUmm[1]
        feature_prop["vel_U"] = vel_ENUmm[2]
        feature_prop["s_date"] = _s_date
        feature_prop["e_date"] = _e_date
        feature_prop["n_day"] = n_day
        feature_geom["type"] = "Point"
        feature_geom["coordinates"] = [pos.lon, pos.lat]

        feature_out = {"type": "Feature", "properties": feature_prop,
                       "geometry": feature_geom}

        features_out_list.append(feature_out)

    f1.close()

    # %% Create GeoJSON
    jsonout_dict = {'type': 'FeatureCollection',
                    'features': features_out_list}

    with open(displ_geojson, 'w') as f:
        json.dump(jsonout_dict, f, indent=2)

    # %% Finish
    elapsed_time = time.time()-start
    hour = int(elapsed_time/3600)
    minite = int(np.mod((elapsed_time/60), 60))
    sec = int(np.mod(elapsed_time, 60))
    print("\nElapsed time: {0:02}h {1:02}m {2:02}s".format(hour, minite, sec))

    print('\n{} Successfully finished!!\n'.format(os.path.basename(argv[0])))
    print('Output: {}'.format(displ_txt))
    print('        {}\n'.format(displ_geojson))


# %% main
if __name__ == "__main__":
    sys.exit(main())
