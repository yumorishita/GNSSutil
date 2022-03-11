#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
v1.0.0  20220311  Yu Morishita

Update pos files from terras.gsi.go.jp.

=======================
Output files in $POSDIR
=======================
- stations.txt
- corrf5o.dat
- yyyy/*.pos
- POS/*.pos (merged)
- POSCOR/*.pos (merged and offset corrected)

=====
Usage
=====
update_pos.py [--full] [--POSDIR POSDIR] [--netrc netrc] [--userid str]
 [--passwd str]

--full    Full update instead of delta update
--POSDIR  Output directory (Default: $POSDIR)
--netrc   netrc file containing login information for ftp (Default: ~/.netrc)
          Format (with permission 600):
              machine  hostname
              login    userid
              password passwd
--userid  User ID of the terras FTP server
--passwd  Password of the terras FTP server

Note: If no netrc is specified, no ~/.netrc exists, or --userid and --passwd
      are not specified, the userid and passwd will be interactively asked
      during the script.

"""
# %% Change log
'''
v1.0.0  20220311  Yu Morishita
- Original implementation
'''
"""
To do:

"""


# %% Import
import getopt
import os
import sys
import time
import netrc
import shutil
from pathlib import Path
from ftplib import FTP
from getpass import getpass
import urllib.request
import datetime as dt
import numpy as np
import pandas as pd

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
    date = 20220311
    author = "Y. Morishita"
    print("\n{} ver{} {} {}".format(os.path.basename(
        argv[0]), ver, date, author), flush=True)
    print("{} {}".format(os.path.basename(
        argv[0]), ' '.join(argv[1:])), flush=True)

    # %% Set default
    isdelta = True
    POSDIR = os.environ['POSDIR']
    netrc_file = None
    userid = None
    passwd = None

    url_terras = 'terras.gsi.go.jp'
    dir_f5 = '/data/coordinates_F5/GPS/'
    url_corr = 'https://mekira.gsi.go.jp/JAPANESE/corrf5o.dat'
    st_year = 1996
    isupdate = False


    # %% Read options
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "h", ["help", "full",
                                   "POSDIR=", "netrc=", "userid=", "passwd="])
        except getopt.error as msg:
            raise Usage(msg)
        for o, a in opts:
            if o == '-h' or o == '--help':
                print(__doc__)
                return 0
            elif o == '--full':
                isdelta = False
            elif o == '--POSDIR':
                POSDIR = a
            elif o == '--netrc':
                netrc_file = a
            elif o == '--userid':
                userid = a
            elif o == '--passwd':
                passwd = a

        # Read login information
        if netrc is None and not os.path.exists('~/.netrc'):
            pass
        elif not (userid is None and passwd is None):
            pass
        else:
            if netrc is not None:
                netrc_dict = netrc.netrc(netrc_file).hosts
            else:
                netrc_dict = netrc.netrc().hosts

            if url_terras in netrc_dict:
                print('\nRead login information for {} from netrc'.format(url_terras))
                userid = netrc_dict[url_terras][0]
                passwd = netrc_dict[url_terras][2]

        if userid is None:
            userid = input("User ID of terras: ")
        if passwd is None:
            passwd = getpass("Passwd of terras (not shown): ")

    except Usage as err:
        print("\nERROR:", file=sys.stderr, end='')
        print("  "+str(err.msg), file=sys.stderr)
        print("\nFor help, use -h or --help.\n", file=sys.stderr)
        return 2

    # %% Directory and file setting
    POSDIR = os.path.abspath(POSDIR)
    pos_dir = os.path.join(POSDIR, 'POS')
    poscor_dir = os.path.join(POSDIR, 'POSCOR')
    station_file = os.path.join(POSDIR, 'stations.txt')
    corr_file = os.path.join(POSDIR, 'corrf5o.dat')

    if not os.path.exists(POSDIR):
        os.mkdir(POSDIR)
    if not os.path.exists(pos_dir):
        os.mkdir(pos_dir)
    if not os.path.exists(poscor_dir):
        os.mkdir(poscor_dir)


    # %% Download pos files by ftp
    print("\nDownload pos files for each year from {}".format(url_terras))
    print("useid: {}".format(userid))
    if isdelta:
        print("\nDelta download:")
        print("If mdate of the remote file is older than the local file and "
              "the size is identical, skip downloading.\n")
    now_dt = dt.datetime.utcnow()

    ftp = FTP(url_terras)
    ftp.login(user=userid, passwd=passwd)

    for year in range(st_year, now_dt.year+1):
        print("Year: {}".format(year))

        pos_dir_y = os.path.join(POSDIR, str(year))
        if not os.path.exists(pos_dir_y):
            os.mkdir(pos_dir_y)

        ftp.cwd(os.path.join(dir_f5, str(year)))
        lines = []
        ftp.dir(lines.append)

        # Get modified date
        for line in lines:
            dir_str = line.split()
            file_name = dir_str[-1]
            out_pos = os.path.join(pos_dir_y, file_name)

            if isdelta and os.path.exists(out_pos): # Delta update if already exist
                local_dt = dt.datetime.fromtimestamp(Path(out_pos).stat().st_mtime)
                local_size = os.path.getsize(out_pos)

                # cf.https://dev.classmethod.jp/articles/python_ftp_modified_dt/
                if ":" in dir_str[7]:
                    remote_dt = dt.datetime.strptime("{} {}".format(
                        now_dt.year, " ".join(dir_str[5:8])), "%Y %b %d %H:%M")
                    if remote_dt > now_dt:
                        remote_dt = dt.strptime("{} {}".format(
                            now_dt.year - 1, " ".join(dir_str[5:8])),
                                                        "%Y %b %d %H:%M")
                else:
                    remote_dt = dt.datetime.strptime(" ".join(dir_str[5:8]),
                                                       "%b %d %Y")
                remote_size = int(dir_str[4])

                # Compare with local file time
                if remote_dt < local_dt and remote_size == local_size:
                    continue # no update

            # Download
            print("Download {}".format(file_name))
            isupdate = True
            with open(out_pos, 'wb') as fp:
                try:
                    ftp.retrbinary('RETR {}'.format(file_name), fp.write)
                except TimeoutError:
                    print('Connetction timed out. Connect again.')
                    # Log in and download again
                    ftp = FTP(url_terras)
                    ftp.login(user=userid, passwd=passwd)
                    ftp.cwd(os.path.join(dir_f5, str(year)))
                    ftp.retrbinary('RETR {}'.format(file_name), fp.write)

    ftp.close()


    # %% Check update
    if not isupdate: # No pos update. Finish
        print("\nNo pos files are updated. Finish")
        return 0


    # %% Get station list
    paths = list(Path(POSDIR).glob(r'[12][09][0-9][0-9]/*.pos'))
    stcodes = sorted(list(set([ path.name.split('.')[0] for path in paths ])))


    # %% Merge pos files
    print("")
    for stcode in stcodes:
        print("Merge {}".format(stcode))
        paths = list(Path(POSDIR).glob(
            r'[12][09][0-9][0-9]/{}.??.pos'.format(stcode)))
        paths = sorted(paths, key = lambda x: x.parent)

        merge_pos = os.path.join(pos_dir, '{}.pos'.format(stcode))
        if os.path.exists(merge_pos):
            os.remove(merge_pos)

        # Merge each year
        for path in paths:
            year = path.parent.name
            with open(path, 'r') as f:
                lines = f.readlines()
            lines_data = [line for line in lines
                          if line.startswith(' {}'.format(year))]
            with open(merge_pos, 'a') as f:
                f.writelines(lines_data)


    # %% Download and read corrf5o.dat
    urllib.request.urlretrieve(url_corr, corr_file)

    names = [ 'site', 'day', 'dx', 'dy', 'dz', 'db', 'dl', 'dh', 'comment' ]
    df_corr = pd.read_csv(corr_file, names=names, comment='#',
                          encoding='shift_jis')

    # %% Offset correction and make stations.txt
    names_pos = ['yyyy', 'mm', 'dd', 'time', 'x', 'y', 'z', 'lat', 'lon', 'h']
    if os.path.exists(station_file):
        os.remove(station_file)

    for stcode in stcodes:
        print("Offset correction for {}".format(stcode))
        merge_pos = os.path.join(pos_dir, '{}.pos'.format(stcode))
        cor_pos = os.path.join(poscor_dir, '{}.pos'.format(stcode))

        _df_corr = df_corr[df_corr['site']==stcode]
        if len(_df_corr) == 0: # No offset, just copy
            shutil.copyfile(merge_pos, cor_pos)
            continue

        df_pos = pd.read_table(merge_pos, sep='\s+', header=None,
                               names=names_pos, index_col=0,
                               parse_dates=[['yyyy','mm','dd']])

        for _corr in _df_corr.itertuples():
            if _corr.day == '0001/01/01': # Initial offset
                offset_dt = df_pos.index[0].to_pydatetime()
            else:
                offset_dt = dt.datetime.strptime(_corr.day, '%Y/%m/%d')

            # Add offset
            for i in range(6): #x,y,z,b,l,h
                df_pos.loc[offset_dt:, names_pos[i+4]] = \
                df_pos.loc[offset_dt:, names_pos[i+4]] + _corr[i+3]

        # Write corrected pos
        df_pos.to_csv(cor_pos, header=False, sep=' ', float_format='%.10e',
                      date_format='%Y %m %d', quotechar=' ')

        # Make stations.txt with corrdinates at first date
        with open(station_file, 'a') as f:
            print('{:>6s}  {:10.7f}  {:11.7f}  {:8.3f}'.format(stcode,
                   df_pos['lat'][0], df_pos['lon'][0], df_pos['h'][0]), file=f)


    # %% Finish
    elapsed_time = time.time()-start
    hour = int(elapsed_time/3600)
    minite = int(np.mod((elapsed_time/60), 60))
    sec = int(np.mod(elapsed_time, 60))
    print("\nElapsed time: {0:02}h {1:02}m {2:02}s".format(hour, minite, sec))

    print('\n{} Successfully finished!!\n'.format(os.path.basename(argv[0])))
    print('Output directory: {}\n'.format(POSDIR))


# %% main
if __name__ == "__main__":
    sys.exit(main())
