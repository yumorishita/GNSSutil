#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions for GNSSutil.

Change log:
v1.1.0 20220315 Yu Morishita
- # as comment in read_pos and read_stations
- Add make_3im_png
v1.0.0 20220215 Yu Morishita
- Original implementation

"""
'''
To do:
- Estimate uncertainty in calc_vel
  (cf. https://www1.kaiho.mlit.go.jp/GIJUTSUKOKUSAI/KENKYU/report/rhr53/rhr53-TR04.pdf)

'''

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt


# %% BLH2XYZ
def BLH2XYZ(lat, lon, h):
    """
    Convert BLH to XYZ (ECEF) coordinates.

    Parameters
    ----------
    lat, lon : float
        Latidude and longitude in degree.
    h : float
        Height in meter.

    Returns
    -------
    X, Y, Z : float
        XYZ coordinates in meter.

    See Also
    --------
    cf. https://www.enri.go.jp/~fks442/K_MUSEN/1st/1st060428rev2.pdf

    """
    a = 6378137.0
    f = 1/298.257222101  # GRS80
#    f = 1/298.257223563 # WGS84
    e2 = f*(2-f)
    N = a/np.sqrt(1-e2*np.sin(np.deg2rad(lat))**2)

    X = (N+h)*np.cos(np.deg2rad(lat))*np.cos(np.deg2rad(lon))
    Y = (N+h)*np.cos(np.deg2rad(lat))*np.sin(np.deg2rad(lon))
    Z = (N*(1-e2)+h)*np.sin(np.deg2rad(lat))

    return X, Y, Z


# %% calc_displacement
def calc_displacement(df_pos, s_date1, s_date2, e_date1, e_date2):
    """
    Calculate displacements in ENU directions.

    Parameters
    ----------
    df_pos : DataFrame
        Pos data with 'yyyy_mm_dd' (index), 'lat', 'lon', 'h' fields.
    s_date1, s_date2 : datetime, 'yyyymmdd' or any date format
        Start date period for displacement calculation.
    e_date1, e_date2 : datetime, 'yyyymmdd' or any date format
        End date period for displacement calculation.

    Returns
    -------
    dENU : array (size: (3,))
        Displacements in ENU directions in meter.
    dENUstd : array (size: (3,))
        Standard deviations of the displacements in ENU directions in meter.
    s_n_day, e_n_day : int
        Number of used date for each period

    """
    df_pos_s = df_pos.loc[s_date1:s_date2]
    df_pos_e = df_pos.loc[e_date1:e_date2]
    s_n_day = len(df_pos_s)
    e_n_day = len(df_pos_e)

    lat_s_mean = df_pos_s['lat'].mean()
    lat_s_std = df_pos_s['lat'].std()
    lat_e_mean = df_pos_e['lat'].mean()
    lat_e_std = df_pos_e['lat'].std()
    lon_s_mean = df_pos_s['lon'].mean()
    lon_s_std = df_pos_s['lon'].std()
    lon_e_mean = df_pos_e['lon'].mean()
    lon_e_std = df_pos_e['lon'].std()
    h_s_mean = df_pos_s['h'].mean()
    h_s_std = df_pos_s['h'].std()
    h_e_mean = df_pos_e['h'].mean()
    h_e_std = df_pos_e['h'].std()

    dlat = lat_e_mean - lat_s_mean
    dlat_std = np.sqrt(lat_e_std**2+lat_s_std**2)
    dlon = lon_e_mean - lon_s_mean
    dlon_std = np.sqrt(lon_e_std**2+lon_s_std**2)
    dh = h_e_mean - h_s_mean
    dh_std = np.sqrt(h_e_std**2+h_s_std**2)

    dE, dN, dU = dBLH2dENU(lat_s_mean, lon_s_mean, h_s_mean, dlat, dlon, dh)
    dEstd, dNstd, dUstd = dBLH2dENU(lat_s_mean, lon_s_mean, h_s_mean,
                                    dlat_std, dlon_std, dh_std)

    dENU = np.array([dE, dN, dU])
    dENUstd = np.array([dEstd, dNstd, dUstd])

    return dENU, dENUstd, s_n_day, e_n_day


# %% calc_vel
def calc_vel(df_pos, s_date=None, e_date=None):
    """
    Calculate velocities in ENU directions.

    Parameters
    ----------
    df_pos : DataFrame
        Pos data with 'yyyy_mm_dd' (index), 'lat', 'lon', 'h' fields.
    s_date, e_date (optional) : datetime, 'yyyymmdd' or any date format
        Start and end date for velocity calculation. The default is None.

    Returns
    -------
    vel_ENU : array (size: (3,))
        Velocities in ENU directions in m/yr.
    s_date, e_date : datetime.date
        Start and end date used for velocity calculation.
    n_day : int
        Number of days used for velocity calculation.

    """
    df_pos = df_pos.loc[s_date:e_date]
    dday = (df_pos.index-df_pos.index[0]).days.to_numpy()
    n_day = len(dday)
    s_date = df_pos.index[0].to_pydatetime().date()
    e_date = df_pos.index[-1].to_pydatetime().date()

    dE, dN, dU = dBLH2dENU(df_pos['lat'][0], df_pos['lon'][0], df_pos['h'][0],
                           (df_pos['lat']-df_pos['lat'][0]).to_numpy(),
                           (df_pos['lon']-df_pos['lon'][0]).to_numpy(),
                           (df_pos['h']-df_pos['h'][0]).to_numpy())

    G = np.stack((np.ones_like(dday), dday), axis=1)

    vel_ENU = np.linalg.lstsq(G, np.array([dE, dN, dU]).T, rcond=None)[0][1]
    vel_ENU = vel_ENU*365.24 # m/day -> m/yr

    return vel_ENU, s_date, e_date, n_day


# %% dBLH2dENU
def dBLH2dENU(lat, lon, h, dlat, dlon, dh):
    """
    Convert BLH dispalcement to ENU displacement

    Parameters
    ----------
    lat, lon : float
        Latitude and longitude in degree.
    h : float
        Height in meter.
    dlat, dlon : float or 1-D array
        Displacement along latitude and longitude in degree.
    dh : float or 1-D array
        Displacement of height in meter.

    Returns
    -------
    dE, dN, dU : float or 1-D array
        Eastward, northward, and upward displacement in meter.

    """
    dX, dY, dZ = dBLH2dXYZ(lat, lon, h, dlat, dlon, dh)
    dE, dN, dU = dXYZ2dENU(lat, lon, dX, dY, dZ)

    return dE, dN, dU


# %% dBLH2dXYZ
def dBLH2dXYZ(lat, lon, h, dlat, dlon, dh):
    """
    Convert BLH dispalcement to XYZ displacement

    Parameters
    ----------
    lat, lon : float
        Latitude and longitude in degree.
    h : float
        Height in meter.
    dlat, dlon : float or 1-D array
        Displacement along latitude and longitude in degree.
    dh : float or 1-D array
        Displacement of height in meter.

    Returns
    -------
    dX, dY, dZ : float or 1-D array
        XYZ (ECEF) displacement in meter.

    """

    X0, Y0, Z0 = BLH2XYZ(lat, lon, h)
    X1, Y1, Z1 = BLH2XYZ(lat+dlat, lon+dlon, h+dh)
    dX = X1 - X0
    dY = Y1 - Y0
    dZ = Z1 - Z0

    return dX, dY, dZ

# %% dENU2dXYZ
def dENU2dXYZ(lat, lon, dE, dN, dU):
    """
    Convert ENU displacement to XYZ displacement

    Parameters
    ----------
    lat, lon : float
        Latitude and longitude in degree.
    dE, dN, dU : float or 1-D array
        Eastward, northward, and upward displacement in meter.

    Returns
    -------
    dX, dY, dZ : float or 1-D array
        XYZ (ECEF) displacement in meter.

    See Also
    --------
    cf. http://sapt.sakura.ne.jp/saptweb/wp-content/uploads/2017/01/enu2xyz.pdf

    """

    Brad = np.deg2rad(lat)
    Lrad = np.deg2rad(lon)

    A = np.zeros((3, 3), dtype=np.float64)
    A[0, 0] = -np.sin(Lrad)
    A[0, 1] = -np.sin(Brad)*np.cos(Lrad)
    A[0, 2] = np.cos(Brad)*np.cos(Lrad)
    A[1, 0] = np.cos(Lrad)
    A[1, 1] = -np.sin(Brad)*np.sin(Lrad)
    A[1, 2] = np.cos(Brad)*np.sin(Lrad)
    A[2, 1] = np.cos(Brad)
    A[2, 2] = np.sin(Brad)

    dX, dY, dZ = np.matmul(A, np.array([dE, dN, dU]))

    return dX, dY, dZ


# %% dENU_sub_ref
def dENU_sub_ref(lat, lon, dE, dN, dU, lat_ref, lon_ref, dE_ref, dN_ref,
                 dU_ref):
    """
    Subtract reference ENU displacement

    Parameters
    ----------
    lat, lon : float
        Latitude and longitude in degree.
    dE, dN, dU : float or 1-D array
        Eastward, northward, and upward displacement in meter.
    lat_ref, lon_ref : float
        Latitude and longitude at the reference in degree.
    dE_ref, dN_ref, dU_ref : float or 1-D array
        Eastward, northward, and upward displacement at the reference in meter.

    Returns
    -------
    dErel, dNrel, dUrel : float or 1-D array
        Relative ENU displacement to the reference in meter.

    """

    dX, dY, dZ = dENU2dXYZ(lat, lon, dE, dN, dU)
    dX_ref, dY_ref, dZ_ref = dENU2dXYZ(lat_ref, lon_ref, dE_ref, dN_ref, dU_ref)
    dXrel = dX - dX_ref
    dYrel = dY - dY_ref
    dZrel = dZ - dZ_ref
    dErel, dNrel, dUrel = dXYZ2dENU(lat, lon, dXrel, dYrel, dZrel)

    return dErel, dNrel, dUrel


# %% dXYZ2dENU
def dXYZ2dENU(lat, lon, dX, dY, dZ):
    """
    Convert XYZ displacement to ENU displacement

    Parameters
    ----------
    lat, lon : float
        Latitude and longitude in degree.
    dX, dY, dZ : float or 1-D array
        XYZ (ECEF) displacement in meter.

    Returns
    -------
    dE, dN, dU : float or 1-D array
        Eastward, northward, and upward displacement in mm.

    See Also
    --------
    cf. http://memo--randum.blogspot.com/2010/06/gps.html

    """

    Brad = np.deg2rad(lat)
    Lrad = np.deg2rad(lon)

    A = np.zeros((3, 3), dtype=np.float64)
    A[0, 0] = -np.sin(Lrad)
    A[0, 1] = np.cos(Lrad)
    A[1, 0] = -np.sin(Brad)*np.cos(Lrad)
    A[1, 1] = -np.sin(Brad)*np.sin(Lrad)
    A[1, 2] = np.cos(Brad)
    A[2, 0] = np.cos(Brad)*np.cos(Lrad)
    A[2, 1] = np.cos(Brad)*np.sin(Lrad)
    A[2, 2] = np.sin(Brad)

    dE, dN, dU = np.matmul(A, np.array([dX, dY, dZ]))

    return dE, dN, dU


#%% make_3im_png
def make_3im_png(data3, pngfile, cmap, title3, vmin=None, vmax=None, cbar=True):
    """
    Make png with 3 images for comparison.
    data3 and title3 must be list with 3 elements.
    """
    ### Plot setting
    interp = 'nearest'

    length, width = data3[0].shape
    figsizex = 12
    xmergin = 4 if cbar else 0
    figsizey = int((figsizex-xmergin)/3*length/width)+2

    fig = plt.figure(figsize = (figsizex, figsizey))

    for i in range(3):
        ax = fig.add_subplot(1, 3, i+1) #index start from 1
        im = ax.imshow(data3[i], vmin=vmin, vmax=vmax, cmap=cmap,
                       interpolation=interp)
        ax.set_title(title3[i])
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        if cbar: fig.colorbar(im, ax=ax)

    plt.tight_layout()
    plt.savefig(pngfile)
    plt.close()

    return


# %% read_dENU
def read_dENU(displ_file):
    """
    Read displacement txt file.

    Parameters
    ----------
    displ_file : str
        Path to a displacement file.
        Format: stcode  lat(deg)  lon(deg)  dE  dN  dU

    Returns
    -------
    df_dENU : DataFrame
        dENU data with 'stcode' (index), 'lat', 'lon', 'dE', 'dN', 'dU' columns.

    """
    names = ['stcode', 'lat', 'lon', 'dE', 'dN', 'dU']
    df_dENU = pd.read_table(displ_file, sep='\s+', header=None,
                               names=names, index_col=0, comment='#')

    return df_dENU


# %% read_pos
def read_pos(pos_file, s_date=None, e_date=None):
    """
    Read *.pos file.

    Parameters
    ----------
    pos_file : str
        Path to a *.pos file.
    s_date (optional) : datetime, 'yyyymmdd' or any date format
        Start date to read. The default is None.
    e_date (optional) : datetime, 'yyyymmdd' or any date format, optional
        End date to read. The default is None.

    Returns
    -------
    df_pos : DataFrame
        Pos data with 'yyyy_mm_dd' (index), 'lat', 'lon', 'h' columns.

    """
    names = ['yyyy', 'mm', 'dd', 'time', 'x', 'y', 'z', 'lat', 'lon', 'h']
    usecols = ['yyyy', 'mm', 'dd', 'lat', 'lon', 'h']
    df_pos = pd.read_table(pos_file, sep='\s+', header=None, names=names,
                           index_col=0, usecols=usecols, comment='#',
                           parse_dates=[['yyyy','mm','dd']])
    df_pos = df_pos.loc[s_date:e_date]

    return df_pos


# %% read_stations
def read_stations(station_file, lat_s=-90, lat_n=90, lon_w=-180, lon_e=180):
    """
    Read station file.

    Parameters
    ----------
    station_file : str
        Path to a station file.
        Format: stcode  lat(deg)  lon(deg)  h(m)
    lat_s, lat_n, lon_w, lon_e (optional) : float
        Area to read. The default is all.

    Returns
    -------
    df_station : DataFrame
        Station data with 'stcode' (index), 'lat', 'lon', 'h' columns.

    """
    names = ['stcode', 'lat', 'lon', 'h']
    df_station = pd.read_table(station_file, sep='\s+', header=None,
                               names=names, index_col=0, comment='#')

    df_station = df_station.query('{}<=lon<={}'.format(lon_w, lon_e)).query(
        '{}<=lat<={}'.format(lat_s, lat_n))

    return df_station


# %% XYZ2BLH
def XYZ2BLH(X, Y, Z):
    """
    Convert XYZ (ECEF) coordinates to BLH.

    Parameters
    ----------
    X, Y, Z : float or 1-D array
        XYZ coordinates in meter.

    Returns
    -------
    lat, lon : float or 1-D array
        Latidude and longitude in degree.
    h : float or 1-D array
        Height in meter.

    See Also
    --------
    cf. https://www.enri.go.jp/~fks442/K_MUSEN/1st/1st060428rev2.pdf

    """
    a = 6378137.0
    f = 1/298.257222101  # GRS80
#    f = 1/298.257223563 # WGS84
    b = a*(1-f)
    p = np.sqrt(X*X+Y*Y)
    theta = np.arctan2(Z*a, p*b)
    e2 = f*(2-f)
    ep2 = e2*a**2/b**2

    lat_rad = np.arctan2(Z+ep2*b*np.sin(theta)**3,
                         p-e2*a*np.cos(theta)**3)
    lat = np.rad2deg(lat_rad)
    lon = np.rad2deg(np.arctan2(Y, X))
    h = p/np.cos(lat_rad)-a/np.sqrt(1-e2*np.sin(lat_rad)**2)

    return lat, lon, h
