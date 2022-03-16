#!/usr/bin/env python3
"""
v1.0 20220316 Yu Morishita

This script filters GeoTIFF image.

=====
Usage
=====
filter_img.py -i in_GeoTIFF [-m method] [-s win_size] [-d deg]

 -i  Input GeoTIFF image
 -m  Filtering method (Default: g)
     g: 2D Gaussian filter
     m: 2D median filter
     r: Deramp
 -s  Window size in km for g and m methods (Default: 2 km)
 -d  Degree of deramp for r method (Default: 1)
     0:  a (Only bias)
     1:  a+bx+cy
     bl: a+bx+cy+dxy
     2:  a+bx+cy+dxy+ex**2+fy**2

Output: ${in_GeoTIFF%.tif}_${method}${win_size|deg}{LP,HP}.tif
        ${in_GeoTIFF%.tif}_${method}${win_size|deg}.png

"""
# %% Change log
'''
v1.0 20220316 Yu Morishita
 - Original implementation
'''


# %% Import




import getopt
import os
import sys
import numpy as np
import time
from osgeo import gdal, osr
from astropy.convolution import Gaussian2DKernel, convolve_fft
import GNSSlib
import filter_lib
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
    ver = '1.0'
    date = 20220316
    author = "Y. Morishita"
    prog_name = os.path.basename(argv[0])
    print(f"\n{prog_name} ver{ver} {date} {author}")
    print(f"{prog_name} {' '.join(argv[1:])}")

    # %% Set default
    geotiff_file = []
    method = 'g'
    win_size_km = 2
    deg = '1'

    compress_option = ['COMPRESS=DEFLATE', 'PREDICTOR=3']

    # %% Read options
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "hi:m:s:d:",
                                       ["help"])
        except getopt.error as msg:
            raise Usage(msg)
        for o, a in opts:
            if o == '-h' or o == '--help':
                print(__doc__)
                return 0
            elif o == '-i':
                geotiff_file = a
            elif o == '-m':
                method = a
            elif o == '-s':
                win_size_km = float(a)
            elif o == '-d':
                deg = a

        if not geotiff_file:
            raise Usage('No GeoTIFF given, -i is not optional!')
        elif not os.path.exists(geotiff_file):
            raise Usage(f'{geotiff_file} does not exist!')
        if not method in ['g', 'm', 'r']:
            raise Usage('method must be either g, m, or r!')
        if not deg in ['0', '1', 'bl', '2']:
            raise Usage('deg must be either 0, 1, bl, or 2!')

    except Usage as err:
        print("\nERROR:", file=sys.stderr, end='')
        print("  "+str(err.msg), file=sys.stderr)
        print("\nFor help, use -h or --help.\n", file=sys.stderr)
        return 2

    # %%  Settings
    if method == 'r':
        suffix = f'{method}{deg}'
    else:
        suffix = f'{method}{win_size_km}'

    geotiff_LP_file = geotiff_file.replace('.tif', f'_{suffix}_LP.tif')
    geotiff_HP_file = geotiff_file.replace('.tif', f'_{suffix}_HP.tif')
    png_file = geotiff_file.replace('.tif', f'_{suffix}.png')

    # %% Read GeoTIFF file
    geotiff = gdal.Open(geotiff_file)
    data = geotiff.ReadAsArray()
    width = geotiff.RasterXSize
    length = geotiff.RasterYSize
    lon_w_p, dlon, _, lat_n_p, _, dlat = geotiff.GetGeoTransform()
    # lat lon are in pixel registration. dlat is negative
    lon_w = lon_w_p + dlon/2  # grid registration
    lat_n = lat_n_p + dlat/2

    nodata = geotiff.GetRasterBand(1).GetNoDataValue() # None if blank
    data[data == nodata] = np.nan


    # %% Filter
    if method == 'r':
        print(f'\nDeramp with {deg} degree')
        data_LP, _ = filter_lib.fit2d(data, deg=deg)
        data_LP[np.isnan(data)] = np.nan

    elif method in ['g', 'm']:

        lat_c = lat_n + dlat*length/2
        lon_c = lon_w + dlon*width/2
        dx_grid_m, dy_grid_m, _ = GNSSlib.dBLH2dENU(lat_c, lon_c, 0,
                                                    np.abs(dlat), dlon, 0)
        print(f'\nX Grid size of GeoTIFF: {dx_grid_m:.1f} m ({dlon} deg)')
        print(f'Y Grid size of GeoTIFF: {dy_grid_m:.1f} m ({np.abs(dlat)} deg)')
        win_size_x_px = win_size_km/dx_grid_m*1000
        win_size_y_px = win_size_km/dy_grid_m*1000

        if method == 'g':
            print(f'\nGaussian filter with win_size of {win_size_km} km')
            print(f'X kernel width: {win_size_x_px:.1f} pixel')
            print(f'Y kernel width: {win_size_y_px:.1f} pixel')

            kernel = Gaussian2DKernel(win_size_x_px, win_size_y_px)
            data_LP = convolve_fft(data, kernel, preserve_nan=True,
                                   fill_value=np.nan, allow_huge=True)

        elif method == 'm':
            win_size_px = int(np.round((win_size_x_px+win_size_y_px)/2)/2)*2+1
            # must be odd integer and isotropic
            print(f'\nMedian filter with win_size of {win_size_km} km')
            print(f'Kernel width: {win_size_px} pixel')
            data_LP = filter_lib.nanmedian2d(data, win_sz=win_size_px)

    data_HP = data - data_LP


    # %% Output png
    data3 = [data, data_LP, data_HP]
    title3 = ['Original', f'LP ({suffix})', 'HP']
    GNSSlib.make_3im_png(data3, png_file, 'coolwarm', title3, vmin=None,
                         vmax=None, cbar=True)


    # %% Output GeoTIFF
    make_geotiff(data_LP, lat_n_p, lon_w_p, dlat, dlon, geotiff_LP_file,
                 compress_option, nodata=nodata)
    make_geotiff(data_HP, lat_n_p, lon_w_p, dlat, dlon, geotiff_HP_file,
                 compress_option, nodata=nodata)


    # %% Finish
    elapsed_time = time.time()-start
    hour = int(elapsed_time/3600)
    minite = int(np.mod((elapsed_time/60), 60))
    sec = int(np.mod(elapsed_time, 60))
    print("\nElapsed time: {0:02}h {1:02}m {2:02}s".format(hour, minite, sec))
    print('\n{} Successfully finished!!\n'.format(os.path.basename(argv[0])))
    print('Output: {}'.format(geotiff_LP_file))
    print('        {}'.format(geotiff_HP_file))
    print('        {}\n'.format(png_file), flush=True)

# %% make_geotiff
def make_geotiff(data, latn_p, lonw_p, dlat, dlon, outfile, compress_option, nodata=None):
    length, width = data.shape
    if data.dtype == np.float32:
        dtype = gdal.GDT_Float32
    elif data.dtype == np.uint8:
        dtype = gdal.GDT_Byte
    elif data.dtype == np.float64:
        data = np.float32(data)
        dtype = gdal.GDT_Float32

    driver = gdal.GetDriverByName('GTiff')
    outRaster = driver.Create(outfile, width, length, 1, dtype, options=compress_option)
    outRaster.SetGeoTransform((lonw_p, dlon, 0, latn_p, 0, dlat))
    outband = outRaster.GetRasterBand(1)
    outband.WriteArray(data)
    if nodata is not None: outband.SetNoDataValue(nodata)
    outRaster.SetMetadataItem('AREA_OR_POINT', 'Area')
    outRasterSRS = osr.SpatialReference()
    outRasterSRS.ImportFromEPSG(4326)
    outRaster.SetProjection(outRasterSRS.ExportToWkt())
    outband.FlushCache()



# %% main
if __name__ == "__main__":
    sys.exit(main())
