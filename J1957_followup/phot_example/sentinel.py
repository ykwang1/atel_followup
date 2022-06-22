#!/usr/bin/env python
#
# @file    sentinel.py
# @brief   Code to support a sentinel service to photometrically monitor sources
# @author  Matthew Graham, Przemek Mroz
#
# <!---------------------------------------------------------------------------
# Copyright (c) 2020 by the California Institute of Technology.
# This software is part of the ZTF software system.
# ---------------------------------------------------------------------------->

import numpy as np
import pandas as pd
import concurrent.futures
import json
import argparse
import os.path
# from sqlalchemy import create_engine, update
# import pymysql
import glob

from astropy.io import fits
from astropy import time
from time import gmtime, strftime

# from penquins import Kowalski
# from ztfutil import load_file
from calibrate_ZTF_phot import *

#------------------------------------------------------------------
# GLOBAL PARAMETERS
#------------------------------------------------------------------

ZTF_BASEURL = "https://irsa.ipac.caltech.edu/ibe/data/ztf/products"
DATA_DIR = '../images'
MATCH_RADIUS =  3 # 1.5
# ZTF_SOURCE_CATALOG = 'ZTF_sources_20200401'
FILTERS = {'g': 1, 'r': 2, 'i': 3}
# DB_ENGINE = create_engine('mysql+pymysql://alcuin:@localhost:3306/sentinel', echo = False)

# with open('passwd.json','r') as f:
#   passwd = json.load(f)

# KOW = Kowalski(username = passwd['kowalski']['username'], password = passwd['kowalski']['password'])

#------------------------------------------------------------------
# METHODS
#------------------------------------------------------------------


def process_photometry(id, ra, dec, frames, filters, field, ccd, quad):
    """ Process a set of photometry files

    Parameters:
    -----------
    id: str
        The ID of the source
    ra: float
        The RA of the source
    dec: float
        The Declination of the source
    frames: list  
        The photometry FITS files retrieved from IRSA
    filters: list
        A list of filters covered from ['g', 'r', 'i']
    field: int
        The field id
    ccd: int
        The ccd
    quad: int
        The quadrant
    """
    
    #----Read the calibration data-----#
    a = np.empty(2)
    b = np.empty(2)
    c = np.empty(2)
    field_cor = np.empty(2)
    coeff = np.empty((2, 441))
    x_knot = np.empty((2, 24))
    y_knot = np.empty((2, 24))

    rcid = 4 * (ccd - 1) + quad - 1
    deg = np.pi / 180.0
    for i in [ccd, quad]:
        print(type(i))
    for i, flt in enumerate(['g','r']):
        a[i], b[i], c[i] = get_magnitude_bias_coefficients(ccd, quad, flt)
        field_cor[i] = get_field_offset_corrections(field, rcid, flt)
        coeff[i], x_knot[i], y_knot[i] = \
            read_spline_correction_coefficients(ccd, quad, flt)

    #-----Extract and calibrate photometry-----#
    df = pd.DataFrame(columns = ['id', 'mjd', 'mag', 'magerr', 'fid', 'pid'])
        
    for frame in frames:

        # Opening FITS file
        filename = '%s/%s_c%02d_o_q%1d_psfcat.fits' % (DATA_DIR, frame, ccd, quad)
        if not os.path.exists(filename): continue
        
        with fits.open(filename) as hdul:
            data = hdul[1].data
            hdr = hdul[0].header

        magzp = hdr['MAGZP'] # magnitude zero point
        fid = hdr['FILTERID'] - 1 # filter ID
        programid = hdr['PROGRMID'] # program ID
        mjd = hdr['OBSMJD'] # + 2400000.5 # JD

        # Source matching
        f = (abs(data['ra'] - ra) < 0.002) & (abs(data['dec'] - dec) < 0.002)
        if sum(f) == 0: continue

        adist = ((data['ra'][f] - ra) * np.cos(dec * deg)) ** 2 
        adist += (data['dec'][f] - dec) ** 2
        adist = np.sqrt(adist) * 3600.0 # angular distance in arcsec
        idx = np.argmin(adist)

        if adist[idx] > MATCH_RADIUS: continue

        x = data['xpos'][f][idx]
        y = data['ypos'][f][idx]
        mag = data['mag'][f][idx] + magzp
        sigmag = data['sigmag'][f][idx]
        flg = data['flags'][f][idx]
        
        mag_cal = calibrate(mag, x, y, a[fid], b[fid], c[fid], field_cor[fid],\
            coeff[fid], x_knot[fid], y_knot[fid])

        df = df.append({'id': id, 'mjd': mjd, 'mag': mag_cal, 'magerr': sigmag, 'fid': fid, 'pid': programid}, ignore_index = True)
        df.to_csv('photdata.csv', index=False)

#     df.to_sql('photdata', DB_ENGINE, index = False, if_exists = 'append', schema = 'sentinel')

#     conn = DB_ENGINE.connect()
#     sql = f'update watchlist set last_updated = {get_time()}'
#     res = conn.execute(sql)
    

    

def test():
    """ Test the sentinel setup and update
    
    Parameters:
    -----------
    """
#     load_watchlist('binary.csv')
#     update()
    frames = glob.glob(f'{DATA_DIR}/*.fits')
#     frames = [f[7:-21] for f in frames]
    process_photometry(1, 299.16300398, 3.4454103, frames, ['g', 'r'], 1483, 13, 2)
if __name__ == "__main__":
    test()
    

