"""
Calibrations of the ZTF PSF raw photometry

P. Mroz @ Caltech, 29 Apr 2020
"""

import numpy as np
import pandas as pd
import warnings

MAGNITUDE_BIAS_COEFFICIENTS = 'calib_data/Mag_Bias.txt'
FIELD_OFFSET_CORRECTIONS = 'calib_data/Field_Corrections.txt'
SPLINE_CORRECTION_DIR = 'calib_data'

def read_magnitude_bias_coefficients():
    """Load the magnitude bias fit coefficients from file Mag_Bias.txt
    """
    names = ['ccd','quad','filter','a','b','c']
    try:
#         d = np.genfromtxt(MAGNITUDE_BIAS_COEFFICIENTS,unpack=True,
#             names=names,dtype=None,encoding='utf-8')
        d = pd.read_csv(MAGNITUDE_BIAS_COEFFICIENTS, delim_whitespace=True,
            names=names,dtype=None,encoding='utf-8', comment='#')
    except IOError:
        exit('Error: Cannot open file %s'%MAGNITUDE_BIAS_COEFFICIENTS)
    return d

def get_magnitude_bias_coefficients(ccd, quad, filter):
    """Get the magnitude bias fit coefficients for a given ccd, quadrant
    and filter
    
    Parameters
    ----------
    ccd : int
        ZTF CCD ID
    quad : int
        Quadrant number
    filter : str
        The 1-letter filter code ('g', 'r', 'i') 
    """
    d = read_magnitude_bias_coefficients()
    f = (d['ccd'] == ccd) & (d['quad'] == quad) & (d['filter'] == filter)
    if sum(f) != 1: 
        exit('Error while loading magnitude bias coefficients.')
    return d['a'][f].values[0], d['b'][f].values[0], d['c'][f].values[0]
    
def read_field_offset_corrections():
    """Load the field offset corrections from file Field_Corrections.txt
    """
    names = ['field','rc','filter','correction']
    try:
#         d = np.genfromtxt(FIELD_OFFSET_CORRECTIONS,unpack=True,
#             names=names,dtype=None,encoding='utf-8')
        d = pd.read_csv(FIELD_OFFSET_CORRECTIONS, delim_whitespace=True,
            names=names,dtype=None,encoding='utf-8', comment='#')
    except IOError:
        exit('Error: Cannot open file %s'%FIELD_OFFSET_CORRECTIONS)
    return d
    
def get_field_offset_corrections(field, rc, filter):
    """Get the field offset correction for a given field, RC, and filter
    
    Parameters
    ----------
    field : int
        ZTF field number
    RC : int
        ZTF RC number
    filter : str
        The 1-letter filter code ('g', 'r', 'i') 
    """
    d = read_field_offset_corrections()
    f = (d['field'] == field) & (d['rc'] == rc) & (d['filter'] == filter)
    if sum(f) > 1:
        exit('Error while loading field offset corrections.')
    elif sum(f) == 0:
        warnings.warn('No field offset correction for field %d'%field)
        return 0.0
    return d['correction'][f].values[0]
    
def read_spline_correction_coefficients(ccd, quad, filter):
    """Load the spline fit correction coefficients from a text file.
    
    Parameters
    ----------
    ccd : int
        ZTF CCD number
    quad : int
        quadrant number
    filter : str
        The 1-letter filter code ('g', 'r', 'i') 
    """
    coeff, x_knot, y_knot = [], [], []
    fname = '%s/ztf_z%1s_c%02d_q%1d.out.txt.spline'%\
        (SPLINE_CORRECTION_DIR,filter,ccd,quad)
    with open(fname,'r') as fp:
        for line in fp:
            tmp = line.split()
            if len(tmp) == 0: continue
            if len(tmp) == 0: continue
            if tmp[0] == 'coeff': coeff.append(float(tmp[3]))
            elif tmp[0] == 'X': x_knot.append(float(tmp[5]))
            elif tmp[0] == 'Y': y_knot.append(float(tmp[5]))
    return np.array(coeff), np.array(x_knot), np.array(y_knot)

def calibrate(magztf, x, y, a, b, c, field_cor, coeff, x_knot, y_knot):
    """Calibrate a single photometric data point
    
    Parameters
    ----------
    magztf : float
        uncalibrated ZTF PSF magnitude
    x : float
        x coordinate
    y : float
        y coordinate
    a : float
        magnitude bias fit coefficient a
    b : float
        magnitude bias fit coefficient b
    c : float
        magnitude bias fit coefficient c
    field_cor : float
        field offset correction
    coeff : array
        spline fit coefficients
    x_knot : array
        spline fit x knots
    y_knot : array
        spline fit y knots
    """
    if magztf < 19.25:
        magbias_cor = a + b*magztf + c*magztf*magztf
    else:
        magbias_cor = a + b*19.25 + c*19.25*19.25
    spline_cor = spline_model(x,y,x_knot,y_knot,coeff)
    mag_calib = magztf + field_cor + magbias_cor + spline_cor
    return mag_calib

def calibrate_all(magztf, x, y, a, b, c, field_cor, coeff, x_knot, y_knot):
    """Calibrate an array of photometric data points
    
    Parameters
    ----------
    magztf : np.array
        uncalibrated ZTF PSF magnitudes
    x : np.array
        x coordinates
    y : np.array
        y coordinates
    a : float
        magnitude bias fit coefficient a
    b : float
        magnitude bias fit coefficient b
    c : float
        magnitude bias fit coefficient c
    field_cor : float
        field offset correction
    coeff : array
        spline fit coefficients
    x_knot : array
        spline fit x knots
    y_knot : array
        spline fit y knots
    """
    magbias_cor = np.zeros(len(magztf))
    f = (magztf < 19.25)
    magbias_cor[f] = a + b*magztf[f] + c*magztf[f]*magztf[f]
    magbias_cor[~f] = a + b*19.25 + c*19.25*19.25
    spline_cor = np.array([spline_model(_x,_y,x_knot,y_knot,coeff) \
        for _x,_y in zip(x,y)])
    mag_calib = magztf + field_cor + magbias_cor + spline_cor
    print(field_cor + magbias_cor + spline_cor)
    return mag_calib

def spline_model(x_in, y_in, tx, ty, coeff):
    """Spline model (by Andrew Drake)
    
    Parameters
    ----------
    x_in : float
        x coordinate
    y_in : float
        y coordinate
    coeff : array
        spline fit coefficients
    x_knot : array
        spline fit x knots
    y_knot : array
        spline fit y knots
    """
    nx = 24
    ny = 24
    kx = 2
    ky = 2

    h = [0.0] * 25
    hh = [0.0] * 25
    w_x = [0.0] * 25
    w_y = [0.0] * 25

    kx1 = kx+1
    nkx1 = nx-kx1
    l = kx1
    l1 = l+1

    while x_in >= tx[l1-1] and l != nkx1:
        l = l1
        l1 = l+1
        
    h[0] = 1.0
    for j in range(1, kx+1):
        for i in range(j):
            hh[i] = h[i]
        h[0] = 0.0
        for i in range(j):
            li = l+i
            lj = li-j
            if tx[li] != tx[lj]:
                f = hh[i] / (tx[li] - tx[lj])
                h[i] = h[i] + f * (tx[li] - x_in)
                h[i+1] = f * (x_in - tx[lj])
            else:
                h[i+1-1] = 0.0
                
    lx = l-kx1
    for j in range(kx1):
        w_x[j] = h[j]

    ky1 = ky+1
    nky1 = ny-ky1
    l = ky1
    l1 = l+1

    while y_in >= ty[l1-1] and l != nky1:
        l = l1
        l1 = l+1
        
    h[0] = 1.0
    for j in range(1, ky+1):
        for i in range(j):
            hh[i] = h[i]
        h[0] = 0.0
        for i in range(j):
            li = l+i
            lj = li-j
            if ty[li] != ty[lj]:
                f = hh[i] / (ty[li] - ty[lj])
                h[i] = h[i] + f * (ty[li] - y_in)
                h[i+1] = f * (y_in - ty[lj])
            else:
                h[i+1-1] = 0.0

    ly = l-ky1
    for j in range(ky1):
        w_y[j] = h[j]

    l = lx*nky1
    for i1 in range(kx1):
        h[i1] = w_x[i1]
        
    l1 = l+ly
    temp = 0.0
    for i1 in range(kx1):
        l2 = l1
        for j1 in range(ky1):
            l2 = l2+1
            temp = temp + coeff[l2-1] * h[i1] * w_y[j1]
        l1 = l1+nky1
        
    return temp
