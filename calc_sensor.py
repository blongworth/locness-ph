#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions to covert raw BGC sensor output to scientific units


@author: dnicholson
"""
import numpy as np
from numpy.polynomial.polynomial import polyval #from low to high

# ************************************************************************
#  SET SOME CONSTANTS
# ************************************************************************
K0 = 273.15 #kelvin
R = 8.3144626 # Gas constant J/(mol K)
F = 96485.33  # Faraday constant Coulomb / mol

# O2 solubility coefficients from Garcia and Gordon (1992) for Kraus and Benson
# fit.
VO2 = 22.3916 # O2 molar volume at STP
o2coef = {
    'a': (2.00907,  3.22014,  4.0501,  4.94457, -0.256847,  3.88767), #what are these?
    'b': (-6.24523e-3, -7.37614e-3, -1.03410e-2, -8.17083e-3),
    'c0':  -4.88682e-7,
    }


def o2_sbe(cal,phase,otv,pres,tem,sal,sref=0,B15=True):
    """
    Parameters
    ----------
    cal :
        ToDo
    phase : numpy array
        phase measurements from sensor
    otv : numpy array
        temperature voltage from sensor
    pres : numpy array
        seawater pressure (dbar) from ctd
    tem : numpy array
        pressure from ctd #temp??
    sal : numpy array
        salinity from ctd
    sref : scalar, optional
        reference salinity setting for sensor. The default is 0.
    B15 : boolean, optional
        Flag to use Bittig 2015 pressure dependence. The default is True.

    Returns
    -------
    O2 : numpy array
        Oxygen concentration in uM
    t_o : numpy array
        Oxygen sensor temperature (deg C)

    References
    ----------
    Bittig et al. 2015:  https://doi.org/10.1175/JTECH-D-15-0108.1
    Garcia and Gordon 1992: https://doi.org/10.4319/lo.1992.37.6.1307
    Benson and Krause 1984: https://doi.org/10.4319/lo.1984.29.3.0620

    Example
    -------
    cal = ToDo

    """
    #Convert cal data into arrays
    cal_a = np.array([cal['A0'], cal['A1'], cal['A2']])
    cal_b = np.array([cal['B0'], cal['B1']])
    cal_c = np.array([cal['C0'], cal['C1'], cal['C2']])
    cal_ta = np.array([cal['TA0'], cal['TA1'], cal['TA2'], cal['TA3']])

    dPhdV = 39.4570707

    L = np.log(1e5*otv / (3.3-otv))
    t_o =  1 / polyval(L,cal_ta ) - K0 #oxygen sensor temp
    t_s = np.log((25+K0 - t_o)/(K0 + t_o)) #Different from mbari and aoml

    # Bittig pressure correction (2015) Part 1, T & O2 independent
    if B15:
        phase = phase + (0.115 * pres / 1000)
        Pcorr = (0.00022* tem + 0.0419) * pres / 1000 + 1
    else:
        Pcorr = np.exp(0.011 * np.maximum(pres,0) / (tem + K0))

    #Salt correction
    Scorr = np.exp((sal-sref) * polyval(t_s, o2coef['b']) + o2coef['c0'] * (sal**2 - sref**2))
    V = phase / dPhdV
    Ksv = polyval(tem, cal_c)
    O2_mlL = (((cal['A0'] + tem * cal['A1'] +  cal['A2']*V**2) /
               (cal['B0'] + V * cal['B1']) - 1.0) / Ksv) * Scorr * Pcorr
    O2 = (1000 / VO2) * O2_mlL #converts from ml/L to uM

    return (O2,t_o,O2_mlL)

# def o2_aa_svu(cal,phase,pres,tem,sal,sref=0,B15=True):
#     if 'b' not in cal:
#         cal['b'] = (-4.29155e-3, -6.90358e-3, -6.93498e-3, -6.24097e-3)
#         cal['c0'] = -3.11680e-7

#     t_s = np.log((25+K0 - tem)/(K0 + tem))
#     # pressure correction
#     if B15:
#         phase    = phase + (0.115 * pres / 1000)
#         Pcorr = (0.00022* tem + 0.0419) * pres / 1000 + 1
#     else:
#         Pcorr = np.exp(0.011 * np.maximum(pres,0) / (tem + K0))
#     # salinity correction
#     Scorr = np.exp((sal-sref) * polyval(t_s,cal['b']) +
#                    cal['c0'] * (sal**2 - sref**2))
#     Ksv = polyval(tem,cal['svu'][:3])
#     P0 = polyval(tem,cal['svu'][3:5])
#     PC = polyval(phase,cal['svu'][5:])
#     optode_uM = ((P0 / PC)-1) / Ksv
#     O2 = optode_uM * Scorr * Pcorr
#     if 'conc' in cal:
#         O2 = polyval(O2,cal['conc'])
#     return (O2)

# def o2sol(tem,sal):
#     """
#     Calculate equilibrium solubility of dissolved oxygen

#     Parameters
#     ----------
#     tem : array_like
#         temperature [deg C]
#     sal : array_like
#         salinity from ctd [pss]

#     Returns
#     -------
#     O2sol
#         Equilibrium solubility (1atm, saturated water vapor) of oxygen [uM]

#     References
#     ----------
#     Garcia and Gordon 1992: https://doi.org/10.4319/lo.1992.37.6.1307
#     """
#     t_s = np.log((25+K0 - tem)/(K0 + tem))
#     lnC = (polyval(t_s,o2coef['a'])+ sal * (polyval(t_s, o2coef['b']) +
#            sal * o2coef['c0']))
#     O2sol = 1000 * np.exp(lnC) / VO2
#     return O2sol

def mcoms(scale,dark,counts):
    out = scale * (counts-dark)
    return out

def ocr(cal,ch1,ch2,ch3,ch4):
    # DOWN_IRRADIANCE380 = 0.01*A1_380*(RAW_DOWNWELLING_IRRADIANCE380 -
    # A0_380) * lm_380
    # ToDo, needs some work

    I1 = 0.01 * cal['a1'][0] * (ch1 - cal['a0'][0]) * cal['Im'][0]
    I2 = 0.01 * cal['a1'][1] * (ch1 - cal['a0'][1]) * cal['Im'][1]
    I3 = 0.01 * cal['a1'][2] * (ch1 - cal['a0'][2]) * cal['Im'][2]
    I4 = 0.01 * cal['a1'][3] * (ch1 - cal['a0'][3]) * cal['Im'][3]
    return (I1,I2,I3,I4)

def pH_sbe(cal,Vrs,pres,tem,sal):
    """
    Calculate pH from raw Seabird sensor output

    Parameters
    ----------
    cal : dict
        calibration coefficients including the following fields
            k0      = Sensor reference potential (intercept at tem = 0C)
            k2      = linear temerature coefficient (slope)
            fp1, fp2, etc.  = sensor dependent pressure coefficients

    Vrs : array_like
        Voltage bewteen reference electrode and ISFET source [V]
    pres : TYPE
        Pressures [decibars]
    tem : TYPE
        Temperature [deg C]
    sal : TYPE
        Salinity (usually CTD salinity) [PSS]

    Returns
    -------
    phfree : TYPE
        mol/kg-seawater scale
    phtot : TYPE
        total proton scale

    """


    Tk   = K0 + tem # degrees Kelvin
    ln10 = np.log(10) # natural log of 10
    pres_bar = pres / 10
    # some newer sensors arriving with 9th order polynomials
    if "F9" in cal.keys():
        Pcoefs = np.array([0,cal['F1'],cal['F2'],cal['F3'],cal['F4'],cal['F5'],cal['F6'],cal['F7'],cal['F8'],cal['F9']]) #Cal coefs to array
    else:   
        Pcoefs = np.array([0,cal['F1'],cal['F2'],cal['F3'],cal['F4'],cal['F5'],cal['F6']]) #Cal coefs to array

    # ************************************************************************
    # CALCULATE PHYSICAL AND THERMODYNAMIC DATA
    # Dickson, A. G., Sabine, C. L., & Christian, J. R. (2007). Guide to best
    # practices for ocean CO2 measurements.
    # ************************************************************************

    # IONIC STRENGTH OF SEAWATER (mol / kg H2O)
    # Varified units by comparing to Dickson et al. 2007: Chap 5, p10 Table 2
    # Dickson et al. 2007: Chap 5, p13 Eq 34
    IonS = 19.924 * sal / (1000 - 1.005 * sal)

    # MEAN SEAWATER SULFATE CONCENTRATION (mol / kg solution)
    # This wants to be mol/kg-seawater  as KHSO4 is on that scale
    # Dickson et al. 2007: Chap 5, p10 Table 2
    Stotal = (0.14 / 96.062) * (sal / 1.80655)

    # MEAN SEAWATER CHLORIDE CONCENTRATION  (mol / kg H20)
    # this wants to be mol/kg H2O as activity is on mol/kg H2O scale
    # Dickson et al. 2007: Chap 5, p10 Table 2
    Cltotal = 0.99889 / 35.453 * sal / 1.80655 #(mol / kg solution)
    Cltotal = Cltotal /(1 - 0.001005 * sal)  # (mol / kg H20)

    # BISULFIDE DISSCIATION CONSTANT AT T,S AND IONIC STRENGTH(mol/kg solution)
    # Dickson et al. 2007: Chap 5, p12 Eq 33
    Khso4 = np.exp(-4276.1 / Tk + 141.328 - 23.093 * np.log(Tk) + \
        (-13856 / Tk + 324.57 - 47.986 * np.log(Tk)) * IonS ** 0.5 + \
        (35474 / Tk - 771.54 + 114.723 * np.log(Tk)) * IonS - \
        2698 / Tk * IonS ** 1.5 + 1776 / Tk * IonS ** 2 + \
        np.log(1 - 0.001005 * sal))
    # Millero 1983 Chemical Oceanography vol 8
    # partial molar volume and compressibility of HSO4 in seawater.
    deltaVHSO4 = -18.03 + 0.0466 * tem + 0.000316 * tem ** 2
    KappaHSO4 = (-4.53 + 0.09 * tem) / 1000

    #######  Press changed from dbar to bar here by / 10
    lnKhso4fac = (-deltaVHSO4 + 0.5 * KappaHSO4 * (pres_bar)) * \
             (pres_bar) / (R * 10 * Tk)

    #  bisulfate association constant at T, S, P
    Khso4TPS = Khso4 * np.exp(lnKhso4fac)

    # GAMMA +/- HCl, activity coefficient of HCl at T/S, P=1
    # ADH is the Debye Huckel constant, calcualted as a polynomial
    # fit to data in Khoo et al. 1977, doi:10.1021/ac50009a016
    # See Martz et al. 2010, DOI 10.4319/lom.2010.8.172, p175
    # Typo in paper 2nd term should be e-4 not e-6
    #
    ADH = polyval(tem,[0.49172143, 6.7524e-4, 3.4286e-6])
    log10gammaHCl = -ADH * np.sqrt(IonS) / (1 + 1.394 * np.sqrt(IonS)) + \
                (0.08885 - 0.000111 * tem) * IonS
    # Millero 1983 partial molar volume of HCl in seawater
    deltaVHCl = 17.85 + 0.1044 * tem - 0.001316 * tem ** 2

    # effect of pressure on activity coefficient of HCl, divide by 2 because
    # its a mean activity coefficient, divide by 10 for units in the cm3 to F
    # conversion.

    log10gammaHCLtP = log10gammaHCl + deltaVHCl*(pres_bar)/(R*Tk*ln10)/2/10
    #  Sensor reference potential

    # ************************************************************************
    if 'K2f3' in cal.keys():
        K2coefs = np.array([cal['K2f0'],cal['K2f1'],cal['K2f2'],cal['K2f3']])
    elif 'K2f2' in cal.keys():
        K2coefs = np.array([cal['K2f0'],cal['K2f1'],cal['K2f2']])
    elif 'K2f1' in cal.keys():
        K2coefs = np.array([cal['K2f0'],cal['K2f1']])
    else:
        K2coefs = np.array([cal['K2']])
    
    K2 = polyval(pres,K2coefs)
    
    k0T = cal['K0'] + K2 * tem # tem  in deg C

    # CALCULATE PRESSURE CORRECTION (POLYNOMIAL FUNCTION OF PRESSURE)
    # ALL SENSORS HAVE A PRESSURE RESPONSE WHICH IS DETERMINED IN THE LAB
    # AND CONTAINED IN THE POLYNOMIAL Pcoefs
    #pc    = [0].append(cal['Pcoefs']) #descending powers & n+1 (add 0)
    pcorr = polyval(pres,Pcoefs)
    k0TP  = k0T + pcorr

    #  pH on free scale then corrected to get to pH total on mol/kg-seawater scale
    #    pHinsituFree = (Vrs - k0TP) / (R * Tk / F * ln10) + \
    #                   log(Cltotal) / ln10 + 2 * log10gammaHCLtP
    #  this will be mol kg H2O  need to convert to mol/kg-seawater
    phfree = (Vrs - k0TP) / (R * Tk / F * ln10) + \
    np.log(Cltotal) / ln10 + 2 * log10gammaHCLtP #mol/kg-H2O scale
    # CONVERT TO mol/kg-seawater scale - JP 2/4/16
    phfree = phfree - np.log10(1 - 0.001005 * sal) #mol/kg-seawater scale

    # convert to total proton scale
    phtot = phfree - np.log10(1 + Stotal / Khso4TPS)
    ###############################

    return (phfree,phtot)

def ph_mfet(vrs, tempc, sal, k0, k2 = -0.001048):
    """
    Calculate pH from raw MFET sensor output
    
    Parameters
    vrs"""

    tempk = tempc + 273.15
    s_t = (R*tempk)/F * np.log(10)
    z = 19.924*sal/(1000-1.005*sal)
    so4_tot = (0.14/96.062)*(sal/1.80655)
    ccl = 0.99889/35.453*sal/1.80655
    mcl = ccl*1000/(1000-sal*35.165/35)

    # BISULFIDE DISSOCIATION CONSTANT AT T,S AND IONIC STRENGTH(mol/kg solution)
    # Dickson et al. 2007: Chap 5, p12 Eq 33
    Khso4 = np.exp(-4276.1 / tempk + 141.328 - 23.093 * np.log(tempk) + \
        (-13856 / tempk + 324.57 - 47.986 * np.log(tempk)) * z ** 0.5 + \
        (35474 / tempk - 771.54 + 114.723 * np.log(tempk)) * z - \
        2698 / tempk * z ** 1.5 + 1776 / tempk * z ** 2 + \
        np.log(1 - 0.001005 * sal))

    dh_const = 0.00000343 * tempc ** 2 + 0.00067524 * tempc + 0.49172143
    log_gamma_hcl = 2 * (-dh_const * np.sqrt(z) / (1 + 1.394 * np.sqrt(z)) + \
        (0.08885 - 0.000111 * tempc) * z)
    phext_free = -(((k0 + k2 * tempc) - vrs) - s_t * (np.log10(mcl) + log_gamma_hcl )) / s_t
    phext_free = phext_free - np.log10((1000 - sal * 35.165 / 35) / 1000)
    phext_tot = phext_free - np.log10(1 + so4_tot / Khso4)

    return (phext_free,phext_tot)

