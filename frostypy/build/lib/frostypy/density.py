# -*- coding: utf-8 -*-
"""
.. contents:: :local:

Density Correlations
--------------------
.. autofunction:: Hayashi
.. autofunction:: modified_Hayashi
.. autofunction:: Hermes_2009
.. autofunction:: Hermes_2014
.. autofunction:: Hosoda
.. autofunction:: Kandula
.. autofunction:: Leoni
.. autofunction:: Sommers
.. autofunction:: Yang

"""
import numpy

__all__ = ['Hayashi', 'modified_Hayashi', 'Hermes_2009', 'Hermes_2014', 'Hosoda',
           'Kandula', 'Leoni', 'Yang']

def Hayashi(Tf):    
    r'''Calculates the density of frost according to the modified Hayashi correlation
    for a flat plate, which is valid for: 
        -15 < Tw < -5 [C],
        16 < Ta < 22 [C], 
        0.5 < psi < 0.8 [-],
        u_a = 0.7 [m/s]

    ..math::
        \rho_f = 207*exp(0.266*T_f - 0.0615*T_w)

    Parameters
    ----------
    Tf : float 
        Frost surface temperature [C]
    
    Returns
    -------
    rho : Density of frost (kg/m**3)

    Examples
    --------
    >>>
    
    References
    ----------
    ..[1] (insert here)

    '''
    return 650*numpy.exp(0.277*Tf)


def modified_Hayashi(Tf, Tw):
    r'''Calculates the density of frost according to the modified Hayashi correlation
    for a flat plate, which is valid for: 
        -15 < Tw < -5 [C],
        16 < Ta < 22 [C], 
        0.5 < psi < 0.8 [-],
        u_a = 0.7 [m/s]

    ..math::
        \rho_f = 207*exp(0.266*T_f - 0.0615*T_w)

    Parameters
    ----------
    Tf : float 
        Frost surface temperature [C]
    Tw : float
        Plate surface temperature [C]
    
    Returns
    -------
    rho : Density of frost (kg/m**3)

    Examples
    --------
    >>>
    
    References
    ----------
    ..[1] (insert here)

    '''
    return numpy.exp(0.266*Tf - 0.0615*Tw)


def Hermes_2009(lamb, t):
    r'''Calculates the density of frost according to the Hermes correlation
    for a flat plate, which is valid for: 
        -15 < Tw < -5 [C],
        16 < Ta < 22 [C], 
        0.5 < psi < 0.8 [-],
        u_a = 0.7 [m/s]

    ..math::
        \rho_f = 207*exp(0.266*T_f - 0.0615*T_w)

    Parameters
    ----------
    lamb : float 
        modified Jakob number [-]
    t : float
        time [s]
    
    Returns
    -------
    rho : Density of frost (kg/m**3)

    Examples
    --------
    >>>
    
    References
    ----------
    ..[1] (insert here)

    '''
    return 2.2 * (lamb ** (3/2)) * numpy.sqrt(t)


def Hermes_2014(Cp, Tsata, Tw, hsub, wa, wsatw, t):
    r'''Calculates the density of frost according to the Hermes correlation
    for a flat plate, which is valid for: 
        -15 < Tw < -5 [C],
        16 < Ta < 22 [C], 
        0.05 < psi < 0.8 [%],
        0 < t < 120 [min]
        u_a = 0.7 [m/s]

    ..math::
        \rho_f = 2.2 * ((Cp*(Tsata - Tw)) / (hsub*(wa - wsatw)))**-1.5* t**0.5

    Parameters
    ----------
    Cp : float 
        Specific Heat [J/kg*K]
    Tsata : float
        Saturation Temperature of Air [K]
    Tw : float
        
    hsub : float
        
    wa : float
        Relative Humidity of air [%]
    wsatw : float
        [%]
    t : float
        time [s]
    
    Returns
    -------
    rho : Frost density [kg/m**3]

    Examples
    --------
    >>>
    
    References
    ----------
    ..[1] (insert here)

    '''    
    return 2.2 * ((Cp*(Tsata - Tw)) / (hsub*(wa - wsatw)))**-1.5* t**0.5


def Hosoda(Tw, Va):
    r'''Calculates the frost density according to the Hosoda and Uzuhashi correlation
    for a flat plate, which is valid for: 
        -15 < Tw < -5 [C],
       
    ..math::
        \rho_f = 340*|Tw|**-0.445 + 85*Va

    Parameters
    ----------
    Va : float 
        Air Velocity [m/s]
    Tw : float
        Plate surface temperature [C]
    
    Returns
    -------
    rho : Frost density [kg/m**3]

    Examples
    --------
    >>>
    
    References
    ----------
    ..[1] (insert here)

    '''    
    return 340*numpy.absolute(Tw)**-0.445 + 85*Va


def Kandula(rho_i, Tf, Tw, Tm, Re, Rec):
    r'''Calculates frost density according to the Kandula correlation (2012)
    for a flat plate, which is valid for: 
        -20 < Tw < -5 [C],
        10 < Ta < 22 [C], 
        0.7 < v < 2.5 [m/s]
        50 < RH < 80 [%]
        
    ..math::
        \rho_f = (rho_i * (0.5*((Tf-Tw)/(Tm-Tw)) 
                     * numpy.exp(-(0.375+1.5*((Tf-Tw)/(Tm-Tw)))) 
                     * (1 - numpy.sqrt(Re/Rec))) )

    Parameters
    ----------
    rho_i : float 
        Ice Density [kg/m^3]
    Tf : float
        Temperature of frost  [K]     
    Tw : float
        Temperature of vapor [K]
    Tm : float
        Temperature of  [K]        
    Re : float
        Reynolds number [-]
    Rec : float
        Critical Reynolds number [-]
        
    Returns
    -------
    rho : Frost density [kg/m^3]

    Examples
    --------
    >>>
    
    References
    ----------
    ..[1] (insert here)

    '''    
    return (rho_i * (0.5*((Tf-Tw)/(Tm-Tw)) 
                     * numpy.exp(-(0.375+1.5*((Tf-Tw)/(Tm-Tw)))) 
                     * (1 - numpy.sqrt(Re/Rec))) )
    

def Leoni(Re, Ja, ww, wa, t):
    r'''Calculates the density of frost according to the Leoni et. al (2020) correlation
    for a flat plate, which is valid for: ???
        -35 < Tw < -15 [C],
        -5 < Ta < 15 [C], 
        1.0 < v < 2.5 [m/s]
        0.00322 < w < 0.00847 [kg_v/kg_a]
        
    ..math::
        \rho_f = 5.47*Re**0.16 * Ja**0.29 * (ww/wa)**0.61 * t**0.34

    Parameters
    ----------
    Re : float
        Reynolds number [-]
    Ja : float
        Jakob number [-]        
    ww : float
        Relative Humidity of [%]
    wa : float
        Relative Humidity of air [%]
    t : float
        Time [s]
    
    Returns
    -------
    rho : Frost density [kg/m^3]

    Examples
    --------
    >>>
    
    References
    ----------
    ..[1] (insert here)

    '''    
    return 5.47*Re**0.16 * Ja**0.29 * (ww/wa)**0.61 * t**0.34
  
    
def Sommers(lamb, sigma, RH, t):
    r'''Calculates the density of frost according to the Sommers et al. correlation
    for a flat plate, which is valid for: 
        -15 < Tw < -5 [C],
        16 < Ta < 22 [C], 
        0.5 < psi < 0.8 [-],
        u_a = 0.7 [m/s]

    ..math::
        \rho_f = (2.16*lamb**(-3/2) * t**0.5 * (45/numpy.log(sigma))**0.95 * (0.8/RH)**0.60)

    Parameters
    ----------
    lamb : float 
        modified Jakob number [-]
    sigma : float
        Contact Angle [degrees]
    RH : float
        relative humidity [%]
    t : float
        time [s]
    
    Returns
    -------
    rho : Frost density (kg/m**3)

    Examples
    --------
    >>>
    
    References
    ----------
    ..[1] (insert here)

    '''    
    return (2.16*lamb**(-3/2) * t**0.5 * (45/numpy.log(sigma))**0.95 * (0.8/RH)**0.60)
   
    
def Yang(rho_i, Re, Fo, wa, Ta, Ttp, Tw):
    r'''Calculates the density of frost according to the Yang and Lee correlation
    for a flat plate, which is valid for: 
        -35 < Tw < -15 [C],
        -5 < Ta < 15 [C], 
        1.0 < v < 2.5 [m/s]
        0.00322 < w < 0.00847 [kg_v/kg_a]
        
    ..math::
        \rho_f = 2.2 * ((Cp*(Tsata - Tw)) / (hsub*(wa - wsatw)))**-1.5* t**0.5

    Parameters
    ----------
    rho_i : float 
        Ice Density [kg/m^3]
    Re : float
        Reynolds number [-]
    Fo : float
        Fourier number [-]        
    wa : float
        Relative Humidity of air [%]
    Ta : float
        Temperature of air  [K]
    Ttp : float
        Triple-point of water
    Tw : float
        Temperature of vapor [K]
    
    Returns
    -------
    rho : Frost density [kg/m**3]

    Examples
    --------
    >>>
    
    References
    ----------
    ..[1] (insert here)

    '''
    return (rho_i * (1.54*10**-4 * Re**0.351 * Fo**0.311 * wa**-0.413 * 
                     (Ta-Ttp)/(Ta-Tw)))