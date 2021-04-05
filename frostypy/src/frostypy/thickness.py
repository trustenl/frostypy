# -*- coding: utf-8 -*-
"""
.. contents:: :local:

Frost Thickness
---------------
.. autofunction:: Cremers_free
.. autofunction:: Cremers
.. autofunction:: Lee
.. autofunction:: Okoroafor
.. autofunction:: Schneider

"""

__all__ = ['Cremers_free', 'Cremers', 'Hermes_2012', 'Okoroafor', 'Schneider']

def Cremers_free(Tf, Tp, t):
    r'''Calculates the frost thickness according to the Cremers & Mehra correlation
    for a tube, which is valid for: 
        Ta = 24 [C]
        Free convection
        9.2 < w < 15 [g/kg_a],
        45 < psi < 90 [%], 

    ..math::
        y_f = 0.20*t*(Tf - Tp)**0.4

    Parameters
    ----------
    Tf : float 
        Frost Temperature [K]
    Tp : float
        Cold Plate Surface [K]
    t : float
        time [s]
    
    Returns
    -------
    y_f : Frost Thickness [m]
    
    Examples
    --------
    >>>
    
    References
    ----------
    ..[1] (insert here)

    '''    
    return 0.20*t*(Tf - Tp)**0.4


def Cremers(Tf, Tp, t):
    r'''Calculates the frost thickness according to the Cremers & Mehra correlation
    for a tube, which is valid for: 
        -90 < Tp < -60 [C]
        0 < t < 600 [min]
        34 < psi < 53 [%], 

    ..math::
        y_f = 0.12*t*(Tf - Tp)**0.43

    Parameters
    ----------
    Tf : float 
        Frost Temperature [K]
    Tp : float
        Cold Plate Surface [K]
    t : float
        time [s]
    
    Returns
    -------
    y_f : Frost Thickness [m]
    
    Examples
    --------
    >>>
    
    References
    ----------
    ..[1] (insert here)

    '''    
    return 0.12*t*(Tf - Tp)**0.43


def Hermes_2012():
    return 

def Okoroafor(b0, b1, t, n):
    return b0 + b1*t**n

def Schneider(Tf, Tm, Tw, k_i, hsub, rho_i, t, p, pfsat, psat):
    r'''Calculates the frost density according to Schneider et. al (2020) correlation
    for a flat plate, which is valid for: ???
        -35 < Tw < -15 [C],
        -5 < Ta < 15 [C], 
        1.0 < v < 2.5 [m/s]
        0.00322 < w < 0.00847 [kg_v/kg_a]
        
    ..math::
        \rho_f = 5.47*Re**0.16 * Ja**0.29 * (ww/wa)**0.61 * t**0.34

    Parameters
    ----------
    Ta : float
        Air Temperature [C]
    Tw : float
        Plate Temperature [C]
    Tm : float
        Melting-point Temperature of water-ice [C]
    Tf : float
        Frost Temperature [C]
    pv : float
        Partial vapor pressure of the air [-]        
    psat : float
        Vapor pressure of saturated air [%]
    psatf : float
        Vapor pressure of saturated air at the frost surface temperature [-]
    kice : float
        Ice thermal conductivity (kJ/kg*K)
    pice : float
        Ice density (kg/m^3)
    isv : float
        
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
    
    Ft = 1 + 0.052 * ((Tf - Tm) / (Tm - Tw))
    return (0.465 * ((k_i * (Tf - Tw))/ (hsub * rho_i))**0.5 * (t/3600)**-0.03 
         * (Tf - Tm)**0.01 * ((p-pfsat) / (psat - pfsat))**0.25 * Ft)


def Sengupta(d, Re, Pr, w, Fo):
    r'''Calculates the frost thickness according to the Sengupta et al. correlation
    for a tube, which is valid for: 
        Ta = 29 [C]
        1.5 < V < 4.4 [m/s]
        10.0 < w < 20.0 [g/kg_a]
        
    ..math::
        y_f = 0.84 * d * Re**-0.15 * Pr**0.65 * (1+w)**0.71 * Fo**0.11

    Parameters
    ----------
    d : float 
        diameter [m]
    Re : float
        Reynolds number [-]
    Pr : float
        Prandtl number [-]
    w : float
        Relative Humidity [%]
    Fo : float
        Fourier number [-]
    
    Returns
    -------
    y_f : Frost Thickness [m]
    
    Examples
    --------
    >>>
    
    References
    ----------
    ..[1] (insert here)

    '''    
    return 0.84 * d * Re**-0.15 * Pr**0.65 * (1+w)**0.71 * Fo**0.11