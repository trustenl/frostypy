# -*- coding: utf-8 -*-
"""
.. contents:: :local:

Thermal Conductivity
--------------------
.. autofunction:: Yonko

.. autofunction::Ostin

"""

__all__ = ['Yonko', 'Ostin']

def Yonko(rho_f):
    r'''Calculates the frost thermal conductivity according to the Yonko & Sepsy correlation
    for a flat plate, which is valid for: 
        Ta = 21 [C],
        1.3 < V < 5.3 [m/s],
        7.5 < w < 15 [g/kg_a],
        -30 < Tp < -5 [C],
        0 < t < 200 [min]

    ..math::
        k_f = 0.014 + 0.00668*rho_f + 0.0001759*rho_f**2

    Parameters
    ----------
    rho_f : float 
        Frost density [kg/m^3]
    
    Returns
    -------
    k_f : Frost Thermal Conductivity [J/kg*K]
    
    Examples
    --------
    >>>
    
    References
    ----------
    ..[1] (insert here)

    '''     
    return 0.014 + 0.00668*rho_f + 0.0001759*rho_f**2


def Ostin(rho_f):
    r'''Calculates the frost thermal conductivity according to the Ostin & Andersson correlation
    for parallel plates, which is valid for: 
        20 <= Ta <= 21 [C],
        V = 3 [m/s],
        4.6 <= w <= 10.50 [g/kg_a],
        -20 <= Tp <= -7 [C],
        0 <= t <= 300 [min]

    ..math::
        k_f = -8.71*10**-3 +4.39*10**-4 * rho_f + 1.05*10**-6*rho_f**2

    Parameters
    ----------
    rho_f : float 
        Frost density [kg/m^3]
    
    Returns
    -------
    k_f : Frost Thermal Conductivity [J/kg*K]
    
    Examples
    --------
    >>>
    
    References
    ----------
    ..[1] (insert here)

    '''     
    return -8.71*10**-3 +4.39*10**-4 * rho_f + 1.05*10**-6*rho_f**2