# -*- coding: utf-8 -*-
"""
.. contents:: :local:

Porosity
---------------
.. autofunction:: Hermes

"""

__all__ = ['Hermes']

# Frost Porosity
def Hermes(rho_f, rho_i, rho_a):
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
    return (rho_f - rho_i)/(rho_a - rho_i)
