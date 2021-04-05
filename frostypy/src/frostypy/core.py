# -*- coding: utf-8 -*-
"""
.. contents:: :local:

Dimensionless Numbers
---------------------
.. autofunction:: Hatta
.. autofunction:: Jakob
.. autofunction:: Fourier
.. autofunction:: Sherwood
.. autofunction:: Nusselt
.. autofunction:: Reynolds

Dimensions
----------
.. autofunction:: thermal_conductivity

"""

import numpy

__all__ = ['Jakob', 'Hatta', 'Fourier', 'Sherwood', 'Nusselt', 'Reynolds']

def Fourier(alpha_a, t, L):
    r'''Calculates the Fourier number 

    ..math::
        Fo = \alpha_a*t\frac(L^2)

    Parameters
    ----------
    alpha_a : float
        Thermal Diffusivity [m^2/s]
    t   : float
        Time [s]
    L : float
        plate length [m]
    
    Returns
    -------
    Fo : Fourier number [-]

    Examples
    --------
    >>>
    
    References
    ----------
    ..[1] (insert here)

    '''   
    return (alpha_a * t) / (L*L)


def Jakob(cp, i_sv, Tsat_a, T_w, omega_a, omega_sat_w):
    r'''Calculates the Jakob number 

    ..math::
        \lambda = \frac(Cp/isv) * (Tsat_a - T_w)/(omega_a - omega_sat_w)

    Parameters
    ----------
    i_sv : float
        Latent heat of sublimation [J/kg]
    Cp   : float
        Specific heat of dry air [J/kg*K] 
    Tsat_a : float
        Saturation temperature of air [C]
    T_w  : float
        Plate surface temperature [C]
    omega_a  : float
        Humidity Ratio [kg_vapor/kg_total]
    omega_sat_w : float
        Humidity ratio @ Tsat
    
    Returns
    -------
    rho : Density of frost (kg/m^3)

    Examples
    --------
    >>>
    
    References
    ----------
    ..[1] (insert here)

    '''        
    return (cp/i_sv) * (Tsat_a - T_w)/(omega_a - omega_sat_w)


def Hatta(x_s, lamb, t, epsilon, D):
    r'''Calculates the Hatta number
    ..math::
        \rho_f = 207*exp(0.266*T_f - 0.0615*T_w)

    Parameters
    ---------   
    x_s     : float
           Frost thickness (m)
    lamb    : float
        Desublimation Coefficient (1/s)
    t       : float
        Time (s)
    epsilon : float
        Porosity (-)
    D       : float
        Diffusivity of water vapor in air (m^2/s)
    
    Returns
    -------
    Ha : Hatta number (-)


    Examples
    --------
    >>>
    
    References
    ----------
    ..[1] (insert here)

    '''  
    return x_s* numpy.sqrt((lamb * t) / (epsilon * D))
    

def Nusselt(h, L, k_a):
    r'''Calculates the Nusselt number
    ..math::
        Nu = h * L / k_a

    Parameters
    ---------   
    h   : float
        Heat transfer coefficient [W/m^2*s]
    L   : float
        Plate length [m]
    k_a : float
        Thermal conductivity of Air [W/mK]
    
    Returns
    -------
    Nu : Nusselt number (-)

    Examples
    --------
    >>>
    
    References
    ----------
    ..[1] (insert here)

    '''      
    return h*L / k_a
    

def Reynolds(V, L, nu):
    r'''Calculates the Reynolds number
    ..math::
        Re = V*L/nu

    Parameters
    ---------   
    V   : float
        Air Velocity [m/s]
    L   : float
        Plate length [m]
    nu : float
        Dynamic Viscosity [m^2/s]
    
    Returns
    -------
    Re : Reynolds number (-)

    Examples
    --------
    >>>
    
    References
    ----------
    ..[1] (insert here)

    ''' 
    return (V*L) / nu
    

def Sherwood(h_m, L, rho_a, D_a):
    r'''Calculates the Sherwood number
    ..math::
        Sh = h_m * L / rho_a * D_a

    Parameters
    ---------   
    h_m  : float
        Mass transfer coefficient [kg/m^2*s]
    L    : float
        Plate length [m]
    rho_a  : float
        Density of Air [kg/m^3]
    D_a : float
        Mass Diffusionn Coefficient [m^2/s]
    
    Returns
    -------
    Sh : Sherwood number (-)


    Examples
    --------
    >>>
    
    References
    ----------
    ..[1] (insert here)

    '''  
    return (h_m * L) / (rho_a * D_a)

def dimensionless_time(k, t, rho, Cp, l):
    r'''Calculates the dimensionless time constant
    ..math::
        \tau = (k*t) / (rho*Cp*l**2)

    Parameters
    ---------   
    k  : float
        Thermal conductivity [W/m^2*K]
    t   : float
        Time [s]
    rho : float
        Density [kg/m^3]
    Cp : float
        Specific heat [J/kg*K]
    l  : float
        Plate length [m]
    
    Returns
    -------
    tau : dimensionless time constant (-)


    Examples
    --------
    >>>
    
    References
    ----------
    ..[1] (insert here)

    '''      
    return (k*t) / (rho*Cp*l**2)


# Dimensional
def thermal_diffusivity(k, rho, Cp):
    r'''Calculates thermal diffusivity or `alpha` for a fluid with the given
    parameters.

    .. math::
        \alpha = \frac{k}{\rho Cp}

    Parameters
    ----------
    k : float
        Thermal conductivity, [W/m/K]
    rho : float
        Density, [kg/m^3]
    Cp : float
        Heat capacity, [J/kg/K]

    Returns
    -------
    alpha : float
        Thermal diffusivity, [m^2/s]

    Notes
    -----

    Examples
    --------
    >>> thermal_diffusivity(k=0.02, rho=1., Cp=1000.)
    2e-05

    References
    ----------
    .. [1] Blevins, Robert D. Applied Fluid Dynamics Handbook. New York, N.Y.:
       Van Nostrand Reinhold Co., 1984.
    '''
    return k/(rho*Cp)