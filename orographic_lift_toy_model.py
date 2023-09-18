"""
Orographic Lift Toy Model

Primary Programmers: 
    Cameron Cummins
    Shokoufeh Khojeh
Contributors: 
    Trevor Brooks
    Sasanka Talukdar
    Constance Marshall
    Chieh-Chen Lee

9/17/23
"""

import numpy as np
import matplotlib.pyplot as plt


def calc_potential_temperature(p0: float, p1: float, T: float, R=287.05, cp=1000) -> float:
    """
    Comptues the potential temperature of an air parcel moving adiabatically from pressure p0 to p1
    
    Keyword arguments:
    p0 -- initial pressure in hPa
    p1 -- new pressure in hPa
    T -- temperature of parcel  in Kelvins from p0 to p1
    R -- gas constant (default 287.05)
    cp -- heat capacity of air parcel (default 1000)
    
    """
    return T*(p0/p1)**(R/cp)


def calc_saturated_vapor_pressure(T: float) -> float:
    """
    Calculates the saturated vapor pressure of an air parcel in hPa
    
    Keyword arguments:
    T -- temperature of air parcel in Kelvins
    """
    return 6.11 * np.exp(5423*((1/273.15) - (1/T)))


def calc_latent_heat_vaporization(T: float):
    """
    Uses linear interpolation to determine change in latent heat of vaporization with temperature.
    Maximum range of 0 C to 100 C
    
    Keyword arguments:
    T -- temperature of air parcel in Kelvins
    """
    return 2.5*(10**6) - T * (2.5*(10**6) - 2.25*(10**6)) / 100


def calc_vapor_pressure(q: float, p: float) -> float:
    """
    Computes the vapor pressure of an air parcel in hPa
    
    Keyword arguments:
    q -- specific humidity, dimensionless kilogram over kilogram
    p -- pressure of air parcel in hPa
    """
    return q * p / (0.622 + 0.378 * q)


def calc_specific_humidity(e: float, p: float) -> float:
    """
    Computes specific humidity as a function of pressure and vapor pressure
    
    Keyword arguments:
    e -- vapor pressure of air parcel in hPa
    p -- pressure of air parcel in hPa
    """
    return 0.622*e / (p - 0.378*e)


def calc_pressure(tv: float, p: float, dz: float, g=9.81, R=287.05):
    """
    Computes pressure of air parcel in hPa
    
    Keyword arguments:
    tv -- virtual temperature of air parcel in Kelvins
    p -- pressure of air parcel in hPa
    dz -- air parcel change in altitude in meters
    g -- graviational constant (default 9.81)
    R -- gas constant for dry air (default 287.05)
    """
    exponential_term = np.exp((g / (tv*R)) * (dz))
    
    return p / exponential_term


def calc_virtual_temperature(t: float, p: float, pv: float, epsilon=0.622) -> float:
    """
    Computes virtual temperature of air parcel in Kelvins
    
    Keyword arguments:
    t -- temperature of air parcel in Kelvins
    p -- pressure of air parcel in hPa
    pv -- vapor pressure of air parcel in hPa
    epsilon -- ratio of molar mass of water vapor to molar mass of dry air (default 0.622)
    """
    # Epsilon is molar mass of water vapor over molar mass of dry air
    return t / (1 - (pv/p)*(1 - epsilon))


def calc_dry_lapse_rate(g=9.81, cp=1000) -> float:
    """
    Computes the adiabatic lapse rate of a dry air parcel in Kelvins per meter
    
    Keyword arguments:
    g -- graviational constant (default 9.81)
    cp -- heat capacity of air parcel (default 1000)
    """
    return -1*(g/cp)


def calc_wet_lapse_rate_coeff(w: float, l: float, cp=1000) -> float:
    """
    Computes the coefficient to multiply the adiabatic lapse rate for a dry parcel of air by
    
    Keyword arguments:
    w -- mixing ratio, dimensionless
    l -- latent heat of condensation in Joules per kilogram
    cp -- heat capacity of dry air (default 1000)
    """
    return 1 / (1 + (l/cp)*w)


def calc_mixing_ratio(q: float) -> float:
    """
    Estimates the mixing ratio of an air parcel using specific humidity
    
    Keyword arguments:
    q -- specific humidity of air parcel, dimensionless kilogram over kilogram
    """
    return q / (1 - q)


def altitude_function(x: float, peak_x=0, peak_height=3000, base_width=1000) -> float:
    """
    Returns relative altitude at the specified displacement from the mountain peak formed by am isosceles triangle.
    
    Keyword arguments:
    x -- displacement from mountain peak in meters
    peak_x -- x position of mountain peak in meters
    peak_height -- vertical height of peak at peak_x in meters
    base_width -- horizontal distance from z = 0 to mountain peak in meters
    """
    slope = peak_height / base_width
    
    if x > peak_x:
        if x < peak_x + base_width:
            return peak_height - slope*(x - peak_x)
    else:
        if peak_x - 1*base_width < x:
            return slope*(base_width + x - peak_x)
    return 0


if __name__ == '__main__':
    # Initial conditions @ x_initial
    x_coords = np.arange(-1000, 1010, 10)
    z_coords = np.zeros(x_coords.shape)

    temperatures = np.zeros(x_coords.shape)
    pressures = np.zeros(x_coords.shape)
    specific_humidities = np.zeros(x_coords.shape)
    relative_humidities = np.zeros(x_coords.shape)
    potential_temperatures = np.zeros(x_coords.shape)
    vapor_pressures = np.zeros(x_coords.shape)

    temperatures[0] = 20 + 273.15
    pressures[0] = 1000
    specific_humidities[0] = 0.01
    vapor_pressures[0] = calc_vapor_pressure(specific_humidities[0], pressures[0])
    relative_humidities[0] = vapor_pressures[0] / calc_saturated_vapor_pressure(temperatures[0])
    potential_temperatures[0] = calc_potential_temperature(pressures[0], pressures[0], temperatures[0])
    z_coords[0] = altitude_function(x_coords[0])

    # Run toy model
    for index in range(1, x_coords.shape[0]):
        # Get state
        temperature = temperatures[index-1]
        pressure = pressures[index-1]
        specific_humidity = specific_humidities[index-1]
        vapor_pressure = vapor_pressures[index-1]

        # Compute derivative variables
        saturated_vapor_pressure = calc_saturated_vapor_pressure(temperature)
        latent_heat = calc_latent_heat_vaporization(temperature)
        lapse_rate = calc_dry_lapse_rate()
        mixing_ratio = calc_mixing_ratio(specific_humidity)

        # Impose constraint
        specific_humidities[index] = specific_humidities[index-1]
        if vapor_pressure > saturated_vapor_pressure:
            vapor_pressure = saturated_vapor_pressure

            # Modify lapse rate
            lapse_rate = lapse_rate*calc_wet_lapse_rate_coeff(mixing_ratio, latent_heat)

            # Compute new specific humidity
            specific_humidities[index] = calc_specific_humidity(vapor_pressure, pressure)

        # Update vapor pressure
        vapor_pressures[index] = calc_vapor_pressure(specific_humidity, pressure)

        # Update altitude
        z_coords[index] = altitude_function(x_coords[index])
        z_delta = z_coords[index] - z_coords[index-1]

        # Update temperature
        temperatures[index] = temperature + lapse_rate*z_delta

        # Update pressure
        virtual_temperature = calc_virtual_temperature(temperature, pressure, vapor_pressure)
        pressures[index] = calc_pressure(virtual_temperature, pressure, z_delta)

        # Update relative humidity
        relative_humidities[index] = calc_vapor_pressure(specific_humidity, pressure) / calc_saturated_vapor_pressure(temperature)

        # Update potential temperature
        potential_temperatures[index] = calc_potential_temperature(pressure, pressures[index], temperature)

    # =================================== Plotting ===================================

    f, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2, figsize=(20, 18), facecolor='w')
    f.suptitle("Air Parcel Simulation State Variables", fontsize=30)

    lw = 2

    ax1.plot(x_coords, temperatures - 273.15, color="Red", label="Temperature")
    ax1.plot(x_coords, potential_temperatures - 273.15, color="Blue", label="Potential Temperature")
    ax2.plot(x_coords, pressures, color="Red", label="Pressure")
    ax3.plot(x_coords, specific_humidities, color="Red", label="Specific Humidity")
    ax4.plot(x_coords, relative_humidities*100, color="Red", label="Relative Humidity")
    ax5.plot(x_coords, vapor_pressures, color="Red", label="Altitude")
    ax6.plot(x_coords, z_coords, color="Red", label="Altitude")

    title_s = 20
    label_s = 12
    legend_s = 15

    ax1.set_title("Temperature", fontsize=title_s)
    ax2.set_title("Pressure", fontsize=title_s)
    ax3.set_title("Specific Humidity", fontsize=title_s)
    ax4.set_title("Relative Humidity", fontsize=title_s)
    ax5.set_title("Vapor Pressure", fontsize=title_s)
    ax6.set_title("Altitude", fontsize=title_s)

    ax1.set_ylabel("Temperature (°C)", fontsize=label_s)
    ax2.set_ylabel("Pressure (hPa)", fontsize=label_s)
    ax3.set_ylabel("Specific Humidity (g/kg)", fontsize=label_s)
    ax4.set_ylabel("Relative Humidity (%)", fontsize=label_s)
    ax5.set_ylabel("Vapor Pressure (hPa)", fontsize=label_s)
    ax6.set_ylabel("Altitude (m)", fontsize=label_s)

    ax1.set_xlabel("Horizontal Distance to Peak (m)", fontsize=label_s)
    ax2.set_xlabel("Horizontal Distance to Peak (m)", fontsize=label_s)
    ax3.set_xlabel("Horizontal Distance to Peak (m)", fontsize=label_s)
    ax4.set_xlabel("Horizontal Distance to Peak (m)", fontsize=label_s)
    ax5.set_xlabel("Horizontal Distance to Peak (m)", fontsize=label_s)
    ax6.set_xlabel("Horizontal Distance to Peak (m)", fontsize=label_s)

    ax1.legend(fontsize=legend_s)
    ax2.legend(fontsize=legend_s)
    ax3.legend(fontsize=legend_s)
    ax4.legend(fontsize=legend_s)
    ax5.legend(fontsize=legend_s)
    ax6.legend(fontsize=legend_s)

    ax1.grid()
    ax2.grid()
    ax3.grid()
    ax4.grid()
    ax5.grid()
    ax6.grid()

    f.savefig("Simulation_Output.png")