"""
Orographic Lift Toy Model Equations

Primary Programmers: 
    Cameron Cummins

11/10/23
"""
import numpy as np
import matplotlib.pyplot as plt


class SaturationVaporPressure():
    """
    Collection of equations approximating saturation vapor pressure for a given temperature T given in Kelvins.
    All equations return saturation vapor pressure in hPa.
    
    Source: https://en.wikipedia.org/wiki/Vapour_pressure_of_water
    """
    def __init__(self, T):
        self.T = T
        self.kPa_to_hPa = 10
        self.mmHg_to_hPa = 75.01
        
    
    def lookup_table(self):
        with open("VaporPressureLookup.json") as f:
            table = json.load(f)

        temps = np.array(table["Temperature (K)"])
        vp = np.array(table["Pressure (hPa)"])

        dt = temps - self.T
        return vp[dt == np.min(dt)][0]
    
    
    def august_eq(self) -> float:
        return np.exp(20.386 - (5132/self.T)) * self.mmHg_to_hPa
    
    
    def august_roche_magnus_eq(self) -> float:
        cT = self.T + 273.15
        return 0.61094*np.exp(17.625*cT/(cT + 243.04))*self.kPa_to_hPa
    
    
    def tetens_eq(self) -> float:
        cT = self.T + 273.15
        return 0.61078*np.exp(17.27*cT/(Ct + 237.3))*self.kPa_to_hPa
    
    
    def buck_eq(self) -> float:
        cT = self.T + 273.15
        return 0.61121*np.exp((18.678 - (cT)/234.5)*(cT / (257.14 + cT)))*self.kPa_to_hPa
    
    
    def goff_gratch_eq(self) -> float:
        steamT = 373.15
        steamPtPressure = 1013.25
        return np.exp(-7.90298*(steamT/(self.T - 1)) + 5.02808*np.log(steamT/self.T) - 1.3816 * (10**(-7)) * (10 ** (11.344 * (1 - self.T / steamT)) - 1) + 8.1328 * 10**(-3) * (10 ** (-3.49149 * (steamT / (self.T - 1))) - 1) + np.log(steamPtPressure))*self.kPa_to_hPa
    
    
def LatentHeatVaporization():
    """
    Collection of equations approximating latent heat of vaporization for a given temperature T given in Kelvins.
    All equations return latent heat of vaporization in kJ/kg.
    
    Lookup Table Source: https://www.engineeringtoolbox.com/water-properties-d_1573.html
    """
    def __init__(self, T):
        self.T = T
            
            
    def lookup_table(self):
        with open("LatentHeatLookup.json") as f:
            table = json.load(f)

        temps = np.array(table["Temperature (K)"])
        latent_heat = np.array(table["Heat of Vaporization (kJ/kg)"])

        dt = temps - self.T
        return latent_heat[dt == np.min(dt)][0]
    
    
    def linear_approx(self):
        return 2.5*(10**6) + self.T * (2.25*(10**6) - 2.5*(10**6)) / 100


class Model():
    """
    Orographic Lift Model
    """
    def __init__(self, z0: float, p0: float, T0: float, q0: float, latent_heat_func, sat_vapor_pressure_func):
        self.z = z0
        self.p = p0
        self.T = T0
        self.q = q0
        self.LH_func = latent_heat_func
        self.SVP_func = sat_vapor_pressure_func
        
        self.__g = -9.81
        self.__R = 287.05
        self.__epsilon = 0.622
    
    
    def update_pressure(self, dz: float):
        g = self.__g
        R = self.__R
        epsilon = self.__epsilon
        
        Tv = self.T / (1 - (pv/self.p)*(1 - epsilon))
        
        self.p = self.p / np.exp((dz * g) / (R * Tv))
    
    
    def timestep(self, dz: float):
        self.z += dz
        
        update_pressure(dz)
        
        
        
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