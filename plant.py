# -*- coding: utf-8 -*-
"""

@author: Aleks
"""

import numpy as np

class PneumCylinder(object):
    """
    Class for the pneumatic cylinder.
    The cylinder is defined by a second order equation for the piston rod
    and two first order equations for the chambers pressures.
    
    -- References --
    [1] 'A High Performance Pneumatic Force Actuator System Part I',
        Richer & Hurmuzlu, 2000
    """
    def __init__(self):
        # Set all parameters to False (unset)
        # Cylinder
        self.A1 = False # Inside piston area [m2]
        self.A2 = False # Outside piston area [m2]
        self.Ar = False # Rod area [m2]
        self.S = False # Stroke length
        self.V01 = False # Dead volume, inside [m3]
        self.V02 = False # Dead volume, outside [m3]
        self.mp = False # Piston mass [kg]
        
        # Load
        self.ml = False # Load mass [kg]
        self.beta = False # Viscous damping coeff. [N.s/m]
        self.Fdf = False # 
        self.Fl = False
        
        # Air chambers
        self.Pa = False # Surround pressure, most likely around 10^5 [Pa]
        self.R = False # Specific gas constant of air [J/kg.K]
        self.T = False # Temperature  [K]
        self.k = False # Specific heat ratio
        self.alpha_in = False 
        self.alpha_out = False
        self.alpha = False
    
    def loadRicherParam(self):
        # Cylinder geom.
        D_pist = 0.75*25.4e-3
        D_rod = 0.25*25.4e-3
        self.A1 = np.pi/4 *D_pist**2 # Inside piston area [m2]
        self.A2 = np.pi/4 *(D_pist**2 -D_rod**2) # Outside piston area [m2]
        self.Ar = np.pi/4 *D_rod**2 # Rod area [m2]
        self.S = 3*25.4e-3 # Stroke length [m]
        self.V01 = self.A1*0.5e-3 # Dead volume, inside [m3], guess
        self.V02 = self.A2*0.5e-3 # Dead volume, outside [m3], guess
        self.mp = 0.1 # Piston mass [kg], guess
        
        # Load
        self.ml = 0 # Load mass [kg]
        self.beta = 4.47 # Viscous damping coeff. [N.s/m] or [Kg/s]
        self.Fdf = 0.486 # Dynamic Coulomb friction [N]
        self.Fl = 0 # Load force
        
        # Air chambers
        self.Pa = 1e5 # Surround pressure, most likely around 10^5 [Pa]
        self.k = 1.4
        self.R = 8.31446261815324/28.9647e-3 # Specific gas constant of air [J/kg.K]
        self.T = 20 +273.15
        self.alpha_in = self.k 
        self.alpha_out = 1
        self.alpha = 1.2
    
    
    # Functions for air volumes
    def V1(self,x):
        return self.V01 +self.A1*(self.S/2 +x)
    def V2(self,x):
        return self.V02 +self.A2*(self.S/2 -x)
    
    
    def derivative(self, y, mdot11, mdot12, mdot21, mdot22):
        """
        y = [x, v, p1, p2]
        
        The inputs are through
        mdot11 : Flow into chamber 1
        mdot12 : Flow out of chamber 1
        mdot21 : Flow into chamber 2
        mdot22 : Flow out of chamber 2
        All positive numbers.
        """
        yd1 = y[1]
        yd2 = (y[2]*self.A1 -y[3]*self.A2 -self.Pa*self.Ar -self.Fl -self.Fdf*np.sign(y[1]) -self.beta*y[1])/(self.ml +self.mp)
        yd3 = (self.R*self.T*(self.alpha_in*mdot11 -self.alpha_out*mdot12) -self.alpha*y[2]*self.A1*y[1])/self.V1(y[0])
        yd4 = (self.R*self.T*(self.alpha_in*mdot21 -self.alpha_out*mdot22) +self.alpha*y[3]*self.A2*y[1])/self.V2(y[0])

        return np.array([yd1,yd2,yd3,yd4])
    
#class Tube(object):
#    def __init__(self):
#        self.L = False # Tube length
#        self.D = False # Tube diameter
#        self.R = 8.31446261815324/28.9647e3 # Specific gas constant of air [J/kg.K]
#        self.T = False # Temperature  [K]
#        self.k = 1.4 # Specific heat ratio
#        self.c =
        

class Valve(object):
    def __init__(self):
        self.k = 1.4
        self.R = 8.31446261815324/28.9647e-3 # Specific gas constant of air [J/kg.K]
        self.C1 = 0.040418
        self.C2 = 0.156174
        self.Cf = False
        self.T = False
        self.n = False
        self.Rh = False
        self.p = False
        self.Ps = False # Supply pressure [Pa]
        self.Pa = False # Exhaust pressure
        self.Pcr = (2/(self.k +1))**(self.k/(self.k -1))
        
    def loadRicherParam(self):
        self.Cf = 0.25
        self.C1 = 0.040418
        self.C2 = 0.156174
        self.T = 20 +273.15
        self.n = 14
        self.Rh = 0.8e-3 
        self.p = 0.9e-3 
        self.Ps = 5e5 
        self.Pa = 1e5 
    
    def Avin(self, x):
        # xs = float, not an array
        if x <= (self.p -self.Rh):
            # Orifice holes covered
            A = 0
        elif x > (self.p -self.Rh) and x < (self.p +self.Rh):
            # Orifice holes partially covered
            A = self.Rh**2 *np.arccos((self.p -x)/self.Rh) -(self.p -x)*np.sqrt(self.Rh**2 -(self.p -x)**2)
        elif x >= (self.p +self.Rh):
            # Holes fully uncovered 
            A = np.pi*self.Rh**2
        else:
            print(x)
            raise ValueError('Something happened in the orifice area equation')
        return A*self.n
    def Avex(self, x):
        x = -x
        # xs = float, not an array
        if x <= (self.p -self.Rh):
            # Orifice holes covered
            A = 0
        elif x > (self.p -self.Rh) and x < (self.p +self.Rh):
            # Orifice holes partially covered
            A = self.Rh**2 *np.arccos((self.p - x)/self.Rh) -(self.p - x)*np.sqrt(self.Rh**2 -(self.p - x)**2)
        elif x >= (self.p +self.R):
            # Holes fully uncovered 
            A = np.pi*self.Rh**2
        else:
            raise ValueError('Something happened in the orifice area equation')
        return A*self.n
    def mdot(self, Pu, Pd, A):
        """
        Mass flow rate [kg/s] as defined by eq. 39-40 in [1]
        
        Pd :: Downstream pressure
        Pu :: Upstream pressure
        A  :: Valve orifice area
        """
#        self.C1 = np.sqrt(self.k/self.R*(2/(self.k +1))**((self.k +1)/(self.k -1)))
#        self.C2 = np.sqrt(2*self.k/self.R/(self.k -1))
#        self.Pcr = (2/(self.k +1))**(self.k/(self.k -1))
        if Pu > Pd:
            if Pd/Pu <= self.Pcr:
                mdv = self.Cf*A*self.C1*Pu/np.sqrt(self.T)
            else:
                mdv = self.Cf*A*self.C2*Pu/np.sqrt(self.T) *(Pd/Pu)**(1/self.k) *np.sqrt(np.abs(1-(Pd/Pu)**((self.k -1)/self.k)))
        else:
            mdv = 0
#        print("Pu = {0:.2f} bar - Pd = {1:.2f} bar - A = {2:.2f} mm2 - Flow: {3}".format(Pu*1e-5, Pd*1e-5, A*1e6, mdv))
        return mdv





