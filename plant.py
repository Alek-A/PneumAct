import numpy as np

class PneumCylinder(object):
    """
    Class for the pneumatic cylinder.
    The cylinder is defined by a second order equation for the piston rod
    and two first order equations for the chambers pressures.
    
    -- References --
    'A High Performance Pneumatic Force Actuator System Part I',
        Richer & Hurmuzlu, 2000
    """
    def __init__(self):
        # Set all parameters to False (unset)
        # Cylinder parameters
        self.ml = False # Load mass [kg]
        self.mp = False # Piston mass [kg]
        self.A1 = False # Inside piston area [m2]
        self.A2 = False # Outside piston area [m2]
        self.Ar = False # Rod area [m2]
        self.Pa = False # Surround pressure, most likely around 10^5 [Pa]
        self.beta = False # Viscous damping coeff. [N.s/m]
        self.Ff = False
        self.Fl = False
        self.S = False # Stroke length
        self.V01 = False # Dead volume, inside [m3]
        self.V02 = False # Dead volume, outside [m3]
        
        # Air chambers
        self.R = 8.31446261815324/28.9647e3 # Specific gas constant of air [J/kg.K]
        self.T = False # Temperature  [K]
        self.k = 1.4 # Specific heat ratio
        self.alpha_in = k 
        self.alpha_out = 1
        self.alpha = 1.2
    
    # Functions for air volumes
    def V1(self,x):
        return self.V01 +self.A1*(self.S/2 +x)
    def V2(self,x):
        return self.V02 +self.A2*(self.S/2 +x)
        
    def updateParameters(self, **kwargs)
        self.__dict__.update(kwargs)
    
    def checkParameters(self):
        checklist = ['ml', 
                     'mp',
                     'A1',
                     'A2',
                     'Ar',
                     'Pa',
                     'beta',
                     'Ff',
                     'Fl',
                     'S',
                     'V01',
                     'V02']
    
    def derivative(self,y,t):
        """
        y = [x, v, p1, p2]
        
        The inputs are through
        mdot11, mdot12, mdot21, mdot22
        """
        yd1 = y[1]
        yd2 = (y[2]*self.A1 -y[3]*self.A2 -self.Pa*self.Ar -self.Fl -self.Ff -self.beta*y[1])/(self.ml +self.mp)
        yd3 = (self.R*self.T*(self.alpha_in*self.mdot11 -self.alpha_out*self.mdot12) -self.alpha*y[2]*self.A1*y[1])/self.V1(y[0])
        yd4 = (self.R*self.T*(self.alpha_in*self.mdot21 -self.alpha_out*self.mdot22) -self.alpha*y[2]*self.A1*y[1])/self.V2(y[0])
        
        return np.array([yd1,yd2,yd3,yd4])
    
class Tube(object):
    def __init__(self):
        self.L = False # Tube length
        self.D = False # Tube diameter
        self.R = 8.31446261815324/28.9647e3 # Specific gas constant of air [J/kg.K]
        self.T = False # Temperature  [K]
        self.k = 1.4 # Specific heat ratio
        self.c = 
class Valve(object):
    def __init__(self):
        self.C1 = 0.040418
        self.C2 = 0.156174
        self.Pcr
        self.nh = False
        self.Rh = False
        self.pw = False
        
    def Avin(self, xs):
        # xs = float, not an array
        if xs <= (self.pw -self.Rh):
            # Orifice holes covered
            return 0
        elif xs > (self.pw -self.Rh) and xs < (self.pw +self.Rh):
            # Orifice holes partially covered
        elif xs > (self.pw +self.Rh):
            # Holes fully uncovered 
            return np.pi*self.nh*self.Rh**2
        else:
            raise ValueError('Something happened in the orifice area equation')
    def Avex(self, xs):
        # xs = float, not an array
        if xs <= -(self.pw +self.Rh):
            # Orifice holes covered
            return 0
        elif xs > -(self.pw +self.Rh) and xs < (self.Rh -self.pw):
            # Orifice holes partially covered
        elif xs > (self.Rh -self.pw):
            # Holes fully uncovered 
            return np.pi*self.nh*self.Rh**2
        else:
            raise ValueError('Something happened in the orifice area equation')
        