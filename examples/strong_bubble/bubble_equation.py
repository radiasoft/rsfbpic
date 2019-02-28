from scipy.integrate import odeint
import numpy as np


class bubble_integrator:
    
    def __init__(self, delta, sigma_xi, N_b):
        """
        ::float:: delta -- the width of the sheath
        ::float:: sigma_xi -- the rms length of the drive bunch
        ::int:: N_b -- the number of electrons in the drive bunch
        """

        self.derivative = bubble_derivative(delta, sigma_xi, N_b)
    
    
    def compute_bubble(self, xi, r0):
        """
        ::array:: xi -- the list of values for xi we want to compute for
        ::float:: r0 -- initial radius
        """
        
        bubble_values = odeint(self.derivative.derivative, [r0, 0.], xi)
        
        return bubble_values
        

class bubble_derivative:
    
    def __init__(self, delta, sigma_xi, N_b):
        """
        ::float:: delta -- the width of the sheath
        ::float:: sigma_xi -- the rms length of the drive bunch
        ::int:: N_b -- the number of electrons in the drive bunch
        """
        
        self.delta = delta
        self.sigma_xi = sigma_xi
        self.N_b = N_b
    
    
    def A(self, r_b):
        
        beta = self.beta(r_b)
        beta_prime = self.beta_prime(r_b)
        
        A = 0.25 + 0.5*beta + 0.125*r_b*beta_prime
        A *= r_b*r_b
        A += 1.
        
        return A
    
    
    def B(self, r_b):
        
        beta = self.beta(r_b)
        beta_prime = self.beta_prime(r_b)
        beta_dbl_prime = self.beta_dbl_prime(r_b)
        
        B = 0.5 + .75*(beta + r_b*beta_prime) + 0.125*r_b*r_b*beta_dbl_prime
        
        return B
    
    
    def C(self, r_b):
        
        beta = self.beta(r_b)
        
        C = 1./(1. + 0.25*beta*r_b*r_b)**2
        C += 1.
        C *= 0.25
        
        return C
    
    
    def lambda_xi(self, xi):
        # the linear charge density
        
        lambda_xi = self.N_b/np.sqrt(2.*np.pi*self.sigma_xi**2)
        lambda_xi *= np.exp(-xi**2/(2.*self.sigma_xi**2))
        
        return lambda_xi
        
        
    def beta(self, r_b):
        
        beta = self.f(r_b)*self.g(r_b) - 1.
        return beta
        
        
    def beta_prime(self, r_b):
        
        # beta' = (d beta/d alpha) X (d alpha/d r)
        beta_prime = (self.df(r_b)*self.g(r_b) + self.f(r_b)*self.dg(r_b))*self.d_alpha(r_b)
        return beta_prime
    
    
    def beta_dbl_prime(self, r_b):
        
        # beta'' = (d / dr) ((d beta/d alpha) X (d alpha/d r))
        #        = [(d2 beta/ d alpha^2) X (d alpha/ d r)^2] + [(d beta/ d alpha) X (d^2 alpha/ d r^2)]
        beta_dbl_prime = self.d2f(r_b) * self.g(r_b) + 2.*self.df(r_b)*self.dg(r_b) + self.f(r_b)*self.d2g(r_b)
        beta_prime = (self.df(r_b)*self.g(r_b) + self.f(r_b)*self.dg(r_b))
        d_alpha = self.d_alpha(r_b)
        d2_alpha = self.d2_alpha(r_b)
        
        d2beta_dr2 = beta_dbl_prime*d_alpha*d_alpha + beta_prime*d2_alpha
        
        return d2beta_dr2
    
    
    def a(self, r_b):
        
        return 1. + self.delta/r_b
    
    
    def d_alpha(self, r_b):
        
        return -self.delta / (r_b**2)
    
    def d2_alpha(self, r_b):
        
        return 2.*self.delta / (r_b**3)
        
    
    def f(self, r_b):
        
        f = (self.a(r_b)*np.log(self.a(r_b)))**2
        return f
        
    def g(self, r_b):
        
        g = 1./(self.a(r_b)**2 - 1)
        return g
    
    
    def df(self, r_b):

        df = 2.*(np.sqrt(self.f(r_b)) + self.f(r_b)/self.a(r_b))
        return df
    
    def dg(self, r_b):
        
        dg = -2.*self.a(r_b)*(self.g(r_b))**2
        return dg
    
    
    def d2f(self, r_b):
        
        d2f = self.df(r_b)*(1./np.sqrt(self.f(r_b)) + 2./self.a(r_b)) - 2.*self.f(r_b)/self.a(r_b)**2
        return d2f
        
        
    def d2g(self, r_b):
        
        d2g = 2.*self.g(r_b)*self.g(r_b)*(4*self.a(r_b)*self.g(r_b) - 1.)
        return 0
        
    
    def derivative(self, coords, xi):
        
        r = coords[0]
        # this is r'
        u = coords[1]
        
        r_prime = u
        u_prime = (self.lambda_xi(xi)/r - self.C(r)*r - self.B(r)*r*u*u)/self.A(r)
        
        return [r_prime, u_prime]