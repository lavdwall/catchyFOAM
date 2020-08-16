import cantera as ct
import numpy as np
import scipy.integrate

class PFR:
    def __init__(self, gas, mdot, diam, isothermal = False, pressuredrop = False, U = 0.0, Tr = 273.15, eps = 0.0):
        self.gas = gas
        self.mdot = mdot
        self.diam = diam
        self.A = np.pi*diam**2/4
        self.isothermal = isothermal
        self.deltap = pressuredrop
        self.eps = eps
        self.P = np.pi*diam
        self.U = U
        self.Tr = Tr
        self.multiplier = 1.0

    def __call__(self, z, y):
        """the ODE function, y' = f(z,y) """
        # State vector is [T, p, Y_1, Y_2, ... Y_K]
        self.gas.set_unnormalized_mass_fractions(y[2:])
        self.gas.TP = y[0], y[1]

        wdot = self.gas.net_production_rates * self.multiplier
        dYdz = wdot * self.gas.molecular_weights / self.mdot * self.A
        
        if (self.isothermal):
            dTdz = 0.0
        else:
            dTdz = (-np.dot(self.gas.partial_molar_enthalpies, wdot)*self.A + self.U*self.P*(self.Tr-y[0])) / (self.mdot*self.gas.cp_mass)
        
        if (self.deltap):
            rho = self.gas.density
            u = self.mdot/rho/self.A
            mu = self.gas.viscosity
            Re = rho*self.diam*u/mu
            psi1 = (-2.457*np.log((7/Re)**0.9+0.27*self.eps/self.diam))**16
            psi2 = (37530/Re)**16
            fD = 8*((8/Re)**12+1/(psi1+psi2)**1.5)**(1.0/12.0)
            dpdz = -fD/self.diam*rho/2*u**2
        else:
            dpdz = 0.0

        return np.hstack((dTdz, dpdz, dYdz))
    
    def set_multiplier(self, m):
        self.multiplier = m
    
    def state(self, y):
        self.gas.set_unnormalized_mass_fractions(y[2:])
        self.gas.TP = y[0], y[1]

    def solve(self, y0, length, nsteps):
        dz = length/nsteps

        # Set up solver
        solver = scipy.integrate.ode(self)
        solver.set_integrator('vode', method='bdf', with_jacobian=True, nsteps=10000, atol=1e-15, rtol=1e-9, order=5, min_step = 1e-15)
        solver.set_initial_value(y0, 0.0)

        # Integrate the equations
        restime = 0.0
        states = ct.SolutionArray(self.gas, 1, extra={'z': [0.0], 'tau':[restime], 'MW':[self.gas.mean_molecular_weight]})
        surfstates = ct.SolutionArray(self.surf, 1, extra={'z': [0.0]})
        while solver.successful() and solver.t < length:
            solver.integrate(solver.t + dz)
            self.state(solver.y)
            restime=restime+dz/(self.mdot/self.gas.density/self.A)
            states.append(self.state(solver.y), z=solver.t, tau=restime, MW=self.gas.mean_molecular_weight)
            surfstates.append(self.surf.state, z=solver.t)
        
        if states.z[-1] < length:
            print("PFR simulation did not converge up to specified length")

        return states, surfstates
    
class surfacePFR:
    def __init__(self, gas, surf, mdot, diam, epsB, av, epsC = 0, isothermal = False, pressuredrop = False, U = 0.0, Tr = 273.15):
        self.gas = gas
        self.surf = surf
        self.mdot = mdot
        self.diam = diam
        self.A = np.pi*diam**2/4
        self.epsB = epsB
        self.epsC = epsC
        self.av = av
        self.isothermal = isothermal
        self.deltap = pressuredrop
        self.P = np.pi*diam
        self.U = U
        self.Tr = Tr
        self.multiplier = 1.0

    def __call__(self, z, y):
        """the ODE function, y' = f(z,y) """
        # State vector is [T, p, Y_1, Y_2, ..., Y_K]
        #self.gas.set_unnormalized_mass_fractions(y[2:])
        self.gas.TPY = y[0], y[1], y[2:]
        self.surf.TP = y[0], y[1]

        wdot = self.gas.net_production_rates
        try:
            self.surf.advance_coverages_to_steady_state()
        except:
            pass
        sdot = self.surf.get_net_production_rates(self.gas)
        rdot = ((self.epsB + (1.0-self.epsB)*self.epsC) * wdot + (1.0 - self.epsB) * self.av * sdot) * self.multiplier
        dYdz = rdot * self.gas.molecular_weights / self.mdot * self.A
        
        if (self.isothermal):
            dTdz = 0.0
        else:
            dTdz = (-np.dot(self.gas.partial_molar_enthalpies, rdot)*self.A + self.U*self.P*(self.Tr-y[0])) / (self.mdot*self.gas.cp_mass)  
        
        if (self.deltap):
            rho = self.gas.density
            u = self.mdot/rho/self.A
            mu = self.gas.viscosity
            dpdz = - 150*mu*u*(1-self.epsB)**2/self.epsB**3 - 1.75*rho*u**2*(1-self.epsB)/self.epsB**3
        else:
            dpdz = 0.0

        return np.hstack((dTdz, dpdz, dYdz))
    
    def coverages(self):
        return self.surf.coverages
    
    def set_multiplier(self, m):
        self.multiplier = m
    
    def state(self, y):
        self.gas.set_unnormalized_mass_fractions(y[2:])
        self.gas.TP = y[0], y[1]
        return self.gas.state
        return self.gas.state

    def solve(self, y0, length, nsteps):
        dz = length/nsteps

        # Set up solver
        solver = scipy.integrate.ode(self)
        solver.set_integrator('vode', method='bdf', with_jacobian=True, nsteps=10000, atol=1e-15, rtol=1e-9, order=5, min_step = 1e-15)
        solver.set_initial_value(y0, 0.0)

        # Integrate the equations
        restime = 0.0
        states = ct.SolutionArray(self.gas, 1, extra={'z': [0.0], 'tau':[restime], 'MW':[self.gas.mean_molecular_weight], 'epsB':self.epsB})
        surfstates = ct.SolutionArray(self.surf, 1, extra={'z': [0.0]})
        while solver.successful() and solver.t < length:
            solver.integrate(solver.t + dz)
            self.state(solver.y)
            restime=restime+dz/(self.mdot/self.gas.density/self.A)
            states.append(self.state(solver.y), z=solver.t, tau=restime, MW=self.gas.mean_molecular_weight, epsB=self.epsB)
            surfstates.append(self.surf.state, z=solver.t)
        
        if states.z[-1] < length:
            print("PFR simulation did not converge up to specified length")

        return states, surfstates
    
class surfacePFR_epsB:
    def __init__(self, gas, surf, mdot, diam, epsB, av, epsC = 0, isothermal = False, pressuredrop = False, U = 0.0, Tr = 273.15):
        self.gas = gas
        self.surf = surf
        self.mdot = mdot
        self.diam = diam
        self.A = np.pi*diam**2/4
        self.epsB_func = epsB
        self.epsB = self.epsB_func(0.0)
        self.epsC = epsC
        self.av = av
        self.isothermal = isothermal
        self.deltap = pressuredrop
        self.P = np.pi*diam
        self.U = U
        self.Tr = Tr
        self.multiplier = 1.0

    def __call__(self, z, y):
        """the ODE function, y' = f(z,y) """
        # State vector is [T, p, Y_1, Y_2, ..., Y_K]
        #self.gas.set_unnormalized_mass_fractions(y[2:])
        self.gas.TPY = y[0], y[1], y[2:]
        self.surf.TP = y[0], y[1]
        self.epsB = self.epsB_func(z)

        wdot = self.gas.net_production_rates
        try:
            self.surf.advance_coverages_to_steady_state()
        except:
            pass
        sdot = self.surf.get_net_production_rates(self.gas)
        rdot = ((self.epsB + (1.0-self.epsB)*self.epsC) * wdot + (1.0 - self.epsB) * self.av * sdot) * self.multiplier
        dYdz = rdot * self.gas.molecular_weights / self.mdot * self.A
        
        if (self.isothermal):
            dTdz = 0.0
        else:
            dTdz = (-np.dot(self.gas.partial_molar_enthalpies, rdot)*self.A + self.U*self.P*(self.Tr-y[0])) / (self.mdot*self.gas.cp_mass)  
        
        if (self.deltap):
            rho = self.gas.density
            u = self.mdot/rho/self.A
            mu = self.gas.viscosity
            dpdz = - 150*mu*u*(1-self.epsB)**2/self.epsB**3 - 1.75*rho*u**2*(1-self.epsB)/self.epsB**3
        else:
            dpdz = 0.0

        return np.hstack((dTdz, dpdz, dYdz))
    
    def coverages(self):
        return self.surf.coverages
    
    def set_multiplier(self, m):
        self.multiplier = m
    
    def state(self, y):
        self.gas.set_unnormalized_mass_fractions(y[2:])
        self.gas.TP = y[0], y[1]
        return self.gas.state

    def solve(self, y0, length, nsteps):
        dz = length/nsteps

        # Set up solver
        solver = scipy.integrate.ode(self)
        solver.set_integrator('vode', method='bdf', with_jacobian=True, nsteps=10000, atol=1e-15, rtol=1e-9, order=5, min_step = 1e-15)
        solver.set_initial_value(y0, 0.0)

        # Integrate the equations
        restime = 0.0
        states = ct.SolutionArray(self.gas, 1, extra={'z': [0.0], 'tau':[restime], 'MW':[self.gas.mean_molecular_weight], 'epsB':self.epsB_func(solver.t)})
        surfstates = ct.SolutionArray(self.surf, 1, extra={'z': [0.0]})
        while solver.successful() and solver.t < length:
            solver.integrate(solver.t + dz)
            self.state(solver.y)
            restime=restime+dz/(self.mdot/self.gas.density/self.A)
            states.append(self.state(solver.y), z=solver.t, tau=restime, MW=self.gas.mean_molecular_weight, epsB=self.epsB_func(solver.t))
            surfstates.append(self.surf.state, z=solver.t)
        
        if states.z[-1] < length:
            print("PFR simulation did not converge up to specified length")

        return states, surfstates
