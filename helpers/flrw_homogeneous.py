from scipy.integrate import ode
import numpy as np
from matplotlib import pyplot as plt

m = 1.0

def V(phi):
    return m**2 * phi**2 / 2.0

def Vprime(phi):
    return m**2 * phi

def Vprimeprime(phi):
    return m**2

def eps(phi):
    return 0.5 * (Vprime(phi) / V(phi))**2

def eta(phi):
    return Vprimeprime(phi) / V(phi)

def Hubble(phi, phi_dot):
    return ((0.5 * phi_dot**2 + V(phi)) / 3.0)**0.5

def f(t,y):
    phi = y[0]
    phi_dot = y[1]
    a = y[2]
    H = Hubble(phi, phi_dot)
    a_dot = H*a
    phi_dot_dot = -3.0 * H * phi_dot - Vprime(phi)
    return [phi_dot, phi_dot_dot, a_dot]

class inflation:
    def __init__(self, phi0, phi0_dot, dt=1.0):
        a0 = 1.0
        y = [phi0, phi0_dot, a0]
        self.ODE = ode(f).set_integrator("dopri5", rtol=1e-8, atol=1e-10, nsteps=100000)
        self.ODE.set_initial_value(y, 0.0)

    def integrate(self, D_t, N):
        dt = D_t * 1.0 / N
        phi = np.zeros(N+1, float)
        phi_dot = np.zeros(N+1, float)
        a = np.zeros(N+1, float)
        t = np.zeros(N+1, float)
        phi[0] = self.ODE.y[0]
        phi_dot[0] = self.ODE.y[1]
        a[0] = self.ODE.y[2]
        t[0] = self.ODE.t
        i = 0
        while (self.ODE.successful() and i < N):
            self.ODE.integrate(self.ODE.t + dt)
            i += 1
            phi[i] = self.ODE.y[0]
            phi_dot[i] = self.ODE.y[1]
            a[i] = self.ODE.y[2]
            t[i] = self.ODE.t
        if (i != N):
            t = t[0:i]
            a = a[0:i]
            phi = phi[0:i]
            phi_dot = phi_dot[0:i]
        H = np.vectorize(Hubble)(phi, phi_dot)
        return [t, phi, phi_dot, a, H]

# If the program is called from the shell, the following is executed:
if __name__ == '__main__':
    efolds = 60.0
    phi0 = np.sqrt((efolds - 0.5) * 4.0)
    phi0_dot = 0.0
    u = inflation(phi0, phi0_dot)
    t_end = 50
    N_plot = 500
    [t, phi, phi_dot, a, H] = u.integrate(t_end, N_plot)
    v = V(phi)
    rho = 0.5 * phi_dot**2 + v

    H_dot = np.gradient(H)
    slow_roll = np.abs(-H_dot / H**2 - 1)
    end_inflation = slow_roll.argmin()

    print('end of inflation at index {}'.format(epsilon))

    np.savetxt("out_N.csv", np.column_stack((t, np.log(a))), delimiter=", ")
    np.savetxt("out_phi.csv", np.column_stack((t, phi)), delimiter=", ")
    np.savetxt("out_a.csv", np.column_stack((t, a)), delimiter=", ")
    np.savetxt("out_H.csv", np.column_stack((t, H)), delimiter=", ")
    np.savetxt("out_V.csv", np.column_stack((phi, m * m * phi * phi / 2.0 )), delimiter=", ")

    plt.subplot(4, 1, 1)
    plt.plot(t, phi)
    plt.ylabel(r"$\phi$")
    plt.xlabel(r"$t$")
    plt.subplot(4, 1, 2)
    plt.plot(t, a)
    plt.ylabel(r"$a$")
    plt.xlabel(r"$t$")
    plt.subplot(4, 1, 3)
    plt.plot(t, np.log(H))
    plt.ylabel(r"$l_H$")
    plt.xlabel(r"$t$")
    plt.subplot(4, 1, 4)
    plt.plot(t, np.log(a))
    plt.ylabel(r"$N$")
    plt.xlabel(r"$t$")
    plt.show()
