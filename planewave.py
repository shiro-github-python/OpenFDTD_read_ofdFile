

import numpy as np
import math

# Constants
C0 = 2.99792458e8

class Planewave_t:
    """
    A class to hold planewave parameters, equivalent to the C struct planewave_t.
    """
    def __init__(self, theta=0, phi=0, ei=None, hi=None, ri=None, r0=None, ai=0, pol=0):
        self.theta = theta
        self.phi = phi
        self.ei = np.zeros(3, dtype=np.float64) if ei is None else np.array(ei, dtype=np.float64)
        self.hi = np.zeros(3, dtype=np.float64) if hi is None else np.array(hi, dtype=np.float64)
        self.ri = np.zeros(3, dtype=np.float64) if ri is None else np.array(ri, dtype=np.float64)
        self.r0 = np.zeros(3, dtype=np.float64) if r0 is None else np.array(r0, dtype=np.float64)
        self.ai = ai
        self.pol = pol

def setup_planewave(theta_deg, phi_deg, pol, Xn, Yn, Zn, Freq2):
    """
    Calculates the detailed parameters for the Planewave_t object.
    This is a Python translation of the setup_planewave function from sol/input2.c.

    Args:
        theta_deg (float): The angle theta in degrees.
        phi_deg (float): The angle phi in degrees.
        pol (int): Polarization (1 for V-pol, 2 for H-pol).
        Xn (np.ndarray): Array of node positions along the x-axis.
        Yn (np.ndarray): Array of node positions along the y-axis.
        Zn (np.ndarray): Array of node positions along the z-axis.
        Freq2 (np.ndarray): Array of frequencies.

    Returns:
        Planewave_t: A fully populated Planewave_t object.
    """
    p = Planewave_t(theta=theta_deg, phi=phi_deg, pol=pol)

    cost = math.cos(math.radians(theta_deg))
    sint = math.sin(math.radians(theta_deg))
    cosp = math.cos(math.radians(phi_deg))
    sinp = math.sin(math.radians(phi_deg))

    # Unit vectors in (r, theta, phi) spherical coordinates
    r1 = np.array([+sint * cosp, +sint * sinp, +cost])
    t1 = np.array([+cost * cosp, +cost * sinp, -sint])
    p1 = np.array([-sinp, +cosp, 0])

    # Propagation vector (ri) is opposite to the radial unit vector (r1)
    p.ri = -r1

    # Electric field unit vector (ei) based on polarization
    if pol == 1:  # V-pol
        p.ei = -t1
    elif pol == 2:  # H-pol
        p.ei = +p1
    
    # Magnetic field unit vector (hi) is the cross product of E and r (H = E x r)
    # Note: In the source, it's calculated as (E x r1), which is equivalent to k x E / eta
    p.hi = np.cross(p.ei, r1)

    # Initial position (r0)
    f0 = (Freq2[0] + Freq2[-1]) / 2
    if f0 == 0:
        # Avoid division by zero if frequencies are not set
        f0 = 1e9 # default to 1GHz
        
    Nx, Ny, Nz = len(Xn) - 1, len(Yn) - 1, len(Zn) - 1
    r_diag = math.sqrt((Xn[0] - Xn[Nx])**2 + (Yn[0] - Yn[Ny])**2 + (Zn[0] - Zn[Nz])**2) / 2
    r = r_diag + (0.5 * C0 / f0)
    
    center = np.array([(Xn[0] + Xn[Nx]) / 2, (Yn[0] + Yn[Ny]) / 2, (Zn[0] + Zn[Nz]) / 2])
    p.r0 = center - (r * p.ri)

    # Waveform parameter (ai)
    p.ai = 4 / (1.27 / f0) if f0 > 0 else 0

    return p


def planewave(freq, x, y, z, i_planewave, planewave_data):
    """
    Calculates the electric and magnetic fields of a plane wave at a given point.

    Args:
        freq (float): Frequency of the plane wave in Hz.
        x (float): x-coordinate of the point.
        y (float): y-coordinate of the point.
        z (float): z-coordinate of the point.
        i_planewave (int): Flag to indicate if planewave is active (1 for active).
        planewave_data (Planewave_t): An object containing the planewave parameters.

    Returns:
        tuple[np.ndarray, np.ndarray]: A tuple containing the electric field vector (e)
                                       and magnetic field vector (h) as NumPy arrays.
                                       Returns (None, None) if planewave is not active.
    """
    if not i_planewave:
        return None, None

    k0 = 2 * np.pi * freq / C0
    x0, y0, z0 = planewave_data.r0

    rri = (x - x0) * planewave_data.ri[0] +
          (y - y0) * planewave_data.ri[1] +
          (z - z0) * planewave_data.ri[2]

    phs = np.exp(1j * k0 * rri)

    # Note: The ei and hi vectors are unit vectors for the field components.
    # They should be treated as complex values if they represent phase at the origin.
    # The C code implicitly handles them as doubles, which we replicate here.
    e = planewave_data.ei.astype(np.complex128) * phs
    h = planewave_data.hi.astype(np.complex128) * phs

    return e, h

if __name__ == '__main__':
    # Example usage:
    # This is a demonstration of how to use the setup and planewave functions.

    # 1. Define simulation parameters (grid, frequency)
    # These would normally be loaded from your simulation setup.
    Xn = np.linspace(-0.5, 0.5, 101)  # 1 meter analysis domain, 100 cells
    Yn = np.linspace(-0.5, 0.5, 101)
    Zn = np.linspace(-0.5, 0.5, 101)
    Freq2 = np.linspace(0.5e9, 1.5e9, 10) # 0.5 GHz to 1.5 GHz

    # 2. Setup Planewave data using the new function
    # Inputs: theta=45 deg, phi=90 deg, Vertical Polarization
    planewave_params = setup_planewave(45, 90, 1, Xn, Yn, Zn, Freq2)
    i_planewave_active = 1

    print("--- Calculated Planewave Parameters ---")
    print(f"Theta: {planewave_params.theta}, Phi: {planewave_params.phi}, Pol: {planewave_params.pol}")
    print(f"Incidence vector (ri): {planewave_params.ri}")
    print(f"E-field vector (ei):   {planewave_params.ei}")
    print(f"H-field vector (hi):   {planewave_params.hi}")
    print(f"Reference point (r0):  {planewave_params.r0}")
    print(f"Waveform param (ai):   {planewave_params.ai}")
    print("-" * 37)


    # 3. Define point and frequency for field calculation
    frequency = 1e9  # 1 GHz
    point_x, point_y, point_z = (0.1, 0.2, 0.3)

    # 4. Calculate planewave fields at the specified point
    e_field, h_field = planewave(
        frequency,
        point_x,
        point_y,
        point_z,
        i_planewave_active,
        planewave_params
    )

    # 5. Print results
    if e_field is not None and h_field is not None:
        print(f"\nCalculation at point ({point_x}, {point_y}, {point_z}) and frequency {frequency/1e9} GHz:")
        print(f"  Electric Field (E):")
        print(f"    Ex: {e_field[0]}")
        print(f"    Ey: {e_field[1]}")
        print(f"    Ez: {e_field[2]}")
        print(f"  Magnetic Field (H):")
        print(f"    Hx: {h_field[0]}")
        print(f"    Hy: {h_field[1]}")
        print(f"    Hz: {h_field[2]}")
    else:
        print("\nPlanewave incidence is not active.")
