import struct
import numpy as np
from Planewave_t import *

def cal_planewave(od, i, j, k):
    """
    od.data["Planewave"]
    typedef struct {
    	double theta, phi;        // direction
    	double ei[3], hi[3];      // E and H unit vector
    	double ri[3], r0[3], ai;  // incidence vector and factor
    	int    pol;               // polarization : 1=V, 2=H
    } planewave_t;                // plane wave incidence
    
    Calculates the electric and magnetic fields of a plane wave at a given point.

    Args:
        freq (float): Frequency of the plane wave in Hz.
        x (float): x-coordinate of the point.
        y (float): y-coordinate of the point.
        z (float): z-coordinate of the point.
        planewave_data (Planewave_t): An object containing the planewave parameters.

    Returns:
        tuple[np.ndarray, np.ndarray]: A tuple containing the electric field vector (e)
                                       and magnetic field vector (h) as NumPy arrays.
                                       Returns (None, None) if planewave is not active.
    """
    freq = 0.0
    for idx in range(od.NFreq2):
        freq = freq + od.data["Freq2"][idx]
    freq = freq / od.NFreq2
    Xn = od.data["Xc"]
    Yn = od.data["Yc"]
    Zn = od.data["Zc"]
    
    C0 = 3e8
    k0 = 2 * np.pi * freq / C0
    
    pw = Planewave_t(od.data["Planewave"])
    x0 = pw.r0[0]
    y0 = pw.r0[1]
    z0 = pw.r0[2]
    
    x = Xn[i]    
    y = Yn[j]    
    z = Zn[k]    

    rri = ((x - x0) * pw.ri[0] +    # planewave_data.ri[0] +
           (y - y0) * pw.ri[1] +    # planewave_data.ri[1] +
           (z - z0) * pw.ri[2])     # planewave_data.ri[2])

    phs = np.exp(-1j * k0 * rri)

    # Note: The ei and hi vectors are unit vectors for the field components.
    # They should be treated as complex values if they represent phase at the origin.
    # The C code implicitly handles them as doubles, which we replicate here.
    # e = planewave_data.ei.astype(np.complex128) * phs
    # h = planewave_data.hi.astype(np.complex128) * phs
    e = pw.ei * phs
    h = pw.hi * phs

    return e, h
    