import struct
import numpy as np
from Planewave_t import *

def cal_planewave(od, i, j, k):
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
    e = pw.ei * phs
    h = pw.hi * phs

    return e, h
    