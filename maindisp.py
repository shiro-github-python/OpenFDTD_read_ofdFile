from ReadOutData import *
from Planewave_t import *
from cal_planewave import *
from cal_Intensity import *
from NA import *
from display_v2 import *
import matplotlib.pyplot as plt

# read data file
filename = 'ofd.out'
od = ReadOutData(filename)
od.readoutdata()
d = od.data

# extract complex field
de = d['cfield'][0] 
# Ex = np.zeros((od.Nx, od.Ny, od.Nz), dtype=np.complex128)
# Ey = np.zeros((od.Nx, od.Ny, od.Nz), dtype=np.complex128)
# Ez = np.zeros((od.Nx, od.Ny, od.Nz), dtype=np.complex128)
E = np.zeros((od.Nx, od.Ny, od.Nz), dtype=np.float64)
for ix in range(od.Nx):
    for iy in range(od.Ny):
        for iz in range(od.Nz):
            e, h = cal_planewave(od, ix, iy, iz)
            
            Ex = complex(de['Ex_r'][NA(od, ix, iy, iz)], de['Ex_i'][NA(od, ix, iy, iz)]) + e[0]
            Ey = complex(de['Ey_r'][NA(od, ix, iy, iz)], de['Ey_i'][NA(od, ix, iy, iz)]) + e[1]
            Ez = complex(de['Ez_r'][NA(od, ix, iy, iz)], de['Ez_i'][NA(od, ix, iy, iz)]) + e[2]

            E[ix][iy][iz] = cal_Intensity(Ex, Ey, Ez)

# disp data
display_v2(E, display_res=512)
