import numpy as np

def cal_Intensity(Cx, Cy, Cz):
    intensity = np.real(Cx*np.conj(Cx) + Cy*np.conj(Cy) + Cz*np.conj(Cz))
    return intensity