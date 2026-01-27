import struct
import numpy as np

class Planewave_t:
    """
    A class to hold planewave parameters, equivalent to the C struct planewave_t.
    """

    def __init__(self, pw):
        self.theta = pw[0]
        self.phi = pw[1]
        self.ei = np.array([pw[2], pw[3], pw[4]])
        self.hi = np.array([pw[5], pw[6], pw[7]])
        self.ri = np.array([pw[8], pw[9], pw[10]])
        self.r0 = np.array([pw[11], pw[12], pw[13]])
        self.ai = pw[14]
        self.pol = pw[15]