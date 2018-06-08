import math
from cmassN import cmassN

class cmassN_coarse(cmassN):

    def binningZ(self): return {"bins":150, "range":(0.43,0.7)}
    def binningRA(self): return {"bins":800 , "range":(105,265)}
    def binningDec(self): return {"bins":305, "range":(-4,57)}
    def binningTheta(self): return {"bins":500, "range":(0,math.pi/2)}

    def maxDeltaRA(self): return None
    def maxDeltaDec(self): return None
    def maxDeltaZ(self): return None
