import math
from cmassN import cmassN

class cmassN_coarse(cmassN):

    def binningZ(self): return {"bins":621, "range":(0.43,0.7)}
    def binningRA(self): return {"bins":3810 , "range":(125,245)}
    def binningDec(self): return {"bins":381, "range":(-4,8)}
    def binningTheta(self): return {"bins":1904, "range":(0,1.04)}
    def binningS(self): return {"bins":100, "range":(0, 200)}

    def maxDeltaRA(self): return None
    def maxDeltaDec(self): return None
    def maxDeltaZ(self): return None
