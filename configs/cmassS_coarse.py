import lasspia as La
import math
from cmassS import cmassS

class cmassS_coarse(cmassS):

    def binningZ(self): return {"bins":1404, "range":(0.43,0.7)}
    def binningRA(self): return {"bins": 2000, "range":(-50,50)} # 6350
    def binningDec(self): return {"bins":1270, "range":(-10,10)}
    def binningTheta(self): return {"bins":1206, "range":(0,19)}

    def binningS(self): return {"bins": 100, "range": (0, 200)}
    def maxDeltaRA(self): return None
    def maxDeltaDec(self): return None
    def maxDeltaZ(self): return None
