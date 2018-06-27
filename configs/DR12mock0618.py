import lasspia as La
import math
from cmassS import cmassS

class DR12mock0618(cmassS):

    def binningZ(self): return {"bins":860, "range":(0.43,0.7)}
    def binningRA(self): return {"bins": 6667, "range":(-50,50)}
    def binningDec(self): return {"bins":2000, "range":(-10,20)}
    def binningTheta(self): return {"bins":12000, "range":(0,math.pi)}

    def binningS(self): return {"bins":100, "range":(0,200)}
    def maxDeltaRA(self): return None
    def maxDeltaDec(self): return None
    def maxDeltaZ(self): return None

