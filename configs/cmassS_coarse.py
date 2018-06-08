import lasspia as La
import math
from cmassS import cmassS

class cmassS_coarse(cmassS):

    def binningZ(self): return {"bins":1404, "range":(0.43,0.7)}
    def binningRA(self): return {"bins": 6350, "range":(-50,50)} # 6350

    def binningS(self): return {"bins": 100, "range": (0, 200)}
    def maxDeltaRA(self): return None
    def maxDeltaDec(self): return None
    def maxDeltaZ(self): return None
