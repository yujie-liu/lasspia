# On 5/24/18 ZB changed the theta binning to make the argument type correct
import lasspia as La
import math
from cmassS import cmassS

class DR12mock(cmassS):

    def inputFilesObserved(self):
        return [self.dataDir() + "PREPPED_mock_DR12_ens000.fits"]
    def inputFilesRandom(self):
        return [self.dataDir() + "randoms_mocksS_yl0618.fits"]

    def binningZ(self): return {"bins":500, "range":(0.43,0.7)}
    def binningRA(self): return {"bins":1000, "range":(-50,50)}
    def binningDec(self): return {"bins":300, "range":(-10,20)}
    def binningTheta(self): return {"bins":1000, "range":(0,math.pi)}
    def binningS(self): return {"bins":100, "range":(0,200)}

    def maxDeltaRA(self): return None
    def maxDeltaDec(self): return None
    def maxDeltaZ(self): return None
