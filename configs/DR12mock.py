# On 5/24/18 ZB changed the theta binning to make the argument type correct
import lasspia as La
import math

class DR12mock(La.configuration):

    def dataDir(self):
        """Directory of catalog files."""
        return '../data/'

    def outputLocation(self):
        return self.dataDir()

    def inputFilesRandom(self):
        return [self.dataDir() + "DR12_mock_modified.fits"]

    def inputFilesObserved(self):
        return [self.dataDir() + "DR12_mock_modified.fits"]

    def catalogRandom(self):
        return La.wrapRandomSDSS(self.inputFilesRandom(), shiftRA=True)

    def catalogObserved(self):
        return La.wrapObservedSDSS(self.inputFilesObserved(), shiftRA=True)

    def binningZ(self): return {"bins":1311, "range":(0.43,1)}
    def binningRA(self): return {"bins":1864 , "range":(125, 245)}
    def binningDec(self): return {"bins":1315, "range":(-10,25)}
    def binningTheta(self): return {"bins":4160, "range":(0,91.7)}

    def chunkSize(self): return 2000

    def maxDeltaRA(self): return None
    def maxDeltaDec(self): return None
    def maxDeltaZ(self): return None

    def binningS(self): return {"bins": 100, "range": (0, 200)}

    '''Configuration affecting only the "integration" routine.'''

    def omegasMKL(self):
        '''Cosmology parameters (\Omega_M, \Omega_K, \Omega_\Lambda).'''
        return (0.274, 0, 0.726)

    def H0(self):
        '''Hubble constant in (h km/s / Mpc)'''
        return 100.

    def lightspeed(self):
        '''Speed of light in km/s'''
        return 299792

