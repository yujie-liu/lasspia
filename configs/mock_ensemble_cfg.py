import lasspia as La
import math

class mock_ensemble_cfg(La.configuration):

    def dataDir(self):
        """Directory of catalog files."""
        return '../data/'

    def outputLocation(self):
        return self.dataDir()

    def inputFilesRandom(self):
        return [self.dataDir() + "randoms_DR12v5_CMASS_North.fits"]

    def inputFilesObserved(self):
        return [self.dataDir() + "PREPPED_mock_DR12_ens009.fits"]

    def catalogRandom(self):
        return La.wrapRandomSDSS(self.inputFilesRandom(), shiftRA=True)

    def catalogObserved(self):
        return La.wrapObservedSDSS(self.inputFilesObserved(), shiftRA=True)

    def binningZ(self): return {"bins":2357, "range":(0.31,0.99)}
    def binningRA(self): return {"bins": 9340, "range":(96,252)}
    def binningDec(self): return {"bins":4430, "range":(-4,70)}
    def binningTheta(self): return {"bins":1627, "range":(0,1.)}

    def chunkSize(self): return 2000


    '''Parameters for avoiding unnecessary combinatorial calculations at large s.
    Galaxies farther apart than these parameters may not be included in result.'''

    def maxDeltaRA(self): return 19
    def maxDeltaDec(self): return 19
    def maxDeltaZ(self): return 0.13


    '''Configuration affecting only the "integration" routine.'''

    def omegasMKL(self):
        '''Cosmology parameters (\Omega_M, \Omega_K, \Omega_\Lambda).'''
        return (0.31, 0, 0.69)

    def H0(self):
        '''Hubble constant in (h km/s / Mpc)'''
        return 70.

    def lightspeed(self):
        '''Speed of light in km/s'''
        return 299792

    def binningS(self): return {"bins":100, "range":(0,300)}
