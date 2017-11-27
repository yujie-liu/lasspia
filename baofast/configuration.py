import numpy as np

class configuration(object):

    def outputLocation(self):
        return "."

    def inputFilesRandom(self):
        '''List of random catalog file names.'''
        pass
    
    def inputFilesObserved(self):
        '''List of observed catalog file names.'''
        pass

    def binningZ(self): pass
    def binningRA(self): pass
    def binningDec(self): pass
    def binningTheta(self): pass

    '''Parameters for avoiding unnecessary combinatorial calculations at large s.
    Galaxies farther apart than these parameters may not be included in result.'''
    def maxDeltaRA(self): return None
    def maxDeltaDec(self): return None
    def maxDeltaZ(self): return None

    def edgesZ(self): return self.edgesFromBinning(self.binningZ())
    def edgesRA(self): return self.edgesFromBinning(self.binningRA())
    def edgesDec(self): return self.edgesFromBinning(self.binningDec())
    def edgesTheta(self): return self.edgesFromBinning(self.binningTheta())



    def __init__(self):
        pass
    
    @property
    def name(self) : return self.__class__.__name__
        
    def stageFileName(self, stage):
        return '/'.join([self.outputLocation().rstrip('/'),
                         "%s_%s.fits" % (self.name, str(stage))])

    @staticmethod
    def binningDD(binnings):
        binning = {"bins": tuple([b['bins'] for b in binnings])}
        if "range" in binnings[0]:
            binning["range"] = [b["range"] for b in binnings]
        return binning

    @staticmethod
    def edgesFromBinning(binning):
        _,edges = np.histogram([], **binning)
        return edges
