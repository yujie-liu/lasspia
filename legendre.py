import numpy as np
import sys
import os.path as osp
from astropy.io import fits
from scipy.special import legendre
import time


# From lasspia.py
def getInstance(argFile, args=(), kwargs={}):
    path = argFile[0].split('/')
    name = path[-1].split('.')[0]
    sys.path.append('/'.join(path[:-1]))
    exec("from %s import %s " % (name, name))
    return eval(name)(*args, **kwargs)


# From lasspia.py
def getKWs(args):
    if args.nCores or args.iJob:
        n = args.nJobs[0] if args.nJobs else 1
        jobs = args.iJob if args.iJob else range(args.nJobs[0])
        return [{"nJobs": n, "iJob": i} for i in jobs]
    if args.nJobs: return {"nJobs": args.nJobs[0]}
    return {}


# From lasspia.py
def getCfgArgs(args):
    cfgArgs = {'txtToFile': args.txtToFile}
    if args.iSliceZ: cfgArgs['iSliceZ'] = args.iSliceZ[0]
    return cfgArgs


def get_tpcf(fits_file, config):
    '''
    Calculate tpcf(sigma, pi) from DD, RR and DR and convert to tpcf(s, mu)
    Inputs:
    + fits: the fits file
    + config: configuration object
    Outputs:
    + s_vec: ndarray
        the bin centers of s as a vector
    + tpcf: ndarray
        the tpcf of s and mu as a 2-d array
    '''
    hdul = fits.open(fits_file)
    centerS = hdul[2].data['bincenter']
    # normalization factors from integration.py
    nRR, nDR, nDD = (lambda h:
                     (h['normrr'],
                      h['normdr'],
                      h['normdd']))(fits.getheader(fits_file, 'TPCF'))
    sigma_vals = hdul[3].data['binCenter']
    print(hdul[2].data['binCenter'])
    pi_vals = hdul[4].data['binCenter']
    sigma_indices = hdul[5].data['iSigma']
    pi_indices = hdul[5].data['iPi']
    s_vals = (np.square(sigma_vals) + np.square(pi_vals)) ** .5
    mu_vals = pi_vals / s_vals
    # test1 = [sigma_vals[i] for i in sigma_indices]
    # test2 = [pi_vals[i] for i in pi_indices]
    # test_s = (np.square(test1) + np.square(test2)) ** .5
    # test_mu = test2/test_s
    # s_index = (test_s/2)
    # mu_index = test_mu * 500/2 + 249
    bin_theta, range_theta = config.binningTheta()['bins'], config.binningTheta()['range']
    bin_theta = int(bin_theta)
    bin_s, range_s = config.binningS()['bins'], config.binningS()['range']
    bin_s = int(bin_s)
    s_vec = np.arange(1, int(range_s[1]), int(range_s[1] / bin_s))
    tpcf = np.empty((bin_s, bin_theta))
    for tuple in hdul[5].data:
        isigma, ipi, rr, dr, dd, dde2 = tuple
        s = (sigma_vals[isigma] ** 2 + pi_vals[ipi] ** 2) ** .5
        mu = pi_vals[ipi] / s
        i_s = min(int(s / (range_s[1] / bin_s)), int(range_s[1] / 2 - 1))
        i_mu = min(int((mu + 1) * bin_theta / 2), bin_theta - 1)
        tpcf[i_s, i_mu] = (dd / nDD - 2 * dr / nDR + rr / nRR) / (rr / nRR)
    return s_vec, tpcf


def legendre_coef(tpcf, l, config):
    """
    Calculate the coefficient of tpcf_l(s) in Legendre expansion for all s
    Inputs:
    + tpcf: ndarray
        the tpcf as a matrix dimensioned by mu and s
    + l: the order of desired Legendre polynomial
    Outputs:
    + coef:
        Return the coefficient of tpcf_l(s) as an array indexed by s values
    + write to the *_legendre.fits file
        The .fits file has columns (s, tpcf0, etpcf0, tpcf2, etpcf2, ...)
    """
    if l % 2 == 1:
        raise ValueError("The order of Legendre polynomials to be even")
    delta_mu = 2 / config.binningTheta()['bins']
    temp = np.copy(tpcf)
    mu_vec = np.linspace(-1, 1, int(config.binningTheta()['bins']))[:, None]  # the vector of mu values
    p = legendre(l)  # Legendre polynomial at order l
    mu_vec = np.array([p(mu) for mu in mu_vec])  # P-_l(mu)
    temp = temp * mu_vec.transpose() * delta_mu  # P_l(mu) * tpcf(s, mu) * delta_mu
    temp = temp * (2 * l + 1) / 2
    return temp.sum(axis=1)


def legendre_to_file():
    from argparse import ArgumentParser
    parser = ArgumentParser(description="Large Scale Structure Probability Integration Algorithm")
    parser.add_argument('configFile', metavar='configFile', type=str, nargs=1,
                        help='A python file containing a subclass of lasspia.configuration of the same name.')
    parser.add_argument('fitsFile', metavar='fitsFile', type=str, nargs=1,
                        help='FITS file.')
    args = parser.parse_args()
    config = getInstance(args.configFile)
    fits_file = osp.join(osp.abspath('../data'), args.fitsFile[0])
    s_vec, tpcf = get_tpcf(fits_file, config)
    hdu = fits.BinTableHDU.from_columns([
        fits.Column(name='s', array=s_vec, format='I'),
        fits.Column(name='tpcf0', array=legendre_coef(tpcf, 0, config), format='D'),
        fits.Column(name='tpcf2', array=legendre_coef(tpcf, 2, config), format='D'),
        fits.Column(name='tpcf4', array=legendre_coef(tpcf, 4, config), format='D'),
        fits.Column(name='tpcf6', array=legendre_coef(tpcf, 6, config), format='D')
    ],
        name='legendre')
    hdu.writeto(config.__class__.__name__ + '_legendre.fits')


if __name__ == '__main__':
    start = time.time()
    legendre_to_file()
    end = time.time()
    print('---------------------------------')
    print(end - start)