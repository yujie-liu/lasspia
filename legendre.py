import numpy as np
import sys
import os.path as osp
from astropy.io import fits
from scipy.special import legendre
import matplotlib.pyplot as plt
from matplotlib import rc
import legendre_util


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


def get_tpcf_sigma_pi(fits_file, config):
    """
    Calculate tpcf(sigma, pi) from integration.fits file for plotting purpose
    :param fits_file: filename
    :param config: configuration object
    :return:
        bin_space:
        bin_s: number of s bins
        bin_theta: number of theta bins
        tpcf: tpcf(sigma, pi) matrix
    """
    hdul = fits.open(fits_file)
    # normalization factors from integration.py
    nRR, nDR, nDD = (lambda h:
                     (h['normrr'],
                      h['normdr'],
                      h['normdd']))(fits.getheader(fits_file, 'TPCF2D'))
    try:
        sigma_vals = hdul[3].data['binCenter']
        pi_vals = hdul[4].data['binCenter']
    except:
        sigma_vals = hdul[2].data['binCenter']
        pi_vals = hdul[3].data['binCenter']
    bin_theta, range_theta = config.binningTheta()['bins'], config.binningTheta()['range']
    bin_theta = int(bin_theta)
    bin_s, range_s = config.binningS()['bins'], config.binningS()['range']
    bin_s = int(bin_s)
    tpcf = np.zeros((len(pi_vals), len(sigma_vals)))
    tpcf_s2 = np.zeros((len(pi_vals), len(sigma_vals)))
    dd_matrix = np.zeros((len(pi_vals), len(sigma_vals)))
    dr_matrix = np.zeros((len(pi_vals), len(sigma_vals)))
    rr_matrix = np.zeros((len(pi_vals), len(sigma_vals)))
    bin_space = pi_vals[1] - pi_vals[0]
    s_bins = int(len(pi_vals) / 2)
    '''
    put tpcf values into (s, mu) grid
    need vectorization for efficiency
    '''
    for tuple in hdul[5].data:
        isigma, ipi, rr, dr, dd, dde2 = tuple
        dd = dd / 2
        sigma = isigma * bin_space - s_bins * bin_space
        pi = ipi * bin_space - s_bins * bin_space
        s = max(1e-6, np.sqrt(sigma ** 2 + pi ** 2))
        if (0 < isigma < bin_s * 2) and (0 <= ipi < bin_s * 2):
            dd_matrix[isigma, ipi] = dd / nDD
            dr_matrix[isigma, ipi] = dr / nDR
            rr_matrix[isigma, ipi] = rr / nRR

    leng = tpcf.shape

    dd_matrix = dd_matrix[int(leng[0] / 2):leng[0], int(leng[1] / 2):leng[1]]
    dd_matrix = np.lib.pad(dd_matrix, ((int(leng[0] / 2), 0), (int(leng[1] / 2), 0)), 'reflect')
    dr_matrix = dr_matrix[int(leng[0] / 2):leng[0], int(leng[1] / 2):leng[1]]
    dr_matrix = np.lib.pad(dr_matrix, ((int(leng[0] / 2), 0), (int(leng[1] / 2), 0)), 'reflect')
    rr_matrix = rr_matrix[int(leng[0] / 2):leng[0], int(leng[1] / 2):leng[1]]
    rr_matrix = np.lib.pad(rr_matrix, ((int(leng[0] / 2), 0), (int(leng[1] / 2), 0)), 'reflect')

    dd_matrix = dd_matrix * 2

    for i_s in range(bin_s * 2):
        for i_p in range(bin_s * 2):
            sigma = i_s * bin_space - s_bins * bin_space
            pi = i_p * bin_space - s_bins * bin_space
            s = max(1e-6, np.sqrt(sigma ** 2 + pi ** 2))
            if s > 300:
                dd_matrix[i_s, i_p] = 0
                dr_matrix[i_s, i_p] = 0
                rr_matrix[i_s, i_p] = 0

    tpcf1 = (dd_matrix - 2 * dr_matrix + rr_matrix) / rr_matrix
    tpcf1[np.isnan(tpcf1)] = 0
    for i_s in range(bin_s * 2):
        for i_p in range(bin_s * 2):
            sigma = i_s * bin_space - s_bins * bin_space
            pi = i_p * bin_space - s_bins * bin_space
            s = max(1e-6, np.sqrt(sigma ** 2 + pi ** 2))
            if s < 300:
                tpcf_s2[i_s, i_p] = s ** 2 * tpcf1[i_s, i_p]

    return bin_space, bin_s, bin_theta, tpcf1, tpcf_s2


def tpcf_to_s_mu(bin_space, s_bins, mu_bins, tpcf, ratio=1.0):
    """
    Convert tpcf(sigma, pi) to tpcf(s, mu)
    :param bin_space:
    :param s_bins: number of s bins
    :param mu_bins: number of mu bins
    :param tpcf: tpcf(sigma, pi) matrix
    :param s_sqr: whether to multiply by s^2
    :return:
        s_vec: s values
        tpcf_s_mu: tpcf(s, mu)
    """
    tpcf_s_mu = np.zeros((s_bins, mu_bins))
    s_vec = np.arange(0, s_bins * bin_space, bin_space)[:, None]
    for (isigma, ipi), val in np.ndenumerate(tpcf):
        sigma = isigma * bin_space - s_bins * bin_space
        pi = ipi * bin_space - s_bins * bin_space
        s = max(1e-6, np.sqrt((ratio * sigma) ** 2 + pi ** 2))
        mu = pi / s
        i_s = int(s / bin_space)
        i_mu = int((mu + 1) * mu_bins / 2)
        if 1 < i_s < s_bins - 1 and i_mu < mu_bins:
            tpcf_s_mu[i_s, i_mu] = val
    tpcf_s_mu = legendre_util.interpolate(tpcf_s_mu)
    return s_vec, tpcf_s_mu


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

    # normalization factors from integration.py
    nRR, nDR, nDD = (lambda h:
                     (h['normrr'],
                      h['normdr'],
                      h['normdd']))(fits.getheader(fits_file, 'TPCF'))
    sigma_vals = hdul[2].data['binCenter']
    pi_vals = hdul[3].data['binCenter']
    s_vec = hdul[1].data['binCenter']
    bin_theta, range_theta = config.binningTheta()['bins'], config.binningTheta()['range']
    bin_theta = int(bin_theta)
    bin_s, range_s = config.binningS()['bins'], config.binningS()['range']
    bin_s = int(bin_s)
    tpcf = np.zeros((bin_s, bin_theta))
    repeat = 0
    repeat_arr = []
    '''
    put tpcf values into (s, mu) grid
    need vectorization for efficiency
    '''
    for tuple in hdul[5].data:
        isigma, ipi, rr, dr, dd, dde2 = tuple
        s = (sigma_vals[isigma] ** 2 + pi_vals[ipi] ** 2) ** .5
        mu = pi_vals[ipi] / s
        i_s = int(s / (range_s[1] / bin_s))
        i_mu = int((mu + 1) * bin_theta / 2)
        val = (dd / nDD - 2 * dr / nDR + rr / nRR) / (rr / nRR)
        if 2 < s < int(range_s[1] - 2) and i_mu <= bin_theta - 1:
            # record number of overwritten entries
            if tpcf[i_s, i_mu] != 0:
                repeat = repeat + 1
                repeat_arr.append(abs(val - tpcf[i_s, i_mu]))
            tpcf[i_s, i_mu] = val

    # interpolation for tpcf(s, mu)
    tpcf1 = np.copy(tpcf)
    for i in range(1, bin_s):
        temp = tpcf[i, :]
        zeros = np.nonzero(temp == 0)[0]
        nonzeros = np.nonzero(temp != 0)[0]
        nonzero_vals = temp[nonzeros]
        if len(nonzeros) >= 10:
            temp[temp == 0] = np.interp(zeros, nonzeros, nonzero_vals)
            tpcf1[i, :] = temp

    tpcf2 = legendre_util.interpolate(tpcf)
    tpcf = .5 * (tpcf1 + tpcf2)
    tpcf[np.where(tpcf1 * tpcf2 == 0)] = 0  # get rid of entries where tpcf1 or tpcf2 are zero
    return s_vec, tpcf


def legendre_coef(tpcf, l, config=None, mu_bins=0):
    """
    Calculate the coefficient of tpcf_l(s) in Legendre expansion for all s
    Inputs:
    + tpcf: ndarray
        the tpcf as a matrix dimensioned by mu and s
    + l: the order of desired Legendre polynomial
    + config: Configuration object
    + mu_bins: if no configuration specified, the number of mu bins
    Outputs:
    + coef:
        Return the coefficient of tpcf_l(s) as an array indexed by s values
    + write to the *_legendre.fits file
        The .fits file has columns (s, tpcf0, etpcf0, tpcf2, etpcf2, ...)
    """
    if l % 2 == 1:
        raise ValueError("\'l\' has to be even")
    if config:
        mu_bins = int(config.binningTheta()['bins'])
    elif mu_bins == 0:
        raise ValueError('Either \'config\' or \'mu_bins\' has to be specified')
    delta_mu = 2 / mu_bins
    temp = np.copy(tpcf)
    nonzeros = [np.count_nonzero(e) for e in temp]
    nonzeros[nonzeros == 0] = 1
    mu_vec = np.linspace(-1, 1, mu_bins)[:, None]  # the vector of mu values
    p = legendre(l)  # Legendre polynomial at order l
    mu_vec = np.array([p(mu) for mu in mu_vec])  # P_l(mu)
    temp = temp * mu_vec.transpose() * delta_mu  # P_l(mu) * tpcf(s, mu) * delta_mu
    temp = temp * (2 * l + 1) / 2
    val = temp.sum(axis=1)
    return val


def legendre_coef_u(tpcf, l, config=None, mu_bins=0):
    '''
    Calculate the Legendre expansion results for uncertainty
    :param tpcf:
    :param l: the order of desired Legendre polynomial
    :param config: Configuration object
    :param mu_bins: if not configuration specified, the number of mu bins
    :return:
        Return the coefficient of tpcf_l(s) as an array indexed by s values
        Write results to *_legendre.fits file
    '''
    if l % 2 == 1:
        raise ValueError("\'l\' has to be even")
    if config:
        mu_bins = int(config.binningTheta()['bins'])
    elif mu_bins == 0:
        raise ValueError('Either \'config\' or \'mu_bins\' has to be specified')
    delta_mu = 2 / mu_bins
    temp = np.copy(tpcf)
    mu_vec = np.linspace(-1, 1, mu_bins)[:, None]  # the vector of mu values
    p = legendre(l)  # Legendre polynomial at order l
    mu_vec = np.array([p(mu) for mu in mu_vec])  # P_l(mu)
    temp = temp * mu_vec.transpose() * delta_mu  # P_l(mu) * tpcf(s, mu) * delta_mu
    temp = temp ** 2
    val = np.sqrt(temp.sum(axis=1))
    val = val * (2 * l + 1) / 2
    return val


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
    rc('font', family='serif')
    rc('font', size=16)
    plt.tight_layout()

    bin_space, bin_s, bin_theta, tpcf_pi_sigma, tpcf_pi_sigma_s2 = get_tpcf_sigma_pi(fits_file, config)
    s_vec, tpcf_s2 = tpcf_to_s_mu(bin_space, bin_s, bin_theta, tpcf_pi_sigma_s2)
    s_vec, tpcf = tpcf_to_s_mu(bin_space, bin_s, bin_theta, tpcf_pi_sigma)
    plt.ylabel(r"$\pi $", fontsize=16)
    plt.xlabel(r"$\sigma $", fontsize=16)
    hdu = fits.BinTableHDU.from_columns([
        fits.Column(name='s', array=s_vec, format='I'),
        fits.Column(name='tpcf0', array=legendre_coef(tpcf_s2, 0, config=config), format='D'),
        fits.Column(name='tpcf2', array=legendre_coef(tpcf_s2, 2, config=config), format='D'),
        fits.Column(name='tpcf4', array=legendre_coef(tpcf_s2, 4, config=config), format='D'),
        fits.Column(name='tpcf6', array=legendre_coef(tpcf_s2, 6, config=config), format='D')
    ],
        name='legendre')
    hdu.writeto(config.__class__.__name__ + '_legendre.fits')


if __name__ == '__main__':
    legendre_to_file()
