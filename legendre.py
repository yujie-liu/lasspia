import numpy as np
import sys
import os.path as osp
from astropy.io import fits
from scipy.special import legendre
import legendre_util as util


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
    print(hdul.info())
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
    tpcf = np.zeros((len(pi_vals), len(sigma_vals)))
    bin_space = pi_vals[1] - pi_vals[0]
    s_bins = int(len(pi_vals) / 2)
    # put tpcf values into (s, mu) grid
    # need vectorization for efficiency
    for tuple in hdul[5].data:

        isigma, ipi, rr, dr, dd, dde2 = tuple
        sigma = isigma * bin_space - s_bins * bin_space
        pi = ipi * bin_space - s_bins * bin_space
        s = max(1e-6, np.sqrt(sigma ** 2 + pi ** 2))
        mu = pi / s
        val = (dd / nDD - 2 * dr / nDR + rr / nRR) / (rr / nRR)
        val = s ** 2 * val
        if isigma > 5 and isigma < bin_s * 2 - 5:
            tpcf[isigma, ipi] = val
    leng = tpcf.shape
    tpcf = tpcf[int(leng[0] / 2):leng[0], 0:int(leng[1] / 2)]
    tpcf = np.fliplr(tpcf)
    tpcf = np.lib.pad(tpcf, ((int(leng[0] / 2), 0), (int(leng[1] / 2), 0)), 'reflect')
    print(tpcf.min(), tpcf.max(), tpcf.mean())
    return bin_space, bin_s, bin_theta, tpcf


def tpcf_to_s_mu(bin_space, s_bins, mu_bins, tpcf):
    """
    Convert tpcf(sigma, pi) to tpcf(s, mu)
    :param bin_space:
    :param s_bins: number of s bins
    :param mu_bins: number of mu bins
    :param tpcf: tpcf(sigma, pi) matrix
    :return:
        s_vec: s values
        tpcf_s_mu: tpcf(s, mu)
    """
    tpcf_s_mu = np.zeros((s_bins, mu_bins))
    s_vec = np.arange(0, s_bins * bin_space, bin_space)[:, None]
    for (isigma, ipi), val in np.ndenumerate(tpcf):
        sigma = isigma * bin_space - s_bins * bin_space
        pi = ipi * bin_space - s_bins * bin_space
        s = max(1e-6, np.sqrt(sigma ** 2 + pi ** 2))
        mu = pi / s
        i_s = int(s / bin_space)
        i_mu = int((mu + 1) * mu_bins / 2)
        if 5 < i_s < s_bins - 5 and i_mu < mu_bins:
            # val = val * s **2
            tpcf_s_mu[i_s, i_mu] = val
    tpcf_s_mu = interpolate(tpcf_s_mu)
    print(tpcf_s_mu.max(), tpcf_s_mu.min(), tpcf_s_mu.mean())
    print(s_vec.min(), s_vec.max(), s_vec.mean())
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
    print(hdul.info())

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
    # put tpcf values into (s, mu) grid
    # need vectorization for efficiency
    for tuple in hdul[5].data:
        isigma, ipi, rr, dr, dd, dde2 = tuple
        s = (sigma_vals[isigma] ** 2 + pi_vals[ipi] ** 2) ** .5
        mu = pi_vals[ipi] / s
        i_s = int(s / (range_s[1] / bin_s))
        i_mu = int((mu + 1) * bin_theta / 2)
        val = (dd / nDD - 2 * dr / nDR + rr / nRR) / (rr / nRR)
        # val = val * s**2
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

    tpcf2 = np.copy(tpcf)
    for i in range(1, bin_theta):
        temp = tpcf[:, i]
        zeros = np.nonzero(temp == 0)[0]
        nonzeros = np.nonzero(temp != 0)[0]
        nonzero_vals = temp[nonzeros]
        if len(nonzeros) >= 10:
            temp[temp == 0] = np.interp(zeros, nonzeros, nonzero_vals)
            tpcf2[:, i] = temp
    tpcf = .5 * (tpcf1 + tpcf2)
    # tpcf[tpcf > 2] = 2
    print(tpcf.min(), tpcf.max(), tpcf.mean())
    tpcf[np.where(tpcf1 * tpcf2 == 0)] = 0  # get rid of entries where tpcf1 or tpcf2 are zero
    return s_vec, tpcf


def interpolate(tpcf):
    """
    Using interpolation to fill in the empty entries in the tpcf matrix
    :param tpcf: ndarray
    :return: filled tpcf matrix
    """
    tpcf2 = np.copy(tpcf)
    for i in range(1, tpcf.shape[1]):
        temp = tpcf[:, i]
        zeros = np.nonzero(temp == 0)[0]
        nonzeros = np.nonzero(temp != 0)[0]
        nonzero_vals = temp[nonzeros]
        if len(nonzeros) >= 10:
            temp[temp == 0] = np.interp(zeros, nonzeros, nonzero_vals)
            tpcf2[:, i] = temp
    return tpcf2


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
        raise ValueError("\'l\' has to be even")
    delta_mu = 2 / config.binningTheta()['bins']
    temp = np.copy(tpcf)
    nonzeros = [np.count_nonzero(e) for e in temp]
    nonzeros[nonzeros == 0] = 1
    mu_vec = np.linspace(-1, 1, int(config.binningTheta()['bins']))[:, None]  # the vector of mu values
    p = legendre(l)  # Legendre polynomial at order l
    mu_vec = np.array([p(mu) for mu in mu_vec])  # P_l(mu)
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
    bin_space, bin_s, bin_theta, tpcf_pi_sigma = get_tpcf_sigma_pi(fits_file, config)
    s_vec, tpcf = tpcf_to_s_mu(bin_space, bin_s, bin_theta, tpcf_pi_sigma)
    tpcf_pi_sigma = np.transpose(tpcf_pi_sigma)

    util.plot_contourf(tpcf_pi_sigma, 1, 1.2, 10, [-300, 300, -300, 300])

    # plt.figure()
    # plt.contourf(tpcf, 50, extent=[-1, 1, 0, 300])
    # plt.colorbar()
    # plt.ylabel(r"$s $", fontsize=16)
    # plt.xlabel(r"$\mu $", fontsize=16)
    # plt.tight_layout()
    # plt.show()

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
    legendre_to_file()
