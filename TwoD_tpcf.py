from astropy.io import fits
import legendre as leg
from scipy.special import legendre
import os
import legendre_util
import matplotlib.pyplot as plt
from matplotlib import rc
import re
import numpy as np


def transform(tpcf, ratio):
    '''
    Transform the tpcf matrix with sigma multiplied by a scaler
    :param tpcf: tpcf matrix
    :param ratio: the number by which sigma is scaled
    :return: the new tpcf
    '''
    tpcf2 = np.zeros(tpcf.shape)
    for (isigma, ipi), val in np.ndenumerate(tpcf):
        new_isigma = int(isigma * ratio)
        if new_isigma < tpcf.shape[0]:
            tpcf2[new_isigma, ipi] = val
    tpcf2 = legendre_util.interpolate(tpcf2, axis=1)
    return tpcf2


def expand(tpcf, title='', uncertainty=False):
    '''
    :param tpcf: tpcf matrix
    :param title: one of 'CR', 'RR', 'LS'
    :param uncertainty: True: to expand the uncertainty terms; False: only to expand normal terms
    :return: a tuple of expansion results (s, order 0, 2, 4, 6)
    '''
    rc('font', family='serif')
    rc('font', size=16)
    s_vec, tpcf2 = leg.tpcf_to_s_mu(2, int(tpcf.shape[0] / 2), 10000, tpcf)
    tpcf2[np.isnan(tpcf2)] = 0
    if uncertainty:
        filename = 'tpcf_u_legendre.fits'
        try:
            os.remove(filename)
        except FileNotFoundError:
            pass

        hdu = fits.BinTableHDU.from_columns([
            fits.Column(name='s', array=s_vec, format='I'),
            fits.Column(name='tpcf0', array=leg.legendre_coef_u(tpcf2, 0, mu_bins=10000), format='D'),
            fits.Column(name='tpcf2', array=leg.legendre_coef_u(tpcf2, 2, mu_bins=10000), format='D'),
            fits.Column(name='tpcf4', array=leg.legendre_coef_u(tpcf2, 4, mu_bins=10000), format='D'),
            fits.Column(name='tpcf6', array=leg.legendre_coef_u(tpcf2, 6, mu_bins=10000), format='D')
        ],
            name='legendre')
    else:
        filename = 'tpcf_legendre.fits'
        try:
            os.remove(filename)
        except FileNotFoundError:
            pass
        hdu = fits.BinTableHDU.from_columns([
            fits.Column(name='s', array=s_vec, format='I'),
            fits.Column(name='tpcf0', array=leg.legendre_coef(tpcf2, 0, mu_bins=10000), format='D'),
            fits.Column(name='tpcf2', array=leg.legendre_coef(tpcf2, 2, mu_bins=10000), format='D'),
            fits.Column(name='tpcf4', array=leg.legendre_coef(tpcf2, 4, mu_bins=10000), format='D'),
            fits.Column(name='tpcf6', array=leg.legendre_coef(tpcf2, 6, mu_bins=10000), format='D')
        ],
            name='legendre')
    hdu.writeto(filename)
    hdul = fits.open(filename)
    coef = hdul[1].data
    return coef['s'], coef['tpcf0'], coef['tpcf2'], coef['tpcf4'], coef['tpcf6']


def find_min2(tpcf, r=1, l=2, step=.1):
    leng = tpcf.shape
    current = r
    l_index = int(l / 2) + 1
    for j in range(7):
        vals = []
        for i in np.arange(- step * 5, step * 5, step):
            tpcf1 = transform(tpcf, r + i)
            tpcf1 = np.lib.pad(tpcf1, ((leng[0], 0), (leng[1], 0)), 'reflect')
            coef = expand(tpcf1)
            vals.append(max(abs(coef[l_index].min()), abs(coef[l_index].max())))
        index = np.argmin(vals)
        print(index)
        r = r - step * 5 + step * (index)
        step = step / 5
        print(vals)
        print(r)


def find_min3(tpcf, r=1, delta_x=.1):
    leng = tpcf.shape
    vals = []
    step = .05
    for i in np.arange(-.3, .3, step):
        tpcf1 = transform(tpcf, r + i)
        tpcf1 = np.lib.pad(tpcf1, ((leng[0], 0), (leng[1], 0)), 'reflect')
        s_1, o0_1, o2_1, o4_1, o6_1 = expand(tpcf1)
        vals.append(max(abs(o2_1.min()), abs(o2_1.max())))
    index = np.argmin(vals)
    print(index)
    r = r - step * 6 + step * (index)
    print(vals)
    print(r)


def find_min(tpcf, r=1, delta_x=.1):
    '''
    Find the sigma-scaler that results in the smallest l=2 term
    Uses Newton's method of optimization
    When the calculated r is within 1e-4 range of the last calculated r, return r
    :param tpcf: tpcf matrix
    :param r: scaler to start with
    :param delta_x:
    :return:
    '''
    leng = tpcf.shape
    while True:
        tpcf1 = transform(tpcf, r)
        tpcf1 = np.lib.pad(tpcf1, ((leng[0], 0), (leng[1], 0)), 'reflect')
        tpcf2 = transform(tpcf, r + delta_x)
        tpcf2 = np.lib.pad(tpcf2, ((leng[0], 0), (leng[1], 0)), 'reflect')
        tpcf3 = transform(tpcf, r - delta_x)
        tpcf3 = np.lib.pad(tpcf3, ((leng[0], 0), (leng[1], 0)), 'reflect')
        s_1, o0_1, o2_1, o4_1, o6_1 = expand(tpcf1)
        v1 = max(abs(o2_1.min()), abs(o2_1.max()))
        s_2, o0_2, o2_2, o4_2, o6_2 = expand(tpcf2)
        v2 = max(abs(o2_2.min()), abs(o2_2.max()))
        s_3, o0_3, o2_3, o4_3, o6_3 = expand(tpcf3)
        v3 = max(abs(o2_3.min()), abs(o2_3.max()))

        diff1 = (v2 - v3) / delta_x / 2
        if abs(diff1) < 1e-5:
            break
        snd_1 = (v1 - v3) / delta_x
        snd_2 = (v2 - v1) / delta_x
        diff2 = (snd_2 - snd_1) / delta_x
        new_r = r - diff1 / diff2
        if abs(new_r - r) < 1e-3:
            break
        r = new_r
        print(diff1, diff2)
        print('Current r: ', r)
    print('R that minimizes l=2 term: ', r)
    return r


def find_min_plot(tpcf, name, r=1, step=.1):
    vals_arr = []
    r_range = np.arange(.5, 1.5, step)
    for l in range(0, 7, 2):
        index = int(l / 2) + 1
        leng = tpcf.shape
        vals = []
        for i in np.arange(-.5, .5, step):
            tpcf1 = transform(tpcf, r + i)
            tpcf1 = np.lib.pad(tpcf1, ((leng[0], 0), (leng[1], 0)), 'reflect')
            coef = expand(tpcf1)
            vals.append(max(abs(coef[index].min()), abs(coef[index].max())))
        vals_arr.append(vals)
    plt.figure(figsize=[8, 6])
    plt.plot(r_range, vals_arr[3])
    plt.plot(r_range, vals_arr[2])
    plt.plot(r_range, vals_arr[1])
    plt.plot(r_range, vals_arr[0])
    plt.legend(['l = 6', 'l = 4', 'l = 2', 'l = 0'], loc='upper left')
    plt.ylabel(r"Minimum of expansion results", fontsize=16)
    plt.xlabel(r"$\sigma $ ratio", fontsize=16)
    # plt.title(r'Minimum of Legendre expansion results with $s^2$ (' + name +')')
    plt.tight_layout()
    plt.savefig('/home/yujie/Desktop/Figures/Figure_0921-_CR.png')
    plt.show()


def plot():
    '''
    :return:
    '''
    names = ['CR', 'RR', 'LS']
    for r in np.arange(r_start, r_end, r_step):
        for name in names:
            if name == 'CR':
                tpcfCR = transform(tpcfCR1, r)
                tpcfCR_u = transform(tpcfCR_u1, r)
                tpcf = np.lib.pad(tpcfCR, ((leng[0], 0), (leng[1], 0)), 'reflect')
                tpcf_u = np.lib.pad(tpcfCR_u, ((leng[0], 0), (leng[1], 0)), 'reflect')
            elif name == 'RR':
                tpcfRR = transform(tpcfRR1, r)
                tpcfRR_u = transform(tpcfRR_u1, r)
                tpcf = np.lib.pad(tpcfRR, ((leng[0], 0), (leng[1], 0)), 'reflect')
                tpcf_u = np.lib.pad(tpcfRR_u, ((leng[0], 0), (leng[1], 0)), 'reflect')
            elif name == 'LS':
                tpcfLS = transform(tpcfLS1, r)
                tpcfLS_u = transform(tpcfLS_u1, r)
                tpcf = np.lib.pad(tpcfLS, ((leng[0], 0), (leng[1], 0)), 'reflect')
                tpcf_u = np.lib.pad(tpcfLS_u, ((leng[0], 0), (leng[1], 0)), 'reflect')
            else:
                raise ValueError('Must be one of CR, RR or LS')

            s, o0, o2, o4, o6 = expand(tpcf, title=name)
            o0min, o2min, o4min, o6min = abs(o0.min()), abs(o2.min()), abs(o4.min()), abs(o6.min())
            min_vals[name] = np.append(min_vals[name], [o0min, o2min, o4min, o6min])
            s, uo0, uo2, uo4, uo6 = expand(tpcf_u, title=name, uncertainty=True)

            extent = [-400, 400, -400, 400]
            if name == 'CR':
                max = 0.032
            elif name == 'RR':
                max = 0.016
            elif name == 'LS':
                max = 0.8
            plt.figure(figsize=[14, 6])
            plt.subplot(1, 2, 2)
            plt.plot(s, o6)
            plt.fill_between(s, o6 - uo6, o6 + uo6, alpha=0.2)
            plt.plot(s, o4)
            plt.fill_between(s, o4 - uo4, o4 + uo4, alpha=0.2)
            plt.plot(s, o2)
            plt.fill_between(s, o2 - uo2, o2 + uo2, alpha=0.2)
            plt.plot(s, o0)
            plt.fill_between(s, o0 - uo0, o0 + uo0, alpha=0.2)
            plt.tight_layout()
            if name == 'CR':
                plt.ylim(-300, 400)
            elif name == 'RR':
                plt.ylim(-80, 100)
            elif name == 'LS':
                plt.ylim(-100, 125)
            plt.xlabel('s')
            plt.ylabel(r'$\widetilde{\xi}(s) s^2$')
            r_format = str(float("{0:.2f}".format(r)))
            plt.title('Legendre expansion of tpcf' + name + '\n with uncertainties \n ' r'$\sigma$ ratio = ' + r_format)
            plt.legend(['l = 6', 'l = 4', 'l = 2', 'l = 0'], loc='upper left')

            plt.subplot(1, 2, 1)
            plt.ylabel(r"$\pi $", fontsize=16)
            plt.xlabel(r"$\sigma $", fontsize=16)
            plt.contourf(tpcf.T, origin='lower', extent=extent, levels=np.arange(0, max, max / 10))
            plt.colorbar()
            plt.tight_layout()
            plt.title('tpcf' + name)
            plt.close('all')

    for name in ['CR', 'RR', 'LS']:
        v = min_vals[name].reshape(-1, 4)
        plt.figure(figsize=[7, 5])
        plt.plot(np.arange(r_start, r_end, r_step), v[:, 3])
        plt.plot(np.arange(r_start, r_end, r_step), v[:, 2])
        plt.plot(np.arange(r_start, r_end, r_step), v[:, 1])
        plt.plot(np.arange(r_start, r_end, r_step), v[:, 0])
        plt.legend(['l = 6', 'l = 4', 'l = 2', 'l = 0'], loc='upper left')
        plt.ylabel(r"Minimum of expansion results", fontsize=16)
        plt.xlabel(r"$\sigma $ ratio", fontsize=16)
        plt.title('Minimum of Legendre expansion results (' + name + ')')
        plt.tight_layout()
        plt.show()


def tpcf_results(fits_file, bin_space=2, bin_s=100, dd_scale=False, dr_scale=False, rr_scale=False):
    hdul = fits.open(fits_file)
    tpcf = np.zeros((bin_s, bin_s))
    tpcf_s2 = np.zeros((bin_s, bin_s))
    dd_matrix = np.zeros((bin_s, bin_s))
    dr_matrix = np.zeros((bin_s, bin_s))
    rr_matrix = np.zeros((bin_s, bin_s))
    nDD = 0
    nRR = 0
    nDR = 0
    for tuple in hdul[1].data:
        isigma, ipi, rr, dr, dd, dde2 = tuple
        dd = dd / 2
        sigma = isigma * bin_space - bin_s * bin_space
        pi = ipi * bin_space - bin_s * bin_space
        s = max(1e-6, np.sqrt(sigma ** 2 + pi ** 2))
        nDD = nDD + dd
        nDR = nDR + dr
        nRR = nRR + rr

    leng = tpcf.shape

    dd_matrix = dd_matrix[int(leng[0] / 2):leng[0], int(leng[1] / 2):leng[1]] / sum(dd_matrix)
    dd_matrix = np.lib.pad(dd_matrix, ((int(leng[0] / 2), 0), (int(leng[1] / 2), 0)), 'reflect')
    dr_matrix = dr_matrix[int(leng[0] / 2):leng[0], int(leng[1] / 2):leng[1]] / sum(dr_matrix)
    dr_matrix = np.lib.pad(dr_matrix, ((int(leng[0] / 2), 0), (int(leng[1] / 2), 0)), 'reflect')
    rr_matrix = rr_matrix[int(leng[0] / 2):leng[0], int(leng[1] / 2):leng[1]] / sum(rr_matrix)
    rr_matrix = np.lib.pad(rr_matrix, ((int(leng[0] / 2), 0), (int(leng[1] / 2), 0)), 'reflect')

    if dd_scale:
        dd_matrix = dd_matrix * 2
    if dr_scale:
        dr_matrix = dr_matrix * 2
    if rr_scale:
        rr_matrix = rr_matrix * 2
    for i_s in range(bin_s * 2):
        for i_p in range(bin_s * 2):
            sigma = i_s * bin_space - bin_s * bin_space
            pi = i_p * bin_space - bin_s * bin_space
            s = max(1e-6, np.sqrt(sigma ** 2 + pi ** 2))
            if s > 300:
                dd_matrix[i_s, i_p] = 0
                dr_matrix[i_s, i_p] = 0
                rr_matrix[i_s, i_p] = 0

    tpcf1 = (dd_matrix - 2 * dr_matrix + rr_matrix) / rr_matrix
    tpcf1[np.isnan(tpcf1)] = 0
    for i_s in range(bin_s * 2):
        for i_p in range(bin_s * 2):
            sigma = i_s * bin_space - bin_s * bin_space
            pi = i_p * bin_space - bin_s * bin_space
            s = max(1e-6, np.sqrt(sigma ** 2 + pi ** 2))
            if s < 300:
                tpcf_s2[i_s, i_p] = s ** 2 * tpcf1[i_s, i_p]
    return tpcf1, tpcf_s2


if __name__ == '__main__':
    from argparse import ArgumentParser
    import os.path as osp

    parser = ArgumentParser(description="Large Scale Structure Probability Integration Algorithm")
    parser.add_argument('txtFile', metavar='txtFile', type=str, nargs=1,
                        help='the txt file for mocks.')
    parser.add_argument('type', metavar='type', type=str, nargs=1,
                        help='one of \'CR\', \'RR\' and \'LS\' ')
    parser.add_argument('-p', action='store_true', help='Set this flag to show tpcf plot')
    parser.add_argument('-dd', action='store_true',
                        help='Set this flag to scale DD by a factor of 2 in calculation of tpcf')
    parser.add_argument('-dr', action='store_true',
                        help='Set this flag to scale DR by a factor of 2 in calculation of tpcf')
    parser.add_argument('-rr', action='store_true',
                        help='Set this flag to scale RR by a factor of 2 in calculation of tpcf')
    args = parser.parse_args()
    data_file = osp.join(osp.abspath('../data'), args.txtFile[0])
    # tpcf = tpcf_results(data_file)
    # find_min_plot(tpcf)
    with open(data_file) as f:
        content = f.read().splitlines()
    max_sigma = int(content[len(content) - 1].split()[0])
    min_sigma = int(content[1].split()[0])
    i = len(content) - 1
    while int(content[i].split()[0]) == max_sigma:
        i = i - 1
    bin_size = max_sigma - int(content[i].split()[0])
    size = int((max_sigma - min_sigma) / bin_size + 1) * 2
    content = content[1:]
    tpcf = np.zeros((size, size))
    tpcfCR1 = np.zeros((size, size))
    tpcfRR1 = np.zeros((size, size))
    tpcfLS1 = np.zeros((size, size))
    tpcfCR_u1 = np.zeros((size, size))
    tpcfRR_u1 = np.zeros((size, size))
    tpcfLS_u1 = np.zeros((size, size))
    for e in content:
        tuple = e.split()
        i_sigma = int(tuple[0])
        i_sigma = int(i_sigma / 2)
        i_pi = int(tuple[1])
        i_pi = int(i_pi / 2)
        tpcfCR1[i_sigma, i_pi] = float(tuple[2])
        tpcfRR1[i_sigma, i_pi] = float(tuple[4])
        tpcfLS1[i_sigma, i_pi] = float(tuple[6])
        tpcfCR_u1[i_sigma, i_pi] = float(tuple[3])
        tpcfRR_u1[i_sigma, i_pi] = float(tuple[5])
        tpcfLS_u1[i_sigma, i_pi] = float(tuple[7])

    leng = tpcfCR_u1.shape
    r_start, r_end, r_step = 0.4, 2.0, 0.05
    min_vals = {
        'CR': [],
        'RR': [],
        'LS': []
    }
    if args.type[0] == 'CR':
        tpcf = tpcfCR1
    elif args.type[0] == 'RR':
        tpcf = tpcfRR1
    elif args.type[0] == 'LS':
        tpcf = tpcfLS1

    # find_min2(tpcf)
    find_min_plot(tpcf, args.type[0])
