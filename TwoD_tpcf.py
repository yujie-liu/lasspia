from astropy.io import fits
import legendre as leg
from scipy.special import legendre
import os
import legendre_util as util
import matplotlib.pyplot as plt
from matplotlib import rc
import re
import numpy as np


def transform(tpcf, ratio):
    tpcf2 = np.zeros(tpcf.shape)
    for (isigma, ipi), val in np.ndenumerate(tpcf):
        new_isigma = int(isigma * ratio)
        if new_isigma < tpcf.shape[0]:
            tpcf2[new_isigma, ipi] = val
    tpcf2 = leg.interpolate(tpcf2, axis=1)
    return tpcf2


def expand(tpcf, title, uncertainty=False, ratio=1):
    rc('font', family='serif')
    rc('font', size=16)
    s_vec, tpcf2 = leg.tpcf_to_s_mu(2, 200, 10000, tpcf)
    tpcf2[np.isnan(tpcf2)] = 0
    extent = [-400, 400, -400, 400]
    if title == 'CR':
        max=0.032
    elif title == 'RR':
        max=0.016
    elif title == 'LS':
        max=0.8
    if uncertainty:
        pass
        #plt.title('tpcf' + title + ' with uncertainties')
        #plt.savefig('/home/yujie/Desktop/Figures/Figure_0714_tpcfv14_'+title+'_u_matrix.png')
    else:
        # plt.ylabel(r"$\pi $", fontsize=16)
        # plt.xlabel(r"$\sigma $", fontsize=16)
        # plt.contourf(tpcf.T, origin='lower', extent=extent, vmax=max)
        # plt.colorbar()
        # plt.tight_layout()
        # plt.title('tpcf' + title)
        # plt.savefig('/home/yujie/Desktop/TPCF_Testing/' + title + '/0714_tpcfv14matrix_' + name + str(r).replace(".", "") + '.png')
        #plt.show()
        #plt.close()
        pass
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
    print(hdul.info())
    coef = hdul[1].data
    # plt.figure(figsize=[7, 6])
    #
    # plt.plot(coef['s'], coef['tpcf6'])
    # plt.plot(coef['s'], coef['tpcf4'])
    # plt.plot(coef['s'], coef['tpcf2'])
    # plt.plot(coef['s'], coef['tpcf0'], color="r")
    # # plt.legend(['l = 0'], loc='upper left')
    # plt.legend(['l = 6', 'l = 4', 'l = 2', 'l = 0'], loc='upper left')
    # # plt.title(r'$\widetilde{\xi_l}(s) = \frac{2l+1}{2}\sum_j{\xi(s, \mu_j)P_l(\mu_j)\Delta \mu_j}$')
    #
    # if uncertainty:
    #     plt.xlabel('s')
    #     plt.ylabel(r'$\Delta \widetilde{\xi}(s) s^2$')
    #     plt.title('Legendre expansion of tpcf'+title+' uncertainties')
    #     plt.tight_layout()
    #     #plt.savefig('/home/yujie/Desktop/Figures/Figure_0714_tpcfv14_'+title+'_u.png')
    #     pass
    # else:
    #     plt.xlabel('s')
    #     plt.ylabel(r'$\widetilde{\xi}(s) s^2$')
    #     plt.title('Legendre expansion of tpcf'+title)
    #     plt.tight_layout()
    #     plt.savefig('/home/yujie/Desktop/Figures/Figure_0714_tpcfv14_'+title+'.png')
    #     plt.show()
    return coef['s'], coef['tpcf0'], coef['tpcf2'], coef['tpcf4'], coef['tpcf6']


if __name__ == '__main__':
    with open('TwoD_tpcf_v14.txt') as f:
        content = f.read().splitlines()

    content = content[1:]
    tpcfCR1 = np.zeros((200, 200))
    tpcfRR1 = np.zeros((200, 200))
    tpcfLS1 = np.zeros((200, 200))
    tpcfCR_u1 = np.zeros((200, 200))
    tpcfRR_u1 = np.zeros((200, 200))
    tpcfLS_u1 = np.zeros((200, 200))
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
        'CR':[],
        'RR':[],
        'LS':[]
    }
    names = ['CR', 'RR', 'LS']
    for r in np.arange(r_start, r_end, r_step):
        for name in ['CR', 'RR', 'LS']:
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

            s, o0, o2, o4, o6 = expand(tpcf, name, ratio=r)
            o0min, o2min, o4min, o6min = abs(o0.min()), abs(o2.min()), abs(o4.min()), abs(o6.min())
            min_vals[name]=np.append(min_vals[name], [o0min, o2min, o4min, o6min])
            s, uo0, uo2, uo4, uo6 = expand(tpcf_u, name, uncertainty=True, ratio=r)

            extent = [-400, 400, -400, 400]
            if name == 'CR':
                max = 0.032
            elif name == 'RR':
                max = 0.016
            elif name == 'LS':
                max = 0.8
            plt.figure(figsize=[14,6])
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
            plt.contourf(tpcf.T, origin='lower', extent=extent, levels=np.arange(0, max, max/10))
            plt.colorbar()
            plt.tight_layout()
            plt.title('tpcf' + name)


            #plt.savefig('/home/yujie/Desktop/TPCF_Testing/v14/'+name+'/0714_tpcfv14_'+name+str(r).replace(".", "")+'.png')
            #plt.show()
            plt.close('all')
    for name in ['CR', 'RR', 'LS']:
        v = min_vals[name].reshape(-1, 4)
        print(v)
        plt.figure(figsize=[7, 5])
        plt.plot(np.arange(r_start, r_end, r_step), v[:,3])
        plt.plot(np.arange(r_start, r_end, r_step), v[:,2])
        plt.plot(np.arange(r_start, r_end, r_step), v[:,1])
        plt.plot(np.arange(r_start, r_end, r_step), v[:,0])
        plt.legend(['l = 6', 'l = 4', 'l = 2', 'l = 0'], loc='upper left')
        plt.ylabel(r"Minimum of expansion results", fontsize=16)
        plt.xlabel(r"$\sigma $ ratio", fontsize=16)
        plt.title('Minimum of Legendre expansion results (' + name+')')
        plt.tight_layout()
        plt.savefig('/home/yujie/Desktop/TPCF_Testing/v14/0717_v14_min_'+name+'.png')
        plt.show()