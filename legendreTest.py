from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib import rc

if __name__ == '__main__':
    hdul = fits.open('cute_legendre.fits')
    coef = hdul[1].data
    rc('font', family='serif')
    rc('font', size=16)
    plt.figure(1)
    # plt.subplot(311)
    # plt.xlim([0, 200])
    # plt.ylim([-0.05, 0.2])

    plt.plot(coef['s'], coef['tpcf6'])
    plt.plot(coef['s'], coef['tpcf4'])
    plt.plot(coef['s'], coef['tpcf2'])
    plt.plot(coef['s'], coef['tpcf0'], color="r")
    # plt.legend(['l = 0'], loc='upper left')
    plt.legend(['l = 6', 'l = 4', 'l = 2', 'l = 0'], loc='upper left')

    plt.xlabel('s')
    plt.ylabel(r'$\widetilde{\xi}(s)$')
    plt.title(r'$\widetilde{\xi_l}(s) = \frac{2l+1}{2}\sum_j{\xi(s, \mu_j)P_l(\mu_j)\Delta \mu_j}$')

    plt.show()
