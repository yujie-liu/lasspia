from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np

if __name__ == '__main__':
    hdul = fits.open('cmassS_coarse_legendre.fits')
    coef = hdul[1].data
    print(coef.columns)
    print(coef)
    plt.plot(coef['s'], coef['tpcf6'])
    plt.plot(coef['s'], coef['tpcf4'])
    plt.plot(coef['s'], coef['tpcf2'])
    plt.plot(coef['s'], coef['tpcf0'])
    plt.legend(['6', '4', '2', '0'], loc='upper left')
    plt.show()
