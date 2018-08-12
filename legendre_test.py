from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib import rc

'''
Plot the *_legendre.fits file for testing
'''

def plot():
    from argparse import ArgumentParser
    parser = ArgumentParser(description="Large Scale Structure Probability Integration Algorithm")
    parser.add_argument('fitsFile', metavar='fitsFile', type=str, nargs=1,
                        help='FITS file.')
    args = parser.parse_args()
    fits_file = args.fitsFile[0]
    hdul = fits.open(fits_file)
    print(hdul.info())
    coef = hdul[1].data
    rc('font', family='serif')
    rc('font', size=16)

    plt.figure(figsize=[7, 6])
    plt.plot(coef['s'], coef['tpcf6'])
    plt.plot(coef['s'], coef['tpcf4'])
    plt.plot(coef['s'], coef['tpcf2'])
    plt.plot(coef['s'], coef['tpcf0'], color="r")
    plt.legend(['l = 6', 'l = 4', 'l = 2', 'l = 0'], loc='upper left')
    plt.xlabel('s')
    plt.ylabel(r'$\widetilde{\xi}(s) s^2$')

    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    plot()