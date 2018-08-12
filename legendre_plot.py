import sys
import os.path as osp
from matplotlib import rc
import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
from argparse import ArgumentParser
import skimage.measure


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


def tpcf_matrix(fits_file):
    hdul = fits.open(fits_file)
    # normalization factors from integration.py
    nRR, nDR, nDD = (lambda h:
                     (h['normrr'],
                      h['normdr'],
                      h['normdd']))(fits.getheader(fits_file, 'TPCF'))
    bins = len(hdul[1].data['binCenter']) * 2
    tpcf1d = np.zeros((bins, bins))
    for tuple in hdul[5].data:
        isigma, ipi, rr, dr, dd, dde2 = tuple
        val = (dd / nDD - 2 * dr / nDR + rr / nRR) / (rr / nRR)
        tpcf1d[isigma - 1, ipi - 1] = val
    return tpcf1d


def plot():
    parser = ArgumentParser(description="Large Scale Structure Probability Integration Algorithm")
    parser.add_argument('configFile', metavar='configFile', type=str, nargs=1,
                        help='A python file containing a subclass of lasspia.configuration of the same name.')
    parser.add_argument('fitsFile', metavar='fitsFile', type=str, nargs=1,
                        help='FITS file.')
    args = parser.parse_args()
    config = getInstance(args.configFile)
    fits_file = osp.join(osp.abspath('../data'), args.fitsFile[0])
    tpcf = tpcf_matrix(fits_file)
    tpcf = skimage.measure.block_reduce(tpcf, (4, 4), np.mean)
    leng = tpcf.shape
    tpcf = tpcf[int(leng[0] / 2):leng[0], 0:int(leng[1] / 2)]
    tpcf1 = np.arcsinh(tpcf)
    tpcf2 = np.log(tpcf)
    rc('font', family='serif')
    rc('font', size=16)
    extent = [-200, 0, 0, 200]
    plt.ylabel(r"$\pi $", fontsize=16)
    plt.xlabel(r"$\sigma $", fontsize=16)
    levels = [0, .6, .8, .9, 0.95, 1, 1.1, 1.2]
    plt.gca().invert_xaxis()
    # plt.contour(tpcf1, levels = levels, colors=('k',), extent=extent)
    plt.contourf(tpcf1, levels=levels, extent=extent)
    plt.title(r'$\xi (\pi, \sigma)$ in $ sinh^{-1} $ scale')
    plt.colorbar()
    plt.show()

    plt.ylabel(r"$\pi $", fontsize=16)
    plt.xlabel(r"$\sigma $", fontsize=16)
    levels = [-1, -.5, -.1, 0, .05, 0.4, 0.7, 1, 1.3, 1.6, 1.9, 2.2]
    plt.gca().invert_xaxis()
    # plt.contour(tpcf2, colors=('k',), levels=levels, extent=extent)
    plt.contourf(tpcf2, levels=levels, extent=extent)
    plt.title(r'$\xi (\pi, \sigma)$ in log scale')
    plt.colorbar()
    plt.show()

    plt.ylabel(r"$\pi $", fontsize=16)
    plt.xlabel(r"$\sigma $", fontsize=16)
    plt.imshow(tpcf, extent=[-200, 0, 0, 200], origin='lower')
    plt.title(r'$\xi (\pi, \sigma)$ unscaled')
    plt.colorbar()
    plt.show()

    levels = [.5, .8, .9, 1, 1.05, 1.1, 1.2, 2, 2.5, 3]
    plt.gca().invert_xaxis()
    plt.ylabel(r"$\pi $", fontsize=16)
    plt.xlabel(r"$\sigma $", fontsize=16)
    # plt.contour(tpcf, colors=('k',), levels=levels, extent=[-200, 0, 0, 200])
    plt.contourf(tpcf, levels=levels, extent=[-200, 0, 0, 200])
    plt.title(r'$\xi (\pi, \sigma)$ unscaled')
    plt.colorbar()
    plt.show()


if __name__ == '__main__':
    plot()
