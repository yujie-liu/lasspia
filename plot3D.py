import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import os.path as osp
from astropy.io import fits
from scipy.special import legendre


def rebin(a, *args):
    shape = a.shape
    lenShape = len(shape)
    factor = np.asarray(shape) / np.asarray(args)
    evList = ['a.reshape('] + \
             ['args[%d],factor[%d],' % (i, i) for i in range(lenShape)] + \
             [')'] + ['.mean(%d)' % (i + 1) for i in range(lenShape)]
    print(''.join(evList))
    eval(''.join(evList))


def legendre_coef(tpcf, l, s_vec):
    if l % 2 == 1:
        raise ValueError("The order of Legendre polynomials to be even")
    delta_mu = 2 / 32
    temp = np.copy(tpcf)
    nonzeros = [np.count_nonzero(e) for e in temp]
    nonzeros[nonzeros == 0] = 1
    nonzeros = np.asarray(nonzeros)[:, None]
    mu_vec = np.linspace(-1, 1, 32)[:, None]  # the vector of mu values
    p = legendre(l)  # Legendre polynomial at order l
    mu_vec = np.array([p(mu) for mu in mu_vec])  # P-_l(mu)
    temp = temp * mu_vec.transpose() * delta_mu  # P_l(mu) * tpcf(s, mu) * delta_mu
    temp = temp * (2 * l + 1) / 2
    print(temp.sum(axis=1).shape)
    return temp.sum(axis=1)


mpl.rc("font", family="serif", size=14)

file = osp.join('/home', 'yujie', 'Documents', 'data', 'corr_North_dr12_3D.dat')
data = np.genfromtxt(file)
pi, sigma, xi, dxi = [data[:, i] for i in range(0, 4)]
DD, DR, RR = [data[:, i] for i in range(4, 7)]

Pi, Sigma = np.meshgrid(pi[:32], pi[:32])
S = np.sqrt(Pi ** 2 + Sigma ** 2)
s_vec = np.sqrt(pi[:32]**2+pi[:32]**2)
print(s_vec)

Xi = xi.reshape(32, 32)

fig, [ax0, ax1] = plt.subplots(1, 2, figsize=(10, 4))
im = ax0.imshow((S ** 2 * Xi).T, origin='lower', interpolation='nearest',
                extent=[0, 200, 0, 200],
                vmin=0, vmax=200,
                cmap='magma')
ax0.set(aspect='equal',
        xlabel='$\sigma$',
        ylabel='$\pi$')
cb = fig.colorbar(im, ax=ax0, label=r'$s^2 \xi(s)$', fraction=0.05)

cf = ax1.contourf(Sigma, Pi, S ** 2 * Xi, 20,
                  vmin=0, vmax=200,
                  cmap='magma')
ax1.set(aspect='equal',
        xlabel='$\sigma$',
        ylabel='$\pi$',
        xlim=(0, 200),
        ylim=(0, 200))
cb = fig.colorbar(im, ax=ax1, label=r'$s^2 \xi(s)$', fraction=0.05)

fig.tight_layout()

plt.show()
Xi = S ** 2 * Xi
hdu = fits.BinTableHDU.from_columns([
    fits.Column(name='s', array=s_vec, format='I'),
    fits.Column(name='tpcf0', array=legendre_coef(Xi, 0, s_vec), format='D'),
    fits.Column(name='tpcf2', array=legendre_coef(Xi, 2, s_vec), format='D'),
    fits.Column(name='tpcf4', array=legendre_coef(Xi, 4, s_vec), format='D'),
    fits.Column(name='tpcf6', array=legendre_coef(Xi, 6, s_vec), format='D')
],
    name='legendre')
hdu.writeto('cute_legendre.fits')
