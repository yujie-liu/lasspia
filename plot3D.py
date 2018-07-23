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


def legendre_coef(Xi, l):
    if l % 2 == 1:
        raise ValueError("\'l\' has to be even")
    delta_mu = 2 / 32
    temp = np.copy(Xi)
    nonzeros = [np.count_nonzero(e) for e in temp]
    nonzeros[nonzeros == 0] = 1
    mu_vec = np.linspace(-1, 1, 32)[:, None]  # the vector of mu values
    p = legendre(l)  # Legendre polynomial at order l
    mu_vec = np.array([p(mu) for mu in mu_vec])  # P_l(mu)
    temp = temp * mu_vec.transpose() * delta_mu  # P_l(mu) * Xi(s, mu) * delta_mu
    temp = temp * (2 * l + 1) / 2
    return temp.sum(axis=1)


mpl.rc("font", family="serif", size=14)

file = osp.join('/home', 'yujie', 'Documents', 'data', 'corr_North_dr12_3D.dat')
data = np.genfromtxt(file)
pi, sigma, xi, dxi = [data[:, i] for i in range(0, 4)]
DD, DR, RR = [data[:, i] for i in range(4, 7)]

Pi, Sigma = np.meshgrid(pi[:32], pi[:32])
S = np.sqrt(Pi ** 2 + Sigma ** 2)
s_vec = np.sqrt(pi[:32]**2+pi[:32]**2)
mu_vec = pi[:32]/s_vec

Xi = xi.reshape(32, 32)

# leng = Xi.shape
# Xi = np.lib.pad(Xi, ((leng[0], 0), (leng[1], 0)), 'reflect')

print(Xi.min(), Xi.mean(), Xi.max())

levels = np.arange(-0.01, 0.03, 0.002)
plt.contourf(Xi.T, 20, origin='lower', interpolation='nearest', extent=[-200, 200, -200, 200], levels=levels)
plt.xlabel(r'$\sigma$')
plt.ylabel(r'$\pi$')
plt.colorbar()
plt.tight_layout()
plt.show()
fig, [ax0, ax1] = plt.subplots(1, 2, figsize=(10, 4))
im = ax0.imshow((S ** 2 * Xi).T, origin='lower', interpolation='nearest',
                extent=[0, 200, 0, 200],
                vmin=0, vmax=200,
                cmap='magma')
im = ax0.imshow((1 * Xi).T, origin='lower', interpolation='nearest',
                extent=[0, 200, 0, 200], vmax=0.03,
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


Xi = np.zeros((32,32))
for tuple in data:
    pi, sigma, xi, dxi =tuple[0:4]
    s = np.sqrt(pi ** 2 + sigma ** 2)
    mu=pi/s
    s_index=int(s*32/280)
    mu_index=int((mu+1)*32/2)
    Xi[s_index,mu_index] = xi

# Xi1 = np.copy(Xi)
# for i in range(1, 32):
#     temp = Xi[i, :]
#     zeros = np.nonzero(temp == 0)[0]
#     nonzeros = np.nonzero(temp != 0)[0]
#     nonzero_vals = temp[nonzeros]
#     if len(nonzeros) >= 10:
#         temp[temp == 0] = np.interp(zeros, nonzeros, nonzero_vals)
#         Xi1[i, :] = temp
# Xi = Xi1

a = np.arange(0, 2.5, 0.5)
b = np.arange(10).reshape(-1, 5)
plt.plot(a, b[0])
plt.show()

hdu = fits.BinTableHDU.from_columns([
    fits.Column(name='s', array=s_vec, format='I'),
    fits.Column(name='tpcf0', array=legendre_coef(Xi, 0)*s_vec**2, format='D'),
    fits.Column(name='tpcf2', array=legendre_coef(Xi, 2)*s_vec**2, format='D'),
    fits.Column(name='tpcf4', array=legendre_coef(Xi, 4)*s_vec**2, format='D'),
    fits.Column(name='tpcf6', array=legendre_coef(Xi, 6)*s_vec**2, format='D')
],
    name='legendre')
hdu.writeto('cute_legendre.fits')

