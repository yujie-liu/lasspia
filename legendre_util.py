import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc


def gaussian(x, mu, sig):
    """
    Helper function of the bao_signal
    :param x:
    :param mu: mean
    :param sig: standard deviation
    :return:
    """
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))


def bao_signal(tpcf, space, radius=100, noise=5, e=0):
    """
    Create a circular/elliptical bao signal in the given tpcf(sigma, pi) matrix
    :param tpcf:
    :param space: bin space
    :param radius: signal radius
    :param noise: standard deviation of the signal
    :param e: eccentricity of the signal
    :return: tpcf with the specified bao signal
    """
    radius_new = int(radius / space)
    if radius_new > tpcf.shape[0]:
        raise ValueError('Radius too large')
    b = np.sqrt(1 / (2 - e ** 2))
    a = np.sqrt(1 - b ** 2)
    print(a, b)
    center = int(tpcf.shape[0] / 2)
    for radius in range(radius_new - 50, radius_new + 50):
        magnitude = gaussian(radius, radius_new, noise)
        for x in range(int(center - radius * a), int(center + radius * a)):
            y = int(max(0, np.sqrt(radius ** 2 - ((x - center) ** 2 / (a ** 2)))) * b)
            tpcf[x, center + y] = magnitude
            tpcf[x, center - y] = magnitude
        for y in range(int(center - radius * b), int(center + radius * b)):
            x = int(max(0, np.sqrt(radius ** 2 - (y - center) ** 2 / (b ** 2))) * a)
            tpcf[center + x, y] = magnitude
            tpcf[center - x, y] = magnitude
    return tpcf


def gaussian_noise(tpcf, sigma):
    """
    Add Gaussian noise to a tpcf(sigma, pi) matrix
    :param tpcf:
    :param sigma: standard deviation of the Gaussian noise
    :return: tpcf with the Gaussian noise
    """
    for (x, y), val in np.ndenumerate(tpcf):
        tpcf[x, y] = np.random.normal(val, sigma)
    return tpcf


def plot_contourf(tpcf, vmin, vmax, levels, extent):
    """
    Self-defined contour-plotting
    :param tpcf:
    :param vmin: min value shown
    :param vmax: max value shown
    :param levels: number of levels
    :param extent: extent
    :return:
    """
    rc('font', family='serif')
    rc('font', size=16)
    plt.tight_layout()
    plt.ylabel(r"$\pi $", fontsize=16)
    plt.xlabel(r"$\sigma $", fontsize=16)
    levels = np.arange(vmin, vmax, (vmax - vmin) / levels)
    plt.contourf(tpcf, levels=levels, origin='lower', extent=extent)
    plt.colorbar()
    plt.tight_layout()
    plt.show()

def plot_threshold(tpcf, extent, title):
    tpcf_sorted = tpcf.flatten()
    tpcf_sorted = sorted(tpcf_sorted)
    first = int(len(tpcf_sorted)/20)
    last = first * 19
    rc('font', family='serif')
    rc('font', size=16)
    plt.tight_layout()
    plt.ylabel(r"$\pi $", fontsize=16)
    plt.xlabel(r"$\sigma $", fontsize=16)
    plt.imshow(tpcf, extent=extent, vmax=tpcf_sorted[last], vmin=tpcf_sorted[first])
    plt.title(title)
    plt.colorbar()
    plt.show()



def interpolate(tpcf, axis=0, threshold=10):
    """
    Using interpolation to fill in the empty entries in the tpcf matrix
    :param tpcf: ndarray
    :param axis: the axis to interpolate along
    :param threshold: for each array, do interpolation only if the number of non-zero entries is above threshold
    :return: filled tpcf matrix
    """
    tpcf2 = np.copy(tpcf)
    if axis == 0:
        for i in range(1, tpcf.shape[0]):
            temp = tpcf[i, :]
            zeros = np.nonzero(temp == 0)[0]
            nonzeros = np.nonzero(temp != 0)[0]
            nonzero_vals = temp[nonzeros]
            if len(nonzeros) >= threshold:
                temp[temp == 0] = np.interp(zeros, nonzeros, nonzero_vals)
                tpcf2[i, :] = temp
    elif axis == 1:
        for i in range(1, tpcf.shape[1]):
            temp = tpcf[:, i]
            zeros = np.nonzero(temp == 0)[0]
            nonzeros = np.nonzero(temp != 0)[0]
            nonzero_vals = temp[nonzeros]
            if len(nonzeros) >= threshold:
                temp[temp == 0] = np.interp(zeros, nonzeros, nonzero_vals)
                tpcf2[:, i] = temp
    return tpcf2
