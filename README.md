# Legendre Expansion of Two-Point Correrlation Function


<!-- toc -->

- [Prerequisites](#prerequisites)
- [Installing](#installing)
- [Introduction](#introduction)
- [Usage](#usage)
  * [Legendre expansion](#legendre-expansion)
  * [Plotting](#plotting)
  * [Test with mocks](#test-with-mocks)
  

<!-- tocstop -->

## Prerequisites

This implementation depends on the following Python packages:
* [AstroPy](http://www.astropy.org)
* [SciPy library](https://github.com/scipy/scipy)
* [NumPy](http://www.numpy.org)
* [Matplotlib](http://matplotlib.org) (optional)

The included configurations rely on DR9 data from the Sloan Digital Sky Survey:
* [SDSS archive](https://data.sdss.org/sas/dr9/boss/lss/)

## Installing

Clone the repository to your local system:
```
git clone https://github.com/yujie-liu/lasspia.git
```

If desired, use virtualenv to install the prerequisites:
```
virtualenv venv
source venv/bin/activate
pip install numpy astropy scipy
```

Download the SDSS galaxy survey catalogs:
```
cd lasspia/
./downloadDR9.py
```

## Introduction

The Legendre expansion mainly consists of the following steps: 
1. Read from `*_integration.fits` file and generate the 2D TPCF matrix: tpcf(sigma, pi)
2. Convert tpcf(sigma, pi) to tpcf(s, mu)
3. Expand tpcf(s, mu) on Legendre polynomials

It is also possible to generate the tpcf(s, mu) matrix directly without the intermediate 
tpcf(sigma, pi) matrix. However, for debugging purposes tpcf(sigma, pi) is still calculated
and can be plotted.

## Usage

### Legendre expansion

* First run Lasspia to get the *_integration.fits file. See [instructions from Lasspia](https://github.com/betchart/lasspia#a-complete-example).
* Make sure the correct configuration file is in the configs folder.
* Run the expansion routine with specified config and *_integration.fits file:
```
./legendre.py DR12_finebins.py DR12_finebins_integration.py
```
* Instead using `./` to execute, depending on the python environment,
 `python` or `python3` will also do.
* After the running is complete, a new fits file 
named *_legendre.fits will be generated in the current directory.

Flags: 
* To see the TPCF matrix as a heatmap, run with the `-p` flag:
```
./legendre.py DR12_finebins.py DR12_finebins_integration.py -p
```
* There has been confusion in the scale factors in DD, DR and RR matrices; 
for instance, the `TwoD_tpcf_*.txt` formatted mocks need to double the DD term. 
The `-dd`, `-dr`, `-rr` flags are created to help with this: 
using those flags will scale the corresponding term by a factor of 2 
and they can be used simultaneously. An example:
```
./legendre.py DR12_finebins.py DR12_finebins_integration.py -dd -dr
```

Notice:
* It is possible to see the warning `RuntimeWarning: invalid value encountered in true_divide`,
because the calculation of tpcf is `(DD-2DR+RR)/RR`, where if any value in RR is zero, 
this warning will be raised. 
* For a file `*_integration.fits`, if the file `*_legendre.fits` already exists, 
running `legendre.py` on it again will overwrite `*_legendre.fits`.
* All the input fits files are expected to be found in the `../Data` folder.

### Plotting
With the newly generated *_legendre.fits file, we can plot the expansion results:
```
./legendre_test DR12_finebins_legendre.py
```

### Test with mocks
For mocks formatted the same as `TwoD_tpcf_*.txt`, run `TwoD_tpcf.py`to generate Legendre expansion and plot the results.
For example:
```
./TwoD_tpcf.py TwoD_tpcf_v11.txt
```
This will calculate the sigma scaler that minimizies the l=2 term's distance to 0.
It also generates a `tpcf_legendre.fits` file. Feed this file to `legendre_test.py`,
we can get the plot of the expansion results of the mock:
```
./legendre_test tpcf_legendre.py
```


Usage syntax and a full list of options can be obtained via the `-h` or `--help` flag:
```
./lasspia.py --help
```
