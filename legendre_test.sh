rm *_legendre.fits
python3 legendre.py '/home/yujie/Documents/lasspia/configs/cmassS_coarse.py' 'cmassS_coarse_integration.fits' >> log.txt

venv/bin/python legendreTest.py >> log.txt
