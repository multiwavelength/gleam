from astropy.io import fits 
import numpy as np
import os, sys
from scipy import signal
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema


filepath = sys.argv[1]

data = fits.open(filepath)[1].data
hdr = fits.open(filepath)[1].header

mask_optical = (data.field('Lambda_air') < 10000)

wl = data.field('Lambda_air')[mask_optical]
line = data.field('OH_SKY')[mask_optical]

data1 = data[mask_optical]
peak_indices = signal.find_peaks_cwt(line, np.arange(10,60))
#peak_indices = argrelextrema(line, np.greater)

plt.plot(wl, line, 'k-')
bright = line > 0.01
#plt.plot(wl[bright], line[bright],  marker='o', c='b', linestyle='None')
plt.plot(wl[peak_indices], line[peak_indices], marker='o', c='b', linestyle='None')
plt.show()
plt.close()


fits.writeto('sky_lines.fits', data1[peak_indices], hdr)

