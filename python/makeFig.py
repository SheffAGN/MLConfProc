import numpy as np
from astropy.io import fits
from astropy.io import ascii
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from astropy.visualization import ImageNormalize, ZScaleInterval

#Image data:
im = fits.open('../data/diffexp-01-201708030001.fits')
imdat = im[1].data
w = WCS(im[1].header)

#Detected candidates:
c = fits.open('../data/diaSrc-01-201708030001.fits')
cdat = c[1].data
cra = 180. * cdat['coord_ra'] / np.pi
cdec = 180. * cdat['coord_dec'] / np.pi
cc = SkyCoord(ra=cra * u.degree, dec=cdec * u.degree)

#Real sources:
var = ascii.read("../data/GOTO_01_20170803_0001.var")
cv = SkyCoord(ra=var['RA'] * u.degree, dec=var['Dec'] * u.degree)

tra = ascii.read("../data/GOTO_01_20170803_0001.trans")
ct = SkyCoord(ra=tra['RA'] * u.degree, dec=tra['Dec'] * u.degree)

#Exclude any candidates near real sources:
idx, d2d, d3d = cc.match_to_catalog_sky(cv)
o = np.where((d2d.to(u.arcsec)).value >= 10.)
cc = cc[o]

idx, d2d, d3d = cc.match_to_catalog_sky(ct)
o = np.where((d2d.to(u.arcsec)).value >= 10.)
cc = cc[o]

#Convert candidate positions to pixels, and exclude any near edges:
cxpx, cypx = cc.to_pixel(w)
xlim = np.logical_and(cxpx >= 200, cxpx <= 7700)
ylim = np.logical_and(cypx >= 200, cypx <=5800)
o = np.where(np.logical_and(xlim, ylim))
cxpx, cypx = cxpx[o].astype(np.int_), cypx[o].astype(np.int_)

#Extract regions around those pixels:
inds = np.array([110, 79, 397, 768, 554, 1010])

f = plt.figure(figsize=(8, 3), dpi=320)
gs1 = gridspec.GridSpec(2,3)
gs1.update(left=0.05, right=0.48, wspace=0.05, hspace=0.05)
i = 0
for ind in inds:
    subim = imdat[cypx[ind]-25:cypx[ind]+25,cxpx[ind]-25:cxpx[ind]+25]
    norm = ImageNormalize(subim, interval=ZScaleInterval())
    ax = plt.subplot(gs1[i / 3, i % 3])
    ax.tick_params(axis='both', \
                   bottom='off', labelbottom='off', \
                   left='off', labelleft='off')
    ax.imshow(subim, norm=norm, cmap='binary')
    if i / 3 == 0 and i % 3 == 1:
        ax.set_title('(a) Examples of bogus thumbnails') 
    i += 1

#Do the same for real sources:
cvxpx, cvypx = cv.to_pixel(w)
ctxpx, ctypx = ct.to_pixel(w)

cxpx = np.append(ctxpx, cvxpx)
cypx = np.append(ctypx, cvypx)

xlim = np.logical_and(cxpx >= 200, cxpx <= 7700)
ylim = np.logical_and(cypx >= 200, cypx <=5800)
o = np.where(np.logical_and(xlim, ylim))
cxpx, cypx = cxpx[o].astype(np.int_), cypx[o].astype(np.int_)

inds = np.array([1, 2, 3, 75, 120, 121])

gs1 = gridspec.GridSpec(2,3)
gs1.update(left=0.52, right=0.95, wspace=0.05, hspace=0.05)
i = 0
for ind in inds:
    subim = imdat[cypx[ind]-25:cypx[ind]+25,cxpx[ind]-25:cxpx[ind]+25]
    norm = ImageNormalize(subim, interval=ZScaleInterval())
    ax = plt.subplot(gs1[i / 3, i % 3])
    ax.tick_params(axis='both', \
                   bottom='off', labelbottom='off', \
                   left='off', labelleft='off')
    ax.imshow(subim, norm=norm, cmap='binary_r')
    if i / 3 == 0 and i % 3 == 1:
        ax.set_title('(b) Examples of real thumbnails') 
    i += 1

plt.savefig('../figs/plot.png')
