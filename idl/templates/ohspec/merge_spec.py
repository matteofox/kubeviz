import numpy as np
import pyfits as pf

iz = pf.open('kmos_oh_spec_iz.fits')
iz_spec = iz[1].data
iz_head = iz[1].header
iz.close()

yj = pf.open('kmos_oh_spec_yj.fits')
yj_spec = yj[1].data
yj_head = yj[1].header
yj.close()

h = pf.open('kmos_oh_spec_h.fits')
h_spec = h[1].data
h_head = h[1].header
h.close()

k = pf.open('kmos_oh_spec_k.fits')
k_spec = k[1].data
k_head = k[1].header
k.close()

iz_wave = iz_head['CRVAL1'] + iz_head['CDELT1']* (np.array(range(1,len(iz_spec)+1))-iz_head['CRPIX1'])
yj_wave = yj_head['CRVAL1'] + yj_head['CDELT1']* (np.array(range(1,len(yj_spec)+1))-yj_head['CRPIX1'])
h_wave = h_head['CRVAL1'] + h_head['CDELT1']* (np.array(range(1,len(h_spec)+1))-h_head['CRPIX1'])
k_wave = k_head['CRVAL1'] + k_head['CDELT1']* (np.array(range(1,len(k_spec)+1))-k_head['CRPIX1'])

minlam = iz_wave[0]
maxlam = k_wave[-1]

lenall = int((maxlam-minlam)/iz_head['CDELT1'])+1

specall = np.zeros((lenall))
waveall = iz_head['CRVAL1'] + iz_head['CDELT1']* (np.array(range(1,len(specall)+1))-iz_head['CRPIX1'])

iz_first = np.where(abs(waveall - iz_wave[0])<1E-6)[0]
yj_first = np.where(abs(waveall - yj_wave[0])<1E-6)[0]
h_first =  np.where(abs(waveall - h_wave[0])<1E-6)[0]
k_first =  np.where(abs(waveall - k_wave[0])<1E-6)[0]

specall[iz_first:iz_first+len(iz_spec)] = iz_spec
specall[yj_first:yj_first+len(yj_spec)] = yj_spec
specall[h_first:h_first+len(h_spec)] = h_spec
specall[k_first:k_first+len(k_spec)] = k_spec

newhead = iz_head
newhead['NAXIS1'] = lenall
newhead['INSTRUME'] = 'KMOS'
newhead['COMMENT'] = 'This file is distributed with Kubeviz and is used to identify the sky lines \
to derive the instrumental resolution. It has not been tested for any other use. The flux is in arb. \
units and normalized to 100 in each atmospheric band'

pf.writeto('test.fits',specall,newhead, clobber=True)
