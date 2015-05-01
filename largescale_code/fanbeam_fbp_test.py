# approximate fbp alg.

from pylab import *
import numpy

ion()

convolve = numpy.convolve

import tomo2D as tomo
import time
from phantoms_tomo2D import generate_shepp_logan

phshepp = generate_shepp_logan()

print '''This program is designed to run in an ipython shell with matplotlib.
From terminal command line run:
ipython --pylab

From within the ipython shell execute:
run fanbeam_fbp_test.py


'''

nx0 = 512
ny0 = 512
xlen0 = 2.
ylen0 = 2.
x00 = -1.
y00 = -1.
fbp_image = tomo.image2D(shape=(nx0,ny0),x0=x00,xlen=xlen0,y0=y00,ylen=ylen0)
image = tomo.image2D(shape=(nx0,ny0),x0=x00,xlen=xlen0,y0=y00,ylen=ylen0)

print "embedding phantom in image"
phshepp.collapse_to(image)

ns0 = 512
nu0 = 512
slen0 = 2.*pi
ulen0 = 4.
s00 = 0.
u00 = -2.


source_rad = 5.
source_det = 8.

du = ulen0/nu0
dup = du*source_rad/source_det


# from kak-slaney
def ramp_kernel(np,du):

   it_is_even = 0
   if mod(np,2)==0:
      it_is_even = 1

   if it_is_even:

      filter = arange(-np/2., np/2., 1.)
      filter = ((-1.)**filter)/(2.*du*du*pi*filter+ du*du*pi)- 1./((8.*du*du)*(pi*filter/2. +pi/4.)**2)

   else:

      filter = arange(-(np-1)/2., (np-1)/2. + 1., 1.)
      filter[(np-1)/2+1::2] = -1./(pi*du*filter[(np-1)/2+1::2])**2
      filter[(np-1)/2-1::-2] = -1./(pi*du*filter[(np-1)/2-1::-2])**2
      filter[(np-1)/2::2] = 0.
      filter[(np-1)/2::-2] = 0.
      filter[(np-1)/2] = 1./(4.*du*du)

   return filter/2.


sino = tomo.sinogram2D(\
   config_name='circular_fan',\
   parms={"source_to_detector":source_det,\
          "radius"            :source_rad},\
   shape=(ns0,nu0),\
   s0=s00,slen=slen0, u0 = u00, ulen = ulen0)

work_sino = tomo.sinogram2D(\
   config_name='circular_fan',\
   parms={"source_to_detector":source_det,\
          "radius"            :source_rad},\
   shape=(ns0,nu0),\
   s0=s00,slen=slen0, u0 = u00, ulen = ulen0)

phshepp.project_to(sino)

rk = ramp_kernel(2*nu0-1,dup) *dup

data_weight = arange( u00 + du/2., u00+ du/2. +ulen0, du)
data_weight *= source_rad/source_det
data_weight = source_rad/sqrt(source_rad**2 + data_weight**2)

print "filtering data..."
t0 = time.time()
for i in range(ns0):
   work_sino.mat[i,:] = convolve(data_weight*sino.mat[i,:],rk,0)
print "filtering time: ",time.time()-t0

print"back projecting ..."
t0 = time.time()
work_sino.weighted_back_project_to(fbp_image,fov=1.)
print "backprojection time: ",time.time()-t0

print "fbp image available in fbp_image.mat, and reference discretized phantom is in image.mat"

print '''
example imshow command:
imshow(fbp_image.mat.transpose()[::-1],vmin = 1.0,vmax=1.05,cmap=cm.gray,interpolation='nearest')'''


