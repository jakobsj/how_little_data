

import tomo2D as tomo
from numpy import *
import time
from phantoms_tomo2D import generate_shepp_logan_HC

phshepp = generate_shepp_logan_HC()

#
print '''This program is designed to run in an ipython shell with matplotlib.
From terminal command line run:
ipython --pylab

From within the ipython shell execute:
run testCUDA.py


'''


print tomo.get_GPU_dev()
tomo.set_GPU_dev(0)
print 'GPU dev now: ',tomo.get_GPU_dev()


#image parameters
nx0=1024
ny0=1024
xlen0=2.0
ylen0=2.0
x00=-1.0
y00=-1.0


# scan configuration parameters
srad = 4.
sd = 8.

#sinogram parameters
ns0=200
slen0=2.*pi
s00=0.
nu0= 2048
alp = arcsin((xlen0/2.)/srad)
ulen0 = 2.*tan(alp)*sd
u00 = -ulen0/2.


nblocks_projection = 512
nthreads_projection = 4

nblocks_backprojection = 128
nthreads_backprojection = 4

fsino = tomo.sinogram2D(\
   config_name='circular_fan',\
   parms={"source_to_detector":sd,\
          "radius"            :srad},\
   shape=(ns0,nu0),\
   s0=s00,slen=slen0, u0 = u00, ulen = ulen0)


cuda_sino =fsino.duplicate()




image = tomo.image2D(shape=(nx0,ny0),x0=x00,xlen=xlen0,y0=y00,ylen=ylen0)
fimage = image.duplicate()
cuda_image = image.duplicate()

print "Embedding phantom to an image..."
phshepp.collapse_to(image)

print('starting fortran projection')
t0 = time.time()
image.project_to(fsino)
print 'finished fortran projection in: ', time.time() - t0
print 'result in fsino'


print('starting CUDA projection')
t0 = time.time()
image.cproject_to(cuda_sino,nblocks = nblocks_projection, blocksize = nthreads_projection)
print 'finished CUDA projection in: ', time.time() - t0
print 'result in cuda_sino'


t0 = time.time()
print('starting fortran backprojection')
fsino.back_project_rays_to(fimage)
print 'finished fortran bp in: ', time.time() - t0
print 'result in fimage'

t0 = time.time()
print('starting CUDA backprojection')
fsino.cback_project_rays_to(cuda_image,nblocks = nblocks_backprojection, blocksize = nthreads_backprojection)
print 'finished CUDA bp in: ', time.time() - t0
print 'result in cuda_image'

print '''
example imshow command:
imshow(cuda_image.mat.transpose()[::-1],cmap=cm.gray,interpolation='nearest')'''



print '''
example imshow command for difference image:
imshow((fimage.mat-cuda_image.mat).transpose()[::-1],vmin = -0.01,vmax=0.01,cmap=cm.gray,interpolation='nearest')'''


