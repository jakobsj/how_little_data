# complex phantom. ideal recovery

import tomo2D as tomo
import numpy as np



using_CUDA = 0


#Full image size is 1024x1024, problem can be reduced by scalefactor
scalefactor = 8   # choose to be a power of 2


if using_CUDA:
   npbl=512/scalefactor    #number of blocks for projection
   npth = 4                #threads per block for projection
   nbpbl = 128/scalefactor #number of blocks for backprojection
   nbpth = 4               #threads per block for backprojection


# create this directory before running
basedir = 'data/study2/'

# MAIN PARAMETER BLOCK
imagefile = "l1min"+str(scalefactor)+"_lam1e-4_"
phantomfile = "structureWalnut.npy"
lam = 1.e-4
itermax = 100001 
imagedumpset=[100,1000,10000]
nviews = [10,18,19,20,21,22,24,26,28,30,40,60]
# end of MAIN PARAMETER BLOCK





if using_CUDA:
   tomo.set_GPU_dev(0)
   print "GPU = ",tomo.get_GPU_dev()




# Image parameters
#    numbers of pixels in each dimension of the image array
nx0 = 1024/scalefactor
ny0 =  1024/scalefactor
#    physical lengths of image array
xlen0 = 2.0             
ylen0 = 2.0             
#    physical corner of image array
x00 = -1.
y00 = -1.
# for the present theoretical studies physical dimension has little meaning
# because we do not attempt to match actual X-ray transmission intensity.
dx = xlen0/nx0
dy = ylen0/ny0

#generate an indicator mask for the largest inscribed circle of the image
mask_rad = 1.0*xlen0*0.5
xar = np.arange(x00+dx/2.,x00+xlen0,dx)[:,np.newaxis]*np.ones([ny0])
yar = np.ones([nx0])[:,np.newaxis]*np.arange(y00+dy/2.,y00+ylen0,dy)
mask = np.ones([nx0,ny0],dtype='float64')
mask[xar**2 + yar**2>mask_rad**2] = 0.




image2Dsl = tomo.image2D(\
        shape=(nx0,ny0),\
        xlen= xlen0, ylen=ylen0,\
        x0 = x00, y0= y00)
q_image=image2Dsl.duplicate()
temp_image=image2Dsl.duplicate()

xim= tomo.image2D(shape=(nx0,ny0),x0=x00,xlen=xlen0,y0=y00,ylen=ylen0)
wimq= tomo.image2D(shape=(nx0,ny0),x0=x00,xlen=xlen0,y0=y00,ylen=ylen0)
yid = tomo.image2D(shape=(nx0,ny0),x0=x00,xlen=xlen0,y0=y00,ylen=ylen0)
wimp= tomo.image2D(shape=(nx0,ny0),x0=x00,xlen=xlen0,y0=y00,ylen=ylen0)
xbarim= tomo.image2D(shape=(nx0,ny0),x0=x00,xlen=xlen0,y0=y00,ylen=ylen0)



image2Dsl.mat = np.load(phantomfile)[::scalefactor,::scalefactor]
image2Dsl.mat *= mask


l1c0 = abs(image2Dsl.mat).sum(dtype='float64')
nsp = (image2Dsl.mat>0.000000001).sum(dtype='float64')
print "sparsity: ",nsp," out of ", mask.sum(dtype='float64')
print "base L1: ", l1c0




# sinogram parameters
nu0  = 2048/scalefactor           # number of detector bins on a linear detector
slen0 = 2.*np.pi                  # scanning arc length in radians
s00 = 0.                          # starting angle
srad = 4.                         # X-ray source to iso-center distance
sd = 8.                           # X-ray source to detector center distance

#compute detector length so that it captures the largest inscribed circle
alp = np.arcsin((xlen0/2.)/srad)
ulen0 = 2.*np.tan(alp)*sd
u00 = -ulen0/2.




#main loop over projection view number
for ns0 in nviews:
   print "views: ", ns0



   sino = tomo.sinogram2D(\
      config_name='circular_fan',\
      parms={"source_to_detector":sd,\
             "radius"            :srad},\
      shape=(ns0,nu0),\
      s0=s00,slen=slen0, u0 = u00, ulen = ulen0)


   work_sino = sino.duplicate()
   ysino = sino.duplicate()



   if using_CUDA:
      image2Dsl.cproject_to(sino,npbl,npth)
   else:
      image2Dsl.project_to(sino)

   


# power method computation
   poweriter=20
   temp_image.mat = (np.random.randn(nx0,ny0).astype('float64'))*mask

   #projection spectral norm
   for iter in range(1,poweriter):
      temp_image.mat  *= mask
      if using_CUDA:
         temp_image.cproject_to(work_sino,npbl,npth)
      else:
         temp_image.project_to(work_sino)
      work_sino.mat /= work_sino.mag()

      q_image.mat.fill(0.)
      if using_CUDA:
         work_sino.cback_project_rays_to(q_image,nbpbl,nbpth)
      else:
         work_sino.back_project_rays_to(q_image)
      q_image.mat  *= mask
      lprojection = q_image.mag()
      q_image.mat /=lprojection
      temp_image.mat = 1.*q_image.mat


   
   lid = 1.0    #for L1 minimization the gradient operator is replaced by the identity

   nu = lprojection/lid
#   print "nu: ",nu


   poweriter=20
   temp_image.mat = (np.random.randn(nx0,ny0).astype('float64'))*mask
   for iter in range(1,poweriter):
      temp_image.mat  *= mask
      if using_CUDA:
         temp_image.cproject_to(work_sino,npbl,npth)
      else:
         temp_image.project_to(work_sino)
      yid.mat = nu*temp_image.mat
      ltotal = np.sqrt(yid.mag()**2 + work_sino.mag()**2)
      yid.mat /=ltotal
      work_sino.mat /= ltotal


      qtemp=nu*yid.mat
      q_image.mat.fill(0.)
      if using_CUDA:
         work_sino.cback_project_rays_to(q_image,nbpbl,nbpth)
      else:
         work_sino.back_project_rays_to(q_image)
      q_image.mat += qtemp
      q_image.mat  *= mask
      ltotal = q_image.mag()
      q_image.mat /=ltotal
      temp_image.mat = 1.*q_image.mat

# end of power method computations




   tau = 1./ltotal
   sig = 1./ltotal
   theta = 1.0
   xbarim.mat.fill(0.)
   xim.mat.fill(0.)
   ysino.mat.fill(0.)
   yid.mat.fill(0.)
  

   npix = mask.sum(dtype='float64')


   l1s = []
   derrs = []
   ierrs = []


# main CP primal-dual loop
   for iter in range(itermax):

      xbarim.mat *= mask
      if using_CUDA:
         xbarim.cproject_to(work_sino,npbl,npth)
      else:
         xbarim.project_to(work_sino)
      ddist=sino.dist_to(work_sino)
      ddistn= ddist/np.sqrt(1.*ns0*nu0)
      derrs.append(ddistn)       #store data error


      ysino.mat = (ysino.mat + sig*(work_sino.mat - sino.mat))



      currentl1 = abs(xbarim.mat).sum(dtype='float64')
      l1s.append(currentl1/l1c0)    #store ratio of current TV to phantom TV

      tg = nu*xbarim.mat
      idist = xbarim.dist_to(image2Dsl)
      ierr=idist/np.sqrt(1.*npix)
      ierrs.append(ierr)            #store image RMSE


      ptil= yid.mat + sig*tg
      ptilmag = abs(ptil)
      ptilmag[ptilmag<(lam/nu)]=(lam/nu)
      yid.mat = (lam/nu)*ptil/ptilmag



      wimq.mat.fill(0.)
      if using_CUDA:
         ysino.cback_project_rays_to(wimq,nbpbl,nbpth)
      else:
         ysino.back_project_rays_to(wimq)
      wimq.mat *= mask

      wimp.mat = nu*yid.mat
      wimp.mat *= mask

      xold = xim.mat*1.
      xim.mat = xim.mat -tau*(wimq.mat+wimp.mat)
      xim.mat *= mask

      xbarim.mat = xim.mat + theta*(xim.mat  -xold)
      xbarim.mat *= mask


      if iter in imagedumpset:
         print "dumping image at iteration: ",iter
         np.save(basedir+imagefile+"nv"+str(ns0)+"_iter"+str(iter)+".npy",xim.mat)



   np.save(basedir+imagefile+"nv"+str(ns0)+"_derrs.npy",np.array(derrs))
   np.save(basedir+imagefile+"nv"+str(ns0)+"_ierrs.npy",np.array(ierrs))
   np.save(basedir+imagefile+"nv"+str(ns0)+"_l1s.npy",np.array(l1s))
   np.save(basedir+imagefile+"nv"+str(ns0)+"_image.npy",xim.mat)
         
