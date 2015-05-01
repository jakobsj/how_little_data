

import numpy
import configs2D
import sino_proj2D



pi = numpy.pi
zeros = numpy.zeros
array = numpy.array
ones = numpy.ones
sin = numpy.sin
cos = numpy.cos
tan = numpy.tan
sqrt = numpy.sqrt
float32 = numpy.float32
float64 = numpy.float64
indicator_type = numpy.int8
#default_type = numpy.float32
default_type = numpy.float64

default_data_shape = 120,128
default_image_shape = 128,128





using_CUDA = 0

if using_CUDA:
   import ctypes

   csinolib = numpy.ctypeslib.load_library('sino_proj2D_CUDA','.')
   cbackproject = csinolib.backproject
   crayproj = csinolib.rayproj
   setGPUdev = csinolib.setGPUdevice
   getGPUdev = csinolib.getGPUdevice


   def set_GPU_dev(devnum = 0):
      cdevnum = ctypes.c_int()
      cdevnum.value = devnum
      setGPUdev(cdevnum)

   def get_GPU_dev():
      cdevnum = ctypes.c_int()
      cdevnum.value = getGPUdev()
      return cdevnum.value






class sinogram2D(object):
   '''2D sinogram array for tomography'''

   def __init__(self,\
      config_name='circular_fan',\
      parms={"source_to_detector":8.,\
             "radius":4.},\
      shape= default_data_shape,\
      slen=pi, ulen=2.0,\
      s0=0., u0=-1.0,\
      s_include_endpoints=0):
      '''Initializes the sinogram array'''

      self.ns=shape[0]
      self.nu=shape[1]
      

      self.mat=zeros(shape,default_type)
      self.indicator=ones(shape,indicator_type)
      self.frame_vectors=zeros((self.ns,8) , default_type)
      self.config_name=config_name
      self.parms=parms

      self.frame=configs2D.configs[config_name]

      self.slen=slen
      self.ulen=ulen
     
      self.s_include_endpoints=s_include_endpoints
      self.s0=s0
      self.u0=u0
      
      self.ds=self.get_ds()
      self.du=self.get_du()

# The following "for" loop creates an array of the configurations from all
# different source location along the trajectory. This array is used later as
# an input for the reconstruction algorithms done by "sinoproj2D.f"
# It uses file config2D.py.
     
      for ip in range(self.ns):
         s=self.s0+ip*self.ds
         self.frame_vectors[ip]=self.frame(s,self.parms)

   
   def __str__(self):
      '''Gives sinogram matrix parameters'''


      nset=sum(sum(self.indicator))
      frac=(100.0 * nset)/(self.ns * self.nu)
      lulu =  '''
      The sinogram is %f %% full.
      ns = %d  
      nu = %d
      
      ds = %f
      du = %f
     
      s0 = %f
      u0 = %f''' % (frac,self.ns,self.nu,self.ds,self.du,self.s0,self.u0)

      return lulu

   def mag(self):
      return sqrt( ((self.indicator*self.mat)**2).sum(dtype='float64') )

   def dist_to(self,sino):
      '''Calculates L2 distance to sino'''
      
      dist = sqrt( (self.indicator*(self.mat-sino.mat)**2).sum(dtype = "float64") )
      return dist



   def duplicate(self):
      
      new_sino = sinogram2D(\
         config_name=self.config_name,\
         parms=self.parms,\
         shape= (self.ns,self.nu),\
         slen=self.slen, ulen=self.ulen,\
         s0=self.s0, u0=self.u0,\
         s_include_endpoints=self.s_include_endpoints)
      new_sino.mat = self.mat.copy()
      new_sino.indicator = self.indicator.copy()
      new_sino.frame_vectors = self.frame_vectors.copy()

      return new_sino

   
   def get_ds(self):
      if self.s_include_endpoints:
         return self.slen/(self.ns-1.)
      else:
	 return self.slen/self.ns

   def get_du(self):
      return self.ulen/self.nu

   def get_parms(self):
      return self.parms

   
   def frame_vecs(self,s):
      '''Returns the config vectors for the trajectory parm s'''

      clist = self.frame(s,self.parms)
      xloc = array(clist[0:2])
      detc = array(clist[2:4])
      euhat = array(clist[4:6])
      ewhat = array(clist[6:8])
     
      return xloc, detc, euhat, ewhat



 
   def back_project_rays_to(self,image):
      '''Backproject onto image (ray driven)'''

      ns=self.ns
      nu=self.nu
      ds=self.ds
      du=self.du      
      s0=self.s0
      u0=self.u0
      config=self.config_name
     
      dx=image.dx
      dy=image.dy       
      x0=image.x0
      y0=image.y0
      nx=image.nx
      ny=image.ny

      image.mat=array(\
                   sino_proj2D.backproject(\
                      self.mat,\
                      self.indicator,\
                      self.frame_vectors,\
                      ns,ds,nu,du,u0,\
                      dx,dy,x0,y0,nx,ny),\
                   order='C')


 
   if using_CUDA:
      def cback_project_rays_to(self,image,nblocks=64,blocksize=4):
         '''Ray-driven backprojection with CUDA'''

         tempsino = self.mat.astype('float32')
         tempimage = image.mat.astype('float32')
         tempfv = self.frame_vectors.astype('float32')

         ns = ctypes.c_int()
         ns.value = self.mat.shape[0]
         nu = ctypes.c_int()
         nu.value = self.mat.shape[1]
         ds = ctypes.c_float()
         ds.value = float32(self.get_ds())
         du = ctypes.c_float()
         du.value = float32(self.get_du())
         s0 = ctypes.c_float()
         s0.value = float32(self.s0)
         u0 = ctypes.c_float()
         u0.value = float32(self.u0)
      
         nx = ctypes.c_int()
         nx.value=image.mat.shape[0]
         ny = ctypes.c_int()
         ny.value=image.mat.shape[1]
         dx = ctypes.c_float()
         dx.value=float32(image.get_dx())
         dy = ctypes.c_float()
         dy.value=float32(image.get_dy())
         x0 = ctypes.c_float()
         x0.value=float32(image.x0)
         y0 = ctypes.c_float()
         y0.value=float32(image.y0)

         nblk = ctypes.c_uint()
         nblk.value=nblocks
         blksz = ctypes.c_uint()
         blksz.value=blocksize

         cbackproject(\
                 ctypes.c_void_p(tempsino.ctypes.data),\
                 ctypes.c_void_p(self.indicator.ctypes.data),\
                 ctypes.c_void_p(tempfv.ctypes.data),\
                 ns,nu,du,u0,\
                 ctypes.c_void_p(tempimage.ctypes.data),\
                 dx,dy,x0,y0,nx,ny,\
                 nblk,blksz)

         image.mat = tempimage.astype(default_type)



   def back_project_pix_to(self,image,fov=1.0,xc=0.,yc=0.):
      '''Pixel-driven backprojection onto image'''

      ns=self.ns
      nu=self.nu
      ds=self.ds
      du=self.du      
      s0=self.s0
      u0=self.u0
     
      dx=image.dx
      dy=image.dy       
      x0=image.x0
      y0=image.y0
      nx=image.nx
      ny=image.ny

      image.mat=array(\
                   sino_proj2D.pixel_driven_backproj(\
                      self.mat,\
                      self.indicator,\
                      self.frame_vectors,\
                      ds,ns,nu,du,u0,\
                      image.mat,dx,dy,x0,y0,nx,ny,\
                      fov,xc,yc),\
                   order = 'C')


   def weighted_back_project_to(self,image,fov=1.0,xc=0.,yc=0.):
      '''Backproject onto image with fan-beam weighting for linear-detector'''

      ns=self.ns
      nu=self.nu
      ds=self.ds
      du=self.du      
      s0=self.s0
      u0=self.u0
      config=self.config_name
     
      dx=image.dx
      dy=image.dy       
      x0=image.x0
      y0=image.y0
      nx=image.nx
      ny=image.ny

      image.mat=array(\
                   sino_proj2D.weighted_pixel_driven_backproj(\
                      self.mat,\
                      self.indicator,\
                      self.frame_vectors,\
                      ds,ns,nu,du,u0,\
                      image.mat,dx,dy,x0,y0,nx,ny,\
                      fov,xc,yc),\
                   order = 'C')






##############################################################

class image2D(object):
   '''2D image array for tomography'''

   def __init__(self,\
                shape=default_image_shape,\
                xlen=2.0,ylen=2.0,\
		x0=-1.0,y0=-1.0):
      '''Initializes the image array'''

      self.mat=zeros(shape,default_type)
      self.nx=shape[0]
      self.ny=shape[1]      
      self.x0=x0
      self.y0=y0
      self.xlen=xlen
      self.ylen=ylen
      self.dx=self.get_dx()
      self.dy=self.get_dy()
     

   def __str__(self):
      '''Gives image matrix parameters'''
      lulu =  '''
nx = %d  
ny = %d

dx = %f
dy = %f

x0 = %f
y0 = %f

xlen = %f
ylen = %f''' % (self.nx,self.ny,\
              self.dx,self.dy,\
	      self.x0,self.y0,\
              self.xlen,self.ylen)

      return lulu

   def mag(self):
      return sqrt((self.mat**2).sum(dtype='float64'))

   def dist_to(self,image):
      '''Calculates L2 distance to image'''
      
      dist = sqrt( ((self.mat-image.mat)**2).sum(dtype = "float64") )
      return dist



   def duplicate(self):
      
      new_image = image2D(\
         shape= (self.nx,self.ny),\
         xlen=self.xlen, ylen=self.ylen,\
         x0=self.x0, y0=self.y0)
      new_image.mat = self.mat.copy()

      return new_image


   def get_dx(self):
      return self.xlen/self.nx

   def get_dy(self):
      return self.ylen/self.ny



   def add_shape(self,shape):
      '''Puts a 2D shape in the image array
	
The attenuation value for the shape is added to each pixel
whose center is in the shape'''

# int() in the following might not be correct for centers
# off the image array
      dx=self.get_dx()
      dy=self.get_dy()
      nx=self.mat.shape[0]
      ny=self.mat.shape[1]
      ncenterx=int((shape.x0-self.x0)/dx)
      ncentery=int((shape.y0-self.y0)/dy)

#      half_square_len=max(ell.ax,ell.ay)
      nhslenx=int(shape.size/dx)+1
      nhsleny=int(shape.size/dy)+1

      nxi=max(0,ncenterx-nhslenx)
      nyi=max(0,ncentery-nhsleny)

      nxf=min(nx,ncenterx+nhslenx)
      nyf=min(ny,ncentery+nhsleny)


      for i in range(nxi,nxf):
         for j in range(nyi,nyf):
            x=self.x0+(i+0.5)*dx
            y=self.y0+(j+0.5)*dy
            self.mat[i,j]=self.mat[i,j]+shape.pixval(x,y)


    
    
   def project_to(self,sino):
      '''Ray-driven projection onto sinogram'''
     
      ns=sino.ns
      nu=sino.nu     
      ds=sino.ds
      du=sino.du     
      s0=sino.s0
      u0=sino.u0
      config=sino.config_name
      
      dx=self.dx
      dy=self.dy     
      x0=self.x0
      y0=self.y0    
      nx=self.nx
      ny=self.ny

      sino.mat=array(\
                  sino_proj2D.rayproj(\
                     sino.indicator,\
                     sino.frame_vectors,\
                     ns,nu,du,u0,\
                     self.mat,dx,dy,x0,y0,nx,ny),\
                  order='C')




   if using_CUDA:    
      def cproject_to(self,sino,nblocks=64,blocksize=4):
         '''CUDA ray-driven projection onto sinogram'''
     
         tempsino = sino.mat.astype('float32')
         tempimage = self.mat.astype('float32')
         tempfv = sino.frame_vectors.astype('float32')
         ns = ctypes.c_int()
         ns.value = sino.mat.shape[0]
         nu = ctypes.c_int()
         nu.value = sino.mat.shape[1]
         ds = ctypes.c_float()
         ds.value = float32(sino.get_ds())
         du = ctypes.c_float()
         du.value = float32(sino.get_du())
         s0 = ctypes.c_float()
         s0.value = float32(sino.s0)
         u0 = ctypes.c_float()
         u0.value = float32(sino.u0)
      
         nx = ctypes.c_int()
         nx.value=self.mat.shape[0]
         ny = ctypes.c_int()
         ny.value=self.mat.shape[1]
         dx = ctypes.c_float()
         dx.value=float32(self.get_dx())
         dy = ctypes.c_float()
         dy.value=float32(self.get_dy())
         x0 = ctypes.c_float()
         x0.value=float32(self.x0)
         y0 = ctypes.c_float()
         y0.value=float32(self.y0)


         nblk = ctypes.c_uint()
         nblk.value=nblocks
         blksz = ctypes.c_uint()
         blksz.value=blocksize

         crayproj(ctypes.c_void_p(tempsino.ctypes.data),\
                  ctypes.c_void_p(sino.indicator.ctypes.data),\
                  ctypes.c_void_p(tempfv.ctypes.data),\
                  ns,nu,du,u0,\
                  ctypes.c_void_p(tempimage.ctypes.data),\
                  dx,dy,x0,y0,nx,ny,\
                  nblk,blksz)

         sino.mat = tempsino.astype(default_type)

    
     
 
#######################################################################

class phantom2D(object):
   '''A collection of 2D shapes
'''
   
   def __init__(self):
      self.num_components = 0
      self.components=[]

   def __str__(self):
      lulu=""
      for i in range(1,self.num_components+1):

	      lulu+="component %d: "%(i)+ self.components[i-1].__str__()+"\n"

      return lulu

   def clear(self):
      self.__init__()

   def add_component(self,component):
      '''adds either a shape to the phantom'''
      self.num_components+=1
      self.components.append(component)

   def collapse_to(self,image):
      '''collapse the 2D phantom to an image array

Note: this function is additive, so be sure to clear the image array first'''

      for i in range(0,self.num_components):
	 image.add_shape(self.components[i])


   def project_to(self,sino):
      '''Projects phantom onto a 2D sinogram.'''

      for i in range(0,self.num_components):
         self.components[i].project_to(sino)



#######################################################################


class shape2D(object):
   '''geometric shape class for building 2D phantoms'''

   def __init__(self,x0=0.,y0=0.,size=0.,att=1.0):
      '''(x0,y0) is the center of the shape, and
size is the radius that encompasses the shape'''
      self.x0=x0
      self.y0=y0
      self.size=size
      self.att=att

   def __str__(self):
      lulu="""Shape parms:
attenuation = %f
center = (%f,%f)
size = %f \n""" % (self.att,self.x0,self.y0,self.size)
      return lulu

   def pixval(self,x,y):
      pass

   def project_to(self,sino):
      pass



########################################################################

class ellipse(shape2D):

   def __init__(self,x0=0.,y0=0.,att=1.0,ax=1.0,ay=1.0,gamma=0.0):
      self.ax=ax
      self.ay=ay
      self.gamma=gamma
      size=max(ax,ay)
      shape2D.__init__(self,x0,y0,size,att)

   def __str__(self):
      lulu = "An ellipse!\n"+shape2D.__str__(self)+ '''ellipse parms:
ax = %f  
ay = %f
gamma = %f \n''' % (self.ax,self.ay,self.gamma)

      return lulu


   def pixval(self,x,y):
      '''Returns attenuation value if x,y is in the ellipse'''

      rel_x = x - self.x0
      rel_y = y - self.y0

      mu=self.att

      cg=cos(self.gamma)
      sg=sin(self.gamma)
      ax=self.ax
      ay=self.ay

      ellipse_lhs=((rel_x*cg+rel_y*sg)/ax)**2+((rel_y*cg-rel_x*sg)/ay)**2
      if ellipse_lhs<=1.0:
         return mu
      else:
         return 0.0


   def project_to(self,sino):
      '''Project ellipse onto predefined sino object'''

      ns=sino.mat.shape[0]
      nu=sino.mat.shape[1]
      ds=sino.get_ds()
      du=sino.get_du()
      s0=sino.s0
      u0=sino.u0
      config=sino.config_name

      ax=self.ax
      ay=self.ay
      x0=self.x0
      y0=self.y0
      gamma=self.gamma
      att=self.att

      sino.mat += array(\
                     sino_proj2D.ellproj(\
                        sino.indicator,sino.frame_vectors,\
                        ns,nu,du,u0,\
                        ax,ay,x0,y0,gamma,att),\
                     order = 'C')



########################################################
