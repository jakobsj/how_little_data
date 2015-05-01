import numpy
from tomo2D import phantom2D
from tomo2D import ellipse
import random

pi=numpy.pi
sin=numpy.sin
cos=numpy.cos
sqrt=numpy.sqrt




def bar_object():
   ph=phantom2D()
   e=ellipse(0.,0.,1.0,5.0,0.5,0.)
   ph.add_component(e)

   e=ellipse(0.1,0.,.5,.05,0.1,0.4)
   ph.add_component(e)

   e=ellipse(-0.5,-0.1,-0.2,.1,0.2,0.)
   ph.add_component(e)

   e=ellipse(0.5,0.2,.2,0.14,0.2,-1.)
   ph.add_component(e)

   return ph



def generate_random_spots(my_seed=0.3,\
                          num_ellipse=30,\
                          x0_range=[0.2,0.6],\
                          y0_range=[-0.3,0.3],\
                          half_axis_range=[0.01,0.04],\
	                  atten_range=[0.9 ,1.1],\
		          angle_range=[0.,pi/2.]):


   ph=phantom2D()
   e=ellipse(0.1,0.,1.0,0.8,0.6,0.)
   ph.add_component(e)

   e=ellipse(0.4,0.,-1.0,0.3,0.4,0.)
   ph.add_component(e)

   e=ellipse(-0.4,0.1,0.05,0.25,0.15,0.)
   ph.add_component(e)
   random.seed(my_seed)

   for i in range(num_ellipse):
      x0=random.uniform(x0_range[0],x0_range[1])
      y0=random.uniform(y0_range[0],y0_range[1])
      atten=random.uniform(atten_range[0],atten_range[1])
      ax=random.uniform(half_axis_range[0],half_axis_range[1])
      ay=random.uniform(half_axis_range[0],half_axis_range[1])
      angle=random.uniform(angle_range[0],angle_range[1])

      e=ellipse(x0,y0,atten,ax,ay,angle)
      ph.add_component(e)

   return ph


def generate_random_ellipses(my_seed=0.77,\
                             num_ellipse=10,\
                             center_range=[0.3,0.6],\
                             half_axis_range=[0.02,0.3],\
                             atten_range=[0.01,0.1],\
		             angle_range=[0.,pi/2.]):


   ph=phantom2D() 
   e=ellipse(0,0,1.0,0.9,0.9,0.)
   ph.add_component(e)
   random.seed(my_seed)

   for i in range(num_ellipse):
      r0=random.uniform(center_range[0],center_range[1])
      a0=random.uniform(0.,2.*pi)
      x0=r0*cos(a0)
      y0=r0*sin(a0)
      atten=random.uniform(atten_range[0],atten_range[1])
      ax=random.uniform(half_axis_range[0],half_axis_range[1])
      ay=random.uniform(half_axis_range[0],half_axis_range[1])
      angle=random.uniform(angle_range[0],angle_range[1])

      e=ellipse(x0,y0,atten,ax,ay,angle)
      ph.add_component(e)

   return ph


def generate_PET_discs():

   ph=phantom2D()      
   e=ellipse(-75.,0.,1.0,25.,25.,0.)
   ph.add_component(e)

   e=ellipse(-75.,0.,10.,5.,5.,0.)
   ph.add_component(e)
      
      
   e=ellipse(0.,0.,1.0,25.,25.,0.)
   ph.add_component(e)
      
   e=ellipse(0.,0.,10.,5.,5.,0.)
   ph.add_component(e)
      
      
   e=ellipse(75.,0.,1.0,25.,25.,0.)
   ph.add_component(e)

   e=ellipse(75.,0.,10.,5.,5.,0.)
   ph.add_component(e)

   return ph
      
      
def generate_test_ellipse():

   r=sqrt(0.5)      
   ph=phantom2D()
   e=ellipse(-0.1,0.0,1.0,r,r,0.)
   ph.add_component(e)
  
   return ph



def generate_shepp_logan():
   ph=phantom2D()
   e=ellipse(0.,0.,2.0,0.92,0.69,pi/2.)
   ph.add_component(e)
      
   e=ellipse(0,-0.0184,-0.98,0.874,0.6624,pi/2.)
   ph.add_component(e)

   e=ellipse(0.22,0.,-0.02,0.31,0.11,pi*72./180.)
   ph.add_component(e)

   e=ellipse(-0.22,0.,-0.02,0.41,0.16,pi*108./180.)
   ph.add_component(e)

   e=ellipse(0.,0.35,0.01,0.25,0.21,pi/2.)
   ph.add_component(e)

   e=ellipse(0.,0.1,0.01,0.046,0.046,0.)
   ph.add_component(e)

   e=ellipse(0.,-0.1,0.01,0.046,0.046,0.)
   ph.add_component(e)

   e=ellipse(-0.08,-0.605,0.01,0.046,0.023,0.)
   ph.add_component(e)

   e=ellipse(0.0,-0.605,0.01,0.023,0.023,0.)
   ph.add_component(e)

   e=ellipse(0.06,-0.605,0.01,0.046,0.023,pi/2.)
   ph.add_component(e)

   return ph




def generate_shepp_logan_HC():
   ph=phantom2D()
   e=ellipse(0.,0.,2.0,0.92,0.69,pi/2.)
   ph.add_component(e)
   
   e=ellipse(0,-0.0184,-0.98,0.874,0.6624,pi/2.)
   ph.add_component(e)

   e=ellipse(0.22,0.,-0.08,0.31,0.11,pi*72./180.)
   ph.add_component(e)

   e=ellipse(-0.22,0.,-0.08,0.41,0.16,pi*108./180.)
   ph.add_component(e)

   e=ellipse(0.,0.35,0.04,0.25,0.21,pi/2.)
   ph.add_component(e)

   e=ellipse(0.,0.1,0.04,0.046,0.046,0.)
   ph.add_component(e)

   e=ellipse(0.,-0.1,0.04,0.046,0.046,0.)
   ph.add_component(e)

   e=ellipse(-0.08,-0.605,0.04,0.046,0.023,0.)
   ph.add_component(e)

   e=ellipse(0.0,-0.605,0.04,0.023,0.023,0.)
   ph.add_component(e)

   e=ellipse(0.06,-0.605,0.04,0.046,0.023,pi/2.)
   ph.add_component(e)

   return ph

