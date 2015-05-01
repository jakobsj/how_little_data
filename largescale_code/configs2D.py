
import numpy

sin=numpy.sin
sqrt=numpy.sqrt
cos=numpy.cos
tan=numpy.tan
pi=numpy.pi

def config_circular_fan(s,parms):

   radius=parms["radius"]
   source_to_detector=parms["source_to_detector"]

# Location of the source
   x_source=radius * cos(s)
   y_source=radius *sin(s)

   
# detector center
   x_det_center=(radius - source_to_detector )*cos(s)
   y_det_center=(radius - source_to_detector )*sin(s)

# unit vector in the direction of the detector line
   eux = -sin(s)
   euy =  cos(s)

# Unit vector in the direction perpendicular to the detector line
   ewx = cos(s)
   ewy = sin(s)

   return\
      x_source,y_source,\
      x_det_center,y_det_center,\
      eux,euy,\
      ewx,ewy

def config_tomosynthesis(s,parms):

   radius=parms["radius"]
   rotation_center_to_detector=parms["rotation_center_to_detector"]

# Location of the source
   x_source=radius * sin(s)
   y_source=radius *cos(s)

   
# detector center
   x_det_center= 0.
   y_det_center= - rotation_center_to_detector

# unit vector in the direction of the detector line
   eux = 1.
   euy = 0.

# Unit vector in the direction perpendicular to the detector line
   ewx =  0.
   ewy =  1.

   return\
      x_source,y_source,\
      x_det_center,y_det_center,\
      eux,euy,\
      ewx,ewy


configs={\
"circular_fan"                : config_circular_fan,\
"tomosynthesis"               : config_tomosynthesis\
}
