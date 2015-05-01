sino_proj2D_CUDA:
	nvcc --compiler-options '-O3 -lpthread -fpic' -o sino_proj2D_CUDA.so -shared sino_proj2D_CUDA_wholeSino.cu
sino_proj2D:
	f2py -c sino_proj2D.f -m sino_proj2D

all: sino_proj2D sino_proj2D_CUDA
