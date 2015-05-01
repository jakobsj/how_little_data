// repelacing sys. mat. storage
// loads in whole sino for POCS and bp
#ifndef max
	#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
//const unsigned int blocksize = 256;
//const unsigned int nblocks = 4;
//unsigned int nsimulrays = blocksize*nblocks;
const int nframecomp=8;

__global__ void multiraybp(int* intparms, float* floatparms,
		   float* sinovals, char* sinoinds, float* smat);


__global__ void multirayproj(int* intparms, float* floatparms,
                   float* sinovals, char* sinoinds, float* smat);



extern "C" void backproject(
   float* sinomat, char* indsino,
   float* frame_vectors,
   int ns, int nu,
   float du,float u0,
   float* smat,
   float dx, float dy,
   float x0,float y0,
   int nx,int ny,
   unsigned int nblocks, unsigned int blocksize);


extern "C" void  rayproj(
   float* sinomat, char* indsino,
   float* frame_vectors,
   int ns, int nu,
   float du,float u0,
   float* smat,
   float dx, float dy,
   float x0,float y0,
   int nx,int ny,
   unsigned int nblocks, unsigned int blocksize);





extern "C" void setGPUdevice(int devnum){
  cudaSetDevice(devnum);
}
extern "C" int getGPUdevice(){
  int devnum;
  cudaGetDevice(&devnum);
  return devnum;
}

/////////////////////////////////////
// error checking routine
/////////////////////////////////////
void checkErrors(char *label)
{
  // we need to synchronise first to catch errors due to
  // asynchroneous operations that would otherwise
  // potentially go unnoticed

  cudaError_t err;

  err = cudaThreadSynchronize();
  if (err != cudaSuccess)
  {
    char *e = (char*) cudaGetErrorString(err);
    fprintf(stderr, "CUDA Error: %s (at %s)", e, label);
  }

  err = cudaGetLastError();
  if (err != cudaSuccess)
  {
    char *e = (char*) cudaGetErrorString(err);
    fprintf(stderr, "CUDA Error: %s (at %s)", e, label);
  }
}





void backproject(
   float* sinomat, char* indsino,
   float* frame_vectors,
   int ns, int nu,
   float du,float u0,
   float* smat,
   float dx, float dy,
   float x0,float y0,
   int nx,int ny,
   unsigned int nblocks, unsigned int blocksize) {


  unsigned int nsimulrays = nblocks*blocksize;

  float xSourceCenter,ySourceCenter,xDetCenter, yDetCenter, eux, euy;

  //  float val;

  int ip;

  int ngroups = (nu+nsimulrays-1)/nsimulrays;



  float* floatparms;
  int* intparms;

  floatparms = (float *)malloc(13*sizeof(float));
  intparms = (int *)malloc(5*sizeof(int));


  size_t imsize = nx*ny*sizeof(float);
  float* d_smat;
  cudaMalloc((void**)&d_smat,imsize);

  size_t sinosize = nu*sizeof(float);
  float* d_sinovals;
  cudaMalloc((void**)&d_sinovals,sinosize);

  size_t indsize = nu*sizeof(char);
  char* d_sinoinds;
  cudaMalloc((void**)&d_sinoinds,indsize);

  size_t ipsize = 5*sizeof(int);
  int* d_intparms;
  cudaMalloc((void**)&d_intparms,ipsize);

  size_t dpsize = 13*sizeof(float);
  float* d_floatparms;
  cudaMalloc((void**)&d_floatparms,dpsize);


  checkErrors("memory allocation");

  intparms[0] = 0;
  intparms[1] = ngroups;
  intparms[2] = nx; intparms[3] = ny; intparms[4] = nu;

  floatparms[0] = 0.;
  floatparms[1] = x0; floatparms[2] = y0; floatparms[3] = dx; floatparms[4] = dy;
  floatparms[5] = du; floatparms[6] = u0;

  cudaMemcpy(d_smat,smat,imsize,cudaMemcpyHostToDevice);
  for (ip=0; ip<ns; ip++){

     xSourceCenter=frame_vectors[ip*nframecomp+ 0];
     ySourceCenter=frame_vectors[ip*nframecomp+ 1];

     xDetCenter=frame_vectors[ip*nframecomp+ 2];
     yDetCenter=frame_vectors[ip*nframecomp+ 3];

     eux=frame_vectors[ip*nframecomp+ 4];
     euy=frame_vectors[ip*nframecomp+ 5];

     floatparms[7] = xSourceCenter; floatparms[8] = ySourceCenter;
     floatparms[9] = xDetCenter; floatparms[10] = yDetCenter;
     floatparms[11] = eux; floatparms[12] = euy;



     cudaMemcpy(d_floatparms,floatparms,dpsize,cudaMemcpyHostToDevice);

       cudaMemcpy(d_intparms,intparms,ipsize,cudaMemcpyHostToDevice);
       cudaMemcpy(d_sinovals,sinomat+ip*nu,sinosize,cudaMemcpyHostToDevice);
       cudaMemcpy(d_sinoinds,indsino+ip*nu,indsize,cudaMemcpyHostToDevice);
       checkErrors("copy data to device");

       multiraybp<<<nblocks,blocksize>>>(d_intparms,d_floatparms,
	            d_sinovals,d_sinoinds,d_smat);
       checkErrors("compute on device");
  }
  cudaMemcpy(smat,d_smat,imsize,cudaMemcpyDeviceToHost);
  checkErrors("copy data from device");
  free(intparms);
  free(floatparms);
  cudaFree(d_smat);
  cudaFree(d_sinovals);
  cudaFree(d_sinoinds);
  cudaFree(d_intparms);
  cudaFree(d_floatparms);
}





void  rayproj(
   float* sinomat, char* indsino,
   float* frame_vectors,
   int ns, int nu,
   float du,float u0,
   float* smat,
   float dx, float dy,
   float x0,float y0,
   int nx,int ny,
   unsigned int nblocks, unsigned int blocksize) {


  float xSourceCenter,ySourceCenter,xDetCenter, yDetCenter, eux, euy;
  unsigned int nsimulrays = nblocks*blocksize;

  //  float val;

  int ip;

  int ngroups = (nu+nsimulrays-1)/nsimulrays;



  float* floatparms;
  int* intparms;

  floatparms = (float *)malloc(13*sizeof(float));
  intparms = (int *)malloc(5*sizeof(int));


  size_t imsize = nx*ny*sizeof(float);
  float* d_smat;
  cudaMalloc((void**)&d_smat,imsize);

  size_t sinosize = nu*sizeof(float);
  float* d_sinovals;
  cudaMalloc((void**)&d_sinovals,sinosize);

  size_t indsize = nu*sizeof(char);
  char* d_sinoinds;
  cudaMalloc((void**)&d_sinoinds,indsize);

  size_t ipsize = 5*sizeof(int);
  int* d_intparms;
  cudaMalloc((void**)&d_intparms,ipsize);

  size_t dpsize = 13*sizeof(float);
  float* d_floatparms;
  cudaMalloc((void**)&d_floatparms,dpsize);



  checkErrors("memory allocation");

  intparms[0] = 0;
  intparms[1] = ngroups;
  intparms[2] = nx; intparms[3] = ny; intparms[4] = nu;

  floatparms[0] = 0.;
  floatparms[1] = x0; floatparms[2] = y0; floatparms[3] = dx; floatparms[4] = dy;
  floatparms[5] = du; floatparms[6] = u0;

  cudaMemcpy(d_smat,smat,imsize,cudaMemcpyHostToDevice);
  for (ip=0; ip<ns; ip++){

     xSourceCenter=frame_vectors[ip*nframecomp+ 0];
     ySourceCenter=frame_vectors[ip*nframecomp+ 1];

     xDetCenter=frame_vectors[ip*nframecomp+ 2];
     yDetCenter=frame_vectors[ip*nframecomp+ 3];

     eux=frame_vectors[ip*nframecomp+ 4];
     euy=frame_vectors[ip*nframecomp+ 5];

     floatparms[7] = xSourceCenter; floatparms[8] = ySourceCenter;
     floatparms[9] = xDetCenter; floatparms[10] = yDetCenter;
     floatparms[11] = eux; floatparms[12] = euy;


  

     cudaMemcpy(d_floatparms,floatparms,dpsize,cudaMemcpyHostToDevice);

       cudaMemcpy(d_intparms,intparms,ipsize,cudaMemcpyHostToDevice);
       cudaMemcpy(d_sinovals,sinomat+ip*nu,sinosize,cudaMemcpyHostToDevice);
       cudaMemcpy(d_sinoinds,indsino+ip*nu,indsize,cudaMemcpyHostToDevice);
       checkErrors("copy data to device");

       multirayproj<<<nblocks,blocksize>>>(d_intparms,d_floatparms,
	            d_sinovals,d_sinoinds,d_smat);
       checkErrors("compute on device");

       cudaMemcpy(sinomat+ip*nu,d_sinovals,sinosize,cudaMemcpyDeviceToHost);
       checkErrors("copy data from device");
      


  }

  free(intparms);
  free(floatparms);
  cudaFree(d_smat);
  cudaFree(d_sinovals);
  cudaFree(d_sinoinds);
  cudaFree(d_intparms);
  cudaFree(d_floatparms);
}





__global__ void multiraybp(int* intparms, float* floatparms,
                   float* sinovals, char* sinoinds, float* smat){




__shared__ int ngroups, nx,ny,nu;
__shared__ float x0,y0,dx,dy,du,u0;
__shared__ float xSourceCenter,ySourceCenter,xDetCenter,yDetCenter,eux,euy;




ngroups=intparms[1];
nx=intparms[2]; ny=intparms[3]; nu=intparms[4];

x0=floatparms[1]; y0=floatparms[2]; dx=floatparms[3]; dy=floatparms[4];
du=floatparms[5]; u0=floatparms[6];

xSourceCenter=floatparms[7]; ySourceCenter=floatparms[8];
xDetCenter=floatparms[9]; yDetCenter=floatparms[10];
eux=floatparms[11]; euy=floatparms[12];


int ng, iyOld, ixOld, ix, iy, jp, jpp;
float u, xbin, ybin, xsource, ysource, xl, yl, xdiff, ydiff, xad, yad, slope, slopeinv, x, y;
float travPixlen, yMid, xMid;
float yIntercept, yIntOld, ydist1, ydist2, frac1, frac2;
float xIntercept, xIntOld, xdist1, xdist2;
float val0;

jpp =  blockDim.x*blockIdx.x + threadIdx.x;
for(ng=0;ng<ngroups;ng++){
jp = ng + ngroups*jpp;
if ((jp<nu) && (sinoinds[jp] == 1)){
   val0 = sinovals[jp];

   u = u0+(jp+0.5f)*du;
   xbin = xDetCenter + eux*u;
   ybin = yDetCenter + euy*u;
   xsource = xSourceCenter;
   ysource = ySourceCenter;


   xl=x0;
   yl=y0;

   xdiff=xbin-xsource;
   ydiff=ybin-ysource;
   xad=fabs(xdiff)*dy;
   yad=fabs(ydiff)*dx;

   if (xad>yad) {
      slope=ydiff/xdiff;
      travPixlen=dx*sqrt(1.0f+slope*slope);
      yIntOld=ysource+slope*(xl-xsource);
      iyOld=floor((yIntOld-y0)/dy);
           for (ix=0; ix<nx; ix++) {
	     x=xl+dx*(ix + 1.0f);
	     yIntercept=ysource+slope*(x-xsource);
	     iy=floor((yIntercept-y0)/dy);
	     if (iy == iyOld) {
	       if ((iy >= 0) && (iy < ny)) {
                 smat[ix*ny + iy] = smat[ix*ny + iy]+
	            val0*travPixlen;
		       }
		    }else{
	       yMid=dy*max(iy,iyOld)+yl;
	       ydist1=fabs(yMid-yIntOld);
	       ydist2=fabs(yIntercept-yMid);
	       frac1=ydist1/(ydist1+ydist2);
	       frac2=1.0f-frac1;
	       if ((iyOld >= 0) && (iyOld < ny)) {
                 smat[ix*ny + iyOld] = smat[ix*ny + iyOld]+
	            val0*frac1*travPixlen;
		       }
	       if ((iy>=0) && (iy<ny)) {
                 smat[ix*ny + iy] = smat[ix*ny + iy]+
	            val0*frac2*travPixlen;
		       }
		    }
	     iyOld=iy;
	     yIntOld=yIntercept;
	    }
		    }else{


	   slopeinv=xdiff/ydiff;
	   travPixlen=dy*sqrt(1.0f+slopeinv*slopeinv);
	   xIntOld=xsource+slopeinv*(yl-ysource);
	   ixOld=floor((xIntOld-x0)/dx);
           for (iy=0; iy<ny; iy++) {
	     y=yl+dy*(iy + 1.0f);
	     xIntercept=xsource+slopeinv*(y-ysource);
	     ix=floor((xIntercept-x0)/dx);
	     if (ix == ixOld) {
	       if ((ix >= 0) && (ix < nx)) {
                 smat[ix*ny + iy] = smat[ix*ny + iy]+
	            val0*travPixlen;
		       }
		    }else{
	       xMid=dx*max(ix,ixOld)+xl;
	       xdist1=fabs(xMid-xIntOld);
	       xdist2=fabs(xIntercept-xMid);
	       frac1=xdist1/(xdist1+xdist2);
	       frac2=1.0f-frac1;
	       if ((ixOld >= 0) && (ixOld < nx)) {
                 smat[ixOld*ny + iy] = smat[ixOld*ny + iy]+
	            val0*frac1*travPixlen;
		       }
	       if ((ix>=0) && (ix<nx)) {
                 smat[ix*ny + iy] = smat[ix*ny + iy]+
	            val0*frac2*travPixlen;
		       }
		    }
	     ixOld=ix;
	     xIntOld=xIntercept;
	    }


		    }




	}
__syncthreads();
}

}




__global__ void multirayproj(int* intparms, float* floatparms,
                   float* sinovals, char* sinoinds, float* smat){




__shared__ int ngroups, nx,ny,nu;
__shared__ float x0,y0,dx,dy,du,u0;
__shared__ float xSourceCenter,ySourceCenter,xDetCenter,yDetCenter,eux,euy;




ngroups=intparms[1];
nx=intparms[2]; ny=intparms[3]; nu=intparms[4];

x0=floatparms[1]; y0=floatparms[2]; dx=floatparms[3]; dy=floatparms[4];
du=floatparms[5]; u0=floatparms[6];

xSourceCenter=floatparms[7]; ySourceCenter=floatparms[8];
xDetCenter=floatparms[9]; yDetCenter=floatparms[10];
eux=floatparms[11]; euy=floatparms[12];


int ng,iyOld, ixOld, ix, iy, jp, jpp;
float u, xbin, ybin, xsource, ysource, xl, yl, xdiff, ydiff, xad, yad, slope, slopeinv, x, y;
float travPixlen, raysum, yMid, xMid;
float yIntercept, yIntOld, ydist1, ydist2, frac1, frac2;
float xIntercept, xIntOld, xdist1, xdist2;

jpp =  blockDim.x*blockIdx.x + threadIdx.x;

for(ng=0; ng<ngroups; ng++){
jp = ng + ngroups*jpp;
if ((jp<nu) && (sinoinds[jp] == 1)){

   u = u0+(jp+0.5f)*du;
   xbin = xDetCenter + eux*u;
   ybin = yDetCenter + euy*u;
   xsource = xSourceCenter;
   ysource = ySourceCenter;


   xl=x0;
   yl=y0;

   xdiff=xbin-xsource;
   ydiff=ybin-ysource;
   xad=fabs(xdiff)*dy;
   yad=fabs(ydiff)*dx;

   if (xad>yad) {
      slope=ydiff/xdiff;
      travPixlen=dx*sqrt(1.0f+slope*slope);
      yIntOld=ysource+slope*(xl-xsource);
      iyOld=floor((yIntOld-y0)/dy);
      raysum=0.;
           for (ix=0; ix<nx; ix++) {
	     x=xl+dx*(ix + 1.0f);
	     yIntercept=ysource+slope*(x-xsource);
	     iy=floor((yIntercept-y0)/dy);
	     if (iy == iyOld) {
	       if ((iy >= 0) && (iy < ny)) {
		 raysum=raysum+travPixlen*smat[ix*ny + iy];
		       }
		    }else{
	       yMid=dy*max(iy,iyOld)+yl;
	       ydist1=fabs(yMid-yIntOld);
	       ydist2=fabs(yIntercept-yMid);
	       frac1=ydist1/(ydist1+ydist2);
	       frac2=1.0f-frac1;
	       if ((iyOld >= 0) && (iyOld < ny)) {
		 raysum=raysum+frac1*travPixlen*smat[ix*ny + iyOld];
		       }
	       if ((iy>=0) && (iy<ny)) {
		 raysum=raysum+frac2*travPixlen*smat[ix*ny + iy];
		       }
		    }
	     iyOld=iy;
	     yIntOld=yIntercept;
	    }
		    }else{


	   slopeinv=xdiff/ydiff;
	   travPixlen=dy*sqrt(1.0f+slopeinv*slopeinv);
	   xIntOld=xsource+slopeinv*(yl-ysource);
	   ixOld=floor((xIntOld-x0)/dx);
	   raysum=0.f;
           for (iy=0; iy<ny; iy++) {
	     y=yl+dy*(iy + 1.0f);
	     xIntercept=xsource+slopeinv*(y-ysource);
	     ix=floor((xIntercept-x0)/dx);
	     if (ix == ixOld) {
	       if ((ix >= 0) && (ix < nx)) {
		 raysum=raysum+travPixlen*smat[ix*ny + iy];
		       }
		    }else{
	       xMid=dx*max(ix,ixOld)+xl;
	       xdist1=fabs(xMid-xIntOld);
	       xdist2=fabs(xIntercept-xMid);
	       frac1=xdist1/(xdist1+xdist2);
	       frac2=1.0f-frac1;
	       if ((ixOld >= 0) && (ixOld < nx)) {
		 raysum=raysum+frac1*travPixlen*smat[ixOld*ny + iy];
		       }
	       if ((ix>=0) && (ix<nx)) {
                 raysum=raysum+frac2*travPixlen*smat[ix*ny + iy];
		       }
		    }
	     ixOld=ix;
	     xIntOld=xIntercept;
	    }


		    }

     sinovals[jp]=raysum;

	}
__syncthreads();
}

}



