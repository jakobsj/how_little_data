c This set of functions are written for proj. and back_proj.
c using a list of rays.

c the tomo.py library attaches a configuration specific
c function that converts a sino index to beam properties used in
c these routines.







       subroutine get_frames_param(du,xsource,ysource, 
     +           xDetCenter, yDetCenter, eux, euy,
     +           ewx, ewy, tcx, tcy, sd_dist, frame,
     +           sourcewidth, binwidth)

c This subroutine extracts the configuration data from the frame vector
      implicit real*8(a-h,o-z)

      Dimension frame(8)

         xSource=frame(1)
         ySource=frame(2)
        
         xDetCenter=frame(3)
         yDetCenter=frame(4)
        
         eux=frame(5)
         euy=frame(6)
        
         ewx=frame(7)
         ewy=frame(8)

c Calculate detector true coordinate center

         udr =(xDetCenter-xSource)*eux+
     +        (yDetCenter-ySource)*euy
       
         tcx = xDetCenter - eux*udr 
         tcy = yDetCenter - euy*udr 
        
c compute source-detector distance

         sd_dist=Sqrt((xSource-tcx)**2+
     +                (ySource-tcy)**2)

         sourcewidth=0
         binwidth=du

        
      return
      end


c ###########################################################


      subroutine ellproj(sinomat,indsino,frame_vectors,
     +           ns,nu,du,u0,ax,ay,x0,y0,gam,att)

c perfoming a projection of a ellips to a sinogram.
c This is a projection of a continueus object

      Implicit Real*8(a-h,o-z)

      Parameter (ilarge=10000)
      Dimension frame_vectors(ns,8)
      Dimension indsino(ns,nu)
      Dimension sinomat(ns,nu)

Cf2py intent(in) ns
Cf2py intent(in) nu
cf2py intent(in) frame_vectors
cf2py depend(ns) frame_vectors
Cf2py intent(out) sinomat
cf2py depend(ns) sinomat
cf2py depend(nu) sinomat
Cf2py intent(in) indsino
cf2py depend(ns) indsino
cf2py depend(nu) indsino

c loop thru sinograms indicator function
      Do ip=1,ns

       call  get_frames_param(du,xsource,ysource,
     +           xDetCenter, yDetCenter, eux, euy,
     +           ewx, ewy, tcx, tcy, sd_dist, frame_vectors(ip,:),
     +           sourcewidth, binwidth)

      Do jp=1,nu

      if (indsino(ip,jp).eq.1) then

         u = u0+(jp-0.5)*du 
         xbin = xDetCenter + eux*u
         ybin = yDetCenter + euy*u

      
c put source and detector points in ellipse frame
      xsp=x0 + (-x0 + xsource)*Cos(gam) - (y0 - ysource)*Sin(gam)
      ysp=y0 + (-y0 + ysource)*Cos(gam) + (x0 - xsource)*Sin(gam)

      xbp=x0 + (-x0 + xbin)*Cos(gam) - (y0 - ybin)*Sin(gam)
      ybp=y0 + (-y0 + ybin)*Cos(gam) + (x0 - xbin)*Sin(gam)


      pathLength2=(4*ax**2*ay**2*(ay**2*(xbp - xsp)**2 - 
     -      (xsp*(-y0 + ybp) + xbp*(y0 - ysp) + (ax - x0)*(ybp - ysp))*
     -       (xsp*(-y0 + ybp) + xbp*(y0 - ysp) - (ax + x0)*(ybp -
     -       ysp)))*((xbp - xsp)**2 + (ybp - ysp)**2))/
     -  (ay**2*(xbp - xsp)**2 + ax**2*(ybp - ysp)**2)**2

      if (pathLength2.lt.0.) pathLength2=0.

      sinomat(ip,jp)=sqrt(pathLength2)*att

      endif
      enddo
      enddo

      return
      end



     
      subroutine rayproj(sinomat,indsino,frame_vectors,ns,nu,du,u0,
     +       smat,dx,dy,x0,y0,nx,ny)

c Projecting an image to a sinogram.
c Note that this is a projection of a discreetized object

      Implicit Real*8(a-h,o-z)

      parameter (ilarge=10000)

      Dimension frame_vectors(ns,8)
      Dimension indsino(ns,nu)
      Dimension sinomat(ns,nu),smat(nx,ny)

Cf2py intent(in) ns
Cf2py intent(in) nu
Cf2py intent(in) nx
Cf2py intent(in) ny
cf2py intent(in) frame_vectors
cf2py depend(ns) frame_vectors
Cf2py intent(out) sinomat
cf2py depend(ns) sinomat
cf2py depend(nu) sinomat
Cf2py intent(in) indsino
cf2py depend(ns) indsino
cf2py depend(nu) indsino
Cf2py intent(in) smat
cf2py depend(nx) smat
cf2py depend(ny) smat


c loop thru sinograms indicator function
      Do ip=1,ns
         
        call  get_frames_param(du,xsource,ysource,
     +           xDetCenter, yDetCenter, eux, euy,
     +           ewx, ewy, tcx, tcy, sd_dist, frame_vectors(ip,:),
     +           sourcewidth, binwidth)

      Do jp=1,nu

      if (indsino(ip,jp).eq.1) then

         u = u0+(jp-0.5)*du 
         xbin = xDetCenter + eux*u
         ybin = yDetCenter + euy*u

      
     
      xl=x0
      yl=y0

      xdiff=xbin-xsource
      ydiff=ybin-ysource
      xad=abs(xdiff)*dy
      yad=abs(ydiff)*dx

      if (xad.gt.yad) then
            slope=ydiff/xdiff
            travPixlen=dx*sqrt(1.0+slope**2)
            yIntOld=ysource+slope*(xl-xsource)
            iyOld=int(ilarge+(yIntOld-y0)/dy)-ilarge
            total=0.
            do ix=1, nx
               x=xl+dx*ix
               yIntercept=ysource+slope*(x-xsource)
               iy=int(ilarge+(yIntercept-y0)/dy)-ilarge

            if (iy.eq.iyOld) then
                  if ((iy.ge.0).and.(iy.lt.ny))
     +               total=total+travPixlen*smat(ix,iy+1)
               else
                  yMid=dy*max(iy,iyOld)+yl
                  ydist1=abs(yMid-yIntOld)
                  ydist2=abs(yIntercept-yMid)
                  frac1=ydist1/(ydist1+ydist2)
                  frac2=1.0-frac1
                  if ((iyOld.ge.0).and.(iyOld.lt.ny)) 
     +                total=total+frac1*travPixlen*smat(ix,iyOld+1)
                  if ((iy.ge.0).and.(iy.lt.ny))
     +               total=total+frac2*travPixlen*smat(ix,iy+1)
               endif
               iyOld=iy
               yIntOld=yIntercept
            enddo
      else
            slopeinv=xdiff/ydiff
            travPixlen=dy*sqrt(1.0+slopeinv**2)
            xIntOld=xsource+slopeinv*(yl-ysource)
            ixOld=int(ilarge+(xIntOld-x0)/dx)-ilarge
            total=0.
            do iy=1, ny
               y=yl+dy*iy
               xIntercept=xsource+slopeinv*(y-ysource)
               ix=int(ilarge+(xIntercept-x0)/dx)-ilarge

               if (ix.eq.ixOld) then
                  if ((ix.ge.0).and.(ix.lt.nx))
     +               total=total+travPixlen*smat(ix+1,iy)
               else
                  xMid=dx*max(ix,ixOld)+xl
                  xdist1=abs(xMid-xIntOld)
                  xdist2=abs(xIntercept-xMid)
                  frac1=xdist1/(xdist1+xdist2)
                  frac2=1.0-frac1
                  if ((ixOld.ge.0).and.(ixOld.lt.nx)) 
     +                total=total+frac1*travPixlen*smat(ixOld+1,iy)
                  if ((ix.ge.0).and.(ix.lt.nx))
     +               total=total+frac2*travPixlen*smat(ix+1,iy)
               endif
               ixOld=ix
               xIntOld=xIntercept
            enddo
      endif

      sinomat(ip,jp)=total

      end if
      end do
      enddo

      return
      end

     


c ###############################################################     

      subroutine backproject(
     +       sinomat,indsino,frame_vectors,ns,ds,nu,du,u0,
     +       smat,dx,dy,x0,y0,nx,ny)
     

c Backprojecting from sinogram to an image. This algorithm is ray deriven.

      Implicit Real*8(a-h,o-z)

      Parameter (ilarge=10000,nraypix=10000)
      Dimension frame_vectors(ns,8)
      Dimension indsino(ns,nu)
      Dimension sinomat(ns,nu),smat(nx,ny)
      Dimension raypix(nraypix),iraypix(nraypix),jraypix(nraypix)
Cf2py intent(in) ns
Cf2py intent(in) ds
Cf2py intent(in) nu
Cf2py intent(in) nx
Cf2py intent(in) ny
cf2py intent(in) frame_vectors
cf2py depend(ns) frame_vectors
Cf2py intent(in) sinomat
cf2py depend(ns) sinomat
cf2py depend(nu) sinomat
Cf2py intent(in) indsino
cf2py depend(ns) indsino
cf2py depend(nu) indsino
Cf2py intent(out) smat
cf2py depend(nx) smat
cf2py depend(ny) smat

c fnr is a normalization factor of the transformation
c between sinogram and image
c      fnr=((ds*du)/(dx*dy))
      fnr=1.0
      
      Do ix=1,nx
         Do iy =1, ny
            smat(ix,iy) = 0.
         End Do
      End Do
  
c loop thru sinograms indicator function
      Do ip=1,ns
     
      call get_frames_param(du,xsource,ysource, 
     +           xDetCenter, yDetCenter, eux, euy,
     +           ewx, ewy, tcx, tcy, sd_dist, frame_vectors(ip,:),
     +           sourcewidth, binwidth)

      Do jp=1,nu

      if (indsino(ip,jp).eq.1) then

         u = u0+(jp-0.5)*du
         xbin = xDetCenter + eux*u
         ybin = yDetCenter + euy*u
         val=sinomat(ip,jp)

      
      xl=x0
      yl=y0

      xdiff=xbin-xsource
      ydiff=ybin-ysource
      xad=abs(xdiff)*dy
      yad=abs(ydiff)*dx

      if (xad.gt.yad) then
            slope=ydiff/xdiff
            travPixlen=dx*sqrt(1.0+slope**2)
            yIntOld=ysource+slope*(xl-xsource)
            iyOld=int(ilarge+(yIntOld-y0)/dy)-ilarge
c            total=0.
            raysum=0.
            irp=0
            do ix=1, nx
               x=xl+dx*ix
               yIntercept=ysource+slope*(x-xsource)
               iy=int(ilarge+(yIntercept-y0)/dy)-ilarge

               if (iy.eq.iyOld) then
                  if ((iy.ge.0).and.(iy.lt.ny)) then
                     irp=irp+1
                     raypix(irp)=travPixlen
                     iraypix(irp)=ix
                     jraypix(irp)=iy+1
                     raysum=raysum+travPixlen*smat(ix,iy+1)
c                     total=total+travPixlen*travPixlen
                  endif
               else
                  yMid=dy*max(iy,iyOld)+yl
                  ydist1=abs(yMid-yIntOld)
                  ydist2=abs(yIntercept-yMid)
                  frac1=ydist1/(ydist1+ydist2)
                  frac2=1.0-frac1
                  if ((iyOld.ge.0).and.(iyOld.lt.ny)) then
                     irp=irp+1
                     raypix(irp)=frac1*travPixlen
                     iraypix(irp)=ix
                     jraypix(irp)=iyOld+1
                     raysum=raysum+frac1*travPixlen*smat(ix,iyOld+1)
c                     total=total+frac1*frac1*travPixlen*travPixlen
                  endif
                  if ((iy.ge.0).and.(iy.lt.ny)) then
                     irp=irp+1
                     raypix(irp)=frac2*travPixlen
                     iraypix(irp)=ix
                     jraypix(irp)=iy+1
                     raysum=raysum+frac2*travPixlen*smat(ix,iy+1)
c                     total=total+frac2*frac2*travPixlen*travPixlen
                  endif
               endif
               iyOld=iy
               yIntOld=yIntercept
            enddo
      else
            slopeinv=xdiff/ydiff
            travPixlen=dy*sqrt(1.0+slopeinv**2)
            xIntOld=xsource+slopeinv*(yl-ysource)
            ixOld=int(ilarge+(xIntOld-x0)/dx)-ilarge
c            total=0.
            raysum=0.
            irp=0
            do iy=1, ny
               y=yl+dy*iy
               xIntercept=xsource+slopeinv*(y-ysource)
               ix=int(ilarge+(xIntercept-x0)/dx)-ilarge

               if (ix.eq.ixOld) then
                  if ((ix.ge.0).and.(ix.lt.nx)) then
                     irp=irp+1
                     raypix(irp)=travPixlen
                     iraypix(irp)=ix+1
                     jraypix(irp)=iy
                     raysum=raysum+travPixlen*smat(ix+1,iy)
c                     total=total+travPixlen*travPixlen
                  endif
               else
                  xMid=dx*max(ix,ixOld)+xl
                  xdist1=abs(xMid-xIntOld)
                  xdist2=abs(xIntercept-xMid)
                  frac1=xdist1/(xdist1+xdist2)
                  frac2=1.0-frac1
                  if ((ixOld.ge.0).and.(ixOld.lt.nx)) then
                     irp=irp+1
                     raypix(irp)=frac1*travPixlen
                     iraypix(irp)=ixOld+1
                     jraypix(irp)=iy
                     raysum=raysum+frac1*travPixlen*smat(ixOld+1,iy)
c                     total=total+frac1*frac1*travPixlen*travPixlen
                  endif
                  if ((ix.ge.0).and.(ix.lt.nx)) then
                     irp=irp+1
                     raypix(irp)=frac2*travPixlen
                     iraypix(irp)=ix+1
                     jraypix(irp)=iy
                     raysum=raysum+frac2*travPixlen*smat(ix+1,iy)
c                     total=total+frac2*frac2*travPixlen*travPixlen
                  endif
               endif
               ixOld=ix
               xIntOld=xIntercept
            enddo
      endif

      nrp=irp
      do irp=1,nrp
         smat(iraypix(irp),jraypix(irp))=
     +       smat(iraypix(irp),jraypix(irp))+ (val*raypix(irp))*fnr
      enddo

      endif
      enddo
      enddo

      return
      end



C CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c Here starts the subroutins of the pixel_driven_backproj
      
      subroutine find_ray_elongation_const(xsource,ysource,
     +     pix_x,pix_y,
     +     xDetCenter,yDetCenter,ewx,ewy,elongation_const)

        implicit real*8(a-h,o-z)
        real*8 xsource,ysource,pix_x,pix_y, diff_x_source_det
        real*8 diff_y_source_det,diff_x_source_pix,diff_y_source_pix
        real*8 xDetCenter,yDetCenter,ewx,ewy,elongation_const
        
        diff_x_source_det = (xDetCenter-xsource)*ewx
        diff_y_source_det = (yDetCenter-ysource)*ewy
        diff_x_source_pix = (pix_x-xsource)*ewx
        diff_y_source_pix = (pix_y-ysource)*ewy        
        elongation_const = ((diff_x_source_det+diff_y_source_det)/
     +       (diff_x_source_pix+diff_y_source_pix))

        return
        end


       subroutine find_interception(xsource,ysource,pix_x,pix_y,
     +          elongation_const, Det_int_x,Det_int_y)

        real*8 xsource,ysource,pix_x,pix_y
        real*8 elongation_const
        real*8 Det_int_x,Det_int_y
              
        Det_int_x = xsource+elongation_const*(pix_x-xsource)
        Det_int_y = ysource+elongation_const*(pix_y-ysource)

        return
        end


       subroutine is_intercepting(Det_int_x,Det_int_y,xDetCenter,
     +             yDetCenter,det_length,eux,euy,
     +             intercepting,upos)

        real*8 Det_int_x,Det_int_y,xDetCenter,yDetCenter,det_length
        real*8 eux,euy,upos
        logical intercepting

         
         upos = ((Det_int_x-xDetCenter)*eux
     +       +(Det_int_y-yDetCenter)*euy)
        if ((det_length/2).lt.abs(upos)) then
           
           intercepting = .FALSE.

        else
           intercepting = .TRUE.

        endif
        return
        end


       subroutine compute_det_value(sinomat,indsino,ip,
     +     upos,u0,nu,ns,binwidth,
     +     det_value)

        implicit real*8(a-h,o-z)
        real*8 u0,binwidth
        real*8 det_value,frac,bin_loc
        integer nbin1,nbin2,ip,nu,ns
        Dimension sinomat(ns,nu)
        Dimension indsino(ns,nu)

        bin_loc = (upos-u0)/binwidth +0.5
        nbin1 = int(bin_loc)
        nbin2 = nbin1+1 

        if ((nbin1.ge.1).and.(nbin1.le.nu-1)) then
              frac= bin_loc - nbin1
              det_value=frac*sinomat(ip,nbin2)*indsino(ip,nbin2)
     +             +(1.-frac)*sinomat(ip,nbin1)*indsino(ip,nbin1)
        else
           det_value = 0.0
        endif
        return
        end

       subroutine put_det_value(sinomat,indsino,ip,
     +     upos,u0,nu,ns,binwidth,
     +     pix_value)

        implicit real*8(a-h,o-z)
        real*8 u0,binwidth
        real*8 pix_value,frac,bin_loc
        integer nbin1,nbin2,ip,nu,ns
        Dimension sinomat(ns,nu)
        Dimension indsino(ns,nu)

        bin_loc = (upos-u0)/binwidth +0.5
        nbin1 = int(bin_loc)
        nbin2 = nbin1+1 

        if ((nbin1.ge.1).and.(nbin1.le.nu-1)) then
              frac= bin_loc - nbin1
              sinomat(ip,nbin2) = sinomat(ip,nbin2) +
     +           frac*pix_value*indsino(ip,nbin2)
              sinomat(ip,nbin1) = sinomat(ip,nbin1) +
     +           (1.-frac)*pix_value*indsino(ip,nbin1)
        endif
        return
        end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C Here starts the main function pixel_driven_backproj

       subroutine pixel_driven_backproj(
     +       sinomat,indsino,frame_vectors,ds,ns,nu,du,u0,
     +       smat,dx,dy,x0,y0,nx,ny,fov,xc,yc)

c This routine does image reconstruction by computing the total contribution 
c from the sinogram to every pixel. It uses linear interpolation to determine
c the value contributed by each ray.

      implicit real*8(a-h,o-z)
      Parameter (ilarge=10000)
      Dimension frame_vectors(ns,8)
      Dimension indsino(ns,nu)
      Dimension sinomat(ns,nu),smat(nx,ny)
      real*8 ds,du,u0,dx,dy,x0,y0,det_length
      real*8 fov,xc,yc,frad, upos
      integer ns,nu,nx,ny,ip,ix,iy
      logical intercepting

Cf2py intent(in) ns
Cf2py intent(in) nu
Cf2py intent(in) nx
Cf2py intent(in) ny
cf2py intent(in) frame_vectors
cf2py depend(ns) frame_vectors
Cf2py intent(in) sinomat
cf2py depend(ns) sinomat
cf2py depend(nu) sinomat
Cf2py intent(in) indsino
cf2py depend(ns) indsino
cf2py depend(nu) indsino
Cf2py intent(in,out) smat
cf2py depend(nx) smat
cf2py depend(ny) smat
      
      det_length = du*nu
      

        Do iy=1,ny          
           pix_y=y0+dy*(iy-0.5)

           Do ix=1,nx
              pix_x=x0+dx*(ix-0.5)
c     Check whether the pixel is within the FOV
            frad=sqrt((pix_x-xc)**2+(pix_y-yc)**2)
            if (frad.lt.fov) then
      
      Do ip=1,ns
       
       call get_frames_param(du,xsource,ysource,
     +        xDetCenter, yDetCenter, eux, euy,
     +        ewx, ewy, tcx, tcy, sd_dist, frame_vectors(ip,:),
     +        sourcewidth, binwidth)
      

c    find source-pixel vector elongation constant in order
c    to find the ray vector
              call find_ray_elongation_const(xsource,ysource,
     +     pix_x,pix_y,
     +     xDetCenter,yDetCenter,ewx,ewy,elongation_const)

c    find interseption point of the ray and the detector line
              call find_interception(xsource,ysource,pix_x,pix_y,
     +          elongation_const, Det_int_x,Det_int_y)

c    Check whether the interseption point is within the detector physical size
              call is_intercepting(Det_int_x,Det_int_y,xDetCenter,
     +             yDetCenter,det_length,eux,euy,
     +             intercepting,upos)

              if (intercepting) then
                 
c    compute the value of the sinogram by interpolation
              call compute_det_value(sinomat,indsino,ip,
     +     upos,u0,nu,ns,binwidth,
     +     det_value)
                
                smat(ix,iy)=smat(ix,iy)+det_value*ds

              endif
      enddo
            endif
           enddo
        enddo

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C  weighted fbp for linear-detector fan-beam recon. based on the pixel_bp above.

       subroutine weighted_pixel_driven_backproj(
     +       sinomat,indsino,frame_vectors,ds,ns,nu,du,u0,
     +       smat,dx,dy,x0,y0,nx,ny,fov,xc,yc)

c This routine does image reconstruction by computing the total contribution 
c from the sinogram to every pixel. It uses linear interpolation to determine
c the value contributed by each ray.

      implicit real*8(a-h,o-z)
      Parameter (ilarge=10000)
      Dimension frame_vectors(ns,8)
      Dimension indsino(ns,nu)
      Dimension sinomat(ns,nu),smat(nx,ny)
      real*8 ds,du,u0,dx,dy,x0,y0,det_length, s, upos
      real*8 fov,xc,yc,frad, fphi, bigd,bigu,fbweight, pio2
      integer ns,nu,nx,ny,ip,ix,iy
      logical intercepting

Cf2py intent(in) ns
Cf2py intent(in) nu
Cf2py intent(in) nx
Cf2py intent(in) ny
cf2py intent(in) frame_vectors
cf2py depend(ns) frame_vectors
Cf2py intent(in) sinomat
cf2py depend(ns) sinomat
cf2py depend(nu) sinomat
Cf2py intent(in) indsino
cf2py depend(ns) indsino
cf2py depend(nu) indsino
Cf2py intent(in,out) smat
cf2py depend(nx) smat
cf2py depend(ny) smat
      
      det_length = du*nu
      pio2=Atan(1.0)*2.0
      
      Do ip=1,ns
c assumming that first view is at 0.
c       print *, "w-bp view: ",ip
       s = (ip-1.)*ds
       
       call get_frames_param(du,xsource,ysource,
     +        xDetCenter, yDetCenter, eux, euy,
     +        ewx, ewy, tcx, tcy, sd_dist, frame_vectors(ip,:),
     +        sourcewidth, binwidth)

        bigd = sqrt(xsource*xsource + ysource*ysource)
      

        Do iy=1,ny          
           pix_y=y0+dy*(iy-0.5)

           Do ix=1,nx
              pix_x=x0+dx*(ix-0.5)
c     Check whether the pixel is within the FOV
            frad=sqrt((pix_x-xc)**2+(pix_y-yc)**2)
            fphi= atan2(pix_y,pix_x)
            if (frad.lt.fov) then

c eq. (106) kak-slaney
              bigu = (bigd+frad*sin(s-fphi - pio2))/bigd
              fbweight = 1./(bigu*bigu)


c    find source-pixel vector elongation constant in order
c    to find the ray vector
              call find_ray_elongation_const(xsource,ysource,
     +     pix_x,pix_y,
     +     xDetCenter,yDetCenter,ewx,ewy,elongation_const)

c    find interseption point of the ray and the detector line
              call find_interception(xsource,ysource,pix_x,pix_y,
     +          elongation_const, Det_int_x,Det_int_y)

c    Check whether the interseption point is within the detector physical size
              call is_intercepting(Det_int_x,Det_int_y,xDetCenter,
     +             yDetCenter,det_length,eux,euy,
     +             intercepting,upos)

              if (intercepting) then
                 
c    compute the value of the sinogram by interpolation
              call compute_det_value(sinomat,indsino,ip,
     +     upos,u0,nu,ns,binwidth,
     +     det_value)
                
                smat(ix,iy)=smat(ix,iy)+fbweight*det_value*ds

              endif
            endif
           enddo

        enddo
      enddo

      return
      end

