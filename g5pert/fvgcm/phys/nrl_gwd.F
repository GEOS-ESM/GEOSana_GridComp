      subroutine nrl_gwd(im,km,kbot,dta,cp,grav,pl,dpl,
     &                   pk,tt,hprime,ut,vt)                                               
 
! Gravity wave drag parameterization based on Palmer, T. N., G. J. Shutts, and
!    R. Swinbank, 1986: Alleviation of a systematic westerly bias in general
!    circulation and numerical weather prediction models through an orographic
!    gravity wave drag parameterization. Quart. J. Roy. Met. Soc, 112, 1001-1039.

! The procedure is that vector wind field (ut,vt) is adjusted (changed) in response
!    to the computed gravity wave drag.
 
! See the latest edition (really out of date) of
!                         The NOGAPS Forecast Model:
!                         A Technical Description
!                         by T.F. Hogan, T.E. Rosmond, and R. Gelaro
! PARAMETERS:
!      Name          Type      Usage            Description
!   ----------    ----------  -------  ----------------------------
!      im          Integer     Input    x-vector length,
!      km          Integer     Input    number of model layers,
!      kbot        Integer     Input    level where gravity wave is induced,
!      dta         Real        Input    leapfrog time step (seconds),
!      cp          Real        Input    specific heat of air (J/K kg),
!      grav        Real        Input    acceleration of gravity (m/s**2),
!      pl          Real        Input    full level pressures levels in pascal (im,km),
!      dpl         Real        Input    pressure thickness in pascal (im,km),
!      tt          Real        Input    virtual potential temperature in K (im,km), 
!      pk          Real        Input    Exner function (p/100000)**capa (im,km),
!      hprime      Real        Input    standard deviation of the terrain height in m (im),
!      ut          Real        In/Out   east-west wind in m/s (im,km), 
!      vt          Real        In/Out   north-south wind in m/s (im,km), 
                                                                                          
      implicit none                                                                       
 
! Input
      integer kbot,km,im
      real pl(im,km)
      real tt(im,km)
      real pk(im,km)
      real dpl(im,km)
      real hprime(im)
      real cp,dta,grav

! Input/Output
      real ut(im,km)
      real vt(im,km)

! Local
      real uhat(im,kbot),bvf(im,kbot),uvlm(im)
      real robvuh(im,kbot),gwtau(im)
      integer i,k

      real ak, csat, dh, eps, epssq
      real ri, tem,temx                                                   
      real uz0, u0               ! De-singularization const by SJL
      parameter (uz0 = 1.E-4, u0 = 1.e-4)
      data eps/1.0E-6/, epssq/1.0E-4/, ak/2.5E-5/                                           

!  Compute magnitude of even level wind component parallel to surface                     
!  wind.
      do i=1,im
         uvlm(i) = 0.5 / max( u0, 
     &           sqrt(ut(i,kbot)*ut(i,kbot)+vt(i,kbot)*vt(i,kbot)))
      enddo

      do 130 k=1,kbot                                                                     
        do i=1,im
          uhat(i,k)=max(eps,uvlm(i)*((ut(i,k)+ut(i,k+1))*ut(i,km)+                  
     &              (vt(i,k)+vt(i,k+1))*vt(i,km)))                                  
!  Compute Brunt-Vaisala frequency and Richardson number at half levels.
          csat=0.5*(tt(i,k)+tt(i,k+1))                                                
          dh = cp*(pk(i,k+1)-pk(i,k))*csat                                              
          bvf(i,k) =grav*sqrt(max(eps,(tt(i,k)-tt(i,k+1))/(dh*csat)))                  
          bvf(i,k) =max(bvf(i,k), epssq)                                                    
          robvuh(i,k)=(pl(i,k+1)-pl(i,k))*bvf(i,k)*uhat(i,k) / dh
        enddo
 130  continue                                                                            
                                                                                          
!  Displacement height and gravity wave stress for kbot surface                           
!  also add gw drag contribution to tridiagonal coeffs for kbot layer                     
                                                                                          
!  Displacement height is limited by requiring froude number in                           
!  kbot layer not exceed one.                                                             

      temx=1.0/grav                                                                       
      do 140 i=1,im
        dh = min(hprime(i),uhat(i,kbot)/bvf(i,kbot),400.0)                           
        gwtau(i) = ak*robvuh(i,kbot)*dh*dh                                                  
 140  continue                                                                            

!  Starting at even level above lowest layer where gw drag occurs, scan                   
!  upward looking for layers with positive ri no. lt 0.25 - where they                    
!  occur, apply Lindzen's saturation hypothosis to find displacement                      
!  height which brings ri no. up to 0.25.  Use this adjusted                              
!  displacement height to compute gw momentum flux on that even level.                    
!  drag is proportional to momentum flux divergence.                                      

! Limit the amount of gravity wave breaking below 450 mb - engineering
! fix to try to prevent too much breaking in the lower atmosphere.

      do 160 k=kbot-1,1,-1
        do 150 i=1,im
          tem =0.5*(tt(i,k)+tt(i,k+1))                                                 
          tem =(cp/grav)*(pk(i,k+1)-pk(i,k))*tem                                       
! SJL: ri could be zero
          ri=(ut(i,k)-ut(i,k+1))**2+(vt(i,k)-vt(i,k+1))**2 
          csat = (bvf(i,k)*tem) / max (uz0, sqrt(ri) )
          dh = gwtau(i) / (ak*robvuh(i,k))
          dh=sqrt(dh) 
          ri=csat*csat
          tem =bvf(i,k)*dh/uhat(i,k)
          tem =ri*(1.0-tem)/(1.0+csat*tem)**2 
          csat=2.0+1.0/csat              
          csat=2.0*sqrt(csat)-csat      
          if(tem.lt.0.25)dh=csat*uhat(i,k)/bvf(i,k) 
          dh=max(dh,0.0)                           
          if((pl(i,k) >= 45000.0).and.(uhat(i,k).le. eps)) dh=0.5*dh 
          if((pl(i,k) <  45000.0).and.(uhat(i,k).le. eps)) dh=0.0   
          tem = ak*robvuh(i,k)*dh*dh  
          csat=grav/dpl(i,k+1)     
          robvuh(i,k+1) = gwtau(i)-tem 
          gwtau(i)= tem               
          tem =dta*robvuh(i,k+1)*csat
          csat=ut(i,k+1)*ut(i,k+1)+vt(i,k+1)*vt(i,k+1)
          csat=sqrt(csat)                            
          tem =uvlm(i)*min(tem, csat)                 
          ut(i,k+1)=ut(i,k+1)-2.0*(ut(i,km)*tem)
          vt(i,k+1)=vt(i,k+1)-2.0*(vt(i,km)*tem)
 150    continue
 160  continue 
      return 
      end
