! -------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      
      module m_eigen 
!
!  Call sequence of routines to compute ordered, normalized, modified
!  EOFs and associated eigenvalues.  Also, for symmetric matrices allows
!  for computation of inverse. 
!
      use precision

      implicit none


      Private

      public eigen_main

      real(r8), parameter  :: zero=0.d0
      real(r8), parameter  :: one=1.d0
      real(r8), parameter  :: two=2.d0
      real(r8), parameter  :: three=3.d0

! !REVISION HISTORY:
!
!  04Jun2003  R. Errico  Initial algorithm

          
      CONTAINS

! -------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------

      subroutine eigen_main (evectrs, evalues, matrix, weights, nk, &
                            testvalue, purpose, name, matrixR)

      external ssyev 
      external dsyev 
      
      character(len=*), intent(in) :: purpose
      character(len=*), intent(in) :: name
      integer,  intent(in)   :: nk
      real(r8), intent(in)   :: matrix(nk,nk)
      real(r8), intent(in)   :: weights(nk)
      real(r8), intent(in)   :: testvalue
      real(r8), intent(out)  :: evectrs(nk,nk)
      real(r8), intent(out)  :: evalues(nk)
      real(r8), intent(out)  :: matrixR(nk,nk)

!
      integer     :: nwork1  
      integer     :: info 
      real(r8)    :: mcopy(nk,nk) 
      real(r8)    :: work1(nk*nk*2)

      nwork1=nk*nk*2

      if (purpose(2:2) == 'S') then ! symmetric input matrix
        evectrs=matrix 
        if (r8 == 4) then
          call ssyev ( purpose(1:1), 'U', nk, evectrs, nk, evalues, & 
                        work1, nwork1, info )
        elseif (r8 == 8) then
          call dsyev ( purpose(1:1), 'U', nk, evectrs, nk, evalues, & 
                        work1, nwork1, info )
        else
          print *,'Invalid precision specified in EigenMain; r8=',r8
          stop
        endif
      else
        print *,'Invalid purpose in EigenMain; purpose=',purpose
        stop
      endif
 
      call Reorder (evectrs, evalues, nk, purpose(1:1) )

      call Normalize (evectrs, weights, nk)
      
      if (testvalue > 0.d0) then 
         call CheckEigen (matrix, evectrs, evalues, nk, name, testvalue )
         if ( purpose(2:2) == 'S' ) then
           call CheckOrtho (evectrs, weights, nk, name, testvalue)
         endif
      endif

      if ( ( purpose(2:2) == 'S') .and. (purpose(3:3) == 'I' ) ) then
        call Inverse ( matrixR, evectrs, evalues, nk, name )
        if (testvalue > 0.d0) then 
          call CheckInverse (matrixR, matrix, nk, name, testvalue)
        endif
      endif
!

      end subroutine eigen_main
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      Subroutine Reorder (evectrs,evalues,nk,cv)
!
!  Reorder eigenvalues and corresponding eigenvectors from largest 
!  eigenvalue to smallest. 
!
      character(len=*), intent(in)  :: cv
      integer, intent(in)  :: nk 
      real(r8), intent(inout)  :: evectrs(nk,nk)
      real(r8), intent(inout)  :: evalues(nk)
!
      integer  :: k, kx,i
      real(r8) :: copy, emax     
!     
      do k=1,nk-1
        emax=-1.e-30 
        do i=k,nk
          if (evalues(i).gt.emax) then
            emax=evalues(i)
            kx=i
          endif 
        enddo
        copy=evalues(kx)
        evalues(kx)=evalues(k)
        evalues(k)=copy
        if (cv(1:1) == 'V') then
          do i=1,nk
            copy=evectrs(i,kx)
            evectrs(i,kx)=evectrs(i,k)
            evectrs(i,k)=copy
          enddo
        endif
      enddo
!
      end subroutine reorder

!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      Subroutine CheckOrtho (evectrs, weights, nk, name, testvalue )
!
! Check whether vectors are orthogonal 
!
      character(len=*), intent(in) :: name
      integer,  intent(in) :: nk
      real(r8), intent(in) :: evectrs(nk,nk)
      real(r8), intent(in) :: weights(nk)
      real(r8), intent(in) :: testvalue
!
      integer   :: i, j, k
      real(r8)  :: sum
!        
      do i=1,nk
        do j=1,nk
          sum=0.d0
          if ( i==j ) sum=-1.d0
          do k=1,nk
            sum=sum+evectrs(k,i)*evectrs(k,j)*weights(k)
          enddo
          if ( abs(sum) > testvalue ) then
              print *,' CheckOrtho fails for ',name,' i,j,sum = ',i,j,sum
          endif
        enddo  
      enddo 

      end subroutine CheckOrtho
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      Subroutine Normalize (evectrs, weights, nk)
!
!  Normalize each eigenvector such that the sum of its squared 
!  components times corresonding weight is 1
!
      integer,  intent(in)    :: nk
      real(r8), intent(in)    :: weights(nk)
      real(r8), intent(inout) :: evectrs(nk,nk)
!
      integer  :: k, i
      real(r8) :: sum  
     
      do k=1,nk
        sum=0.d0   
        do i=1,nk
          sum=sum+evectrs(i,k)*evectrs(i,k)
        enddo
        sum=1./sqrt(sum)
        evectrs(:,k)=evectrs(:,k)*sum
      enddo
!
      end subroutine normalize
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      Subroutine CheckEigen (matrix, evectrs, evalues, nk, name, &
                             testvalue)
 
      character(len=*), intent(in) :: name
      integer,  intent(in) :: nk
      real(r8), intent(in) :: matrix(nk,nk)
      real(r8), intent(in) :: evectrs(nk,nk)
      real(r8), intent(in) :: evalues(nk)
      real(r8), intent(in) :: testvalue
!
      integer  :: i, j, k
      real(r8) :: vmax, dmax, work 

!
      do k=1,nk
        vmax=0.
        dmax=0.
        do j=1,nk
          work=-evectrs(j,k)*evalues(k)
          if (abs(work).gt.vmax) vmax=abs(work)
          do i=1,nk
            work=work+matrix(j,i)*evectrs(i,k)
          enddo
          if (abs(work).gt.dmax) dmax=abs(work)
        enddo
        if (dmax.gt.vmax*testvalue) print *,'CheckEigen', &
                                    name,' ',k,dmax,vmax
      enddo
!
      end subroutine CheckEigen
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      Subroutine Inverse ( matrixR, evectrs, evalues, nk, name )      

      character(len=*), intent(in) :: name
      integer,  intent(in)  :: nk
      real(r8), intent(in)  :: evectrs(nk,nk)
      real(r8), intent(in)  :: evalues(nk)
      real(r8), intent(out) :: matrixR(nk,nk)

      integer               :: i, j, k
      real                  :: minvalue 
!
! Invert symmetric matrix using eigenvalues
!

      minvalue=1.e30
      do i=1,nk      
        if (abs(evalues(i)) < minvalue) minvalue=abs(evalues(i))
      enddo
      if (minvalue < 1.d-30) then 
        print *,'matrix has eigenvalues as abs values as low as ',minvalue
        print *,'inverse not computed'
        stop
      endif   

      do k=1,nk
        do j=1,nk
          matrixR(k,j)=0.d0
          do i=1,nk 
            matrixR(k,j)=matrixR(k,j)+evectrs(k,i)*evectrs(j,i) &
                         /evalues(i)
          enddo
        enddo
      enddo
!
      end subroutine Inverse
!
!  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!
      Subroutine CheckInverse (matrixR, matrix, nk, name, testvalue)
      character(len=*), intent(in) :: name
      integer,  intent(in)  :: nk
      real(r8), intent(in)  :: testvalue
      real(r8), intent(in)  :: matrixR(nk,nk)
      real(r8), intent(in)  :: matrix(nk,nk)

      integer               :: i, j, k
      real(r8)              :: cmax1, cmax2, cmax, s  
!
! check inverse of matrix
!
      cmax1=0.
      cmax2=0.
      do k=1,nk
        do j=1,nk
          if (cmax1.lt.abs( matrix(j,k))) cmax1=abs( matrix(j,k))
          if (cmax2.lt.abs(matrixR(j,k))) cmax2=abs(matrixR(j,k))
        enddo
      enddo
      cmax=cmax1*cmax2*nk*testvalue
!
      do k=1,nk
        do j=1,nk
          if (k.eq.j) then
            s=-1.d0
          else
            s=0.d0
          endif
          do i=1,nk
            s=s+matrixR(k,i)*matrix(i,j)
          enddo
          if (abs(s) > cmax) print *,'CheckInverse Failed ',name,' ' &
                            ,k,j,s,cmax
        enddo       
      enddo
!
      end subroutine CheckInverse

      end module m_eigen




