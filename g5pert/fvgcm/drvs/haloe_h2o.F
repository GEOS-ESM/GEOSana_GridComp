!_RT #include <params.h>
      subroutine haloe_h2o(q_unit, q_flnm, nymd, nhms, kw_top, kw_bot,
     &                     press, q_ref, jfirst, jlast)

# if defined ( SPMD ) 
      use mod_comm, only: gid, mp_barrier, mp_bcst_real, mp_bcst_r2d
      use mod_comm, only: jnp, nl
#endif

      implicit none

      integer   l_data
      parameter (l_data = 25)

      real vol2m
      parameter ( vol2m = 0.625E-6 )

! input
      integer jfirst, jlast
      character*(*) q_flnm
      integer q_unit
      integer nymd
      integer nhms
      integer kw_top, kw_bot
      real    press(nl)

! output
      real q_ref(nl, jfirst:jlast)

! local
      real p_data(l_data)

      real q1(jnp, l_data)
      real q2(jnp, l_data)
 
! Save vertically interpolated data
!_RT  real qm1(nl, jnp)
!_RT  real qm2(nl, jnp)
!_RT  save qm1, qm2
      real, allocatable, save :: qm1(:,:), qm2(:,:)

      integer  md(12)                             ! mn_day
      integer hms(12)                             ! hr_min_sec

      real    d_time(12)                          ! data time
      real    m_time                              ! model time
      integer m, m2
      integer n, nsecf


# if defined ( SPMD ) 
#else
      integer gid
      parameter ( gid = 0 )
#endif

      real*4 qt(jnp, l_data)
      integer j, k
      real lnp_m(nl)
      real lnp_d(l_data)
      real vmax, qmax,  qmin

      real res, fac, pi
      logical first

      data p_data / 1., 1.5, 2.2, 3.2, 4.6, 6.8, 10., 15., 22.,
     &             32.,  46., 68., 100., 150., 220., 320., 460.,
     &             680., 1000., 1470., 2150., 3160., 4630.,
     &             6810., 10000./
! Date for  non-leap-year
      data md  / 19990116,   19990215,  19990316,  19990416,   19990516,
     &           19990616,   19990716,  19990816,  19990916,   19991016,
     &           19991116,   19991216 /

      data hms /120000,      0, 120000,     0, 120000, 
     &               0, 120000, 120000,     0, 120000,
     &               0, 120000                 /

      data first / .true./
      save m2
      save d_time
      
      nsecf(n)  = n/10000*3600 + mod(n,10000)/100* 60 + mod(n,100)

      if ( first ) then
           allocate( qm1(nl, jnp) )
           allocate( qm2(nl, jnp) )
      endif

      if ( gid == 0 ) then

      call mcalday(nymd, nhms, m_time)

! Correct for no leap yr climatology
      if ( m_time .lt. 0. ) then                          ! leap year
           m_time = abs( m_time )
           if( mod(nymd,10000) .gt. 228 ) m_time = m_time - 1.
      endif

      if(first) then

         write(6,*) 'Unit to read H2O data is ', q_unit
         open (unit=q_unit,file=q_flnm,status='old',
     &         form='unformatted',access='direct',recl=jnp*l_data*4)
         write(6,*) 'File for H2O is opened.'

         read (q_unit, rec=1) qt

                 do k=1,l_data
                    do j=1,jnp
                       q2(j,k) = qt(j,k)
                    enddo
                 enddo
         m = 1

! Compute d_time(12)
         
      do m=1,12
         call mcalday(md(m), hms(m), d_time(m))
      enddo

      if( m_time .le. d_time(1) ) then

! Model time is early January.

           read (q_unit, rec=12) qt

                 do k=1,l_data
                    do j=1,jnp
                       q1(j,k) = qt(j,k)
                    enddo
                 enddo
             m2 = 1
             write(6,*) 'Dec and Jan H2O data will be used'
      else

! Model time is beyond 12Z,0116

          do m=2,12

                 do k=1,l_data
                    do j=1,jnp
                       q1(j,k) = q2(j,k)
                    enddo
                 enddo

             read(q_unit, rec=m) qt

                 do k=1,l_data
                    do j=1,jnp
                       q2(j,k) = qt(j,k)
                    enddo
                 enddo

             if( m_time .le. d_time(m) ) then
                 m2 = m
                 write(6,*) 'm2 = ', m2
                 go to 123
             endif
          enddo

! Model time is late Dec
          write(6,*) 'Dec and Jan H2O data will be used'
          m2 = 1
                 do k=1,l_data
                    do j=1,jnp
                       q1(j,k) = q2(j,k)
                    enddo
                 enddo

          read(q_unit, rec=1) qt

                 do k=1,l_data
                    do j=1,jnp
                       q2(j,k) = qt(j,k)
                    enddo
                 enddo

         m = 1
123      write(6,*) ' '
      endif

! Original data is in grads format with k=1 the lowest level

      call zflip(q1, jnp, l_data, 1)
      call zflip(q2, jnp, l_data, 1)

! Interpolate to model levels
      do k=1,kw_bot
         lnp_m(k) = log(press(k))
      enddo

      do k=1,l_data
         lnp_d(k) = log(p_data(k))
      enddo

         call intpz(q1, qm1, jnp, l_data, lnp_d, lnp_m,
     &              kw_top, kw_bot, nl)
         call intpz(q2, qm2, jnp, l_data, lnp_d, lnp_m,
     &              kw_top, kw_bot, nl)

      first = .false.
      endif
!----- End init

      if(m2 .eq. 1 .and. m_time .gt. d_time(12))  go to 3333
      if(m_time .gt. d_time(m2)) then

          if(m2 .eq. 12) then
              m2 = 1
          else
              m2 = m2+1
          endif

          write(6,*) 'Reading ',m2,'th month of H2O data'
          read(q_unit, rec=m2) qt

                 do k=1,l_data
                    do j=1,jnp
                       q2(j,k) = qt(j,k)
                    enddo
                 enddo

!        qmax = vmax(q2, qmin, jnp*l_data)
!        write(6,*) m2,' qmax=',qmax,' qmin=',qmin

          call zflip(q2, jnp, l_data, 1)

          do j=1, jnp
            do k=1,kw_bot
               qm1(k,j) = qm2(k,j)
            enddo
          enddo

      do k=1,kw_bot
         lnp_m(k) = log(press(k))
      enddo

      do k=1,l_data
         lnp_d(k) = log(p_data(k))
      enddo

      call intpz(q2, qm2, jnp, l_data, lnp_d, lnp_m,
     &           kw_top, kw_bot, nl)

      endif

3333  continue
      if( m2 .eq. 1) Then
              res = 366. - d_time(12)
          if(m_time .le. d_time(12) )then
!            early Jan
              fac = (m_time + res) / (d_time(1) + res)
          else
!            
              fac = (m_time - d_time(12)) / (d_time(1) + res)
          endif
      else
          fac = (m_time - d_time(m2-1)) / (d_time(m2) - d_time(m2-1))
      endif

! Interpolate H2O data in time and
! convert to mass mixing ratio.

      pi = 4.*atan(1.)
      fac = 0.5 * ( 1. - cos(fac*pi) )
      endif                  ! End master id section

# if defined ( SPMD )
! broadcast fac, qm1, and qm2
      call mp_barrier
      call mp_bcst_real(fac)
      call mp_bcst_r2d(nl, jnp, jfirst, jlast, qm1, 0)
      call mp_bcst_r2d(nl, jnp, jfirst, jlast, qm2, 0)
#endif

!$omp parallel do private(j,k)
      do j=jfirst, jlast
         do k=1, kw_bot
            q_ref(k,j) = vol2m*(qm1(k,j) + fac*(qm2(k,j)-qm1(k,j)))
         enddo
      enddo

      return
      end
