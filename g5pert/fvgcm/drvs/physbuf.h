! NCAR CCM3 Physics Buffer

      real fsns(imr,jnp)      ! surface solar absorbed flux
      real precst(imr,jnp)    ! total snow precip. rate
      common /physbf3/ fsns, precst
 
! surface flag
      real icemask(imr,jnp)
      common/icemask/icemask

