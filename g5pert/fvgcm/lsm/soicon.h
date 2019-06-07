* ------------------------ code history ---------------------------
* source file:       soicon.h
* purpose:           "soil type" constants
* date last revised: March 1996 - lsm version 1
* author:            Gordon Bonan
* standardized:      J. Truesdale, Feb 1996
* reviewed:          G. Bonan, Feb 1996
* -----------------------------------------------------------------

      common /soicon_i/ istice,istdlak,istslak,istwet

      integer istice     !land ice     "soil" type
      integer istdlak    !deep lake    "soil" type
      integer istslak    !shallow lake "soil" type
      integer istwet     !wetland      "soil" type

      common /soicon_r/ eg,rlsoi

      real eg(mst)       !ground emissivity: varies with "soil" type
      real rlsoi(mst)    !ground roughness length (m): varies with "soil" type 

* ------------------------ end soicon.h ---------------------------
 
