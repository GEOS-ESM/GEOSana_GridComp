c
c $Id$
c $Author$
c
C
C Subgrid scale stationary orographic gravity wave drag isotropic
C standard deviations.  Computed from a global, 10' resolution
C topographic height dataset from the US Navy.  Standard deviations
C about the mean topographic height were computed on a spatial scale
C dependent on the model resolution; for T42 and lower resolutions, 2.0
C X 2.0 degrees was used; for T63, 1.67 X 1.67 degrees; and for T106 1.0
C X 1.0 degrees.  The isotropic standard deviations were averaged
C together for each model grid box, then smoothed twice by a 1-2-1
C spatial filter in both latitude and longitude.  Values over the oceans
C are zero.
C
      real sgh     ! Standard deviations for gravity wave drag
C
      common /comsgh/ sgh(plond,plat)
C
 
