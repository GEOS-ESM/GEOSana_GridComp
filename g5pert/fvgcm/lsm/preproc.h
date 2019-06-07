#include <fvgcm.h>
#ifndef PREPROC_SET
#define PREPROC_SET
#define LSMLAT FVGCM_LAT
#define LSMLON FVGCM_LON
#define NCPREC NCDOUBLE
#define COUP_CCM
#define SHELL_MSS
#if defined( CRAY_T3E )
#define ISIZE INTEGER8
#else
#define ISIZE INTEGER4
#endif
#if ( defined RS6K )
#define LOGLEN 4
#else
#define LOGLEN 8
#endif
#endif
 
