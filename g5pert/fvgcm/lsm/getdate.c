/*
 $Id$
 $Author$
*/
#include <errno.h>
#include <stdio.h>
#include <time.h>
/* 
*** FORTRAN-callable routine to return current date and time information on
*** various UNIX platforms
*/

#undef FORTRANCAPS
#undef FORTRANUNDERSCORE
#if (defined( SUN ) || defined( __sgi ) || defined( DEC ) || defined( linux ))
#define FORTRANUNDERSCORE
#elif (defined( CRAY ) || defined( CRAY_T3E ))
#define FORTRANCAPS
#endif

#if ( defined FORTRANCAPS )
#define get_date_ GET_DATE
#elif ( ! defined FORTRANUNDERSCORE )
#define get_date_ get_date
#endif

void get_date_(int *mon, int *day, int *yr, int *hr, int *min, int *sec)
{
  extern int sys_nerr;        /* # of system error codes */

#ifdef linux
  extern const char *const sys_errlist[]; /* english translation of error codes */
#else
  extern char *sys_errlist[]; /* english translation of error codes */
#endif

  extern int errno;           /* error code # */

  time_t clock;
  struct tm *tmptr;

  if ( ( clock = time( 0 ) ) == -1 ) {
    fprintf( stderr, "time  failed: %d, %s\n",
	    errno, sys_errlist[errno] );
    exit( -1  );
  }
  tmptr = localtime(&clock);
  *mon = tmptr->tm_mon + 1; /* Change 0-11 range to 1-12 */
  *day = tmptr->tm_mday;
  *yr  = tmptr->tm_year;
  *hr = tmptr->tm_hour;
  *min = tmptr->tm_min;
  *sec = tmptr->tm_sec;
  return;
}
