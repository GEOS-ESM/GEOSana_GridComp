#ifdef _DARPACK_
#define _SSCAL		DSCAL
#define _SAUPD		DSAUPD
#define _PSAUPD		PDSAUPD
#define _SEUPD		DSEUPD
#define _PSEUPD		PDSEUPD
#define _MOUT		DMOUT
#define _PMOUT		PDMOUT
#define _VOUT		DVOUT
#define _SSUM		DSUM
#define _PSNRM2		PDNORM2
#else
#define _SSCAL		SSCAL
#define _SAUPD		SSAUPD
#define _PSAUPD		PSSAUPD
#define _SEUPD		SSEUPD
#define _PSEUPD		PSSEUPD
#define _MOUT		SMOUT
#define _PMOUT		PSMOUT
#define _VOUT		SVOUT
#define _SSUM		SSUM
#define _PSNRM2		PSNORM2
#endif
