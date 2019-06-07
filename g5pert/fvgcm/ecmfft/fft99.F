c This is the FFT package from ECMWF
c set99 has been renamed FFTFAX and FFT991 renamed as RFFTMLT
c
      subroutine fft99(a,work,trigs,ifax,inc,jump,n,lot,isign)
c
c PURPOSE      PERFORMS MULTIPLE FAST FOURIER TRANSFORMS.  THIS PACKAGE
c              WILL PERFORM A NUMBER OF SIMULTANEOUS REAL/HALF-COMPLEX
c              PERIODIC FOURIER TRANSFORMS OR CORRESPONDING INVERSE
c              TRANSFORMS, I.E.  GIVEN A SET OF REAL DATA VECTORS, THE
c              PACKAGE RETURNS A SET OF 'HALF-COMPLEX' FOURIER
c              COEFFICIENT VECTORS, OR VICE VERSA.  THE LENGTH OF THE
c              TRANSFORMS MUST BE AN EVEN NUMBER GREATER THAN 4 THAT HAS
c              NO OTHER FACTORS EXCEPT POSSIBLY POWERS OF 2, 3, AND 5.
c              THIS IS AN ALL FORTRAN VERSION OF THE CRAYLIB PACKAGE
c              THAT IS MOSTLY WRITTEN IN CAL.
c
c              THE PACKAGE FFT99F CONTAINS SEVERAL USER-LEVEL ROUTINES:
c
c            SUBROUTINE SET99
c                AN INITIALIZATION ROUTINE THAT MUST BE CALLED ONCE
c                BEFORE A SEQUENCE OF CALLS TO THE FFT ROUTINES
c                (PROVIDED THAT N IS NOT CHANGED).
c
c            SUBROUTINES FFT99 AND FFT991
c                TWO FFT ROUTINES THAT RETURN SLIGHTLY DIFFERENT
c                ARRANGEMENTS OF THE DATA IN GRIDPOINT SPACE.
c
c
c ACCESS       THIS FORTRAN VERSION MAY BE ACCESSED WITH
c
c                   *FORTRAN,P=XLIB,SN=FFT99F
c
c              TO ACCESS THE CRAY OBJECT CODE, CALLING THE USER ENTRY
c              POINTS FROM A CRAY PROGRAM IS SUFFICIENT.  THE SOURCE
c              FORTRAN AND CAL CODE FOR THE CRAYLIB VERSION MAY BE
c              ACCESSED USING
c
c                   FETCH P=CRAYLIB,SN=FFT99
c                   FETCH P=CRAYLIB,SN=CAL99
c
c USAGE        LET N BE OF THE FORM 2**P * 3**Q * 5**R, WHERE P .GE. 1,
c              Q .GE. 0, AND R .GE. 0.  THEN A TYPICAL SEQUENCE OF
c              CALLS TO TRANSFORM A GIVEN SET OF REAL VECTORS OF LENGTH
c              N TO A SET OF 'HALF-COMPLEX' FOURIER COEFFICIENT VECTORS
c              OF LENGTH N IS
c
c                   DIMENSION IFAX(13),TRIGS(3*N/2+1),A(M*(N+2)),
c                  +          WORK(M*(N+1))
c
c                   CALL SET99 (TRIGS, IFAX, N)
c                   CALL FFT99 (A,WORK,TRIGS,IFAX,INC,JUMP,N,M,ISIGN)
c
c              SEE THE INDIVIDUAL WRITE-UPS FOR SET99, FFT99, AND
c              FFT991 BELOW, FOR A DETAILED DESCRIPTION OF THE
c              ARGUMENTS.
c
c HISTORY      THE PACKAGE WAS WRITTEN BY CLIVE TEMPERTON AT ECMWF IN
c              NOVEMBER, 1978.  IT WAS MODIFIED, DOCUMENTED, AND TESTED
c              FOR NCAR BY RUSS REW IN SEPTEMBER, 1980.
c
c-----------------------------------------------------------------------
c
c SUBROUTINE SET99 (TRIGS, IFAX, N)
c
c PURPOSE      A SET-UP ROUTINE FOR FFT99 AND FFT991.  IT NEED ONLY BE
c              CALLED ONCE BEFORE A SEQUENCE OF CALLS TO THE FFT
c              ROUTINES (PROVIDED THAT N IS NOT CHANGED).
c
c ARGUMENT     IFAX(13),TRIGS(3*N/2+1)
c DIMENSIONS
c
c ARGUMENTS
c
c ON INPUT     TRIGS
c               A FLOATING POINT ARRAY OF DIMENSION 3*N/2 IF N/2 IS
c               EVEN, OR 3*N/2+1 IF N/2 IS ODD.
c
c              IFAX
c               AN INTEGER ARRAY.  THE NUMBER OF ELEMENTS ACTUALLY USED
c               WILL DEPEND ON THE FACTORIZATION OF N.  DIMENSIONING
c               IFAX FOR 13 SUFFICES FOR ALL N LESS THAN A MILLION.
c
c              N
c               AN EVEN NUMBER GREATER THAN 4 THAT HAS NO PRIME FACTOR
c               GREATER THAN 5.  N IS THE LENGTH OF THE TRANSFORMS (SEE
c               THE DOCUMENTATION FOR FFT99 AND FFT991 FOR THE
c               DEFINITIONS OF THE TRANSFORMS).
c
c ON OUTPUT    IFAX
c               CONTAINS THE FACTORIZATION OF N/2.  IFAX(1) IS THE
c               NUMBER OF FACTORS, AND THE FACTORS THEMSELVES ARE STORED
c               IN IFAX(2),IFAX(3),...  IF SET99 IS CALLED WITH N ODD,
c               OR IF N HAS ANY PRIME FACTORS GREATER THAN 5, IFAX(1)
c               IS SET TO -99.
c
c              TRIGS
c               AN ARRAY OF TRIGONOMETRIC FUNCTION VALUES SUBSEQUENTLY
c               USED BY THE FFT ROUTINES.
c
c-----------------------------------------------------------------------
c
c SUBROUTINE FFT991 (A,WORK,TRIGS,IFAX,INC,JUMP,N,M,ISIGN)
c                       AND
c SUBROUTINE FFT99 (A,WORK,TRIGS,IFAX,INC,JUMP,N,M,ISIGN)
c
c PURPOSE      PERFORM A NUMBER OF SIMULTANEOUS REAL/HALF-COMPLEX
c              PERIODIC FOURIER TRANSFORMS OR CORRESPONDING INVERSE
c              TRANSFORMS, USING ORDINARY SPATIAL ORDER OF GRIDPOINT
c              VALUES (FFT991) OR EXPLICIT CYCLIC CONTINUITY IN THE
c              GRIDPOINT VALUES (FFT99).  GIVEN A SET
c              OF REAL DATA VECTORS, THE PACKAGE RETURNS A SET OF
c              'HALF-COMPLEX' FOURIER COEFFICIENT VECTORS, OR VICE
c              VERSA.  THE LENGTH OF THE TRANSFORMS MUST BE AN EVEN
c              NUMBER THAT HAS NO OTHER FACTORS EXCEPT POSSIBLY POWERS
c              OF 2, 3, AND 5.  THESE VERSION OF FFT991 AND FFT99 ARE
c              OPTIMIZED FOR USE ON THE CRAY-1.
c
c ARGUMENT     A(M*(N+2)), WORK(M*(N+1)), TRIGS(3*N/2+1), IFAX(13)
c DIMENSIONS
c
c ARGUMENTS
c
c ON INPUT     A
c               AN ARRAY OF LENGTH M*(N+2) CONTAINING THE INPUT DATA
c               OR COEFFICIENT VECTORS.  THIS ARRAY IS OVERWRITTEN BY
c               THE RESULTS.
c
c              WORK
c               A WORK ARRAY OF DIMENSION M*(N+1)
c
c              TRIGS
c               AN ARRAY SET UP BY SET99, WHICH MUST BE CALLED FIRST.
c
c              IFAX
c               AN ARRAY SET UP BY SET99, WHICH MUST BE CALLED FIRST.
c
c              INC
c               THE INCREMENT (IN WORDS) BETWEEN SUCCESSIVE ELEMENTS OF
c               EACH DATA OR COEFFICIENT VECTOR (E.G.  INC=1 FOR
c               CONSECUTIVELY STORED DATA).
c
c              JUMP
c               THE INCREMENT (IN WORDS) BETWEEN THE FIRST ELEMENTS OF
c               SUCCESSIVE DATA OR COEFFICIENT VECTORS.  ON THE CRAY-1,
c               TRY TO ARRANGE DATA SO THAT JUMP IS NOT A MULTIPLE OF 8
c               (TO AVOID MEMORY BANK CONFLICTS).  FOR CLARIFICATION OF
c               INC AND JUMP, SEE THE EXAMPLES BELOW.
c
c              N
c               THE LENGTH OF EACH TRANSFORM (SEE DEFINITION OF
c               TRANSFORMS, BELOW).
c
c              M
c               THE NUMBER OF TRANSFORMS TO BE DONE SIMULTANEOUSLY.
c
c              ISIGN
c               = +1 FOR A TRANSFORM FROM FOURIER COEFFICIENTS TO
c                    GRIDPOINT VALUES.
c               = -1 FOR A TRANSFORM FROM GRIDPOINT VALUES TO FOURIER
c                    COEFFICIENTS.
c
c ON OUTPUT    A
c               IF ISIGN = +1, AND M COEFFICIENT VECTORS ARE SUPPLIED
c               EACH CONTAINING THE SEQUENCE:
c
c               A(0),B(0),A(1),B(1),...,A(N/2),B(N/2)  (N+2 VALUES)
c
c               THEN THE RESULT CONSISTS OF M DATA VECTORS EACH
c               CONTAINING THE CORRESPONDING N+2 GRIDPOINT VALUES:
c
c               FOR FFT991, X(0), X(1), X(2),...,X(N-1),0,0.
c               FOR FFT99, X(N-1),X(0),X(1),X(2),...,X(N-1),X(0).
c                   (EXPLICIT CYCLIC CONTINUITY)
c
c               WHEN ISIGN = +1, THE TRANSFORM IS DEFINED BY:
c                 X(J)=SUM(K=0,...,N-1)(C(K)*EXP(2*I*J*K*PI/N))
c                 WHERE C(K)=A(K)+I*B(K) AND C(N-K)=A(K)-I*B(K)
c                 AND I=SQRT (-1)
c
c               IF ISIGN = -1, AND M DATA VECTORS ARE SUPPLIED EACH
c               CONTAINING A SEQUENCE OF GRIDPOINT VALUES X(J) AS
c               DEFINED ABOVE, THEN THE RESULT CONSISTS OF M VECTORS
c               EACH CONTAINING THE CORRESPONDING FOURIER COFFICIENTS
c               A(K), B(K), 0 .LE. K .LE N/2.
c
c               WHEN ISIGN = -1, THE INVERSE TRANSFORM IS DEFINED BY:
c                 C(K)=(1/N)*SUM(J=0,...,N-1)(X(J)*EXP(-2*I*J*K*PI/N))
c                 WHERE C(K)=A(K)+I*B(K) AND I=SQRT(-1)
c
c               A CALL WITH ISIGN=+1 FOLLOWED BY A CALL WITH ISIGN=-1
c               (OR VICE VERSA) RETURNS THE ORIGINAL DATA.
c
c               NOTE: THE FACT THAT THE GRIDPOINT VALUES X(J) ARE REAL
c               IMPLIES THAT B(0)=B(N/2)=0.  FOR A CALL WITH ISIGN=+1,
c               IT IS NOT ACTUALLY NECESSARY TO SUPPLY THESE ZEROS.
c
c EXAMPLES      GIVEN 19 DATA VECTORS EACH OF LENGTH 64 (+2 FOR EXPLICIT
c               CYCLIC CONTINUITY), COMPUTE THE CORRESPONDING VECTORS OF
c               FOURIER COEFFICIENTS.  THE DATA MAY, FOR EXAMPLE, BE
c               ARRANGED LIKE THIS:
c
c FIRST DATA   A(1)=    . . .                A(66)=             A(70)
c VECTOR       X(63) X(0) X(1) X(2) ... X(63) X(0)  (4 EMPTY LOCATIONS)
c
c SECOND DATA  A(71)=   . . .                                  A(140)
c VECTOR       X(63) X(0) X(1) X(2) ... X(63) X(0)  (4 EMPTY LOCATIONS)
c
c               AND SO ON.  HERE INC=1, JUMP=70, N=64, M=19, ISIGN=-1,
c               AND FFT99 SHOULD BE USED (BECAUSE OF THE EXPLICIT CYCLIC
c               CONTINUITY).
c
c               ALTERNATIVELY THE DATA MAY BE ARRANGED LIKE THIS:
c
c                FIRST         SECOND                          LAST
c                DATA          DATA                            DATA
c                VECTOR        VECTOR                          VECTOR
c
c                 A(1)=         A(2)=                           A(19)=
c
c                 X(63)         X(63)       . . .               X(63)
c        A(20)=   X(0)          X(0)        . . .               X(0)
c        A(39)=   X(1)          X(1)        . . .               X(1)
c                  .             .                               .
c                  .             .                               .
c                  .             .                               .
c
c               IN WHICH CASE WE HAVE INC=19, JUMP=1, AND THE REMAINING
c               PARAMETERS ARE THE SAME AS BEFORE.  IN EITHER CASE, EACH
c               COEFFICIENT VECTOR OVERWRITES THE CORRESPONDING INPUT
c               DATA VECTOR.
c
c-----------------------------------------------------------------------
      dimension a(jump*lot),work((n+1)*lot),trigs(3*n/2+1),ifax(13)
c
c     SUBROUTINE "FFT99" - MULTIPLE FAST REAL PERIODIC TRANSFORM
c     CORRESPONDING TO OLD SCALAR ROUTINE FFT9
c     PROCEDURE USED TO CONVERT TO HALF-LENGTH COMPLEX TRANSFORM
c     IS GIVEN BY COOLEY, LEWIS AND WELCH (J. SOUND VIB., VOL. 12
c     (1970), 315-337)
c
c     A IS THE ARRAY CONTAINING INPUT AND OUTPUT DATA
c     WORK IS AN AREA OF SIZE (N+1)*LOT
c     TRIGS IS A PREVIOUSLY PREPARED LIST OF TRIG FUNCTION VALUES
c     IFAX IS A PREVIOUSLY PREPARED LIST OF FACTORS OF N/2
c     INC IS THE INCREMENT WITHIN EACH DATA 'VECTOR'
c         (E.G. INC=1 FOR CONSECUTIVELY STORED DATA)
c     JUMP IS THE INCREMENT BETWEEN THE START OF EACH DATA VECTOR
c     N IS THE LENGTH OF THE DATA VECTORS
c     LOT IS THE NUMBER OF DATA VECTORS
c     ISIGN = +1 FOR TRANSFORM FROM SPECTRAL TO GRIDPOINT
c           = -1 FOR TRANSFORM FROM GRIDPOINT TO SPECTRAL
c
c     ORDERING OF COEFFICIENTS:
c         A(0),B(0),A(1),B(1),A(2),B(2),...,A(N/2),B(N/2)
c         WHERE B(0)=B(N/2)=0; (N+2) LOCATIONS REQUIRED
c
c     ORDERING OF DATA:
c         X(N-1),X(0),X(1),X(2),...,X(N),X(0)
c         I.E. EXPLICIT CYCLIC CONTINUITY; (N+2) LOCATIONS REQUIRED
c
c     VECTORIZATION IS ACHIEVED ON CRAY BY DOING THE TRANSFORMS IN
c     PARALLEL
c
c     *** N.B. N IS ASSUMED TO BE AN EVEN NUMBER
c
c     DEFINITION OF TRANSFORMS:
c     -------------------------
c
c     ISIGN=+1: X(J)=SUM(K=0,...,N-1)(C(K)*EXP(2*I*J*K*PI/N))
c         WHERE C(K)=A(K)+I*B(K) AND C(N-K)=A(K)-I*B(K)
c
c     ISIGN=-1: A(K)=(1/N)*SUM(J=0,...,N-1)(X(J)*COS(2*J*K*PI/N))
c               B(K)=-(1/N)*SUM(J=0,...,N-1)(X(J)*SIN(2*J*K*PI/N))
c
c
c
c
      nfax=ifax(1)
      nx=n+1
      nh=n/2
      ink=inc+inc
      if (isign.eq.+1) go to 30
c
c     IF NECESSARY, TRANSFER DATA TO WORK AREA
      igo=50
      if (mod(nfax,2).eq.1) goto 40
      ibase=inc+1
      jbase=1
      do 20 l=1,lot
      i=ibase
      j=jbase
cDIR$ IVDEP
      do 10 m=1,n
      work(j)=a(i)
      i=i+inc
      j=j+1
   10 continue
      ibase=ibase+jump
      jbase=jbase+nx
   20 continue
c
      igo=60
      go to 40
c
c     PREPROCESSING (ISIGN=+1)
c     ------------------------
c
   30 continue
      call fft99a(a,work,trigs,inc,jump,n,lot)
      igo=60
c
c     COMPLEX TRANSFORM
c     -----------------
c
   40 continue
      ia=inc+1
      la=1
      do 80 k=1,nfax
      if (igo.eq.60) go to 60
   50 continue
      call vpassm(a(ia),a(ia+inc),work(1),work(2),trigs,
     *   ink,2,jump,nx,lot,nh,ifax(k+1),la)
      igo=60
      go to 70
   60 continue
      call vpassm(work(1),work(2),a(ia),a(ia+inc),trigs,
     *    2,ink,nx,jump,lot,nh,ifax(k+1),la)
      igo=50
   70 continue
      la=la*ifax(k+1)
   80 continue
c
      if (isign.eq.-1) go to 130
c
c     IF NECESSARY, TRANSFER DATA FROM WORK AREA
      if (mod(nfax,2).eq.1) go to 110
      ibase=1
      jbase=ia
      do 100 l=1,lot
      i=ibase
      j=jbase
cDIR$ IVDEP
      do 90 m=1,n
      a(j)=work(i)
      i=i+1
      j=j+inc
   90 continue
      ibase=ibase+nx
      jbase=jbase+jump
  100 continue
c
c     FILL IN CYCLIC BOUNDARY POINTS
  110 continue
      ia=1
      ib=n*inc+1
cDIR$ IVDEP
      do 120 l=1,lot
      a(ia)=a(ib)
      a(ib+inc)=a(ia+inc)
      ia=ia+jump
      ib=ib+jump
  120 continue
      go to 140
c
c     POSTPROCESSING (ISIGN=-1):
c     --------------------------
c
  130 continue
      call fft99b(work,a,trigs,inc,jump,n,lot)
c
  140 continue
      return
      end
