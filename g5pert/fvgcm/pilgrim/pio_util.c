#if (defined SPMD) && (defined PIO)
#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <errno.h>
#include <signal.h>


#define CACHE_LINE_SIZE 128
#define SECTION_ROUND (1024*1024)
#define BARRIER_AREA_ROUND (16*1024)
#define round_up(n, round)  (((n + (round - 1)) / round) * round)


#define ROUND (1024*1024)
#define maxpe  200

int pid_of_doutd;     /* I/O's pid for refout  */
int pid_of_goutd;     /* I/O's pid for diagout */


char *mmap_addr_4gout;       //each PE has its own copy
unsigned long size_of_pe;    //each PE has its own copy

int fd4gout[maxpe];        //# of PEs


/**********************************************************************/
void u_exit(){
  exit(0);
}

/***********************************************************/
void getfd_4gout_c_(int *numpro) {
int inumpro = *numpro;
int i;
char buf[100];

for (i=0; i<inumpro; ++i) {
    sprintf(buf, "./.diagtmp.%d", i);
    unlink(buf);                  
    fd4gout[i] = open(buf, O_RDWR|O_CREAT|O_NONBLOCK, 0600);
    if (fd4gout[i] < 0) {
       perror("open of the shm file in getfd_c");
       exit(1);
    }

}

 return;

}

/**********************************************************/ 
/* please call allocate_shm first, and then call deallocate_shm.
It is not good for each PE to call allocate_shm twice before call deallocate  */
void deallocate_shm_4gout_c_()
{
   munmap(mmap_addr_4gout, size_of_pe);

}
/**********************************************************/ 
void allocate_shm_4gout_c_ (long *pointer, long *size, int *gid, int *tr) {

  char *mmap_addr;
  char buf[100];
  long isize = *size;
  int  igid  = *gid;
  int  itr   = *tr ;
  int io=1;
  int fd;
  struct stat statbuf;

  sprintf(buf, "./.diagtmp.%d", igid);

//  printf("buf %s itr %d\n", buf,  itr);


  if ( isize == 0 ) {          
       io=0;
       fd = open(buf, O_RDWR|O_CREAT|O_NONBLOCK, 0600);
  } else {                          //read and write
    if (igid >= 0 ) {               // called by each PE
      fd = dup(fd4gout[igid]);
    } else {                        // called from pio_init, 
                                    // for "peinfo" and other shared variables
                                    // or a read-only file
//      fd4gout[maxpe+igid] = open(buf, O_RDWR|O_CREAT|O_NONBLOCK, 0600);
      fd = open(buf, O_RDWR|O_CREAT|O_NONBLOCK, 0600);
      if (fd < 0) {
          perror("open of the shm file in allocate_shm");
          exit(1);
      }
//      fd=dup(fd4gout[maxpe+igid]);
    }
  }

//  fd1 = open("./info", O_RDWR|O_CREAT|O_NONBLOCK, 0600);

  if (isize == 0 ) {                //called by the PIO daemon
     if (fstat(fd, &statbuf) < 0) {
         perror("fstat error");
     }
     isize= (long) statbuf.st_size;           //mask?
//     sprintf(buf, "%d", isize);
//     write(fd1, buf, strlen(buf));
  }



  if ( io !=0 && itr == 1 ) {                //for a new r/w file
     /* Grow the file */
     if (ftruncate(fd, (off_t)isize) < 0) {
         perror("ftruncate(2) of the shm file");
         exit(1);
     }
//     if (igid < 0 ){
//        unlink(buf);    //4debug,           
//     }
  }


//  printf("call mmap \n");
  /* Map the file "shared" */
  mmap_addr = (char *)mmap(0, (size_t) isize,
                        PROT_READ|PROT_WRITE, MAP_SHARED, fd, 0);
  if (mmap_addr == (char *)MAP_FAILED) {
    perror("mmap(2) of the shm file");
    exit(1);
  }

  *pointer = (long) mmap_addr;

  mmap_addr_4gout = mmap_addr;     //used in deallocate
  size_of_pe      = isize    ;     //used in deallocate




     if (igid < 0 ){
//        close(fd);         //important when fd is used twice
     }


 return;
}



/**********************************************************************/

void do_io_4gout (){

  collect_io_4gout_f_();           //call a Fortran subroutine

  return;

}


/*****************************/
void pio_handler_4gout(){

//printf("waiting for SIGUSR1 in pio_handler_4gout \n");
while (1) {
//  signal(SIGUSR2, do_io_4gout);   /* must be inside the loop! */
  signal(SIGUSR1, do_io_4gout);   /* must be inside the loop! */
  usleep(10);   /*keep it! */

}


}



/**********************************************************************/
void spawn_iod_4gout_c_(){
int pid;
int retid;
// int status;
int pfd[2];
// long prtmp;

signal(SIGTERM, u_exit);
signal(SIGINT, u_exit);
signal(SIGQUIT, u_exit);

  if (pipe(pfd) == -1) {
     perror("pipe");
     exit(1);
  }

if ((pid=fork())<0)
     exit(-1);
if (pid == 0){  /* child */

  switch(pid=fork())
  {

   case -1:
      perror("fork");
      break;
   case 0:
      pio_handler_4gout();
      break;
   default:
      write(pfd[1], (char *)&pid, sizeof(int));
      exit(0);     /* return grandchild's pid, not good */
      break;
  }


 }
retid=waitpid(pid, NULL, 0);
read(pfd[0], (char *) &pid_of_goutd, sizeof( int));         /*? length? */

// printf("pid_of_goutd %d \n", pid_of_goutd);
if ( retid != pid)
 {
   perror("waitpid");
   exit(3);
 }

}

/**********************************************************************/
void awake_goutd_c_() {
kill (pid_of_goutd, SIGUSR1 );
return;

}
#endif
