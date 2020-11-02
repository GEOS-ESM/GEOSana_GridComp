#!/usr/bin/env perl
#--------------------------------------------------
#
# Purpose: create info files at given date/time
#          according to information in database. 
#
# Usage:
#
#  gsiinfo.pl [options] EXPID NYMD NHMS
#
# !REVISION HISTORY:
#
#   14Apr2012 Todling  Initial code
#                      (based on Jing's satinfo analyzer code)
#
#--------------------------------------------------
use Env;                 # make env vars readily available
use File::Basename;      # for basename(), dirname()
use File::Copy "cp";     # for cp()
use Getopt::Long;        # load module with GetOptions function
use Time::Local;         # time functions
use FindBin;             # so we can find where this script resides

# look for perl packages in the following locations
#--------------------------------------------------
use lib ( "$FindBin::Bin", "$FVROOT/bin", "$ESMADIR/$ARCH/bin" );

  GetOptions ( "log",
               "debug",
               "h" );

# Initialize variables
# --------------------
  init();

# Run program for each case
# -------------------------
# info_run("conv"); # convinfo: do link below until this is mature
  Assignfn( "$fvwork/gmao_global_convinfo.rc", "convinfo");
  info_run("sat");
  info_run("oz");
 
  exit(0);

#..................................................................

sub info_run {

   my ( $this ) = @_;
   my ( $iuse_flags );

   if ( $this eq "conv" ) {
       $exec_def = "make_convinfo.x";
       $db_def   = "gmao_convinfo.db";
       $tmpl_def = "gmao_global_convinfo.rc";
       $gsiname  = "convinfo";
       $iuse_flags = "notused_flag = -2";
   }
   if ( $this eq "oz" ) {
       $exec_def = "make_ozinfo.x";
       $db_def   = "gmao_ozinfo.db";
       $tmpl_def = "gmao_global_ozinfo.rc";
       $gsiname  = "ozinfo";
       $iuse_flags = "notused_flag = -2";
   }
   if ( $this eq "sat" ) {
       $exec_def = "make_satinfo.x";
       $db_def   = "gmao_satinfo.db";
       $tmpl_def = "gmao_global_satinfo.rc";
       $gsiname  = "satinfo";
       $iuse_flags = "notused_flag = -2";
   }

#
#  Create Satellite info ( ozinfo | gmao_global_ozinfo.yyyymmdd_hhz.txt );
#  /-----------------------------------------------------------------------
   # Setup make_gsiinfo.x() through a simple NAMELIST, such that
   # all configurable parameters are explicit to the developers
   # of analyzer().

   $mksi_sidb = "$FVROOT/etc/$db_def";
     if ( ! -r $mksi_sidb ) { $mksi_sidb = "$FVROOT/etc/$db_def" };
     if ( $this eq "conv" ) {
        if ( $ENV{MKSICN_SIDB} ) { $mksi_sidb = $ENV{MKSICN_SIDB} };
     }
     if ( $this eq "oz" ) {
        if ( $ENV{MKSIOZ_SIDB} ) { $mksi_sidb = $ENV{MKSIOZ_SIDB} };
     }
     if ( $this eq "sat" ) {
        if ( $ENV{MKSI_SIDB} ) { $mksi_sidb = $ENV{MKSI_SIDB} };
     }

   $mksi_tmpl = "$fvwork/$tmpl_def";
     if ( ! -r $mksi_tmpl ) { $mksi_tmpl = "$FVROOT/etc/$tmpl_def" };
     if ( $this eq "conv" ) {
        if ( $ENV{MKSICN_TMPL} ) { $mksi_tmpl = $ENV{MKSICN_TMPL} };
     }
     if ( $this eq "oz" ) {
        if ( $ENV{MKSIOZ_TMPL} ) { $mksi_tmpl = $ENV{MKSIOZ_TMPL} };
     }
     if ( $this eq "sat" ) {
        if ( $ENV{MKSI_TMPL} ) { $mksi_tmpl = $ENV{MKSI_TMPL} };
     }

   if ( $this eq "sat" ) { # Jing: can we reconcile these?
      $cmd = "echo '&setup
	   nymd = $nymd
	   nhms = $nhms
	   $iuse_flags
	   satinfo_tmpl = \"$mksi_tmpl\"
	   satinfo_outf = \"./$gsiname\"
	   dbname    = \"$mksi_sidb\"
	   nowarn = .true.
	   /' | $FVROOT/bin/$exec_def " ;
   } else {
      $cmd = "echo '&setup
	   nymd = $nymd
	   nhms = $nhms
	   $iuse_flags
	   info_tmpl = \"$mksi_tmpl\"
	   info_outf = \"./$gsiname\"
	   dbname    = \"$mksi_sidb\"
	   nowarn = .true.
	   /' | $FVROOT/bin/$exec_def " ;
   }

   print " $cmd\n";
#  $rc = System($cmd, "$log_ana","$exec_def");
   $rc = system($cmd);
   die ">>>> ERROR <<< running $exec_def" if ( $rc );
   print " $0: $exec_def \$rc =  $rc\n";

#  make a timely archive of the cached "ozinfo" file.
   cp( "./$gsiname","$fvwork/$expid.gmao_global_${gsiname}.${nymd}_${hh}z.txt" );

#  end-of-make_ozinfo()
#  \------------------------------------------------------------------------

}
#....................................................................................
sub init {

   if ( $#ARGV  <  2 || $opt_h ) {
     print STDERR " Improper input parameters ; see usage:\n";
     usage();
   } else {              # required command line args
     $expid  = $ARGV[0];
     $nymd   = $ARGV[1];
     $nhms   = sprintf("%6.6d",$ARGV[2]);
   }
   $yyyy = substr($nymd,0,4);
   $mm   = substr($nymd,4,2);
   $dd   = substr($nymd,6,2);
   $hh   = substr($nhms,0,2);

   
   $fvwork = $ENV{FVWORK};
   if ( ! $fvwork ) {
       print "env var FVWORK must be defined. Aborting ...\n";
       exit(1) unless ($opt_debug);
   }

   $log = "$expid.gsiinfo.log.${nymd}_${hh}z.txt";
   if ($opt_log) {
       $log = $opt_log;
   }

}

#....................................................................................
sub System {

    my ( $cmd, $logfile, $xname ) = @_;
    my ( @zname );

    open SAVEOUT, ">&STDOUT";  # save stdout
    open SAVEERR, ">&STDERR";  # save stderr

    open STDOUT, ">>$logfile" or die "can't redirect stdout";
    open STDERR, ">>$logfile" or die "can't redirect stderr";

    select STDERR; $| = 1;     # make it unbuffered
    select STDOUT; $| = 1;     # make it unbuffered

    @zname = split(" ", $cmd);
    if ( "$zname[0]" eq "mpirun" || "$zname[0]" eq "prun" ) {
      $rc1 = system( "zeit_ci.x -r $fvwork/.zeit $xname");
    } else {
      $rc1 = system( "zeit_ci.x -r $fvwork/.zeit $zname[0]");
    }

    $rc = system ( $cmd );     # run the shell command

    if ( "$zname[0]" eq "mpirun" || "$zname[0]" eq "prun" ) {
      $rc2 = system( "zeit_co.x -r $fvwork/.zeit $xname");
    } else {
      $rc2 = system( "zeit_co.x -r $fvwork/.zeit $zname[0]");
    }

    # Bitwise shift returns actual UNIX return code
    $exit_code = $rc >> 8;

    close STDOUT;
    close STDERR;

    open STDOUT, ">&SAVEOUT" ;  # restore stdout
    open STDERR, ">&SAVEERR" ;  # restore stdout

    return $exit_code;

  }

#....................................................................................
#
# Tick - advance date/time by nsecs seconds
#
sub tick {
    my ( $nymd, $nhms, $nsecs ) = @_;

    if("$nsecs" == "0" ) {
        return ($nymd, $nhms);
    }

    $yyyy1  = substr($nymd,0,4);
    $mm1    = substr($nymd,4,2);
    $dd1    = substr($nymd,6,2);

    $hh1 = 0 unless ( $hh1 = substr($nhms,0,2));
    $mn1 = 0 unless ( $mn1 = substr($nhms,2,2));
    $ss1 = 0 unless ( $ss1 = substr($nhms,4,2));
    $time1 = timegm($ss1,$mn1,$hh1,$dd1,$mm1-1,$yyyy1) + $nsecs;
    ($ss1,$mn1,$hh1,$dd1,$mm1,$yyyy1,undef,undef,undef) = gmtime($time1);

    $nymd = (1900+$yyyy1)*10000 + ($mm1+1)*100 + $dd1;
    $nhms = sprintf("%6.6d",$hh1*10000 + $mn1*100 + $ss1);
    return ($nymd, $nhms);

}

#....................................................................................
sub Assignfn {

# Assignfn - assigns fn to given file name fname.
# fname = old file
# fn = new file (links to old)
  my ( $fname, $fn ) = @_;
  unlink($fn) if ( -e $fn ) ;
  symlink("$fname","$fn");

}

#....................................................................................
sub usage {

   print <<"EOF";

NAME
     gsiinfo.pl - generates info files for GSI according to database

SYNOPSIS

     gsiinfo.pl [options] EXPID NYMD NHMS

DESCRIPTION

     This script involves creates info files for a given
     date and time using a database

     expid   experiment id
     nymd    date, as in YYYYMMDD
     nhms    time, as in HHMMSS

    Optinal Arguents:

    -h      help (echoes this usage)
    -log    name of log files
             (default: expid.sac.log.nymd_hhz.txt)
    -debug  going through it bug do not actually run it

    Optinal Env Vars:
      MKSI_SIDB    - name/location of satellite-database
      MKSI_TMPL    - name/location of satellite-template file
      MKSICN_SIDB  - name/location of conventional-obs-database
      MKSICN_TMPL  - name/location of conventional-obs-template file
      MKSIOZ_SIDB  - name/location of oz-database
      MKSIOZ_TMPL  - name/location of oz-template file

EOF

  exit(1);
}

