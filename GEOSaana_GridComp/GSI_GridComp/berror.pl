#!/usr/bin/env perl
#--------------------------------------------------
#
# Purpose: shuffle members of an ensemble kept for
#          the purpose of generating a dynamic
#          representation of a background error 
#          covariance for GSI.
#
# Usage:
#
#  berror.pl [options] expid nymd nhms what
#
# !TO DO:
#   1. need a way to recover from crashes
#
# !REVISION HISTORY:
#
#   23Nov2013 Todling  Initial code
#
#--------------------------------------------------
use Env;                 # make env vars readily available
use File::Basename;      # for basename(), dirname()
use File::Copy "cp";     # for cp()
use Getopt::Long;        # load module with GetOptions function
use Time::Local;         # time functions
use FindBin;             # so we can find where this script resides
use lib ( "$FindBin::Bin" );
#use POSIX qw(strftime);
use Manipulate_time;

#use strict;
#use warnings;

# Main ...
#{

my($expid,$nymd,$nhms,$freq,$debug,$what);

# FVROOT is where the binaries have been installed
# ------------------------------------------------
$fvroot  = dirname($FindBin::Bin);
$fvroot  =~ s|/u/.realmounts/share|/share|;   # for portability across
                                              # NAS machines

# Initialize variables
# --------------------
  init();

# Make ensemble available to analysis
# -----------------------------------
  if ("$what" eq "set") {
     print " Setting ensemble for background error covariance estimation ... \n";
     link_ensemble2analysis();
  }

# Update ensemble
# ---------------
  if ("$what" eq "update") {
     print " Updating ensemble used to estimate background error covariance  ... \n";
     update_ensemble();
  }

  print "berror.pl: successfully completed.\n\n";
  exit(0);
#}
#.........................................................
sub init {

   GetOptions ( "adir=s",
                "edir=s",
                "rc=s",
                "debug",
                "h" );

   if ( $#ARGV  <  4 || $opt_h ) {
     usage();
   } else {              # required command line args
     $expid = $ARGV[0];  # experiment identifier
     $nymd  = $ARGV[1];  # analysis date
     $nhms  = $ARGV[2];  # analysis time
     $freq  = $ARGV[3];  # analysis frequency
     $what  = $ARGV[4];  # set or update
   }
   $hh = `echo $nhms | cut -c1-2`; $hh = sprintf("%2.2d",$hh);
   $freq_sc = $freq * 3600; 
   $nymdhhz_now = "${nymd}_${hh}z";

   ($nymd_fro, $nhms_fro) = tick ($nymd, $nhms, $freq_sc);
   $hh_fro = `echo $nhms_fro | cut -c1-2`; $hh_fro = sprintf("%2.2d",$hh_fro);
   $nymdhhz_fro = "${nymd_fro}_${hh_fro}z";

   die "berror.pl, ERROR: need env(FVHOME) to be defined" unless($ENV{FVHOME});

   if($opt_debug) {$debug = 1};

   if($opt_adir) {
      $ANADIR = $opt_adir;
   }

   if($opt_edir) {
      $ENSDIR = $opt_edir;
   } else {
      $ENSDIR = "$FVHOME/aens4berror";
   }
   if (! -d $ENSDIR) {
      print "berror.pl: nothing to do.\n\n";
      exit(0);
   }
   
   $NCSUFFIX = "nc4"; # wired for now

   if ( $opt_rc ) {
      $fnmember_now = `echorc.x -template $expid $nymd $nhms -rc $gcrc upper-air_bkg_filename`;
      $fnmember_fro = `echorc.x -template $expid $nymd_fro $nhms_fro -rc $gcrc upper-air_bkg_filename`;
   } else {
      $fnmember_now = "$expid.bkg.eta.$nymdhhz_now.$NCSUFFIX";
      $fnmember_fro = "$expid.bkg.eta.$nymdhhz_fro.$NCSUFFIX";
   }

}
#.........................................................
sub update_ensemble {

   my($nmem,$memtag,$mem_p1,$wc);

   $wc = `ls -1d $ENSDIR/mem* | wc | cut -c1-10`; chomp($wc);  # don't quite like this
   $nmem = $wc;
   die "berror.pl, ERROR: cannot find ensemble" unless($nmem > 0);
   print "Found $nmem ensemble members\n";

   # make sure files in the ensemble directory are never looked at by the archiving procedure 
   if( ! -e "$ENSDIR/.no_archiving" ) { system("touch $ENSDIR/.no_archiving") };

   # store for a while ...
   if( -e "$ENSDIR/hold" ) {
      $rc = system("/bin/rm -r $ENSDIR/hold");
   }
   $rc = system("/bin/mkdir -p $ENSDIR/hold");
   print   "mv $ENSDIR/mem001 $ENSDIR/hold \n";
   if ( ! $debug ){
     rename("$ENSDIR/mem001","$ENSDIR/hold");
   }
   $rc = system("/bin/mkdir -p $ENSDIR/mem001");

# slide all nmem-1 members, that is, 2->1, 3->2, etc
   $nc = 0;
   while( $nc < $nmem - 1 ) {
     $nc = $nc + 1;
     $memtag = $nc    ; $memtag = sprintf("%3.3d",$memtag);
     $mem_p1 = $nc + 1; $mem_p1 = sprintf("%3.3d",$mem_p1);
     print "mv $ENSDIR/mem$mem_p1/$fnmember_now   $ENSDIR/mem$memtag/$fnmember_fro \n";
     print "reset_time.x $ENSDIR/mem$memtag/$fnmember_fro $nymd_fro $nhms_fro -9 \n";
     if ( ! $debug ) {
       # shift ensemble member  
       rename("$ENSDIR/mem$mem_p1/$fnmember_now","$ENSDIR/mem$memtag/$fnmember_fro");
       # update date/time of ensemble making it suitable to current analysis date/time
       $rc = system("reset_time.x $ENSDIR/mem$memtag/$fnmember_fro $nymd_fro $nhms_fro -9" );
       die ">>> ERROR <<< cannot reset time of member mem$memtag " if ( $rc );
     }
   }

   return unless ( $ANADIR );

# add new member as last one
   $memtag = sprintf("%3.3d",$nmem);
   print "cp $ANADIR/$fnmember_now $ENSDIR/mem$memtag/$fnmember_fro \n";
   print "reset_time.x $ENSDIR/mem$memtag/$fnmember_fro $nymd_fro $nhms_fro -9 \n";
   if ( ! $debug ){
     cp("$ANADIR/$fnmember_now","$ENSDIR/mem$memtag/$fnmember_fro");
     $rc =  system("reset_time.x $ENSDIR/mem$memtag/$fnmember_fro $nymd_fro $nhms_fro -9" );
     die ">>> ERROR <<< cannot reset time of member mem$memtag " if ( $rc );
   }

}

#.......................................
sub link_ensemble2analysis {

  return unless ($ANADIR);

  $local = `pwd`; chomp($local);
  print "cd $ANADIR\n"; 
  chdir("$ANADIR");

  foreach $dirmem (`ls -d $ENSDIR/mem*`) {
    chomp($dirmem);
    unlink($dirmem) if ( -d "$dirmem" ); # need this here since Assignfn is really for files (not dirs)
    $cmd = "ln -sf $dirmem . "; 
    print "$cmd \n"; 
    if (! $debug) {
      system($cmd);
    }
  }
  print "cd $local\n"; 
  chdir("$local");

}

#.......................................
sub usage {

   print <<"EOF";

NAME
     berror.pl - update ensemble of fields used to construct background error covariance

SYNOPSIS

     berror.pl [options] what

DESCRIPTION

    This procedure is to be used while running GSI. It updates an existing ensemble of
    fields that are used by GSI to construct its background error covariance. Typically,
    experiments using this feature should begin by placing an ensemble of fields under
    FVHOME/aens4berror; at each analysis time, the present procedure updates the ensemble
    by tossing one of its members and adding a new one from a recent GCM forecast. 

    At present, the ensemble is made up of 6-hour background fields. That is, the background
    error covariance matrix created with GCM is
     
          B = 1/(m-1) sum_m ( xb(m) - xb_bar ) ( xb(m) - xb_bar )

    where xb_bar is the mean of the background xb(m), with m representing the ensemble
    members. The members are made up of consecutive 6-hour background fields, kept for 
    a time lag of M*6 hours, where M is the size of the ensemble.

    Required Arguments when Updating Ensemble:
   
     expid  - experiment identifier
     nymd   - date of analysis (as in YYYYMMDD)
     nhms   - time of analysis (as in HHMMSS)
     freq   - analysis frequency (in hours)
     what   -    set: before GSI runs - make existing ensemble available to GSI 
              update: updates ensemble by tossing oldest member and add new one
                      (this requires more entries)

    Optional Arguments:

    -h      help (echoes this usage)
    -adir   name of diretory where analysis is running
            (NO DEFAULT; necessary to establish link with analysis)
    -edir   name of diretory where ensemble resides
            (default, FVHOME/aens4berror)   
    -rc     name of resource file defining background
            file name template (default template: expid.bkg.eta.nymd_hhz.nc4)

    Optional Env Vars:
      FVHOME  - experiment home directory

EOF

  exit(1);
}

