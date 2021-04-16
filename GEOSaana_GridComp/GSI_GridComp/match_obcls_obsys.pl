#!/usr/bin/env perl

# match_obcls_obsys.pl
#
#  03oct2016 Todling  Created to overcome MPT mpit problem with GSI system call
#
use Env;                 # make env vars readily available
use File::Basename;      # for basename(), dirname()
use File::Copy "cp";     # for cp()
use Getopt::Long;        # load module with GetOptions function
use Time::Local;         # time functions
use FindBin;             # so we can find where this script resides
use List::MoreUtils qw(uniq);

# look for perl packages in the following locations
#--------------------------------------------------
use lib ( "$FindBin::Bin", "$FVROOT/bin", "$ESMADIR/$ARCH/bin" );

GetOptions ( "debug",
             "dir=s",
             "h" );

# FVROOT is where the binaries have been installed
# ------------------------------------------------
$fvroot  = dirname($FindBin::Bin);
$fvroot  =~ s|/u/.realmounts/share|/share|;   # for portability across
                                              # NAS machines

  init();

  run();

#..................................................................
sub init {

  if ( $#ARGV < 3 || $opt_h ) {
       print STDERR "missing entries in argument list; see usage";
       usage();
  } else {              # required command lile args
       $nymd  = $ARGV[0];
       $nhms  = sprintf("%6.6d",$ARGV[1]);
       $gridcomprc = $ARGV[2];
       $gsirc      = $ARGV[3];
  }
  chomp($gridcomprc);
  chomp($gsirc);

  $obsdir = "Obsloc"; 
  if ( $opt_dir ) {
     $obsdir = $opt_dir;
  }

  if ($opt_debug) {
     print "$nymd \n";
     print "$nhms \n";
     print "$gridcomprc \n";
     print "$gsirc \n";
  }

}
#..................................................................
sub run {

my @want_obscls = `$fvroot/bin/echorc.x -rc $gridcomprc -ncol 1 observation_files`;
my @fname_tmpl  = `$fvroot/bin/echorc.x -template null $nymd $nhms -rc $gridcomprc -ncol 3 observation_files`;

my @gsirc_lnname = `$fvroot/bin/echorc.x -rc $gsirc -ncol 1 OBS_INPUT`;
my @gsirc_obscls = `$fvroot/bin/echorc.x -rc $gsirc -ncol 8 OBS_INPUT`;

$ii=-1;
foreach $var ( @fname_tmpl ) {
   chomp($var);
   $ii = $ii + 1;
   @array_tmpl[$ii] = $var; 
}
$ii=-1;
foreach $var ( @gsirc_lnname ) {
   chomp($var);
   $ii = $ii + 1;
   @array_lnname[$ii] = $var; 
}

my $iln = -1; my $ll = -1;
# loop over obs-classes in gsi.rc
foreach $var1 ( @want_obscls ) {
   chomp($var1);
   $iln = $iln + 1;
   # find a match with gridcomp
   $jln = -1;
   foreach $var2 ( @gsirc_obscls ) {
        chomp($var2);
        $jln = $jln + 1;
        if ($var1 eq $var2) {
           $ll = $ll + 1;
           $tmlist[$ll] = "$array_tmpl[$iln]";
           $lnlist[$ll] = "$array_lnname[$jln]";
           last:
        }
   }
}
@lnlist = uniq @lnlist;
@tmlist = uniq @tmlist;
$ii=-1;
foreach $var ( @lnlist ){
  chomp($var); $var =~ s/^\s+|\s+$//g;
  $ii = $ii + 1;
  if ( $opt_debug ){
     print "ln -s $obsdir/$tmlist[$ii] $var \n";
  } else {
     print "ln -s $obsdir/$tmlist[$ii]   $var \n";
        Assignfn("$obsdir/$tmlist[$ii]","$var");
  }
}

}
#............................................................................. 
sub Assign {

  my ( $fname, $lu ) = @_;

  $f77name = "fort.$lu";
  unlink($f77name) if ( -e $f77name ) ;
  symlink("$fname","$f77name");

}

sub Assignfn {

# Assignfn - assigns fn to given file name fname.
# fname = old file
# fn = new file (links to old)
  my ( $fname, $fn ) = @_;
  unlink($fn) if ( -e $fn ) ;
  symlink("$fname","$fn");

 }

sub fullpath {

    my ( $fname ) = @_;

    $dirn = dirname("$fname");
    chomp($dirn = `pwd`) if ( "$dirn" eq "." ) ;
    $name = basename("$fname");
    $full = "$dirn/$name";
    return ($full);

  }
#............................................................................. 

sub usage {

   print <<"EOF";

NAME
     match_obcls_obsys.pl - appends obs table to GSI_GridComp.rc

SYNOPSIS

     match_obcls_obsys.pl [options] nymd nhms GSI_GridComp.rc gsi.rc

DESCRIPTION

     Establish links between templated observation files names
     and GSI dfile names.

    Optinal Arguents:

    -h             help (echoes this usage)
    -debug         echo results, but does not do anything
    -dir   DIR     location where observations files to be found
                   (default: Obsloc)

    Required Env Vars:

    Optional Env Vars:

EOF

  exit(1);
}

