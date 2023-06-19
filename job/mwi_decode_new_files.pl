#!/usr/bin/env perl
#
# ##################################################################
#
# NAME:
# mwi_decode_new_files.pl
#
# PURPOSE:
# Script for finding new EUMETSAT MWI files and send to decoding, one by one.
#
# REQUIREMENTS:
#
# INPUT:
#
# OUTPUT:
#
# NOTES:
#
# BUGS:
#
# AUTHOR:
# Steinar Eastwood, met.no, 16.06.2023


sub usage {
    my($progname) = @_;

    print "\nSYNTAX: $progname (-i <inpdir> -o <outdir> -z <satid> -r -a -v -s -d <rdate>  -h)\n";
    print "<inpdir> : Optional input dir, if not using default value (SSMI_L2PG=$SSMI_L2PG).\n";
    print "<outdir> : Optional output dir, if not using default value (SSMI_L2PG=$SSMI_L2PG).\n";
    print "<satid>  : Using specific DMSP satellite, e.g. f18. Default: Uses everything in inpdir.\n";
    print "-r       : Also reprocess files that have edge/type probabilities already\n";
    print "-a       : Process all files, not check UPDATED file\n";
    print "-v       : Run in verbose mode, more log output\n";
    print "-h       : Print this help message.\n";
    print "-s       : Process all files with suffix _sn*.nc.\n";
    print "<rdate>  : Optional, specific date to process 'YYYYMMDD' (for reproc).\n";
    exit(1);
}


use strict;
use Time::Local;
use File::Path qw(make_path);
use Getopt::Std;
use vars qw($opt_r $opt_a $opt_i $opt_v $opt_D $opt_h $opt_p $opt_o $opt_x $opt_z $opt_s $opt_d);
sub printMSG($$);

my($ret,$err,$lpfmtime,$arglist,@files,$numf,$i,$pdf_coeff_fname_arc,$pdf_coeff_fname_ant);
my($fname,$fname_new,$fmtime,$newlpfmtime,$fdate,$dc,$protime,$prodate);
my($procflg,@fnsp,$reproc,$inpdir,$outdir,$proc_all,$verbose,$fname_arc,$fname_ant);
my($instname,$satname,@Ts);

#
my $command_decode = "python3 ./process_mwi.py";

#Defaults
my $inpdir = "";
my $outdir = "";
my $satid = 'sgb1';
my $period = '015';
my $daysahead = 2;
my $suffix = "snd"; 
my $rdate = '';
my $verbose = 0;
my $proc_all = 0;

# Get commandoline parameters (if any)
$opt_r = $opt_a = $opt_i = $opt_v = $opt_h = $opt_o = $opt_z = $opt_s = $opt_d = 0;
getopts('ai:vEho:z:sd:');
if ($opt_a) {
    printMSG "Process all files, not check UPDATED file",1;
    $proc_all = 1;
}
if ($opt_i) {
    $inpdir = $opt_i;
    $outdir = $opt_i;
    printMSG "Process files in this directory instead: $inpdir",1;
}
if ($opt_z) {
    $satid = lc($opt_z);
    printMSG "Uses only data from DMSP-$satid",1;
}
if ($opt_v) {
    printMSG "Running in verbose mode",1;
    $verbose = 1;
}
if ($opt_h) {
    usage($0);
    exit(1);
}
if ($opt_o) {
    $outdir = $opt_o;
    if ( !(-d $outdir) ) {
	printMSG "error: Could not find outdir: $outdir",0;
	exit(1);
    }
    printMSG "Place processed files in this directory instead: $outdir",1;
}
if ( $opt_s ) {
    $suffix = "sn";
    printMSG "Process all files with suffix _sn*.nc.",1;
}

if ( $opt_d ) {
    if (length($opt_d) != 8) {
        printMSG "Wrong date-format. Should be YYYYMMDD",1;
        exit(1);
    }
    $rdate = $opt_d;
    printMSG "Process product files for this date: $rdate",1;
}


my $updated = $outdir."/UPDATED_prob";

# Find mtime of last processed file
if (!$proc_all) {
    if (!(-e $updated)) {
        printMSG "Touching UPDATED since do not exist: $updated",1;
        $ret = system("touch $updated");
        if ($ret) {
            printMSG "ERROR($0): Could not touch file $updated",0;
            exit(2);
        }
        $lpfmtime = 1;
    }
    else {
        $lpfmtime = (stat($updated))[9];
    }
}
else {
    $lpfmtime = 1;
}
$newlpfmtime = $lpfmtime;

if ($verbose) {
    printMSG "Updated file: $updated",1;
    printMSG "lpfmtime: $lpfmtime",1;
}


# Loop through SSMIMWI NetCDF dir and find new files to process
opendir(SRC,$inpdir);
@files = grep /.*SGB1-MWI.*\.nc$/, readdir(SRC);
closedir(SRC);
@files = sort(@files);
$numf = $#files + 1;
printMSG "number of files: $numf",1;

for ($i=0;$i<$numf;$i++) {

    $fname   = $inpdir."/".$files[$i];
    $fmtime  = (stat($fname))[9];
    if ($verbose) {  printf "\t$fname $fmtime $lpfmtime\n";  }

    if ($fmtime <= $lpfmtime && !$proc_all) {
        if ($verbose) {  printf "\tskip since too old\n";  }
        next;
    }

    $arglist = "-i $fname -o $outdir";
    $ret = system("$command_decode $arglist");
    if ($ret) {
        printMSG "ERROR($0): in running $command_decode $arglist",0;
        $err ++;
        next;
    }
    else {
        if ($fmtime > $newlpfmtime) {
            $newlpfmtime = $fmtime;
        }
        if ($fname ne $fname_new) {
            rename $fname, $fname_new;
        }
    }
}

# Touch updated file with time of last processed file
if (!$proc_all) {
    printMSG "Touching $updated",1;
    #$ret = system("touch $updated");
    utime $newlpfmtime, $newlpfmtime, $updated;
}





#
# PURPOSE: 
# Write message to stdout or stderr with time stamp.
#
# SYNTAX:
# printMSG <message>,<destination>
#       <destination> == 0 or 2 => message goes to only STDERR
#       <destination> == 1 or others => message goes to STDOUT
#
# Format of output: "$logdate  $message\n"
#
# If first two characters of $message is "\n", a new line is printed 
# before $logdate.
#
sub printMSG($$) {
  my $message = shift;
  my $dest = shift;
  my($logdate,$ptime,@T,$prefix,$prefix2);

  $ptime = time;
  @T = gmtime(time);
  $logdate = sprintf("%04d.%02d.%02d %02d:%02d:%02d",
    $T[5]+1900,$T[4]+1,$T[3],$T[2],$T[1],$T[0]);
  $prefix = substr($message,0,1);
  if ($prefix eq "\n") {
    $message = substr($message,1);
  }
  else { 
    $prefix  = "";
  }

  if ($dest == 0 || $dest == 2) {
    print STDERR "$prefix$logdate  $message\n";
  }
  else {
    print STDOUT "$prefix$logdate  $message\n";
  }
}
