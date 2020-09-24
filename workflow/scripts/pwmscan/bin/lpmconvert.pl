#!/usr/bin/perl -w
# FILE: pfmconvert
# CREATE DATE: 3/07/2014
# AUTHOR: Giovanna Ambrosini 
#
#
# Part of the code is based on an implemetation by
# William Stafford Noble and Timothy L. Bailey
# Created in 1999
# ORIG: transfac2meme.pl
#
# DESCRIPTION: Convert a Letter Probability Matrix (LPM) file to a Position Weihgt Matrix (PWM) with log odds weights.

use Math::Round;

# Set up global variables. Assume uniform.
my @bases = ("A", "C", "G", "T");
my $num_bases = 4;
my %bg = ();
$bg{"A"} = 0.25;
$bg{"C"} = 0.25;
$bg{"G"} = 0.25;
$bg{"T"} = 0.25;
#my $c = 0.0001;			# default pseudo-weight
my $c = 0;				# default pseudo-weight
my $logscl = 100;
my $minscore = -10000;

my $usage = "USAGE: lpmconvert.pl [options] <matrix file>

  Options: 
	   -bg <background file>	set of f_a
	   -c <pseudo weight>		add pseudo weight fraction <c> distributed according to residue priors
					default: $c
           -m <low value score>         set lowest value score (if set, default pseudo-weight=0)
                                        default: $minscore
           -n <log scaling factor>      set log scaling factor (int)
                                        default: 100
           -noheader                    write raw matrix (without header)
           -o <outfile>                 output file name
                                        default: no output file
  
  Convert a Letter Probability Matrix (LPM) file to a Position Weihgt Matrix (PWM) with integer log odds weights.
  A letter probability matrix (LPM) records the position-dependent probability normalized to 1 of each residue or nucleotide.\n\n";


my $out_file = "";
my $ofile = 0;


# Process command line arguments.
if (scalar(@ARGV) == 0) {
  printf(STDERR $usage);
  exit(1);
}
my $species = "";
my $id_list = "";
my $header = 1;
while (scalar(@ARGV) > 1) {
  $next_arg = shift(@ARGV);
  if ($next_arg eq "-bg") {
    $bg_file = shift(@ARGV);
  } elsif ($next_arg eq "-c") {
    $c = shift(@ARGV);
  } elsif ($next_arg eq "-m") {
    $minscore = shift(@ARGV);
  } elsif ($next_arg eq "-n") {
    $logscl = shift(@ARGV);
  } elsif ($next_arg eq "-noheader") {
    $header = 0;
  } elsif ($next_arg eq "-o") {
    $out_file = shift(@ARGV);
    $ofile = 1;
  } else {
    print(STDERR "Illegal argument ($next_arg)\n");
    exit(1);
  }
}
($matrix_file) = @ARGV;

# read the background file
if (defined($bg_file)) {
  open($bg_file, "<$bg_file") || die("Can't open $bg_file.\n");
  $total_bg = 0;
  while (<$bg_file>) {
    next if (/^#/);			# skip comments
    ($a, $f) = split;
    if ($a eq "A" || $a eq "a") {
      $bg{"A"} = $f; 
      $total_bg += $f;
    } elsif ($a eq "C" || $a eq "c") {
      $bg{"C"} = $f; 
      $total_bg += $f;
    } elsif ($a eq "G" || $a eq "g") {
      $bg{"G"} = $f; 
      $total_bg += $f;
    } elsif ($a eq "T" || $a eq "t") {
      $bg{"T"} = $f; 
      $total_bg += $f;
    }
  }
  # make sure they sum to 1
  foreach $key (keys %bg) {
    $bg{$key} /= $total_bg;
    #printf STDERR "$key $bg{$key}\n";
  }
}  # background file

# Open the matrix file for reading.
open(MF, "<$matrix_file") || die("Can't open $matrix_file.\n");
# Open the output file
if ($ofile) {
  open(OF, ">$out_file") || die("Can't open $out_file.\n");
}

# Print the MEME header.
print("ALPHABET= ACGT\n\n");
print("strands: + -\n\n");
print("Background letter frequencies (from dataset with add-one prior applied):\n");
printf("A %f C %f G %f T %f\n\n",  $bg{"A"}, $bg{"C"}, $bg{"G"}, $bg{"T"});

# Read the input file.
$num_motifs = 0;
my $i_motif = 0;
my $len = 0;
my $i_base = 0;
my $num_seqs = 0;
my $matrix_name = "";
my $curspos = 0;
my $prev_curspos = 0;
my $first = 1;
while ($line = <MF>) {
  # It is a Count Matrix: nucleotides are columns and position-dependent occurrencies are rows
  $i_motif = 0;
  # Check header (new matrix)
  if ($line  =~ /^#/ or $line  =~ /^>/ or $first) {
    # Have we reached a new matrix?
    $first = 0;
    if ($line  =~ /^#/ or $line  =~ /^>/) {
      if ($line  =~ /^>letter-probability matrix (.*):.*/) {
        $matrix_name = $1;
      } else {
        ($matrix_name) = $line =~/^[\#\>]\s*(\S+).*/;
      }
      if ($matrix_name eq "") {
        $matrix_name = "Unknown";
      }
    } else {
      $matrix_name = "Unknown";
      @counts = split(' ', $line);
      # Store the contents of this row.
      for ($i_base = 0; $i_base < $num_bases; $i_base++) {
        $motif{$i_base, $i_motif} = shift(@counts);
        chomp $motif{$i_base, $i_motif};
      }
      $i_motif++;
    }
    # Read the motif.
    while ($line = <MF>) {
      $curspos = tell(MF);

      if (! defined $line ) {
        print "\n";
        last;
      }
      #print $line;
      @counts = split(' ', $line);

      if ($line =~ />/) {
        seek(MF, $prev_curspos, 0);
        #print("Moving cursor to POS $prev_curspos\n");
        last;
      }
      
      # Store the contents of this row.
      for ($i_base = 0; $i_base < $num_bases; $i_base++) {
	$motif{$i_base, $i_motif} = shift(@counts);
        chomp $motif{$i_base, $i_motif};
      }
      $i_motif++;
      if ($curspos ne $prev_curspos) {
       $prev_curspos = $curspos;
       #print "PREV CURSOR POS : $prev_curspos\n";
      }
    } # END OF WHILE ON SINGLE MATRIX
    $len = $i_motif;
    print "ID $matrix_name len=$len\n";
    for ($i_motif = 0; $i_motif < $len; $i_motif++) {
      for ($i_base = 0; $i_base < $num_bases; $i_base++) {
        printf ("%9.6f ", $motif{$i_base, $i_motif});
      }
      print "\n";
    }
    print "\n";

    # Print the PWM.
    $num_motifs++;
    print(STDERR "Printing motif $matrix_name.\n");
    print(STDERR "MOTIF $num_motifs $matrix_name\n\n");
    print(STDERR "BL   MOTIF $num_motifs width=$len seqs=$num_seqs\n");

    # PSSM
    if ($header and $ofile) {
      print(OF ">log-odds matrix $matrix_name: alength= $num_bases w= $len n= 0 bayes= 0 E= 0\n");
    }
    print("log-odds matrix $matrix_name: alength= $num_bases w= $len n= 0 bayes= 0 E= 0\n");
    for ($i_motif = 0; $i_motif < $len; $i_motif++) {
      for ($i_base = 0; $i_base < $num_bases; $i_base++) {
        if ($ofile) {
          if ($c > 0) {
            printf(OF "%7d ", round((log( ($motif{$i_base, $i_motif} + $bg{$bases[$i_base]} * $c) / ($bg{$bases[$i_base]} * (1 + $c)) )/log(2.0))*$logscl) );
          } else {
            if ($motif{$i_base, $i_motif} > 0) {
              printf(OF "%7d ", round((log( $motif{$i_base, $i_motif} / $bg{$bases[$i_base]} )/log(2.0))*$logscl) );
            } else {
              printf(OF "%7d ", $minscore);
            }
          }
        }
        if ($c > 0) {
            printf("%7d ", round((log( ($motif{$i_base, $i_motif} + $bg{$bases[$i_base]} * $c) / ($bg{$bases[$i_base]} * (1 + $c)) )/log(2.0))*$logscl) );
        } else {
          if ($motif{$i_base, $i_motif} > 0) {
            printf("%7d ", round((log( $motif{$i_base, $i_motif} / $bg{$bases[$i_base]} )/log(2.0))*$logscl) );
          } else {
            printf("%7d ", $minscore);
          }
        }
      }
      if ($ofile) {
        print(OF "\n");
      }
      print("\n");
    }
    print("\n");
  }
}
print(STDERR "Converted $num_motifs motifs.\n");

close(MF);
if ($ofile) {
  close(OF);
}
