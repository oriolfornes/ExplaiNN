#!/usr/bin/perl -w
# FILE: pfmconvert
# CREATE DATE: 3/07/2014
# AUTHOR: Giovanna Ambrosini 
#
# Part of the code is based on an implemetation by
# William Stafford Noble and Timothy L. Bailey
# Created in 1999
# ORIG: transfac2meme.pl

# DESCRIPTION: Convert a Position Frequency Matrix (PFM) file to a Position Weihgt Matrix (PWM) with log odds weights.

use Scalar::Util::Numeric qw(isnum isint isfloat);
use Math::Round;

# Set up global variables. Assume uniform.
my @bases = ("A", "C", "G", "T");
my $num_bases = 4;
my %bg = ();
$bg{"A"} = 0.25;
$bg{"C"} = 0.25;
$bg{"G"} = 0.25;
$bg{"T"} = 0.25;

#my $c = 0.0001;			# default pseudocount fraction
my $c = 0;				# default pseudocount fraction
my $logscl = 100;
my $minscore = -10000;

my $usage = "USAGE: pfmconvert.pl [options] <matrix file>

  Options: 
	   -bg <background file>	set of f_a
           -c <pseudo weight>           add an arbitrary pseudo weight fraction <c> to each freq
					default: $c
           -m <low value score>         set lowest value score
                                        default: -10000
           -mswap                       by default colums represent nucleotides A,C,G,T
                                        if set matrix raws represent nucleotides A,C,G,T
           -n <log scaling factor>      set log scaling factor (int)
                                        default: 100
           -noheader                    write raw matrix (without header)
           -o <outfile>                 output file name
                                        default: no output file
           -w <l>|<p>                   <l> select log-odds matrix format for written output
                                        <p> select letter-probability matrix format for written output
                                        default format: log-odds matrix
  
  Convert a Position Frequency Matrix (PFM) file to a Position Weihgt Matrix (PWM) with either integer log odds weights or letter probabilities.
  A position frequency matrix (PFM) records the position-dependent frequency (or counts) of each residue or nucleotide.\n\n";


my $logodds = 1;
my $letprob = 0;
my $out_file = "";
my $ofile = 0;
my $m_swap = 0;


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
  } elsif ($next_arg eq "-mswap") {
    $m_swap = 1;
  } elsif ($next_arg eq "-n") {
    $logscl = shift(@ARGV);
  } elsif ($next_arg eq "-noheader") {
    $header = 0;
  } elsif ($next_arg eq "-o") {
    $out_file = shift(@ARGV);
    $ofile = 1;
  } elsif ($next_arg eq "-w") {
    my $format = shift(@ARGV);
    if ($format eq "l") {
      $logodds = 1;
    } elsif ($format eq "p") {
      $logodds = 0;
      $letprob = 1;
    }
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
my $float_flag = 0;
while ($line = <MF>) {
  # It is a Count Matrix: nucleotides are columns and position-dependent occurrencies are rows
  $i_motif = 0;
  $i_pos = 0;
  $i_base = 0;
  # Check header (new matrix)
  if ($line  =~ /^#/ or $line  =~ /^>/ or $first) {
    # Have we reached a new matrix?
    #print "$line\n";
    $first = 0;
    if ($line  =~ /^#/ or $line  =~ /^>/) {
      ($matrix_name) = $line =~/^[\#\>]\s*(\S+\s*\S*).*/;
      $matrix_name =~ s/\r|\n//g;
      if ($matrix_name eq "") {
        $matrix_name = "Unknown";
      }
    } else {
      $matrix_name = "Unknown";
      @counts = split(' ', $line);
      if ($m_swap) { # Raws are nucleotides ACGT, colums are positions
        $i_motif = scalar(@counts);
        # Store the contents of this row.
        for ($i_pos = 0; $i_pos < $i_motif;  $i_pos++) {
          $motif{$i_base, $i_pos} = shift(@counts);
          #print "  motif{$i_base, $i_motif} : $motif{$i_base, $i_motif}";
          if (isfloat($motif{$i_base, $i_pos})) {
            $float_flag = 1;
          }
        }
        $i_base++;
      } else {  # Colums are nucleotides ACGT, rwas are positions
        # Store the contents of this row.
        for ($i_base = 0; $i_base < $num_bases; $i_base++) {
	  $motif{$i_base, $i_motif} = shift(@counts);
          if (isfloat($motif{$i_base, $i_motif})) {
            $float_flag = 1;
          }
        }
        $i_motif++;
      }
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
      
      if ($m_swap) {  # Raws are nucleotides ACGT, colums are positions
        $i_motif = scalar(@counts);
        # Store the contents of this row.
        for ($i_pos = 0; $i_pos < $i_motif;  $i_pos++) {
          $motif{$i_base, $i_pos} = shift(@counts);
          #print "  motif{$i_base, $i_motif} : $motif{$i_base, $i_motif}";
          if (isfloat($motif{$i_base, $i_pos})) {
            $float_flag = 1;
          }
        }
        $i_base++;
      } else { # Colums are nucleotides ACGT, rwas are positions
      # Store the contents of this row.
        for ($i_base = 0; $i_base < $num_bases; $i_base++) {
    	  $motif{$i_base, $i_motif} = shift(@counts);
          if (isfloat($motif{$i_base, $i_motif})) {
            $float_flag = 1;
          }
        }
        $i_motif++;
      }
      if ($curspos ne $prev_curspos) {
       $prev_curspos = $curspos;
       #print "PREV CURSOR POS : $prev_curspos\n";
      }
    } # END OF WHILE ON SINGLE MATRIX
    $len = $i_motif;
    print "ID $matrix_name len=$len\n";
    for ($i_motif = 0; $i_motif < $len; $i_motif++) {
      for ($i_base = 0; $i_base < $num_bases; $i_base++) {
        if ($float_flag) {
          printf ("%7.2f ", $motif{$i_base, $i_motif});
        } else {
          printf ("%7d ", $motif{$i_base, $i_motif});
        }
      }
      print "\n";
    }
    print "\n";
    # Convert the motif to frequencies.
    # If $c != 0, add pseudocount fraction $c
    for ($i_motif = 0; $i_motif < $len; $i_motif++) {
      # motif columns may have different counts
      $num_seqs = 0;
      for ($i_base = 0; $i_base < $num_bases; $i_base++) {
        $num_seqs += $motif{$i_base, $i_motif};
      }
      for ($i_base = 0; $i_base < $num_bases; $i_base++) {
         $motif{$i_base, $i_motif} = 
         ($motif{$i_base, $i_motif} + ($bg{$bases[$i_base]} * $num_seqs * $c) ) / 
         ($num_seqs * (1 + $c));
      }
    }
    # Print the Letter-probability matrix (PSFM)
    if ($header and $ofile and $letprob) {
        print(OF ">letter-probability matrix $matrix_name: ");
        print(OF "alength= $num_bases w= $len nsites= $num_seqs E= 0\n");
    }
    print("letter-probability matrix $matrix_name: ");
    print("alength= $num_bases w= $len nsites= $num_seqs E= 0\n");
    for ($i_motif = 0; $i_motif < $len; $i_motif++) {
      for ($i_base = 0; $i_base < $num_bases; $i_base++) {
        if ($ofile and $letprob) {
          printf(OF "  %8.6f\t", $motif{$i_base, $i_motif});
        }
        printf("  %8.6f\t", $motif{$i_base, $i_motif});
      }
      if ($ofile and $letprob) {
        print(OF "\n");
      }
      print("\n");
    }
    print("\n");

    # Print the PWM.
    $num_motifs++;
    print(STDERR "Printing motif $matrix_name.\n");
    print(STDERR "MOTIF $num_motifs $matrix_name\n\n");
    print(STDERR "BL   MOTIF $num_motifs width=$len seqs=$num_seqs\n");

    # PSSM
    if ($header and $ofile and $logodds) {
      print(OF ">log-odds matrix $matrix_name: alength= $num_bases w= $len n= 0 bayes= 0 E= 0\n");
    }
    print("log-odds matrix $matrix_name: alength= $num_bases w= $len n= 0 bayes= 0 E= 0\n");
    for ($i_motif = 0; $i_motif < $len; $i_motif++) {
      for ($i_base = 0; $i_base < $num_bases; $i_base++) {
        if ($ofile and $logodds) {
          if ($motif{$i_base, $i_motif} > 0) {
            printf(OF "%7d ", round((log( $motif{$i_base, $i_motif} / $bg{$bases[$i_base]} )/log(2.0))*$logscl) );
          } else {
            printf(OF "%7d ", $minscore);
          }
        }
        if ($motif{$i_base, $i_motif} > 0) {
          printf("%7d ", round((log( $motif{$i_base, $i_motif} / $bg{$bases[$i_base]} )/log(2.0))*$logscl) );
        } else {
          printf("%7d ", $minscore);
        }
      }
      if ($ofile and $logodds) {
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
