#!/usr/bin/perl -w
# FILE: pwm2lpmconvert.pl 
# CREATE DATE: 21/02/2018
# AUTHOR: Giovanna Ambrosini 
#

# DESCRIPTION: Convert a position weigh matrix (PWM) to a letter probability matrix
#              assuming, for the background model, an independent succession of 
#              nucleotides (Bernoulli model). Use the natural logarithm.

use List::Util qw[min max sum];
use Scalar::Util::Numeric qw(isnum isint isfloat);
use Math::Round;

# Set up global variables. Assume uniform.
my @bases = ("A", "C", "G", "T");
my $num_bases = 4;

my $logscl = 1;

my $usage = "USAGE: pwm2pfmconvert.pl [options] <matrix file>

  Options:
           -n <log scaling factor>      set log scaling factor (int) default: 1
           -noheader                    write raw matrix (without header)
           -o <outfile>                 output file name
                                        default: no output file
           -ssa                         input PWM matrix is in SSA-like format 

  Convert a Position Weight Matrix (PWM) to letter probability format. 
  The conversion formula assumes, for the background model, an independent
  succession of nucleotides (Bernoulli model), and the natural logarithm base.
  PWM weights can be normalized by a logarithm scaling factor (by default=1).
\n";


my $out_file = "";
my $ofile = 0;
my $inp_ssa_flag = 0;
# Process command line arguments.
if (scalar(@ARGV) == 0) {
  printf(STDERR $usage);
  exit(1);
}
my $header = 1;
while (scalar(@ARGV) > 1) {
  $next_arg = shift(@ARGV);
  if ($next_arg eq "-n") {
    $logscl = shift(@ARGV);
  } elsif ($next_arg eq "-noheader") {
    $header = 0;
  } elsif ($next_arg eq "-o") {
    $out_file = shift(@ARGV);
    $ofile = 1;
  } elsif ($next_arg eq "-ssa") {
    $inp_ssa_flag = 1;
  } else {
    print(STDERR "Illegal argument ($next_arg)\n");
    exit(1);
  }
}
($matrix_file) = @ARGV;

# Open the matrix file for reading.
open(MF, "<$matrix_file") || die("Can't open $matrix_file.\n");
# Open the output file
if ($ofile) {
  open(OF, ">$out_file") || die("Can't open $out_file.\n");
}

# Print the MEME-like header.
print("ALPHABET= ACGT\n\n");
print("strands: + -\n\n");

# Read the input file.
$num_motifs = 0;
$num_skipped = 0;
my $i_motif = 0;
my $width = 0;
my $i_base = 0;
my $coff = 0;
my $matrix_name = "";
my $print_it = 0;

my $min = 0;
my $max = 0;
my $tmp_max = 0;
my @mat = ();
my $k = 1;
my $float_flag = 1;
my @norm_sum = ();
while ($line = <MF>) {
  if (!$inp_ssa_flag) {
    # Matrix has only position score columns for each base
    print STDERR "RAW PWM format.\n"; 
    $i_motif = 0;
    if ($line  =~ /^#/ or $line  =~ /^>/) {
      # Have we reached a new matrix?
      ($matrix_name) = $line =~/^[\#\>]\s*(\S+).*/;
      # Read the motif.
      while (<MF>) {
        last if (/\/\//);
        (@scores) = split(' ', $_);
        # Store the contents of this row.
        for ($i_base = 0; $i_base < $num_bases; $i_base++) {
          $motif{$i_base, $i_motif} = shift(@scores);
          $mat[$k][$i_base] = $motif{$i_base, $i_motif};
          if (($mat[$k][$i_base] != 0) and isint($mat[$k][$i_base])) {
            $float_flag = 0;
          }
        }
        $i_motif++;
        $k++;
      }
    } else {
      $matrix_name = "Unknown motif";
      (@scores) = split(' ', $line);
      for ($i_base = 0; $i_base < $num_bases; $i_base++) {
        $motif{$i_base, $i_motif} = shift(@scores);
        $mat[$k][$i_base] = $motif{$i_base, $i_motif};
        if (($mat[$k][$i_base] != 0) and (isint $mat[$k][$i_base])) {
          $float_flag = 0;
        }
      }
      $i_motif++;
      $k++;
      # Read the motif.
      while (<MF>) {
        last if (/\/\//);
        (@scores) = split(' ', $_);
        # Store the contents of this row.
        for ($i_base = 0; $i_base < $num_bases; $i_base++) {
          $motif{$i_base, $i_motif} = shift(@scores);
          $mat[$k][$i_base] = $motif{$i_base, $i_motif};
          if (($mat[$k][$i_base] != 0) and (isint $mat[$k][$i_base])) {
            $float_flag = 0;
          }
        }
        $i_motif++;
        $k++;
      }
    }
    $width = $i_motif;
    print "> $matrix_name len=$width\n";
    for ($i_motif = 0; $i_motif < $width; $i_motif++) {
      for ($i_base = 0; $i_base < $num_bases; $i_base++) {
        if ($float_flag) {
          printf ("%7.2f ", $motif{$i_base, $i_motif});
        } else {
          printf ("%7d ", $motif{$i_base, $i_motif});
        }
      }
      print "\n";
    }
    print "//\n";
    print "\n";
    # Convert the motif to letter-probabilities and print it.
    @norm_sum = ();
    for ($i_motif = 0; $i_motif < $width; $i_motif++) {
      for ($i_base = 0; $i_base < $num_bases; $i_base++) {
        $motif{$i_base, $i_motif} = exp($motif{$i_base, $i_motif}/$logscl);
        $norm_sum[$i_motif] += $motif{$i_base, $i_motif};
      }
    }
    $print_it = 1;
    print("letter-probability matrix $matrix_name: alength=$num_bases w=$width n= 0 bayes= 0 E= 0\n");
    for ($i_motif = 0; $i_motif < $width; $i_motif++) {
      # motif columns may have different counts
      for ($i_base = 0; $i_base < $num_bases; $i_base++) {
        printf ("%8.6f ", $motif{$i_base, $i_motif}/$norm_sum[$i_motif]); 
      }
      print("\n");
    }
    print("\n");
    # Print the motif to file.
    if ($print_it) {
      $num_motifs++;
      if ($ofile) {
        print(STDERR "Printing motif $matrix_name to file $out_file.\n");
        print(STDERR "MOTIF $num_motifs $matrix_name width=$width\n\n");

        # PSSM
        if ($header) {
          print(OF ">letter-probability matrix $matrix_name: alength= $num_bases w= $width n= 0 bayes= 0 E= 0\n");
        }
        for ($i_motif = 0; $i_motif < $width; $i_motif++) {
          for ($i_base = 0; $i_base < $num_bases; $i_base++) {
            printf(OF "%8.6f ", $motif{$i_base, $i_motif}/$norm_sum[$i_motif]);
          }
          print(OF "\n");
        }
      }
    } else {
      $num_skipped++;
      #print(STDERR "Skipping motif $matrix_name.\n");
    }
  } else {
    # SSA format Matrix (inlcuding header and additional columns
    # Split the line into identifier and everything else.
    print STDERR "Standard SSA PWM format.\n"; 
    ($id, @data) = split(' ', $line);

    # Have we reached a new matrix?
    if ($id eq "TI") {
      print STDERR "ID = $id\n";
      $matrix_name = shift(@data);
      # Read to the beginning of the motif.
      while (($id ne "XX")) {
        $line = <MF>;

        if (! defined($line)) {
          die ("Can't find first XX line for SSA matrix $matrix_name.\n");
        }
        ($id, @data) = split(' ', $line);
      } 
      $id = "";
      while (($id ne "XX")) {
        $line = <MF>;

        if (! defined($line)) {
  	  die ("Can't find second XX line for SSA matrix $matrix_name.\n");
        }
        ($id, @data) = split(' ', $line);
        # Store the coff line.
        if ($id eq "CO") {
          $coff = shift(@data);
        }
      }
      # Read the motif.
      $i_motif = 0;
      while () {
        $line = <MF>;

        if (! defined $line ) {
  	  die ("Can't find `//' line for SSA matrix $matrix_name.\n");
        }

        ($id, @counts) = split(' ', $line);

        if ($id eq "//") {
  	  last;
        }
        # Store the contents of this row.
        $pos{$i_motif} = shift(@counts);
        for ($i_base = 0; $i_base < $num_bases; $i_base++) {
	  $motif{$i_base, $i_motif} = shift(@counts);
          $mat[$k][$i_base] = $motif{$i_base, $i_motif};
          if (($mat[$k][$i_base] != 0) and isint($mat[$k][$i_base])) {
            $float_flag = 0;
          }
        }
        $i_motif++;
        $k++;
      } # END OF WHILE ON SINGLE MATRIX
      $width = $i_motif;
      print "TI $matrix_name len=$width\n";
      print "XX\n";
      print "CO    $coff\n";
      print "XX\n";
      for ($i_motif = 0; $i_motif < $width; $i_motif++) {
        for ($i_base = 0; $i_base < $num_bases; $i_base++) {
          if ($i_base == 0) {
            printf ("%7d", $pos{$i_motif});
          }
          if ($float_flag) {
            printf ("%7.2f", $motif{$i_base, $i_motif});
          } else {
            printf ("%7d ", $motif{$i_base, $i_motif});
          }
        }
        print "\n";
      }
      print "//\n";
      print "\n";
      # Convert the motif to letter probabilities and print it.
      @norm_sum = ();
      for ($i_motif = 0; $i_motif < $width; $i_motif++) {
        for ($i_base = 0; $i_base < $num_bases; $i_base++) {
          $motif{$i_base, $i_motif} = exp($motif{$i_base, $i_motif}/$logscl);
          $norm_sum[$i_motif] += $motif{$i_base, $i_motif};
        }
      }
      $print_it = 1;
      print("letter-probability matrix $matrix_name: alength=$num_bases w=$width n= 0 bayes= 0 E= 0\n");
      for ($i_motif = 0; $i_motif < $width; $i_motif++) {
        # motif columns may have different counts
        for ($i_base = 0; $i_base < $num_bases; $i_base++) {
  	  printf ("%8.6f ", $motif{$i_base, $i_motif}/$norm_sum[$i_motif]); 
        }
        print("\n");
      }
      print("\n");
    
      # Print the motif to file.
      if ($print_it) {
        $num_motifs++;
        if ($ofile) {
          print(STDERR "Printing motif $matrix_name to file $out_file.\n");
          print(STDERR "MOTIF $num_motifs $matrix_name width=$width\n\n");

          # PSSM
          if ($header) {
            print(OF ">letter-probability matrix $matrix_name: alength= $num_bases w= $width n= 0 bayes= 0 E= 0\n");
          }
          for ($i_motif = 0; $i_motif < $width; $i_motif++) {
            for ($i_base = 0; $i_base < $num_bases; $i_base++) {
              printf(OF "%8.6f ", $motif{$i_base, $i_motif}/$norm_sum[$i_motif]);
            }
            print(OF "\n");
          }
        }
      } else {
        $num_skipped++;
        #print(STDERR "Skipping motif $matrix_name.\n");
      }
    } # END ON MATRIX HEADER TI
  } # END ON MATRIX FORMAT : Raw vs Standard SSA
}
print(STDERR "Converted $num_motifs motifs.\n");
print(STDERR "Skipped $num_skipped motifs.\n");

close(MF);
if ($ofile) {
  close(OF);
}
