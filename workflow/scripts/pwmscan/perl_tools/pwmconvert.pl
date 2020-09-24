#!/usr/bin/perl -w
# FILE: pwmconvert
# CREATE DATE: 10/07/2015
# AUTHOR: Giovanna Ambrosini 
#
# Part of the code is based on an implemetation by
# William Stafford Noble and Timothy L. Bailey
# Created in 1999
# ORIG: transfac2meme.pl

# DESCRIPTION: Convert a SSA-formatted or plain-text PWM to an integer PWM.

use List::Util qw[min max sum];
use Scalar::Util::Numeric qw(isnum isint isfloat);
use Math::Round;

# Set up global variables. Assume uniform.
my @bases = ("A", "C", "G", "T");
my $num_bases = 4;

my $usage = "USAGE: pwmconvert.pl [options] <matrix file>

  Options: 
           -c                           check format : check whether PWM scores are Real or Integer 
                                        and print out the corresponding Real/Integer flag
           -n <scaling factor>          set scaling factor (int) for real format only
                                        If not set, the program estimates its value 
           -noheader                    write raw matrix (without header)
           -o <outfile>                 output file name
                                        default: no output file
           -ssa                         input PWM matrix is in SSA-like format 

  Convert a real or SSA-formatted Position Weight Matrix (PWM) to integer plain-text format. 
  Optionally, the program performs a check to see whether the matrix scores are Real or Integer
  numbers, and prints the corresponding Real/Integer flag.
\n";


my $out_file = "";
my $ofile = 0;
my $inp_ssa_flag = 0;
my $check_flag = 0;
my $scale = 1;
my $setscale_flag = 0;
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
  if ($next_arg eq "-c") {
    $check_flag = 1;
  } elsif ($next_arg eq "-n") {
    $scale = shift(@ARGV);
    if ($scale ne "") {
       $setscale_flag = 1; 
    } else {
       $scale = 1;
    }
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
    # Check whether the PWM is integer or real
    if ($check_flag) {
      if ($float_flag) {
        print "Real PWM\n";
      } else {
        print "Integer PWM\n";
      }
      exit (0);
    }
    # Estimate scaling factor
    if ($float_flag and ($setscale_flag == 0)) {
      for ($k=1; $k<=$width; $k++) {
        $min=min(@{$mat[$k]});
        $tmp_max=max(@{$mat[$k]});
        $max += $tmp_max-$min;
      }
      #print "MAX range : $max\n";
      if ($max <= 50 ) {
        $scale = 1000;
      } elsif ($max <= 500) {
        $scale = 100;
      } elsif ($max <= 10000) {
        $scale = 10;
      } elsif ($max > 10000) {
        $scale = 1;
      }
    }
    # Convert the motif to integer scores and print it.
    $print_it = 1;
    print("matrix $matrix_name: alength=$num_bases w=$width n= 0 bayes= 0 E= 0\n");
    for ($i_motif = 0; $i_motif < $width; $i_motif++) {
      # motif columns may have different counts
      for ($i_base = 0; $i_base < $num_bases; $i_base++) {
        printf ("%7d ", round($motif{$i_base, $i_motif} * $scale)); 
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
          print(OF ">matrix $matrix_name: alength= $num_bases w= $width n= 0 bayes= 0 E= 0\n");
        }
        for ($i_motif = 0; $i_motif < $width; $i_motif++) {
          for ($i_base = 0; $i_base < $num_bases; $i_base++) {
            printf(OF "%7d ", round($motif{$i_base, $i_motif} * $scale));
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
      if ($check_flag) {
        if ($float_flag) {
          print "Real PWM\n";
        } else {
          print "Integer PWM\n";
        }
        exit (0);
      }
      #print "Estimate scaling factor...\n";
      #print "Scale = $scale\n";
      #print "Float Flag = $float_flag\n";
      if ($float_flag and ($setscale_flag == 0)) {
        for ($k=1; $k<=$width; $k++) {
          $min=min(@{$mat[$k]});
          $tmp_max=max(@{$mat[$k]});
          $max += $tmp_max-$min;
        }
        #print "MAX range : $max\n";
        if ($max <= 50 ) {
          $scale = 100;
        } elsif ($max <= 1000) {
          $scale = 10;
        } elsif ($max > 1000) {
          $scale = 1;
        }
      }
      # Convert the motif to integer scores and print it.
      $print_it = 1;
      print("matrix $matrix_name: alength=$num_bases w=$width n= 0 bayes= 0 E= 0\n");
      for ($i_motif = 0; $i_motif < $width; $i_motif++) {
        # motif columns may have different counts
        for ($i_base = 0; $i_base < $num_bases; $i_base++) {
  	  printf ("%7d ", round($motif{$i_base, $i_motif} * $scale)); 
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
            print(OF ">matrix $matrix_name: alength= $num_bases w= $width n= 0 bayes= 0 E= 0\n");
          }
          for ($i_motif = 0; $i_motif < $width; $i_motif++) {
            for ($i_base = 0; $i_base < $num_bases; $i_base++) {
              printf(OF "%7d ", round($motif{$i_base, $i_motif} * $scale));
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
print(STDOUT "Scaling :\n");
print(STDOUT "$scale\n");

close(MF);
if ($ofile) {
  close(OF);
}
