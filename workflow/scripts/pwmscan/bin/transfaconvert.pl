#!/usr/bin/perl -w
# FILE: transfaconvert
# CREATE DATE: 8/31/2013
# AUTHOR: Giovanna Ambrosini 
#
# Part of the code is based on an implemetation by
# William Stafford Noble and Timothy L. Bailey
# Created in 1999
# ORIG: transfac2meme.pl
# DESCRIPTION: Convert a Transfac matrix file to MEME output format.

# Giovanna Ambrosini 18/10/2017 
# Add pseudo weight fraction to correct frequencies
# 

use Math::Round;

# Set up global variables. Assume uniform.
my @bases = ("A", "C", "G", "T");
my $num_bases = 4;
my %bg = ();
$bg{"A"} = 0.25;
$bg{"C"} = 0.25;
$bg{"G"} = 0.25;
$bg{"T"} = 0.25;
my $c = 0;				# default total pseudocounts
my $logscl = 100;
my $minscore = -10000;

my $usage = "USAGE: transfaconvert.pl [options] <matrix file>

  Options: -species <name>
           -skip <transfac ID> (may be repeated)
           -ids <file containing list of transfac IDs>
	   -bg <background file>	set of f_a
           -c <pseudo weight>           add an arbitrary pseudo weight fraction <c> to each freq
                                        default: $c
           -m <low value score>         set low value score
                                        default: $minscore
           -n <log scaling factor>      set log scaling factor (int)
                                        default: 100
           -noheader                    write raw matrix (without header)
           -o <outfile>                 output file name
                                        default: no output file
           -w <l>|<p>                   <l> select log-odds matrix format for written output
                                        <p> select letter-probability matrix format for written output
                                        default format: original TRANSFAC frequency matrix
  
  Convert a Transfac matrix file to MEME output format (i.e. integer log likelihoods or letter-probability matrix). 
  N.B. Dollar signs in TRANSFAC IDs are converted to underscores.\n";


my $logodds = 0;
my $letprob = 0;
my $defout  = 0;
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
  if ($next_arg eq "-species") {
    $species = shift(@ARGV);
  } elsif ($next_arg eq "-skip") {
    $skips{shift(@ARGV)} = 1;
  } elsif ($next_arg eq "-ids") {
    $id_list = shift(@ARGV);
  } elsif ($next_arg eq "-bg") {
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
  } elsif ($next_arg eq "-w") {
    my $f = shift(@ARGV);
    if ($f eq "l") {
      $logodds = 1;
    } elsif ($f eq "p") {
      $letprob = 1;
    } else {
      print(STDERR "Illegal argument ($next_arg) $f\n");
      exit(1);
    }
  } else {
    print(STDERR "Illegal argument ($next_arg)\n");
    exit(1);
  }
}
($matrix_file) = @ARGV;

if (($logodds or $letprob) and $ofile == 0) {
  print(STDERR "Please, give a filename for your output file (-o <outfile>)\n");
  exit(1);
}
if ($logodds == 0 and $letprob == 0 and $ofile == 1) {
  $defout = 1;
}
# Store the target IDs.
%id_list = &read_list_from_file($id_list);

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
$num_skipped = 0;
my $i_motif = 0;
my $width = 0;
my $i_base = 0;
my $num_seqs = 0;
while ($line = <MF>) {

  # Split the line into identifier and everything else.
  ($id, @data) = split(' ', $line);

  # Have we reached a new matrix?
  if ($id eq "ID") {
    $matrix_name = shift(@data);
    $matrix_name =~ tr/\$/_/;

    # Read to the beginning of the motif.
    $this_species = "";
    # Old versions of TRANSFAC use pee-zero; new use pee-oh.
    while (($id ne "PO") && ($id ne "P0")) {
      $line = <MF>;

      if (! defined($line)) {
	die ("Can't find PO line for TRANSFAC matrix $matrix_name.\n");
      }
      ($id, @data) = split(' ', $line);

      # Store the species line.
      if ($id eq "BF") {
	$this_species .= $line;
      }
    }

    # Read the motif.
    $i_motif = 0;
    while () {
      $line = <MF>;

      if (! defined $line ) {
	die ("Can't find `XX' line for TRANSFAC matrix $matrix_name.\n");
      }

      ($id, @counts) = split(' ', $line);

      if ($id eq "XX") {
	last;
      }

      # Store the contents of this row.
      for ($i_base = 0; $i_base < $num_bases; $i_base++) {
	$motif{$i_base, $i_motif} = shift(@counts);
      }
      
      $i_motif++;
    } # END OF WHILE ON SINGLE MATRIX
    $width = $i_motif;
    if ($defout) {
      print OF "ID $matrix_name len=$width\n";
      print OF "BF\n";
      print OF "PO      A       C       G       T\n";
    }
    print "ID $matrix_name len=$width\n";
    print "BF\n";
    print "PO      A       C       G       T\n";
    for ($i_motif = 0; $i_motif < $width; $i_motif++) {
      for ($i_base = 0; $i_base < $num_bases; $i_base++) {
        if ($defout) {
          if ($i_base == 0) {
            printf (OF "%02d", $i_motif);
          } 
          printf (OF "%7d ", $motif{$i_base, $i_motif});
        }
        if ($i_base == 0) {
          printf ("%02d", $i_motif);
        } 
        printf ("%7d ", $motif{$i_base, $i_motif});
      }
      if ($defout) {
        print OF "\n";
      }
      print "\n";
    }
    if ($defout) {
      print OF "XX\n";
    }
    print "XX\n";
    print "\n";
    # If $c != 0, add pseudocount fraction $c
    for ($i_motif = 0; $i_motif < $width; $i_motif++) {
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
    # Convert the motif to frequencies.
#   Old stuff
#    for ($i_motif = 0; $i_motif < $width; $i_motif++) {
#      # motif columns may have different counts
#      $num_seqs = 0;
#      for ($i_base = 0; $i_base < $num_bases; $i_base++) {
#	$num_seqs += $motif{$i_base, $i_motif};
#      }
#      for ($i_base = 0; $i_base < $num_bases; $i_base++) {
#	$motif{$i_base, $i_motif} = 
#          ($motif{$i_base, $i_motif} + ($b * $bg{$bases[$i_base]}) ) / 
#          ($num_seqs + $b);
#      }
#    }

    ###### Decide whether to print the motif.

    # If no criteria are given, then print it.
    $print_it = 1;

    # If we were given a list,
    if ($id_list ne "") {
      #  was this matrix in the list?
      $print_it = defined($id_list{$matrix_name});
    }

    # If we were given a species.
    elsif ($species ne "") {
      # is this the right species?
      $print_it = ($this_species =~ m/$species/);
      if ($this_species eq "") {
	print(STDERR "Warning: No species given for $matrix_name.\n");
      }
    }

    # Were we explicitly asked to skip this one?
    if (defined($skips{$matrix_name})) {
      $print_it = 0;
    } 
    
    # Print the motif.
    if ($print_it) {
      $num_motifs++;
      print(STDERR "Printing motif $matrix_name.\n");
      print(STDERR "MOTIF $num_motifs $matrix_name\n\n");
      print(STDERR "BL   MOTIF $num_motifs width=$width seqs=$num_seqs\n");

      # PSSM
      if ($logodds) {
        if ($header) {
          print(OF ">log-odds matrix $matrix_name: alength= $num_bases w= $width n= 0 bayes= 0 E= 0\n");
        }
      }
      print("log-odds matrix $matrix_name: alength= $num_bases w= $width n= 0 bayes= 0 E= 0\n");
      for ($i_motif = 0; $i_motif < $width; $i_motif++) {
        for ($i_base = 0; $i_base < $num_bases; $i_base++) {
          if ($logodds) {
            if ($motif{$i_base, $i_motif}) {
              printf(OF "%7d ", round((log( $motif{$i_base, $i_motif} / $bg{$bases[$i_base]} )/log(2.0))*$logscl) );
            } else {
              printf(OF "%7d ", $minscore);
            }
          }
          if ($motif{$i_base, $i_motif}) {
            printf("%7d ", round((log( $motif{$i_base, $i_motif} / $bg{$bases[$i_base]} )/log(2.0))*$logscl) );
          } else {
              printf("%7d ", $minscore);
          }
        }
        if ($logodds) {
          print(OF "\n");
        }
        print("\n");
      }
      print("\n");

      # PSFM
      if ($letprob) {
        if ($header) {
          print(OF ">letter-probability matrix $matrix_name: ");
          print(OF "alength= $num_bases w= $width nsites= $num_seqs E= 0\n");
        }
      }
      print("letter-probability matrix $matrix_name: ");
      print("alength= $num_bases w= $width nsites= $num_seqs E= 0\n");
      for ($i_motif = 0; $i_motif < $width; $i_motif++) {
        for ($i_base = 0; $i_base < $num_bases; $i_base++) {
          if ($letprob) {
            printf(OF "  %8.6f\t", $motif{$i_base, $i_motif});
          }
          printf("  %8.6f\t", $motif{$i_base, $i_motif});
        }
        if ($letprob) {
          print(OF "\n");
        }
        print("\n");
      }
      print("\n");
    } else {
      $num_skipped++;
      #print(STDERR "Skipping motif $matrix_name.\n");
    }
  }
}
print(STDERR "Converted $num_motifs motifs.\n");
print(STDERR "Skipped $num_skipped motifs.\n");

close(MF);
if ($ofile) {
  close(OF);
}

sub read_list_from_file {
  my($id_list) = @_;
  my($line, %return_value);

  if ($id_list eq "") {
    return(%return_value);
  }

  open($id_list, "<$id_list") || die("Can't open $id_list.");
  while ($line = <$id_list>) {
    chomp($line);
    @words = split(' ', $line);
    foreach $word (@words) {
      $return_value{$word} = 1;
    }
  }
  close($id_list);

  return(%return_value);
}
