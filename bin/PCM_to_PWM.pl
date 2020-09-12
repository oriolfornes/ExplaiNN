#!/usr/bin/perl -w
#
################################################################################
#
################################################################################

use strict;

################################################################################
# Main
################################################################################

#usage: perl PFM_to_PWM_JASPAR.pl /path/to/the/PFM/bundle/

use Bio::SeqIO;
use TFBS::SiteSet;
use Getopt::Long;
use TFBS::DB::JASPAR;

my $matrix_file = $ARGV[0];
GetOptions(
  '-f=s' => \$matrix_file,
);

my $matrix_set;
if ($matrix_file){
  $matrix_set = read_PFMs($matrix_file);
}

chdir './PWM';
my $pwm;
my $miter = $matrix_set->Iterator();
while (my $matrix = $miter->next){
  $pwm = $matrix->to_PWM();
  my $ID = $pwm->ID;
  my $name = $pwm->name;
  my $out_file = "${ID}_${name}.pfm";
  my $str = $pwm->rawprint();

  print $str or die "Cannot write to $out_file";
}

exit;


sub read_PFMs
{
    my ($file) = @_;

    open(FH, $file) || die "Error opening PFM file $file\n";

    my $matrix_set = TFBS::MatrixSet->new();

    my $name          = '';
    my $ID            = '';
    my $matrix_string = '';
    my $line_count    = 0;
    while (my $line = <FH>) {
        chomp $line;
        next if !$line;
        if ($line =~ /^>\s*(\S+)\t(\S+)/) {
            $ID = $1;
            $name = $2;
        } else {
            if ($line =~ /^\s*[ACGT]\s*\[\s*(.*)\s*\]/) {
                # line of the form: A [ # # # ... # ]
                $matrix_string .= "$1\n";
            } elsif ($line =~ /^\s*\d+/) {
                # line of the form: # # # ... #
                $matrix_string .= "$line\n";
            } else {
                next;
            }
            $line_count++;

            if ($line_count == 4) {
                my $pfm = TFBS::Matrix::PFM->new(
                    -matrixstring => $matrix_string,
                    -name         => $name,
                    -ID           => $ID
                );
                $matrix_set->add_Matrix($pfm);

                $line_count    = 0;
                $name          = '';
                $ID            = '';
                $matrix_string = '';
            }
        }
    }
    close(FH);

    return $matrix_set;
}
