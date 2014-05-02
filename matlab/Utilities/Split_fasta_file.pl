#!/usr/bin/perl -w

#####################################################################################################################
# Written by Eran Elhaik
# Ver 1.10
#####################################################################################################################
# The program accept input of a Fasta file, that contains many files
# and create a directory Fasta with all the fasta files inside
# Activation: perl divide_file.pl filename output_dir
# perl Split_fasta_file.pl /amber1/archive/chaklab/users/eelhaik/Genomes/Robin/original/Lepmu1_AssemblyScaffolds.fa /amber1/archive/chaklab/users/eelhaik/Genomes/Robin/scaffolds_with_ns/
#####################################################################################################################
# Ver 1.10: I added an option for output_dir
#####################################################################################################################


# input file
my $fileName = $ARGV[0];
my $output_dir = $ARGV[1];

# input file handle
my $IN = *IN;

my @data_file = ();
my $line = '';
my $temp_seq = '';
my $length;
my $fasta_extension='.fa';

# open file
open (IN, "<".$fileName ) || die "can't open file: $!";
@data_file = <IN>;
close (IN) || die "can't close file: $!";

my $fasta_header = shift @data_file;
# my $s =  substr($fasta_header,2);

($fasta_header,$length) = split(/\s+/, substr($fasta_header,1),2);
print $fasta_header;

foreach $line (@data_file){

	# discard blank line
	if($line =~ /^\s*$/){
		next;
	chomp $line;
	# new gene
	}elsif($line =~ /^>/){
		# write the sequence to file 
		open (OUT,">".$output_dir.$fasta_header.$fasta_extension ) || die "can't open file: $!";
                print OUT ">".$fasta_header."\n".$temp_seq; #Print the header
                close (OUT) || die "can't open file: $!";
		
		($fasta_header,$length) = split(/\s+/, substr($line,1),2);
		$temp_seq = '';		
		next;
	}else{
		# concat the sequence
		$temp_seq .= $line;
 	}
}

open (OUT,">".$output_dir.$fasta_header.$fasta_extension ) || die "can't open file: $!";
 print OUT ">".$fasta_header."\n".$temp_seq; #Print the header
close (OUT) || die "can't open file: $!";
