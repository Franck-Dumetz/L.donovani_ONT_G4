#!/usr/bin/perl -w 
# creates a gff from Blast output

use strict;

open (IN, "/local/projects-t3/SerreDLab-3/fdumetz/Leishmania/Ld1S_genome/Flye_scaffold/blast_annotation/tblastn.out"); 

open (OUT, ">/local/projects-t3/SerreDLab-3/fdumetz/Leishmania/Ld1S_genome/Flye_scaffold/blast_annotation/Ld1S_blastout_p89_m98.gff");
#print OUT "GeneID\tChr\tStart\tEnd\tStrand\tDefinition\n" ### this line writes a header

my $last = "";


while (<IN>) {
	my $line = $_;
	chomp $line;
#	print "$line\n";
	
	my @blast = split "\t", $line;
	
	if ($blast[0] ne $last) {
		if ($blast[8] < $blast[9]) {print OUT "$blast[1]\tDumetzDB\tprotein_coding_gene\t$blast[8]\t$blast[9]\t.\t+\t.\t$blast[0]\n"}
			else {print OUT "$blast[1]\tDumetzDB\tprotein_coding_gene\t$blast[9]\t$blast[8]\t.\t-\t.\t$blast[0]\n"}
		}
	
	$last = $blast[0];
	
}

close IN;
close OUT;

exit;
