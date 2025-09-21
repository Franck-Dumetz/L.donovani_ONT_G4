#!/usr/bin/perl

### Rename stringtie output based on an another stringtie

use strict;
use warnings;

open (PRO, "/PATH/Ld1S_stg_filtered_final_renamed.gtf");
open (AMA, "/PATH/Ama_monocistron_SL_FiltCov.gtf");

open (OUT, ">/PATH/Ama_monocistron_SL_FiltCov_renamed.gtf");

my %gene;

while (<AMA>) {
    chomp;
    my $gene_split = $_;
	my @gff_line = split /\t/, $gene_split;
    my $ama_chr = $gff_line[0];
    my $ama_db = $gff_line[1];
    my $ama_qual = $gff_line[2];
	my $ama_start = $gff_line[3];
	my $ama_end = $gff_line[4];
    my $ama_1000 = $gff_line[5];
    my $ama_strand = $gff_line[6];
    my $ama_dot = $gff_line[7];
    my $toprint = "$ama_chr\t$ama_db\t$ama_qual\t$ama_start\t$ama_end\t$ama_1000\t$ama_strand\t$ama_dot";
    for (my $bp = $ama_start; $bp < $ama_end; $bp++) {					
			my $position = join("_",$ama_chr,$bp);
			$gene{$position} = $toprint;
            #print "$position\tama\n";
    }
}

close AMA;

while (<PRO>) {
    chomp;
    my $gene_split = $_;
	my @gff_line = split /\t/, $gene_split;
    my $pro_chr = $gff_line[0];
	my $pro_start = $gff_line[3];
	my $pro_end = $gff_line[4];
    my $pro_comment = $gff_line[8];
    my $midpoint = $pro_start + int(($pro_end - $pro_start)/2);
    my $position = join("_",$pro_chr,$midpoint);
    if (defined $gene{$position}) {print OUT "$gene{$position}\t$pro_comment\n"} 
        #print "$position\tpro\n";
}

close PRO;
close OUT;

exit;
