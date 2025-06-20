#!/usr/bin/perl

### finding isoforms that are at least 10% of the overall coverage OR a minimum of 5 reads per CDS.
### The isoform that are detected will have at least 100 nucleotides difference in 3' or in 5' UTR length (or both)

use strict;
use warnings;

open (SGT, "/local/projects-t3/SerreDLab-3/fdumetz/Leishmania/Ld_annotation/Step2_Annotation_curation/Step5_Transdecoder/Stringtie_renamed/Ld1S_stg_filtered_final_renamed_transcriptONLY.txt");
open (MONO, "/local/projects-t3/SerreDLab-3/fdumetz/Leishmania/Ld_annotation/Step2_Annotation_curation/Step5_Transdecoder/Stringtie_renamed/Ld1S_transdecoder_monocistron_coord.txt");
open (BAM, "samtools view /local/projects-t3/SerreDLab-3/fdumetz/Leishmania/Ld_ONT/Annotation/Amastigotes/Ld1S_AxAma_ONT_induro.bam |");

open (OUT, ">/local/projects-t3/SerreDLab-3/fdumetz/Leishmania/Ld_annotation/Step2_Annotation_curation/Step5_Transdecoder/Stringtie_renamed/Ama_isoform_wiggle100.txt");

my %transcript;

while (<SGT>) {
    chomp;
    my $line = $_;
	my @sgt_line = split /\s/, $line;
    my $chr = $sgt_line[0];
	my $sgt_id = $sgt_line[11];
    my $sgt_start = $sgt_line[3];
    my $sgt_end = $sgt_line[4]; 
    $transcript{$sgt_id}[0] = $chr;
    $transcript{$sgt_id}[1] = $sgt_start;
    $transcript{$sgt_id}[2] = $sgt_end;
    #print STDERR "$sgt_id\n";
}

close SGT;


while (<MONO>) {
    chomp;
    my $line = $_;
	my @mono_line = split /\s/, $line;
	my $mono_id = $mono_line[0];
    my $mono_start = $mono_line[10];
    my $mono_end = $mono_line[11];   
    if (defined $transcript{$mono_id}) {
        $transcript{$mono_id}[3] = $mono_start + $transcript{$mono_id}[1];
        $transcript{$mono_id}[4] = $mono_end + $transcript{$mono_id}[1];
    }
}

close MONO;

my %isoform;
my %coverage;

while (<BAM>){
    chomp;
    my $line = $_;
    my @bam_line = split /\t/, $line;
    my $bam_chr = $bam_line[2];
    my $bam_start = $bam_line[3];
    my $size = length($bam_line[9]);
    my $bam_end = ($bam_start); #### remove short clipping check with Janne
                    #Define softclipping numbers from cigar string:
    my $cigar = $bam_line[5];
    my $softclip_start = 0;
    my $softclip_end = 0;
    if ($cigar =~ /^(\d+)S/) {
        $softclip_start = $1;
    }
    if ($cigar =~ /(\d+)S$/) {
        $softclip_end = $1;
    }
    my $size_clipped = $size - $softclip_start - $softclip_end;
    $bam_end = $bam_end + $size_clipped;
    #Check if the correct numbers are extracted:
    #print "$softclip_start\t$softclip_end\t$cigar\t$size\t$size_clipped\n";
    #sleep (1);
    my $s = $bam_start;  ### wiggle room of 10 nucleotides from the start of each read
    my $e = $bam_end;
    #print STDERR "$bam_start test\t$size test\t$bam_end\n";
    #print STDERR "$s\n";
    #sleep (1);
    #print STDERR "$bam_chr\t$s\t$e\n";
    #sleep (1);
    foreach my $t (keys %transcript) {
        if (defined $transcript{$t}[3]) {
           if (($bam_chr eq $transcript{$t}[0])&&($transcript{$t}[1] <= $s)&&($s<= $transcript{$t}[3])&&($transcript{$t}[4] <= $e)&&($e <= $transcript{$t}[2])) {
            $s = int($s/100);  ### wiggle room of 10 nucleotides from the start of each read
            $e = int($e/100);
            $isoform{"$t\_$s\_$e"}++; 
            $coverage{$t}++;
            #print STDERR "$transcript{$t}[0]\t$transcript{$t}[1]\t$transcript{$t}[2]\t$transcript{$t}[3]\t$transcript{$t}[4]\n";
            }
        }
    } 
}
close MONO;
close BAM;
foreach my $i (keys %isoform) {
    if ($isoform{$i}>= 5) {
        my @name = split "_", $i;
        if ($isoform{$i} >= (0.1 * $coverage{$name[0]})) {
            print OUT "$i\t$isoform{$i}\t$coverage{$name[0]}\n";
        }
    }
}

exit;
