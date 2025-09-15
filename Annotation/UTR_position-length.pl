#!/usr/bin/perl

### Length and position of the UTR

use strict;
use warnings;


open (SGT, "PATH/Ld1S_stg_filtered_final_renamed_transcriptONLY.txt");
open (MONO, "PATH/Ld1S_transdecoder_monocistron_coord.txt");

open (OUT, ">PATH/UTR.txt");

my %transcript;
print OUT "Transcript_id\tTranscript_length\t5'_length\t3'_start\t3'_length\n";
while (<SGT>) {
    chomp;
    my $line = $_;
	my @sgt_line = split /\s/, $line;
	my $sgt_id = $sgt_line[11];
    my $sgt_length = $sgt_line[4] - $sgt_line[3]; 
    $transcript{$sgt_id} = $sgt_length;
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
    #print STDERR "$mono_id\n";
    if (defined $transcript{$mono_id}) {
        my $UTR = $transcript{$mono_id} - $mono_end;
        print OUT "$mono_id\t$transcript{$mono_id}\t$mono_start\t$mono_end\t$UTR\n";
    }
}

close MONO;
close OUT;
exit;
