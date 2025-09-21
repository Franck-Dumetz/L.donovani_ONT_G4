#!/usr/bin/perl

### counting 5' and 3' length isoforms per gene

use strict;
use warnings;

open (ISO, "/PATH/Pro_iso.txt");

open (OUT, ">/PATH/Pro_iso_count.txt");

my $last = "";

my $UTRL = 0;
my $UTRR = 0;
my $iso = 0;
my $strand = "";

my %endL;
my %endR;

my $scrap = <ISO>; # remove first line 
print OUT "Gene_ID\tNumber_of_isoform\tNumber_of_5'UTR_polymorphism\tNumber_of_3'UTR_polymorphism\n";

while (<ISO>) {
    chomp;
    my $line = $_;
	my @iso_line = split /\s/, $line;
    my $iso_ID = $iso_line[0];
    my $iso_strand = $iso_line[1];
    my $iso_number = $iso_line[2];
    my $iso_right = $iso_line[3];
    my $iso_left = $iso_line[4];
    if ($iso_ID eq $last) {
        $iso++;
        if (defined $endR{$iso_right}) {}
        else {
            $UTRR++;
            $endR{$iso_right} = 1;
        }
        if (defined $endL{$iso_left}) {}
        else {
            $UTRL++;
            $endL{$iso_left} = 1;
        }
    }

    else {
        if ($strand eq "+") {
            print OUT "$last\t$iso\t$UTRL\t$UTRR\n";
        }
        elsif ($strand eq "-") {
            print OUT "$last\t$iso\t$UTRR\t$UTRL\n";
        }
        $last = $iso_ID;
        $iso = 1;
        $UTRL = 1;
        $UTRR = 1;
        $strand = $iso_strand;
        undef %endR;
        undef %endL;
        $endR{$iso_right} = 1;
        $endL{$iso_left} = 1;
    }
}

close ISO;
exit;
