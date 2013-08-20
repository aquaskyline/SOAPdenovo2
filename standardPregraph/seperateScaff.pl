#!/usr/bin/perl -w
use strict;
use Getopt::Long;

=head1 Description

Program: seperateScaff.pl

Author: BGI-Shenzhen

Version: 1.0

Contact: soap@genomics.cn

This script is used to seperate singletons from scaffolds in *.scafSeq file 
output by SOAPdenovo. There are two new output files as follow:

1) *.scafSeq.scaffold, containing scaffolds only.

2) *.scafSeq.singleton, containing singletons only.

User can set the minimun length of sequence to output.

=head1 Usage

perl seperateScaff.pl [options] >seperateScaff.log

	--scafSeq	*.scafSeq file output by SOAPdenovo, required

	--scafLen	minumum length of scaffold to output, default 0

	--sigtLen	minumum length of singleton to output, default 0

	-h		output help information

=cut

die `pod2text $0` if (@ARGV < 1);
my ($scafSeq, $help);
my ($scafLen, $sigtLen) = (0, 0);

GetOptions (
	"scafSeq:s"=>\$scafSeq,
	"scafLen:i"=>\$scafLen,
	"sigtLen:i"=>\$sigtLen,
	"h"=>\$help
);

die `pod2text $0` if ($help);
if (! -f $scafSeq) {
	print "\n[Error] $scafSeq does not exist.\n\n";
	die `pod2text $0`;
}

if ($scafLen < 0) {
	print "\n[Warning] Set minimum length of scaffold is negative, program will use 0 instead.\n\n";
	$scafLen = 0;
}

if ($sigtLen < 0) {
	print "\n[Warning] Set minimum length of singleton is negative, program will use 0 instead.\n\n";
	$sigtLen = 0;
}

print "Command:\n  perl $0 --scafSeq $scafSeq --scafLen $scafLen --sigtLen $sigtLen\n\n";

my ($gotScaf, $gotSigt, $id, $seq, $len) = (0, 0, "", "", 0);
my ($totalScafNum, $totalScafLen, $totalSigtNum, $totalSigtLen) = (0, 0, 0, 0);
my ($outputScafNum, $outputScafLen, $outputSigtNum, $outputSigtLen) = (0, 0, 0, 0);
open SCAF, ">$scafSeq.scaffold" or die "Can't open output file: $scafSeq.scaffold\n";
open SIGT, ">$scafSeq.singleton" or die "Can't open output file: $scafSeq.singleton\n";
open IN, "$scafSeq" or die "Can't open input file: $scafSeq\n";
while (<IN>) {
	if (/^>/) {
		handleOneSeq();

		$id = $_;
		$seq = "";
		$len = 0;
		if ($_ =~ /^>scaff/) {
			$gotScaf = 1;
			$gotSigt = 0;
		}
		elsif ($_ =~ /^>C/) {
			$gotSigt = 1;
			$gotScaf = 0;
		}
	}
	else {
		$seq .= $_;
		$len += length ($_) - 1;
	}
}
handleOneSeq();
close IN;
close SCAF;
close SIGT;

print "Total scaffold number: $totalScafNum\nTotal scaffold length: $totalScafLen\nNumber of scaffold no shorter than $scafLen: $outputScafNum\nLength of scaffold no shorter than $scafLen: $outputScafLen\n\n";
print "Total singleton number: $totalSigtNum\nTotal singleton length: $totalSigtLen\nNumber of singleton no shorter than $sigtLen: $outputSigtNum\nLength of singleton no shorter than $sigtLen: $outputSigtLen\n";

#################################
########## Sub Routine ##########
#################################

sub handleOneSeq {
	if ($gotScaf == 1) {
		$totalScafNum++;
		$totalScafLen += $len;
		if ($len >= $scafLen) {
			$outputScafNum++;
			$outputScafLen += $len;
			print SCAF "$id$seq";
		}
	}
	elsif ($gotSigt == 1) {
		$totalSigtNum++;
		$totalSigtLen += $len;
		if ($len >= $sigtLen) {
			$outputSigtNum++;
			$outputSigtLen += $len;
			print SIGT "$id$seq";
		}
	}
	1;
}
