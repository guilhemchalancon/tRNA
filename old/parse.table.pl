#!/usr/bin/perl
my $outfile = "$ARGV[0]";

while(<>){
chomp;
my @f= split/\t/;
$f[0]=~s/.*(Y[A-Z][LR]\d{3}[WC](-[A-Z])?|Q\d{4}|T[A-Z]\([A-Z]{3}\)[A-Z]\d?)[_-]?.*/$1/g;
print "$f[0]\n";
push @output, join("\t",@f);
}

# Careful: to run from current directory where input file is
$outfile=~s/(.*)/cleaned.$1/;

open (OUT, ">> $outfile"); 

print OUT join("\n",@output);
close OUT
