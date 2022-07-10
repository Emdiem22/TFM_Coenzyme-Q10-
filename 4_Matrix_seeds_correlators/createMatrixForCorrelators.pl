#!/usr/bin/perl
use strict;

my @s;
my %n;
my %corr;
my $F = $ARGV[0] || "./common_targets.tsv";

open in, "./seeds.tsv";
while (<in>) {
  chomp;

  my ($id) = split/\t/;
  push @s, $id;
}
close in;

open in, $F;
while (<in>) {
  chomp;

  my ($i1, $i2) = split/\t/;
  $n{$i1}{$i2}++;
  $corr{$i2} = 1;
}
close in;

print "Identifier\t";
print join "\t", sort @s;
print "\n";
foreach my $s1 (sort keys %corr) {
  print $s1;
  foreach my $s2 (sort @s) {
    my $x = 0;
    $x = $n{$s2}{$s1} if $n{$s2}{$s1};
    print "\t$x";
  }
  print "\n";
}

exit;

