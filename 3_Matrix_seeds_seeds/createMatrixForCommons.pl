#!/usr/bin/perl
use strict;

my @s;
my %n;
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
  push @{$n{$i1}}, $i2;
}
close in;

print "Identifier\t";
print join "\t", sort @s;
print "\n";
foreach my $s1 (sort @s) {
  print $s1;
  foreach my $s2 (sort @s) {
    my $nn = 0;
    my %rep;

    for my $s (@{$n{$s1}}, @{$n{$s2}}) {
      $nn++ if ($rep{$s});
      $rep{$s} = 1;
    }
    print "\t$nn";
  }
  print "\n";
}

exit;
