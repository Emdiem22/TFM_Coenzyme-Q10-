#!/usr/bin/perl
use strict;

# Annotations
my @GT;
my %go;
my %te;
my %re;
my %rr; 
my %rt;
my %kw;
my %gn;
my %ph;
my %ac;
my %ac2; 
my %tt;
my %ort;
my %ol2;

my $F = $ARGV[0] || ".";
my $II = $ARGV[1] || 0; 

# Gather Genename
open in, "./biomart_human_gn.tsv";
while (<in>) {
  chomp;

  next if /^Gene stable ID/;
  next if /\t$/;
  my ($id, $c) = split/\t/;
  push @{$gn{$id}}, $c;
  @{$kw{$id}} = ();
  @{$go{$id}} = ();
  @{$te{$id}} = ();
  @{$re{$id}} = ();
  @{$rt{$id}} = ();
  @{$ph{$id}} = ();
  @{$ort{$id}} = ();
  @{$ol2{$id}} = ();
  @{$tt{$c}} = ();
}
close in;

# Gather TRRUST
open in, "./trrust_rawdata.human.tsv";
while (<in>) {
  chomp;

  my ($a1, $a2, $a3) = split/\t/;
  next unless $a3 eq "Activation";
  push @{$tt{$a1}}, $a2;
}
close in;

# Gather SPROT
open in, "./biomart_human_sprot.tsv";
while (<in>) {
  chomp;

  next if /^Gene stable ID/;
  next if /\t$/;
  my ($id, $c) = split/\t/;
  push @{$ac{$id}}, $c;
  push @{$ac2{$c}}, $id;
}
close in;

# Gather GO
open in, "./biomart_human_go.tsv";
while (<in>) {
  chomp;

  next if /^Gene stable ID/;
  next if /\t$/;
  my ($id, $c) = split/\t/;
  push @{$go{$id}}, $c;
}
close in;

# Gather GO terms
open in, "./biomart_human_terms.tsv";
while (<in>) {
  chomp;

  next if /^Gene stable ID/;
  next if /\t$/;
  my ($id, $c) = split/\t/;
  push @{$te{$id}}, $c;
}
close in;

# Gather Reactome terms
open in, "./ReactomePathways.txt";
while (<in>) {
  chomp;

  next unless /Homo sapiens$/;
  my ($id, $c) = split/\t/;
  $rr{$id} = $c;
}
close in;

# Gather Reactome
open in, "./biomart_human_reactome.tsv";
while (<in>) {
  chomp;

  next if /^Gene stable ID/;
  next if /\t$/;
  my ($id, $c) = split/\t/;
  push @{$re{$id}}, $c;
  push @{$rt{$id}}, $rr{$c};
}
close in;

# Gather Phenotype
open in, "./biomart_human_phenotype.tsv";
while (<in>) {
  chomp;

  next if /^Gene stable ID/;
  next if /\t$/;
  my ($id, $c) = split/\t/;
  push @{$ph{$id}}, $c;
}
close in;

# Gather KW
open in, "./uniprot_human_kw.tsv";
while (<in>) {
  chomp;

  next if /\t$/;
  my ($id2, $c) = split/\t/;
  foreach my $id (@{$ac2{$id2}}) {
    push @{$kw{$id}}, $c;
  }
}
close in;

# Orthologous
open in, "./biomart_human_celegans.tsv";
while (<in>) {
  chomp;

  next unless /\t[01]$/;
  my ($a1, $a2, $a3) = split/\t/;
  push @{$ort{$a1}}, "$a2($a3)";
}
close in;

# Ortholist
open in, "./ortholist_master.tsv";
while (<in>) {
  chomp;

  my ($a2, $a1) = (split/\t/)[0, 4];
  push @{$ol2{$a1}}, "$a2";
}
close in;


# Results
if ($II == 1) {
  &INITIAL; 
  exit;
}

# Gather repeated genes
my %GEGE;
my (@G) = `ls $F/*_selgenes_direct.csv $F/*_selgenes_inverse.csv`;
foreach my $g (@G) {
  chomp $g;

  my ($gg) = split/_/, $g;
  $gg =~ s/[^\/]+\///g;

  # Gather correlation p-value
  my %pp = ();
  open in, $g;
  while (<in>) {
    chomp;
    my $ge = $_;
    push @{$GEGE{$ge}}, $gg;
  }
}

# Gather repeated genes II
my %GEGE2;
my (@G) = `ls $F/*_selgenes_direct.csv $F/*_selgenes_inverse.csv`;
foreach my $g (@G) {
  chomp $g;

  my ($gg) = split/_/, $g;
  $gg =~ s/[^\/]+\///g;

  # Gather correlation p-value
  my %pp = ();
  open in, $g;
  while (<in>) {
    chomp;
    my $ge = $_;
    push @{$GEGE2{$ge}}, $gg;
  }
}

# Gene targets
print "Seed\tSign\tGene_target\tEnsembl\tRepeated\tN_repeated\tRepeated2\tN_repeated2\tP-value\tOrtholist2\tOrthologs\tTRRUST(+)\tPhenotype\tGO\tGOterm\tKeyword\tReactome\tReactome_terms\n";
my (@G) = `ls $F/*_selgenes_direct_scores.csv $F/*_selgenes_inverse_scores.csv`;
foreach my $g (@G) {
  chomp $g;
  my $f = $g;

  # Gather correlation p-value
  my %pp = ();
  my $filescore = $f;
  $f =~ s/_scores//;
  open in, $filescore;
  while (<in>) {
    chomp;
    
    my $c = 6;
    $c = 7 if $filescore =~ /_inverse_/;
    my ($ensg, $p) = (split/,/)[0, $c];
    $pp{$ensg} = $p;
  }
  close in;

  #$g =~ s/results\///;
  #my ($g, $s) = split/_/, $g;
  #if ($s eq "0.7") { $s = "positive" } else { $s = "negative" }
  $g =~ s/$F\///;
  my ($g, $s) = (split/_/, $g)[0, 2];
  if ($s eq "direct") { $s = "positive" } else { $s = "negative" }
  
  # results
  open gene, $f;
  while (<gene>) {
    chomp;
    
    my $ge = $_;
    my $gg = join ";", @{$gn{$ge}};
    
    print "$g\t$s\t$gg\t$ge\t";
    print join ";", @{$GEGE{$ge}};
    print "\t";
    print $#{$GEGE{$ge}} + 1;
    print "\t";
    print join ";", @{$GEGE2{$ge}};
    print "\t";
    print $#{$GEGE2{$ge}} + 1;
    print "\t";
    print "$pp{$ge}\t";
    print join ";", @{$ol2{$ge}};
    print "\t";
    print join ";", @{$ort{$ge}};
    print "\t";
    print join ";", @{$tt{$gg}};
    print "\t";
    print join ";", @{$ph{$ge}};
    print "\t";
    print join ";", @{$go{$ge}};
    print "\t";
    print join ";", @{$te{$ge}};
    print "\t";
    print join ";", @{$kw{$ge}};
    print "\t";
    print join ";", @{$re{$ge}};
    print "\t";
    print join ";", @{$rt{$ge}};
    print "\n";
  }
  close gene;
}

# Subrutines
sub INITIAL () {
  print "Genename\tEnsembl\tPhenotype\tGO\tGOterm\tKeyword\tReactome\tReactome_terms\n";
  my (@G) = `ls experiments/*.tsv`;
  foreach my $g (@G) {
    chomp $g;

    open gene, $g;
    my $ge = <gene>;
    $ge = <gene>;
    ($ge) = split/\t/, $ge;
    close gene;
    $g =~ s/experiments\///;
    $g =~ s/\.tsv//;

    # Output
    print "$g\t$ge\t";
    print join ";", @{$ph{$ge}};
    print "\t";
    print join ";", @{$go{$ge}};
    print "\t";
    print join ";", @{$te{$ge}};
    print "\t";
    print join ";", @{$kw{$ge}};
    print "\t";
    print join ";", @{$re{$ge}};
    print "\t";
    print join ";", @{$rt{$ge}};
    print "\n";
  }
  return;
}
