#!/usr/bin/env perl
use strict;
use warnings;

# Convert NCBI Genomes 'feature_table' (e.g.  KSHV
# ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/838/265/GCF_000838265.1_ViralProj14158/)
# to Circos text file format (ID START END NAME)

# Usage: $0 NCBI_Genomes_Feature_Table_File

my %labels;
while (<>) {
  my @fields = split (/\t/, $_);
  
  my $feature = $fields[0];
  next if ($feature =~ /^gene/);
  
  my $id = $fields[6];
  my $start = $fields[7];
  my $end = $fields[8];
  my $name = $fields[14] || $fields[10] || $fields[12];  # symbol or product_accession or related_accession or name (name is last because it can be too long for visualization purposes)

  if (!$name) {
    $name = $feature . ++$labels{$feature};
  }
  
  $name =~ s/ /_/g;

  print join(" ", $id, $start, $end, $name), $/;
}


