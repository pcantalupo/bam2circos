#!/usr/bin/env perl

=head1 NAME

bam2circos.pl - Convert BAM file to Circos plot

=head1 SYNOPSIS

bam2circos.pl [optional] -k FILE -b FILE

 Required:
   -b|bam               BAM file
   -k|karyotype         Karyotype file

 Optional:
   -d|depth FILE        Depth file
   -s|slidingpid FILE   Sliding window percent identity file
   -w|window INT        Window size (default: 50) for sliding window pid calculation
   -o|overlap INT       Overlap size (default: 10) for sliding window pid calculation
   -a|annotations FILE  Annotation file
   --keep               Keep the generated circos configuration and text files
   -t|test              Generate circos configuration files but don't execute circos
   -h|help             	Full documentation

=head1 OPTIONS

=over 4


=item B<-b|bam>

    BAM file containing aligments to one or more reference sequences. 
    Reference IDs in RNAME column must match those found in the karyotype,
    annotation, depth and slidingpid files.

=item B<-k|karyotype>

    The karyotype file needs to conform to the Circos karyotype file format:
    http://circos.ca/tutorials/lessons/ideograms/karyotypes/.  CAREFUL:
    There must only be 1 space between the fields.  For most viruses, the
    file will contain one row to define the full length of the viral genome. 
    You can add rows to define segmented viruses.  The 3rd column (ID) must
    match the RNAME column (reference column) of the BAM file.  I found that
    a | (vertical bar) is not allowed in the ID.

    Example: Bunyaviridae have 3 segments
    chr - NC_014397 L 0 6404 chr1
    chr - NC_014396 M 0 3885 chr2
    chr - NC_014395 S 0 1690 chr3

=item B<-a|annot>
    
    The annotation file needs to conform to the Circos text track file
    format: http://circos.ca/documentation/tutorials/2d_tracks/text_1/
    
    Example (two genes named FOO and BAR)
    NC_001669 450 650 FOO
    NC_001669 600 1200 BAR

=item B<-d|depth>

    The file needs to conform to the Circos data file format:
    http://circos.ca/documentation/tutorials/2d_tracks/histograms/
    
    Example (base 1 has depth of 5 and base 2 depth of 10)
    NC_001669 0 1 5
    NC_001669 1 2 10

=item B<-s|slidingpid>

    The file needs to conform to the Circos data file format above except
    that value of the fourth column is the sliding window percent identity
    
    Example (50bp sliding window where 1st 50bp is 88.5%, 2nd 50bp is 93.4%)
    NC_001669 0 50 88.5
    NC_001669 51 100 93.4

=item B<-help>

    Print full documentation and exits.

=back

=head1 DESCRIPTION

B<This program> will do blah blah blah.

=head1 REQUIREMENTS

=over 4

=item circos

=item samtools

=item awk

=back

=cut

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use List::Util qw/sum/;

my %conf = (circos   => "virus.circos.conf",
            depth    => "virus.depth.conf",
            ideogram => "virus.ideogram.conf",
            ticks    => "virus.ticks.conf",
            slidingpid => "virus.slidingpid.conf",
            );
my ($bam, $depth, $pid, $annot, $karyo);
my ($window, $overlap) = (50, 10);
my ($help, $man, $keep, $test) = (0, 0, 0, 0);
GetOptions ('b|bam=s'        => \$bam,
            'd|depth=s'      => \$depth,
            's|slidingpid=s' => \$pid,
            'w|window=i'     => \$window,
            'o|overlap=i'    => \$overlap,
            'a|annotations=s'=> \$annot,
            'k|karyotype=s'  => \$karyo,
            'h|help'         => \$help,
            't|test'         => \$test,
            'keep'           => \$keep,
            ) || pod2usage (-verbose => 0);

pod2usage(-verbose => 2) if ($help);
pod2usage(-verbose => 1) if (!$bam || !$karyo);

initialize (\%conf, $karyo, $annot, $pid, $depth, $bam, $window, $overlap);

run($conf{circos}) if ($test == 0);

cleanup (\%conf, $keep);
   
exit 0;




################################################################################

sub initialize {
  my ($conf, $karyotype, $annot, $pid, $depth, $bam, $window, $overlap) = @_;

  my $annottemp = "temp.annot.txt";
  $annot //= $annottemp;
  if ($annot eq $annottemp) {
    create_annot_file ($annottemp, $karyotype);
  }
  
  my $depthtemp = "temp.depth.txt";
  $depth //= $depthtemp;
  if ($depth eq $depthtemp) {
    create_depth_file ($depthtemp, $bam, $karyotype);
  }

  my $slidingpidtemp = "temp.slidingpid.txt";
  $pid //= $slidingpidtemp;
  if ($pid eq $slidingpidtemp) {
    create_slidingpid_file($slidingpidtemp, $window, $overlap, $bam, $karyotype);
  }
   
  slidingpid_conf($conf->{slidingpid}, $pid);  
  circos_conf($conf, $karyotype, $annot, $pid);
  depth_conf($conf->{depth}, $depth);
  ticks_conf($conf->{ticks});
  ideogram_conf($conf->{ideogram});

}

sub run {
  my ($circosconf) = @_;
  my $command = "circos -conf $circosconf";
  print "Running circos: $command\n";
  `$command > circos.log 2>&1`;
}

sub cleanup {
  my ($conf, $keep) = @_;
  
  if ($keep) {
    print STDERR "Not deleting generated circos conf files\n";
  }
  else {
    foreach my $conffile (%$conf) {
      unlink ($conffile);
    }
    
    unlink ("temp.annot.txt", "temp.depth.txt", "temp.slidingpid.txt");
  } 
}

sub create_slidingpid_file {
  my ($slidingpidfile, $window, $overlap, $bam, $karyotype, $random) = @_;

  my %vlen;
  open (my $kin, "<", $karyotype) or die ("Can't open $karyotype: $!\n");
  while (my $line = <$kin>) {
    my (undef, undef, $id, $name, $start, $end) = split (/ /, $line);
    $vlen{$id} = $end;
  }
  close $kin;
   
  my $text = "";
  if ($random) {
    # calculate a random window percent identity
#    my $w = 50;
#    my $iter = $length / $w;
#    my $s = 0;
#    my $e = $w - 1;
#    for (my $i = 0; $i < $iter; $i++) {
#      my $line_s = $s + $w * $i;
#      my $line_e = $e + $w * $i;
#      my $x = 75 + int(rand(100 - 75));
#      $text .= "$id $line_s $line_e $x\n";
#    }
  }
  else {
    # calculate sliding window percent identity
    open (my $bam_in, "samtools view $bam | ") or die ("Can't open pipe from samtools view of a bam file: $!\n");
    my $bases;   # 0th index in array will not be used.
    while (<$bam_in>) {
      my ($q, $f, $r, $p, $m, $c, undef, undef, undef, $seq, undef, @rest) = split(/\t/, $_);
      my ($nm) = $_ =~ /NM:i:(\d+)/;
      my $rlen = length($seq);
      my $avg_nm_read = ($rlen - $nm) / $rlen;
      my $end = $rlen + $p - 1;
      for (; $p <= $end; $p++) {
        push (@{$bases->{$r}->[$p]}, $avg_nm_read);
      }    
    }
    
    # calculate the mean of the average NM / base values for each base
    foreach my $id (keys %$bases) {
      for (my $i = 1; $i <= $#{$bases->{$id}}; $i++) {
        my $mean = mean( @{ $bases->{$id}->[$i] } );
        #print "$i $mean",$/;
        $bases->{$id}->[$i] = $mean;
      }
    }

    foreach my $id (keys %$bases) {
      my $length = $vlen{$id};
      for (my $i = 1; $i <= $length - $window - 1; $i += ($window - $overlap) ) {
        my $end = $i + $window - 1;
        $text .= "$id $i $end " . 100 * mean(@{ $bases->{$id} }[$i .. $end]) . "\n";
      }
    }
  }    

  open (my $sout, ">", $slidingpidfile) or die ("Can't open $slidingpidfile for writing: $!\n");
  print $sout $text;
  close $sout;    
} 

 

sub create_annot_file {
  my ($annotfile, $karyotype) = @_;
  
  # get ID, start and end for each chromosome and create an annotation for each chromosome
  open (my $kin, "<", $karyotype) or die ("Can't open $karyotype: $!\n");
  open (my $aout, ">", $annotfile) or die ("Can't open $annotfile for writing: $!\n");
  while (my $line = <$kin>) {
    my (undef, undef, $id, $name, $start, $end) = split (/ /, $line);
    print $aout "$id $start $end $id\n";
  }
  close $kin;
  close $aout;
}


sub create_depth_file {
  my ($depthfile, $bam, $karyotype) = @_;

  my %vid;
  open (my $kin, "<", $karyotype) or die ("Can't open $karyotype: $!\n");
  while (my $line = <$kin>) {
    my (undef, undef, $id) = split (/ /, $line);
    $vid{$id}++;
  }
  close $kin;

  my $awk = q{awk -F "\t" '{print $1,$2-1,$2,$3}'};
  my $command = join("", "samtools depth $bam | ", $awk);
  print "Running samtools depth: $command\n";
  my @output = `$command`;

  my $text = "";
  foreach (@output) {
    my ($id) = split;
    $text .= $_ if (exists $vid{$id});
  }

  open (my $tout, ">", $depthfile) or die ("Can't open $depthfile for writing: $!\n");
  print $tout $text;
  close $tout;
}


sub circos_conf {
  my ($conf, $karyotype, $annot, $pid) = @_;

  my $slidingpid_text = "";
  if ($pid) {
    $slidingpid_text = <<HEREDOC;
###################
### SLIDING PID ###
<plot>
<<include ./$conf->{slidingpid}>>
</plot>
HEREDOC
  }

  my $text = <<HEREDOC;
<<include etc/colors_fonts_patterns.conf>>
<<include ./$conf->{ticks}>>
<<include ./$conf->{ideogram}>>

karyotype = ./$karyotype

<image>
<<include etc/image.conf>>
</image>

<<include etc/housekeeping.conf>>

#chromosomes_units           = 1000
#chromosomes_display_default = yes

<plots>
############
### ORFS ###
<plot>
type      = tile
file      = $annot
r1        = 1r
r0        = 0.85r
layers    = 3
margin    = 0.01u
thickness = 50
padding   = 4
layers_overflow  = hide
orientation      = out
stroke_thickness = 1
stroke_color     = grey
color            = vlorange

#<rules>
#<rule>
#condition = var(size) < 1kb
#color     = lgrey
#</rule>
#</rules>

</plot>

<plot>
type = text
file = $annot
r1 = 1r
r0 = 0.85r
show_links     = yes
link_dims      = 2p,2p,2p,2p,2p
link_thickness = 4p
link_color     = red
#label_rotate = no
label_parallel = yes
label_size = 30p
label_font = bold
</plot>


################
### DEPTH ###
<plot>
<<include ./$conf->{depth}>>
# radius 70 to 80
</plot>

$slidingpid_text

</plots>
HEREDOC

  writefile ($text, $conf->{circos});
  return $text;
}

sub ideogram_conf {
  my ($conffile) = @_;
   
  my $text = <<HEREDOC;
### IDEOGRAM ###
<ideogram>

<spacing>
default = 0.01r
break   = 0.5r
</spacing>

# ideogram position, fill and outline
radius    = 0.9r
thickness = 50p  # default = 20p
fill      = yes
stroke_color = dgrey
stroke_thickness = 2p

# minimum definition for ideogram labels
show_label = yes
label_font = default
label_radius = 1r + 75p
label_size = 50
label_parallel = yes

</ideogram>
HEREDOC

  writefile ($text, $conffile);
  return $text;
}


sub ticks_conf {
  my ($conffile) = @_;
  
  my $text = <<HEREDOC;
### TICKS ###
show_ticks = yes
show_tick_labels = yes

<ticks>
radius           = 1r
color            = black
thickness        = 2p

# the tick label is derived by multiplying the tick position
# by 'multiplier' and casting it in 'format':
#
# sprintf(format,position*multiplier)
multiplier       = 1

# %d   - integer
# %f   - float
# %.1f - float with one decimal
# %.2f - float with two decimals
#
# for other formats, see http://perldoc.perl.org/functions/sprintf.html
format           = %d

<tick>
spacing        = 100u
size           = 10p
</tick>

<tick>
spacing        = 1000u
size           = 25p
show_label     = yes
label_size     = 36p
label_offset   = 10p
format         = %d
</tick>

</ticks>
HEREDOC
  
  writefile($text, $conffile);
  return $text;
}


sub depth_conf {
  my ($conffile, $depthfile) = @_;

  my $text = <<HEREDOC;
### COVERAGE ###
type = line
file = $depthfile
r1 = 0.80r
r0 = 0.70r
thickness = 0p
fill_color=red
max_gap=1b
HEREDOC

  writefile ($text, $conffile);
  return $text;
}


sub slidingpid_conf {
  my ($conffile, $slidingpidfile) = @_;
  
  my $text = <<HEREDOC;
### SLIDING PID ###
type      = line
thickness = 10   # 2

max_gap = 1u
file    = $slidingpidfile
color   = chr15
#fill_color = chr17
min     = 0
max     = 100
r1      = 0.65r
r0      = 0.55r

<axes>
<axis>
color     = lgrey
thickness = 1
spacing   = 0.1r
</axis>
</axes>
HEREDOC
  
  writefile($text, $conffile);
  return $text;
}




sub writefile {
  my ($text, $outfile) = @_;
  open (my $out, ">", $outfile) or die ("Can't open $outfile for writing: $!\n");
  print $out $text;
  print STDERR "Wrote $outfile file\n";
  close ($out);  
}



sub mean {
  if (@_) {
    return sum(@_)/@_;
  }
  else {
    return 0;
  }
}
