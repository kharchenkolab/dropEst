#!/usr/bin/perl

for (my $i = 2; $i <= 27; $i++) {
  my $cmd = "bsub -q short -W 1:00 -R 'rusage[mem=3000]' ../../cp/scripts/analyze_preseqr.r $i";
  print "$cmd\n";
  system($cmd);
}
