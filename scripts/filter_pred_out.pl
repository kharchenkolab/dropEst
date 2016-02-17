#!/usr/bin/perl
my @files = `ls *.rds.pred_out`;

open my $file, ">all.pred_out";
print $file "BamCount\tReads\tUmigs\n";
for my $name (@files) {
  open my $in_file, $name;
  my $line = <$in_file>;
  close $in_file;

  my ($size, $id) = $name =~ /.*_(\d+)_(\d+)\.rds\.pred_out/;
  my ($reads, $umigs) = $line =~ /(\d+): (\d+)/;
  
  print $file "$size\t$reads\t$umigs\n"; 
}

close $file;
