package misc;

use strict;
use warnings;

sub format_block {
  my ($block, $nindent) = @_;

  my $ind = ' ' x $nindent;
  $block =~ s/^/${ind}/mg;

  return $block;
}

sub format_group {
  my ($lines, $nindent) = @_;

  my $output = "";
  foreach my $line (@$lines) {
    $output .= ' ' x $nindent . $line . "\n";
  }

  return $output;
}

1;
