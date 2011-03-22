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

sub break_terms {
  my ($argstr) = @_;

  my $string = $argstr;
  my @params;
  while ($string ne "") {
    if ($string =~ /^([^"']*?),/) {
      push(@params, $1);
      $string =~ s/^\Q$1\E\s*,?\s*//;
    } else {
      my $tmp = $string;
      my $dq = ($tmp =~ tr/"/X/);
      my $sq = ($tmp =~ tr/'/X/);
      if (($dq >= 2) || ($sq >= 2)) {



        my ($ext, $new_string, $pre) = extract_quotelike($string,'[^"\']*');
        $string = $new_string;
        push(@params, "$pre$ext");
        $string =~ s/^\s*,\s*//;
      } else {
        push(@params, "$string");
        $string = "";
      }
    }
  }
  return @params;
}

sub wrap_line {
  my ($line, $nindent, $max_len) = @_;

  my $ind = ' ' x $nindent;

  my $out_block = "";
  my @lines = split(/(?<=\n)/, $line);
  foreach my $newline (@lines) {
    my @wrap = split(/(?<=\s)/, $newline);

    my $out_line = "";
    foreach my $w (@wrap) {
      if ((length($out_line) > 0) and ((length($out_line) + length($w)) > $max_len)) {
        $out_block .= $out_line;
        $out_line = "&\n" . $ind;
      }
      $out_line .= $w;
    }
    $out_block .= $out_line;
  }

  return $out_block;
}

1;
