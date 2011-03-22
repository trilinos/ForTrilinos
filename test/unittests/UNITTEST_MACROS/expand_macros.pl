#!/usr/bin/perl

use strict;
use warnings;

use macroexp;
use misc;
use parser;

sub process_line {
  my ($input, $nindent) = @_;
  my $result = parser::expand_line($input);
  if (defined($result)) {
    my $wrapped = misc::wrap_line($result, 4, 100);
    my $output = misc::format_block($wrapped, $nindent);
    return $output;
  }
  return "";
}

macroexp::load_macros();

my $nindent = 0;

my $prev_line = "";
while (<>) {
  my $new_line = $_;
  if (($prev_line =~ /\&\s*$/) || ($new_line =~ /^\s*\&/)) {
    $prev_line =~ s/\s*\&?\s*$//;
    $new_line =~ s/^\s*\&?\s*//;
    $prev_line .= " " . $new_line;
  } else {
    my $output = process_line($prev_line, $nindent);
    print $output;
    $prev_line = $new_line;
  }
}

my $output = process_line($prev_line, $nindent);
print $output;

print "\n";
