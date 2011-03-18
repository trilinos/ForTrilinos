#!/usr/bin/perl

use strict;
use warnings;

use macroexp;
use misc;
use parser;

macroexp::load_macros();

my $nindent = 0;

while (<>) {
  my $input = $_;
  my $result = parser::expand_line($input);
  if (defined($result)) {
    my $output = misc::format_block($result, $nindent);
    print $output;
  }
}

print "\n";
