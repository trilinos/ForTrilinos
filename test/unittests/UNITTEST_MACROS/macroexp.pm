package macroexp;

#use strict;
use warnings;

use parser;

sub strfy { return "\"" . $_[0] . "\""; }
sub printany { return "print *," . $_[0]; }
sub printlit { return "print *," . strfy($_[0]); }
sub echo { return printlit($_[0]) . "\n" . doit($_[0]); }

sub printany_list {
  my @terms = @_;
  my $str = "print *";
  foreach my $t (@terms) {
    $str .= "," . $t;
  }
  return $str;
}

sub printlit_list {
  my @terms = @_;
  my $str = "print *";
  foreach my $t (@terms) {
    $str .= "," . strfy($t);
  }
  return $str;
}

sub srnd { return "(" . $_[0] . ")"; }
sub doit { return "$_[0]"; }
sub assert_fail { return $_[0] . " = .FALSE.\n" . printlit_list("Assertion failed on line ", $_[1]); }

sub print_test {
  return printany_list(strfy("TEST: $_[0] = "), $_[0],
                    strfy(" ?$_[1]? $_[2] = "), $_[2]); }

sub binary_test {
  # pass string with comma-separated args, comparison op string, fail operator
  my @args = parser::break_args($_[0]);
  my $ret = print_test($args[0], $_[1], $args[1]);
  $ret .= "\nif (" . srnd($args[0]) . " " . $_[2] . " " . srnd($args[1]) . ") then";
  $ret .= "\n" . assert_fail("success", "?");
  $ret .= "\nendif";
  return $ret;
}

sub test_equality { return binary_test($_[0], "==", ".NE."); }
sub test_inequality { return binary_test($_[0], "!=", ".EQ."); }
sub test_equiv { return binary_test($_[0], ".eqv.", ".NEQV."); }
sub test_lessequal { return binary_test($_[0], "<=", ".GT."); }


my %macro_map;

sub load_macros {
  $macro_map{'ECHO'} = "macroexp::echo";
  $macro_map{'STRINGIFY'} = "macroexp::strfy";
  $macro_map{'PRINTANY'} = "macroexp::printany";
  $macro_map{'PRINTLIT'} = "macroexp::printlit";
  $macro_map{'TEST_EQUALITY'} = "macroexp::test_equality";
  $macro_map{'TEST_INEQUALITY'} = "macroexp::test_inequality";
  $macro_map{'TEST_EQUIV'} = "macroexp::test_equiv";
  $macro_map{'TEST_LESSEQUAL'} = "macroexp::test_lessequal";
}

sub expand_macro {
  my ($key, $args) = @_;

  if (defined($macro_map{$key})) {
    my $fn = $macro_map{$key};
    return &{$fn}($args);
  }

  return undef;
}

1;
