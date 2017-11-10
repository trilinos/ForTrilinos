%module forbelos

%include <std_string.i>

%include "ForTrilinosBelos_config.hpp"

%ignore Belos::toString;

#define BELOS_DEPRECATED
%include "Belos_Types.i"
