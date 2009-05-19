// This C++ file is just here to trick autotools into using the C++ linker
// when linking the *.f03 programs.  Using the C++ Linker is less involved than
// using the Fortran linker because the number libraries that must be explicitly
// linked is lower.  This idea was copied from CTrilinos.
// 
