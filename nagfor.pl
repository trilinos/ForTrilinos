#!/usr/bin/perl 

# 
# Wrapper for compilers having troubles with preprocessor.
#
# Assumes the compiler is in your PATH
#
# Works for only 1 input file at a time.
# 
# Written for a version of NAG that could not preprocess
# lines longer than 132 chars.
# 
$compiler="nagfor";

#
# Scan through the argument list separating:
#   1. input file needing preprocessing:
#      .Fnn .ffnn  where nn are two digits (covers
#              f77 f90 f95 f03)
#      .F   .ff   
#   2. input file not (necessarily) needing preprocessing:
#      .fnn .f   SHOULD WE ADD .FOR?? 
#
#   3. Options affecting the preprocessor:
#      -Idir -Dmacro
#      
#   4. Output clause: -o outfile or -ooutfile
#   
#   5. Option to force preprocessing even for files in 2.
#
#   6. -c to have compile but not link (affects output file name)
#   
#   7. Any other option, irrelevant to the preprocessor
#      and intended for the base fortran compiler
#
#   8. Script-specific options: currently only -debug
#      They must be removed from @ARGV, lest they get
#      passed to the base compiler. 
#

MAIN: for ($i=0; (($i<$#ARGV), $_=$ARGV[$i]) ; $i++) {
#    print "debug: $_ \n";
    SWITCH: {
	if (/.*\.F[0-9][0-9]$|.*\.ff[0-9][0-9]$|.*\.F$|.*\.ff$/) {
		$infile=$_;
		$preproc=1;
		last SWITCH;	    }
	    
	if (/.*\.f[0-9][0-9]$|.*\.f$/) {
		$infile=$_;
		last SWITCH;	   }
	if (/^-I|^-D/) {
		push (@ccopt,$_);
		last SWITCH;   }
	    
	if (/^-debug$/) { 
		$debug=1;
		splice(@ARGV,$i,1); 
		last SWITCH;   }
	
	    
	if (/^-o$/) {
		$outfile=$ARGV[$i+1];
		$i++;
		last SWITCH;   }

	if (/^-o(.*)$/) {
		$outfile=$1;
		last SWITCH;   }
		
	if (/^-fpp$/) {
		$preproc=1;
		last SWITCH;   }
	if (/^-c$/) {
		$componly=1;
		push (@fopt,$_);
		last SWITCH;   }

        push (@fopt,$_);
		
	}
}

#
# Ok, so: we have to preprocess if 
#  1. the flag is set and
#  2. there is an input file
#   remember we might just have been called to link
#   an existing .o
#

if (($infile)&&(-e $infile)&&($preproc)) {
    #
    #
    # Build a name for the temporary file
    # holding preprocessed source. Should be the same basename
    # with single lowercase "f" plus digits, if any.
    # Take out prefix to write in current directory, just
    # in case source dir is not writable, and make sure not
    # to conflict with existing files.
    # 
    #
    $pfile=$infile;
    $pfile =~ s/.*\///;
    $pfile =~ s/\.[Ff][Ff]*([0-9]*)$/\.f\1/;
    while ( -e $pfile ) {
	$pfile="XtmpX_$pfile";
    }
    
    #
    # Do the preprocessing; adjust return code if needed. 
    #
    push (@cpp,"gcc","-traditional-cpp","-undef","-x", "c","-E","-P");
    push (@cpp,@ccopt);
    push (@cpp,"-o",$pfile,$infile);
    print "@cpp \n";
    $rc=system(@cpp);
    ($rc == 0) or exit $rc>>8;

    
    #
    # Now for compilation: if outfile was not forced, and if
    # we are supposed to produce a .o file, we have to
    # build the correct name.
    # Need to push the -I options, if any; they are often used
    # to specify search path for module files. 
    #
    push (@fcomp,$compiler);
    push (@fcomp,@ccopt,@fopt,$pfile);
    if ((!$outfile)&&($componly)) {
	$outfile=$infile;
	$outfile =~ s/.*\///;
	$outfile =~ s/\.[Ff][Ff]*([0-9]*)$/\.o/;
    }
    if ($outfile) {
	push (@fcomp,"-o",$outfile);
    } 

    $rc=system(@fcomp);
    unlink $pfile;
    exit $rc>>8;    
    
} else {
    #
    # No input file was recognized, and/or no preprocessing is necessary.
    # Perhaps  we were just liking?
    # In any case, invoke the base compiler
    # with the arglist (minus script specific options if any)
    # 
    push (@args,"nagfor");
    push (@args, @ARGV);
    exec @args;
}
