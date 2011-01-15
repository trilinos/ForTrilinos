#!/bin/sh

# 
# Wrapper for compilers having troubles with preprocessor.
#
# Works for only 1 input file at a time.
# To make it work for more than 1 file, I should
# probably rewrite it in perl. 
# 
# Assumes:
#  1. /bin/sh is a Bourne-like shell
#  2. sed is available
#  3. gcc is available
#  4. Current directory is writable
#

#
# Written for a version of NAG that could not preprocess
# lines longer than 132 chars.
# Only these two variables depend on the compiler:
#   the compiler executable name
#   the option (if any) to force preprocessing.
compiler="nagfor"
forcepre="-fpp"

#
# The script itself
# 
origarg="$*"
preproc="no"
infile=
outfile=
copt=
fopt=

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
#      -Idir -I dir -Dmacro
#      I can never remember if -I dir is allowed or not..
#      
#   4. Output clause: -o outfile
#   
#   5. Option to force preprocessing even for files in 2.
#   
#   6. Any other option, irrelevant to the preprocessor
#      and intended for the base fortran compiler
#
while [ $# -gt 0 ]
do
    case "$1" in 
	*.F[0-9][0-9]|*.ff[0-9][0-9]) infile="$1"
	    preproc="yes"
	    ;;
	*.F|*.ff) infile="$1"
	    preproc="yes"
	    ;;
	*.f[0-9][0-9]|*.f) infile="$1"
	    ;;
	-I) copt="$copt $1 $2"
	    shift
	    ;;
	-I*|-D*) copt="$copt $1"
	    ;;
	-o) shift
	    outfile=$1
	    ;; 
	$forcepre) preproc="yes"
	    ;;
	*) fopt="$fopt $1"
	    ;;
    esac
    shift
done


if test "x$preproc" = "xyes"; then 
    #
    # Build a name for the temporary file
    # holding preprocessed source. Should be the same basename
    # with single lowercase "f" plus digits, if any.
    # Take out prefix to write in current directory, just
    # in case source dir is not writable.
    # 
    pfile=`echo $infile | sed 's/\.[Ff][Ff]*\([^\.]*\)$/.f\1/
s|.*/||'`
    #
    # make sure intermediate file name is not used already
    # 
    while [ -f $pfile ]
    do
	pfile="XtmpX_$pfile"
    done
    #
    # Preprocess:
    #   -E    only runs the preprocessor
    #   -undef does not predefine GNU specific macros
    #          (but keeps system macros)
    #   -P do not add line number directives
    #      (just in case the Fortran compiler has troubles)
    #   Make sure to report the return code if gcc -E fails
    # 
    gcc -undef -E -P $copt -o $pfile $infile || exit $?
    if test "x$outfile" = "x"; then
	#
	# Make sure outfile has the right name if not already forced
	# Write in current directory, the compiler would do it anyway
	# but potentially with the wrong name. 
	# 
	outfile=`echo $infile | sed 's/\.[Ff][Ff]*[^\.]*$/.o/
s|.*/||'`
    fi
    $compiler $copt $fopt $pfile -o $outfile 
    #
    # Make sure to exit with the compiler's return code
    # 
    rc=$?
    /bin/rm -f $pfile
else
    $compiler $origarg
    rc=$?
fi

exit $rc