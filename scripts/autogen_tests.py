#!/usr/bin/env python
#Copyright 2017-2018, UT-Battelle, LLC
#
#SPDX-License-Identifier: BSD-3-Clause
#License-Filename: LICENSE
import os
import re
import sys
from argparse import ArgumentParser

"""
This script is a total hack that starting by reading fortpetra.f90 to create
skeletons for unit tests for Tpetra objects.  I've since modified it a bit to
work for other packages, but I am sure that it will have to be modified when
working with a new package to remove any implicit assumptions from its origins
with Tpetra
"""

class NamespaceTestGenerator:
    def __init__(self, namespace, procs, filename, f90_module):
        self.namespace = namespace
        self.procs = procs
        self.filename = filename
        self.f90_module = f90_module

    def write(self, stream):
        if stream is None:
            stream = sys.stdout

        # Write the preamble to the test suite
        stream.write("""\
!Copyright 2017-2018, UT-Battelle, LLC
!
!SPDX-License-Identifier: BSD-3-Clause
!License-Filename: LICENSE
program test_{0}
#include "<CONFIG_FILE>"
#include "FortranTestUtilities.h"
  use iso_fortran_env
  use, intrinsic :: iso_c_binding
  use forteuchos
  use {3}

  implicit none
  type(TeuchosComm) :: comm
  character(len={1}), parameter :: FILENAME="{2}"

  SETUP_TEST()

#if FORTRILINOS_USE_MPI
  call comm%create(MPI_COMM_WORLD); FORTRILINOS_CHECK_IERR()
#else
  call comm%create()
#endif

""".format(self.namespace, len(self.filename)+2, self.filename, self.f90_module))

        # Write call to each individual test
        for proc in self.procs:
            stream.write('  ADD_SUBTEST_AND_RUN({0}_{1})\n'.format(self.namespace, proc.name))
            continue

        # Error reporting, clean up, MPI shutdown
        stream.write('\n  call comm%release()\n\n')
        stream.write('  TEARDOWN_TEST()\n\ncontains\n\n')

        # Now write each individual test
        for proc in self.procs:
            stream.write('  ! {0} !\n'.format(proc.name.center(74, '-')))
            stream.write(proc.generic_test())
            stream.write('\n')

        # Finish it all off
        stream.write('end program test_{0}'.format(self.namespace))


class Procedure:
    def __init__(self, namespace, name, args, proctype, instructions):
        self.namespace = namespace
        self.name = name
        self.args = args
        self.type = proctype
        self.instructions = instructions
        self.argtypes = self.parse_instructions_for_argtypes()

    def parse_instructions_for_argtypes(self):
        def keep(x):
            if not x.split():
                return False
            if 'intent(' in x:
                return False
            if x.strip() == 'target':
                return False
            return True
        regex = '(?P<decl>.*) :: (?P<names>.*)'
        argtypes = [None] * len(self.args)
        for line in self.instructions.split('\n'):
            if not line.split():
                continue
            match = re.search(regex, line)
            if not match:
                continue
            decl = [x.strip() for x in match.group('decl').split(',') if keep(x)]
            names = match.group('names')
            assert '(' not in names
            for (i, arg) in enumerate(self.args):
                if re.search(r'\b{0}\b'.format(arg), names):
                    assert argtypes[i] is None
                    argtypes[i] = decl
        if any([x is None for x in argtypes]):
            assert False, 'Unfound argtype'
        return argtypes

    def generic_test(self):
        testname = '{0}_{1}'.format(self.namespace, self.name)
        s = '''\
  FORTRILINOS_UNIT_TEST({0})
    {1}
    OUT0("Starting {0}!")

    success = .false.

    {2}
    {3}{4}

    write(*,*) '{0}: Test not yet implemented'

    OUT0("Finished {0}!")

  END_FORTRILINOS_UNIT_TEST({0})
'''
        # Form declarations and initializations
        decl = []
        init = []
        fini = []
        obj_init = None
        obj_fini = None
        the_args = []
        for (i, arg) in enumerate(self.args):
            argtype = ', '.join(self.argtypes[i])
            name = 'Obj' if arg == 'self' else arg
            if name == 'comm':
                continue
            if argtype.lower().startswith(('class', 'type')):
                argtype = argtype.replace('class', 'type')
                i = ['!call {0}%create(); TEST_IERR()'.format(name),]
                f = '!call {0}%release(); TEST_IERR()'.format(name)
                if name == 'Obj':
                    obj_init = i
                    obj_fini = f
                else:
                    init.extend(i)
                    fini.append(f)
            elif 'dimension(:)' in argtype:
                assert name != 'Obj'
                argtype = argtype.replace('dimension(:)', 'allocatable')
                name = '{0}(:)'.format(name)
                init.append('!allocate({0}(0))'.format(name))
                fini.append('!deallocate({0})'.format(name))
            elif name != 'fresult':
                assert name != 'Obj'
                init.append('{0} = 0'.format(name))
            if name != 'Obj':
                the_args.append(name)
            decl.append('{0} :: {1}'.format(argtype, name))
        assert obj_init is not None
        assert obj_fini is not None
        init.extend(obj_init)
        fini.append(obj_fini)
        if self.type == 'subroutine':
            call = '!call Obj%{0}({1}); TEST_IERR()\n'.format(
                self.name, ', '.join(the_args))
        elif self.type == 'function':
            call = '!fresult = Obj%{0}({1}); TEST_IERR()\n'.format(
                self.name, ', '.join(the_args))
        else:
            assert False, 'Wrong type'
        decl = '\n    '.join(decl)
        init = '\n    '.join(init)
        fini = '' if not fini else '\n    ' + '\n    '.join(fini)
        return s.format(testname, decl, init, call, fini)



def main():
    p = ArgumentParser()
    p.add_argument('filename')
    p.add_argument('--package',
                   help='Package name.  Inferred from filename if not give')
    args = p.parse_args()

    filename = args.filename
    assert os.path.isfile(filename)
    package = args.package
    if package is None:
        # Infer package from filename:
        # for<name>.f90 => package = Name
        root, ext = os.path.splitext(os.path.basename(filename))
        package = re.sub(r'^for', '', root).title()
    contents = open(filename).read()

    f90_module = re.search(r'\Wmodule\s+(?P<name>\w+)', contents)
    if f90_module is None:
        raise Exception('Unable to determine module name')
    f90_module = f90_module.group('name').strip()

    # Find public procedures of the fortpetra module
    namespaces = {}
    re1 = re.compile('procedure\s+::\s+(?P<proc>\w+)\s+=\>\s+swigf_(?P<proc2>\w+)\W')
    re2 = '(?ims){2} swigf_{0}_{1}\(' \
          '(?P<args>.*?)\)(?P<instructions>.*?)' \
          'end {2}'
    for (proc, proc2) in re1.findall(contents):

        if proc in ('create', 'release'):
            continue

        #if proc.startswith(('set_', 'get_')):
        #    # Not sure about RowInfo right now
        #    assert 'RowInfo' in proc2
        #    continue

        namespace, proc2 = proc2.split('_', 1)

        # I'm not sure this assertion is right, but it should warn me of a situation
        # I have not anticipated
        assert proc == proc2, 'Mismatching procedure names: {0}, {1}'.format(proc, proc2)

        for proctype in ('subroutine', 'function'):
            regex = re2.format(namespace, proc, proctype)
            match = re.search(regex, contents)
            if match:
                break
        else:
            raise Exception('Unable to find procedure {0} for {1}'.format(proc2, namespace))


        args = [x.strip() for x in match.group('args').split(',') if x.split()]
        theproc = Procedure(namespace, proc2, args, proctype,
                            match.group('instructions'))
        namespaces.setdefault(namespace, []).append(theproc)

    for (namespace, procs) in namespaces.items():
        filename = 'Test{0}_autogen.f90'.format(namespace)
        s = namespace.lower().replace(package, package+'_')
        write_filename = 'test_{0}.f90'.format(s)
        doc = NamespaceTestGenerator(namespace, procs, write_filename, f90_module)
        with open(filename, 'w') as fh:
            doc.write(fh)

if __name__ == '__main__':
    main()
