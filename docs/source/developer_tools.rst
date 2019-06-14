Developer Tools
===============

Run ForTrilinos development environment in a Docker container
-------------------------------------------------------------

To start a container from the ForTrilinos pre-built Docker image that is used in the
automated build on Jenkins, run:

.. code:: bash

    [host]$ cd docker
    [host]$ echo COMPOSE_PROJECT_NAME=$USER > .env # [optional] specify a project name
    [host]$ docker-compose pull # pull the most up-to-date version of the ForTrilinos base image
    [host]$ docker-compose up -d # start the container

This will mount the local ForTrilinos source directory into the container at
``$TRILINOS_DIR/packages/ForTrilinos``. The environment variable ``TRILINOS_DIR``
is already defined and contains the path to a version of Trilinos that
has been downloaded into the ForTrilinos base image. The exact version of the
Trilinos version is given by a hash in ``docker/trilinos_version``. We recommend
you use a ``.env`` file to specify an alternate project name (the default being
the directory name, i.e. ``docker``).  This will let you run multiple isolated
environments on a single host.  Here the service name will be prefixed by your
username which will prevent interferences with other developers on the same
system.

Then to launch an interactive Bash session inside that container, do:

.. code:: bash

    [host]$ docker-compose exec fortrilinos_dev bash

Patch, configure, build, and test as you would usually do:

.. code:: bash

    [container]$ cd $TRILINOS_DIR
    [container]$ cd $TRILINOS_DIR/packages/ForTrilinos
    [container]$ mkdir build && cd build
    [container]$ ../scripts/docker_cmake
    [container]$ make -j<N>
    [container]$ ctest -j<N>

Do not forget to cleanup after yourself:

.. code:: bash

    [container]$ exit
    [host]$ docker-compose down # stop and remove the container

Auto-generating source files using SWIG
---------------------------------------

ForTrilinos carries generated ``.F90`` and ``_wrap.cxx`` files with its source. Thus, introducing new ``.i`` files
requires a developer to generate the corresponding wrapper and proxy files. To do this automatically, the CMake
configuration of the build must include

.. code-block:: bash

    -D ForTrilinos_ENABLE_DeveloperMode=ON

It also requires a SWIG installation in the ``$PATH`` with Fortran generator enabled. It is available at
`swig-fortran/swig <https://github.com/swig-fortran/swig>`_.

If one does simultaneous development in both SWIG and ForTrilinos, it is convenient to skip the ``make install`` step.
This can be done by adding the following configuration options to the script (assuming you build swig in-source):

.. code-block:: bash

    -D SWIG_EXECUTABLE="$SWIG_DIR/swig"
    -D SWIG_DIR="$SWIG_DIR/Lib"

.. warning::

    ForTrilinos does not automatically pick up the changes in the files in the SWIG library, such as
    ``fortypemaps.swg``. Therefore, to regenerate all wrapper files after changes in SWIG one must touch all ``.i``
    files in ForTrilinos.


.. .. warning::
.. 
    .. When using the developer mode, all of the patches ``scripts/patches`` directory of the ForTrilinos source tree must be applied.  If you previously used the ``scripts/patches/apply-patches`` script, this was done for you. Otherwise, be sure that all of the patches are applied.
