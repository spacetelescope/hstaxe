.. _installing:

Installing hstaxe
=================

.. note::
   As of current testing, HSTaXe does not support Windows. Please use our
   legacy instructions for attempting to use HSTaXe on Windows

Preparing Your Local Environment
--------------------------------
We recommend using Anaconda to manage your ``hstaxe`` environment.

You may want to consider installing ``hstaxe`` in a new virtual or conda
environment to avoid version conflicts with other packages you may have
installed.

Start by creating an empty conda environment:

.. code-block:: bash

    conda create --name hstaxe-env "python>=3.8, <3.11"
    conda activate hstaxe-env

Build Prerequisites
^^^^^^^^^^^^^^^^^^^
Because the core modules of ``hstaxe`` are written in C, we require some
additional prerequisites to be installed before attempting to install ``hstaxe``:

.. code-block:: bash

    conda install gsl cfitsio make automake autoconf libtool pkg-config -y
    conda install wcstools -c https://conda.anaconda.org/conda-forge/ --override-channels -y

Some architectures (for example systems with Apple M1 processors) may not have a
wheel available for ``tables``. This can manifest as a failure to install ``tables``
during the next ``hstaxe`` installation step and can be solved by running the
following (which allows ``tables`` to build from source) before repeating the
failed install step:

.. code-block:: bash

    conda install hdf5 -c conda-forge

Some users have also reported needing to install ``openblas`` in order to get a successful
install working. If you have followed the rest of the ``hstaxe`` installation instructions
and are still having problems, it may be worth adding the following command at this step:

.. code-block:: bash

    conda install openblas -c conda-forge

Installing HSTaXe
-----------------
``hstaxe`` is distributed through PyPI. To install the latest release of ``hstaxe``:

.. code-block:: bash

    pip install hstaxe --no-cache-dir

The ``--no-cache-dir`` flag is optional, but it ensures that you are getting the most
up-to-date versions of dependencies, rather than using anything cached on your local system.

To instead install the latest development version, you can either install from our
GitHub repository directly:

.. code-block:: bash

    pip install git+https://github.com/spacetelescope/hstaxe.git

or alternatively, clone the repository locally and install the clone:

.. code-block:: bash

    git clone https://github.com/spacetelescope/hstaxe.git
    cd hstaxe
    pip install .

If you want to run the example notebooks, you will also need to install Jupyter:

.. code-block:: bash

    pip install jupyter

Legacy Astroconda Installation
------------------------------
For historical preservation, we provide the original installation instructions
for installing ``hstaxe`` via Astroconda. We preserved a premade conda
environment yaml in the repository for reproducability:

.. code-block:: bash

    conda create --name hstaxe-env --file legacy_astroconda_environment.yml
    conda activate hstaxe-env
    conda install hstaxe -c https://ssb.stsci.edu/astroconda --override-channels


Package Structure
-----------------

The ``hstaxe`` software is composed of a combination of routines written in
ANSI-C and python. Many of the python modules use the C executables to
do their work, while some perform all operations within the python
module itself. The C executables reside in the cextern directory,
while the python source routines reside in hstaxe tree.


Validating the aXe installation
-------------------------------

Jupyter notebooks for validation (and as examples of using HSTaXe) are available
at https://github.com/spacetelescope/hstaxe/tree/main/cookbooks. We recommend cloning the
repository as described above, and creating a new environment by running the following in
the directory with the ``cookbook_env.yml`` file:

.. code-block:: bash

    conda env create -f cookbook_env.yml

This will set up an environment called ``hstaxe_cookbook_acs_wfc3`` that has some packages
installed that are needed to run the notebooks but are not strictly necessary for ``hstaxe``
itself, such as ``hstcal`` and ``acstools``. Once this environment is activated, you will need
to download some data and configuration files as detailed in the ``Introduction`` section of
each notebook, after which the notebooks should run all the way through without any other
changes if ``hstaxe`` has been installed successfully.
