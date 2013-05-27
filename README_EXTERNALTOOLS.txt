
The symlink externaltools links to ../externaltools/externaltools

Checkout externaltools in the same directory containing this higgstautau
repository:

svn co svn+ssh://${USER}@svn.cern.ch/reps/atlasphys/Physics/Higgs/HSG4/software/common/externaltools/trunk externaltools

See the README in externaltools for further instructions on how to checkout
and build all of the tools.

externaltools is a RootCore alternative using the waf build system.
externaltools creates a Python package under which all tools are subpackages.
