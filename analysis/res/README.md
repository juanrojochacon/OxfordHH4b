Merging (and Plotting) README
==============

Key Changes
-----------

### Merging
- `yodamerge` rewritten in C++ to avoid needing to load Python and `yoda` module
- New `yoda-merge` also writes out a flat file, removing need to run `yoda2flat`
- Merging scripts ported to ZSH. Compiling with `zcompile` speeds these up
- In total, merging is sped up by a factor of 4 to 6

### Plotting
- New `plotting` Python modules has been added providing two functions:
  - `plot1D` for plotting (multiple) 1D histograms, normalized by area.
  - `plot2D` for plotting a single 2D histogram

Usage
------------

### Environment Setup
The code requires some environment setup. Much of the Python related setup has
been moved into a virtualenv using Python 2.7.13 built with GCC 4.9.3, and with
Numpy, Matplotlib, and Seaborn (needed for the plotting module) installed.

All environment setup can be done on the PPLXINT machines by sourcing
`/home/stanislaus/PUBLIC/setupEnv.sh` (this file is world readable).

### Merging
Compile the appropriate merging script (`mergeSubSamples.zsh` or `mergeBackground.zsh`)
with `zcompile`, then run it from ZSH. If running on the batch system, ensure the
submitted script uses ZSH: the shebang line should be
`#!/usr/bin/env zsh`. An example script is below
~~~zsh
#!/usr/bin/env zsh
#PBS -l cput=01:59:00
#PBS -l walltime=01:59:00
#PBS -o /home/stanislaus/oxfordhh4b/analysis/batchlog/
#PBS -e /home/stanislaus/oxfordhh4b/analysis/batchlog/
cd $PBS_O_WORKDIR
. /home/stanislaus/PUBLIC/setupEnv.sh
#tar -cJf baseline_noPU_atlas_qcd.tar.xz baseline_noPU_atlas_qcd
cd baseline_noPU_atlas_qcd
../mergeSubSamples.zsh diHiggs
../mergeSubSamples.zsh SHERPA
~~~

Unfortunately, `yoda-merge.cpp` doesn't currently build on PPLXINT. Either use the binary in
the git repo, or build it with `gcc -o yoda-merge yoda-merge.cpp -O2 -lboost_program_options`
on an x86_64 (i.e. any recent) Linux system with a recent GCC and Boost. I've only tested with
GCC 6.3 and Boost 1.63, but it *should* work with slightly older versions (likely GCC 5 and above,
though adding `-std=gnu++14` may be advisable on versions older than 6.1).

### Plotting
The main factor causing plotting to be slow was the time taken to load Python, and the Numpy and
Matplotlib modules. This can be mitigated by doing all plotting in one Python script. The `plotting`
module helps with this. Import it with `from plotting import *`, then use the `plot1D` and `plot2D`
functions to do all plotting in one script.

Assuming the environment has been setup as described above (i.e. you are using the virtualenv), the
shebang line for any Python script should be `#!/usr/bin/env python2`. Do *not* use `#!/usr/bin/python`
as this will use the system Python interpreter. `# !/usr/bin/env python` (no 2) should suffice, but may
cause problems if / when the naming scheme for Python interpreters switches to `python` calling Python
3 by default.
