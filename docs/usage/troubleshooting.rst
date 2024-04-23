Troubleshooting
===============

Conda Env Building Issues
-----------------------------
Building the conda environment can take a lot of resources, especially if you're not using `mamba`. Make sure the process isn't being killed 
(e.g because it's taking up too many resources in a shared environment). If it is, make sure you try building the environment with more resources. 

If you're on MacOSX, make sure to remove some of the requirements in abcenv.yml. Those requirements to delete are flagged in the file.

If there are incompatibility issues, try building off the 'release.yml' conda environment.


malloc: Heap corruption detected
--------------------------------
We've seen this happen when running on MacOSX during the prediction rule. It's an error thrown by the hicstraw library and happens the first time you use it. 
Re-running the pipeline should fix it.
