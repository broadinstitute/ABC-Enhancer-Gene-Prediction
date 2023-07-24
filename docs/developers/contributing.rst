Contributing
============

To start contributing to the codebase, request Anthony or Jesse to be added as 
a collaborator in the github repo.

What we should incorporate

* Bug fixes
* New features
* Code improvements

What we shouldn't incorporate

* Experimental code
* User specific configs (e.g my version of `biosamples.tsv``)
* Reference files not used in the generic version of the code

Process
------------------

1. Start from the **dev branch**: ``git checkout dev``
2. Create a new branch: ``git checkout -b update_dev_docs``
3. Make the appropriate changes to the files you care about
	- Update the docs if applicable
	- Run/Add/Update tests. See :doc:`testing`
	- If you made significant changes, do a full ABC run. Run `CRISPRi benchmarking <https://github.com/EngreitzLab/CRISPR_comparison>`_ to ensure no regressions
4. Git add, commit, and push your changes to Github (You may have to rebase on top of the latest changes to **dev**)
5. Create a PR (Pull Request) to merge your branch into **dev**
6. Add atancoder (Anthony) as a reviewer
7. Merge PR after approval and all checks pass