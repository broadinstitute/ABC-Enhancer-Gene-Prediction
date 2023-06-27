Developers
==========

To start contributing to the codebase, request Anthony or Jesse to be added as 
a collaborator in the github repo.

We should incorporate bug fixes, new features, and refactoring 
continuously. We don't need to integrate experimental stuff. 

1. Start from the **dev branch**: ``git checkout dev``
2. Create a new branch: ``git checkout -b update_dev_docs``
3. Make the appropriate changes to the files you care about
	- Update the docs if applicable
	- If you made significant changes, do a full ABC run. Run CRISPRi benchmarking to ensure no regressions
4. Git add, commit, and push your changes to Github (You may have to rebase on top of the latest changes to **dev**)
5. Create a PR (Pull Request) to merge your branch into **dev**
6. Add atancoder (Anthony) as a reviewer
7. Merge PR after approval and all checks pass


Testing
-------

Series of tests that we automatically run in your PR

#. Ability to build conda environment
#. Full ABC run on chr22 matches the expected output
#. Python Linter (`black <https://pypi.org/project/black/>`_)

Testing Infrastructure
----------------------

Tests are run and managed by `CircleCI <https://app.circleci.com/pipelines/github/broadinstitute>`_. Code config can be found in ``.cicleci`` folder