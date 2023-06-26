Developers
==========

To start contributing to the codebase, request Anthony or Jesse to be added as 
a collaborator in the github repo.

We should incorporate bug fixes, new features, and refactoring 
continuously. We don't need to integrate experimental stuff. 

1. Start from the **dev branch**: ``git checkout dev``
2. Create a new branch: ``git checkout -b update_dev_docs``
3. Make the appropriate changes to the files you care about
4. Git add, commit, and push your changes to Github (You may have to rebase on top of the latest changes to **dev**)
5. Create a PR (Pull Request) to merge your branch into **dev**
6. Add atancoder (Anthony) as a reviewer
7. Merge PR after approval and tests pass


Testing
-------

Series of tests that we automatically run in your PR

1. Full ABC run on chr22 matches the expected output
2. Python Linter (`black <https://pypi.org/project/black/>`_)

Testing Infrastructure
----------------------

Tests are run and managed by CircleCI. Code config can be found here (insert link)

Everytime a PR is created, we run the tests. All tests must pass in order to get merged into dev
