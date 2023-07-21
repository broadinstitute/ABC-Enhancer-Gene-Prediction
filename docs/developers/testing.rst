Testing
*******

Tests
=====

Series of checks that we continuously run on your PR and the **dev**
branch

#. Ability to build the abc conda environment
#. Code is formatted properly. We use (`Python black <https://pypi.org/project/black/>`_)
#. ABC Tests

ABC Tests
=========
These are integration tests, where we execute ABC from start to finish.
We compare the results from the test with an expected output. 

These tests 
aim to prevent any regressions that may come from updating the code or the 
conda environment.

We should be continuously adding more test to increase coverage of different 
paths of the ABC code. Here's the coverage so far

* DNase-seq + H3K27ac + K562 Hi-C run on chr22 K562 cell type

You can find all the tests and code in the `/tests` directory.

To run the ABC test suite locally

#. cd into the base repo (`ABC-Enhancer-Gene-Prediction`)
#. activate the abc-env (`conda activate abc-env`)
#. run `pytest`

Writing an ABC Test
-------------------
Simple
^^^^^^

If you want to test on certain inputs, but don't need to modify any special 
params, this option is for you.

#. Open up `tests/config/test_biosamples.tsv`
#. Move your input files into `tests/resources`
#. Add a new line indicating with the biosample name and input files references (BAM, Hi-C, ...)

	a. Try to limit the biosample to only 1 chromosome (e.g chr22). You'll need to provide an alt_TSS and alt_genes file with lines only for that chromosome.  
#. Put the expected results under `tests/expected_output/generic/{biosample_name}` 

When you run `pytest`, the test suite will automatically pick up your new biosample

Custom
^^^^^^

If you need to modify certain params to ABC scripts, like `threshold` or `pval`, you'll need to do this

#. Copy `tests/config/generic_config.yml` into your own separate yml file
#. In the new yml file, change `TEST_CONFIG_NAME` and any params
#. Create your own `biosamples.tsv` file and have your new config reference that for `biosamplesTable`
#. Put your expected results into `tests/expected_output/{new_TEST_CONFIG_NAME}/{biosample_name}`
#. Add a new test function in `tests/test_full_abc_run.py` and run the test specifying your new config file

Testing Infrastructure
======================

All the checks/tests are run and managed by `CircleCI <https://app.circleci.com/pipelines/github/broadinstitute>`_. 

CircleCI will run all the checks 

* On every pull request
* Everyday on the `dev` branch

Continuously testing ensures that the codebase is working as intended and allows 
us to catch issues quickly. 

The configuration for CircleCI can be found in the `.circleci/config.yml` file. You 
can find the CircleCI docs `here <https://circleci.com/docs/>`_
