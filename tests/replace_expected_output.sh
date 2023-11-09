#!/bin/bash
# Used when we want to change the expected output to match the latest test output run
set -e

rm -rf expected_output/*
cp -R test_output/* expected_output/