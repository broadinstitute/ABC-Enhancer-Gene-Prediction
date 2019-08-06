#!/bin/bash

## Copy the source files here
cp -r ../../src .
## Build the docker image
docker build -t quay.io/nbarkas/abc-general-container .
## Remove the copy of the source files
rm -rf src
## Let user know they can now push the image
echo You can now push the image with "docker push quay.io/nbarkas/abc-general-container:latest"