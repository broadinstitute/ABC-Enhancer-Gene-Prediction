#!/bin/bash

cp -r ../../src .
docker build -t quay.io/nbarkas/abc-general-container .
rm -rf src
echo You can now push the image with "docker push quay.io/nbarkas/abc-general-container:latest"