#!/bin/bash

## Run the tests on the cloud. This requires a local cromwell that has access to google api
## See instructions here to configure this: https://cromwell.readthedocs.io/en/develop/tutorials/PipelinesApi101/
## and add the google.conf file to the current directory

CROMWELL_JAR_PATH=/Users/nbarkas/jars/cromwell-44.jar

java -Dconfig.file=google.conf -jar $CROMWELL_JAR_PATH run test_ABC.wdl -i test_inputs_cloud.json
