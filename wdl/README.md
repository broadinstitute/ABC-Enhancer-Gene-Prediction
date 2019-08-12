# Running ABC pipeline as a WDL

## Initial Setup
1. Install Cromwell and WOMTools locally
   - Download the cromwell and womtools .jar files from [here](https://github.com/broadinstitute/cromwell/releases)
   - For ease of use, put them in the same place somewhere on your machine. In this example, we'll be using `~/jars/`
2. Create your `input.json` based on the one below. If you're not running this in terra, you may want to create a script to automate the creation of this input file.

3. Create the `bedgraph.tar.gz` file used in `ABCpipeline.HiCdirTar` by tarring the `/example/input_data/HiC/bedgraph` directory
4. Download and install docker, we use [Docker Desktop](https://www.docker.com/products/docker-desktop)

## Running the Pipeline
1. Start cromwell by running `java -jar ~/jars/cromwell.jar server` from the command line
2. Run the script by running `java -jar ~/jars/cromwell.jar run /path/to/abc-pipeline.wdl --inputs /path/to/inputs.json --options /path/to/options.json`
