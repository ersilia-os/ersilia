# Build Instructions

1. From the `base` directory build the base image `docker build -t ersiliaio/base:latest .`
1. From the `model` directory in the `Dockerfile`:
    * Change the `FROM` line to match `ersiliaio/base:latest .`
    * Update the `ARG` line for the model you want to build
    * Build the image `docker build -t ersiliaio/eos3b5e:latest .`
1. Run `docker run -p 80:8080 ersiliaio/eos3b5e`
