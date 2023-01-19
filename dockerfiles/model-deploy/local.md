# Tips

1. From the `base` directory build the base image `docker build -t ersiliaos/base:latest .`
1. From the `model` directory in the `Dockerfile`:
    * Change the `FROM` line to match `ersiliaos/base:latest .`
    * Update the `ARG` line for the model you want to build
    * Build the image `docker build -t ersiliaos/eos4e40:latest .`

# Troubleshooting
1. Enter the base image in interactive mode: `docker run -i -t --entrypoint bash ersiliaos/base:latest`
1. Inside docker, activate ersilia: `source activate; conda activate ersilia`
1. Fetch the model again and inspect as usual: `ersilia -v fetch eos4e40`