# Build Instructions

1. From the `base` directory build the base image `docker build -t ersiliaio/base:latest .`
1. From the `model` directory in the `Dockerfile`:
    * Change the `FROM` line to match `ersiliaio/base:latest .`
    * Update the `ARG` line for the model you want to build
    * Build the image `docker build -t ersiliaio/eos3b5e:latest .`
1. These images should be pushed to a registry using `docker push <IMAGE>`
1. Create an app according to documentation - [https://docs.digitalocean.com/products/app-platform/how-to/create-apps/](https://docs.digitalocean.com/products/app-platform/how-to/create-apps/)
    * Be sure to set HTTP port to 80
    * Add a health check to `/healthz` endpoint
    * Change the name and plan of the app for your use case
