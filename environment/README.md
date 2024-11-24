# Docker Environment Setup for CodeOcean

This directory contains Docker configuration for running the pipeline in CodeOcean environment.

## Structure

- `Dockerfile`: Defines the container image with all required dependencies
- `postInstall`: Contains additional setup scripts run after container creation (optional)

## Requirements

- Docker version 20.10 or higher
- At least 8GB RAM allocated to Docker
- 20GB free disk space

## Usage

Read the [CodeOcean documentation](https://codeocean.com/docs/codeocean-environment) for information on how to set up and run the environment.
   ```

## Note

This Docker configuration is specifically optimized for CodeOcean execution. For local development or other cloud environments, please refer to the main installation guide in the project root.
