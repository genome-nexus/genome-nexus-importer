# This base image starts up mongo
# This version needs to correspond with the helm chart version
FROM bitnami/mongodb:4.0.12

# Use .dockerignore file to ignore unwanted files
# These files are used by import_mongo.sh to initialize mongo
# Creating directories as root
# Set user back to the one in base image
USER root
RUN mkdir -p /data
COPY data/ /data/
USER 1001

# Import data into mongodb
COPY scripts/import_mongo.sh /docker-entrypoint-initdb.d/

