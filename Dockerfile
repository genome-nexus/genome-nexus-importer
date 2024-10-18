# This base image starts up mongo
# This version needs to correspond with the helm chart version
ARG MONGODBVERSION=4.0.12
FROM bitnami/mongodb:${MONGODBVERSION} as build

# Use .dockerignore file to ignore unwanted files
# These files are used by import_mongo.sh to initialize mongo
# Creating directories as root
# Set user back to the one in base image
USER root
RUN mkdir -p /data
COPY data/ /data/

ARG ARG_REF_ENSEMBL_VERSION
ENV REF_ENSEMBL_VERSION=${ARG_REF_ENSEMBL_VERSION}
ARG SPECIES=homo_sapiens
ARG MUTATIONASSESSOR=false

# Import data into mongodb
COPY scripts/import_mongo.sh /docker-entrypoint-initdb.d/
COPY scripts/setup.sh /setup.sh
RUN /setup.sh

FROM bitnami/mongodb:${MONGODBVERSION}
COPY --from=build /bitnami/mongodb /bitnami/seed
COPY /scripts/startup.sh /startup.sh

USER root
RUN chown -R 1001 /bitnami/seed
USER 1001

CMD [ "/startup.sh" ]
