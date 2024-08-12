# This base image starts up mongo
# This version needs to correspond with the helm chart version
ARG MONGODBVERSION=4.0.12
FROM docker.io/mongo:${MONGODBVERSION} as build

# Use .dockerignore file to ignore unwanted files
# These files are used by import_mongo.sh to initialize mongo
# Creating directories as root
# Set user back to the one in base image
USER root
RUN mkdir -p /data
RUN mkdir -p /seed && chmod 777 /seed
COPY data/ /data/

ARG ARG_REF_ENSEMBL_VERSION
ENV REF_ENSEMBL_VERSION=${ARG_REF_ENSEMBL_VERSION}
ARG SPECIES=homo_sapiens
ARG MUTATIONASSESSOR=false

# Import data into mongodb
COPY scripts/import_mongo.sh /docker-entrypoint-initdb.d/import_mongo.sh
COPY scripts/init.sh /init.sh
RUN /init.sh mongod

RUN ls -la /data/db

FROM docker.io/mongo:${MONGODBVERSION}
COPY --from=build /seed /seed
COPY /scripts/startup.sh /startup.sh

CMD [ "/startup.sh" ]
