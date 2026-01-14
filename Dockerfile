ARG MONGODBVERSION=8.2.3

FROM public.ecr.aws/bitnami/mongodb:${MONGODBVERSION}

# Define build arguments
ARG ARG_REF_ENSEMBL_VERSION
ENV REF_ENSEMBL_VERSION=${ARG_REF_ENSEMBL_VERSION}
ARG SPECIES=homo_sapiens
ENV SPECIES=${SPECIES}

# Define additional annotation resources arguments
ARG MUTATIONASSESSOR=false
ENV MUTATIONASSESSOR=${MUTATIONASSESSOR}
ARG SLIM_MODE=false
ENV SLIM_MODE=${SLIM_MODE}

USER root

# Create directories for scripts and data storage and copy data to directory inside the container
RUN mkdir -p /scripts /data
COPY data/ /data/
COPY scripts/startup.sh /scripts/

# Make all scripts in the /scripts directory executable
RUN chmod +x /scripts/*.sh

# Store environment variables in a file
# When deploying using bitnami chart, env might be overridden when starting the container, use a file to persist custom env
RUN echo "export REF_ENSEMBL_VERSION=${REF_ENSEMBL_VERSION}" >> /scripts/persisted_env.sh && \
    echo "export SPECIES=${SPECIES}" >> /scripts/persisted_env.sh && \
    echo "export MUTATIONASSESSOR=${MUTATIONASSESSOR}" >> /scripts/persisted_env.sh && \
    echo "export SLIM_MODE=${SLIM_MODE}" >> /scripts/persisted_env.sh && \
    chmod +x /scripts/persisted_env.sh

# Change ownership of the /data directory and its contents to non-root user 1001
RUN chown -R 1001 /data

# Switch to the non-root user
USER 1001

# Copy the MongoDB initialization script into the /docker-entrypoint-initdb.d/ directory
# This directory is automatically scanned and executed by MongoDB during the first database initialization
COPY scripts/import_mongo.sh /docker-entrypoint-initdb.d/

# Set the default command to execute the custom startup script when the container runs
# The startup script arranges the setup and start of MongoDB
CMD [ "/scripts/startup.sh" ]
