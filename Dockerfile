ARG MONGODBVERSION=8.2.3

FROM mongo:${MONGODBVERSION}

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

# Create directories for scripts and data storage
RUN mkdir -p /scripts /data

# Copy data to directory inside the container
COPY data/ /data/

# Store environment variables in a file for persistence
RUN echo "export REF_ENSEMBL_VERSION=${REF_ENSEMBL_VERSION}" >> /scripts/persisted_env.sh && \
    echo "export SPECIES=${SPECIES}" >> /scripts/persisted_env.sh && \
    echo "export MUTATIONASSESSOR=${MUTATIONASSESSOR}" >> /scripts/persisted_env.sh && \
    echo "export SLIM_MODE=${SLIM_MODE}" >> /scripts/persisted_env.sh && \
    chmod +x /scripts/persisted_env.sh

# Change ownership to mongodb user (official image uses mongodb:mongodb)
RUN chown -R mongodb:mongodb /data /scripts

# Copy the MongoDB initialization script
# This directory is automatically scanned and executed by MongoDB during first database initialization
COPY scripts/import_mongo.sh /docker-entrypoint-initdb.d/
RUN chmod +x /docker-entrypoint-initdb.d/import_mongo.sh

# Use default MongoDB entrypoint - no need to override
