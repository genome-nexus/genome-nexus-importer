ARG MONGODBVERSION=7.0.28

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

# Install required tools for import script (curl for downloads, gzip for decompression)
RUN apt-get update && apt-get install -y curl gzip && rm -rf /var/lib/apt/lists/*

# Create directories for scripts and data storage
RUN mkdir -p /scripts /data /data/common_input

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

# Create log directory for our custom init
RUN mkdir -p /var/log/mongodb && chown mongodb:mongodb /var/log/mongodb

# Copy custom entrypoint wrapper to fix race condition
COPY docker-entrypoint-wrapper.sh /usr/local/bin/docker-entrypoint-wrapper.sh
RUN chmod +x /usr/local/bin/docker-entrypoint-wrapper.sh

# Use custom entrypoint that adds delay after init
ENTRYPOINT ["/usr/local/bin/docker-entrypoint-wrapper.sh"]