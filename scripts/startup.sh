#!/bin/bash
set -e

# Load persisted environment variables
if [ -f /scripts/persisted_env.sh ]; then
    source /scripts/persisted_env.sh
fi

# The new Bitnami MongoDB 8.x uses /opt/bitnami/scripts/mongodb/entrypoint.sh
# Just call it without arguments - it handles everything internally
if [ -f /opt/bitnami/scripts/mongodb/entrypoint.sh ]; then
    exec /opt/bitnami/scripts/mongodb/entrypoint.sh mongod
elif [ -f /setup.sh ] && [ -f /run.sh ]; then
    # Legacy Bitnami structure (4.x)
    /setup.sh
    exec /run.sh
else
    # Fallback
    exec mongod
fi