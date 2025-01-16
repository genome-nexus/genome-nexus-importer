#!/bin/bash
set -e

DIR=/bitnami/mongodb

# Start MongoDB in background
# The setup.sh script (provided by the Bitnami MongoDB image) configures MongoDB, initializes the database if needed, and sets up the environment.
# It triggers import_mongo.sh if MongoDB is being initialized for the first time.
/setup.sh

# This script starts the MongoDB server using the configuration and data prepared by /setup.sh.
/run.sh