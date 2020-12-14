#!/bin/bash

DIR=/bitnami/mongodb

if [ -d "$DIR/data/db" ]; then
    echo "Data directory exists. Skipping."
else
    echo "Data directory not existing. Initializing seed in $DIR"
    rm -rf $DIR/db
    cp -r /bitnami/seed/* $DIR
fi

/setup.sh
/run.sh