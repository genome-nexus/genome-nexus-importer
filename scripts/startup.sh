#!/bin/bash

DIR=/data/db

if [ -d "$DIR/journal" ]; then
    echo "Data directory exists. Skipping."
else
    echo "Data directory not existing. Initializing seed in $DIR"
    rm -rf $DIR/*
    cp -r /seed/db/* $DIR
fi

docker-entrypoint.sh mongod
