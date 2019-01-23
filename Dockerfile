FROM mongo:3.6.2
COPY data/ /data/
COPY scripts/import_mongo.sh /docker-entrypoint-initdb.d/
CMD ["mongod"]
