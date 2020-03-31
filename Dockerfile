FROM mongo:3.6.2

# Easier to copy everything and remove unused directories, as there's no --parent parameter in docker COPY
COPY data/ /data/
RUN rm -rf /data/*/input /data/*/tmp

RUN apt-get update && apt-get install -y wget
RUN mkdir /mongodump
RUN chmod ug+rw -R /mongodump

# Import data into mongodb
COPY scripts/import_mongo.sh /docker-entrypoint-initdb.d/
CMD ["mongod"]
