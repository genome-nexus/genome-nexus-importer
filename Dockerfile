FROM mongo:3.6.2
COPY export/* /import/export/
COPY scripts/* /import/scripts/
RUN mongod --fork --logpath /var/initdb.log && .//import/scripts/import_mongo.sh mongodb://127.0.0.1:27017/annotator && rm -rf /import
CMD ["mongod"]
