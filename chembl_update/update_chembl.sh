#!/bin/sh

VERSION=$1

if [ -z "$VERSION" ]; then
	echo No chembl db version given;
	exit;
fi;

HOST="chembl.cycqd59qnrsj.us-east-2.rds.amazonaws.com"

URL="ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_${VERSION}_postgresql.tar.gz"
FILE="chembl_${VERSION}_postgresql.tar.gz"
ARCHIVE="chembl_${VERSION}/chembl_${VERSION}_postgresql/chembl_${VERSION}_postgresql.dmp"
PASS='chembl1889'

echo updating to version $VERSION, downloading $URL

if [ ! -s $ARCHIVE ]; then
	wget $URL
	tar xzf $FILE
fi


PGPASSWORD=$PASS psql -h $HOST -U chembl -c "DROP DATABASE chembl_loading "
PGPASSWORD=$PASS psql -h $HOST -U chembl -c "CREATE DATABASE chembl_loading"
PGPASSWORD=$PASS pg_restore -h $HOST -U chembl -d chembl_loading --no-owner  --verbose $ARCHIVE
PGPASSWORD=$PASS psql -h $HOST -U chembl -c "ALTER DATABASE chembl RENAME TO chembl_old"
PGPASSWORD=$PASS psql -h $HOST -U chembl -c "ALTER DATABASE chembl_loading RENAME TO chembl"
echo Manually drop database 'chembl_old' after checking update: DROP DATABASE chembl_old
