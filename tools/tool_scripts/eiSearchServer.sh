#!/bin/bash

DIR=$(dirname $0)

$DIR/eiSearchServer.R 2>&1|logger --id --tag eiSearchServer &

PID=$(jobs -p)

#echo pid $PID

echo $PID > /var/run/eiSearchServer/eiSearchServer.pid
