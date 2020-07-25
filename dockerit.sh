#!/bin/sh
docker build -t chapter13 .
docker run -it --entrypoint /bin/bash chapter13
