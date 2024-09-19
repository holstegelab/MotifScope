#!/bin/bash
#drop docker image motifscope if it exists
docker rmi motifscope
#build docker image motifscope
docker build -t motifscope .
