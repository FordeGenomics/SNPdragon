#!/bin/bash

BUILD=snpdragon
TAG=$1

docker image build --no-cache --file ./Dockerfile --tag ${BUILD}:${TAG} .
docker save "${BUILD}:${TAG}" -o "${BUILD}_${TAG}.tar"
#docker image rm "${BUILD}:${TAG}"
singularity build --force "${BUILD}_${TAG}.img" "docker-archive://./${BUILD}_${TAG}.tar"
rm "./${BUILD}_${TAG}.tar"

