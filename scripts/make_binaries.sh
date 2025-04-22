#!/usr/bin/env bash

set -eux

version=$(cat VERSION)

mkdir -p bin

docker build . -f docker/remeta.ubuntu20.dockerfile -t remeta_v${version}_x86_64_ubuntu20
id=$(docker create remeta_v${version}_x86_64_ubuntu20)
docker cp ${id}:/usr/local/bin/remeta ./bin/remeta_v${version}_x86_64_ubuntu20
docker rm ${id}
docker image rm remeta_v${version}_x86_64_ubuntu20

docker build . -f docker/remeta.ubuntu22.dockerfile -t remeta_v${version}_x86_64_ubuntu22
id=$(docker create remeta_v${version}_x86_64_ubuntu22)
docker cp ${id}:/usr/local/bin/remeta ./bin/remeta_v${version}_x86_64_ubuntu22
docker rm ${id}
docker image rm remeta_v${version}_x86_64_ubuntu22