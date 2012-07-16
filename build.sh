#!/bin/sh

cmake ./
make

tar -c -z -f ./results.tar.gz ./bin 