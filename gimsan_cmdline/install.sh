#!/bin/bash

#Install GIMSAN
#See README for instructions

tar -zxvf weblogo.2.8.2.tar.gz

cd gibbsmarkov/code
make clean
make
cd ../../

cd column_dependency_app/code
make clean
make
cd ../../
