#!/bin/bash

IFS='
'
for i in *.fa
do
mv $i "${i// /_}"
done

