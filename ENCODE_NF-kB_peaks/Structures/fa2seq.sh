#!/bin/bash

for i in *.fa
do
  sed -i '1s/^/; comment_line\n/' $i #add a single comment line
  sed -i 's/>//' $i #remove the > symbol making a title
  sed -i '${s/$/1/}' $i #append the 1 to the last line
done

