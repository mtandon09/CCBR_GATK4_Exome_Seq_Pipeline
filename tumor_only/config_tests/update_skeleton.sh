#!/usr/bin/bash

source_dir="/data/MyPART/mt_tools/tumor_normal_pipeline/test_6"

for file in *
do
    if [ -f "$source_dir/$file" ]
    then
        echo "Copying over $source_dir/$file"
        cp "$source_dir/$file" .
    fi
done
