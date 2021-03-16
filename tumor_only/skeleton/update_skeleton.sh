#!/usr/bin/bash

source_dir="/scratch/tandonm/exome_pipeline_dev_mt/tumor_only/skeleton"

for file in *
do
    if [ -f "$source_dir/$file" ]
    then
        echo "Copying over $source_dir/$file"
        cp "$source_dir/$file" .
    fi
done
