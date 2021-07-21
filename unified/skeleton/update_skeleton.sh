#!/usr/bin/bash

source_dir="/scratch/tandonm/exome_pipeline_dev_mt/unified/skeleton"

for file in *
do
    if [ -f "$source_dir/$file" ]
    then
        echo "Copying over $source_dir/$file"
        cp "$source_dir/$file" .
    fi
done


for file in $source_dir/workflow/*
do
    if [ -f "$file" ]
    then
        echo "Copying over $file"
        cp -f "$file" ./workflow
    fi
done

for file in $source_dir/workflow/rules/*
do
    if [ -f "$file" ]
    then
        echo "Copying over $file"
        cp -f "$file" ./workflow/rules
    fi
done

for file in $source_dir/config/*
do
    if [ -f "$file" ]
    then
        echo "Copying over $file"
        cp -f "$file" ./config
    fi
done
