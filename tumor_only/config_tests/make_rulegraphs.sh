#!/usr/bin/bash

configsdir="example_configs"  ## where config jsons are stored
realconfigfile="references_hg38.json"  ## The filename that the run command will expect

pngdir="example_configs/rulegraphs"  ## Where to output rule graph PNGs
if [ ! -d $pngdir ]; then mkdir -p $pngdir; fi

## Loop through jsons, symlink as "real" config, and generate a rule graph for each
for json in $configsdir/*.json; do
    if [ -f $realconfigfile ]; then rm $realconfigfile; fi
    ln -s $json $realconfigfile
    
    jsonlabel=$(echo $(basename $json) | sed -e 's/.json$//')
    pngfile="$pngdir/rules.$jsonlabel.png"
    ./run.sh --rulegraph $pngfile
done
