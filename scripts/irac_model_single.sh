#!/bin/bash

root=$1

cd /GrizliImaging/

date > ${root}.start.txt
aws s3 cp ${root}.start.txt s3://grizli-v1/IRAC/Log/MStart/ --acl public-read

irac_model_single.py ${root}

rm -rf /GrizliImging/${root}/

aws s3 cp /tmp/${root}.finished.txt s3://grizli-v1/IRAC/Log/MFinished/ --acl public-read
