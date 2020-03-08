#!/bin/bash

root=$1

cd /GrizliImaging/

date > ${root}.start.txt
aws s3 cp ${root}.start.txt s3://grizli-v1/IRAC/Log/MStart/ --acl public-read

aws s3 cp s3://grizli-v1/IRAC/Log/MFinished/${root}.finished.txt /tmp/ --acl public-read

irac_model_single.py ${root}

aws s3 cp /tmp/${root}.finished.txt s3://grizli-v1/IRAC/Log/MFinished/ --acl public-read
