#!/bin/bash

### Mount NVME volume on AWS instances, e.g., c5d

sudo mkfs.ext4 -E nodiscard /dev/nvme1n1
sudo mount -o discard /dev/nvme1n1 /GrizliImaging/
sudo chown ec2-user /GrizliImaging
