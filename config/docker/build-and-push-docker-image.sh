#! /bin/bash

DOCKER_DESTINATION='dockerreg.chtc.wisc.edu/dhoconno/mhc-alleles-from-gdna-amplicons:27751'

# build Docker image

docker build -t $DOCKER_DESTINATION -f Dockerfile .

# push Docker image

docker push $DOCKER_DESTINATION
