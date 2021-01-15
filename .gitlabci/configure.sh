#!/bin/bash

## Configure a CI file based on a template file (first and only
## argument to this script).

## Known variables that will be substituted:
## @EIC_TAG@             - docker tag used for the EIC container
## @JUGGLER_TAG@         - output tag for the Juggler version

TEMPLATE_FILE=$1
OUTPUT_FILE=${TEMPLATE_FILE%.*}

echo "Configuring CI file: ${TEMPLATE_FILE}"
echo "Output will be written to: ${OUTPUT_FILE}"

TAG=$EIC_VERSION

sed "s/@EIC_TAG@/$EIC_TAG/g" $TEMPLATE_FILE | \
  sed "s/@JUGGLER_TAG@/$JUGGLER_TAG/g" > ${OUTPUT_FILE}

echo "Done"
