#!/bin/bash
PMMC_BUILD_VERSION=$(awk -F"-" '{ print $2 }' <<< $(ls dist/ | head -1))
PMMC_VERSIONS_STRING=$(pip index versions pmmc | grep versions:)
PMMC_VERSIONS_STRING=${PMMC_VERSIONS_STRING#*:}
UPLOAD_TO_PYPI=1
while IFS=', ' read -ra PMMC_VERSIONS_ARRAY; do
  for VERSION in "${PMMC_VERSIONS_ARRAY[@]}"; do
    if [ "$PMMC_BUILD_VERSION" = "$VERSION" ]; then
      UPLOAD_TO_PYPI=0
    fi
  done;
done <<< "$PMMC_VERSIONS_STRING"
if [ "$UPLOAD_TO_PYPI" = 1 ]; then
  echo "Wheel version wasn't found on PyPi.";
else
  echo "Wheel was found on PyPi.";
fi
echo "perform_pypi_upload=$UPLOAD_TO_PYPI" >> $GITHUB_OUTPUT
