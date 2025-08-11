#!/usr/bin/env bash

SCRIPT_DIR="$(dirname "$(realpath "$0")")"
revisionlog="${SCRIPT_DIR}/../../../revisions.log"
tail -n1 $revisionlog | sed 's/Branch \(.*\) deployed.*/\1/'
