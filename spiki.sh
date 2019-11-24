#!/usr/bin/env bash

SCRIPT_LOCATION_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )" 


export PATH=$SCRIPT_LOCATION_DIR:$PATH

python $SCRIPT_LOCATION_DIR/spiki.py
