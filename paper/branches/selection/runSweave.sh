#!/bin/sh
echo "Sweave(\"selection.rnw\")" | LC_ALL="C" R --no-save --no-restore
