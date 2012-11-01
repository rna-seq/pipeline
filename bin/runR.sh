#! /bin/sh

# Author: Na Li <nali@umn.edu>
# Created: 2005/06/16 18:51:57
# Time-stamp: "Thu Dec 02 09:59:57 CST 2004 (nali@bass)"
# Purpose: run R scripts that allow command line arguments

if [ $# -lt 1 ]; then
    echo "Usage: `basename $0` Rscript.R [ arguments ]"
    exit 2
fi
RCMD="R --vanilla -q"
RSCRIPT=$1
if [ ! -f $RSCRIPT ]; then
    echo "R script '$RSCRIPT' does not exist"
    exit 3
fi
shift
echo "source (\"${RSCRIPT}\", echo = TRUE)" | ${RCMD} --args $@

