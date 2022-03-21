#!/bin/sh

#
# initialize the current LCIO version by sourcing this in the top level directory
#

export LCIO=`pwd`
export PATH=$LCIO/bin:$LCIO/tools:$PATH
export LD_LIBRARY_PATH=$LCIO/lib:$LCIO/lib64:$LD_LIBRARY_PATH
export PYTHONPATH=$LCIO/src/python:$PYTHONPATH
alias pylcio='python $LCIO/src/python/pylcio.py'
