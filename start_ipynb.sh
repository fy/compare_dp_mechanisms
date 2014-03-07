#!/bin/bash
BASEDIR=$(dirname $0)

ipython notebook --profile=nbserver --pylab inline --NotebookManager.save_script=True --NotebookManager.notebook_dir="$BASEDIR/notebooks"
