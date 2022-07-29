#!/bin/sh
## prints all sha of commits that are seen twice or more in the log
## with the same author and the same time stamp.

usage() {
    echo ""
    echo "==========================================================="
    echo "Usage:"
    echo "        $1 existing_local_git_repository"
    echo "==========================================================="
    echo "*** OR if you are already in a git repo tree *** "
    echo "-----------------------------------------------------------"
    echo "        $1 ."
    echo "-----------------------------------------------------------"
    echo ""
}

if [ $# -ne 1 ]
then
    echo "Wrong number of arguments (num_args=$#)"
    usage $0
    exit 1
else
    if [ -d $1 ]
    then
	cd $1
	chk=`git rev-parse --is-inside-work-tree 2>&1|cut -d: -f1`
	if [ $chk = "true" ]
	then
	    git log --date=format:'%s' --pretty=format:'%H:%ad%al' | awk -F: 'a[$2]++{print $0}' | cut -d: -f1
	    exit 0
	else
	    echo "ERR: $1 is not a git repository"
	    usage $0
	    exit 3
	fi
    else
	usage $0
	exit 2
    fi
fi
