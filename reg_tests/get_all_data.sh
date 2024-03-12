#!/bin/bash

downloaddir=("testcases unit_tests")
startdir=$(dirname "$(readlink -f $0)")
codedir=`readlink -f "${startdir}/.."` # Assumption about source code location

for dirname in ${downloaddir[@]}; do
    cd "${startdir}/$dirname"
    for name in *; do
      if [ -d $name ]; then
        echo ""
        echo "===== Getting data for $name ====="       
        sh "get_${dirname%?}_data.sh" $name
      fi
    done
done
exit 0
