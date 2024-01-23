#!/bin/bash

for name in *; do
    if [ -d $name ]; then
        echo ""
        echo "===== Getting data for $name ====="       
        sh "get_unit_test_data.sh" $name
    fi
done
exit 0
