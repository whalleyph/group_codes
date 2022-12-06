#!/bin/bash
for py in *.py; do
    echo Processing $py
    STR=$(head -n 1 $py)
    SUB='#!'
    if [[ "$STR" == *"$SUB"* ]]; then
        echo "It's there."
    fi
    sed -i 1d $py
done