#!/bin/bash
from lib import get_enz

for file in ./seq_a_test/* ; do
	echo "$file"
	./main_multi.py -seq "$file"
	echo "ok"

./get_enz.py
done
