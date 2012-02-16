#!/bin/sh

while true
do
    python test_leaks.py
    if [ $? != 0 ]; then
	break
    fi
done
