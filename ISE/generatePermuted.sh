#!/bin/sh
head -1 cleaned > permuted
head -2 cleaned | tail -1 | ./permuteRow.py >> permuted
tail -n +3 cleaned >> permuted
