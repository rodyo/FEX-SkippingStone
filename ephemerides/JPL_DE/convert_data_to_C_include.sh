#!/bin/bash

# This command serializes the static binary data file JPLEPH20002060Bin.405
# into a C-style include file. This allows one to directly include the data
# into the binary, without the need to read any files, check filenames, etc.

xxd -i JPLEPH20002060Bin.405 JPLEPH20002060Bin.c
xxd -i JPLEPH19402100Bin.405 JPLEPH19402100Bin.c