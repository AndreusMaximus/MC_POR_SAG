#!/bin/bash

input_file="trace_viewer.txt"
output_file="trace_out.txt"
pattern="(11)"

grep -v "$pattern" "$input_file" > "$output_file"