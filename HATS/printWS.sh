#!/bin/sh

now=$(date +%Y-%m-%d)
#ws_file_name = '/homs/HATS/data/aux/hats-'$now'.ws'
ws_file_name='/ogma/data/HATS/aux/hats-'$now'.ws'
last_line=$(tail -1 $ws_file_name)
time_stamp=$(echo $last_line | awk -F',' '{print $1}')
temperature=$(echo $last_line | awk -F',' '{print $3}' | awk -F'=' '{print $2}' | sed 's/C//')
humidity=$(echo $last_line | awk -F',' '{print $4}' | awk -F'=' '{print $2}' | sed 's/P//')
pressure=$(echo $last_line | awk -F',' '{print $5}' | awk -F'=' '{print $2}' | sed 's/H//')

echo $time_stamp , $temperature , $humidity , $pressure




