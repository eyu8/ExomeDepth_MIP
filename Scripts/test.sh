#!/bin/bash


sleep 10 | echo "10" &
sleep 15 | echo "15" &

wait $(jobs -p)
echo "done"
