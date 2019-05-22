#!/bin/sh

echo $PWD
echo "==>Arguments read:"
echo $@

echo $PWD > module_log.txt
echo "==>Arguments read:" >> module_log.txt
echo $@ >> module_log.txt