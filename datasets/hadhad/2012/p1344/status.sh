#!/bin/bash

echo ""
echo "MC:"
echo "======================="
ami list datasets -F -L --fields prodsys_status,events mc.txt

echo ""
echo "DATA:"
echo "======================="
ami list datasets -F -L --fields prodsys_status,events data.txt
