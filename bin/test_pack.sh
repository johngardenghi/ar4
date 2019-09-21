#!/bin/bash

for i in {1000,2000,4000,6000,8000,10000}
do
    time=$(date)
    echo $time "- Running packing problem with N =" $i
    ./pack -N $i -s arp -p 2 -pl 2 > $i"_arp_p2.out"
done