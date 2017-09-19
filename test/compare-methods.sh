#!/bin/bash

if [ $# -ne 3 ]
then
    echo " Usage : $0   dataset   method1   method2"
    echo "     possible methods are :"
    echo "      lsape_bunke,  lsape_rw"
    echo "      lsape_multi_bunke,   lsape_multi_rw,   lsape_multi_greedy"
    echo "      multi_random"
    echo 
    echo "      ipfpe_bunke, ipfpe_rw"
    echo "      ipfpe_multi_bunke,   ipfpe_multi_rw,   ipfpe_multi_random,   ipfpe_multi_greedy"
    echo "      gnccp"
    echo
    echo " For the 'multi' methods, give the number of editpaths k to be considered as : "
    echo "     $0   ...   \"method -p k\" ... "
    exit
fi


dataset=$1
method1=$2
method2=$3

paste <(./chemical-edit-distances $dataset -s -m $method1) <(./chemical-edit-distances $dataset -s -m $method2)
