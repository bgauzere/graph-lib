#!/bin/sh

./xp_PRL2017 ~/work/Datasets/Acyclic/dataset_bps.ds  -o acyclic_metric -s
./xp_PRL2017 ~/work/Datasets/MAO/dataset.ds  -o mao_metric -s

./xp_PRL2017 ~/work/Datasets/Acyclic/dataset_bps.ds  -o acyclic_nonmetric -s -c 6,2,2,3,1,1
./xp_PRL2017 ~/work/Datasets/MAO/dataset.ds  -o mao_nonmetric -s -c 6,2,2,3,1,1
