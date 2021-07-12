#!/bin/bash

rm -fr {0..10}

cp -r 0.original 0

blockMesh
moveDynamicMesh

cp -r 10/polyMesh/ 0
cp 10/polyMesh/points constant/polyMesh/
rm -rf {1..10}

superDeliciousVanilla


# **************************************************************************************************** #
