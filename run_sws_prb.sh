#!/bin/bash

#nams="a b c d e f g h i j k l m n o p q"
#for n in $nams; do
#    name="Ntop_swei1_"$n
#    python main.py $name
#    gnuplot -e "xi=120;yi=60" topomake.p
#    mv topology.eps ./self/${name}_topo.eps
#    mv history.log ./self/${name}_hist.log
#done
#nams="a b c d e f g h i j k l m n o p q"
#for n in $nams; do
#    name="Ntop_swei2_"$n
#    python main.py $name
#    gnuplot -e "xi=120;yi=60" topomake.p
#    mv topology.eps ./self/${name}_topo.eps
#    mv history.log ./self/${name}_hist.log
#done
nams="a b c d e f g h i j k l m n o p q"
for n in $nams; do
    name="Ntop_swei3_"$n
    python main.py $name
    gnuplot -e "xi=240;yi=60" topomake.p
    mv topology.eps ./self/${name}_topo.eps
    mv history.log ./self/${name}_hist.log
done
nams="a b c d e f g h i j k l m n o p q"
for n in $nams; do
    name="Ntop_swei4_"$n
    python main.py $name
    gnuplot -e "xi=240;yi=60" topomake.p
    mv topology.eps ./self/${name}_topo.eps
    mv history.log ./self/${name}_hist.log
done
