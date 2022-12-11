#!/bin/bash

nams="3a 3c 3f 3g"
for n in $nams; do
    name="Ntop_swei_enf"$n
    python main.py $name
    gnuplot -e "xi=120;yi=60" topomake.p
    mv topology.eps ./tesu/${name}_topo.eps
    mv history.log ./tesu/${name}_hist.log
done
nams="2a 2c 2f 2g"
for n in $nams; do
    name="Ntop_swei_enf"$n
    python main.py $name
    gnuplot -e "xi=240;yi=60" topomake.p
    mv topology.eps ./tesu/${name}_topo.eps
    mv history.log ./tesu/${name}_hist.log
done
nams="1a 1c 1f 1g"
for n in $nams; do
    name="Ntop_swei_enf"$n
    python main.py $name
    gnuplot -e "xi=120;yi=60" topomake.p
    mv topology.eps ./tesu/${name}_topo.eps
    mv history.log ./tesu/${name}_hist.log
done
nams="4a 4c 4f 4g"
for n in $nams; do
    name="Ntop_swei_enf"$n
    python main.py $name
    gnuplot -e "xi=240;yi=60" topomake.p
    mv topology.eps ./tesu/${name}_topo.eps
    mv history.log ./tesu/${name}_hist.log
done

