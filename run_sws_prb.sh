#!/bin/bash

nams="3a 3b 3c 3d 3e 3f 3g"
nams="3b"
for n in $nams; do
    name="Ntop_swei_enf"$n
    python main.py $name
    gnuplot -e "xi=120;yi=60" topomake.p
    mv topology.eps ./resu/${name}_topo.eps
    mv history.log ./resu/${name}_hist.log
done
nams="4a 4b 4c 4d 4e 4f 4g"
nams="4b"
for n in $nams; do
    name="Ntop_swei_enf"$n
    python main.py $name
    gnuplot -e "xi=240;yi=60" topomake.p
    mv topology.eps ./resu/${name}_topo.eps
    mv history.log ./resu/${name}_hist.log
done
nams="1a 1b 1c 1d 1e 1f 1g"
nams="1b"
for n in $nams; do
    name="Ntop_swei_enf"$n
    python main.py $name
    gnuplot -e "xi=120;yi=60" topomake.p
    mv topology.eps ./resu/${name}_topo.eps
    mv history.log ./resu/${name}_hist.log
done
nams="2a 2b 2c 2d 2e 2f 2g"
nams="2b"
for n in $nams; do
    name="Ntop_swei_enf"$n
    python main.py $name
    gnuplot -e "xi=240;yi=60" topomake.p
    mv topology.eps ./resu/${name}_topo.eps
    mv history.log ./resu/${name}_hist.log
done
