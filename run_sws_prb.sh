#!/bin/bash

nams="1a 1b 1c 1d 1e 1l 1f 1g 1i 1jj 1kk"
for n in $nams; do
    name="Ntop_swei_enf"$n
    python main.py $name
    gnuplot -e "xi=120;yi=60" topomake.p
    mv topology.eps ./self/${name}_topo.eps
    mv history.log ./self/${name}_hist.log
done
nams="3a 3b 3c 3d 3e 3l 3f 3g 3i 3jj 3kk"
for n in $nams; do
    name="Ntop_swei_enf"$n
    python main.py $name
    gnuplot -e "xi=120;yi=60" topomake.p
    mv topology.eps ./self/${name}_topo.eps
    mv history.log ./self/${name}_hist.log
done
nams="2a 2b 2c 2d 2e 2l 2f 2g 2i 2jj 2kk"
for n in $nams; do
    name="Ntop_swei_enf"$n
    python main.py $name
    gnuplot -e "xi=240;yi=60" topomake.p
    mv topology.eps ./self/${name}_topo.eps
    mv history.log ./self/${name}_hist.log
done
nams="4a 4b 4c 4d 4e 4l 4f 4g 4i 4jj 4kk"
for n in $nams; do
    name="Ntop_swei_enf"$n
    python main.py $name
    gnuplot -e "xi=240;yi=60" topomake.p
    mv topology.eps ./self/${name}_topo.eps
    mv history.log ./self/${name}_hist.log
done
