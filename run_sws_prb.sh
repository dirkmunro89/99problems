#!/bin/bash
nams="2a 2b 2c 2d 2e 2f 2g"
for n in $nams; do
    name="Ntop_swei_enf"$n
    python main.py $name
    gnuplot -e "xi=240;yi=60" topomake.p
    mv topology.eps ./resu/${name}_topo.eps
    mv history.log ./resu/${name}_hist.log
done
nams="1a 1b 1c 1d 1e 1f 1g"
for n in $nams; do
    name="Ntop_swei_enf"$n
    python main.py $name
    gnuplot -e "xi=120;yi=60" topomake.p
    mv topology.eps ./resu/${name}_topo.eps
    mv history.log ./resu/${name}_hist.log
done
nams="1aa 1bb 1cc 1dd 1ee 1ff 1gg"
for n in $nams; do
    name="Ntop_swei_enf"$n
    python main.py $name
    gnuplot -e "xi=120;yi=60" topomake.p
    mv topology.eps ./resu/${name}_topo.eps
    mv history.log ./resu/${name}_hist.log
done
nams="2aa 2bb 2cc 2dd 2ee 2ff 2gg"
for n in $nams; do
    name="Ntop_swei_enf"$n
    python main.py $name
    gnuplot -e "xi=240;yi=60" topomake.p
    mv topology.eps ./resu/${name}_topo.eps
    mv history.log ./resu/${name}_hist.log
done
