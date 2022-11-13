
name="Ntop_swei_enf2a"
python main.py $name
gnuplot -e "xi=240;yi=60" topomake.p
mv topology.eps ./resu/${name}_topo.eps
mv history.log ./resu/${name}_hist.log

name="Ntop_swei_enf2b"
python main.py $name
gnuplot -e "xi=240;yi=60" topomake.p
mv topology.eps ./resu/${name}_topo.eps
mv history.log ./resu/${name}_hist.log

name="Ntop_swei_enf2c"
python main.py $name
gnuplot -e "xi=240;yi=60" topomake.p
mv topology.eps ./resu/${name}_topo.eps
mv history.log ./resu/${name}_hist.log

name="Ntop_swei_enf2d"
python main.py $name
gnuplot -e "xi=240;yi=60" topomake.p
mv topology.eps ./resu/${name}_topo.eps
mv history.log ./resu/${name}_hist.log

name="Ntop_swei_enf2e"
python main.py $name
gnuplot -e "xi=240;yi=60" topomake.p
mv topology.eps ./resu/${name}_topo.eps
mv history.log ./resu/${name}_hist.log

name="Ntop_swei_enf1a"
python main.py $name
gnuplot -e "xi=120;yi=60" topomake.p
mv topology.eps ./resu/${name}_topo.eps
mv history.log ./resu/${name}_hist.log

name="Ntop_swei_enf1b"
python main.py $name
gnuplot -e "xi=120;yi=60" topomake.p
mv topology.eps ./resu/${name}_topo.eps
mv history.log ./resu/${name}_hist.log

name="Ntop_swei_enf1c"
python main.py $name
gnuplot -e "xi=120;yi=60" topomake.p
mv topology.eps ./resu/${name}_topo.eps
mv history.log ./resu/${name}_hist.log

name="Ntop_swei_enf1d"
python main.py $name
gnuplot -e "xi=120;yi=60" topomake.p
mv topology.eps ./resu/${name}_topo.eps
mv history.log ./resu/${name}_hist.log

name="Ntop_swei_enf1e"
python main.py $name
gnuplot -e "xi=120;yi=60" topomake.p
mv topology.eps ./resu/${name}_topo.eps
mv history.log ./resu/${name}_hist.log


