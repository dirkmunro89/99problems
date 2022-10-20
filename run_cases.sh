
name="Ntop_swei_enf1a"
python main.py $name
gnuplot -e "xi=120;yi=60" topomake.p
mv topology.eps ./rsr/${name}_topo.eps
mv history.log ./rsr/${name}_hist.log

name="Ntop_swei_enf1b"
python main.py $name
gnuplot -e "xi=120;yi=60" topomake.p
mv topology.eps ./rsr/${name}_topo.eps
mv history.log ./rsr/${name}_hist.log

name="Ntop_swei_enf1c"
python main.py $name
gnuplot -e "xi=120;yi=60" topomake.p
mv topology.eps ./rsr/${name}_topo.eps
mv history.log ./rsr/${name}_hist.log

name="Ntop_swei_enf3a"
python main.py $name
gnuplot -e "xi=120;yi=60" topomake.p
mv topology.eps ./rsr/${name}_topo.eps
mv history.log ./rsr/${name}_hist.log

name="Ntop_swei_enf3b"
python main.py $name
gnuplot -e "xi=120;yi=60" topomake.p
mv topology.eps ./rsr/${name}_topo.eps
mv history.log ./rsr/${name}_hist.log

name="Ntop_swei_enf3c"
python main.py $name
gnuplot -e "xi=120;yi=60" topomake.p
mv topology.eps ./rsr/${name}_topo.eps
mv history.log ./rsr/${name}_hist.log

#name="Ntop_swei_enf2a"
#python main.py $name
#gnuplot -e "xi=240;yi=60" topomake.p
#mv topology.eps ./res/${name}_topo.eps
#mv history.log ./res/${name}_hist.log

#name="Ntop_swei_enf2b"
#python main.py $name
#gnuplot -e "xi=240;yi=60" topomake.p
#mv topology.eps ./res/${name}_topo.eps
#mv history.log ./res/${name}_hist.log

#name="Ntop_swei_enf2c"
#python main.py $name
#gnuplot -e "xi=240;yi=60" topomake.p
#mv topology.eps ./res/${name}_topo.eps
#mv history.log ./res/${name}_hist.log

#name="Ntop_swei_enf3a"
#python main.py $name
#gnuplot -e "xi=240;yi=60" topomake.p
#mv topology.eps ./res/${name}_topo.eps
#mv history.log ./res/${name}_hist.log

#name="Ntop_swei_enf3b"
#python main.py $name
#gnuplot -e "xi=240;yi=60" topomake.p
#mv topology.eps ./res/${name}_topo.eps
#mv history.log ./res/${name}_hist.log

#name="Ntop_swei_enf3c"
#python main.py $name
#gnuplot -e "xi=240;yi=60" topomake.p
#mv topology.eps ./res/${name}_topo.eps
#mv history.log ./res/${name}_hist.log
