
name="Ntop_swei_enf4a"
python main.py $name
gnuplot -e "xi=240;yi=60" topomake.p
mv topology.eps ./res/${name}_topo.eps
mv history.log ./res/${name}_hist.log

name="Ntop_swei_enf4b"
python main.py $name
gnuplot -e "xi=240;yi=60" topomake.p
mv topology.eps ./res/${name}_topo.eps
mv history.log ./res/${name}_hist.log

name="Ntop_swei_enf4c"
python main.py $name
gnuplot -e "xi=240;yi=60" topomake.p
mv topology.eps ./res/${name}_topo.eps
mv history.log ./res/${name}_hist.log

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
