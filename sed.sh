nams="3a 3b 3c 3d 3e 3f 3g 3i 3j 3k 3jj 3kk"
for n in $nams; do
    sed -i "s/dext=0#int(np.ceil(rmin))/dext=int(np.ceil(rmin))/g" ./prob/pNtop_swei_enf${n}.py
    sed -i "s/dext=0/dext=int(np.ceil(rmin))/g" ./prob/pNtop_swei_enf${n}.py
done
nams="4a 4b 4c 4d 4e 4f 4g 4i 4j 4k 4jj 4kk"
for n in $nams; do
    sed -i "s/dext=0#int(np.ceil(rmin))/dext=int(np.ceil(rmin))/g" ./prob/pNtop_swei_enf${n}.py
    sed -i "s/dext=0/dext=int(np.ceil(rmin))/g" ./prob/pNtop_swei_enf${n}.py
done

