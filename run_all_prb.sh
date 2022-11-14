#
grep 'if' prbs.py | grep -o -P "(?<=').*(?=')" > prb_lst.log
#
rm prb_err.log 2> /dev/null
while read p; do
    echo "Problem: $p"
    echo "Running problem $p. The run is error free if no further output follows." >> prb_err.log
    python main.py $p 2>> prb_err.log
    cat history.log  | cut -c 1-71 > ./logs/"p${p}_hist.log"
done < prb_lst.log
#
