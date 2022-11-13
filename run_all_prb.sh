#
grep 'if' prbs.py | grep -o -P "(?<=').*(?=')" > prb_lst.log
#
rm prb_err.log 2> /dev/null
while read p; do
    echo "Problem: $p"
    python main.py $p 2>> prb_err.log
    cp history.log ./logs/"P${p}_hist.log"
done < prb_lst.log
#
