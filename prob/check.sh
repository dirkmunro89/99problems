nams="a b c d e f g h i j k l m n o p q r s t u v w x y z"
for n in $nams; do
	diff pNtop_swei3_$n.py pNtop_swei4_${n}.py
done
