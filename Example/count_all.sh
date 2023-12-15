for b in Bottle*dat;  do
	echo $b
	perl ../../count_numbers.pl -i $b | sort -nk2 > $b.count
done

