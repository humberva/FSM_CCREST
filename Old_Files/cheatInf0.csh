#!/bin/csh
set file=$1
set file_name=`echo $file`
set fl=`echo $file | cut -c1`
set n=`wc -l "$file" | tr -d "[ ]" | cut -d "$fl" -f1`
sed -e 's/1\.0/1/g' < "$file" > temp
tr -d "[.]" < temp > new_"$file".txt
foreach num (`seq 1 1 $n`)
	set mterms=`sed -n "$num"p new_"$file".txt | tr -dc "[*]" | wc -c`
	if ($mterms > 1) then
		set cvar=`sed -n "$num"p new_"$file".txt | cut -d "=" -f1`
		set a=`sed -n "$num"p new_"$file".txt | cut -d "=" -f2 | cut -d "*" -f1 | tr -d "[;]"`
		set b=`sed -n "$num"p new_"$file".txt | cut -d "=" -f2 | cut -d "*" -f2 | tr -d "[;]"`
		set term="$cvar"_1
		echo "$term = nansum([$a * $b, 0]);" >> good_"$file"
		foreach tn (`seq 2 1 $mterms`)
			set newterm="$cvar"_"$tn"
			set tn_1=`expr $tn + 1`
			set b=`sed -n "$num"p new_"$file".txt | cut -d "=" -f2 | cut -d "*" -f"$tn_1" | tr -d "[;]"`
			echo "$newterm = nansum([$term * $b,0]);" >> good_"$file"
			set term=$newterm
		end
		sed -e s/$newterm/$cvar/g < good_"$file" > temp
		mv temp good_"$file"
	endif

	if ($mterms == 1) then
		set cvar=`sed -n "$num"p new_"$file".txt | cut -d "=" -f1`
                set a=`sed -n "$num"p new_"$file".txt | cut -d "=" -f2 | cut -d "*" -f1 | tr -d "[;]"`
                set b=`sed -n "$num"p new_"$file".txt | cut -d "=" -f2 | cut -d "*" -f2 | tr -d "[;]"`
                echo "$cvar = nansum([$a * $b, 0]);" >> good_"$file"
        endif

	if ($mterms == 0) then
		sed -n "$num"p new_"$file".txt >> good_"$file"
	endif
end
