#!/bin/csh
set file=$1
set nl=`wc -l $file | cut -d " " -f7`
foreach num (`seq 1 1 $nl`)
	set a=`cut -d " " -f1 $file | sed -n "$num"p`
	set b=`cut -d " " -f3 $file | cut -d ";" -f1 | sed -n "$num"p`
	echo "$a = max(min($b,10^308),-10^308);" >> new_"$file"
end
