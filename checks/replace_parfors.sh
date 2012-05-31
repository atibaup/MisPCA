!/bin/bash

# This script replaces the parallelizable "for"s by "parfors" to decrease running time on Matlab

dir="."
search="% Parallelizable for (deactivated for OCTAVE compatibility) \n for"
replace="parfor"
find ${dir} -type f -exec grep -l "${search}" {} \; | while read file
do
	sed "s/${search}/${replace}/ig" "$file" > tmp
	mv tmp "$file"
	echo "Modified: " $file	
done
