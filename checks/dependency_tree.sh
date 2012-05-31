#!/bin/bash

# Find script names and store in file.
ScriptDirs="Experiments TestScripts"

for dir in $ScriptDirs
do 
	find $dir -type f -name "*.m" | xargs  grep -H -v -n "function" | grep ":1:" | cut -d ':' -f 1 > scripts.$dir.txt
	echo $dir
done

cat $(ls scripts.*) > scriptNames.txt
rm scripts.*

# find function names
FunctionDirs="Utilities Experiments"
for dir in $FunctionDirs
do
        find $dir -type f -name "*.m" | xargs  grep -H -n "function" | grep ":1:" | cut -d ':' -f 1 > functions.$dir.txt
        echo $dir
done

cat $(ls functions.*) > functionNames.txt
rm functions.*
# for each script, find invoked functions

# for each function, find their children


#dir="."
#search="% Parallelizable for (deactivated for OCTAVE compatibility) \n for"
#replace="parfor"
#find ${dir} -type f -exec grep -l "${search}" {} \; | while read file
#do
#        sed "s/${search}/${replace}/ig" "$file" > tmp
#        mv tmp "$file"
#        echo "Modified: " $file 
#done

