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

