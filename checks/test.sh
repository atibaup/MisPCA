#!/bin/bash

functionList=$(cut functionNames.txt -d "/" -f 3,4 | cut -d "/" -f 2 | cut -d "." -f 1 )
functionList=${functionList// /|}

Script2Function=""
Invokedfunctions=""

k=1
cat scriptNames.txt | while IFS= read -r file;
do
	echo $k":"$file

	Script2Function="$Script2Function"$file

	for function in $functionList
	do
		found=$(grep $function"(" $file )
		if [ -n "$found" ]
		then
			Script2Function="$Script2Function"" "$function
			
			if [ -z "$(echo $Invokedfunctions |  grep $function )" ]
			then
				Invokedfunctions="$Invokedfunctions""-"$function
			fi
		fi
	done

	Script2Function="$Script2Function""\n"
	
	echo -e $Script2Function > Script2Function.txt
	echo -e $Invokedfunctions > Invokedfunctions.txt

	k=$(($k+1))
done

Function2Function=""

export Invokedfunctions=$(cat Invokedfunctions.txt)

cat functionNames.txt | while IFS= read -r file;
do
        echo $k":"$file

        Function2Function="$Function2Function"$file

        for function in $functionList
        do
                found=$(grep $function"(" $file )
                if [ -n "$found" ]
                then
                        Function2Function="$Function2Function"" "$function

                        if [ -z "$(echo $Invokedfunctions |  grep $function )" ]
                        then
                                Invokedfunctions="$Invokedfunctions""-"$function
                        fi
                fi
        done

        Function2Function="$Function2Function""\n"

        echo -e $Function2Function > Function2Function.txt

        echo -e ${Invokedfunctions//'-'/"\n"} > Invokedfunctions.txt

        k=$(($k+1))
done

# move non-invoked functions to the garbage/old directory

cat functionNames.txt | grep -vEf Invokedfunctions.txt

