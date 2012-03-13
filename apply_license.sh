#!/bin/bash

git grep -H -c -i copyright | cut -d ':' -f1 > .done
ls *.cpp *.h > .cpp_h

cat .cpp_h | jtset -d .done > .todo

while read line; do
(echo '/*'; cat LICENSE; echo '*/'; cat $line) > $line.tmp
mv $line.tmp $line
done < .todo

#sed -ie 's/^/\/\/ /'

ls *.sh *.R Makefile > .sh_R_Makefile
cat .sh_R_Makefile | jtset -d .done > .todo

while read line; do
#sed -ie 's/^/# /' $line
done < .todo
