grep 'total time:'|cut -d' ' -f2,5|sed -e 's/s//'|sort -n|cut -d' ' -f2
