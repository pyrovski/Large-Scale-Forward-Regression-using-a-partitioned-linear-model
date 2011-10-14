grep commu | cut -d' ' -f2,6 | sed -e 's/s//' | sort -n | cut -d' ' -f2
