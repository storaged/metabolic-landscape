solutions_dir=$1
paste $1/program*.sol | grep "^v[0-9]*" | awk '{for(x=1;x<=NF;x++)if(x % 2 == 0)printf "%s", $x (x == NF || x == (NF-1)?"\n":" ")}' > tmp.tmp
paste $(ls $1/program*.sol | head -n 1 ) | grep -o  "^v[0-9]*" > colnames.tmp
ls $1/program*.sol | grep -o program.* | sed 's/\..*//' | tr '\n' ' ' > final.tmp ; echo >> final.tmp ; paste -d' ' colnames.tmp tmp.tmp >> final.tmp
rm tmp.tmp colnames.tmp 
