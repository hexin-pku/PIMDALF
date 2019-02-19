#!/bin/bash

parm_file=./s.par
exe=./mespimd
python3 parm_conv.py ${parm_file}

while read var value
do
    export "$var"="$value"
    echo $var $value
done < ${parm_file}flx

freedom=1
echo freedom ${freedom}
echo

# reach a equilibrium
${exe} -p ${parm_file} -s e -o 0
rm 0.ana 0.trj time.dat list.dat

for(( i=1;i<100;i++ ))
do
    ${exe} -p ${parm_file} -s r -o ${i} >> time.dat
    rm ${i}.trj
    python3 process.py ${i} ${bead} ${freedom} >> list.dat
done

python3 stat_list.py > result.dat
