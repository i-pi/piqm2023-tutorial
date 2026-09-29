#!/bin/bash

git clean -f .

cp ../pimd-therm/simulation.restart_11 therm_checkpoint.chk

python -u $(which i-pi) rpmd.xml &> log.ipi &

sleep 10

for i in {1..32} ; do 
    i-pi-driver -u -a oh-rpmd -m morsedia &> driver.out &
done
