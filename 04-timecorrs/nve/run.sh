#!/bin/bash

git clean -f .

cp ../class-therm/simulation.restart_11 therm.chk

python -u $(which i-pi) nve.xml &> log.ipi &

sleep 10

i-pi-driver -u -a oh-nve -m morsedia &> driver.out &
