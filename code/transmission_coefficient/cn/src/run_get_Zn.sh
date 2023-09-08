#!/bin/bash
for config in {1..100}
do
  for vel in {1..5}
  do
	if [ -d config_${config}_vel_${vel} ]
	then
		cd config_${config}_vel_${vel}
    		echo " working on config_${config}_vel_${vel} "
    		if [ -f md-1.xyz ]
    		then
			cp ../toZn2.py .
			sed '/^     323/,+1 d' md-1.xyz > temp_traj.xyz
			python toZn2.py > Zn_timeseries2.dat
			rm temp_traj.xyz
		fi
		cd ..
	fi
   done
done
