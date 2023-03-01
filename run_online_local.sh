#!/bin/bash

logext="experiment_online"

d=10

for tag in $(seq 1 1 1)
do
	for bandwidth in 500
	do
		python3 throttle.py -i lo -b $bandwidth -l 1

	  for n in $(seq 80 16 80)
	  do
	    rm client_info.txt
	    touch client_info.txt
	    port=12345
	    for i in $(seq 0 1 $(($n - 1)))
	    do
	      echo "$i,127.0.0.1,$(($port + $i))" | tee -a client_info.txt
	    done

	    for s in 100000
	    do

	      for percentage in 0.6 0.7 0.8 0.9
	      do
		
		echo "running: n=$n i s=$s d=$d corruption=$percentage"
		logdir="logs/logs_${logext}_${n}_${s}_${d}_${percentage}_${bandwidth}_${tag}"

		mkdir -p $logdir

		( ./build/ours_online.x $n 0 $s $d $percentage | tee "${logdir}/party_0.log" ) &

		for i in $(seq 1 1 $(($n - 1))); do
		    ( ./build/ours_online.x $n $i $s $d $percentage &>"${logdir}/party_${i}.log" ) &
		    pid=$!
		done

		echo "waiting for OUR experiment to finish ..."
		wait

		sleep 1

	      done
	    done
	  done
  done
done
