#!/bin/bash
d=$1
ip1=$2
ip2=$3
group=$4

usage () {
    echo $@
    echo "usage: $0 [N] [size] [depth] [ip1] [ip2]"
    exit 0
}

#if ! [[ $n =~ ^[0-9]+$ ]]; then
    #usage "not a number:" $n
#fi

logext="experiment"

for tag in $(seq 1 1 5)
do
  for n in $(seq 16 16 32)
  do
    rm client_info.txt
    touch client_info.txt
    port=12345
    for i in $(seq 0 1 $(($n / 2 - 1)))
    do
      echo "$i,$ip1,$(($port + $i))" | tee -a client_info.txt
    done
    for i in $(seq $(($n / 2)) 1 $(($n - 1)))
    do
      echo "$i,$ip2,$(($port + $i))" | tee -a client_info.txt
    done

    for s in 100 1000
    do

      for percentage in 0.6 0.7 0.8 0.9
      do
        
        echo "running: n=$n i s=$s d=$d corruption=$percentage"
        logdir="logs/logs_${logext}_${n}_${s}_${d}_${percentage}_${tag}"

        mkdir -p $logdir

        if [ $group -eq 0 ]
        then
          ( ./build/ours.x $n 0 $s $d $percentage | tee "${logdir}/party_0.log" ) &

          for i in $(seq 1 1 $(($n / 2 - 1))); do
              ( ./build/ours.x $n $i $s $d $percentage &>"${logdir}/party_${i}.log" ) &
              pid=$!
          done
          echo
        else
          ( ./build/ours.x $n $(($n / 2)) $s $d $percentage | tee "${logdir}/party_0.log" ) &

          for i in $(seq $(($n / 2 + 1)) 1 $(($n - 1))); do
              ( ./build/ours.x $n $i $s $d $percentage &>"${logdir}/party_${i}.log" ) &
              pid=$!
          done
          echo
        fi

        echo "waiting for OUR experiment to finish ..."
        wait

      done
    done
  done
done
