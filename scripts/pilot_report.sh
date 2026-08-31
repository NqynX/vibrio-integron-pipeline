#!/bin/bash
JOB=${1:?usage: pilot_report.sh <array job id>}
P=$SCRATCH/vibrio-integron-pipeline
T=$(wc -l < $P/genome_list.txt)
grep -h '^TIMING' $P/logs/if_${JOB}_*.out | awk -v t=$T '
{split($3,a,"=");split($4,b,"=");split($5,c,"=");s+=a[2];m+=b[2];x+=c[2];n++
 if(a[2]>ms)ms=a[2]; if(b[2]>mm)mm=b[2]}
END{if(n==0){print "no TIMING lines yet";exit 1}
printf "pilot genomes   : %d\n",n
printf "mean/max time   : %.1f / %.1f min\n",s/n/60,ms/60
printf "mean/max output : %.1f / %d MB\n",m/n,mm
printf "mean attC       : %.1f\n",x/n
printf "\n--- extrapolated to %d ---\n",t
printf "CPU hours       : %.1f\n",s/n*t/3600
printf "wall @50 concur : %.1f h\n",s/n*t/3600/50
printf "TOTAL DISK      : %.1f GB\n",m/n*t/1024}'
