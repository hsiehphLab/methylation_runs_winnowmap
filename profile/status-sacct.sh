#!/usr/bin/env bash

# Check status of Slurm job

jobid="$1"

if [[ "$jobid" == Submitted ]]
then
  echo smk-simple-slurm: Invalid job ID: "$jobid" >&2
  echo smk-simple-slurm: Did you remember to add the flag --parsable to your sbatch call? >&2
  exit 1
fi

output=`sacct -j "$jobid" --allclusters --format State --noheader | head -n 1 | awk '{print $1}'`

# this script appears to be run in the top level fastcn directory--not in profile
echo $jobid $output >>sacct.log


# $output == "" added DG since I've seen cases in which the job is
# successfully running but sacct returns nothing.  This script used to
# return failed in such cases.  So snakemake would resubmit the job,
# sometimes causing collisions with the original running job.


if [[ $output =~ ^(COMPLETED).* ]]
then
  echo success
elif [[ $output =~ ^(RUNNING|PENDING|COMPLETING|CONFIGURING|SUSPENDED|REQUEUED).* ]]
then
  echo running
elif [ "$output" = "" ]
then
    echo running
else
  echo failed
fi
