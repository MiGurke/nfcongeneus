#!/bin/bash

if [ "$1" == "-h" ]; then
  echo "Usage: `basename $0` -i /PATH/TO/INPUT/FOLDER/ -o /PATH/TO/OUTPUT/FOLDER/ -n number_raxml_jobs_started_in_parallel[INT] -t number_of_threads_used_by_each_raxml_job[INT] -m number_of_ml_searches_done_by_each_raxml_job[INT] "
  exit 0
fi

if [ "$1" == "-i" ]; then
  infolder=$2
else
  echo "Usage: `basename $0` -i /PATH/TO/INPUT/FOLDER/ -o /PATH/TO/OUTPUT/FOLDER/ -n number_raxml_jobs_started_in_parallel[INT] -t number_of_threads_used_by_each_raxml_job[INT] -m number_of_ml_searches_done_by_each_raxml_job[INT] "
  exit 0
fi

if [ "$3" == "-o" ]; then
  outfolder=$4
else
  echo "Usage: `basename $0` -i /PATH/TO/INPUT/FOLDER/ -o /PATH/TO/OUTPUT/FOLDER/ -n number_raxml_jobs_started_in_parallel[INT] -t number_of_threads_used_by_each_raxml_job[INT] -m number_of_ml_searches_done_by_each_raxml_job[INT] "
  exit 0
fi

if [ "$5" == "-n" ]; then
  N=$6
else
  echo "Usage: `basename $0` -i /PATH/TO/INPUT/FOLDER/ -o /PATH/TO/OUTPUT/FOLDER/ -n number_raxml_jobs_started_in_parallel[INT] -t number_of_threads_used_by_each_raxml_job[INT] -m number_of_ml_searches_done_by_each_raxml_job[INT] "
  exit 0
fi

if [ "$7" == "-t" ]; then
  threads=$8
else
  echo "Usage: `basename $0` -i /PATH/TO/INPUT/FOLDER/ -o /PATH/TO/OUTPUT/FOLDER/ -n number_raxml_jobs_started_in_parallel[INT] -t number_of_threads_used_by_each_raxml_job[INT] -m number_of_ml_searches_done_by_each_raxml_job[INT] "
  exit 0
fi

if [ "$9" == "-m" ]; then
  MLsearches=${10}
else
  echo "Usage: `basename $0` -i /PATH/TO/INPUT/FOLDER/ -o /PATH/TO/OUTPUT/FOLDER/ -n number_raxml_jobs_started_in_parallel[INT] -t number_of_threads_used_by_each_raxml_job[INT] -m number_of_ml_searches_done_by_each_raxml_job[INT] "
  exit 0
fi


cd $outfolder
c=0
for fasta in $infolder*.fasta; do
   ((i=i%N)); ((i++==0)) && wait
   sed -i 's/://g' $fasta 
   raxmlHPC-PTHREADS -m GTRGAMMA -T $threads -p 12345 -s $fasta -N $MLsearches -n $(basename $fasta .fasta) > /dev/null &
   c=$((c+1))
   echo "$c fasta files processed."
done
