#!/bin/bash

start_time=`date +%s`
count=1
for label_size in 16 32 64
do
    job_start=`date +%s`
    ./rs input/barbara_noise32.bmp /dev/null $label_size >> result.txt
    job_end=`date +%s`
    time=$((job_end - job_start));
    echo "@zurvlk 処理${count}が完了しましたっ! 所要時間[${time}s]" | tw --pipe --user="zurvlk"
    count=$count+1
done
end_time=`date +%s`
time=$((end_time -start_time));

echo "@zurvlk 全ての処理が完了しましたっ! 総所要時間[${time}s]" | tw --pipe --user="zurvlk"
