#!/bin/bash

start_time=`date +%s`
./ishikawa input/barbara_noise32.bmp /dev/null 16 >> result.txt
./ishikawa input/barbara_noise32.bmp /dev/null 32 >> result.txt
./ishikawa input/barbara_noise32.bmp /dev/null 64 >> result.txt
end_time=`date +%s`
time=$((end_time -start_time));

echo "@zurvlk 処理が完了しましたっ! 所要時間[${time}s]" | tw --pipe --user="zurvlk"
