
# nohup bash run.sh > nohup.out 2>&1 &

#rm nohup.out
#touch nohup.out
rm verbose.out
nbruns=1
for i in $(seq 1 1 $nbruns)
do
  echo "Run $i/$nbruns"
  TZ=":Indian/Reunion" date +'%a, %b %d, %Y  %r'
  nohup Rscript models.r --no-save > verbose.out 2>&1 || true & # the || true is to keep going if there is an error
  wait $!
  e=$(TZ=":Indian/Reunion" date +'%a, %b %d, %Y  %r')
  echo "$e -- done"
  #ls -1tlh /data/sdm_rbq/sdms | head -5
  #echo "Run $i/$nbruns done"
done  


