


cat verbose.out | tail -10
free -m | awk '/Mem/{print $3/$2}'
grep -E 'Approx inference \(total\)' verbose.out
grep -E 'Starting' verbose.out | sort
grep -E 'Running species' verbose.out | sort
cat nohup.out