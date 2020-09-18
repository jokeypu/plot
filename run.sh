let NN=300
export MAX_EVT=0
rmm /Users/jokey/Panda/workspace/plot/doc/Fit_par.txt 
for((i=0;i<NN;i++))
do
	MAX_EVT=$i
	cd /Users/jokey/Panda/workspace/plot
	expect run.expect 
done
