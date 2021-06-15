local: a 

a:  
	./run.pbs 2>&1 | tee log.txt 

batch: b 

b: 
	sbatch run.pbs 

clean:
	rm -rf output/ *.tmp *.out *.tmp *.dat
