SUBDIRS := $(wildcard */.)

local: a 

a:  
	for dir in $(SUBDIRS); do \
		cd $$dir; \
		./run.pbs; \
		cd ..; \
	done

batch: b 

b: 
	for dir in $(SUBDIRS); do \
		cd $$dir; \
		sbatch run.pbs; \
		cd ..; \
	done

clean:
	rm -rf output/ *.tmp *.out *.tmp *.dat
