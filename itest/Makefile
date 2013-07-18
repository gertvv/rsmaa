all: hansen.Rout

.PHONY: all

hansen.Rout: hansen.R\
		hansen-weights-ord.txt hansen-meas.txt\
		hansen-ra-ord.txt hansen-cf-ord.txt 
	R --vanilla --slave -f $< > $@ 2>&1

%.txt: %.txt.gz
	gunzip -c $< > $@
