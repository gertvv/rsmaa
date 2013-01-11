all: smaa.so

clean:
	rm -f *.so *.o

.PHONY: all clean

%.so: %.c
	R CMD SHLIB -o $@ $<
