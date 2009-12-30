PROJLIBFILE=simpmesh.c tettracing.c mcx_utils.c
OBJFILE =simpmesh.o tettracing.o
GCC     =gcc
PROJFILE =tetray
OPT    = -g -lm #-O3
#OPT    = -O3

.o:
	$(GCC) -c $(OPT) $(PROJLIBFILE)
all:
	$(GCC) -c $(OPT) $(PROJLIBFILE)
	$(GCC)    $(OPT) $(PROJFILE).c $(OBJFILE) -o $(PROJFILE)
clean:
	rm -f *.o $(PROJFILE)
install:
	cp -a $(PROJFILE) ~/bin
