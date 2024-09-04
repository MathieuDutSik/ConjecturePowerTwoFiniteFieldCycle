CFILE = $(wildcard *.C)
OFILE = $(wildcard *.o)
PROG = $(patsubst %.o,%.prog, $(OFILE))
OBJS = $(patsubst %.C,%.o, $(CFILE))


%.o: %.C ; g++ -I/home/mathieu/opt/ntl-5.4/include $< -c -o $@

%.prog: %.o ; g++ $< /home/mathieu/opt/ntl-5.4/lib/libntl.a /home/mathieu/opt/gmp/lib/libgmp.a -o $@




all: $(PROG) $(OBJS)

clean:
	rm -f $(PROG) $(OBJS)
