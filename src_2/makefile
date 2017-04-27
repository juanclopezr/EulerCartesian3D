CC = gcc
CC_OPTIONS = -g -O3
EXEC = euler.x
LIBS = -lm
SOD = 0
SEDOV = 1

OBJS = struct.o init.o solver.o io.o main.o problem.o tube.o
INCL = Makefile struct.h io.h init.h solver.h problem.h tube.h

all: $(EXEC)

plotsedov: #*_sedov.dat
	python plotter.py $(SEDOV)

plotshock: tube.dat
	python plotter.py $(SOD)

tube.dat: $(EXEC)
	./$(EXEC) $(SOD)

*_sedov.dat: $(EXEC)
	./$(EXEC) $(SEDOV)

shock: $(EXEC) tube.dat

sedov: $(EXEC) %_sedov.dat

exec: $(EXEC)

.c.o:
	$(CC) $(CC_OPTIONS) -c $<

$(EXEC): $(OBJS)
	$(CC) $(OBJS) $(LIBS) -o $(EXEC)
	rm $(OBJS)

clean:
	rm -f $(OBJS) *~ core* ${EXEC} *.dat
