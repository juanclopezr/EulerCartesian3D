all : tube_plot.py tube.dat
	python tube_plot.py

tube.dat : a.out
	./a.out

a.out : tube.c
	gcc -lm -g tube.c

animate : tube_plot.py tube.dat
	python tube_plot.py animate 

clean :
	rm *.dat a.out
