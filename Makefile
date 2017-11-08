#Makefile

# CXX = icc
# CXXFLAGS = -Wall -parallel -lm -fPIC

CC=gcc
all:ishikawa

ishikawa:graph.o bmp.o ford_fulkerson.o queue_stack.o iskw_drv.o ishikawa.o -lm
	$(CC) graph.o bmp.o ford_fulkerson.o queue_stack.o iskw_drv.o ishikawa.o -o ishikawa -lm -g3
.c.o:
	$(CC) -c $< -Wall -l -g3
remove:
	rm -f *.o
clean:
	rm -f *.o *~ ishikawa
