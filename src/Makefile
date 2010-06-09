# ARA* + AD* author: Leonardo Moura ##
############################# Makefile ##########################
CC=g++
FLAGS=-Wno-deprecated -lglut -lGLU -lGL -lm
BOOST=/usr/local/boost_1_39_0
INCLUDE=-I $(BOOST)
OPENGL=-lglut -lGLU -lGL -lm
all: main.o ara.o ad.o
	  $(CC) $(FLAGS) -o ARA main.o ad.o ara.o drawGraph.o
	  ./ARA
#heap.o: heap.h
#	  $(CC) -o heap.o heap.h -W -Wall -ansi -pedantic
#PriorityQueue.o:	PriorityQueue.h heap.o
#	  $(CC) -o PriorityQueue.o -combine heap.o PriorityQueue.h -W
drawGraph.o:	drawGraph.h drawGraph.cpp
	  $(CC) $(FLAGS) $(INCLUDE) $(OPENGL) -c drawGraph.cpp
ara.o:	Ara.cpp Ara.h PriorityQueue.h
	  $(CC) $(FLAGS) $(INCLUDE) -o ara.o -c Ara.cpp 
ad.o:	Ad.cpp Ad.h PriorityQueue.h drawGraph.o
	  $(CC) $(FLAGS) $(INCLUDE) -o ad.o -c Ad.cpp
main.o: main.cpp ad.o ara.o
	  $(CC) $(FLAGS) $(INCLUDE) -o main.o -c main.cpp
clean:
	  rm -rf *.o
mrproper: clean
	  rm -rf all
