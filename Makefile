SRC=$(shell ls src/*.cpp)
INC=$(shell ls src/*.h)

isoplotter: $(SRC) $(INC)
	g++ -g -std=c++11 $(SRC) -o isoplotter -lseqio

clean:
	rm isoplotter