CC	      = g++
CC_FLAGS        = -g3 -O3 -Wall -D_GLIBCXX_DEBUG -I  /opt/homebrew/Cellar/gsl/2.7.1/include/
LD_FLAGS        = -L/opt/homebrew/Cellar/gsl/2.7.1/lib  -lgsl -lgslcblas -lm -lstdc++ 
SAM		= diffusion.o absorbing.o reflecting.o io.o utilities.o data.o analysis.o

diff: $(SAM)
	$(CC) $(CC_FLAGS) $(SAM) -o run_diffusion  $(LD_FLAGS)
diffusion.o: diffusion.cpp
	$(CC) $(CC_FLAGS) -c diffusion.cpp
absorbing.o: absorbing.cpp
	$(CC) $(CC_FLAGS) -c absorbing.cpp
reflecting.o: reflecting.cpp
	$(CC) $(CC_FLAGS) -c reflecting.cpp
data.o: data.cpp
	$(CC) $(CC_FLAGS) -c data.cpp
io.o: io.cpp
	$(CC) $(CC_FLAGS) -c io.cpp
utilities.o: utilities.cpp
	$(CC) $(CC_FLAGS) -c utilities.cpp
analysis.o: analysis.cpp
	$(CC) $(CC_FLAGS) -c analysis.cpp

