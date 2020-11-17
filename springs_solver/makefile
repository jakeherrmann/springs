CC := g++
CC_FLAGS := -Wall -Wno-int-in-bool-context -std=gnu++17 -O3 -msse -msse2 -msse3 -mfpmath=sse -mtune=native -march=native
# -fopenmp
CPPFILES := $(wildcard *.cpp)
OBJFILES := $(addprefix obj/,$(CPPFILES:.cpp=.o))

springs_solver.exe: $(OBJFILES)
	$(CC) $(CC_FLAGS) $^ -o $@
	
obj/%.o: %.cpp
	$(CC) $(CC_FLAGS) -c $< -o $@
	
clean:
	rm -f obj/*.o

