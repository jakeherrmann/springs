CC := g++
# CC_FLAGS := -Wall -Wno-int-in-bool-context -std=gnu++17 -O3 -msse -msse2 -msse3 -mfpmath=sse -mtune=native -march=native #single-threaded
CC_FLAGS := -Wall -Wno-int-in-bool-context -std=gnu++2a -O3 -msse -msse2 -msse3 -mfpmath=sse -mtune=native -march=native -fopenmp -mvzeroupper -g #multi-threaded
CPPFILES := $(wildcard *.cpp)
OBJDIR := obj
OBJFILES := $(addprefix $(OBJDIR)/,$(CPPFILES:.cpp=.o))
LIB := -lmetis

dir_guard=@mkdir -p $(@D)

springs_solver.exe: $(OBJFILES)
	$(CC) $(CC_FLAGS) $^ -o $@ $(LIB) 
	
$(OBJDIR)/%.o: %.cpp
	$(dir_guard)
	$(CC) $(CC_FLAGS) -c $< -o $@
	
clean:
	rm -f $(OBJDIR)/*.o

