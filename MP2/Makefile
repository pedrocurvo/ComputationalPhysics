CXXFLAGS := $(shell root-config --cflags)
ROOTLIB := $(shell root-config --libs)
EIGEN_INC := $(shell pkg-config --cflags eigen3)

OBJS := bin/lightmap.o bin/MP2s.o

.PRECIOUS: $(OBJS)

VPATH = main:src

#just compilling
bin/%.o: %.cpp
	g++ -c $(CXXFLAGS) $^ -o $@ -I src $(EIGEN_INC)

# just linking 
bin/%.exe : bin/%.o $(OBJS)
	@echo "Compiling"
	g++ $(CXXFLAGS) $^ -o $@ -I src $(ROOTLIB) $(EIGEN_INC)

bin/%.exe : %.C
	@echo "Compiling"
	g++ $(CXXFLAGS) $^ -o $@ -I src $(ROOTLIB) $(EIGEN_INC)

#regra de limpeza
clean :
	rm -f bin/*.o bin/*.exe