CXXFLAGS := $(shell root-config --cflags)
ROOTLIB := $(shell root-config --libs)
EIGEN_INC := $(shell pkg-config --cflags eigen3)


SRC_CPP := $(wildcard src/*.cpp)
SRC_C := $(wildcard src/*.C)

# OBJS := bin/DataPoints.o bin/Interpolator.o
OBJS := $(patsubst src/%.cpp, bin/%.o, $(SRC_CPP))
OBJS += $(patsubst src/%.C, bin/%.o, $(SRC_C))

.PRECIOUS: $(OBJS)

VPATH = main:src main/aulas

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
	clear