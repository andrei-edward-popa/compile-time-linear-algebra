CXX=g++
EXE=test
LD_FLAGS=-lfmt -lgtest
INCLUDE=-Iinclude
OBJ_FILES=$(CPP_FILES:.cpp=.o)
CPP_FILES=$(shell find src/ -name "*.cpp")
D_FILES=$(CPP_FILES:.cpp=.d)
CXX_FLAGS=-DPRECISION=${PRECISION} -std=c++2b -O3 -Wall -Wextra -Werror -Wpedantic -Wdeprecated -Wconversion -flto=auto

$(EXE): $(OBJ_FILES)
	$(CXX) $^ $(INCLUDE) $(CXX_FLAGS) -o $@ $(LD_FLAGS)
	strip $(EXE)

src/%.o: src/%.cpp
	$(CXX) -MMD -c $< $(INCLUDE) $(CXX_FLAGS) -o $@ $(LD_FLAGS)

.PHONY: mrproper
mrproper:
	rm -f test
	find src/ -name "*.o" | xargs rm -f
	find src/ -name "*.d" | xargs rm -f
	
-include $(D_FILES)

ifndef PRECISION
PRECISION=4
endif

