CXX=g++
EXE=test
LD_FLAGS=-lfmt
INCLUDE=-Iinclude
OBJ_FILES=$(CPP_FILES:.cpp=.o)
CPP_FILES=$(shell find src/ -name "*.cpp")
CXX_FLAGS=-std=c++2b -O3 -Wall -Wextra -Werror -Wpedantic -Wdeprecated

$(EXE): $(OBJ_FILES)
	$(CXX) $^ $(INCLUDE) $(CXX_FLAGS) -o $@ $(LD_FLAGS)
	strip $(EXE)

src/%.o: src/%.cpp
	$(CXX) -MMD -c $< $(INCLUDE) $(CXX_FLAGS) -o $@ $(LD_FLAGS)

.PHONY: clean
clean:
	rm -f test
	find src/ -name "*.o" | xargs rm -f
	find src/ -name "*.d" | xargs rm -f

