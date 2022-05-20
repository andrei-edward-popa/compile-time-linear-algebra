CXX=g++
APP=app
TEST=test
LD_FLAGS=-lfmt
INCLUDE=-Iinclude
APP_OBJ_FILES=$(APP_CPP_FILES:.cpp=.o)
TEST_OBJ_FILES=$(TEST_CPP_FILES:.cpp=.o)
APP_CPP_FILES=$(shell find src/ -name "*.cpp")
TEST_CPP_FILES=$(shell find tests/ -name "*.cpp")
D_FILES=$(CPP_FILES:.cpp=.d)
CXX_FLAGS=-DPRECISION=${PRECISION} -std=c++2b -Ofast -Wall -Wextra -Werror -Wpedantic -Wdeprecated -Wconversion -Wshadow -Wnon-virtual-dtor -Wcast-align -Wpointer-arith -Wunused -Woverloaded-virtual -flto=auto

$(APP): $(APP_OBJ_FILES)
	$(CXX) $^ $(INCLUDE) $(CXX_FLAGS) -o $@ $(LD_FLAGS)
	strip $(APP)
	
$(TEST): $(TEST_OBJ_FILES)
	$(CXX) $^ $(INCLUDE) $(CXX_FLAGS) -o $@ $(LD_FLAGS) -lgtest
	strip $(TEST)

src/%.o: src/%.cpp
	$(CXX) -MMD -c $< $(INCLUDE) $(CXX_FLAGS) -o $@ $(LD_FLAGS)
	
tests/%.o: tests/%.cpp
	$(CXX) -MMD -c $< $(INCLUDE) $(CXX_FLAGS) -o $@ $(LD_FLAGS)

.PHONY: mrproper
mrproper:
	rm -f app
	rm -f test
	find src/ -name "*.o" | xargs rm -f
	find src/ -name "*.d" | xargs rm -f
	find tests/ -name "*.o" | xargs rm -f
	find tests/ -name "*.d" | xargs rm -f
	
-include $(D_FILES)

ifndef PRECISION
PRECISION=4
endif

