CXX = g++
CXX_FLAGS +=-g -Wall -O3 -std=c++17
PROTO_FLAGS = `pkg-config --cflags --libs protobuf`

BIN = metric

BUILD_DIR = ../build

CPP = main.cpp $(wildcard *.cpp)

CC = $(wildcard include/*.pb.cc)

OBJ = $(CPP:%.cpp=$(BUILD_DIR)/%.o)

DEP = $(OBJ:%.o=%.d)

$(BIN) : $(BUILD_DIR)/$(BIN)

$(BUILD_DIR)/$(BIN) : $(OBJ)
	mkdir -p $(@D)
	$(CXX) $(CXX_FLAGS) $(CC) $^ -o $@ $(PROTO_FLAGS)

-include $(DEP)

$(BUILD_DIR)/%.o : %.cpp
	mkdir -p $(@D)

	$(CXX) $(CXX_FLAGS) -MMD -c $< -o $@
.PHONY : clean
clean :
	-rm $(BUILD_DIR)/$(BIN) $(OBJ) $(DEP)
