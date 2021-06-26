
CC := g++
CFLAGS := -Wno-sign-compare -g -std=c++17
INCLUDES := -I . -I eigen -I eigen/unsupported
TARGET := data_gen

# $(wildcard *.cpp /xxx/xxx/*.cpp): get all .cpp files from the current directory and dir "/xxx/xxx/"
SRCS := $(wildcard *.cxx)
# $(patsubst %.cpp,%.o,$(SRCS)): substitute all ".cpp" file name strings to ".o" file name strings
OBJS := $(patsubst %.cxx,%.o,$(SRCS))

all: $(TARGET)
$(TARGET): $(OBJS)
	$(CC) -o $@ $^
%.o: %.cxx
	$(CC) $(CFLAGS) $(INCLUDES) -c $<
clean:
	rm -rf $(TARGET) *.o

.PHONY: all clean