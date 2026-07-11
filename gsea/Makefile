# define the C++ compiler to use.
# Default to the portable `c++` driver (resolves to g++ on Linux and clang++
# on macOS). Override on the command line if needed, e.g. `make CC=g++`.
CC = c++

# define any compile-time flags
CFLAGS = -Wall -g -O3

# define any directories containing header files other than /usr/include
# (headers, including the bundled TCLAP argument parser, live under ./src)
INCLUDES = -I./src

# define the C++ source files
MAIN1 = ./src/calculategsea.cpp

# shared support code (random number generators and stats helpers)
APIS = ./src/math_api.cpp ./src/rvgs.cpp ./src/rngs.cpp
API_OBJS = $(APIS:.cpp=.o)
MAIN1_OBJS = $(MAIN1:.cpp=.o)

# define the executable file
MAIN1_APP = bin/mageckGSEA

all:    $(MAIN1_APP)

$(MAIN1_APP): $(API_OBJS) $(MAIN1_OBJS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $(MAIN1_APP) $(API_OBJS) $(MAIN1_OBJS)

# suffix replacement rule for building .o's from .cpp's
.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $<  -o $@

clean:
	$(RM) $(API_OBJS) $(MAIN1_OBJS) $(MAIN1_APP)
