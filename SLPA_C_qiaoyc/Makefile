# modified by Yuchen Qiao
#

# library path
GPERF_LIB_LOCATION=~/user/lib

# compiler 
CC=g++ # general version
#CC=/usr/local/bin/g++-6

# compiling options flags
CFLAGS=-c -O3 -fopenmp
#CFLAGS=-c -Wall
#CFLAGS= -rdynamic -g -pg -c 
#CFLAGS= -static -g -pg -c 
#CFLAGS= -g -pg -c 

# link options flags
LDFLAGS=

# link for tcmalloc
TC_LDFLAGS=
TC_LDFLAGS+= -ltcmalloc -L $(GPERF_LIB_LOCATION) -Wl,-R$(GPERF_LIB_LOCATION)

all: main_SLPA main_SLPA_tc

#multi-thread
main_SLPA: main_SLPA.o classSLPA.o Net.o NODE.o CommonFuns.o fileOpts.o rndNumbers.o
	$(CC) main_SLPA.o classSLPA.o Net.o NODE.o CommonFuns.o fileOpts.o rndNumbers.o -o SLPA -pthread -fopenmp 

main_SLPA_tc: main_SLPA.o classSLPA.o Net.o NODE.o CommonFuns.o fileOpts.o rndNumbers.o
	$(CC) main_SLPA.o classSLPA.o Net.o NODE.o CommonFuns.o fileOpts.o rndNumbers.o -o SLPA_tc -pthread -fopenmp $(TC_LDFLAGS)

#Normal
#main_SLPA: main_SLPA.o classSLPA.o Net.o NODE.o CommonFuns.o fileOpts.o rndNumbers.o
#	$(CC) main_SLPA.o classSLPA.o Net.o NODE.o CommonFuns.o fileOpts.o rndNumbers.o -o SLPA -lpthread

#debug
#main_SLPA: main_SLPA.o classSLPA.o Net.o NODE.o CommonFuns.o fileOpts.o rndNumbers.o
#	$(CC) -g -pg main_SLPA.o classSLPA.o Net.o NODE.o CommonFuns.o fileOpts.o rndNumbers.o -o SLPA


main_SLPA.o: main_SLPA.cpp
	$(CC) $(CFLAGS) main_SLPA.cpp

classSLPA.o: SLPA.cpp
	$(CC) $(CFLAGS) SLPA.cpp -o classSLPA.o

Net.o: Net.cpp
	$(CC) $(CFLAGS) Net.cpp

NODE.o: NODE.cpp
	$(CC) $(CFLAGS) NODE.cpp

CommonFuns.o: CommonFuns.cpp
	$(CC) $(CFLAGS) CommonFuns.cpp

fileOpts.o: fileOpts.cpp
	$(CC) $(CFLAGS) fileOpts.cpp

rndNumbers.o: rndNumbers.cpp
	$(CC) $(CFLAGS) rndNumbers.cpp


clean:
	rm -rf *o SLPA SLPA_tc
