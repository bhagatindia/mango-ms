CXX = g++
MSTOOLKIT = mstoolkit
HASH = hash
PROTOBUF = protobuf
HARDKLOR = hardklor
override CXXFLAGS +=  -g  -std=c++11 -Wall -Wextra -static -Wno-char-subscripts -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -D__LINUX__ -I$(MSTOOLKIT)/include
EXECNAME = mango.exe
OBJS = mango.o mango_Preprocess.o mango_Search.o mango_MassSpecUtils.o  $(HASH)/mango-hash.o $(HASH)/protein_pep_hash.pb.o
DEPS = mango.h Common.h mango_Data.h mango_DataInternal.h mango_Preprocess.h mango_MassSpecUtils.h

LIBS = -L$(MSTOOLKIT) -lmstoolkitlite -lm -pthread -L/usr/local/lib -lprotobuf 
ifdef MSYSTEM
   LIBS += -lws2_32
endif



mango.exe: $(OBJS)
	git submodule init; git submodule update
	cd $(MSTOOLKIT) ; make lite 
	cd $(HASH) ; make; 
	rm MSToolkit; ln -s mstoolkit MSToolkit
	cd $(HARDKLOR) ; make; 
	${CXX} $(CXXFLAGS) $(OBJS) $(LIBS) -o ${EXECNAME}

mango.o: mango.cpp $(DEPS)
	cd $(HASH) ; make; 
	git submodule init; git submodule update
	${CXX} ${CXXFLAGS} mango.cpp -c

mango_Preprocess.o: mango_Preprocess.cpp Common.h mango_Preprocess.h mango.h Common.h mango_Data.h mango_DataInternal.h
	git submodule init; git submodule update
	${CXX} ${CXXFLAGS} mango_Preprocess.cpp -c

mango_Search.o: mango_Search.cpp Common.h mango_Search.h mango.h Common.h mango_Data.h mango_DataInternal.h
	git submodule init; git submodule update
	${CXX} ${CXXFLAGS} mango_Search.cpp -c

mango_MassSpecUtils.o: mango_MassSpecUtils.cpp Common.h mango_MassSpecUtils.h mango.h Common.h mango_Data.h mango_DataInternal.h
	git submodule init; git submodule update
	${CXX} ${CXXFLAGS} mango_MassSpecUtils.cpp -c

clean:
	rm -f *.o ${EXECNAME}
	cd $(MSTOOLKIT) ; make clean
	cd $(HASH) ; make clean
	cd $(PROTOBUF); make clean
