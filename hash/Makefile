CXX = g++
override CXXFLAGS += -g  -std=c++11 -Wall -Wextra -static 

LIBS=lprotobuf

mango-hash: protein_pep_hash.proto mango-hash.cpp
	protoc --cpp_out=./ protein_pep_hash.proto	
	g++ -g -c -std=c++11 mango-hash.cpp protein_pep_hash.pb.cc
	ar rvs mango-hash.a mango-hash.o protein_pep_hash.pb.o

phd-main:
	g++   -g -std=c++11 phd_main.cpp mango-hash.a -L/usr/local/lib -lprotobuf -o phd-main

clean:
	rm protein_pep_hash.pb.cc protein_pep_hash.pb.h mango-hash.o mango-hash.a protein_pep_hash.pb.o phd-main

