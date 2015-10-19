CC=clang++

all: bin/encode bin/usr

bin/encode: obj/encode.o
	${CC} -o $@ $^ -L${RDKIT_ROOT}/lib -lGraphMol -lFileParsers -lSmilesParse -lSubstructMatch -lRDGeneral

bin/usr: obj/main.o
	${CC} -o $@ $^ -pthread -L${BOOST_ROOT}/lib -lboost_system -lboost_thread -lboost_filesystem -lboost_iostreams -lboost_date_time -L${RDKIT_ROOT}/lib -lGraphMol -lFileParsers -lSmilesParse -lSubstructMatch -lRDGeneral -L${MONGODBCXXDRIVER_ROOT}/sharedclient -lmongoclient -L${POCO_ROOT}/lib -lPocoFoundation -lPocoNet

obj/encode.o: src/encode.cpp
	${CC} -o $@ $< -c -std=c++14 -O2 -Wall -I${BOOST_ROOT} -I${RDKIT_ROOT}/include/rdkit

obj/main.o: src/main.cpp
	${CC} -o $@ $< -c -std=c++14 -O2 -mavx -Wall -Wno-unused-local-typedef -Wno-deprecated-declarations -Wno-deprecated-register -I${BOOST_ROOT} -I${RDKIT_ROOT}/include/rdkit -I${MONGODBCXXDRIVER_ROOT}/src -I${POCO_ROOT}/include

clean:
	rm -f bin/* obj/*

