CC=clang++

all: bin/embed bin/encode bin/usr

bin/embed: obj/embed.o
	${CC} -o $@ $^ -L${RDKIT_ROOT}/lib -lGraphMol -lFileParsers -lSmilesParse -lDistGeomHelpers

bin/encode: obj/encode.o
	${CC} -o $@ $^ -L${RDKIT_ROOT}/lib -lGraphMol -lFileParsers -lSmilesParse -lSubstructMatch -lRDGeneral

bin/usr: obj/main.o obj/io_service_pool.o obj/safe_counter.o
	${CC} -o $@ $^ -pthread -L${BOOST_ROOT}/lib -lboost_system -lboost_filesystem -lboost_date_time -L${RDKIT_ROOT}/lib -lRDGeneral -lGraphMol -lFileParsers -lSmilesParse -lFingerprints -lSubstructMatch -lAlignment -lMolTransforms -L${MONGODBCXXDRIVER_ROOT}/sharedclient -lmongoclient

obj/embed.o: src/embed.cpp
	${CC} -o $@ $< -c -std=c++14 -O2 -Wall -I${RDKIT_ROOT}/include/rdkit

obj/encode.o: src/encode.cpp
	${CC} -o $@ $< -c -std=c++14 -O2 -Wall -I${RDKIT_ROOT}/include/rdkit

obj/main.o: src/main.cpp
	${CC} -o $@ $< -c -std=c++14 -O2 -Wall -Wno-unused-local-typedef -Wno-deprecated-declarations -Wno-deprecated-register -I${BOOST_ROOT} -I${RDKIT_ROOT}/include/rdkit -I${MONGODBCXXDRIVER_ROOT}/src

obj/%.o: src/%.cpp
	${CC} -o $@ $< -c -std=c++14 -O2 -Wall -I${BOOST_ROOT}

clean:
	rm -f bin/* obj/*
