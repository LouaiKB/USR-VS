CC=clang++

all: bin/embed_static bin/embed bin/validatesdf bin/encode bin/usr

bin/embed_static: obj/embed.o
	${CC} -o $@ $^ -static -pthread -L${RDKIT_ROOT}/lib -lRDKitDepictor_static -lRDKitDistGeomHelpers_static -lRDKitDistGeometry_static -lRDKitFileParsers_static -lRDKitForceFieldHelpers_static -lRDKitSmilesParse_static -lRDKitSubstructMatch_static -lRDKitGraphMol_static -lRDKitForceField_static -lRDKitEigenSolvers_static -lRDKitAlignment_static -lRDKitRDGeometryLib_static -lRDKitRDGeneral_static -L${BOOST_ROOT}/lib -lboost_thread -lboost_system

bin/embed: obj/embed.o
	${CC} -o $@ $^ -L${RDKIT_ROOT}/lib -lRDKitDistGeomHelpers -lRDKitFileParsers -lRDKitSmilesParse -lRDKitDepictor -lRDKitGraphMol -lRDKitRDGeometryLib -lRDKitRDGeneral

bin/validatesdf: obj/validatesdf.o
	${CC} -o $@ $^ -L${RDKIT_ROOT}/lib -lRDKitFileParsers -lRDKitSmilesParse -lRDKitSubstructMatch -lRDKitGraphMol

bin/encode: obj/encode.o
	${CC} -o $@ $^ -L${RDKIT_ROOT}/lib -lRDKitFileParsers -lRDKitSmilesParse -lRDKitSubstructMatch -lRDKitGraphMol -lRDKitRDGeneral

bin/usr: obj/main.o obj/io_service_pool.o obj/safe_counter.o
	${CC} -o $@ $^ -pthread -L${BOOST_ROOT}/lib -lboost_thread -lboost_system -lboost_filesystem -lboost_date_time -L${RDKIT_ROOT}/lib -lRDKitMolTransforms -lRDKitFingerprints -lRDKitFileParsers -lRDKitSmilesParse -lRDKitSubstructMatch -lRDKitDepictor -lRDKitGraphMol -lRDKitAlignment -lRDKitRDGeometryLib -lRDKitRDGeneral -L${MONGODBCXXDRIVER_ROOT}/sharedclient -lmongoclient

obj/embed.o: src/embed.cpp
	${CC} -o $@ $< -c -std=c++14 -O2 -Wall -I${RDKIT_ROOT}/include/rdkit

obj/validatesdf.o: src/validatesdf.cpp
	${CC} -o $@ $< -c -std=c++14 -O2 -Wall -I${RDKIT_ROOT}/include/rdkit

obj/encode.o: src/encode.cpp
	${CC} -o $@ $< -c -std=c++14 -O2 -Wall -I${RDKIT_ROOT}/include/rdkit

obj/main.o: src/main.cpp
	${CC} -o $@ $< -c -std=c++14 -O2 -Wall -Wno-unused-local-typedef -Wno-deprecated-declarations -Wno-deprecated-register -I${BOOST_ROOT}/include -I${RDKIT_ROOT}/include/rdkit -I${MONGODBCXXDRIVER_ROOT}/src

obj/%.o: src/%.cpp
	${CC} -o $@ $< -c -std=c++14 -O2 -Wall -I${BOOST_ROOT}/include

clean:
	rm -f bin/* obj/*
