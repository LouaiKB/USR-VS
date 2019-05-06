CC=clang++

all: bin/embed_static bin/embed bin/validatesdf bin/encode bin/usrd

bin/embed_static: obj/embed.o
	${CC} -o $@ $^ -static -pthread -L${RDKIT_ROOT}/lib -lRDKitMolDraw2D_static -lRDKitDistGeomHelpers_static -lRDKitChemReactions_static -lRDKitFileParsers_static -lRDKitDistGeometry_static -lRDKitDepictor_static -lRDKitForceFieldHelpers_static -lRDKitMolTransforms_static -lRDKitSubstructMatch_static -lRDKitSmilesParse_static -lRDKitForceField_static -lRDKitEigenSolvers_static -lRDKitAlignment_static -lRDKitcoordgenlib_static -lRDKitmaeparser_static -lRDKitGraphMol_static -lRDKitRDGeometryLib_static -lRDKitRDGeneral_static

bin/embed: obj/embed.o
	${CC} -o $@ $^ -L${RDKIT_ROOT}/lib -lRDKitMolDraw2D -lRDKitDistGeomHelpers -lRDKitFileParsers -lRDKitDepictor -lRDKitGraphMol -lRDKitRDGeneral

bin/validatesdf: obj/validatesdf.o
	${CC} -o $@ $^ -L${RDKIT_ROOT}/lib -lRDKitFileParsers -lRDKitSubstructMatch -lRDKitSmilesParse -lRDKitGraphMol

bin/encode: obj/encode.o
	${CC} -o $@ $^ -L${RDKIT_ROOT}/lib -lRDKitFileParsers -lRDKitSubstructMatch -lRDKitSmilesParse -lRDKitGraphMol -lRDKitRDGeneral

bin/usrd: obj/main.o obj/io_service_pool.o obj/safe_counter.o
	${CC} -o $@ $^ -pthread -L${BOOST_ROOT}/lib -lboost_date_time -L${RDKIT_ROOT}/lib -lRDKitMolDraw2D -lRDKitFingerprints -lRDKitFileParsers -lRDKitDepictor -lRDKitMolTransforms -lRDKitSubstructMatch -lRDKitSmilesParse -lRDKitAlignment -lRDKitGraphMol -lRDKitRDGeometryLib -lRDKitRDGeneral -L${MONGO_CXX_DRIVER_ROOT}/lib64 -lmongocxx -lbsoncxx

obj/embed.o: src/embed.cpp
	${CC} -o $@ $< -c -std=c++17 -O2 -Wall -I${BOOST_ROOT}/include -I${RDKIT_ROOT}/include/rdkit

obj/validatesdf.o: src/validatesdf.cpp
	${CC} -o $@ $< -c -std=c++17 -O2 -Wall -I${BOOST_ROOT}/include -I${RDKIT_ROOT}/include/rdkit

obj/encode.o: src/encode.cpp
	${CC} -o $@ $< -c -std=c++17 -O2 -Wall -I${BOOST_ROOT}/include -I${RDKIT_ROOT}/include/rdkit

obj/main.o: src/main.cpp
	${CC} -o $@ $< -c -std=c++17 -O2 -Wall -Wno-unused-local-typedef -Wno-deprecated-declarations -Wno-deprecated-register -I${BOOST_ROOT}/include -I${RDKIT_ROOT}/include/rdkit -I${MONGO_CXX_DRIVER_ROOT}/include/mongocxx/v_noabi -I${MONGO_CXX_DRIVER_ROOT}/include/bsoncxx/v_noabi

obj/%.o: src/%.cpp
	${CC} -o $@ $< -c -std=c++17 -O2 -Wall -I${BOOST_ROOT}/include

clean:
	rm -f bin/* obj/*
