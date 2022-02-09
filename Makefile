CC=clang++

#  The use of static rdkit libraries for compiling the embed program was due to the consideration of portability: embed was developed on one machine and deployed on another machine
# all: bin/embed_static bin/embed bin/validatesdf bin/encode bin/usr
# bin/embed_static: obj/embed.o
# 	${CC} -o $@ $^ -static -pthread -L${RDKIT_ROOT}/lib -lRDKitDepictor_static -lRDKitDistGeomHelpers_static -lRDKitDistGeometry_static -lRDKitFileParsers_static -lRDKitForceFieldHelpers_static -lRDKitSmilesParse_static -lRDKitSubstructMatch_static -lRDKitGraphMol_static -lRDKitForceField_static -lRDKitEigenSolvers_static -lRDKitAlignment_static -lRDKitRDGeometryLib_static -lRDKitRDGeneral_static -L${BOOST_ROOT}/lib -lboost_thread -lboost_system

all: bin/embed bin/validatesdf bin/encode bin/usr

bin/embed: obj/embed.o
	${CC} -o $@ $^ -L${RDKIT_ROOT}/build/lib -lChemicalFeatures -lInfoTheory -lSimDivPickers -lhc -lSLNParse -lTrajectory -lReducedGraphs -lStructChecker -lMMPA -lMolHash -lFMCS -lMolDraw2D -lMolCatalog -lShapeHelpers -lMolChemicalFeatures -lFragCatalog -lChemReactions -lFingerprints -lDescriptors -lPartialCharges -lFilterCatalog -lSubgraphs -lChemTransforms -lDepictor -lCatalogs -lDistGeomHelpers -lDistGeometry -lMolAlign -lMolTransforms -lFileParsers -lForceFieldHelpers -lSmilesParse -lSubstructMatch -lGraphMol -lForceField -lOptimizer -lEigenSolvers -lAlignment -lRDGeometryLib -lDataStructs -lRDGeneral

bin/validatesdf: obj/validatesdf.o
	${CC} -o $@ $^ -L${RDKIT_ROOT}/build/lib -lChemicalFeatures -lInfoTheory -lSimDivPickers -lhc -lSLNParse -lTrajectory -lReducedGraphs -lStructChecker -lMMPA -lMolHash -lFMCS -lMolDraw2D -lMolCatalog -lShapeHelpers -lMolChemicalFeatures -lFragCatalog -lChemReactions -lFingerprints -lDescriptors -lPartialCharges -lFilterCatalog -lSubgraphs -lChemTransforms -lDepictor -lCatalogs -lDistGeomHelpers -lDistGeometry -lMolAlign -lMolTransforms -lFileParsers -lForceFieldHelpers -lSmilesParse -lSubstructMatch -lGraphMol -lForceField -lOptimizer -lEigenSolvers -lAlignment -lRDGeometryLib -lDataStructs -lRDGeneral

bin/encode: obj/encode.o
	${CC} -o $@ $^ -L${RDKIT_ROOT}/build/lib -lChemicalFeatures -lInfoTheory -lSimDivPickers -lhc -lSLNParse -lTrajectory -lReducedGraphs -lStructChecker -lMMPA -lMolHash -lFMCS -lMolDraw2D -lMolCatalog -lShapeHelpers -lMolChemicalFeatures -lFragCatalog -lChemReactions -lFingerprints -lDescriptors -lPartialCharges -lFilterCatalog -lSubgraphs -lChemTransforms -lDepictor -lCatalogs -lDistGeomHelpers -lDistGeometry -lMolAlign -lMolTransforms -lFileParsers -lForceFieldHelpers -lSmilesParse -lSubstructMatch -lGraphMol -lForceField -lOptimizer -lEigenSolvers -lAlignment -lRDGeometryLib -lDataStructs -lRDGeneral

bin/usr: obj/main.o obj/io_service_pool.o obj/safe_counter.o
	${CC} -o $@ $^ -pthread -L${BOOST_LIBS} -lboost_date_time -lboost_thread -lboost_system -lboost_filesystem -L${RDKIT_ROOT}/build/lib -lChemicalFeatures -lInfoTheory -lSimDivPickers -lhc -lSLNParse -lTrajectory -lReducedGraphs -lStructChecker -lMMPA -lMolHash -lFMCS -lMolDraw2D -lMolCatalog -lShapeHelpers -lMolChemicalFeatures -lFragCatalog -lChemReactions -lFingerprints -lDescriptors -lPartialCharges -lFilterCatalog -lSubgraphs -lChemTransforms -lDepictor -lCatalogs -lDistGeomHelpers -lDistGeometry -lMolAlign -lMolTransforms -lFileParsers -lForceFieldHelpers -lSmilesParse -lSubstructMatch -lGraphMol -lForceField -lOptimizer -lEigenSolvers -lAlignment -lRDGeometryLib -lDataStructs -lRDGeneral -L${MONGODBCXXDRIVER_ROOT}/build/src/mongocxx -lmongocxx -L${MONGODBCXXDRIVER_ROOT}/build/src/bsoncxx -lbsoncxx

obj/embed.o: src/embed.cpp
	${CC} -o $@ $< -c -std=c++14 -O2 -Wall -I${RDKIT_ROOT}/Code

obj/validatesdf.o: src/validatesdf.cpp
	${CC} -o $@ $< -c -std=c++14 -O2 -Wall -I${RDKIT_ROOT}/Code

obj/encode.o: src/encode.cpp
	${CC} -o $@ $< -c -std=c++14 -O2 -Wall -I${RDKIT_ROOT}/Code

obj/main.o: src/main.cpp
	${CC} -o $@ $< -c -std=c++14 -O2 -Wall -Wno-unused-local-typedef -Wno-deprecated-declarations -Wno-deprecated-register -I${BOOST_HEADERS} -I${RDKIT_ROOT}/Code -I${MONGODBCXX_HEADERS} -I${BSONCXX_HEADERS}

obj/%.o: src/%.cpp
	${CC} -o $@ $< -c -std=c++14 -O2 -Wall -I${BOOST_HEADERS}

clean:
	rm -f bin/* obj/*
