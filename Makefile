CC=clang++

# all: bin/embed_static bin/embed bin/validatesdf bin/encode bin/usr
all: bin/embed bin/validatesdf bin/encode bin/usr

# bin/embed_static: obj/embed.o
# 	${CC} -o $@ $^ -static -pthread -L${RDKIT_ROOT}/lib -lRDKitDepictor_static -lRDKitDistGeomHelpers_static -lRDKitDistGeometry_static -lRDKitFileParsers_static -lRDKitForceFieldHelpers_static -lRDKitSmilesParse_static -lRDKitSubstructMatch_static -lRDKitGraphMol_static -lRDKitForceField_static -lRDKitEigenSolvers_static -lRDKitAlignment_static -lRDKitRDGeometryLib_static -lRDKitRDGeneral_static -L${BOOST_ROOT}/lib -lboost_thread -lboost_system

# bin/embed: obj/embed.o
# 	${CC} -o $@ $^ -L${RDKIT_ROOT}/lib -lRDKitavalon_clib -lRDKitcoordgen -lRDKitfreesasa_clib -lRDKitga -lRDKithc -lRDKitmaeparser -lRDKitMolInterchange -lRDKitMolStandardize -lRDKitMolTransforms -lRDKitDistGeomHelpers -lRDKitO3AAlign -lRDKitOptimizer -lRDKitPartialCharges -lRDKitRDBoost -lRDKitRDGeneral -lRDKitRDGeometryLib -lRDKitRDInchiLib -lRDKitRDStreams -lRDKitReducedGraphs -lRDKitRGroupDecomposition -lRDKitRingDecomposerLib -lRDKitScaffoldNetwork -lRDKitShapeHelpers -lRDKitSimDivPickers -lRDKitSLNParse -lRDKitSmilesParse -lRDKitSubgraphs -lRDKitSubstructLibrary -lRDKitSubstructMatch -lRDKitTautomerQuery -lRDKitTrajectory -lRDKitAbbreviations -lRDKitAlignment -lRDKitAvalonLib -lRDKitCatalogs -lRDKitChemicalFeatures -lRDKitChemReactions -lRDKitChemTransforms -lRDKitCIPLabeler -lRDKitDataStructs -lRDKitDepictor -lRDKitDeprotect -lRDKitDescriptors -lRDKitDistGeometry -lRDKitDistGeomHelpers -lRDKitEHTLib -lRDKitEigenSolvers -lRDKitFileParsers -lRDKitFilterCatalog -lRDKitFingerprints -lRDKitFMCS -lRDKitForceFieldHelpers -lRDKitForceField -lRDKitFragCatalog -lRDKitFreeSASALib -lRDKitGraphMol -lRDKitInchi -lRDKitInfoTheory -lRDKitMMPA -lRDKitMolAlign -lRDKitMolCatalog -lRDKitMolChemicalFeatures -lRDKitMolDraw2D -lRDKitMolEnumerator -lRDKitMolHash 
bin/embed: obj/embed.o
	${CC} -o $@ $^ -L ~/RDKit/RDKIT_2016/build/lib -lChemicalFeatures -lInfoTheory -lSimDivPickers -lhc -lSLNParse -lTrajectory -lReducedGraphs -lStructChecker -lMMPA -lMolHash -lFMCS -lMolDraw2D -lMolCatalog -lShapeHelpers -lMolChemicalFeatures -lFragCatalog -lChemReactions -lFingerprints -lDescriptors -lPartialCharges -lFilterCatalog -lSubgraphs -lChemTransforms -lDepictor -lCatalogs -lDistGeomHelpers -lDistGeometry -lMolAlign -lMolTransforms -lFileParsers -lForceFieldHelpers -lSmilesParse -lSubstructMatch -lGraphMol -lForceField -lOptimizer -lEigenSolvers -lAlignment -lRDGeometryLib -lDataStructs -lRDGeneral

# bin/validatesdf: obj/validatesdf.o
# 	${CC} -o $@ $^ -L${RDKIT_ROOT}/lib -lRDKitFileParsers -lRDKitSmilesParse -lRDKitSubstructMatch -lRDKitGraphMol

bin/validatesdf: obj/validatesdf.o
	${CC} -o $@ $^ -L ~/RDKit/RDKIT_2016/build/lib -lChemicalFeatures -lInfoTheory -lSimDivPickers -lhc -lSLNParse -lTrajectory -lReducedGraphs -lStructChecker -lMMPA -lMolHash -lFMCS -lMolDraw2D -lMolCatalog -lShapeHelpers -lMolChemicalFeatures -lFragCatalog -lChemReactions -lFingerprints -lDescriptors -lPartialCharges -lFilterCatalog -lSubgraphs -lChemTransforms -lDepictor -lCatalogs -lDistGeomHelpers -lDistGeometry -lMolAlign -lMolTransforms -lFileParsers -lForceFieldHelpers -lSmilesParse -lSubstructMatch -lGraphMol -lForceField -lOptimizer -lEigenSolvers -lAlignment -lRDGeometryLib -lDataStructs -lRDGeneral

# bin/encode: obj/encode.o
# 	${CC} -o $@ $^ -L${RDKIT_ROOT}/lib -lRDKitFileParsers -lRDKitSmilesParse -lRDKitSubstructMatch -lRDKitGraphMol -lRDKitRDGeneral

bin/encode: obj/encode.o
	${CC} -o $@ $^ -L ~/RDKit/RDKIT_2016/build/lib -lChemicalFeatures -lInfoTheory -lSimDivPickers -lhc -lSLNParse -lTrajectory -lReducedGraphs -lStructChecker -lMMPA -lMolHash -lFMCS -lMolDraw2D -lMolCatalog -lShapeHelpers -lMolChemicalFeatures -lFragCatalog -lChemReactions -lFingerprints -lDescriptors -lPartialCharges -lFilterCatalog -lSubgraphs -lChemTransforms -lDepictor -lCatalogs -lDistGeomHelpers -lDistGeometry -lMolAlign -lMolTransforms -lFileParsers -lForceFieldHelpers -lSmilesParse -lSubstructMatch -lGraphMol -lForceField -lOptimizer -lEigenSolvers -lAlignment -lRDGeometryLib -lDataStructs -lRDGeneral

# bin/usr: obj/main.o obj/io_service_pool.o obj/safe_counter.o
# 	${CC} -o $@ $^ -pthread -L${RDKIT_ROOT}/lib -lboost_thread -lboost_system -lboost_filesystem -lboost_date_time -L${RDKIT_ROOT}/lib -lRDKitMolTransforms -lRDKitFingerprints -lRDKitFileParsers -lRDKitSmilesParse -lRDKitSubstructMatch -lRDKitDepictor -lRDKitGraphMol -lRDKitAlignment -lRDKitRDGeometryLib -lRDKitRDGeneral -L${MONGODBCXXDRIVER_ROOT}/sharedclient -lmongoclient

bin/usr: obj/main.o obj/io_service_pool.o obj/safe_counter.o
	${CC} -o $@ $^ -pthread -L${BOOST_ROOT}/stage/lib -lboost_date_time -lboost_thread -lboost_system -lboost_filesystem -L ~/RDKit/RDKIT_2016/build/lib -lChemicalFeatures -lInfoTheory -lSimDivPickers -lhc -lSLNParse -lTrajectory -lReducedGraphs -lStructChecker -lMMPA -lMolHash -lFMCS -lMolDraw2D -lMolCatalog -lShapeHelpers -lMolChemicalFeatures -lFragCatalog -lChemReactions -lFingerprints -lDescriptors -lPartialCharges -lFilterCatalog -lSubgraphs -lChemTransforms -lDepictor -lCatalogs -lDistGeomHelpers -lDistGeometry -lMolAlign -lMolTransforms -lFileParsers -lForceFieldHelpers -lSmilesParse -lSubstructMatch -lGraphMol -lForceField -lOptimizer -lEigenSolvers -lAlignment -lRDGeometryLib -lDataStructs -lRDGeneral -L ~/Documents/mongo-cxx-driver-r3.6.6/build/src/mongocxx -lmongocxx -L ~/Documents/mongo-cxx-driver-r3.6.6/build/src/bsoncxx -lbsoncxx


obj/embed.o: src/embed.cpp
	${CC} -o $@ $< -c -std=c++14 -O2 -Wall -I${RDKIT_ROOT}/include/rdkit2016

obj/validatesdf.o: src/validatesdf.cpp
	${CC} -o $@ $< -c -std=c++14 -O2 -Wall -I${RDKIT_ROOT}/include/rdkit2016

obj/encode.o: src/encode.cpp
	${CC} -o $@ $< -c -std=c++14 -O2 -Wall -I${RDKIT_ROOT}/include/rdkit2016

obj/main.o: src/main.cpp
	${CC} -o $@ $< -c -std=c++14 -O2 -Wall -Wno-unused-local-typedef -Wno-deprecated-declarations -Wno-deprecated-register -I${BOOST_ROOT} -I${RDKIT_ROOT}/include/rdkit2016 -I /usr/local/include/mongocxx/v_noabi -I /usr/local/include/bsoncxx/v_noabi

obj/%.o: src/%.cpp
	${CC} -o $@ $< -c -std=c++14 -O2 -Wall -I${BOOST_ROOT}

clean:
	rm -f bin/* obj/*


