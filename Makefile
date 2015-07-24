CC=clang++

bin/usr: obj/main.o
	${CC} -o $@ $^ -pthread -L${OPENBABEL_ROOT}/lib -lopenbabel -L${BOOST_ROOT}/lib -lboost_system -lboost_thread -lboost_filesystem -lboost_iostreams -lboost_date_time -L${POCO_ROOT}/lib/Linux/x86_64 -lPocoFoundation -lPocoNet -L${MONGODBCXXDRIVER_ROOT}/sharedclient -lmongoclient

obj/main.o: src/main.cpp
	${CC} -o $@ $< -c -std=c++11 -O2 -mavx -Wall -Wno-unused-local-typedef -Wno-deprecated-declarations -Wno-deprecated-register -I${OPENBABEL_ROOT}/include -I${BOOST_ROOT} -I${MONGODBCXXDRIVER_ROOT}/src -I${POCO_ROOT}/Foundation/include -I${POCO_ROOT}/Net/include

clean:
	rm -f bin/usr obj/*
