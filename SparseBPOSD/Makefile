WARN =-Wall -Wsign-compare -Wextra
CPP=g++ -O -g ${WARN} 
hpcc1= `itpp-config --cflags`
hpcc2= `itpp-config --libs`

hpcc:
	ssh -X xliu261@cluster.hpcc.ucr.edu

111test:111test.cpp
		${CPP}  `pkg-config --cflags itpp` -o test BP.h 111test.cpp  `pkg-config --libs itpp`

gen_dep:BP.h BP.cpp modified_BP.h modified_BP.cpp gen_depol_mat.cpp
		${CPP}  `pkg-config --cflags itpp` -o gen_depol BP.h BP.cpp modified_BP.h modified_BP.cpp gen_depol_mat.cpp `pkg-config --libs itpp`
		
MBP:BP.h BP.cpp modified_BP.h modified_BP.cpp modified_BP_main.cpp
		${CPP}  `pkg-config --cflags itpp` -o MBP BP.h BP.cpp modified_BP.h modified_BP.cpp modified_BP_main.cpp `pkg-config --libs itpp`
		
OSD:BP.h BP.cpp modified_BP.h modified_BP.cpp OSD_main.cpp
		${CPP}  `pkg-config --cflags itpp` -o OSD BP.h BP.cpp modified_BP.h modified_BP.cpp OSD_main.cpp `pkg-config --libs itpp`

newBP:BP.h BP.cpp modified_BP.h modified_BP.cpp new_BP_main.cpp
		${CPP}  `pkg-config --cflags itpp` -o newBP BP.h BP.cpp modified_BP.h modified_BP.cpp new_BP_main.cpp `pkg-config --libs itpp`		

ABP:BP.h BP.cpp modified_BP.h modified_BP.cpp alternative_decoder_main.cpp
		${CPP}  `pkg-config --cflags itpp` -o ABP BP.h BP.cpp modified_BP.h modified_BP.cpp alternative_decoder_main.cpp `pkg-config --libs itpp`

CBP:BP.h BP.cpp modified_BP.h modified_BP.cpp circuit_main.cpp
		${CPP}  `pkg-config --cflags itpp` -o CBP BP.h BP.cpp modified_BP.h modified_BP.cpp circuit_main.cpp `pkg-config --libs itpp`
		
FOSD_CBP:BP.h BP.cpp modified_BP.h modified_BP.cpp circuit_main.cpp Sparse.cpp Sparse.h
		${CPP}  `pkg-config --cflags itpp` -o FOSD_CBP BP.h BP.cpp modified_BP.h modified_BP.cpp circuit_main.cpp Sparse.cpp Sparse.h `pkg-config --libs itpp`
		
FOSD_CBP_hpcc:BP.h BP.cpp modified_BP.h modified_BP.cpp circuit_main.cpp Sparse.cpp Sparse.h
		${CPP}  ${hpcc1}  -o FOSD_CBP_hpcc BP.h BP.cpp modified_BP.h modified_BP.cpp circuit_main.cpp Sparse.cpp Sparse.h ${hpcc2} 
		
speedtest:BP.h BP.cpp modified_BP.h modified_BP.cpp speedtest.cpp Sparse.cpp Sparse.h
		${CPP}  `pkg-config --cflags itpp` -o speedtest BP.h BP.cpp modified_BP.h modified_BP.cpp speedtest.cpp Sparse.cpp Sparse.h `pkg-config --libs itpp`
		
ABP-hpcc:BP.h BP.cpp modified_BP.h modified_BP.cpp alternative_decoder_main.cpp
		${CPP} ${hpcc1} -o ABP BP.h BP.cpp modified_BP.h modified_BP.cpp alternative_decoder_main.cpp ${hpcc2}

FAST:BP.h BP.cpp modified_BP.h modified_BP.cpp Fast_OSD_main.cpp
		${CPP}  `pkg-config --cflags itpp` -o Fast BP.h BP.cpp modified_BP.h modified_BP.cpp Fast_OSD_main.cpp `pkg-config --libs itpp`
		
NBP:BP.h BP.cpp modified_BP.h modified_BP.cpp new_decoders.cpp  new_decoders.h  new_decoders_main.cpp
		${CPP}  `pkg-config --cflags itpp` -o NBP BP.h BP.cpp modified_BP.h modified_BP.cpp new_decoders.cpp  new_decoders.h  new_decoders_main.cpp `pkg-config --libs itpp`

M_test:BP.h BP.cpp modified_BP.h modified_BP.cpp modified_test.cpp
		${CPP}  `pkg-config --cflags itpp` -o M_test BP.h BP.cpp modified_BP.h modified_BP.cpp modified_test.cpp `pkg-config --libs itpp`

PQ:BP.h BP.cpp modified_BP.h modified_BP.cpp gen_P_Q.cpp
		${CPP}  `pkg-config --cflags itpp` -o gen_PQ BP.h BP.cpp modified_BP.h modified_BP.cpp gen_P_Q.cpp `pkg-config --libs itpp`

PQ_noR:BP.h BP.cpp modified_BP.h modified_BP.cpp gen_P_Q_no_R.cpp
		${CPP}  `pkg-config --cflags itpp` -o gen_PQ_no_R BP.h BP.cpp modified_BP.h modified_BP.cpp gen_P_Q_no_R.cpp `pkg-config --libs itpp`


MLSR:BP.h BP.cpp modified_BP.h modified_BP.cpp ML_suc_rate_main.cpp
		${CPP}  `pkg-config --cflags itpp` -o MLSR BP.h BP.cpp modified_BP.h modified_BP.cpp ML_suc_rate_main.cpp `pkg-config --libs itpp`


MLSRwrong:BP.h BP.cpp modified_BP.h modified_BP.cpp ML_suc_rate_main.cpp
		${CPP}  `pkg-config --cflags itpp` -o MLSRwrong BP.h BP.cpp modified_BP.h modified_BP.cpp ML_suc_rate_main.cpp `pkg-config --libs itpp`

check: checkerboard.cpp  BP_quantum_func1.cpp BP_basic_func1.cpp
	${CPP}  `pkg-config --cflags itpp` -o checkerboard checkerboard.cpp BP_quantum_func1.cpp BP_basic_func1.cpp  `pkg-config --libs itpp`

test: test.cpp  BP.h BP.cpp
	${CPP}  `pkg-config --cflags itpp` -o test test.cpp BP.h BP.cpp  `pkg-config --libs itpp`
	
GPTtest: gpttest.cpp  Sparse.h Sparse.cpp
	${CPP}  `pkg-config --cflags itpp` -o Gtest gpttest.cpp Sparse.h Sparse.cpp `pkg-config --libs itpp`

test2: test.cpp 
	${CPP}  `pkg-config --cflags itpp` -o test test.cpp   `pkg-config --libs itpp`

qtest1: quantum_test1.cpp  BP_quantum_func1.cpp BP_basic_func1.cpp
	${CPP}  `pkg-config --cflags itpp` -o qtest1 quantum_test1.cpp BP_quantum_func1.cpp BP_basic_func1.cpp  `pkg-config --libs itpp`

cla_BP:BP.h BP.cpp classical_main.cpp
	${CPP}  `pkg-config --cflags itpp` -o classical_BP BP.h BP.cpp classical_main.cpp  `pkg-config --libs itpp`


quan_BP:BP.h BP.cpp quantum_main.cpp
		${CPP}  `pkg-config --cflags itpp` -o quantum_BP BP.h BP.cpp quantum_main.cpp  `pkg-config --libs itpp`
		

		

