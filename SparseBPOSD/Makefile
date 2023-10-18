WARN =-Wall -Wsign-compare -Wextra
CPP=g++ -O -g ${WARN}
hpcc1= `itpp-config --cflags`
hpcc2= `itpp-config --libs`

hpcc:
	ssh -X xliu261@cluster.hpcc.ucr.edu


		
FOSD_CBP:BP.h BP.cpp modified_BP.h modified_BP.cpp circuit_main.cpp Sparse.cpp Sparse.h
		${CPP}  `pkg-config --cflags itpp` -o FOSD_CBP BP.h BP.cpp modified_BP.h modified_BP.cpp circuit_main.cpp Sparse.cpp Sparse.h `pkg-config --libs itpp`
		
FOSD_CBPO3:BP.h BP.cpp modified_BP.h modified_BP.cpp circuit_main.cpp Sparse.cpp Sparse.h
		${CPP} -O3 `pkg-config --cflags itpp` -o FOSD_CBP_O3 BP.h BP.cpp modified_BP.h modified_BP.cpp circuit_main.cpp Sparse.cpp Sparse.h `pkg-config --libs itpp`
		
MBP1:BP.h BP.cpp modified_BP.h modified_BP.cpp MBP1_del_cols_main.cpp Sparse.cpp Sparse.h
		${CPP}  `pkg-config --cflags itpp` -o MBP1 BP.h BP.cpp modified_BP.h modified_BP.cpp MBP1_del_cols_main.cpp Sparse.cpp Sparse.h `pkg-config --libs itpp`
		
MBP2:BP.h BP.cpp modified_BP.h modified_BP.cpp MBP2_del_cols_right_main.cpp Sparse.cpp Sparse.h
		${CPP}  `pkg-config --cflags itpp` -o MBP2 BP.h BP.cpp modified_BP.h modified_BP.cpp MBP2_del_cols_right_main.cpp Sparse.cpp Sparse.h `pkg-config --libs itpp`
		
FOSD_CBP_hpcc:BP.h BP.cpp modified_BP.h modified_BP.cpp circuit_main.cpp Sparse.cpp Sparse.h
		${CPP}  ${hpcc1}  -o FOSD_CBP_hpcc BP.h BP.cpp modified_BP.h modified_BP.cpp circuit_main.cpp Sparse.cpp Sparse.h ${hpcc2} 
		




		

		

