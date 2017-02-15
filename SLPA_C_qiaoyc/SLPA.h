//============================================================================
// Name        : SLPA.h
// Author      : Jierui Xie (xiej2@rpi.edu)
// Date        : Oct. 2011
// Version     :
// Copyright   : All rights reserved.
// Description : SLPA algorithm for community detection.
// Web Site    : https://sites.google.com/site/communitydetectionslpa/
// Publication:
//             J. Xie, B. K. Szymanski and X. Liu, "SLPA: Uncovering Overlapping Communities in Social Networks via A Speaker-listener Interaction Dynamic Process", IEEE ICDM workshop on DMCCI 2011, Vancouver, CA.
//============================================================================
#ifndef SLPA_H_
#define SLPA_H_

#include "Net.h"
#include "NODE.h"
#include <map>
#include <vector>
#include <utility>
#include <tr1/unordered_map>
#include <unordered_map>
#include "MersenneTwister.h"

#include <hash_map>
#include <hash_set>


//---------------------------
//		Multi-threading
//---------------------------
typedef std::tr1::unordered_map<int, int> UOrderedH_INT_INT;

struct thread_data{
	int  startind;
	int  endind;

	int *pIndicator;

	//expect to do sharro copy of the pointers
	vector<vector<int>* > cpm;
	vector<UOrderedH_INT_INT* > vectHTable;
};


class SLPA {
public:
	//---------------------------
	//		network parameters
	//---------------------------
	Net* net;
	string netName;
	string fileName_net;
	string networkPath;

	bool isUseLargestComp; //***
	//---------------------------
	//		SLPA parameters
	//---------------------------
	vector<double> THRS;      //thr
	vector<int> THRCS; 		  //thr count
	bool isSyn;  			  //is synchronous version?
	int maxT;
	int maxRun;
	int version;
	//---------------------------
	//		more
	//---------------------------
	string outputDir;

	MTRand mtrand1;
	MTRand mtrand2;

	SLPA(string, vector<double>,int ,int ,string ,bool,int );
	virtual ~SLPA();

	void start();

	void pre_initial_THRCS();
	void initWQueue_more();

	//void GLPA_asyn();
	void GLPA_asyn_pointer(); // the original version
	

	int ceateHistogram_selRandMax(const vector<int>& wordsList);
	void post_createWQHistogram_MapEntryList();
	void post_thresholding(vector<pair<int,int> >& pairList, int thrc, vector<int>& WS);

	//need to change
	//vector<vector<int> > post_sameLabelDisconnectedComponents(vector<vector<int> >& cpm);
	//static void show_cpm(vector<vector<int> >& cpm);
	static void sort_cpm(vector<vector<int> >& cpm);

	//cpm pointer function
	vector<vector<int>* > post_removeSubset_UorderedHashTable_cpmpointer(vector<vector<int>* >& cpm);
	static void sort_cpm_pointer(vector<vector<int>* >& cpm);

	void write2txt_CPM_pointer(string fileName,vector<vector<int>* >& cpm);
	void post_threshold_createCPM_pointer(int thrc,string fileName);

	void dothreshold_createCPM(int thrc,vector<vector<int> >& cpm);
	void dothreshold_createCPM_pointer(int thrc,vector<vector<int>* >& cpm);

	static bool isDEBUG;


	//---------------------------
	//		Multi-threading
	//---------------------------
	int numThreads;

	void decomposeTasks(int numTasks,int numThd,int stInds[],int enInds[]);
	static void *removesubset_onethread(void *threadarg);
	vector<vector<int>* > post_removeSubset_UorderedHashTable_cpmpointer_MultiThd(vector<vector<int>* >& cpm);

	//---------------------------
	// code add by qiao_yuchen
	//---------------------------
	vector<MTRand> mtrand1s;
	vector<MTRand> mtrand2s;
	SLPA(string, vector<double>, int, int, string, bool, int, int); // new construction function

	void start_time(); // instrument time recorder into the original code
	void start_qiao_v1(); // modified version by qiao_yuchen

	void GLPA_asyn_pointer_time(); // intrument time recoder into the original code

	void GLPA_asyn_pointer_omp_v1(); // add openmp on original code, actually it is synchronized version, please notice that
	void GLPA_asyn_pointer_omp_v2(); // add openmp on original code, and it is non-synchronized version, please notice that
	void GLPA_asyn_pointer_omp_v3(); // add openmp on original code, generate multi-threads out of for-loop of maxT in synchronized way
	void GLPA_asyn_pointer_omp_v4(); // add openmp on original code, generate multi-threads out of for-loop of maxT in non-synchronized way
	void GLPA_asyn_pointer_omp_v5(); // add openmp on original code, multi-threads out of outer for-loop, shared variables on heap while private on each thread's own stack in synchronized way 
	void GLPA_asyn_pointer_omp_v6(); // add openmp on original code, multi-threads out of outer for-loop, shared variables on heap while private on each thread's own stack in non-synchronized way
	void GLPA_asyn_pointer_omp_v7(); // add openmp on original code, multi-threads out of outer for-loop, shared variables on heap while private on each thread's own stack， reduce memory allocate times for length variation of vector in synchronized way 
	
	void GLPA_asyn_pointer_qiao_v0(); // modified version by qiao_yuchen, serial version
	void GLPA_asyn_pointer_qiao_v1(); // modified version by qiao_yuchen
	void GLPA_asyn_pointer_qiao_v2(); // modified version by qiao_yuchen with openmp

	int selectMostFrequentLabel_v1(map<int, int>& labelsList);
	int selectMostFrequentLabel_v2(map<int, int>& labelsList, vector<int>& mostLabelsList);

	int ceateHistogram_selRandMax_qiao_v1(const vector<int>& wordsList); // use multi random generator
	int ceateHistogram_selRandMax_qiao_v2(const vector<int>& wordsList, MTRand& mtrand1_s); // use multi random generator and put them on heap

	void SetWQueueSize();
};

#endif /* SLPA_H_ */
