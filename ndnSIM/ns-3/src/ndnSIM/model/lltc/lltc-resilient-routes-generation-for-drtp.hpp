/*
 * This work is licensed under CC BY-NC-SA 4.0
 * (https://creativecommons.org/licenses/by-nc-sa/4.0/).
 * Copyright (c) 2021 Boyang Zhou
 *
 * This file is a part of "Disruption Resilient Transport Protocol"
 * (https://github.com/zhouby-zjl/drtp/).
 * Written by Boyang Zhou (zhouby@zhejianglab.com)
 *
 * This software is protected by the patents numbered with PCT/CN2021/075891,
 * ZL202110344405.7 and ZL202110144836.9, as well as the software copyrights
 * numbered with 2020SR1875227 and 2020SR1875228.
 */

#ifndef RESILIENT_ROUTES_GENERATION_HPP_
#define RESILIENT_ROUTES_GENERATION_HPP_

#include "ns3/ndnSIM/model/lltc/lltc-resilient-routes-generation.hpp"


namespace lltc {

class ResilientRouteGenerationForDRTP : public ResilientRouteGenerationGeneric {
public:
	LltcResilientRouteVectors* genResilientRoutesVectors(LltcGraph* g, int srcNodeId, LltcConfiguration* config);
	ResilientRoutes* constructResilientRoutes(LltcGraph* g, LltcResilientRouteVectors* rrv, int dstNodeId,
																LltcConfiguration* config);
	vector<double>* evaluateResilientRouteEEFR(LltcGraph* g, ResilientRoutes* rr,
							LltcConfiguration* config, double linkLossRate_min, double linkLossRate_max,
							double linkLossRate_step, double eta);

	double computeResilientRouteEEFR(LltcGraph* g, ResilientRoutes* rr, LltcConfiguration* config);
	double computeResilientRouteEEDT_N_InStep(LltcGraph* g, ResilientRoutes* rr, LltcConfiguration* config);
	double computeResilientRouteEEDT_N_OutOfStep(LltcGraph* g, ResilientRoutes* rr, LltcConfiguration* config, int m);
	double computeResilientRouteEEDT_R_OutOfStep(LltcGraph* g, ResilientRoutes* rr, LltcConfiguration* config, int m, size_t reseqQueueSize);

private:
	vector<RedundantSubPath*>* genRedundantSubPaths(LltcGraph* g, int srcNodeId, vector<int>* dstNodeIds,
					unordered_set<int>* invalidLinkIds, LltcConfiguration* config, double maxAllowedRetransDelay);
	RedundantSubPath* genPath(LltcGraph* g, int srcNodeId, vector<int>* dstNodeIds,
			unordered_set<int>* expelledLinkIds, unordered_set<int>* invalidLinkIds,
			LltcConfiguration* config, double maxAllowedRetransDelay);

	unordered_set<int>* convertPathNodeIdsArrayToLinkIdsSet(LltcGraph* g, vector<int>* retransDstNodeIds);
	vector<int>* convertPrevArrayToPathNodeIdsArray(int* prevNodeIdx, int srcNodeId, int uNodeId);
	double computePHI(LltcLink* e, LltcNode* uNode, LltcNode* vNode, vector<RedundantSubPath*>* rsps,
							unordered_set<RedundantSubPath*>* rsps_invalid, int nACK);


	double computeMaxRetransDelay(vector<RedundantSubPath*>* rsps, unordered_set<RedundantSubPath*>* rsps_invalid,
							int primRetransDelayInUs, LltcConfiguration* config);
	double computeSubPathRetransDelay(LltcGraph* g, RedundantSubPath* sp, int numRetransRequests);
	bool isInArray(vector<RedundantSubPath*>* rsps, RedundantSubPath* rsp);

	double* computeProbVectorSelectingSubPath(LltcGraph* g, vector<RedundantSubPath*>* X, LltcConfiguration* config);  // return an array of double values


	double pruningResilientRoutes(LltcGraph* g, ResilientRoutes* rr, LltcConfiguration* config);
	//void addFibsFromPath(Path* path, unordered_map<int, FIB*>* fibMap, LltcGraph* g, int pathId);

	double computePmuPiatDriftRange(double freq, double freqDriftRatio, int maxConsecutiveDriftPackets);

	void computeRetranTimes(LltcGraph* g, ResilientRoutes* rr);


	static int _nextPathID;

};


}

#endif /* RESILIENT_ROUTES_GENERATION_HPP_ */
