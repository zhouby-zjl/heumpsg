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

#include <iostream>
#include <algorithm>
#include <unordered_set>
#include "lltc-fibonacci-heap-key-funcs.hpp"
#include "lltc-resilient-routes-generation-for-drtp.hpp"

using namespace std;
using namespace lltc;

int ResilientRouteGenerationForDRTP::_nextPathID = 0;

vector<RedundantSubPath*>* ResilientRoutes::getRSPsByNodeId(int nodeId) {
	if (hrsp->find(nodeId) == hrsp->end()) return NULL;
	return (*hrsp)[nodeId];
}

void ResilientRoutes::setRSPsToNodeId(int nodeId, vector<RedundantSubPath*>* rsp) {
	(*hrsp)[nodeId] = rsp;
}

unordered_map<int, vector<RedundantSubPath*>*>* ResilientRoutes::getHRSP() {
	return hrsp;
}

void ResilientRoutes::removeRSPByNodeId(int nodeId) {
	hrsp->erase(nodeId);
}

LltcResilientRouteVectors* ResilientRouteGenerationForDRTP::genResilientRoutesVectors(LltcGraph* g, int srcNodeId,
														LltcConfiguration* config) {
	int n = g->getNNodes();
	LltcResilientRouteVectors* rrv = new LltcResilientRouteVectors(srcNodeId, n);

	double PSIs[n];
	double PHIs[n];
	double ALPHAs[n];
	double GOALs[n];

	bool isNodeIdInQ[n];
	for (int i = 0; i < n; ++i) {
		isNodeIdInQ[i] = true;
		PSIs[i] = 1.0;
		PHIs[i] = 0;
		ALPHAs[i] = 1.0;

		rrv->prevNodeIdx[i] = -1;
		rrv->expectedLossRates[i] = 1.0;
		rrv->rspsMap[i] = NULL;
		rrv->maxPrimPathDelays[i] = 0.0;
		rrv->maxDelays[i] = 0.0;

		GOALs[i] = 1.0;
	}

	rrv->expectedLossRates[srcNodeId] = 0.0;
	rrv->maxDelays[srcNodeId] = computePmuPiatDriftRange(config->pmuFreq, config->maxDriftRatioForPmuFreq,
														config->maxConsecutiveDriftPackets);
	GOALs[srcNodeId] = 0.0;
	PSIs[srcNodeId] = 1.0;
	PHIs[srcNodeId] = 0.0;
	ALPHAs[srcNodeId] = 1.0;

	//FHKeyFuncLossRate keyFunc(rrv->expectedLossRates);
	FHKeyFuncLossRate keyFunc(GOALs);
	FibonacciHeap q(&keyFunc);

	for (int i = 0; i < n; ++i) {
		q.insert((FHObjectPtr) g->nodes[i]);
	}

	while (q.size() > 0) {
		//cout << q.size() << endl;
		FibonacciHeapNode* heapNode = q.removeMax();
		int uNodeId = *((int*) (heapNode->getNodeObject()));
		LltcNode* uNode = g->getNodeById(uNodeId);
		isNodeIdInQ[uNode->nodeId] = false;

		for (unsigned int i = 0; i < uNode->eLinks->size(); ++i) {
			LltcLink* e = (*(uNode->eLinks))[i];
			LltcNode* vNode = LltcGraph::getLinkOpSide(e, uNode);
			int vNodeId = vNode->nodeId;

			if (!isNodeIdInQ[vNodeId]) continue;

			int primPathLinkDelay = LltcGraph::getLinkDelayTotal(e, uNode, vNode);
			int primPathLinkDelay_up = LltcGraph::getLinkDelayTotal(e, vNode, uNode);
			double primPathDelay = rrv->maxPrimPathDelays[uNodeId] + primPathLinkDelay;
			if (primPathDelay > config->maxPathDelay) continue;

			vector<int>* retransDstNodeIds = convertPrevArrayToPathNodeIdsArray(rrv->prevNodeIdx,
																				srcNodeId, uNode->nodeId);

			unordered_set<int>* invalidLinkIds = convertPathNodeIdsArrayToLinkIdsSet(g, retransDstNodeIds);
			invalidLinkIds->insert(e->linkId);

			double curRrMaxDelay = rrv->maxDelays[uNodeId] + primPathLinkDelay;
			double maxAllowedRetransDelay = config->maxPathDelay - curRrMaxDelay;
			vector<RedundantSubPath*>* rsps = genRedundantSubPaths(g, vNode->nodeId, retransDstNodeIds,
										invalidLinkIds, config, maxAllowedRetransDelay);

			double PHI, PSI, F_overV;

			PHI = computePHI(e, uNode, vNode, rsps, NULL, config->numRetransRequests);
			PSI = PSIs[uNodeId] * (1 - PHIs[uNodeId]);

			F_overV = rrv->expectedLossRates[uNode->nodeId] + PHI * PSI;


			double alpha_overV = ALPHAs[uNodeId] * (1 - pow(1 - LltcGraph::getLinkRelia(e, uNode, vNode, false),
					config->numDataRetransReportsToSend));

			//double all_cur = ALPHAs[vNodeId] * rrv->expectedLossRates[vNode->nodeId] + (1 - ALPHAs[vNodeId]);
			double all_overV = alpha_overV * F_overV + (1 - alpha_overV);

			//if (F_overV < rrv->expectedLossRates[vNode->nodeId]) {
			//if (all_overV < all_cur) {
			if (all_overV < GOALs[vNodeId]) {
				rrv->maxPrimPathDelays[vNodeId] = primPathDelay;
				int primLinkRetransDelayInUs = (primPathLinkDelay_up + primPathLinkDelay) * config->numRetransRequests;

				rrv->maxDelays[vNodeId] = curRrMaxDelay + computeMaxRetransDelay(rsps, NULL, primLinkRetransDelayInUs, config);

				rrv->expectedLossRates[vNodeId] = F_overV;
				ALPHAs[vNodeId] = alpha_overV;
				GOALs[vNodeId] = all_overV;
				q.updateKey((FHObjectPtr) &(vNode->nodeId), (FHObjectPtr) &(vNode->nodeId));

				if (rsps->size() > 0)
					rrv->rspsMap[vNodeId] = rsps;

				PHIs[vNodeId] = PHI;
				PSIs[vNodeId] = PSI;

				rrv->prevNodeIdx[vNodeId] = uNode->nodeId;
			}


		}
	}



	return rrv;
}

ResilientRoutes* ResilientRouteGenerationForDRTP::constructResilientRoutes(LltcGraph* g, LltcResilientRouteVectors* rrv,
															int dstNodeId, LltcConfiguration* config) {
	vector<int>* nodeIds = convertPrevArrayToPathNodeIdsArray(rrv->prevNodeIdx, rrv->srcNodeId, dstNodeId);
	if (nodeIds->size() <= 1) {
		return NULL;
	}
	int base_pathID = 1;
	_nextPathID = base_pathID;
	Path* primaryPath = new Path();
	primaryPath->nodeIds = nodeIds;
	primaryPath->pathID = _nextPathID++;

	unordered_set<int> linkIDInPP;
	for (size_t i = 1; i < nodeIds->size(); ++i) {
		LltcLink* l_pp = g->getLinkByNodeIds((*nodeIds)[i - 1], (*nodeIds)[i]);
		linkIDInPP.insert(l_pp->linkId);
	}


	unordered_map<int, vector<RedundantSubPath*>*>* hrsp = new unordered_map<int, vector<RedundantSubPath*>*>();
	for (unsigned int i = 0; i < nodeIds->size(); ++i) {
		int nodeId = (*nodeIds)[i];
		vector<RedundantSubPath*>* rsps = rrv->rspsMap[nodeId];
		if (rsps == NULL) {
			(*hrsp)[nodeId] = new vector<RedundantSubPath*>();
			continue;
		}

		vector<int> rsps_overlapped;
		vector<int> rsps_removed;
		size_t j = 0;
		for (vector<RedundantSubPath*>::iterator iter = rsps->begin();
				iter != rsps->end(); ++iter) {
			RedundantSubPath* rsp_uNode = *iter;
			vector<int>* nodeIds = rsp_uNode->path->nodeIds;
			bool isLinkInRSP = false;
			for (int m = 1; m < nodeIds->size() - 1; ++m) {
				LltcLink* e_RSP = g->getLinkByNodeIds((*nodeIds)[m - 1], (*nodeIds)[m]);
				if (linkIDInPP.find(e_RSP->linkId) != linkIDInPP.end()) {
					isLinkInRSP = true;
					break;
				}
			}

			if (isLinkInRSP) {
				rsps_overlapped.push_back(j);
			}


			rsp_uNode->path->pathID = _nextPathID++;
			++j;
		}

		for (vector<int>::reverse_iterator iter = rsps_overlapped.rbegin(); iter != rsps_overlapped.rend(); ++iter) {
			//rsps->erase(rsps->begin() + *iter);
		}


		(*hrsp)[nodeId] = rsps;
	}


	ResilientRoutes* rr = new ResilientRoutes(primaryPath, hrsp);

	computeRetranTimes(g, rr);
	rr->maxRrDelay = computeResilientRouteEEDT_N_InStep(g, rr, config);
	rr->expectedRrLossRate = computeResilientRouteEEFR(g, rr, config);

	if (LltcConfig::LLTC_DISABLE_RETRAN) {
		for (unsigned int i = 1; i < nodeIds->size(); ++i) {
			vector<RedundantSubPath*>* rsps = rr->getRSPsByNodeId((*nodeIds)[i]);

		}
	}

	return rr;
}

void ResilientRouteGenerationForDRTP::computeRetranTimes(LltcGraph* g, ResilientRoutes* rr) {
	Path* primaryPath = rr->primaryPath;
	unordered_map<int, vector<RedundantSubPath*>*>* hrsp = rr->getHRSP();
	vector<int>* nodeIds = primaryPath->nodeIds;
	size_t n = nodeIds->size();
	int maxCulQueueTimesInUs[n];
	for (size_t i = 0; i < n; ++i) maxCulQueueTimesInUs[i] = 0;

	(*rr->retransTimeouts)[0] = 0;
	int culQueueTimeInUs = 0;
	rr->primPathDelay = 0;
	for (unsigned int i = 1; i < nodeIds->size(); ++i) {
		LltcLink* l = g->getLinkByNodeIds((*nodeIds)[i - 1], (*nodeIds)[i]);
		culQueueTimeInUs += g->getLinkDelayQueue(l, (*nodeIds)[i - 1], (*nodeIds)[i]);
		rr->primPathDelay += g->getLinkDelayTotal(l, (*nodeIds)[i - 1], (*nodeIds)[i]);
		maxCulQueueTimesInUs[i] = culQueueTimeInUs;

		int nodeId = (*nodeIds)[i];
		vector<RedundantSubPath*>* rsps = (*hrsp)[nodeId];
		if (rsps == NULL || rsps->size() == 0) continue;
		int maxTimeoutForRetrans = -1;
		for (RedundantSubPath* rsp : *rsps) {
			if (rsp->retransDelayInMultiTimes > maxTimeoutForRetrans) {
				maxTimeoutForRetrans = rsp->retransDelayInMultiTimes;
			}
		}
		(*rr->retransTimeouts)[i] = maxTimeoutForRetrans;
	}

	(*rr->waitingTimeoutOnNextPacketInUs)[0] = 0;
	(*rr->waitingTimeoutOnNextPacketForLostReportInUs)[0] = 0;
	int culRetransTimeout = 0;
	int maxPrimPathDelays = 0;
	for (unsigned int i = 1; i < nodeIds->size(); ++i) {
		culRetransTimeout += (*rr->retransTimeouts)[i - 1];

		(*rr->waitingTimeoutOnNextPacketInUs)[i] = maxCulQueueTimesInUs[i];
		LltcLink* e = g->getLinkByNodeIds((*nodeIds)[i - 1], (*nodeIds)[i]);
		maxPrimPathDelays += e->abDelay_total;

		(*rr->waitingTimeoutOnNextPacketForLostReportInUs)[i] = maxPrimPathDelays +
											culRetransTimeout;
								//(double) ((*rr->retransTimeouts)[i]) * (1.0 / (double) config->numRetransRequests);

		int culRetransTimeoutAtNodeJ = 0;
		for (unsigned int j = i - 1; j >= 1; --j) {
			culRetransTimeoutAtNodeJ += (*rr->retransTimeouts)[j];
			(*(rr->waitingTimeoutOnRetransReportInUs[i]))[(*nodeIds)[j]] = culRetransTimeoutAtNodeJ +
								maxCulQueueTimesInUs[i] - maxCulQueueTimesInUs[j];
		}
	}



}


vector<RedundantSubPath*>* ResilientRouteGenerationForDRTP::genRedundantSubPaths(LltcGraph* g, int srcNodeId, vector<int>* dstNodeIds,
				unordered_set<int>* invalidLinkIds, LltcConfiguration* config, double maxAllowedRetransDelay) {
	vector<RedundantSubPath*>* r = new vector<RedundantSubPath*>();
	unordered_set<int> expelledLinkIds;

	for (int i = 0; i < config->beta; ++i) {
		RedundantSubPath* sp = genPath(g, srcNodeId, dstNodeIds, &expelledLinkIds,
					invalidLinkIds, config, maxAllowedRetransDelay);
		if (sp == NULL || isInArray(r, sp)) break;
		r->push_back(sp);

		if (i < config->beta - 1) {
			for (unsigned int j = 1; j < sp->path->nodeIds->size(); ++j) {
				LltcLink* l = g->getLinkByNodeIds((*(sp->path->nodeIds))[j - 1], (*(sp->path->nodeIds))[j]);
				expelledLinkIds.insert(l->linkId);
			}
		}
	}
	return r;
}

RedundantSubPath* ResilientRouteGenerationForDRTP::genPath(LltcGraph* g, int srcNodeId, vector<int>* dstNodeIds,
				unordered_set<int>* expelledLinkIds, unordered_set<int>* invalidLinkIds,
				LltcConfiguration* config, double maxAllowedRetransDelay) {
	size_t nNodes = g->getNNodes();
	LltcNode* virDstNode = g->addVirtNode(dstNodeIds);
	int n = nNodes + 1;
	double logNodeRelias[n];
	double _logNodeRelias[n];
	int upDelays[n];
	int downDelays[n];
	int retransDelaysInMultiTimes[n];
	bool isNodeIdInQ[n];
	int prevNodeIdx[n];

	for (int i = 0; i < n; ++i) {
		logNodeRelias[i] = DOUBLE_NEGATIVE_INFINITY;
		_logNodeRelias[i] = DOUBLE_NEGATIVE_INFINITY;
		isNodeIdInQ[i] = true;
		upDelays[i] = 0;
		downDelays[i] = 0;
		retransDelaysInMultiTimes[i] = 0;
		prevNodeIdx[i] = -1;
	}

	logNodeRelias[srcNodeId] = 0.0;
	_logNodeRelias[srcNodeId] = 0.0;
	upDelays[srcNodeId] = 0;
	downDelays[srcNodeId] = 0;

	FHKeyFuncNodeRelias keyFunc(logNodeRelias);
	FibonacciHeap q(&keyFunc);

	for (int i = 0; i < n; ++i) {
		q.insert((FHObjectPtr) g->nodes[i]);
	}

	double logExpellingRelia = g->getMinLinkRelia(true) - log(3);

	while (q.size() > 0) {
		int uNodeId = *((int*) q.removeMaxObject());

		LltcNode* uNode = g->getNodeById(uNodeId);
		isNodeIdInQ[uNodeId] = false;

		for (LltcLink* e : *uNode->eLinks) {
			if (invalidLinkIds->find(e->linkId) != invalidLinkIds->end()) continue;

			LltcNode* vNode = LltcGraph::getLinkOpSide(e, uNode);
			int vNodeId = vNode->nodeId;
			int* vNodeIdPtr = &(vNode->nodeId);
			if (!isNodeIdInQ[vNodeId]) continue;

			int upDelay = LltcGraph::getLinkDelayTotal(e, uNode, vNode) + upDelays[uNodeId];
			int downDelay = LltcGraph::getLinkDelayTotal(e, vNode, uNode) + downDelays[uNodeId];
			double maxRetransDelay = retransDelaysInMultiTimes[uNodeId] +
					(double) (e->abDelay_total + e->baDelay_total) * (double) (config->numRetransRequests);

			if (maxRetransDelay > maxAllowedRetransDelay) 	continue;

			double logLinkRelia;
			double _logLinkRelia = LltcGraph::getLinkRelia(e, uNode, vNode, true) + LltcGraph::getLinkRelia(e, vNode, uNode, true);
			if (expelledLinkIds->find(e->linkId) == expelledLinkIds->end()) {
				logLinkRelia = _logLinkRelia;
			} else {
				logLinkRelia = logExpellingRelia;
			}

			double logReliaOverV = logLinkRelia + logNodeRelias[uNodeId];
			double prevLogReliaOverV = logNodeRelias[vNodeId];
			if (logReliaOverV > prevLogReliaOverV) {
				prevNodeIdx[vNodeId] = uNode->nodeId;

				logNodeRelias[vNodeId] = logReliaOverV;
				q.updateKey((FHObjectPtr) vNodeIdPtr, (FHObjectPtr) vNodeIdPtr);

				_logNodeRelias[vNodeId] = _logLinkRelia + _logNodeRelias[uNodeId];
				upDelays[vNodeId] = upDelay;
				downDelays[vNodeId] = downDelay;
				retransDelaysInMultiTimes[vNodeId] = maxRetransDelay;
			}
		}
	}

	if (prevNodeIdx[virDstNode->nodeId] == -1) {
		g->removeVirtNode();
		return NULL;
	}

	int pathNodes = 1;
	int v = virDstNode->nodeId;
	while (v != srcNodeId) {
		v = prevNodeIdx[v];
		++pathNodes;
	}

	Path* path = new Path();
	path->setLength(pathNodes - 1);
	(*path->nodeIds)[0] = srcNodeId;
	v = virDstNode->nodeId;
	int p = pathNodes - 1;
	do {
		if (p < pathNodes - 1) {
			(*path->nodeIds)[p] = v;
		}
		v = prevNodeIdx[v];
		--p;
	} while (v != srcNodeId);

	(*path->nodeIds)[0] = srcNodeId;

	RedundantSubPath* rsp = new RedundantSubPath(path, 0, pathNodes - 2);
	rsp->oneWayDataDeliveryLossRate = 1 - pow(exp(1), _logNodeRelias[virDstNode->nodeId]);
	rsp->dataDeliveryLossRate = pow(rsp->oneWayDataDeliveryLossRate, config->numRetransRequests);
	rsp->linkUpDelay = upDelays[virDstNode->nodeId];
	rsp->linkDownDelay = downDelays[virDstNode->nodeId];
	rsp->retransDelayInMultiTimes = retransDelaysInMultiTimes[virDstNode->nodeId];
	rsp->retransDelayInSingleTime = rsp->retransDelayInMultiTimes / config->numRetransRequests;

	g->removeVirtNode();

	return rsp;
}

unordered_set<int>* ResilientRouteGenerationForDRTP::convertPathNodeIdsArrayToLinkIdsSet(LltcGraph* g, vector<int>* retransDstNodeIds) {
	unordered_set<int>* linkIdSet = new unordered_set<int>();
	for (unsigned int i = 1; i < retransDstNodeIds->size(); ++i) {
		int linkId = g->getLinkByNodeIds((*retransDstNodeIds)[i], (*retransDstNodeIds)[i - 1])->linkId;
		linkIdSet->insert(linkId);
	}
	return linkIdSet;
}

vector<int>* ResilientRouteGenerationForDRTP::convertPrevArrayToPathNodeIdsArray(int* prevNodeIdx, int srcNodeId, int uNodeId) {
	int pathNodes = 1;
	int v = uNodeId;
	while (v != srcNodeId) {
		v = prevNodeIdx[v];
		if (v == -1) break;
		++pathNodes;
	}

	vector<int>* path = new vector<int>();
	for (int i = 0; i < pathNodes; ++i) {
		path->push_back(0);
	}

	v = uNodeId;
	int p = pathNodes - 1;
	do {
		(*path)[p] = v;
		v = prevNodeIdx[v];
		--p;
	} while (p >= 0 && v != srcNodeId);

	if (pathNodes >= 2)
		(*path)[0] = srcNodeId;
	return path;
}


double ResilientRouteGenerationForDRTP::computePHI(LltcLink* e, LltcNode* uNode, LltcNode* vNode,
					vector<RedundantSubPath*>* rsps, unordered_set<RedundantSubPath*>* rsps_invalid, int numRetransRequests) {
	double r = 1.0;
	for (RedundantSubPath* rsp : *rsps) {
		if (rsps_invalid != NULL && rsps_invalid->find(rsp) == rsps_invalid->end()) {
			continue;
		}
		r *= rsp->dataDeliveryLossRate;
	}
	double uvRelia = LltcGraph::getLinkRelia(e, uNode, vNode, false);
	double vuRelia = LltcGraph::getLinkRelia(e, vNode, uNode, false);
	r *= (1.0 - uvRelia) * pow(1.0 - uvRelia * vuRelia, numRetransRequests);
	return r;
}


double ResilientRouteGenerationForDRTP::computeResilientRouteEEDT_N_InStep(LltcGraph* g, ResilientRoutes* rr, LltcConfiguration* config) {
	int nHopPrimPath = rr->primaryPath->nodeIds->size();
	double delay = computePmuPiatDriftRange(config->pmuFreq, config->maxDriftRatioForPmuFreq, config->maxConsecutiveDriftPackets);

	for (int i = 1; i < nHopPrimPath; ++i) {
		int nodeA_id = (*(rr->primaryPath->nodeIds))[i - 1];
		int nodeB_id = (*(rr->primaryPath->nodeIds))[i];
		LltcLink* l = g->getLinkByNodeIds(nodeA_id, nodeB_id);

		int primLinkDelay_down = g->getLinkDelayTotal(l, nodeA_id, nodeB_id);
		int primLinkQueueDelay_down = g->getLinkDelayQueue(l, nodeA_id, nodeB_id);

		int maxTimeoutForRetransInUs = (l->abDelay_total + l->baDelay_total) * (config->numRetransRequests - 1);
		vector<RedundantSubPath*>* X = rr->getRSPsByNodeId(nodeB_id);
		if (X != NULL && X->size() > 1) {
			for (unsigned int s = 0; s < X->size(); ++s) {
				RedundantSubPath* sp = (*X)[s];
				if (sp->retransDelayInMultiTimes > maxTimeoutForRetransInUs) {
					maxTimeoutForRetransInUs = sp->retransDelayInMultiTimes;
				}
			}
		}

		delay += primLinkDelay_down + primLinkQueueDelay_down + maxTimeoutForRetransInUs;
	}

	return delay;
}

double ResilientRouteGenerationForDRTP::computeResilientRouteEEDT_N_OutOfStep(LltcGraph* g, ResilientRoutes* rr,
							LltcConfiguration* config, int m) {
	int nHopPrimPath = rr->primaryPath->nodeIds->size();
	double delay = computePmuPiatDriftRange(config->pmuFreq, config->maxDriftRatioForPmuFreq, config->maxConsecutiveDriftPackets);
	double delay_pp = 0.0;
	double delay_theta_prev = 0.0;
	double delay_theta_last = 0.0;
	double delay_pp_link_last_all_up = 0.0;
	double delay_pp_link_last_all_down = 0.0;
	double delay_pp_link_last_queue_down = 0.0;

	for (int i = 1; i < nHopPrimPath; ++i) {
		int nodeA_id = (*(rr->primaryPath->nodeIds))[i - 1];
		int nodeB_id = (*(rr->primaryPath->nodeIds))[i];
		LltcLink* l = g->getLinkByNodeIds(nodeA_id, nodeB_id);

		int primLinkQueueDelay_down = g->getLinkDelayQueue(l, nodeA_id, nodeB_id);

		int maxTimeoutForRetransInUs = (l->abDelay_total + l->baDelay_total) * (config->numRetransRequests - 1);
		vector<RedundantSubPath*>* X = rr->getRSPsByNodeId(nodeB_id);
		if (X != NULL && X->size() > 1) {
			for (unsigned int s = 0; s < X->size(); ++s) {
				RedundantSubPath* sp = (*X)[s];
				if (sp->retransDelayInMultiTimes > maxTimeoutForRetransInUs) {
					maxTimeoutForRetransInUs = sp->retransDelayInMultiTimes;
				}
			}
		}

		delay_pp += primLinkQueueDelay_down;

		if (i < nHopPrimPath - 1) {
			delay_theta_prev += maxTimeoutForRetransInUs;
		} else {
			delay_theta_last = maxTimeoutForRetransInUs;
			delay_pp_link_last_all_up = g->getLinkDelayTotal(l, nodeB_id, nodeA_id);
			delay_pp_link_last_all_down = g->getLinkDelayTotal(l, nodeA_id, nodeB_id);
			delay_pp_link_last_queue_down = primLinkQueueDelay_down;
		}
	}

	delay += m * (delay_pp_link_last_all_up + delay_pp_link_last_all_down + delay_pp + delay_theta_prev);
	delay += (m + 1) * delay_pp_link_last_queue_down  + delay_theta_last;

	return delay;
}

double ResilientRouteGenerationForDRTP::computeResilientRouteEEDT_R_OutOfStep(LltcGraph* g, ResilientRoutes* rr,
						LltcConfiguration* config, int m, size_t reseqQueueSize) {
	double eedt_n = this->computeResilientRouteEEDT_N_OutOfStep(g, rr, config, m);
	double eedt_r = (double) (1000000 * config->pmuFreq) / (double) reseqQueueSize + eedt_n;
	return eedt_r;
}

double ResilientRouteGenerationForDRTP::computePmuPiatDriftRange(double freq, double freqDriftRatio, int maxConsecutiveDriftPackets) {
	return  (double) maxConsecutiveDriftPackets * 2 * freqDriftRatio * (1000000.0 / freq) / (1 - freqDriftRatio * freqDriftRatio);
}

double ResilientRouteGenerationForDRTP::computeMaxRetransDelay(vector<RedundantSubPath*>* rsps,
		unordered_set<RedundantSubPath*>* rsps_invalid,
		int primRetransDelayInUs, LltcConfiguration* config) {
	double maxRetransDelay = (double) primRetransDelayInUs;
	for (RedundantSubPath* sp : *rsps) {
		if (rsps_invalid != NULL && rsps_invalid->find(sp) != rsps_invalid->end()) {
			continue;
		}
		if (sp->retransDelayInMultiTimes > maxRetransDelay) {
			maxRetransDelay = sp->retransDelayInMultiTimes;
		}
	}
	return maxRetransDelay;
}

double ResilientRouteGenerationForDRTP::computeSubPathRetransDelay(LltcGraph* g, RedundantSubPath* sp, int nACK) {
	double retransDelay = 0.0;

	for (int q = 0; q < nACK - 1; ++q) {
		retransDelay += pow(sp->oneWayDataDeliveryLossRate, q) * (1.0 - (double)sp->oneWayDataDeliveryLossRate) *
				(double)(sp->retransDelayInMultiTimes * q + sp->linkUpDelay + sp->linkDownDelay);
	}
	return retransDelay;
}

bool ResilientRouteGenerationForDRTP::isInArray(vector<RedundantSubPath*>* rsps, RedundantSubPath* rsp) {
	for (RedundantSubPath* _rsp : *rsps) {
		if (_rsp->path->nodeIds->size() != rsp->path->nodeIds->size()) continue;
		bool isTheSame = true;
		for (unsigned int i = 0; i < _rsp->path->nodeIds->size(); ++i) {
			if ((*(_rsp->path->nodeIds))[i] != (*(rsp->path->nodeIds))[i]) {
				isTheSame = false;
				break;
			}
		}
		if (isTheSame) return true;
	}
	return false;
}

double ResilientRouteGenerationForDRTP::computeResilientRouteEEFR(LltcGraph* g, ResilientRoutes* rr, LltcConfiguration* config) {
	int nHopPrimPath = rr->primaryPath->nodeIds->size();
	double lossRate = 0.0;
	double PSI = 0.0;
	double PSI_prev = 1.0;
	double PHI = 0.0;
	double PHI_prev = 0.0;

	for (int i = 1; i < nHopPrimPath; ++i) {
		int uNodeId = (*(rr->primaryPath->nodeIds))[i - 1];
		int vNodeId = (*(rr->primaryPath->nodeIds))[i];
		LltcNode* uNode = g->getNodeById(uNodeId);
		LltcNode* vNode = g->getNodeById(vNodeId);
		LltcLink* e = g->getLinkByNodeIds(uNodeId, vNodeId);

		vector<RedundantSubPath*>* rsps = rr->getRSPsByNodeId(vNodeId);
		if (rsps != NULL) {
			PHI = this->computePHI(e, uNode, vNode, rsps, NULL, config->numRetransRequests);
		} else {
			double uvRelia = LltcGraph::getLinkRelia(e, uNode, vNode, false);
			double vuRelia = LltcGraph::getLinkRelia(e, vNode, uNode, false);
			PHI = (1.0 - uvRelia) * pow(1.0 - uvRelia * vuRelia, config->numRetransRequests);
		}

		PSI = PSI_prev * (1 - PHI_prev);
		lossRate += PSI * PHI;

		PSI_prev = PSI;
		PHI_prev = PHI;
	}

	return lossRate;
}


vector<double>* ResilientRouteGenerationForDRTP::evaluateResilientRouteEEFR(LltcGraph* g, ResilientRoutes* rr,
						LltcConfiguration* config, double linkLossRate_min, double linkLossRate_max,
						double linkLossRate_step, double eta) {
	vector<double>* lossRates = new vector<double>();

	int nTests = (linkLossRate_max - linkLossRate_min) / linkLossRate_step;

	for (int k = 0; k < nTests; ++k) {
		double linkLossRate = linkLossRate_min + k * linkLossRate_step;
		int nHopPrimPath = rr->primaryPath->nodeIds->size();
		double lossRate = 0.0;
		double PSI = 0.0;
		double PSI_prev = 1.0;
		double PHI = 0.0;
		double PHI_prev = 0.0;

		for (int i = 1; i < nHopPrimPath; ++i) {
			int vNodeId = (*(rr->primaryPath->nodeIds))[i];

			vector<RedundantSubPath*>* rsps = rr->getRSPsByNodeId(vNodeId);
			PHI = 1.0;
			if (rsps != NULL) {
				for (RedundantSubPath* rsp : *rsps) {
					double dataDeliveryLossRate = pow(1.0 - pow(1.0 - linkLossRate, (rsp->path->nodeIds->size() - 1) * 2),
							eta);
					PHI *= dataDeliveryLossRate;
				}
			}
			double linkRelia = 1.0 - linkLossRate;
			PHI *= linkLossRate * pow(1.0 - linkRelia * linkRelia, eta);

			PSI = PSI_prev * (1 - PHI_prev);
			lossRate += PSI * PHI;

			PSI_prev = PSI;
			PHI_prev = PHI;
		}

		lossRates->push_back(lossRate);
	}

	return lossRates;
}

double* ResilientRouteGenerationForDRTP::computeProbVectorSelectingSubPath(LltcGraph* g, vector<RedundantSubPath*>* X, LltcConfiguration* config) {
	int s = X->size();
	double* prVec = new double[s];
	double lossRateVec[s];

	for (int i = 0; i < s; ++i) {
		RedundantSubPath* sp = (*X)[i];
		lossRateVec[i] = sp->dataDeliveryLossRate;
	}

	double probLossRateForPrevSubPaths = 1.0;
	prVec[0] = 1 - lossRateVec[0];
	for (int i = 1; i < s; ++i) {
		probLossRateForPrevSubPaths *= lossRateVec[i - 1];
		prVec[i] = probLossRateForPrevSubPaths * (1 - lossRateVec[i]);
	}

	return prVec;
}

double ResilientRouteGenerationForDRTP::pruningResilientRoutes(LltcGraph* g, ResilientRoutes* rr, LltcConfiguration* config) {
	FHKeyFuncRSPTimeoutForRetrans keyFunc;
	FibonacciHeap allRSPs(&keyFunc);

	for (unordered_map<int, vector<RedundantSubPath*>*>::iterator iter = rr->getHRSP()->begin(); iter != rr->getHRSP()->end(); ++iter) {
		vector<RedundantSubPath*>* rsps = iter->second;
		for (RedundantSubPath *rsp : *rsps) {
			allRSPs.insert((FHObjectPtr) rsp);
		}
	}

	double maxDelay = 0.0;
	while (allRSPs.size() >= 1) {
		maxDelay = computeResilientRouteEEDT_N_InStep(g, rr, config);
		if (maxDelay <= config->maxPathDelay) break;

		RedundantSubPath* maxRSP = (RedundantSubPath*) allRSPs.removeMax()->getNodeObject();

		int srcNodeId = (*maxRSP->path->nodeIds)[maxRSP->retransSourceNodeIdx];
		vector<RedundantSubPath*>* rsps = rr->getRSPsByNodeId(srcNodeId);
		if (rsps == NULL) continue;

		vector<RedundantSubPath*>::iterator maxRSP_iter = find(rsps->begin(),  rsps->end(), maxRSP);
		rsps->erase(maxRSP_iter);

		if (rsps->size() == 0) {
			rr->removeRSPByNodeId(srcNodeId);
		}
	}

	if (allRSPs.size() == 0) {
		maxDelay = rr->primPathDelay;
	}

	return maxDelay;
}

Path::Path() {
	nodeIds = new vector<int>();
	pathID = 0;
}

void Path::setLength(int n) {
	nodeIds = new vector<int>();

	for (int i = 0; i < n; ++i) {
		nodeIds->push_back(0);
	}
}

string* Path::toString() {
	stringstream ss;
	ss << pathID << ": " << nodeIds->size() << ", (";
	for (int nodeId : *nodeIds) {
		ss << nodeId << " ";
	}
	ss << ")";
	string* r = new string(ss.str());
	return r;
}

RedundantSubPath::RedundantSubPath(Path* path, int retransSourceNodeIdx, int retransDstNodeIdx) {
	this->path = path;
	this->retransSourceNodeIdx = retransSourceNodeIdx;
	this->retransDstNodeIdx = retransDstNodeIdx;
}

string* RedundantSubPath::toString() {
	stringstream ss;

	ss << path->pathID << ": " << path->nodeIds->size() << ", (";
	for (int i = retransSourceNodeIdx; i <= retransDstNodeIdx; ++i) {
		ss << (*path->nodeIds)[i] << " ";
	}
	ss << "), ";

	ss << "linkDownDelay: " << linkDownDelay << ", ";
	ss << "linkUpDelay: " << linkUpDelay << ", ";
	ss << "dataDeliveryLossRate: " << dataDeliveryLossRate << ", ";
	ss << "timeoutForRetrans: " << retransDelayInMultiTimes;

	string* s = new string(ss.str());
	return s;
}

ResilientRoutes::ResilientRoutes(Path* primaryPath, unordered_map<int, vector<RedundantSubPath*>*>* hrsp) {
	this->primaryPath = primaryPath;
	this->hrsp = hrsp;
	size_t n = primaryPath->nodeIds->size();

	this->waitingTimeoutOnNextPacketInUs = new vector<int>();
	this->waitingTimeoutOnNextPacketForLostReportInUs = new vector<int>();
	this->maxPIATInUsUnderFailover = new vector<int>();
	retransTimeouts = new vector<int>();
	waitingTimeoutOnRetransReportInUs = new unordered_map<int, int>*[n];

	for (size_t i = 0; i < n; ++i) {
		waitingTimeoutOnNextPacketInUs->push_back(0);
		waitingTimeoutOnNextPacketForLostReportInUs->push_back(0);
		maxPIATInUsUnderFailover->push_back(0);
		waitingTimeoutOnRetransReportInUs[i] = new unordered_map<int, int>();
		retransTimeouts->push_back(0);
	}
}

string* ResilientRoutes::toString() {
	stringstream ss;
	ss << "Resilient Routes Dump: \n" << *(primaryPath->toString()) + "\n";

	unsigned int n = primaryPath->nodeIds->size();
	for (unsigned int i = 0; i < n; ++i) {
		int nodeId = (*primaryPath->nodeIds)[i];

		if ((*hrsp).find(nodeId) != hrsp->end()) {
			ss << "--> " << nodeId << ": \n";
			vector<RedundantSubPath*>* rsps = (*hrsp)[nodeId];
			for (RedundantSubPath* rsp : *rsps) {
				ss << "\t" << *(rsp->toString()) << "\n";
			}
		}
	}
	/*
	for (unordered_map<int, vector<RedundantSubPath*>*>::iterator iter = hrsp->begin(); iter != hrsp->end(); ++iter) {
		ss << "--> " << iter->first << ": \n";
		vector<RedundantSubPath*>* rsps = iter->second;
		if (rsps == NULL) continue;
		for (RedundantSubPath* rsp : *rsps) {
			ss << "\t" << *(rsp->toString()) << "\n";
		}
	}*/

	ss << "RR loss rate: " << expectedRrLossRate << ", maxPathDelay: " << maxRrDelay << "\n";


	ss << "retransTimeouts: " << endl;
	for (unsigned int i = 0; i < n; ++i) {
		ss << to_string((*retransTimeouts)[i]) << (i < n - 1 ? ", " : "");
	}
	ss << endl;

	ss << "waitingTimeoutOnNextPacketInUs: " << endl;
	for (unsigned int i = 0; i < n; ++i) {
		ss << to_string((*waitingTimeoutOnNextPacketInUs)[i]) << (i < n - 1 ? ", " : "");
	}
	ss << endl;

	ss << "waitingTimeoutOnNextPacketForLostReportInUs: " << endl;
	for (unsigned int i = 0; i < n; ++i) {
		ss << to_string((*waitingTimeoutOnNextPacketForLostReportInUs)[i]) << (i < n - 1 ? ", " : "");
	}
	ss << endl;

	ss << "waitingTimeoutOnRetransReportInUs: " << endl;
	for (unsigned int i = 0; i < n; ++i) {
		unordered_map<int, int>* item = waitingTimeoutOnRetransReportInUs[i];
		ss << "--> " << (*primaryPath->nodeIds)[i] << ": ";
		for (unordered_map<int, int>::iterator iter = item->begin(); iter != item->end(); ++iter) {
			ss << "(" << iter->first << ", " << iter->second << "), ";
		}
		ss << endl;
	}
	ss << endl;

	string *s = new string(ss.str());
	return s;
}



LltcConfiguration::LltcConfiguration(int numRetransRequests, int beta, double pmuFreq, double maxDriftRatioForPmuFreq,
		int maxConsecutiveDriftPackets, double maxPathDelay, int numDataRetransReportsToSend) {
	this->numRetransRequests = numRetransRequests;
	this->beta = beta;
	this->pmuFreq = pmuFreq;
	this->maxDriftRatioForPmuFreq = maxDriftRatioForPmuFreq;
	this->maxConsecutiveDriftPackets = maxConsecutiveDriftPackets;
	this->maxPathDelay = maxPathDelay;
	this->numDataRetransReportsToSend = numDataRetransReportsToSend;
}

LltcResilientRouteVectors::LltcResilientRouteVectors(int srcNodeId, int n) {
	this->srcNodeId = srcNodeId;
	this->n = n;
	this->maxDelays = new double[n];
	this->expectedLossRates = new double[n];
	this->prevNodeIdx = new int[n];
	this->rspsMap = new vector<RedundantSubPath*>*[n];
	this->maxPrimPathDelays = new double[n];
}
