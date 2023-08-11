# Heuristic MPSG Generation Algorithm (HeuMPSG) for Disruption Resilient Transport Protocol (DRTP) 

## HeuMPSG Introduction
DRTP maintains a recursive hop-by-hop retransmission process using the subpaths-based routes resiliently structured as the multi-path subgraph (MPSG), where MPSG consists of a reliable primary path (PP) and its redundant subpaths (RSPs) to connect two hops in PP as much as possible. However, for the recursiveness, it is non-trivial to efficiently generate the resilient MPSGs with high computation complexity being NP-complete. To exploit the resilience, HeuMPSG generate MPSGs less fragile, while for each MPSG, its recursive process meets a constraint on the end-to-end packet delivery time (EEDT) in maximum. It heuristically finds a candidate PP hop-by-hop, and for each hop, searches maximally disjoint RSPs, with a lower computation cost. For real power grids, HueMPSG can significantly improve the resilience with more RSPs averaged at most 171.51%, at 7.07 times faster speed. Under heavy random faults, HueMPSG significantly improves the resilience with lower EEFRs and EEDTs at most in 29.11% and 2.41 s, respectively. The source code of HeuMPSG can be found in:
- ./ndnSIM/ns-3/src/ndnSIM/model/lltc/lltc-resilient-routes-generation-for-drtp.cpp
- ./ndnSIM/ns-3/src/ndnSIM/model/lltc/lltc-resilient-routes-generation-for-drtp.hpp

## How to Run the Source Code of HeuMPSG?
1. Build DRTP based on ndnSIM according to the instrutions available in https://github.com/zhouby-zjl/drtp/README.md
2. Copy the source code of HeuMPSG to the ndnSIM directory of DRTP
3. To test HeuMPSG with DRTP, you can run:  ./waf --run scratch/drtp-sim-static-optmpsg --command-template="%s drtp-config-heumpsg.ini". To only test HeuMPSG, you can run: ./waf --run scratch/drtp-sim-static-optmpsg-cp --command-template="%s drtp-config-heumpsg-cp.ini". Wherein, the ini files are available in the root directory of HeuMPSG that can be flexibly redefined for test purposes. 

## HeuMPSG's Ongoing Papers
Boyang Zhou, Tong Ye, Wenjie Yu, and Chunming Wu. Enhancing Resilience of Phasor Data Acquisition in Industrial Networks of Power Transmission Grids. IEEE Globecom Workshop on SRINetworks. 2023 (under review). 

## Links to Other Related Source Code
- DRTP: https://github.com/zhouby-zjl/drtp/
- ndnSIM: https://github.com/named-data-ndnSIM/

