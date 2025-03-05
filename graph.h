#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <chrono>
#include "tools.h"
#include "LPFunction.h"
#include <list>
#include <stdlib.h>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <stack>
#include <boost/thread/thread.hpp>
#include <boost/functional/hash.hpp>
#include <boost/bind.hpp>
//#include "heap.h"

using namespace std;
using namespace benchmark;

static int threadNum = 1;


typedef struct COMPARENODE
{
    pair<int, int> pif;//ID, minCost
    bool operator()(const struct COMPARENODE &a, const struct COMPARENODE &b) const
    {
        return a.pif.second > b.pif.second;
    }
} compareNode;

struct pair_hash
{
    template<class T1, class T2>
    std::size_t operator()(const std::pair<T1, T2> &p) const
    {
        auto h1 = std::hash<T1>{}(p.first);
        auto h2 = std::hash<T2>{}(p.second);
        return h1 ^ h2;
    }
};

struct Edge
{
    int ID1, ID2, edgeID, length, twinID;
    vector<int> vX, vY; //Store f(t) in the window
    vector<int> vXFull, vYFull;
    LPFunction lpf;
    Edge()
    {
        ID1 = -1;
        ID2 = -1;
        edgeID = -1;
        length = -1;
        twinID = -1;
        vX.clear();
        vY.clear();
        vXFull.clear();
        vYFull.clear();
    }
    Edge(int id1, int id2, LPFunction &lpf)
    {
        ID1 = id1;
        ID2 = id2;
        this->lpf = lpf;
    }
};

struct NeighborInfo
{
    int nid;
    int w;
    int c;
};

static vector<int> GlobalNodeOrder;

struct TreeNode
{
    vector<pair<int, pair<LPFunction, int>>> vertIn, vertOut;
    vector<int> posIn, posOut;
    vector<int> cntIn, cntOut; //the label value and corresponding count number
    vector<LPFunction> disIn, disOut; //the distance value and corresponding count number
    set<int> changedPos;
    unordered_set<int> DisReIn, DisReOut;
    vector<int> vChildren;
    int height, hdepth;//hdepty is the deepest Node that a vertex still exists
    int uniqueVertex;
    int parentNodeID;
    vector<int> ancIDs;

    TreeNode()
    {
        vertIn.clear();
        vertOut.clear();
        posIn.clear();
        posOut.clear();
        disIn.clear();
        disOut.clear();
        vChildren.clear();
        parentNodeID = -1;
        uniqueVertex = -1;
        height = 0;
        hdepth = 0;
        changedPos.clear();
        DisReIn.clear();
        DisReOut.clear();
    }
};

struct subA
{
    int d;
    int a;
    int k1;
    int k2;
    bool bSE;
//	seEnum* se;
};

struct CompareA
{
    bool operator()(subA &a1, subA &a2)
    {
        if (a1.d < a2.d)
            return false;
        else if (a1.d == a2.d)
            return a1.a > a2.a;
        return true;
    }
};

struct labelEnum
{
    int aNum; //access Node number
    int ID1, ID2; //ID2 has higher order
    benchmark::heap<2, int, int> Q;
    vector<subA> vSubA;
    //priority_queue<subA, vector<subA>, CompareA> Q;//d, aID
    vector<int> vDistance; //d of a
    vector<int> vA; //pivot nodes
    unordered_map<int, int> umAPos;
    vector<int> vK; //Top-K Distance
    vector<int> vKA; //a that creates k`
    vector<pair<int, int> > vSubK; //For path retrieval
    vector<bool> vbFinished; //Finish Enumeration of A_i
    vector<pair<bool, bool> > vpbFinished;
    vector<pair<vector<int>, vector<int> > > vpvADistance;
    vector<pair<pair<int, int>, pair<int, int> > > vR;
    bool bFinished;

};

class Graph
{
public:
    Graph() {}

    Graph(string filenameGraph, string filenameMap, string filenameRoad, string filenameSpeed, int slotNum)
    {
        //if(filenameGraph == "./beijingNodeClean")
        if (filenameGraph == "../beijingNodeDirected")
        {
            readBeijingMapDirectedNew(filenameGraph, filenameMap);
            readBeijingTD(filenameRoad, filenameSpeed, slotNum);
        }
    }

    //True: Only adjList and adjListMap are used
    //False: Also use R versions
    int nodeNum;
    int edgeNum;
    int itvLen;

    vector<vector<pair<int, int> > > adjList;    //neighborID, Distance
    vector<vector<pair<int, int> > > adjListR;
    vector<vector<pair<int, int> > > adjEdge;    //neighborID, edgeID
    vector<vector<pair<int, int> > > adjEdgeR;

    vector<unordered_map<int, int> > vmEdge;
//	unordered_map<pair<int, int>, int, pair_hash> umEdge;

    vector<Edge> vEdge;
    vector<Edge> vEdgeR;

    vector<int> vNodeMap;  //From the old to New
    vector<int> vNodeMapR; //From New to the old

    vector<pair<double, double> > vCoor;
    unordered_map<string, pair<int, int> > mCoor;
    double minX, minY, maxX, maxY;

    void timeWinSwitch(int win1, int win2);

    void timeWinSwitchSimpleExtend(int win1, int win2);

    void timeWinSwitchSimpleExtendRange(int win1, int win2);

    int readBeijingMapDirectedNew(string filenameGraph, string filenameMap);

    int readBeijingTD(string filenameRoad, string filenameSpeed, int slotNum);

    void readDirectedGraph(string mapfile, string speedfile, int lowerB, int upperB, int slotNum);

    void readUndirectedGraph(string mapfile, string speedfile, string order, int slotnum);

    static void readODs(string filename, vector<pair<int, int>> &vOD, int nodeNum);

    static void readUpdates(string updateFile, int nodeNum, vector<pair<pair<int, int>, pair<int, double>>> &testdata);

    static void readUpdates2(string updateFile, int nodeNum, int timeSlot, double inc,
                             vector<pair<pair<int, int>, pair<int, double>>> &testdata);

    void readExtension(vector<pair<int, int>> &testdata, int appendPoint);

    void H2HDecBatDiGraph(vector<pair<pair<int, int>, pair<int, double>>> &wBatch, int timeWinId, int timeWinLen);

    void EachNodeDecMainThread(const vector<int> &startSources);

    void EachNodeDecMainRange(int begin, int end, const vector<int> &startSources);

    void EachNodeDecMainRec(int sourceR);

    void EachNodeDecMain(int sourceR);

    bool lpfIsSupportBy(
            LPFunction &f1, const LPFunction &halfLPF,
            const vector<int> &halfLPFChangedPoses, bool isSC) const;

    void CHUpdate(vector<pair<int, int>> &wBatch, int upd);

    LPFunction forwardSearch(int sourceID, int targetID, int &miny);

    void H2HIncBatDiGraph(vector<pair<int, int>> &wBatch, int updX);

    void EachNodeIncMainThread(const vector<int> &uRank);

    void EachNodeIncMain(int uRank);

    void EachNodeIncMainRec(int uRank);

    void H2HExtendBatDiGraph(vector<pair<int, int>> &wBatch);

    void EachNodeExtendMainThread(const vector<int> &uRank);

    void EachNodeExtendMain(int uRank);

    void EachNodeExtendMainRec(int uRank);

    //TDH2H
    int lBound, uBound, deltaT;
    vector<map<int, pair<LPFunction, int>>> E;
    vector<map<int, pair<LPFunction, int>>> ER;
    //the shortcut value and corresponding count number
//    vector<unordered_map<int,int>> Ecnt;

    vector<Semaphore *> vSm;
    vector<int> DD, DD2, DD3, DD4;
    vector<int> vNodeOrder;
    vector<int> NodeOrder;
    //vector<vector<pair<int,pair<int,int>>>> outLabel;
    vector<vector<pair<int, pair<LPFunction, int>>> > NeighborCon;
    vector<vector<pair<int, pair<LPFunction, int>>> > NeighborConR;
    //<nodeID1, nodeID2>, <support nodes>
    //nodeID1 < nodeID2
    //Supportive nodes are incoperated in the LPF now
    vector<map<int, vector<int>>> SCconNodesMT;

    // store every LPF created during node contraction or updates, nodeID1 != nodeID2 != nodeID3
    vector<unordered_map<int, unordered_map<int, CatSupRec>>> intermediateSCs;
    // store every LPF created during label construction, nodeID1 != nodeID2
    vector<unordered_map<int, unordered_map<int, CatSupRec>>> intermediateLBsIn;
    vector<unordered_map<int, unordered_map<int, CatSupRec>>> intermediateLBsOut;

    void loadTD();

    void TDConstruction();

    void detNodeOrder();

    void CHConstruction();

    void TDContractionThread(
            int x, pair<int, int> pID, vector<pair<int, pair<LPFunction, int>>> &Neigh,
            vector<pair<int, pair<LPFunction, int>>> &NeighR);

    void makeTree();

    void makeIndex();

    void makeIndexDFS(int p, vector<int> &ancs);

    void makeIndexDFSThread();

    void makeIndexDFSRec(int p, vector<int> &ancs);

    void extendIndex();

    void lastChangedPosInit(int begin, int end);

    //int deepestAnc(int x, vector<pair<int,pair<int,int>>> &vert);
    int deepestAnc(vector<pair<int, pair<LPFunction, int>>> &vert);

    void insertE(int u, int v, LPFunction &lpf);

    void deleteE(int u, int v);

    vector<int> rank;
    vector<TreeNode> Tree;
    int heightMax;
    vector<vector<int>> VidtoTNid;//one vertex exist in those tree nodes (nodeID--->tree Node rank)

    vector<int> EulerSeq;
    vector<int> toRMQ;
    vector<vector<int>> RMQIndex;

    void makeRMQ();

    void makeRMQDFS(int p, int height);

    LPFunction QueryH2H(int ID1, int ID2);

    int LCAQuery(int _p, int _q);

//	vector<vector<seEnum> > vvSE;
//	vector<vector<seEnum> > vvSEQuery;//Copy from vvSE for each query 
//    vector<map<int, pair<LPFunction, int>> > vmSE, vmSER;
//    vector<unordered_map<int, int> > vmSEPos;
//	vector<unordered_map<int, bool> > vmEFinished;//SE+EInBG
//	vector<unordered_map<int, bool> > vmEFinishedQuery;//SE+EInBG
    //bool seCompare(seEnum se1, seEnum se2);
//    vector<vector<labelEnum> > vvLE;
//    vector<vector<labelEnum> > vvLEQuery;
//    vector<unordered_map<int, int> > vmLEPos;
//    vector<unordered_map<int, bool> > vmLFinished;//SE+EInBG
//    vector<unordered_map<int, bool> > vmLFinishedQuery;//SE+EInBG

    LPFunction forwardSearch(int sourceID, int targetID);

    void decUpdate(vector<pair<pair<int, int>, pair<int, double>>> &testdata, int timeWId, int timeWLen);

    void incUpdate(vector<pair<pair<int, int>, pair<int, double>>> &testdata, int timeWId, int timeWLen);

    void indexSize();

    void indexSizeThread(int begin, int end, int &scSize, int &lbSize);

    int Dijkstra(int ID1, int ID2);

    void Dijkstra(int ID1, int ID2, int t, int &result);

    int DijkstraPath(int ID1, int ID2, vector<int> &vPath, vector<int> &vPathEdge);

    int DijkstraPath2(int ID1, int ID2, unordered_set<int> &sRemovedNode, vector<int> &vPath, vector<int> &vPathEdge);

    int AStar(int ID1, int ID2);

    LPFunction QueryH2H(int ID1, int ID2,  long long int &LCASize);

    int QueryH2HFixedT(int ID1, int ID2, int t,  long long int &LCASize);

    int QueryCHFixedT(int ID1, int ID2, int t);

    LPFunction QueryCHItvT(int ID1, int ID2, int &t);

    int AStarPath(int ID1, int ID2, vector<int> &vPath, vector<int> &vPathEdge, string &city);

    int EuclideanDistance(int ID1, int ID2);

    int EuclideanDistanceAdaptive(int ID1, int ID2, int latU, int lonU);

    // precompute
    vector<TreeNode> precomputeSCForA();

    void makeIndexForA(vector<TreeNode> &TreeForA);

    void makeIndexDFSThreadForA(vector<TreeNode> &TreeForA);

    void makeIndexDFSRecForA(vector<TreeNode> &TreeForA, int p, vector<int> &ancs);

    void makeIndexDFSForA(vector<TreeNode> &TreeForA, int p, vector<int> &ancs);
};


#endif
