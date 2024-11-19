//
// Created by Yehong Xu on 15/2/2024.
//

#ifndef SWTD_TDGTREE_H
#define SWTD_TDGTREE_H

#endif //SWTD_TDGTREE_H

#include "graph.h"
#include "heap.h"

class TDGTree
{
public:

    int PARTITION_PART_NF = 4;
    int LEAF_TAU = 6;
    int acc;
    typedef struct
    {
        bool isborder;
        int treeid;
        vector<int> HGPtreepath;
    } HGPNode;

    typedef struct
    {
        bool isleaf;
        vector<int> leafnodes;
        int father;
        int level;
        int order;
        vector<int> borders;
        vector<int> children;
        vector<unordered_map<int, LPFunction>> dynMat;
        vector<unordered_map<int, LPFunction>> TDIF;
        set<int> X; // union of border nodes in its child nodes
    } HGPTreeNode;

    vector<HGPNode> HGPNodes;
    vector<HGPTreeNode> HGPTree;

    typedef struct
    {
        int rank;
        set<int> nset;
    } Status;

    Graph graph;

    explicit TDGTree(Graph &g, int acc);

    void ReadGTree(string &path);
    void dynMInit();
    void dynMInitTread(int begin, int end);

    void TDFloyd(HGPTreeNode &nd);
//    void TDFloyd(vector<int> &vertices, vector<unordered_map<int, LPFunction>> &dynMat);
    void TDFloydThread(vector<unordered_map<int, LPFunction>> &dynMat, vector<int> &X, int begin, int end);

    void BotUpLocalOpt();
    void BotUpLocalMain(HGPTreeNode &ni);
    void BotUpLocalOpt2();
    void BotUpLocalOptThread(vector<int> &needProcess, int begin, int end);

    void TopDGlobalOpt();
    void TopDGlobalOptMain(HGPTreeNode &ni);
    void TopDGlobalOpt2();
    void TopDGlobalOptThread(vector<int> &needProcess, int begin, int end);

    void TwoPMatrixBuild();
    int LCA(int u, int v);
    void Update(vector<pair<int, int>> &wBatch);
    void ArchiveDynM(string &path) const;
    void ReadDynM(string &path);
    LPFunction TIPQuery(int u, int v);
    int TDSPQuery(int u, int v, int t);
};
