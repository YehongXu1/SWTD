//
// Created by Yehong Xu on 2/2/2024.
//

#ifndef SWTD_DHL_H
#define SWTD_DHL_H

#include "TDGTree.h"

class DHL: public TDGTree
{

public:

    DHL(Graph &g, int acc);

    void TDIFGen();
    void TDIFGenThread(vector<int> &leafNodes, int begin, int end);

    Graph borderGraph;
    vector<int> borderNodes;

    void BorderGraphCon();

    void BuildDHLIndex();

    void updateBorderLabel(vector<pair<int, int>> &wBatch);
    LPFunction DHLQuery(int u, int v);

    int DHLQuery(int u, int v, int t);

    LPFunction disBetweenBorders(int id1, int id2);

    int disBetweenBorders(int id1, int id2, int t);
};

#endif //SWTD_DHL_H
