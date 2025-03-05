//
// Created by Yehong Xu on 2/2/2024.
//

#include "DHL.h"

DHL::DHL(Graph &g) : TDGTree(g)
{
}

void DHL::BuildDHLIndex()
{
    cout << "Start building DHL index" << endl;

    Timer timer;
    timer.tick();
    TwoPMatrixBuild();
    timer.tock();
    cout << "TwoPMatrixBuild time: " << timer.duration().count() << endl;

    timer.tick();
    TDIFGen();
    timer.tock();
    cout << "TDIFGen time: " << timer.duration().count() << endl;

    timer.tick();
    BorderGraphCon();
    timer.tock();
    cout << "BorderGraphCon time: " << timer.duration().count() << endl;

    timer.tick();
    borderGraph.CHConstruction();
    timer.tock();
    cout << "BorderLabelCon time: " << timer.duration().count() << endl;

    borderGraph.makeTree();
}

void DHL::TDIFGen()
{
    vector<int> leafNodes;
    for (int i = 0; i < HGPTree.size(); i++)
    {
        if (HGPTree[i].isleaf)
        {
            leafNodes.push_back(i);
        }
    }

    int thisThreadNum = threadNum;
    int step = leafNodes.size() / thisThreadNum;
    vector<pair<int, int>> vpID;
    for (int i = 0; i < thisThreadNum; i++)
    {
        int begin = i * step;
        int end = (i + 1) * step;
        if (i == thisThreadNum - 1)
        {
            end = leafNodes.size();
        }
        vpID.emplace_back(begin, end);
    }

    boost::thread_group threads;
    for (int i = 0; i < thisThreadNum; i++)
    {
        threads.add_thread(
                new boost::thread(
                        &DHL::TDIFGenThread, this, leafNodes, vpID[i].first, vpID[i].second
                ));
    }
    threads.join_all();
}

void DHL::TDIFGenThread(vector<int> &leafNodes, int begin, int end)
{
    for (int k = begin; k < end; k++)
    {
        HGPTreeNode &nd = HGPTree[leafNodes[k]];

        nd.TDIF.assign(graph.nodeNum, unordered_map < int, LPFunction > ());
        for (int i = 0; i < nd.borders.size(); i++)
        {
            for (int j = 0; j < nd.borders.size(); j++)
            {
                if (i == j)
                    continue;

                int id1 = nd.borders[i];
                int id2 = nd.borders[j];

                LPFunction lpfMin(id1, id2, graph.lBound, graph.uBound);
                for (const auto &vToW: graph.adjEdge[id1])
                {
                    if (vToW.first == id2)
                    {
                        lpfMin = graph.vEdge[vToW.second].lpf;
                    }
                    if (!HGPNodes[vToW.first].isborder)
                    {
                        for (const auto &vToW2: graph.adjEdge[vToW.first])
                        {
                            if (vToW2.first == id2)
                            {
                                LPFunction lpf1 = graph.vEdge[vToW.second].lpf;
                                LPFunction lpf2 = graph.vEdge[vToW2.second].lpf;
                                if (lpfMin.vX.size() < 2)
                                {
                                    lpfMin = lpf1.LPFCatSupport(lpf2, graph.lBound, graph.uBound);
                                } else
                                {
                                    LPFunction lpfTmp = lpf1.LPFCatSupport(lpf2, graph.lBound, graph.uBound);
                                    if (lpfMin.vX.size() > 1 and lpf1.minY + lpf2.minY <= lpfMin.minY and
                                        !lpfMin.dominate(lpfTmp))
                                    {
                                        lpfMin = lpfMin.LPFMinSupForDec(lpfTmp);
                                    }
                                }
                            }
                        }
                    }
                }

                if (lpfMin.lowerBound != graph.lBound or lpfMin.upperBound != graph.uBound)
                {
                    assert(false);
                }
                nd.TDIF[id1][id2] = lpfMin;
            }
        }
    }
}

void DHL::BorderGraphCon()
{
    borderGraph.nodeNum = graph.nodeNum;
    borderGraph.edgeNum = graph.edgeNum;
    borderGraph.lBound = graph.lBound;
    borderGraph.uBound = graph.uBound;
    borderGraph.itvLen = graph.itvLen;
    borderGraph.adjEdge.assign(graph.nodeNum, vector<pair<int, int>>());
    borderGraph.adjEdgeR.assign(graph.nodeNum, vector<pair<int, int>>());
    borderGraph.vEdge.reserve(graph.edgeNum);
    borderGraph.intermediateSCs.clear();
    borderGraph.intermediateSCs.assign(
            borderGraph.nodeNum, unordered_map < int, unordered_map < int, CatSupRec >> ());
    borderGraph.intermediateLBsIn.clear();
    borderGraph.intermediateLBsIn.assign(
            borderGraph.nodeNum, unordered_map < int, unordered_map < int, CatSupRec >> ());
    borderGraph.intermediateLBsOut.clear();
    borderGraph.intermediateLBsOut.assign(
            borderGraph.nodeNum, unordered_map < int, unordered_map < int, CatSupRec >> ());

    for (const auto &node: HGPTree)
    {
        if (node.isleaf and node.borders.size() > 1)
        {
            for (int j = 0; j < node.borders.size(); j++)
            {
                for (int k = 0; k < node.borders.size(); k++)
                {
                    if (j == k)
                        continue;
                    int id1 = node.borders[j];
                    int id2 = node.borders[k];

                    if (node.TDIF[id1].find(id2) != node.TDIF[id1].end())
                    {
                        LPFunction lpf = node.TDIF[id1].find(id2)->second;
                        if (lpf.vX.size() > 1)
                        {
                            if (lpf.lowerBound != graph.lBound or lpf.upperBound != graph.uBound)
                            {
                                lpf.display();
                                assert(false);
                            }
                            borderGraph.adjEdge[id1].emplace_back(id2, (int) borderGraph.vEdge.size());
                            borderGraph.adjEdgeR[id2].emplace_back(id1, (int) borderGraph.vEdge.size());
                            Edge e(id1, id2, lpf);
                            borderGraph.vEdge.push_back(e);
                        }
                    }
                }
            }
        }

    }

    for (int i = 0; i < HGPNodes.size(); i++)
    {
        if (HGPNodes[i].isborder)
        {
            borderNodes.push_back(i);
        }
    }

    // Minimum degree first
    borderGraph.detNodeOrder();

    // Maximum degree first
    vector<int> vNodeOrder = borderGraph.vNodeOrder;
    borderGraph.vNodeOrder.clear();
    vector<int> NodeOrderInBG = borderGraph.NodeOrder;
    borderGraph.NodeOrder.clear();
    borderGraph.NodeOrder.assign(graph.nodeNum, -1);
    for (int i = vNodeOrder.size() - 1; i >= 0; i--)
    {
        borderGraph.vNodeOrder.push_back(vNodeOrder[i]);
    }

    for (int k = 0; k < (int) borderGraph.vNodeOrder.size(); k++)
    {
        borderGraph.NodeOrder[borderGraph.vNodeOrder[k]] = k;
//        cout << vNodeOrder[k] << ": " << k << endl;
    }
}

int DHL::disBetweenBorders(int id1, int id2, int t)
{
//    LPFunction lpf(id1, id2, graph.lBound, graph.uBound);
    int cost = INF;
    if (borderGraph.NodeOrder[id1] < borderGraph.NodeOrder[id2] and
        id1 < borderGraph.rank.size() and borderGraph.rank[id1] >= 0)
    {
        for (const auto &p: borderGraph.Tree[borderGraph.rank[id1]].vertOut)
        {
            int cost2 = -1;
            if (p.first == id2)
            {
                cost2 = p.second.first.getY(t);
            } else
            {
                if (p.first >= borderGraph.rank.size() or borderGraph.rank[p.first] < 0)
                    continue;

                for (const auto &p2: borderGraph.Tree[borderGraph.rank[p.first]].vertOut)
                {
                    if (p2.first == id2)
                    {
                        LPFunction lpf2 = p.second.first, lpf3 = p2.second.first;
                        int cost1 = lpf2.getY(t);
                        if (t + cost1 < graph.uBound)
                        {
                            cost2 = lpf3.getY(t + cost1);
                        }
                    }
                }
            }

            if (cost2 > -1 and cost2 < cost)
                cost = cost2;
        }
    } else if (id2 < borderGraph.rank.size() and borderGraph.rank[id2] >= 0)
    {
        for (const auto &p: borderGraph.Tree[borderGraph.rank[id2]].vertIn)
        {
            int cost2 = -1;
            if (p.first == id1)
            {
                cost2 = p.second.first.getY(t);
            } else
            {
                if (p.first >= borderGraph.rank.size() or borderGraph.rank[p.first] < 0)
                    continue;
                for (const auto &p2: borderGraph.Tree[borderGraph.rank[p.first]].vertOut)
                {
                    if (p2.first == id1)
                    {
                        LPFunction lpf2 = p.second.first, lpf3 = p2.second.first;
                        int cost1 = lpf2.getY(t);
                        if (cost1 + t < graph.uBound)
                            cost2 = lpf3.getY(cost1 + t);
                    }
                }
            }

            if (cost2 > 0 and cost2 < cost)
                cost = cost2;
        }
    }
    return cost;
}

LPFunction DHL::disBetweenBorders(int id1, int id2)
{
    LPFunction lpf(id1, id2, graph.lBound, graph.uBound);
    if (borderGraph.NodeOrder[id1] < borderGraph.NodeOrder[id2] and id1 < borderGraph.rank.size() and
        borderGraph.rank[id1] >= 0)
    {
        for (const auto &p: borderGraph.Tree[borderGraph.rank[id1]].vertOut)
        {
            LPFunction lpfTmp;
            if (p.first == id2)
            {
                lpfTmp = p.second.first;
            } else
            {
                if (p.first >= borderGraph.rank.size() or borderGraph.rank[p.first] < 0)
                    continue;
                for (const auto &p2: borderGraph.Tree[borderGraph.rank[p.first]].vertOut)
                {
                    if (p2.first == id2)
                    {
                        LPFunction lpf2 = p.second.first, lpf3 = p2.second.first;
                        if (lpf2.ID2 == lpf3.ID1)
                            lpfTmp = lpf2.LPFCatSupport(lpf3, graph.lBound, graph.uBound);
                    }
                }
            }

            if (lpf.vX.size() < 2)
            {
                lpf = lpfTmp;
            } else if (lpfTmp.minY < lpf.maxY and !lpf.dominate(lpfTmp))
            {
                lpf = lpf.LPFMinSupForDec(lpfTmp);
            }
        }
    } else if (id2 < borderGraph.rank.size() and borderGraph.rank[id2] >= 0)
    {
        for (const auto &p: borderGraph.Tree[borderGraph.rank[id2]].vertIn)
        {
            LPFunction lpfTmp;
            if (p.first == id1)
            {
                lpfTmp = p.second.first;
            } else
            {
                if (p.first >= borderGraph.rank.size() or borderGraph.rank[p.first] < 0)
                    continue;
                for (const auto &p2: borderGraph.Tree[borderGraph.rank[p.first]].vertOut)
                {
                    if (p2.first == id1)
                    {
                        LPFunction lpf2 = p.second.first, lpf3 = p2.second.first;
                        if (lpf2.ID2 == lpf3.ID1)
                            lpfTmp = lpf2.LPFCatSupport(lpf3, graph.lBound, graph.uBound);
                    }
                }
            }

            if (lpf.vX.size() < 2)
            {
                lpf = lpfTmp;
            } else if (lpfTmp.minY < lpf.maxY and !lpf.dominate(lpfTmp))
            {
                lpf = lpf.LPFMinSupForDec(lpfTmp);
            }
        }
    }
    return lpf;
}

LPFunction DHL::DHLQuery(int u, int v)
{
    LPFunction lpf(u, v, graph.lBound, graph.uBound);
    int LCAId = LCA(u, v);
    if (HGPTree[LCAId].isleaf)
    {
        lpf = HGPTree[LCAId].dynMat[u][v];
    } else
    {
        for (const auto &bid: HGPTree[LCAId].borders)
        {
            for (const auto &bid2: HGPTree[LCAId].borders)
            {
                int test;
                LPFunction sb1 = graph.forwardSearch(u, bid, test);
                LPFunction b1b2 = disBetweenBorders(bid, bid2);
                LPFunction b2t = graph.forwardSearch(bid2, v, test);
                if (sb1.vX.size() < 2 or b1b2.vX.size() < 2 or b2t.vX.size() < 2)
                    continue;
                if (sb1.ID2 != b1b2.ID1 or b1b2.ID2 != b2t.ID1)
                    continue;

                if (lpf.vX.size() <= 1)
                {
                    LPFunction tmp = sb1.LPFCatSupport(b1b2, graph.lBound, graph.uBound);
                    lpf = tmp.LPFCatSupport(b2t, graph.lBound, graph.uBound);
                } else if (sb1.minY + b1b2.minY + b2t.minY < lpf.maxY)
                {
                    LPFunction tmp = sb1.LPFCatSupport(b1b2, graph.lBound, graph.uBound);
                    tmp = tmp.LPFCatSupport(b2t, graph.lBound, graph.uBound);
                    if (!lpf.dominate(tmp) and lpf.ID1 == tmp.ID1 and lpf.ID2 == tmp.ID2)
                        lpf = lpf.LPFMinSupForDec(tmp);
                }
            }
        }
    }
    return lpf;
}

int DHL::DHLQuery(int u, int v, int t)
{
    int cost = INF;
    int LCAId = LCA(u, v);
    if (HGPTree[LCAId].isleaf)
    {
        cost = HGPTree[LCAId].dynMat[u][v].getY(t);
    } else
    {
        for (const auto &bid: HGPTree[LCAId].borders)
        {
            for (const auto &bid2: HGPTree[LCAId].borders)
            {
                int sb1;
                graph.Dijkstra(u, bid, t, sb1);
                int b1b2 = disBetweenBorders(bid, bid2, t + sb1);
                int b2t;
                graph.Dijkstra(bid2, v, t + sb1 + b1b2, b2t);
                if (sb1 + b1b2 + b2t < cost)
                    cost = sb1 + b1b2 + b2t;
            }
        }
    }
    return cost;
}

void DHL::updateBorderLabel(vector<pair<int, int>> &wBatch, int upd)
{
    borderGraph.E.clear();
    borderGraph.ER.clear();
    for (int i = 0; i < borderGraph.nodeNum; i++)
    {
        borderGraph.NeighborCon[i].clear();
        borderGraph.NeighborConR[i].clear();
        borderGraph.Tree[i].vertIn.clear();
        borderGraph.Tree[i].vertOut.clear();
    }
    borderGraph.CHConstruction();
    Update(wBatch, upd);
    TDIFGen();
}




