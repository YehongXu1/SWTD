//
// Created by Yehong Xu on 11/2/2024.
//

#include "DHL.h"

TDGTree::TDGTree(Graph &g, int accInput) : graph(g)
{
    acc = accInput;
}

void TDGTree::ReadGTree(string &path)
{
    ifstream in(path.c_str());
    if (!in)
    {
        cout << "Cannot open file " << path << endl;
        return;
    }

    HGPTree.reserve(graph.nodeNum);
    HGPNode hgpnode;
    hgpnode.isborder = false;
    HGPNodes.assign(graph.nodeNum, hgpnode);

    string line;
    while (getline(in, line))
    {
        vector<string> vs = Tools::split(line, " ");
        vector<int> info;
        for (auto &v: vs)
        {
            info.emplace_back(atoi(v.c_str()));
        }

        int tNId = info[0];
        int tFId = info[1];
        bool isLeaf = info[2];
        int borN = info[3];
        vector<int> borders = vector<int>(info.begin() + 4, info.begin() + 4 + borN);
        int leafN = info[4 + borN];
        vector<int> leafnodes = vector<int>(info.begin() + 5 + borN, info.begin() + 5 + borN + leafN);
        HGPTreeNode treeNode;
        treeNode.order = tNId;
        treeNode.father = tFId;
        treeNode.isleaf = isLeaf;
        treeNode.borders = borders;
        for (const auto &bid: borders)
        {
            HGPNodes[bid].isborder = true;
        }

        treeNode.leafnodes = leafnodes;
        if (tFId >= 0)
        {
            HGPTree[treeNode.father].children.emplace_back(tNId);
            assert(HGPTree[treeNode.father].order < tNId);
            treeNode.level = HGPTree[tFId].level + 1;
        } else
        {
            treeNode.level = 0;
        }
        HGPTree.emplace_back(treeNode);
    }
    in.close();

    unordered_map<int, int> heightCnt;
    for (int tNid = 0; tNid < HGPTree.size(); tNid++)
    {
        if (HGPTree[tNid].isleaf)
        {
            for (const auto &nid: HGPTree[tNid].leafnodes)
            {
                HGPNode &nd = HGPNodes[nid];
                nd.treeid = tNid;
                vector<int> treePath;
                treePath.emplace_back(tNid);
                int f = HGPTree[tNid].father;
                while (f >= 0)
                {
                    treePath.emplace_back(f);
                    f = HGPTree[f].father;
                }

                nd.HGPtreepath.reserve(treePath.size());
                for (int i = treePath.size() - 1; i >= 0; i--)
                {
                    nd.HGPtreepath.emplace_back(treePath[i]);
                }
            }
        }
    }

    for (int i = 0; i < HGPTree.size(); i++)
    {
        Semaphore *sm = new Semaphore(1);
        graph.vSm.push_back(std::move(sm));
    }

//    for (const auto &p: heightCnt)
//    {
//        cout << p.first << " " << p.second << endl;
//    }
    cout << "Finish building GTree" << endl;
}

void TDGTree::dynMInit()
{
    int thisThreadNum = threadNum;
    int step = (int) HGPTree.size() / thisThreadNum;
    vector<pair<int, int>> intervals;
    for (int i = 0; i < thisThreadNum; i++)
    {
        int begin = i * step;
        int end = (i + 1) * step;
        if (i == thisThreadNum - 1)
            end = (int) HGPTree.size();
        intervals.emplace_back(make_pair(begin, end));
    }

    boost::thread_group threads;
    for (int i = 0; i < thisThreadNum; i++)
    {
        threads.create_thread(boost::bind(
                &TDGTree::dynMInitTread, this, intervals[i].first, intervals[i].second));
    }
    threads.join_all();
//    cout << "Finish init dynMat" << endl;
}

void TDGTree::dynMInitTread(int begin, int end)
{
    for (int i = begin; i < end; i++)
    {
        HGPTreeNode &nd = HGPTree[i];
        nd.dynMat.assign(graph.nodeNum, unordered_map < int, LPFunction > ());
        if (!nd.isleaf)
        {
            for (const auto &cid: nd.children)
            {
                for (const auto &bid: HGPTree[cid].borders)
                {
                    nd.X.insert(bid);
                }
            }
        } else
        {
            nd.X = set < int > (nd.leafnodes.begin(), nd.leafnodes.end());
        }

        for (const auto &bid: nd.X)
        {
            for (const auto &bid2: nd.X)
            {
                if (bid == bid2)
                    continue;
                nd.dynMat[bid][bid2] = LPFunction(bid, bid2, graph.lBound, graph.uBound);
                for (const auto &p: graph.adjEdge[bid])
                {
                    if (p.first == bid2)
                    {
                        nd.dynMat[bid][bid2] = graph.vEdge[p.second].lpf;
                        break;
                    }
                }
            }
        }
    }
}

void TDGTree::TDFloyd(HGPTreeNode &nd)
{
    if (nd.X.empty())
        return;
    vector<int> X(nd.X.begin(), nd.X.end());

//    cout << "\t" << nd.order << " " << X.size() << endl;
    for (int i = 0; i < X.size(); i++)
    {
        int k = X[i];
        for (const auto &s: X)
        {
            if (s == k)
                continue;
            for (const auto &e: X)
            {
                if (s == e || e == k || nd.dynMat[s][k].vX.size() < 2 || nd.dynMat[k][e].vX.size() < 2)
                    continue;

                LPFunction lpfTmp = nd.dynMat[s][k].LPFCatSupport(nd.dynMat[k][e], graph.lBound, graph.uBound, acc);
//                if (nd.dynMat[s][k].minY < 0 or nd.dynMat[k][e].minY > 7200)
//                {
////                    nd.dynMat[s][e].display();
//                    nd.dynMat[s][k].display();
//                    nd.dynMat[k][e].display();
//                    assert(false);
//                }
                if (nd.dynMat[s][e].vX.size() < 2)
                {
                    nd.dynMat[s][e] = lpfTmp;
                } else if (nd.dynMat[s][k].minY + nd.dynMat[k][e].minY < nd.dynMat[s][e].maxY and
                           !nd.dynMat[s][e].dominate(lpfTmp, acc))
                    nd.dynMat[s][e] = nd.dynMat[s][e].LPFMinSupForDec(lpfTmp, acc);
            }
        }

    }

}

void TDGTree::BotUpLocalOpt()
{
    int k = (int) HGPTree.size();
    cout << "TDGTree Bottom Up Construction" << endl;

    for (int i = k - 1; i >= 0; i--)
    {
        if (i % 1000 == 0)
            cout << "\t" << i << " " << HGPTree[i].X.size() << endl;

        HGPTreeNode &ni = HGPTree[i];
        BotUpLocalMain(ni);
    }
    cout << "Finish bottom up local opt" << endl;
}

void TDGTree::BotUpLocalMain(HGPTreeNode &ni)
{
    if (ni.X.empty())
        return;
    TDFloyd(ni);

    if (ni.father < 0)
        return;

    HGPTreeNode &nf = HGPTree[ni.father];
    if (nf.X.empty())
        return;
    vector<int> uB;
    for (const auto &bid: nf.X)
    {
        if (ni.X.find(bid) != ni.X.end())
        {
            uB.push_back(bid);
        }
    }

    for (const auto &bid1: uB)
    {
        for (const auto &bid2: uB)
        {
            if (bid1 == bid2)
                continue;

            assert(nf.dynMat.size() > bid1 and ni.dynMat.size() > bid1);
            assert(nf.dynMat[bid1].find(bid2) != nf.dynMat[bid1].end() and
                   ni.dynMat[bid1].find(bid2) != ni.dynMat[bid1].end());
            graph.vSm[nf.order]->wait();
            nf.dynMat[bid1][bid2] = ni.dynMat[bid1][bid2];
            graph.vSm[nf.order]->notify();
        }
    }
}

void TDGTree::BotUpLocalOpt2()
{
    vector<int> needProcess;
    int height = 0;
    for (const auto &node: HGPTree)
    {
        if (node.isleaf)
        {
            height = height < node.level ? node.level : height;
            needProcess.emplace_back(node.order);
        }
    }

    int thisThreadNum = threadNum;
    int depth = 0;
    while (!needProcess.empty())
    {
        cout << "depth: " << height - depth << " needProcess: " << needProcess.size() << endl;
        thisThreadNum = thisThreadNum < needProcess.size() ? thisThreadNum : needProcess.size();
        int step = (int) needProcess.size() / thisThreadNum;
        vector<pair<int, int>> intervals;
        for (int i = 0; i < thisThreadNum; i++)
        {
            int begin = i * step;
            int end = (i + 1) * step;
            if (i == thisThreadNum - 1)
                end = (int) needProcess.size();
            intervals.emplace_back(make_pair(begin, end));
        }

        boost::thread_group threads;
        for (int i = 0; i < thisThreadNum; i++)
        {
            threads.create_thread(boost::bind(
                    &TDGTree::BotUpLocalOptThread, this,
                    boost::ref(needProcess), intervals[i].first, intervals[i].second));
        }
        threads.join_all();

        set<int> needProcess2;
        for (const auto &nid: needProcess)
        {
            if (HGPTree[nid].father >= 0)
                needProcess2.insert(HGPTree[nid].father);
        }
        needProcess.clear();
        needProcess.assign(needProcess2.begin(), needProcess2.end());
        depth += 1;
    }
}

void TDGTree::BotUpLocalOptThread(vector<int> &needProcess, int begin, int end)
{
    for (int i = begin; i < end; i++)
    {
        HGPTreeNode &ni = HGPTree[needProcess[i]];
        BotUpLocalMain(ni);
    }
}

void TDGTree::TopDGlobalOpt()
{
    int i = 0;
    cout << "TDGTree Top Down Construction" << endl;
    for (auto &ni: HGPTree)
    {
        if (i % 1000 == 0)
            cout << "\t" << i << " " << HGPTree[i].X.size() << endl;
        i++;

        TopDGlobalOptMain(ni);
    }
    cout << "Finish top down global opt" << endl;
}

void TDGTree::TopDGlobalOpt2()
{
    vector<int> needProcess;
    for (const auto &ni: HGPTree[0].children)
    {
        needProcess.push_back(ni);
    }

    int thisThreadNum = threadNum;
    int depth = 0;
    while (!needProcess.empty())
    {
        cout << "depth: " << depth << " needProcess: " << needProcess.size() << endl;

        thisThreadNum = thisThreadNum < needProcess.size() ? thisThreadNum : needProcess.size();
        int step = (int) needProcess.size() / thisThreadNum;
        vector<pair<int, int>> intervals;
        for (int i = 0; i < thisThreadNum; i++)
        {
            int begin = i * step;
            int end = (i + 1) * step;
            if (i == thisThreadNum - 1)
                end = (int) needProcess.size();
            intervals.emplace_back(make_pair(begin, end));
        }

        boost::thread_group threads;
        for (int i = 0; i < thisThreadNum; i++)
        {
            threads.create_thread(boost::bind(
                    &TDGTree::TopDGlobalOptThread, this,
                    boost::ref(needProcess), intervals[i].first, intervals[i].second));
        }
        threads.join_all();

        set<int> needProcess2;
        for (const auto &nid: needProcess)
        {
            for (const auto &cid: HGPTree[nid].children)
            {
                needProcess2.insert(cid);
            }
        }
        needProcess.clear();
        needProcess.assign(needProcess2.begin(), needProcess2.end());

        depth += 1;
    }
}

void TDGTree::TopDGlobalOptThread(vector<int> &needProcess, int begin, int end)
{
    for (int i = begin; i < end; i++)
    {
        HGPTreeNode &ni = HGPTree[needProcess[i]];
        TopDGlobalOptMain(ni);
    }
}

void TDGTree::TopDGlobalOptMain(HGPTreeNode &ni)
{
    if (ni.isleaf)
        return;

    for (const auto &gcId: ni.children)
    {
        HGPTreeNode &nc = HGPTree[gcId];
        vector<int> uBorders;
        for (const auto &bid: nc.X)
        {
            if (ni.X.find(bid) != ni.X.end())
            {
                uBorders.push_back(bid);
            }
        }

        bool dirty = false;
        for (const auto &bid1: uBorders)
        {
            for (const auto &bid2: uBorders)
            {
                if (bid1 == bid2) continue;

                if (nc.dynMat[bid1].find(bid2) == nc.dynMat[bid1].end() ||
                    nc.dynMat[bid1][bid2].vX.size() < 2 ||
                    !ni.dynMat[bid1][bid2].equal(nc.dynMat[bid1][bid2], acc))
                {
                    nc.dynMat[bid1][bid2] = ni.dynMat[bid1][bid2];
                    dirty = true;
                }
            }
        }

        if (dirty)
            TDFloyd(nc);
    }
}

void TDGTree::TwoPMatrixBuild()
{
    Timer timer;

    timer.tick();
    dynMInit();
    timer.tock();
    cout << "dynMInit time: " << timer.duration().count() << endl;

    timer.tick();
    BotUpLocalOpt2();
    timer.tock();
    cout << "BotUpLocalOpt time: " << timer.duration().count() << endl;

    timer.tick();
    TopDGlobalOpt2();
    timer.tock();
    cout << "TopDGlobalOpt time: " << timer.duration().count() << endl;
}

void TDGTree::Update(vector<pair<int, int>> &wBatch)
{
    dynMInit();
    benchmark::heap<2, int, int> UQ(graph.nodeNum);
    benchmark::heap<2, int, int> DQ(graph.nodeNum);
    for (const auto &p: wBatch)
    {
        int u = p.first;
        int v = p.second;
        int lId = graph.NodeOrder[u] < graph.NodeOrder[v] ? u : v;
        std::random_device rd; // obtain a random number from hardware
        std::mt19937 gen(rd()); // seed the generator
        std::uniform_int_distribution<> distr(1, 10);

        int vX, oldW, newW;
        for (auto &pp: graph.adjEdge[u])
        {
            if (pp.first == v)
            {
                LPFunction &lpf = graph.vEdge[pp.second].lpf;
                if (distr(gen) <= 7)
                {
                    vX = lpf.vX.back();
                    oldW = lpf.vY.back();
                } else
                {
                    vX = lpf.vX.back() - graph.deltaT;
                    oldW = lpf.getY(vX);
                }

                if (distr(gen) <= 5)
                    newW = ceil((double) oldW * 1.1);
                else
                    newW = ceil((double) oldW * 0.9);

                vector<int>::const_iterator low;
                low = lower_bound(lpf.vX.begin(), lpf.vX.end(), vX);
                int pos = (int) (low - lpf.vX.begin());

                if (lpf.vX[pos] == vX)
                {
                    lpf.vY[pos] = newW;
                } else if (lpf.vX[pos] > vX)
                {
                    lpf.vX.insert(low, vX);
                    lpf.vY.insert(lpf.vY.begin() + pos, newW);
//                    assert(lpf.vX[pos] == vX);
                }

                vector<int> newX = lpf.vX, newY = lpf.vY;
                lpf.setValue(newX, newY, -1, acc);
                for (int j = 0; j < lpf.vX.size() - 1; j++)
                {
                    lpf.vSupportPre.push_back({{lId, {j}}});
                }

                break;
            }
        }

        int LCAId = LCA(u, v);
//        if (HGPTree[LCAId].dynMat[u].find(v) == HGPTree[LCAId].dynMat[u].end())
//        {
//            cout << "LCAId: " << LCAId << " u: " << u << " v: " << v << endl;
//            assert(false);
//        }
        UQ.update(LCAId, HGPTree.size() - LCAId);
    }

    int nid, key;
    while (!UQ.empty())
    {
        UQ.extract_min(nid, key);

        TDFloyd(HGPTree[nid]);

        if (nid == 0)
        {
            DQ.update(nid, 0);
            assert(UQ.empty());
            break;
        }

        set<int> uBorders;
        HGPTreeNode ni = HGPTree[nid];
        HGPTreeNode nf = HGPTree[ni.father];
        for (const auto &bid: ni.X)
        {
            if (nf.X.find(bid) != nf.X.end())
            {
                uBorders.insert(bid);
            }
        }

        bool dirty = false;
        for (const auto &bid1: uBorders)
        {
            for (const auto &bid2: uBorders)
            {
                if (bid1 == bid2)
                    continue;

                if (!nf.dynMat[bid1][bid2].equal(ni.dynMat[bid1][bid2], acc))
                {
                    if (nf.dynMat[bid1][bid2].vX.size() > 1 and
                        nf.dynMat[bid1][bid2].dominate(ni.dynMat[bid1][bid2], acc))
                    {
                        nf.dynMat[bid1][bid2] = nf.dynMat[bid1][bid2].LPFMinSupForDec(ni.dynMat[bid1][bid2], acc);
                    } else
                        nf.dynMat[bid1][bid2] = ni.dynMat[bid1][bid2];
                    dirty = true;
                }
            }
        }

        if (dirty)
        {
            UQ.update(ni.father, HGPTree.size() - ni.order);
            DQ.update(ni.father, ni.order);
        }
    }

    cout << "finish bot-up update" << endl;
    while (!DQ.empty())
    {
        DQ.extract_min(nid, key);
        HGPTreeNode &ni = HGPTree[nid];
        if (ni.isleaf)
        {
            continue;
        }

        for (const auto &gcId: ni.children)
        {
            set<int> uBorders;
            HGPTreeNode &nc = HGPTree[gcId];
            for (const auto &bid: ni.X)
            {
                if (nc.X.find(bid) != nc.X.end())
                {
                    uBorders.insert(bid);
                }
            }

            bool dirty = false;
            for (const auto &bid1: uBorders)
            {
                for (const auto &bid2: uBorders)
                {
                    if (bid1 == bid2)
                        continue;

                    if (!ni.dynMat[bid1][bid2].equal(nc.dynMat[bid1][bid2], acc))
                    {
                        dirty = true;
                        nc.dynMat[bid1][bid2] = ni.dynMat[bid1][bid2];
                    }
                }
            }

            if (dirty)
            {
                TDFloyd(nc);
                DQ.update(gcId, nc.order);
            }
        }
    }
}

int TDGTree::LCA(int u, int v)
{
    if (HGPNodes[u].treeid == HGPNodes[v].treeid)
    {
        return HGPNodes[u].treeid;
    }
    vector<int> trP1 = HGPNodes[u].HGPtreepath;
    vector<int> trP2 = HGPNodes[v].HGPtreepath;

    int pos = 0;
    while (pos < trP1.size() and pos < trP2.size() and trP1[pos] == trP2[pos])
    {
        pos++;
    }

    return trP1[pos - 1];
}

int TDGTree::TDSPQuery(int u, int v, int t)
{
//    t = t % graph.lBound;

    if (HGPNodes[u].treeid == HGPNodes[v].treeid)
    {
        int treeId = HGPNodes[u].treeid;
        if (HGPTree[treeId].isleaf)
        {
            if(HGPTree[treeId].dynMat[u].find(v) == HGPTree[treeId].dynMat[u].end())
                return -1;
            return HGPTree[treeId].dynMat[u][v].getY(t);
        }
    }

    int LCAId = LCA(u, v);
    if (HGPTree[LCAId].isleaf)
    {
        if (HGPTree[LCAId].dynMat[u].find(v) == HGPTree[LCAId].dynMat[u].end())
            return -1;
        return HGPTree[LCAId].dynMat[u][v].getY(t);
    }

    vector<int> path;
    int pos = HGPNodes[u].HGPtreepath.size() - 1;
    while (pos >= 0 and HGPNodes[u].HGPtreepath[pos] != LCAId)
    {
        path.emplace_back(HGPNodes[u].HGPtreepath[pos]);
        pos--;
    }

    pos = 0;
    while (pos < HGPNodes[v].HGPtreepath.size() and HGPNodes[v].HGPtreepath[pos] != LCAId)
    {
        pos++;
    }
    while (pos < HGPNodes[v].HGPtreepath.size())
    {
        path.emplace_back(HGPNodes[v].HGPtreepath[pos]);
        pos++;
    }

    if(path.empty())
        return -1;
    vector<int> vLPFs(graph.nodeNum, -1); // from u bid stored at vLPFs[bid]
    vLPFs[u] = t;
    for (int vb: HGPTree[path[0]].borders)
    {
        if (vb == u or HGPTree[path[0]].dynMat[u].find(vb) == HGPTree[path[0]].dynMat[u].end() or
            HGPTree[path[0]].dynMat[u][vb].vX.size() < 2)
            continue;
        vLPFs[vb] = t + HGPTree[path[0]].dynMat[u][vb].getY(t);
    }

    for (int j = 0; j < path.size() - 1; j++)
    {
        vector<int> X;
        for (const auto vb: HGPTree[path[j]].X)
        {
            if (HGPTree[path[j + 1]].X.find(vb) != HGPTree[path[j + 1]].X.end())
            {
                X.push_back(vb);
            }
        }

        for (const auto vb: HGPTree[path[j + 1]].X)
        {
            if (HGPTree[path[j]].X.find(vb) == HGPTree[path[j]].X.end())
            {
                int cost = INF;
                bool f = false;
                for (int x: X)
                {
                    if (x == u)
                        continue;
                    int cost1 = vLPFs[x];
                    if (cost1 < graph.uBound and cost1 >= graph.lBound)
                    {
                        f = true;
                        int cost2 = cost1 + HGPTree[path[j + 1]].dynMat[x][vb].getY(cost1);
                        if (cost2 < cost)
                            cost = cost2;
                    }
                }
                if (f)
                    vLPFs[vb] = cost;
            }
        }
    }

    int cost = INF;
    bool f = false;
    for (int vb: HGPTree[path.back()].borders)
    {
        int cost1 = vLPFs[vb];
        if (cost1 < graph.uBound and cost1 >= graph.lBound and vb != v and vb != u
            and HGPTree[path.back()].dynMat[vb].find(v) != HGPTree[path.back()].dynMat[vb].end())
        {
            f = true;
            int cost2 = cost1 + HGPTree[path.back()].dynMat[vb][v].getY(cost1);
            if (cost2 < cost)
                cost = cost2;
        }
    }
    if (f)
        cost -= t;
    return cost;
}

LPFunction TDGTree::TIPQuery(int u, int v)
{
    if (HGPNodes[u].treeid == HGPNodes[v].treeid)
    {
        int treeId = HGPNodes[u].treeid;
        if (HGPTree[treeId].isleaf)
        {
            if(HGPTree[treeId].dynMat[u].find(v) == HGPTree[treeId].dynMat[u].end())
                return LPFunction();
            return HGPTree[treeId].dynMat[u][v];
        }
    }

    int LCAId = LCA(u, v);
    if (HGPTree[LCAId].isleaf)
    {
        if(HGPTree[LCAId].dynMat[u].find(v) == HGPTree[LCAId].dynMat[u].end())
            return LPFunction();
        return HGPTree[LCAId].dynMat[u][v];
    }

    vector<int> path;

    int pos = HGPNodes[u].HGPtreepath.size() - 1;
    while (pos >= 0 and HGPNodes[u].HGPtreepath[pos] != LCAId)
    {
        path.emplace_back(HGPNodes[u].HGPtreepath[pos]);
        pos--;
    }

    pos = 0;
    while (pos < HGPNodes[v].HGPtreepath.size() and HGPNodes[v].HGPtreepath[pos] != LCAId)
    {
        pos++;
    }
    while (pos < HGPNodes[v].HGPtreepath.size())
    {
        path.emplace_back(HGPNodes[v].HGPtreepath[pos]);
        pos++;
    }

    if(path.empty())
        return LPFunction();
    vector<LPFunction> vLPFs(graph.nodeNum, LPFunction()); // from u bid stored at vLPFs[bid]
    for (const auto &bid: HGPTree[path[0]].borders)
    {
        if (bid == u)
            continue;

        if (HGPTree[path[0]].dynMat[u].find(bid) == HGPTree[path[0]].dynMat[u].end()
            or HGPTree[path[0]].dynMat[u][bid].vX.size() < 2)
        {
//            cout << "u: " << u << " bid: " << bid << endl;
//            for (const auto &vid: path)
//            {
//                cout << vid << " ";
//            }
//            cout << endl;

            for (const auto &vid: HGPNodes[u].HGPtreepath)
            {
                cout << vid << " ";
            }
            cout << endl;

            for (const auto &vid: HGPNodes[v].HGPtreepath)
            {
                cout << vid << " ";
            }
            cout << endl;

            return LPFunction();
//            assert(false);
        }
        vLPFs[bid] = HGPTree[path[0]].dynMat[u][bid];
    }

    for (int j = 0; j < path.size() - 1; j++)
    {
        vector<int> X;
        for (const auto vb: HGPTree[path[j]].X)
        {
            if (HGPTree[path[j + 1]].X.find(vb) != HGPTree[path[j + 1]].X.end())
            {
                X.push_back(vb);
            }
        }

        for (const auto &bid: HGPTree[path[j + 1]].X)
        {
            if (HGPTree[path[j]].X.find(bid) == HGPTree[path[j]].X.end())
            {
                LPFunction lpfMin(u, bid, graph.lBound, graph.uBound);
                for (const auto &xid: X)
                {
                    LPFunction lpf = vLPFs[xid];
                    LPFunction lpf2 = HGPTree[path[j + 1]].dynMat[xid][bid];
                    if (xid == u)
                    {
                        vLPFs[bid] = lpf2;
                    } else
                    {
                        if (lpf.vX.size() < 2 or lpf2.vX.size() < 2)
                            continue;

                        LPFunction lpfTmp = lpf.LPFCatSupport(lpf2, graph.lBound, graph.uBound, acc);
                        if (lpfMin.vX.size() < 2)
                        {
                            lpfMin = lpfTmp;
                        } else if (lpf.minY + lpf2.minY <= lpfMin.maxY and !lpfMin.dominate(lpfMin, acc))
                        {
                            lpfMin = lpfMin.LPFMinSupForDec(lpfTmp, acc);
                        }
                    }
                }
                if (lpfMin.vX.size() > 1)
                    vLPFs[bid] = lpfMin;
            }
        }
    }

    LPFunction lpfMin(u, v, graph.lBound, graph.uBound);
    for (const auto &bid: HGPTree[path.back()].borders)
    {
        if (bid == u or bid == v)
            continue;
        LPFunction lpf = vLPFs[bid];
        if (lpf.vX.size() < 2 or HGPTree[path.back()].dynMat[bid].find(v) == HGPTree[path.back()].dynMat[bid].end())
        {
            continue;
        }
        LPFunction lpf2 = HGPTree[path.back()].dynMat[bid][v];
        if (lpf2.vX.size() < 2)
        {
            continue;
        }
        LPFunction lpfTmp = lpf.LPFCatSupport(lpf2, graph.lBound, graph.uBound, acc);
        if (lpfMin.vX.size() < 2)
        {
            lpfMin = lpfTmp;
        } else if (lpf.minY + lpf2.minY <= lpfMin.maxY and !lpfMin.dominate(lpfTmp, acc))
        {
            lpfMin = lpfMin.LPFMinSupForDec(lpfTmp, acc);
        }
    }

    return lpfMin;
}

void TDGTree::ArchiveDynM(string &path) const
{
    for (int i = 0; i < HGPTree.size(); i++)
    {
        HGPTreeNode nd = HGPTree[i];
        for (auto &bid1: nd.X)
        {
            for (auto &bid2: nd.X)
            {
                if (bid1 == bid2)
                    continue;
                if (nd.dynMat[bid1].find(bid2) != nd.dynMat[bid1].end() and nd.dynMat[bid1][bid2].vX.size() > 1)
                {
                    vector<int> v = {i};
//                    nd.dynMat[bid1][bid2].writeLPF(path, v);
                }
            }
        }
    }
}

void TDGTree::ReadDynM(string &path)
{
    ifstream in(path.c_str());
    if (!in)
    {
        cout << "Cannot open file " << path << endl;
        return;
    }

    int i = -1;
    vector<string> lpfStrings;

    string line;
    while (getline(in, line))
    {
        if (line.find("LPF ") != string::npos)
        {
            if (!lpfStrings.empty())
            {
                assert(0 <= i < HGPTree.size());
//                LPFunction lpf = LPFunction::readLPFText(lpfStrings);

//                assert(HGPTree[i].X.find(lpf.ID1) != HGPTree[i].X.end()
//                       && HGPTree[i].X.find(lpf.ID2) != HGPTree[i].X.end());
//
//                HGPTree[i].dynMat[lpf.ID1][lpf.ID2] = lpf;
                lpfStrings.clear();
            }

            i = atoi(line.substr(4).c_str());

        } else if (!line.empty() and line != "")
        {
            lpfStrings.emplace_back(line);
        }
    }
    in.close();
    cout << "Finish reading dynMat" << endl;
}