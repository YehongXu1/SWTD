#include "graph.h"

void Graph::loadTD()
{
    TDConstruction();
}

void Graph::TDConstruction()
{
//    cout << "TD Construction" << endl;

    Timer clock;
//    clock.tick();
    CHConstruction();
//    clock.tock();
//    cout << "Contraction time: " << clock.duration().count() << endl;

    makeTree();

//    int max = 0;
//    int total = 0;
//    int ancSize = 0;
//    int maxAncSize = 0;
//    int maxWidth = 0;
//    int totalWidth = 0;
//    for (const auto &xx: Tree)
//    {
//        int s = (int) xx.vChildren.size();
//        int w = (int) xx.vertOut.size();
//        if (w > maxWidth)
//            maxWidth = w;
//        totalWidth += w;
//        if (s > max)
//            max = s;
//        total += s;
//        if (s == 0)
//        {
//            ancSize += (int) xx.ancIDs.size();
//            if (xx.ancIDs.size() > maxAncSize)
//                maxAncSize = (int) xx.ancIDs.size();
//        }
//    }
//
//    cout << "max childSize no:" << max << "\ttotal childSize no:" << total << endl;
//    cout << "max tree width:" << maxWidth << "\ttotal tree width:" << totalWidth << endl;
//    cout << "max anc size:" << maxAncSize << "\ttotal anc size:" << ancSize << endl;

    makeRMQ();

//    clock.tick();
    makeIndex();
//    clock.tock();
//    cout << "Making index time: " << clock.duration().count() << endl;

//    lBound += itvLen;
//    uBound += itvLen;
//    clock.tick();
//    extendIndex();
//    clock.tock();
//    cout << "Extend index time: " << clock.duration().count() << endl;
}

void Graph::extendIndex()
{
    lBound += itvLen;
    uBound += itvLen;
    vector<pair<int, int>> vpID;
    int thisThreadNum = threadNum;
    int step = nodeNum / thisThreadNum;

    vpID.reserve(thisThreadNum);
    for (int i = 0; i < thisThreadNum - 1; i++)
    {
        vpID.emplace_back(i * step, (i + 1) * step);
    }
    vpID.emplace_back((thisThreadNum - 1) * step, nodeNum);


    boost::thread_group threads;
    for (int i = 0; i < thisThreadNum; i++)
    {
        threads.add_thread(
                new boost::thread(
                        &Graph::lastChangedPosInit, this, vpID[i].first, vpID[i].second));
    }
    threads.join_all();
}

void Graph::lastChangedPosInit(int begin, int end)
{
    for (int i = begin; i < end; i++)
    {
        for (auto &j: adjEdge[i])
        {
            vEdge[j.second].lpf.extendFunction(lBound, uBound);
        }

        for (auto &j: Tree[i].vertIn)
        {
            j.second.first.extendFunction(lBound, uBound);
        }

        for (auto &j: Tree[i].vertOut)
        {
            j.second.first.extendFunction(lBound, uBound);
        }

        for (auto &j: Tree[i].disIn)
        {
            j.extendFunction(lBound, uBound);
        }

        for (auto &j: Tree[i].disOut)
        {
            j.extendFunction(lBound, uBound);
        }
    }
}

void Graph::CHConstruction()
{
    intermediateSCs.clear();
    intermediateSCs.assign(nodeNum, unordered_map < int, unordered_map < int, CatSupRec >> ());

    intermediateLBsOut.clear();
    intermediateLBsOut.assign(nodeNum, unordered_map < int, unordered_map < int, CatSupRec >> ());

    intermediateLBsIn.clear();
    intermediateLBsIn.assign(nodeNum, unordered_map < int, unordered_map < int, CatSupRec >> ());


    SCconNodesMT.assign(nodeNum, {});
    int debugID2 = 1, debugID1 = 1;
    //Contracted Graph E
    //Neighbor, Function
    map<int, pair<LPFunction, int>> m;
    E.assign(nodeNum, m);
    ER.assign(nodeNum, m);
    for (int i = 0; i < (int) adjEdge.size(); i++)
    {
        for (int j = 0; j < (int) adjEdge[i].size(); j++)
        {
            int u = i, v = adjEdge[i][j].first;
            LPFunction lpf = vEdge[adjEdge[i][j].second].lpf;
            E[u].insert({v, {lpf, 1}});

            CatSupRec catSupRec;
            catSupRec.vX1 = lpf.vX;
            catSupRec.vX2 = lpf.vX;
            catSupRec.vY = lpf.vY;
//            catSupRec.extendEnd(uBound + itvLen, itvLen);
            if (intermediateSCs[u].find(u) == intermediateSCs[u].end())
                intermediateSCs[u].insert({u, {{v, catSupRec}}});
            else if (intermediateSCs[u][u].find(v) == intermediateSCs[u][u].end())
                intermediateSCs[u][u].insert({v, catSupRec});
            else
                intermediateSCs[u][u][v] = catSupRec;
        }
    }

    for (int i = 0; i < (int) adjEdgeR.size(); i++)
    {
        for (int j = 0; j < (int) adjEdgeR[i].size(); j++)
            // there is an edge (j, i) in the original graph
            ER[i].insert({adjEdgeR[i][j].first, {vEdge[adjEdgeR[i][j].second].lpf, 1}});
    }

    vector<bool> exist(nodeNum, true);
    vector<bool> change(nodeNum, false);

    for (int i = 0; i < nodeNum; i++)
    {
        Semaphore *sm = new Semaphore(1);
        vSm.push_back(std::move(sm));
    }

    //Contracted Neighbors
    //Neighbor, <Distance, Supportive Node Number>
    vector<pair<int, pair<LPFunction, int>>> vect;
    NeighborCon.assign(nodeNum, vect);
    NeighborConR.assign(nodeNum, vect);
    //	SCconNodes.clear();

//    cout << "** begin to contract" << endl;
    int count = 0;
    map<int, pair<LPFunction, int>> mlpf;
//    vmSE.assign(nodeNum, mlpf);
//    vmSER.assign(nodeNum, mlpf);

    for (auto &x: vNodeOrder)
    {
        count++;
//        if (x == 48251)
//            cout << 1;
//        if (count % 50000 == 0)
//            cout << "count: " << count << endl;
        vector<pair<int, pair<LPFunction, int>>> Neigh, NeighR;
        for (auto &it: E[x])
            if (exist[it.first])
                Neigh.emplace_back(it);

        for (auto &it: ER[x])
            if (exist[it.first])
                NeighR.emplace_back(it);

        NeighborCon[x].assign(Neigh.begin(), Neigh.end());
        NeighborConR[x].assign(NeighR.begin(), NeighR.end());

        //Maintain E
        for (auto &i: Neigh)
        {
            int y = i.first;
            deleteE(x, y);
        }

        for (auto &i: NeighR)
        {
            int y = i.first;
            deleteE(y, x);
        }

        int ID1, ID2;
        LPFunction lpfTmp;
//        cout << x << " " << Neigh.size() << endl;
        if (NeighR.size() * Neigh.size() > 500 && NeighR.size() > 30)
        {
            vector<pair<int, int>> vpID;
            vpID.reserve(threadNum);
            int step = NeighR.size() / threadNum;
            for (int i = 0; i < threadNum - 1; i++)
            {
                vpID.emplace_back(i * step, (i + 1) * step);
            }
            vpID.emplace_back((threadNum - 1) * step, NeighR.size());

            for (auto &i: NeighR)
            {
                for (auto &j: Neigh)
                {
                    ID1 = i.first;
                    ID2 = j.first;
                    if (ID1 == ID2)
                        continue;
                    if (ID1 == debugID1 and ID2 == debugID2)
                    {
                        cout << "x:" << x << endl;
                        i.second.first.display();
                        j.second.first.display();
                    }

                    LPFunction lpfDummy;
                    lpfDummy.ID1 = ID1;
                    lpfDummy.ID2 = ID2;
                    if (E[ID1].find(ID2) == E[ID1].end())
                        E[ID1].insert({ID2, {lpfDummy, 1}});

                    if (ER[ID2].find(ID1) == ER[ID2].end())
                        ER[ID2].insert({ID1, {lpfDummy, 1}});

                    SCconNodesMT[ID1][ID2].emplace_back(x);
                    if (intermediateSCs[ID1].find(x) == intermediateSCs[ID1].end())
                    {
                        intermediateSCs[ID1].insert({x, {}});
                    }
                    if (intermediateSCs[ID1][x].find(ID2) == intermediateSCs[ID1][x].end())
                    {
                        CatSupRec catSupRec;
                        intermediateSCs[ID1][x].insert({ID2, catSupRec});
                    }

//                    if (vmSE[ID1].find(ID2) == vmSE[ID1].end())
//                        vmSE[ID1][ID2] = {lpfDummy, 1};
//
//                    if (vmSER[ID2].find(ID1) == vmSER[ID2].end())
//                        vmSER[ID2][ID1] = {lpfDummy, 1};
                }
            }
            boost::thread_group threads;
            for (int i = 0; i < threadNum; i++)
            {
                threads.add_thread(
                        new boost::thread(
                                &Graph::TDContractionThread, this, x, boost::ref(vpID[i]),
                                boost::ref(Neigh), boost::ref(NeighR)
                        ));
            }
            threads.join_all();
        } else
        {
            for (auto &i: NeighR)
            {
                for (auto &j: Neigh)
                {
                    ID1 = i.first;
                    ID2 = j.first;
                    if (ID1 == ID2)
                        continue;

                    if (ID1 == debugID1 and ID2 == debugID2)
                    {
                        cout << "x:" << x << endl;
                        i.second.first.display();
                        j.second.first.display();
                    }

//                    if (j.second.first.vX.size() <= 1)
//                        continue;

                    if (intermediateSCs[ID1].find(x) == intermediateSCs[ID1].end())
                    {
                        intermediateSCs[ID1].insert({x, {}});
                    }
                    CatSupRec catSupRec;
                    lpfTmp = i.second.first.LPFCatSupport(j.second.first, lBound, uBound, catSupRec, acc);
//                    catSupRec.extendEnd(uBound + itvLen, itvLen);
                    intermediateSCs[ID1][x].insert({ID2, catSupRec});

                    SCconNodesMT[ID1][ID2].push_back(x);
                    insertE(ID1, ID2, lpfTmp);

//                    if (vmSE[ID1].find(ID2) == vmSE[ID1].end() or vmSE[ID1][ID2].first.vX.size() <= 1)
//                        vmSE[ID1][ID2] = {lpfTmp, 1};
//                    else
//                        vmSE[ID1][ID2] = {vmSE[ID1][ID2].first.LPFMinSupForDec(lpfTmp), 1};
//
//
//                    if (vmSER[ID2].find(ID1) == vmSER[ID2].end() or vmSER[ID2][ID1].first.vX.size() <= 1)
//                        vmSER[ID2][ID1] = {lpfTmp, 1};
//                    else
//                        vmSER[ID2][ID1] = {vmSER[ID2][ID1].first.LPFMinSupForDec(lpfTmp), 1};
                }
            }
        }
    }

//    cout << "** finish contraction" << endl;
}

void Graph::TDContractionThread(
        int x, pair<int, int> pID, vector<pair<int, pair<LPFunction, int>>> &Neigh,
        vector<pair<int, pair<LPFunction, int>>> &NeighR)
{
    int ID1, ID2;
    int skip = 0;
    int cat = 0;
    int min = 0;
    for (int i = pID.first; i < pID.second; i++)
    {
        ID1 = NeighR[i].first;
        if (NeighR[i].second.first.vX.size() <= 1)
        {
//            cout << "Left Size:" << NeighR[i].second.first.vX.size() << endl;
            skip++;
            continue;
        }

        for (auto &j: Neigh)
        {
            ID2 = j.first;

            if (ID1 == 1 and ID2 == 1)
                cout << 1;

            if (ID1 == ID2)
            {
                skip++;
                continue;
            }

            if (j.second.first.vX.size() <= 1)
            {
                skip++;
                continue;
            }

            //E[ID1][ID2] exists a function and it is better

            if (E[ID1][ID2].first.vX.size() > 1 &&
                NeighR[i].second.first.minY + j.second.first.minY > E[ID1][ID2].first.maxY)
            {
                skip++;
                continue;
            }

            CatSupRec &catSupRec = intermediateSCs[ID1][x][ID2];
            LPFunction lpfTmp = NeighR[i].second.first.LPFCatSupport(j.second.first, lBound, uBound, catSupRec, acc);
//            catSupRec.extendEnd(uBound + itvLen, itvLen);

            cat++;
            if (lpfTmp.vX.size() <= 1)
            {
                skip++;
                continue;
            }

            LPFunction lpfExist = E[ID1][ID2].first;

            if (lpfExist.vX.size() > 1)
            {
                min++;
                if (!lpfExist.dominate(lpfTmp, acc))
                {
                    E[ID1][ID2] = {E[ID1][ID2].first.LPFMinSupByPlaneSweep(lpfTmp, acc), 1};
                }
            } else
            {
                E[ID1][ID2] = {lpfTmp, 1};
            }


            vSm[ID2]->wait();
            ER[ID2][ID1] = E[ID1][ID2];
            vSm[ID2]->notify();

//            vmSER[ID2][ID1] = ER[ID2][ID1];
//            vmSE[ID1][ID2] = vmSER[ID2][ID1];

//            assert(lpfTmp.ID1 != lpfTmp.ID2);
        }
    }
//    int total = NeighR.size() * Neigh.size();
//    cout << "skip:" << skip << "\t" << (double) skip / total << "\t" << total - skip << endl;
//    cout << "cat:" << cat << "\tmin:" << min << endl;
}

void Graph::makeTree()
{
    vector<int> vecemp;
    VidtoTNid.assign(nodeNum, vecemp);
    rank.clear();
    rank.assign(nodeNum, 0);
    int len = vNodeOrder.size() - 1;
    heightMax = 0;

    TreeNode rootNode;
    int x = vNodeOrder[len];
    assert(x >= 0);
    while (x == -1) //to skip those vertices whose ID is -1
    {
        len--;
        x = vNodeOrder[len];
    }
    rootNode.vertOut = NeighborCon[x];
    rootNode.vertIn = NeighborConR[x];
    rootNode.uniqueVertex = x;
    rootNode.parentNodeID = -1;
    rootNode.height = 1;
    rank[x] = 0;
    Tree.emplace_back(rootNode);
    len--;

    int nn;
    for (; len >= 0; len--)
    {
        int y = vNodeOrder[len];
        if (NeighborCon[y].empty() || NeighborConR[y].empty())
            continue;
        TreeNode nodeTmp;
        nodeTmp.vertOut = NeighborCon[y];
        nodeTmp.vertIn = NeighborConR[y];

        nodeTmp.uniqueVertex = y;

        int paRank;
        if (!NeighborCon[y].empty())
            paRank = deepestAnc(NeighborCon[y]);
        else
            paRank = deepestAnc(NeighborConR[y]);

        Tree[paRank].vChildren.emplace_back(Tree.size());
        nodeTmp.parentNodeID = paRank;
        nodeTmp.height = Tree[paRank].height + 1;
        nodeTmp.hdepth = Tree[paRank].height + 1;

        for (auto &i: NeighborCon[y])
        {
            nn = i.first;
            VidtoTNid[nn].emplace_back(Tree.size());
            if (Tree[rank[nn]].hdepth < Tree[paRank].height + 1)
                Tree[rank[nn]].hdepth = Tree[paRank].height + 1;
        }

        if (nodeTmp.height > heightMax)
            heightMax = nodeTmp.height;

        rank[y] = Tree.size();
        Tree.emplace_back(nodeTmp);
    }
//    cout << "** making tree" << ", max height: " << heightMax << endl;
}

void Graph::makeRMQ()
{
    // explanation in https://www.youtube.com/watch?v=sD1IoalFomA
    //EulerSeq.clear();
    toRMQ.assign(nodeNum, 0);
    //RMQIndex.clear();
    makeRMQDFS(0, 1);
    RMQIndex.push_back(EulerSeq);

    int m = EulerSeq.size();
    for (int i = 2, k = 1; i < m; i = i * 2, k++)
    {
        vector<int> tmp;
        //tmp.clear();
        tmp.assign(m, 0);
        for (int j = 0; j < m - i; j++)
        {
            int x = RMQIndex[k - 1][j], y = RMQIndex[k - 1][j + i / 2];
            if (Tree[x].height < Tree[y].height)
                tmp[j] = x;
            else tmp[j] = y;
        }
        RMQIndex.push_back(tmp);
    }
}

int Graph::LCAQuery(int _p, int _q)
{
    int p = toRMQ[_p], q = toRMQ[_q];
    if (p > q)
    {
        int x = p;
        p = q;
        q = x;
    }
    int len = q - p + 1;
    int i = 1, k = 0;
    while (i * 2 < len)
    {
        i *= 2;
        k++;
    }
    q = q - i + 1;
    if (Tree[RMQIndex[k][p]].height < Tree[RMQIndex[k][q]].height)
        return RMQIndex[k][p];
    else return RMQIndex[k][q];
}

LPFunction Graph::QueryH2H(int ID1, int ID2, long long int &LCASize)
{
    if (ID1 == ID2)
        return LPFunction();
    int r1 = rank[ID1], r2 = rank[ID2];
    int LCA = LCAQuery(r1, r2);

    if (LCA == r1)
        return Tree[r2].disIn[Tree[r1].posIn.back()];
    else if (LCA == r2)
        return Tree[r1].disOut[Tree[r2].posOut.back()];
    else
    {
        LCASize += Tree[LCA].disIn.size();

        LPFunction lpfMin;
        for (int i = 0; i < Tree[LCA].posOut.size(); i++)
        {
            int x = Tree[LCA].posOut[i]; // x-th ancestor of LCA
            if (lpfMin.vX.size() < 2)
            {
                // lpfInit is not initialized
                lpfMin = Tree[r1].disOut[x].LPFCatSupport(Tree[r2].disIn[x], lBound, uBound, acc);
            } else
            {
                // lpfInit is initialized
                if (Tree[r1].disOut[x].minY + Tree[r2].disIn[x].minY > lpfMin.maxY)
                    continue;
                LPFunction disTmp = Tree[r1].disOut[x].LPFCatSupport(Tree[r2].disIn[x], lBound, uBound, acc);
                if (!lpfMin.dominate(disTmp, acc))
                {
                    lpfMin = lpfMin.LPFMinSupForDec(disTmp, acc);
                }
            }
//            for (int j = 0; j < Tree[LCA].posIn.size(); j++)
//            {
//                int y = Tree[LCA].posIn[j]; // y-th ancestor of LCA
//                if (x == y)
//                {
//
//                    break;
//                }
//            }
        }
        return lpfMin;
    }
}

int Graph::QueryH2HFixedT(int ID1, int ID2, int t, long long int &LCASize)
{
    if (ID1 == ID2)
        return 0;
    int r1 = rank[ID1], r2 = rank[ID2];
    int LCA = LCAQuery(r1, r2);

    int departT = t;
    if (LCA == r1)
        return Tree[r2].disIn[Tree[r1].posIn.back()].getY(departT);
    else if (LCA == r2)
        return Tree[r1].disOut[Tree[r2].posOut.back()].getY(departT);
    else
    {
        LCASize += Tree[LCA].ancIDs.size();
        int costMin = INF;
        for (int i = 0; i < Tree[LCA].posOut.size(); i++)
        {
            int x = Tree[LCA].posOut[i]; // x-th ancestor of LCA
            int cost1 = Tree[r2].disIn[x].getY(departT);
            int cost = Tree[r1].disOut[x].getY(cost1 + departT);
            if (cost1 + cost < costMin)
                costMin = cost + cost1;

        }
        return costMin;
    }
}

int Graph::QueryCHFixedT(int ID1, int ID2, int t)
{
    benchmark::heap<2, int, int> queueFWD(nodeNum), queueBKD(nodeNum);
    queueFWD.update(ID1, t), queueBKD.update(ID2, 0);

    vector<int> distFWD(nodeNum, INF), distBKD(nodeNum, INF);
    vector<bool> visitedFWD(nodeNum, false), visitedBKD(nodeNum, false);
    int u, distU;
    int v, distV;

    distFWD[ID1] = t;
    distBKD[ID2] = 0;
    int mu = INF;

    vector<int> succVidBKD(nodeNum, -1);
    while (!queueFWD.empty() and !queueBKD.empty())
    {
        queueFWD.extract_min(u, distU), queueBKD.extract_min(v, distV);
        visitedFWD[u] = true, visitedBKD[v] = true;
        if (mu < distV or u == ID2)
            break;
        for (const auto &p: Tree[rank[u]].vertOut)
        {
            int x = p.first;
            int neighborLength = p.second.first.getY(distU);
            assert(neighborLength + distU < uBound);
            if (!visitedFWD[x] and neighborLength + distU < distFWD[x])
            {
                distFWD[x] = neighborLength + distU;
                queueFWD.update(x, distFWD[x]);
            }

            if (visitedBKD[x])
            {
                int vid = x, succVid = succVidBKD[x]; // vid -> succVid, vid contracted later
                int dist = distFWD[x];
                while (succVid >= 0)
                {
                    for (const auto &pp: Tree[rank[succVid]].vertIn)
                    {
                        if (pp.first == vid)
                        {
                            dist += pp.second.first.getY(dist);
                            break;
                        }
                    }

                    vid = succVid;
                    succVid = succVidBKD[succVid];
                }
                if (mu > dist)
                    mu = dist;
            }

        }

        for (const auto &p: Tree[rank[v]].vertIn)
        {
            int x = p.first;
            int neighborLength = p.second.first.minY;
            if (!visitedBKD[x] and !visitedFWD[x] and neighborLength + distV < distBKD[x])
            {
                distBKD[x] = neighborLength + distV;
                succVidBKD[x] = v;
                queueBKD.update(x, distBKD[x]);
            }
        }

    }

    return mu - t;
}

LPFunction Graph::QueryCHItvT(int ID1, int ID2, int &miny)
{
    benchmark::heap<2, int, int> queFWD(nodeNum + 1), queBKD(nodeNum + 1);
    vector<int> vTimeFWD(nodeNum + 1, INF), vTimeBKD(nodeNum + 1, INF);
    vector<bool> inHeapFWD(nodeNum + 1, false), inHeapBKD(nodeNum + 1, false); //True if in Heap
    vector<bool> visitedFWD(nodeNum + 1, false), visitedBKD(nodeNum + 1, false);
    vector<bool> reachedFWD(nodeNum + 1, false), reachedBKD(nodeNum + 1, false); // True if ever reached
    vector<LPFunction> distFWD(nodeNum + 1, LPFunction()), distBKD(nodeNum + 1, LPFunction());
    vector<int> succIDVec(nodeNum + 1, -1);

    int u, distU = 0, v, distV = 0, x;

    bool xUpdated;
    LPFunction mu = LPFunction();

    vTimeFWD[ID1] = lBound, vTimeBKD[ID2] = lBound;
    inHeapFWD[ID1] = true, inHeapBKD[ID2] = true;;
    visitedFWD[ID1] = true, visitedBKD[ID2] = true;
    reachedFWD[ID1] = true, reachedBKD[ID2] = true;
    queFWD.update(ID1, lBound), queBKD.update(ID2, lBound);

    while (!queFWD.empty() and !queBKD.empty())
    {
        queFWD.extract_min(u, distU), queBKD.extract_min(v, distV);
//        if (mu.dominate(distFWD[]))
//            break;
        inHeapFWD[u] = false, inHeapBKD[v] = false;
        visitedFWD[u] = true, visitedBKD[v] = true;


//        if (u == ID2 or mu.minY < distV)
//            break;

        for (const auto &p: Tree[rank[u]].vertOut)
        {
            xUpdated = false;
            x = p.first;

            LPFunction lpf;
            if (u == ID1)
            {
                lpf = p.second.first;
            } else
            {
                LPFunction lpftmp = p.second.first;
                lpf = distFWD[u].LPFCatSupport(lpftmp, lBound, uBound, acc);
            }

            if (lpf.vX.size() == 1)
                continue;

            if (reachedFWD[x])
            {
                if (distFWD[x].vX.size() > 1 and distFWD[x].dominate(lpf, acc))
                    continue;
                distFWD[x] = distFWD[x].LPFMinSupForDec(lpf, acc);
                xUpdated = true;
            } else
            {
                distFWD[x] = lpf;
                reachedFWD[x] = true;
                xUpdated = true;
            }

            int minCost = distFWD[x].minY;

            //Updated and not in Heap
            if ((!inHeapFWD[x] && xUpdated) || !visitedFWD[x])
            {
                vTimeFWD[x] = minCost;
                queFWD.update(x, minCost);
                inHeapFWD[x] = true;
                visitedFWD[x] = true;
            }
                //Updated and in Heap, key changed
            else if (xUpdated && inHeapFWD[x])
            {
                vTimeFWD[x] = minCost;
                queFWD.update(x, minCost);
            }

            if (visitedBKD[x])
            {
                int vid = x, succVid = succIDVec[x]; // vid -> succVid, vid contracted later
                LPFunction dist = distFWD[x];
                while (succVid >= 0)
                {
                    for (const auto &pp: Tree[rank[succVid]].vertIn)
                    {
                        if (pp.first == vid)
                        {
                            LPFunction copy = pp.second.first;
                            dist = dist.LPFCatSupport(copy, lBound, uBound, acc); //pp.second.first.getY(dist);
                            break;
                        }
                    }

                    vid = succVid;
                    succVid = succIDVec[succVid];
                }
                if (mu.vX.size() < 2)
                    mu = dist;
                else if (!mu.dominate(dist, acc))
                {
                    mu = mu.LPFMinSupForDec(dist, acc);
                } else
                {
                    return mu;
                }
            }
        }

        for (const auto &p: Tree[rank[v]].vertIn)
        {
            xUpdated = false;
            x = p.first;
            if (visitedFWD[x])
                continue;
            LPFunction lpf;
            if (v == ID2)
            {
                lpf = p.second.first; // x - > v (ID2)
            } else
            {
                LPFunction lpftmp = p.second.first;
                lpf = lpftmp.LPFCatSupport(distBKD[v], lBound, uBound, acc); // x -> v -> ID2
            }

            if (lpf.vX.size() == 1)
                continue;

            if (reachedBKD[x])
            {
                //distFWD[x].display()
                //lpf.display();
                if (distBKD[x].vX.size() > 1 and distBKD[x].dominate(lpf, acc))
                    continue;
//                distFWD[x] = distFWD[x].LPFMinNew3(lpf);
                distBKD[x] = distBKD[x].LPFMinSupForDec(lpf, acc);
                xUpdated = true;
                succIDVec[x] = v;
            } else
            {
                distBKD[x] = lpf;
                reachedBKD[x] = true;
                xUpdated = true;
                succIDVec[x] = v;
            }

            int minCost = distBKD[x].minY;

            //Updated and not in Heap
            if ((!inHeapBKD[x] && xUpdated) || !visitedBKD[x])
            {
                vTimeBKD[x] = minCost;
                queBKD.update(x, minCost);
                inHeapBKD[x] = true;
                visitedBKD[x] = true;
            }
                //Updated and in Heap, key changed
            else if (xUpdated && inHeapBKD[x])
            {
                vTimeBKD[x] = minCost;
                queBKD.update(x, minCost);
            }
        }
    }
    return mu;
}

void Graph::makeRMQDFS(int p, int height)
{
    toRMQ[p] = EulerSeq.size();
    EulerSeq.push_back(p);
    for (int i = 0; i < Tree[p].vChildren.size(); i++)
    {
        makeRMQDFS(Tree[p].vChildren[i], height + 1);
        EulerSeq.push_back(p);
    }
}

void Graph::makeIndex()
{
    //initialize
    vector<int> list; //list.clear();
    Tree[0].ancIDs = {};
    list.emplace_back(Tree[0].uniqueVertex);
    Tree[0].posIn.clear();
    Tree[0].posOut.clear();
    Tree[0].posIn.emplace_back(0);
    Tree[0].posOut.emplace_back(0);

    Timer clock;
    clock.tick();
    makeIndexDFSThread();
    clock.tock();
//    cout << "** making index " << clock.duration().count() << " ms" << endl;
}

void Graph::makeIndexDFSThread()
{
    // ancs has the same size as the number of ancestors of p (i.e. the height of p)
    //initialize
    vector<pair<int, vector<int>>> S;
    for (int i: Tree[0].vChildren)
        S.push_back({i, {Tree[0].uniqueVertex}});

//    int height = 0;
    Timer clock;
    while (!S.empty() and S.size() < threadNum)
    {
        clock.tick();

        boost::thread_group threads;
        for (auto &i: S)
        {
            threads.add_thread(
                    new boost::thread(
                            &Graph::makeIndexDFS, this, i.first, boost::ref(i.second)
                    ));
        }
        threads.join_all();

//        int thisThreadNum = threadNum;
//        if (S.size() < thisThreadNum)
//        {
//
//        } else
//        {
//            vector<pair<int, int>> vPIDs;
//            vPIDs.reserve(thisThreadNum);
//            int step = (int) S.size() / thisThreadNum;
//            for (int i = 0; i < thisThreadNum - 1; i++)
//                vPIDs.emplace_back(i * step, (i + 1) * step);
//
//            vPIDs.emplace_back((thisThreadNum - 1) * step, S.size());
//
//            boost::thread_group threads;
//            for (int i = 0; i < thisThreadNum; i++)
//            {
//                threads.add_thread(
//                        new boost::thread(
//                                &Graph::makeIndexDFSRange, this,
//                                vPIDs[i].first, vPIDs[i].second, boost::ref(S)
//                        ));
//            }
//            threads.join_all();
//        }
//        cout << "height: " << height << " S: " << S.size() << " time: " << clock.duration().count() << " ms" << endl;
//        height++;

        vector<pair<int, vector<int>>> S2;
        S2.reserve(S.size() * 6);
        for (const auto &p: S)
        {
            vector<int> ancs = p.second;
            ancs.emplace_back(Tree[p.first].uniqueVertex);
            for (int i: Tree[p.first].vChildren)
            {
                S2.emplace_back(i, ancs);
            }
        }
        clock.tock();
        S = S2;
    }

    if (S.empty())
        return;
//    cout << "thread num: " << S.size() << endl;
    boost::thread_group threads;
    for (auto &i: S)
    {
        threads.add_thread(
                new boost::thread(
                        &Graph::makeIndexDFSRec, this, i.first, boost::ref(i.second)
                ));
    }
    threads.join_all();
}

void Graph::makeIndexDFSRec(int p, vector<int> &ancs)
{
    makeIndexDFS(p, ancs);
    ancs.emplace_back(Tree[p].uniqueVertex);
    for (int i: Tree[p].vChildren)
    {
        makeIndexDFSRec(i, ancs);
    }
    ancs.pop_back();
}

void Graph::makeIndexDFSRange(int begin, int end, vector<pair<int, vector<int>>> &pAncs)
{
    for (int i = begin; i < end; i++)
    {
        makeIndexDFS(pAncs[i].first, pAncs[i].second);
    }
}

void Graph::makeIndexDFS(int p, vector<int> &ancs)
{
    // ancs has the same size as the number of ancestors of p (i.e. the height of p)
    //initialize
    int NeiNumIn = Tree[p].vertIn.size(); // vectIn.size() == vectOut.size()
    int NeiNumOut = Tree[p].vertOut.size(); // vectIn.size() == vectOut.size()
//    if (p % 1000 == 0)
//    cout << "p: " << p << " height: " << ancs.size() << " " << NeiNumIn << endl;
    if (NeiNumIn == 0 and NeiNumOut == 0)
        return;
    Tree[p].posIn.assign(NeiNumIn + 1, 0);
    Tree[p].posOut.assign(NeiNumOut + 1, 0);
    Tree[p].ancIDs.assign(ancs.begin(), ancs.end());

    LPFunction lpf;
    Tree[p].disIn.assign(ancs.size(), lpf);
    Tree[p].disOut.assign(ancs.size(), lpf);
    Tree[p].cntIn.assign(ancs.size(), 1);
    Tree[p].cntOut.assign(ancs.size(), 1);

    //pos
    int pId = Tree[p].uniqueVertex;
//    if (pId == 1964)
//        cout << 1;
    intermediateLBsIn[pId] = {{pId, {}}};
    for (int i = 0; i < NeiNumIn; i++)
    {
        int aid = Tree[p].vertIn[i].first;
        for (int j = 0; j < ancs.size(); j++)
        {
            if (aid == ancs[j])
            {
                Tree[p].posIn[i] = j; // posIn[i] stores ancestor j of p
                Tree[p].disIn[j] = Tree[p].vertIn[i].second.first;
                Tree[p].disIn[j].scSupToLbSup(pId);

                CatSupRec catSupRec;
                catSupRec.vX1 = Tree[p].disIn[j].vX;
                catSupRec.vX2 = Tree[p].disIn[j].vX;
                catSupRec.vY = Tree[p].disIn[j].vY;
//                catSupRec.extendEnd(uBound + itvLen, itvLen);

                intermediateLBsIn[pId][pId].insert({aid, catSupRec});
                break;
            }
        }
    }
    Tree[p].posIn[NeiNumIn] = (int) ancs.size(); // ancs.size() == height - 1

    intermediateLBsOut[pId] = {{pId, {}}};
    for (int i = 0; i < NeiNumOut; i++)
    {
        int aid = Tree[p].vertOut[i].first;
        for (int j = 0; j < ancs.size(); j++)
        {
            if (aid == ancs[j])
            {
                Tree[p].posOut[i] = j;
                Tree[p].disOut[j] = Tree[p].vertOut[i].second.first;
                Tree[p].disOut[j].scSupToLbSup(pId);

                CatSupRec catSupRec;
                catSupRec.vX1 = Tree[p].disOut[j].vX;
                catSupRec.vX2 = Tree[p].disOut[j].vX;
                catSupRec.vY = Tree[p].disOut[j].vY;
//                catSupRec.extendEnd(uBound + itvLen, itvLen);

                intermediateLBsOut[pId][pId].insert({aid, catSupRec});
                break;
            }
        }
    }
    Tree[p].posOut[NeiNumOut] = (int) ancs.size();

    /* in-labels */
    for (int i = 0; i < NeiNumIn; i++)
    {
        if (Tree[p].vertIn[i].second.first.vX.size() < 2)
            continue;
        int x = Tree[p].vertIn[i].first;
        LPFunction scXP = Tree[p].vertIn[i].second.first; // LPF from x to p
        int k = Tree[p].posIn[i]; //the kth ancestor is x

        scXP.scSupToLbSup(pId);
        for (int j = 0; j < ancs.size(); j++)
        {
            // p <- x <- y
            int y = ancs[j]; // the jth ancestor is y
            LPFunction disYX; // LPF from y to x
            if (k < j)
            {
                // x is the ancestor of y
                disYX = Tree[rank[y]].disOut[k];
            } else if (k > j)
            {
                // y is the ancestor of x
                disYX = Tree[rank[x]].disIn[j];
            }
            LPFunction disYPExist = Tree[p].disIn[j];
            if (intermediateLBsIn[pId].find(x) == intermediateLBsIn[pId].end())
                intermediateLBsIn[pId].insert({x, {}});
            intermediateLBsIn[pId][x].insert({y, {}});
            CatSupRec &catSupRec = intermediateLBsIn[pId][x][y];
            if (disYPExist.vX.size() <= 1)
            {
                LPFunction disMin = disYX.LPFCatSupport(scXP, lBound, uBound, catSupRec, acc);
                Tree[p].disIn[j] = disMin;
            } else if (k != j and disYX.minY + scXP.minY <= disYPExist.maxY)
            {
                LPFunction disTmp = disYX.LPFCatSupport(scXP, lBound, uBound, catSupRec, acc);
                if (!disYPExist.dominate(disTmp, acc))
                    Tree[p].disIn[j] = disYPExist.LPFMinSupByPlaneSweep(disTmp, acc);
            }
        }
    }

    /* out-labels */
    for (int i = 0; i < NeiNumOut; i++)
    {
        if (Tree[p].vertOut[i].second.first.vX.size() < 2)
            continue;
        int x = Tree[p].vertOut[i].first;
        LPFunction scPX = Tree[p].vertOut[i].second.first; // LPF from p to x
        int k = Tree[p].posOut[i];//the kth ancestor is x

        scPX.scSupToLbSup(pId);
        for (int j = 0; j < ancs.size(); j++)
        {
            int y = ancs[j]; //the jth ancestor is y
            LPFunction disXY; // LPF from x to y
            if (k < j)
                // y is deeper
                disXY = Tree[rank[y]].disIn[k];
            else if (k > j)
                // x is deeper, contract earlier
                disXY = Tree[rank[x]].disOut[j];
            LPFunction disPYExist = Tree[p].disOut[j];
            if (intermediateLBsOut[pId].find(x) == intermediateLBsOut[pId].end())
                intermediateLBsOut[pId].insert({x, {}});
            intermediateLBsOut[pId][x].insert({y, {}});
            CatSupRec &catSupRec = intermediateLBsOut[pId][x][y];
            if (disPYExist.vX.size() <= 1)
            {
                LPFunction disMin = scPX.LPFCatSupport(disXY, lBound, uBound, catSupRec, acc);
                Tree[p].disOut[j] = disMin;
            } else if (k != j and scPX.minY + disXY.minY <= disPYExist.maxY)
            {
                LPFunction disTmp = scPX.LPFCatSupport(disXY, lBound, uBound, catSupRec, acc);
                if (!disPYExist.dominate(disTmp, acc))
                    Tree[p].disOut[j] = disPYExist.LPFMinSupByPlaneSweep(disTmp, acc);
            }
        }
    }
}

int Graph::deepestAnc(vector<pair<int, pair<LPFunction, int>>> &vert)
{
    if (vert.empty())
        cout << "Empty vert!" << endl;

    int nearest = vert[0].first;
    for (int i = 1; i < (int) vert.size(); i++)
    {
        if (rank[vert[i].first] > rank[nearest])
            nearest = vert[i].first;
    }
    int p = rank[nearest];
    return p;
}

void Graph::deleteE(int u, int v)
{
    if (E[u].find(v) != E[u].end())
        E[u].erase(E[u].find(v));

    if (ER[v].find(u) != ER[v].end())
        ER[v].erase(ER[v].find(u));
}

void Graph::insertE(int u, int v, LPFunction &lpf)
{
    LPFunction lpfExist;
    if (E[u].find(v) == E[u].end() or E[u][v].first.vX.size() <= 1)
    {
        E[u].insert({v, {lpf, 1}});
    } else
    {
        lpfExist = E[u][v].first;
        if (lpf.vX.size() > 1 and !lpfExist.dominate(lpf, acc))
        {
            E[u][v] = {lpfExist.LPFMinSupByPlaneSweep(lpf, acc), 1};
        }
    }

    ER[v][u] = E[u][v];
}

vector<int> NodeOrderss;

struct OrderCompp
{//prior to reture the vertex with smaller order
    int x;

    explicit OrderCompp(int _x)
    {
        x = _x;
    }

    bool operator<(const OrderCompp &d) const
    {
        if (x == d.x)
        {//avoid the redundant
            return false;
        } else
        {
            if (x != d.x)
                return NodeOrderss[x] < NodeOrderss[d.x];
        }
    }
};

void Graph::H2HDecBatDiGraph(vector<pair<int, int>> &wBatch, int timeWinId, int timeWinLen)
{
    cout << "Index maintenance for " << wBatch.size() << " decrease updates" << endl;
    assert(!wBatch.empty());
    // update weights of shortcuts
    // find the sources of top-down distance updates

    for (auto &i: Tree)
    {
        i.DisReIn.clear();//record the star weight change (causing the distance change)
        i.DisReOut.clear();//record the star weight change (causing the distance change)
    }

    // NodeOrderss.clear();
    NodeOrderss.assign(NodeOrder.begin(), NodeOrder.end());
    vector<set<int>> SCreOut(nodeNum, set < int > ()), SCreIn(nodeNum, set < int > ()); //SCre.clear();
    set < OrderCompp > OC; //OC.clear();//vertexID in decreasing Node order

    int debugID1 = 1, debugID2 = 1;

    for (const auto &k: wBatch)
    {
        int a, b, vX, newW, oldW;
        a = k.first;
        b = k.second;
        vX = lBound + 7 * 300;

        if (a == debugID1 and b == debugID2)
            cout << 1;
        LPFunction newEdge;
        bool outOfBound = true;
        for (const auto &i: adjEdge[a])
        {
            if (i.first == b)
            {
                LPFunction &lpf = vEdge[i.second].lpf;

                if (vX > lpf.vX.back() or vX < lpf.vX[0])
                {
                    continue;
                }

                outOfBound = false;
                oldW = lpf.getY(vX);
                newW = ceil((double) oldW * 1.0 - oldW * 0.1);

                vector<int>::const_iterator low;
                low = lower_bound(lpf.vX.begin(), lpf.vX.end(), vX);
                int pos = (int) (low - lpf.vX.begin());
                assert(pos >= 0 and pos < lpf.vX.size() and lpf.vX[pos] >= vX);

                if (lpf.vX[pos] == vX)
                {
                    lpf.vY[pos] = newW;
                } else if (lpf.vX[pos] > vX)
                {
                    lpf.vX.insert(low, vX);
                    lpf.vY.insert(lpf.vY.begin() + pos, newW);
                }

                lpf.lastItvSubToChange = pos > 0 ? pos - 1 : 0;

                vector<int> newX = lpf.vX, newY = lpf.vY;
                lpf.setValue(newX, newY, -1, 0);

                int lid = NodeOrder[a] < NodeOrder[b] ? a : b;
                lpf.vSupportPre.clear();
                for (int j = 0; j < lpf.vX.size() - 1; j++)
                {
                    lpf.vSupportPre.push_back({{lid, {j}}});
                }

                newEdge = lpf;
                break;
            }
        }

        if (newEdge.vX.size() < 2 or outOfBound)
            continue;

        if (NodeOrder[a] < NodeOrder[b])
        {
            for (int i = 0; i < Tree[rank[a]].vertOut.size(); i++)
            {
                if (Tree[rank[a]].vertOut[i].first == b)
                {
                    LPFunction lpfExist = Tree[rank[a]].vertOut[i].second.first;
                    if (!lpfExist.dominate(newEdge, acc))
                    {
                        Tree[rank[a]].vertOut[i].second = {lpfExist.LPFMinSupForDec(newEdge, acc), 1};
                        if (a == debugID1 and b == debugID2)
                        {
                            lpfExist.display();
                            newEdge.display();
                            Tree[rank[a]].vertOut[i].second.first.display();
                        }

                        Tree[rank[a]].DisReOut.insert(b);
                        SCreOut[a].insert(b);
                        OC.insert(OrderCompp(a)); // collect the source of updates
                    }
                    break;
                }
            }
        } else
        {
            for (int i = 0; i < Tree[rank[b]].vertIn.size(); i++)
            {
                if (Tree[rank[b]].vertIn[i].first == a)
                {
                    LPFunction lpfExist = Tree[rank[b]].vertIn[i].second.first;
                    if (!lpfExist.dominate(newEdge, acc))
                    {
                        Tree[rank[b]].vertIn[i].second = {lpfExist.LPFMinSupForDec(newEdge, acc), 1};
                        if (a == debugID1 and b == debugID2)
                        {
                            lpfExist.display();
                            newEdge.display();
                            Tree[rank[b]].vertIn[i].second.first.display();
                        }
                        SCreIn[b].insert(a);
                        OC.insert(OrderCompp(b)); // collect the source of updates
                    }
                    break;
                }
            }
        }
    }
//    cout << "** update origin edges" << endl;

    vector<int> ProBeginVertexSet; // source of top-down distance updates
    vector<int> ProBeginVertexSetNew;

    while (!OC.empty())
    {
        int ProID = (*OC.begin()).x; // lid
        OC.erase(OC.begin());
        vector<pair<int, LPFunction>> Vert; // accessories of ProID
        bool ProIDdisChaIn = false, ProIDdisChaOut = false; //to see if the distance labeling of proID change or not

        unordered_map<int, int> scProcThisTurn;
        for (const auto &Cid: SCreIn[ProID])
        {
            if (ProID == debugID2 and Cid == debugID1)
                cout << "ProId: " << ProID << " from Cid: " << Cid << endl;

            map<int, LPFunction> HneiOut;
            vector<pair<int, LPFunction>> LneiOut;
            for (auto &j: Tree[rank[ProID]].vertOut)
            {
                if (NodeOrder[j.first] > NodeOrder[Cid])
                    HneiOut[j.first] = j.second.first;
                else if (NodeOrder[j.first] < NodeOrder[Cid])
                    LneiOut.emplace_back(j.first, j.second.first);
            }

            LPFunction CidToPidFun;
            for (auto &j: Tree[rank[ProID]].vertIn)
            {
                if (j.first == Cid)
                {
                    // j.first == Cid, appear once
                    CidToPidFun = j.second.first;
                    break;
                }
            }

            int cidH = Tree[rank[Cid]].height - 1;
            LPFunction lpfExist = Tree[rank[ProID]].disIn[cidH];
            if (ProID == debugID2 and Cid == debugID1)
            {
                lpfExist.display();
                CidToPidFun.display();
//                lpfExist.LPFMinSupForDec(CidToPidFun).display();
            }
            if (!lpfExist.dominate(CidToPidFun, acc))
            {
                LPFunction CidToPidFun2 = CidToPidFun;
                CidToPidFun2.scSupToLbSup(ProID);

                Tree[rank[ProID]].disIn[cidH] = lpfExist.LPFMinSupForDec(CidToPidFun2, acc);
                ProIDdisChaIn = true;
                // DisRe contains the accessories whose distance labeling has changed
                Tree[rank[ProID]].DisReIn.insert(Cid);
            }

            // propogate the updates of tree edges to the ancestors of ProID
            int hid, lid;
            LPFunction existSC, lpfTmp;
            for (int j = 0; j < Tree[rank[Cid]].vertOut.size(); j++)
            {
                hid = Tree[rank[Cid]].vertOut[j].first;
                if (HneiOut.find(hid) != HneiOut.end())
                {
                    // otherwise, e(hid, cid) is not supported by (a,b)
                    // how to concatenate?
//                    CatSupRec catSupRec = intermediateLBsOut[Cid][hid][ProID];
                    lpfTmp = CidToPidFun.LPFCatSupport(HneiOut[hid], lBound, uBound, acc);
                    existSC = Tree[rank[Cid]].vertOut[j].second.first;

                    if (!existSC.dominate(lpfTmp, acc))
                    {
                        Tree[rank[Cid]].vertOut[j].second = {existSC.LPFMinSupForDec(lpfTmp, acc), 1};
                        SCreOut[Cid].insert(hid);
                        OC.insert(OrderCompp(Cid));
                        scProcThisTurn.insert({Cid, hid});
                        if (Cid == debugID1 and ProID == debugID2)
                        {
                            cout << "Cid: " << Cid << " to hid: " << hid << " via ProId: " << ProID << endl;
                            CidToPidFun.display();
                            HneiOut[hid].display();
                            existSC.display();
                            lpfTmp.display();
                        }
                    }
                }
            }

            for (auto &j: LneiOut)
            {
                lid = j.first;
                for (int k = 0; k < Tree[rank[lid]].vertIn.size(); k++)
                {
                    if (Tree[rank[lid]].vertIn[k].first == Cid)
                    {
                        lpfTmp = CidToPidFun.LPFCatSupport(j.second, lBound, uBound, acc);
                        existSC = Tree[rank[lid]].vertIn[k].second.first;
                        if (!existSC.dominate(lpfTmp, acc))
                        {
                            Tree[rank[lid]].vertIn[k].second = {existSC.LPFMinSupForDec(lpfTmp, acc), 1};
                            SCreIn[lid].insert(Cid);
                            OC.insert(OrderCompp(lid));
                            scProcThisTurn.insert({Cid, lid});
                            if (Cid == debugID1 and ProID == debugID2)
                            {
                                cout << "lid: " << lid << " from Cid: " << Cid << " via ProId: " << ProID << endl;
                                CidToPidFun.display();
                                j.second.display();
                                existSC.display();
                                lpfTmp.display();
                            }
                        }
                        break;
                    }
                }
            }
        }

        for (const auto &Cid: SCreOut[ProID])
        {
            if (ProID == debugID1 and Cid == debugID2)
                cout << "ProId: " << ProID << " to Cid: " << Cid << endl;

            int cidH = Tree[rank[Cid]].height - 1;
            map<int, LPFunction> HneiIn;
            vector<pair<int, LPFunction>> LneiIn;
            for (auto &j: Tree[rank[ProID]].vertIn)
            {
                if (NodeOrder[j.first] > NodeOrder[Cid])
                    HneiIn[j.first] = j.second.first;
                else if (NodeOrder[j.first] < NodeOrder[Cid])
                    LneiIn.emplace_back(j.first, j.second.first);
            }
            LPFunction PidToCidFun;
            for (const auto &j: Tree[rank[ProID]].vertOut)
            {
                // auto &j: Tree[rank[ProID]].vertOut
                if (j.first == Cid)
                {
                    // j.first == Cid, appear once
                    PidToCidFun = j.second.first;
                    break;
                }
            }
            assert(PidToCidFun.vX.size() > 1);
            LPFunction lpfExist = Tree[rank[ProID]].disOut[cidH];
//            if (ProID == debugID1 and Cid == debugID2)
//            {
//                lpfExist.display();
//                PidToCidFun.display();
//            }
            if (!lpfExist.dominate(PidToCidFun, acc))
            {
                LPFunction PidToCidFun2 = PidToCidFun;
                PidToCidFun2.scSupToLbSup(ProID);
                Tree[rank[ProID]].disOut[cidH] = lpfExist.LPFMinSupForDec(PidToCidFun2, acc);
                ProIDdisChaOut = true;
                // DisRe contains the accessories whose distance labeling has changed
                Tree[rank[ProID]].DisReOut.insert(Cid);
            }

            // propogate the updates of tree edges to the ancestors of ProID
            int hid, lid;
            LPFunction existSC, lpfTmp;
            for (int j = 0; j < Tree[rank[Cid]].vertIn.size(); j++)
            {
                hid = Tree[rank[Cid]].vertIn[j].first;
                if (scProcThisTurn.find(hid) != scProcThisTurn.end() and scProcThisTurn[hid] == Cid)
                    continue;

                if (HneiIn.find(hid) != HneiIn.end())
                {
                    // otherwise, e(hid, cid) is not supported by (a,b)
                    // how to concatenate?
                    lpfTmp = HneiIn[hid].LPFCatSupport(PidToCidFun, lBound, uBound, acc);
                    existSC = Tree[rank[Cid]].vertIn[j].second.first;

                    if (!existSC.dominate(lpfTmp, acc))
                    {
                        Tree[rank[Cid]].vertIn[j].second = {existSC.LPFMinSupForDec(lpfTmp, acc), 1};
                        SCreIn[Cid].insert(hid);
                        OC.insert(OrderCompp(Cid));
                        if (ProID == debugID1 and Cid == debugID2)
                        {
                            cout << "Cid: " << Cid << " from hid: " << hid << " via ProId: " << ProID << endl;
                            HneiIn[hid].display();
                            PidToCidFun.display();
                            existSC.display();
                            lpfTmp.display();
                            Tree[rank[Cid]].vertIn[j].second.first.display();
                        }
                    }

                }
            }

            for (auto &j: LneiIn)
            {
                lid = j.first;
                if (scProcThisTurn.find(lid) != scProcThisTurn.end() and scProcThisTurn[lid] == Cid)
                    continue;

                for (int k = 0; k < Tree[rank[lid]].vertOut.size(); k++)
                {
                    if (Tree[rank[lid]].vertOut[k].first == Cid)
                    {
                        existSC = Tree[rank[lid]].vertOut[k].second.first;
                        lpfTmp = j.second.LPFCatSupport(PidToCidFun, lBound, uBound, acc);
                        if (!existSC.dominate(lpfTmp, acc))
                        {
                            Tree[rank[lid]].vertOut[k].second = {existSC.LPFMinSupForDec(lpfTmp, acc), 1};
                            SCreOut[lid].insert(Cid);
                            OC.insert(OrderCompp(lid));
                            if (ProID == debugID1 and Cid == debugID2)
                            {
                                cout << "lid: " << lid << " to Cid: " << Cid << " via ProId: " << ProID << endl;
                                j.second.display();
                                PidToCidFun.display();
                                existSC.display();
                                lpfTmp.display();
                            }
                        }
                        break;
                    }
                }
            }
        }

        if (ProIDdisChaIn or ProIDdisChaOut)
        {
            // if the distance labeling is detected changed
            ProBeginVertexSetNew.clear();
            ProBeginVertexSetNew.reserve(ProBeginVertexSet.size() + 1);
            ProBeginVertexSetNew.emplace_back(ProID);
            int rnew = rank[ProID], r;
            for (int i: ProBeginVertexSet)
            {
                r = rank[i]; // r is the rank of proID from previous iterations

                // if rnew is an ancestor of r, then r is not the source of top-down distance updates
                if (LCAQuery(rnew, r) != rnew)
                {
                    ProBeginVertexSetNew.push_back(i);
                }
            }
            ProBeginVertexSet = ProBeginVertexSetNew;
        }
    }

    cout << "** finish bottom-up refresh" << endl;

    // top-down manner to update distance labels
    Timer clock;
    clock.tick();
    EachNodeDecMainThread(ProBeginVertexSet);
//    for (int ProBeginVertexID: ProBeginVertexSet)
//    {
//        EachNodeDecMain(rank[ProBeginVertexID]);
//    }
    clock.tock();
    cout << "** time for top-down refresh: " << clock.duration().count() << " seconds" << endl;
}

void Graph::EachNodeDecMainThread(const vector<int> &startSources)
{
    vector<int> vRanks;
    vRanks.reserve(startSources.size());
    for (const auto &i: startSources)
    {
        vRanks.emplace_back(rank[i]);
    }

    while (!vRanks.empty() and vRanks.size() < threadNum)
    {
        int sourceR = vRanks.back();
        vRanks.pop_back();

        EachNodeDecMain(sourceR);
        for (int i: Tree[sourceR].vChildren)
            vRanks.emplace_back(i);
    }

    if (vRanks.empty())
        return;

    assert(vRanks.size() >= threadNum);
    vector<pair<int, int>> vpID;
    int step = vRanks.size() / threadNum;
    vpID.reserve(threadNum);

    for (int i = 0; i < threadNum - 1; i++)
    {
        vpID.emplace_back(i * step, (i + 1) * step);
    }
    vpID.emplace_back((threadNum - 1) * step, vRanks.size());

    boost::thread_group threads;
    for (int i = 0; i < threadNum; i++)
    {
        threads.create_thread(
                boost::bind(&Graph::EachNodeDecMainRange,
                            this, vpID[i].first, vpID[i].second, boost::ref(vRanks))
        );
    }
    threads.join_all();
}

void Graph::EachNodeDecMainRange(int begin, int end, const vector<int> &vRanks)
{
    for (int i = begin; i < end; i++)
    {
        EachNodeDecMainRec(vRanks[i]);
    }
}

void Graph::EachNodeDecMainRec(int vRank)
{
    EachNodeDecMain(vRank);
    for (int j: Tree[vRank].vChildren)
        EachNodeDecMainRec(j);
}

void Graph::EachNodeDecMain(int sourceR)
{
    int debugID1 = -1, debugID2 = -1;
    vector<int> line = Tree[sourceR].ancIDs;
    int source = Tree[sourceR].uniqueVertex;
    if (source == debugID1)
        cout << 1;
    for (int k = 0; k < Tree[sourceR].vertIn.size(); k++)
    {
        // b as a support vertex
        int b = Tree[sourceR].vertIn[k].first, bH = Tree[rank[b]].height - 1;
        LPFunction scBS = Tree[sourceR].vertIn[k].second.first; // func: b to sourceR
        scBS.scSupToLbSup(source);
        bool sbLabChanged = Tree[sourceR].DisReIn.find(b) != Tree[sourceR].DisReIn.end();
        for (int i = 0; i < line.size(); i++)
        {
            int a = line[i];
            if (Tree[sourceR].disIn[bH].vSupportContains(source))
            {
                // scBS at least partly supports disBS
                LPFunction lpfExist = Tree[sourceR].disIn[i];
                LPFunction lpfTmp;
                if (i > bH)
                {
                    // a -> b -> source
                    bool cond = Tree[rank[a]].DisReOut.find(b) != Tree[rank[a]].DisReOut.end();
                    if (sbLabChanged or cond)
                    {
                        lpfTmp = Tree[rank[a]].disOut[bH].LPFCatSupport(scBS, lBound, uBound, acc);
                    }
                } else if (i < bH)
                {
                    // a -> b -> source
                    bool cond = Tree[rank[b]].DisReIn.find(a) != Tree[rank[b]].DisReIn.end();
                    if (sbLabChanged or cond)
                    {
                        lpfTmp = Tree[rank[b]].disIn[i].LPFCatSupport(scBS, lBound, uBound, acc);
                    }
                }

                if (lpfTmp.vX.size() > 1 and a == 1 and source == 1)
                {
                    cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
                    lpfExist.display();
                    lpfTmp.display();
                }
                if (lpfTmp.vX.size() > 1 and !lpfExist.dominate(lpfTmp, acc))
                {
                    Tree[sourceR].disIn[i] = lpfExist.LPFMinSupForDec(lpfTmp, acc);
                    Tree[sourceR].DisReIn.insert(a);
                }
            }
        }
    }

    for (int k = 0; k < Tree[sourceR].vertOut.size(); k++)
    {
        // b as a support vertex
        int b = Tree[sourceR].vertOut[k].first, bH = Tree[rank[b]].height - 1;
        LPFunction scSB = Tree[sourceR].vertOut[k].second.first; // func: sourceR to b
        scSB.scSupToLbSup(source);
        bool scLabChanged = Tree[sourceR].DisReOut.find(b) != Tree[sourceR].DisReOut.end();
        for (int i = 0; i < line.size(); i++)
        {
            int a = line[i];
            if (a == debugID2)
                cout << 1;
            if (Tree[sourceR].disOut[bH].vSupportContains(source))
            {
                // scSB at least partly supports disSB
                // shortcut (sourceR, b) is not supported by any other accessory Node of sourceR
                LPFunction lpfExist = Tree[sourceR].disOut[i], lpfTmp;
                if (i > bH)
                {
                    bool cond = Tree[rank[a]].DisReIn.find(b) != Tree[rank[a]].DisReIn.end();
                    if (scLabChanged or cond)
                    {
                        lpfTmp = scSB.LPFCatSupport(Tree[rank[a]].disIn[bH], lBound, uBound, acc);
                    }
                } else if (i < bH)
                {
                    bool cond = Tree[rank[b]].DisReOut.find(a) != Tree[rank[b]].DisReOut.end();
                    if (scLabChanged or cond)
                    {
                        // i contracts after both sourceR and b, so se(b,i) is in Node of b
                        lpfTmp = scSB.LPFCatSupport(Tree[rank[b]].disOut[i], lBound, uBound, acc);
                    }
                }

                if (lpfTmp.vX.size() > 1 and !lpfExist.dominate(lpfTmp, acc))
                {
                    Tree[sourceR].disOut[i] = lpfExist.LPFMinSupForDec(lpfTmp, acc);
                    Tree[sourceR].DisReOut.insert(a);
                }
            }
        }
    }
}

bool Graph::lpfIsSupportBy(
        LPFunction &f1, const LPFunction &halfLPF, const vector<int> &halfLPFChangedPoses, bool isSC) const
{
    // check if the change of halfLPF influences f1, halfLPFChangedPoses records the places where halfLPF is updated
    assert(!halfLPFChangedPoses.empty());
    if (halfLPFChangedPoses.size() == 1 and halfLPFChangedPoses[0] == -1)
    {
        // this case is when halfLPF's old turning points do not change but it is appended with new turning points
        return true;
    }

    if (f1.ID1 == halfLPF.ID1 and f1.ID2 == halfLPF.ID2)
    {
        unordered_map<int, int> f2SupParts;
        int w = NodeOrder[f1.ID1] < NodeOrder[f1.ID2] ? f1.ID1 : f1.ID2;

        bool changed = false;
        for (int i = 0; i < f1.vSupportPre.size(); i++)
        {
            if (f1.vSupportPre[i].find(w) == f1.vSupportPre[i].end())
                continue;

            for (const int &pos: f1.vSupportPre[i].at(w))
            {
                if (find(halfLPFChangedPoses.begin(),
                         halfLPFChangedPoses.end(), pos) != halfLPFChangedPoses.end())
                {
                    int &cnt = f1.cntOfEachInt[i];
                    cnt -= 1;
                    if (cnt <= 0)
                        changed = true;
                }
            }
        }
        return changed;
    }

    // u -> w -> v
    int u = f1.ID1, w = -1, v = f1.ID2;
    bool backHalf;
    if (halfLPF.ID1 == u)
    {
        w = halfLPF.ID2;
        backHalf = false;
    } else
    {
        if (halfLPF.ID2 != v)
        {
            f1.display();
            halfLPF.display();
            assert(false);
        }
        w = halfLPF.ID1;
        backHalf = true;
    }

    if (!f1.vSupportContains(w))
        return false;

    CatSupRec supportInfo;
    if (isSC)
    {
        if (intermediateSCs[u].find(w) != intermediateSCs[u].end() and
            intermediateSCs[u].at(w).find(v) != intermediateSCs[u].at(w).end())
            supportInfo = intermediateSCs[u].at(w).find(v)->second;
        else
            return false;
    } else
    {
        if (NodeOrder[u] < NodeOrder[v])
        {
            if (intermediateLBsIn[u].find(w) != intermediateLBsIn[u].end()
                and intermediateLBsIn[u].at(w).find(v) != intermediateLBsIn[u].at(w).end())
                supportInfo = intermediateLBsIn[u].at(w).find(v)->second;
            else
                return false;
        } else
        {
            if (intermediateLBsOut[u].find(w) != intermediateLBsOut[u].end()
                and intermediateLBsOut[u].at(w).find(v) != intermediateLBsOut[u].at(w).end())
                supportInfo = intermediateLBsOut[u].at(w).find(v)->second;
            else
                return false;
        }
    }
    if (supportInfo.vX1.empty())
        return false;

    // map changed parts of halfLPF to intermediateSCs[w][u][v]
    vector<int> cctLPFvX;
    if (backHalf)
        cctLPFvX = supportInfo.vX2;
    else
        cctLPFvX = supportInfo.vX1;

    unordered_set<int> changedPosInConcatLPF;
    int itv = 0;
//    assert(!cctLPFvX.empty() and halfLPF.vX.front() <= cctLPFvX.front());
    for (const auto &changedPosInHalfLPF: halfLPFChangedPoses)
    {
        // changedPosInConcatLPF is the changed part of intermediate[w][u][v] (i.e., suppInfo.first)
        int t1 = halfLPF.vX[changedPosInHalfLPF], t2 = halfLPF.vX[changedPosInHalfLPF + 1];
        while (itv + 1 < cctLPFvX.size() and cctLPFvX[itv + 1] < t1)
            itv++;

        while (itv < cctLPFvX.size() and cctLPFvX[itv] < t2)
        {
            changedPosInConcatLPF.insert(itv);
            itv++;
        }
    }

    bool changed = false;
    for (int i = 0; i < f1.vSupportPre.size(); i++)
    {
        if (f1.vSupportPre[i].find(w) == f1.vSupportPre[i].end())
            continue;
        for (const int &pos: f1.vSupportPre[i].at(w))
        {
            // cctLPF[pos] supports f1[i]
            if (changedPosInConcatLPF.find(pos) != changedPosInConcatLPF.end())
            {
                int &cnt = f1.cntOfEachInt[i];
                cnt -= 1;
                if (cnt <= 0)
                    changed = true;
            }
        }
    }
    return changed;
}

void Graph::CHUpdate(vector<pair<int, int>> &wBatch)
{
    cout << "Index maintenance for " << wBatch.size() << " increase updates" << endl;
    unordered_map < pair<int, int>, LPFunction, boost::hash<pair<int, int>> > OCdisOld;

    // NodeOrderss.clear();
    NodeOrderss.assign(NodeOrder.begin(), NodeOrder.end());
    vector<set<int>> SCreIn, SCreOut; //SCre.clear();
    set<int> ss; //ss.clear();
    SCreIn.assign(nodeNum, ss), SCreOut.assign(nodeNum, ss); //{vertexID, set<int>}
    set < OrderCompp > OC;
    OC.clear();//vertexID in decreasing Node order

    int debugID1 = 1, debugID2 = 1;
    std::random_device rd; // obtain a random number from hardware
    std::mt19937 gen(rd()); // seed the generator
    std::uniform_int_distribution<> distr(1, 10);

    for (const auto &k: wBatch)
    {
        int a = k.first;
        int b = k.second;
        if (a == debugID1 and b == debugID2)
            cout << 1;
        int lId = NodeOrder[a] < NodeOrder[b] ? a : b;
//        int vX = uBound;
        int newW, oldW, pos, vX;
//        newW = k.second.second;
//        double frac = k.second.second;
        bool outOfBound = true;
        for (const auto &i: adjEdge[a])
        {
            if (i.first == b)
            {
                LPFunction &lpf = vEdge[i.second].lpf;

                outOfBound = false;
                if (distr(gen) <= 7)
                {
                    vX = lpf.vX.back();
                    oldW = lpf.vY.back();
                } else
                {
                    vX = lpf.vX.back() - deltaT;
                    oldW = lpf.getY(vX);
                }

                if (distr(gen) <= 5)
                    newW = ceil((double) oldW * 1.1);
                else
                    newW = ceil((double) oldW * 0.9);

                vector<int>::const_iterator low;
                low = lower_bound(lpf.vX.begin(), lpf.vX.end(), vX);
                pos = (int) (low - lpf.vX.begin());

                if (lpf.vX[pos] == vX)
                {
                    lpf.vY[pos] = newW;
                } else if (lpf.vX[pos] > vX)
                {
                    lpf.vX.insert(low, vX);
                    lpf.vY.insert(lpf.vY.begin() + pos, newW);
                    assert(lpf.vX[pos] == vX);
                }

                vector<int> newX = lpf.vX, newY = lpf.vY;
                lpf.setValue(newX, newY, -1, acc);
                for (int j = 0; j < lpf.vX.size() - 1; j++)
                {
                    lpf.vSupportPre.push_back({{lId, {j}}});
                }

//                low = lower_bound(lpf.vX.begin(), lpf.vX.end(), vX);
//                pos = (int) (low - lpf.vX.begin());
//                lpf.lastItvSubToChange = pos > 0 ? pos - 1 : 0;

                CatSupRec catSupRec;
                catSupRec.vX1 = lpf.vX;
                catSupRec.vX2 = lpf.vX;
                catSupRec.vY = lpf.vY;

                intermediateSCs[a][a][b] = catSupRec;
                break;
            }
        }

        if (outOfBound)
            continue;

        if (lId == a)
        {
            for (const auto &j: Tree[rank[a]].vertOut)
            {
                if (j.first == b and j.second.first.vSupportContains(a))
                {
                    OCdisOld[{a, b}] = j.second.first;
                    SCreOut[a].insert(b);
                    OC.insert(OrderCompp(a));
                    break;
                }
            }
        } else
        {
            for (const auto &j: Tree[rank[b]].vertIn)
            {
                if (j.first == a and j.second.first.vSupportContains(b))
                {
                    OCdisOld[{a, b}] = j.second.first;
                    SCreIn[b].insert(a);
                    OC.insert(OrderCompp(b));
                    break;
                }
            }
        }
    }
    cout << "** finish sc refresh" << endl;

    int ProID;
    bool influence;

    while (!OC.empty())
    {
        ProID = (*OC.begin()).x;
        if (ProID == -1)
            cout << 1;

        OC.erase(OC.begin());
        influence = false;

//        cout << "ProID: " << ProID << " " << Tree[rank[ProID]].disIn.size() << endl;
        for (const auto &Cid: SCreIn[ProID])
        {
            if (Cid == debugID1 and ProID == debugID2)
                cout << 1;
            int CidH = Tree[rank[Cid]].height - 1;
            /* check if sc(cid,pid) is really influenced, if so, update it with the new function */
            LPFunction scCPNew;
            for (const auto &i: adjEdge[Cid])
            {
                if (i.first == ProID)
                {
                    scCPNew = vEdge[i.second].lpf;
                    break;
                }
            }

            for (auto &vs: SCconNodesMT[Cid][ProID])
            {
                LPFunction CidToVS, vsToProId;
                for (const auto &i: Tree[rank[vs]].vertIn)
                {
                    if (i.first == Cid)
                    {
                        CidToVS = i.second.first;
                        break;
                    }
                }
                for (const auto &i: Tree[rank[vs]].vertOut)
                {
                    if (i.first == ProID)
                    {
                        vsToProId = i.second.first;
                        break;
                    }
                }

                LPFunction CPTmp;
                bool empty = scCPNew.vX.size() <= 1;
                CatSupRec &catSupRec = intermediateSCs[Cid][vs][ProID];
                if (empty)
                {
                    CPTmp = CidToVS.LPFCatSupport(vsToProId, lBound, uBound, catSupRec, acc);
                    scCPNew = CPTmp;
                } else if (CidToVS.minY + vsToProId.minY <= scCPNew.maxY)
                {
                    CPTmp = CidToVS.LPFCatSupport(vsToProId, lBound, uBound, catSupRec, acc);
                }

                if (!empty and CPTmp.vX.size() > 1 and !scCPNew.dominate(CPTmp, acc))
                {
                    scCPNew = scCPNew.LPFMinSupByPlaneSweep(CPTmp, acc);
                }
            }

            LPFunction scCPOld = OCdisOld[{Cid, ProID}];

            vector<int> changedPos;
            if (scCPOld.equal(scCPNew, changedPos))
                continue;

            for (auto &j: Tree[rank[ProID]].vertIn)
            {
                if (j.first == Cid)
                {
                    assert(scCPNew.ID1 == Cid and scCPNew.ID2 == ProID);
                    j.second = {scCPNew, 1};
                    break;
                }
            }

            // HneiOut: <accessories contract after Cid, lpf from ProId to hid
            unordered_set<int> HneiOut;
            // LneiOut: accessories contract after ProId and before Cid, lpf from ProId to lid
            vector<int> LneiOut;
            for (auto &j: Tree[rank[ProID]].vertOut)
            {
                // j.first is the accessory of ProID
                if (NodeOrder[j.first] > NodeOrder[Cid])
                    HneiOut.insert(j.first);
                else if (NodeOrder[j.first] < NodeOrder[Cid])
                    LneiOut.emplace_back(j.first);
            }

            // in theory, scCPOld should dominate scCPNew, but due to rounding, domination may not hold
            bool newSCDec = !scCPOld.dominate(scCPNew, acc); // assert(newSCDec == false);

            for (auto &i: Tree[rank[Cid]].vertOut)
            {
                int hid = i.first;
                if (HneiOut.find(hid) != HneiOut.end())
                {
                    bool changed =
                            newSCDec or i.second.second <= 0 or
                            lpfIsSupportBy(i.second.first, scCPOld, changedPos, true);
                    if (i.second.second > 0 and changed)
                    {
                        OCdisOld[{Cid, hid}] = i.second.first;
                        SCreOut[Cid].insert(hid);
                        OC.insert(OrderCompp(Cid));
                        i.second.second = 0;
                    }
                }
            }

            for (auto &i: LneiOut)
            {
                int lid = i;
                for (auto &j: Tree[rank[lid]].vertIn)
                {
                    if (j.first == Cid)
                    {
                        bool changed =
                                newSCDec or j.second.second <= 0 or
                                lpfIsSupportBy(j.second.first, scCPOld, changedPos, true);
                        if (j.second.second > 0 and changed)
                        {
                            OCdisOld[{Cid, lid}] = j.second.first;
                            SCreIn[lid].insert(Cid);
                            OC.insert(OrderCompp(lid));
                            j.second.second = 0;
                        }
                        break;
                    }
                }
            }
        }

        for (const auto &Cid: SCreOut[ProID])
        {
            if (Cid == debugID2 and ProID == debugID1)
                cout << 1;
            int CidH = Tree[rank[Cid]].height - 1;

            /* check if sc(pid,cid) is really influenced, if so, update it with the new function */
            LPFunction scPCNew;
            for (const auto &i: adjEdge[ProID])
            {
                if (i.first == Cid)
                {
                    scPCNew = vEdge[i.second].lpf;
                    break;
                }
            }

            for (auto &vs: SCconNodesMT[ProID][Cid])
            {
                LPFunction ProidToVS, vsToCid;
                for (const auto &i: Tree[rank[vs]].vertIn)
                {
                    if (i.first == ProID)
                    {
                        ProidToVS = i.second.first;
                        break;
                    }
                }

                for (const auto &i: Tree[rank[vs]].vertOut)
                {
                    if (i.first == Cid)
                    {
                        vsToCid = i.second.first;
                        break;
                    }
                }

                assert(ProidToVS.ID1 == ProID and vsToCid.ID2 == Cid);
                LPFunction PCTmp;
                bool empty = scPCNew.vX.size() <= 1;
                CatSupRec &catSupRec = intermediateSCs[ProID][vs][Cid];
                if (empty)
                {
                    PCTmp = ProidToVS.LPFCatSupport(vsToCid, lBound, uBound, catSupRec, acc);
                    scPCNew = PCTmp;
                } else if (ProidToVS.minY + vsToCid.minY <= scPCNew.maxY)
                {
                    PCTmp = ProidToVS.LPFCatSupport(vsToCid, lBound, uBound, catSupRec, acc);
                }

                if (!empty and PCTmp.vX.size() > 1 and !scPCNew.dominate(PCTmp, acc))
                {
                    scPCNew = scPCNew.LPFMinSupByPlaneSweep(PCTmp, acc);
                }
            }

            LPFunction scPCOld = OCdisOld[{ProID, Cid}];
            vector<int> changedPos;
            if (scPCOld.equal(scPCNew, changedPos))
                continue;
            for (auto &j: Tree[rank[ProID]].vertOut)
            {
                if (j.first == Cid)
                {
                    assert(scPCNew.ID1 == ProID and scPCNew.ID2 == Cid);
                    j.second = {scPCNew, 1};
                    break;
                }
            }

            unordered_set<int> HneiIn;
            vector<int> LneiIn;
            for (auto &i: Tree[rank[ProID]].vertIn)
            {
                // j.first is the accessory of ProID
                if (NodeOrder[i.first] > NodeOrder[Cid])
                    HneiIn.insert(i.first);
                else if (NodeOrder[i.first] < NodeOrder[Cid])
                    LneiIn.emplace_back(i.first);
            }

            // in theory, scCPOld should dominate scCPNew, but due to rounding, domination may not hold
            bool newSCDec = !scPCOld.dominate(scPCNew, acc); // assert(newSCDec == false);

            for (auto &i: Tree[rank[Cid]].vertIn)
            {
                int hid = i.first;
                if (HneiIn.find(hid) != HneiIn.end())
                {
                    bool changed =
                            newSCDec or i.second.second <= 0 or
                            lpfIsSupportBy(i.second.first, scPCOld, changedPos, true);
                    if (i.second.second > 0 and changed)
                    {
                        OCdisOld[{hid, Cid}] = i.second.first;
                        SCreIn[Cid].insert(hid);
                        OC.insert(OrderCompp(Cid));
                        i.second.second = 0;
                    }
                }
            }

            for (auto &lid: LneiIn)
            {
                for (auto &i: Tree[rank[lid]].vertOut)
                {
                    if (i.first == Cid)
                    {
                        bool changed =
                                newSCDec or i.second.second <= 0 or
                                lpfIsSupportBy(i.second.first, scPCOld, changedPos, true);
                        if (i.second.second > 0 and changed)
                        {
                            OCdisOld[{lid, Cid}] = i.second.first;
                            SCreOut[lid].insert(Cid);
                            OC.insert(OrderCompp(lid));
                            i.second.second = 0;
                        }
                        break;
                    }
                }
            }
        }
    }

    cout << "** finish bottom-up refresh" << endl;
}

void Graph::H2HIncBatDiGraph(vector<pair<int, int>> &wBatch)
{
    cout << "Index maintenance for " << wBatch.size() << " increase updates" << endl;
    unordered_map < pair<int, int>, LPFunction, boost::hash<pair<int, int>> > OCdisOld;

    // NodeOrderss.clear();
    NodeOrderss.assign(NodeOrder.begin(), NodeOrder.end());
    vector<set<int>> SCreIn, SCreOut; //SCre.clear();
    set<int> ss; //ss.clear();
    SCreIn.assign(nodeNum, ss), SCreOut.assign(nodeNum, ss); //{vertexID, set<int>}
    set < OrderCompp > OC;
    OC.clear();//vertexID in decreasing Node order

    int debugID1 = 1, debugID2 = 1;
    std::random_device rd; // obtain a random number from hardware
    std::mt19937 gen(rd()); // seed the generator
    std::uniform_int_distribution<> distr(1, 10);

    for (const auto &k: wBatch)
    {
        int a = k.first;
        int b = k.second;
        if (a == debugID1 and b == debugID2)
            cout << 1;
        int lId = NodeOrder[a] < NodeOrder[b] ? a : b;

        int newW, oldW, pos, updX;
        bool outOfBound = true;
        for (const auto &i: adjEdge[a])
        {
            if (i.first == b)
            {
                LPFunction &lpf = vEdge[i.second].lpf;

                outOfBound = false;
                if (distr(gen) <= 7)
                {
                    updX = lpf.vX.back();
                    oldW = lpf.vY.back();
                } else
                {
                    updX = lpf.vX.back() - deltaT;
                    oldW = lpf.getY(updX);
                }

                if (distr(gen) <= 5)
                    newW = ceil((double) oldW * 1.1);
                else
                    newW = ceil((double) oldW * 0.9);

                vector<int>::const_iterator low;
                low = lower_bound(lpf.vX.begin(), lpf.vX.end(), updX);
                pos = (int) (low - lpf.vX.begin());

                if (lpf.vX[pos] == updX)
                {
                    lpf.vY[pos] = newW;
                } else if (lpf.vX[pos] > updX)
                {
                    lpf.vX.insert(low, updX);
                    lpf.vY.insert(lpf.vY.begin() + pos, newW);
                    assert(lpf.vX[pos] == updX);
                }

                vector<int> newX = lpf.vX, newY = lpf.vY;
                lpf.setValue(newX, newY, -1, acc);
                for (int j = 0; j < lpf.vX.size() - 1; j++)
                {
                    lpf.vSupportPre.push_back({{lId, {j}}});
                }

//                low = lower_bound(lpf.vX.begin(), lpf.vX.end(), vX);
//                pos = (int) (low - lpf.vX.begin());
//                lpf.lastItvSubToChange = pos > 0 ? pos - 1 : 0;

                CatSupRec catSupRec;
                catSupRec.vX1 = lpf.vX;
                catSupRec.vX2 = lpf.vX;
                catSupRec.vY = lpf.vY;

                intermediateSCs[a][a][b] = catSupRec;
                break;
            }
        }

        if (outOfBound)
            continue;

        if (lId == a)
        {
            for (const auto &j: Tree[rank[a]].vertOut)
            {
                if (j.first == b and j.second.first.vSupportContains(a))
                {
                    OCdisOld[{a, b}] = j.second.first;
                    SCreOut[a].insert(b);
                    OC.insert(OrderCompp(a));
                    break;
                }
            }
        } else
        {
            for (const auto &j: Tree[rank[b]].vertIn)
            {
                if (j.first == a and j.second.first.vSupportContains(b))
                {
                    OCdisOld[{a, b}] = j.second.first;
                    SCreIn[b].insert(a);
                    OC.insert(OrderCompp(b));
                    break;
                }
            }
        }
    }
    cout << "** finish sc refresh" << endl;

    vector<int> ProBeginVertexSet;
    vector<int> ProBeginVertexSetNew;
    int ProID;
    bool influence;

    while (!OC.empty())
    {
        ProID = (*OC.begin()).x;
        if (ProID == -1)
            cout << 1;

        OC.erase(OC.begin());
        influence = false;

//        cout << "ProID: " << ProID << " " << Tree[rank[ProID]].disIn.size() << endl;
        for (const auto &Cid: SCreIn[ProID])
        {
            if (Cid == debugID1 and ProID == debugID2)
                cout << 1;
            int CidH = Tree[rank[Cid]].height - 1;
            /* check if sc(cid,pid) is really influenced, if so, update it with the new function */
            LPFunction scCPNew;
            for (const auto &i: adjEdge[Cid])
            {
                if (i.first == ProID)
                {
                    scCPNew = vEdge[i.second].lpf;
                    break;
                }
            }

            for (auto &vs: SCconNodesMT[Cid][ProID])
            {
                LPFunction CidToVS, vsToProId;
                for (const auto &i: Tree[rank[vs]].vertIn)
                {
                    if (i.first == Cid)
                    {
                        CidToVS = i.second.first;
                        break;
                    }
                }
                for (const auto &i: Tree[rank[vs]].vertOut)
                {
                    if (i.first == ProID)
                    {
                        vsToProId = i.second.first;
                        break;
                    }
                }

                LPFunction CPTmp;
                bool empty = scCPNew.vX.size() <= 1;
                CatSupRec &catSupRec = intermediateSCs[Cid][vs][ProID];
                if (empty)
                {
                    CPTmp = CidToVS.LPFCatSupport(vsToProId, lBound, uBound, catSupRec, acc);
                    scCPNew = CPTmp;
                } else if (CidToVS.minY + vsToProId.minY <= scCPNew.maxY)
                {
                    CPTmp = CidToVS.LPFCatSupport(vsToProId, lBound, uBound, catSupRec, acc);
                }

                if (!empty and CPTmp.vX.size() > 1 and !scCPNew.dominate(CPTmp, acc))
                {
                    scCPNew = scCPNew.LPFMinSupByPlaneSweep(CPTmp, acc);
                }
            }

            LPFunction scCPOld = OCdisOld[{Cid, ProID}];

            vector<int> changedPos;
            if (scCPOld.equal(scCPNew, changedPos))
                continue;

            for (auto &j: Tree[rank[ProID]].vertIn)
            {
                if (j.first == Cid)
                {
                    assert(scCPNew.ID1 == Cid and scCPNew.ID2 == ProID);
                    j.second = {scCPNew, 1};
                    break;
                }
            }

            // HneiOut: <accessories contract after Cid, lpf from ProId to hid
            unordered_set<int> HneiOut;
            // LneiOut: accessories contract after ProId and before Cid, lpf from ProId to lid
            vector<int> LneiOut;
            for (auto &j: Tree[rank[ProID]].vertOut)
            {
                // j.first is the accessory of ProID
                if (NodeOrder[j.first] > NodeOrder[Cid])
                    HneiOut.insert(j.first);
                else if (NodeOrder[j.first] < NodeOrder[Cid])
                    LneiOut.emplace_back(j.first);
            }

            // in theory, scCPOld should dominate scCPNew, but due to rounding, domination may not hold
            bool newSCDec = !scCPOld.dominate(scCPNew, acc); // assert(newSCDec == false);

            for (auto &i: Tree[rank[Cid]].vertOut)
            {
                int hid = i.first;
                if (HneiOut.find(hid) != HneiOut.end())
                {
                    bool changed =
                            newSCDec or i.second.second <= 0 or
                            lpfIsSupportBy(i.second.first, scCPOld, changedPos, true);
                    if (i.second.second > 0 and changed)
                    {
                        OCdisOld[{Cid, hid}] = i.second.first;
                        SCreOut[Cid].insert(hid);
                        OC.insert(OrderCompp(Cid));
                        i.second.second = 0;
                    }
                }
            }

            for (auto &i: LneiOut)
            {
                int lid = i;
                for (auto &j: Tree[rank[lid]].vertIn)
                {
                    if (j.first == Cid)
                    {
                        bool changed =
                                newSCDec or j.second.second <= 0 or
                                lpfIsSupportBy(j.second.first, scCPOld, changedPos, true);
                        if (j.second.second > 0 and changed)
                        {
                            OCdisOld[{Cid, lid}] = j.second.first;
                            SCreIn[lid].insert(Cid);
                            OC.insert(OrderCompp(lid));
                            j.second.second = 0;
                        }
                        break;
                    }
                }
            }

            int pAncs = (int) Tree[rank[ProID]].ancIDs.size();
            for (int i = 0; i < pAncs; i++)
            {
                // i -> cid -> pid
                LPFunction &disIP = Tree[rank[ProID]].disIn[i];
                if (newSCDec)
                {
                    if (i == CidH)
                        influence = true;
                    Tree[rank[ProID]].cntIn[i] = 0;
                } else if (lpfIsSupportBy(disIP, scCPOld, changedPos, false))
                {
                    Tree[rank[ProID]].cntIn[i] = 0;
                    if (i == CidH)
                        influence = true;
                }
            }

//            scCPNew.scSupToLbSup(ProID);
            CatSupRec catSupRec;
            catSupRec.vX1 = scCPNew.vX;
            catSupRec.vX2 = scCPNew.vX;
            catSupRec.vY = scCPNew.vY;
            intermediateLBsIn[ProID][ProID][Cid] = catSupRec;
        }

        for (const auto &Cid: SCreOut[ProID])
        {
            if (Cid == debugID2 and ProID == debugID1)
                cout << 1;
            int CidH = Tree[rank[Cid]].height - 1;

            /* check if sc(pid,cid) is really influenced, if so, update it with the new function */
            LPFunction scPCNew;
            for (const auto &i: adjEdge[ProID])
            {
                if (i.first == Cid)
                {
                    scPCNew = vEdge[i.second].lpf;
                    break;
                }
            }

            for (auto &vs: SCconNodesMT[ProID][Cid])
            {
                LPFunction ProidToVS, vsToCid;
                for (const auto &i: Tree[rank[vs]].vertIn)
                {
                    if (i.first == ProID)
                    {
                        ProidToVS = i.second.first;
                        break;
                    }
                }

                for (const auto &i: Tree[rank[vs]].vertOut)
                {
                    if (i.first == Cid)
                    {
                        vsToCid = i.second.first;
                        break;
                    }
                }

                assert(ProidToVS.ID1 == ProID and vsToCid.ID2 == Cid);
                LPFunction PCTmp;
                bool empty = scPCNew.vX.size() <= 1;
                CatSupRec &catSupRec = intermediateSCs[ProID][vs][Cid];
                if (empty)
                {
                    PCTmp = ProidToVS.LPFCatSupport(vsToCid, lBound, uBound, catSupRec, acc);
                    scPCNew = PCTmp;
                } else if (ProidToVS.minY + vsToCid.minY <= scPCNew.maxY)
                {
                    PCTmp = ProidToVS.LPFCatSupport(vsToCid, lBound, uBound, catSupRec, acc);
                }

                if (!empty and PCTmp.vX.size() > 1 and !scPCNew.dominate(PCTmp, acc))
                {
                    scPCNew = scPCNew.LPFMinSupByPlaneSweep(PCTmp, acc);
                }
            }

            LPFunction scPCOld = OCdisOld[{ProID, Cid}];
            vector<int> changedPos;
            if (scPCOld.equal(scPCNew, changedPos))
                continue;
            for (auto &j: Tree[rank[ProID]].vertOut)
            {
                if (j.first == Cid)
                {
                    assert(scPCNew.ID1 == ProID and scPCNew.ID2 == Cid);
                    j.second = {scPCNew, 1};
                    break;
                }
            }

            unordered_set<int> HneiIn;
            vector<int> LneiIn;
            for (auto &i: Tree[rank[ProID]].vertIn)
            {
                // j.first is the accessory of ProID
                if (NodeOrder[i.first] > NodeOrder[Cid])
                    HneiIn.insert(i.first);
                else if (NodeOrder[i.first] < NodeOrder[Cid])
                    LneiIn.emplace_back(i.first);
            }

            // in theory, scCPOld should dominate scCPNew, but due to rounding, domination may not hold
            bool newSCDec = !scPCOld.dominate(scPCNew, acc); // assert(newSCDec == false);

            for (auto &i: Tree[rank[Cid]].vertIn)
            {
                int hid = i.first;
                if (HneiIn.find(hid) != HneiIn.end())
                {
                    bool changed =
                            newSCDec or i.second.second <= 0 or
                            lpfIsSupportBy(i.second.first, scPCOld, changedPos, true);
                    if (i.second.second > 0 and changed)
                    {
                        OCdisOld[{hid, Cid}] = i.second.first;
                        SCreIn[Cid].insert(hid);
                        OC.insert(OrderCompp(Cid));
                        i.second.second = 0;
                    }
                }
            }

            for (auto &lid: LneiIn)
            {
                for (auto &i: Tree[rank[lid]].vertOut)
                {
                    if (i.first == Cid)
                    {
                        bool changed = newSCDec or i.second.second <= 0 or
                                       lpfIsSupportBy(i.second.first, scPCOld, changedPos, true);
                        if (i.second.second > 0 and changed)
                        {
                            OCdisOld[{lid, Cid}] = i.second.first;
                            SCreOut[lid].insert(Cid);
                            OC.insert(OrderCompp(lid));
                            i.second.second = 0;
                        }
                        break;
                    }
                }
            }

            influence = true;
            int pAncs = (int) Tree[rank[ProID]].ancIDs.size();
            for (int i = 0; i < pAncs; i++)
            {
                LPFunction &disPI = Tree[rank[ProID]].disOut[i], disCI;
                if (newSCDec)
                {
                    if (i == CidH)
                        influence = true;
                    Tree[rank[ProID]].cntOut[i] = 0;
                } else if (lpfIsSupportBy(disPI, scPCOld, changedPos, false))
                {
                    Tree[rank[ProID]].cntOut[i] = 0;
                    if (i == CidH)
                        influence = true;
                }
            }

//            scPCNew.scSupToLbSup(ProID);
            CatSupRec catSupRec;
            catSupRec.vX1 = scPCNew.vX;
            catSupRec.vX2 = scPCNew.vX;
            catSupRec.vY = scPCNew.vY;
            intermediateLBsOut[ProID][ProID][Cid] = catSupRec;
        }

        if (influence)
        {
            ProBeginVertexSetNew.clear();
            ProBeginVertexSetNew.reserve(ProBeginVertexSet.size() + 1);
            ProBeginVertexSetNew.push_back(ProID);
            int rnew = rank[ProID], r;
            for (int i: ProBeginVertexSet)
            {
                r = rank[i];
                if (LCAQuery(rnew, r) != rnew)
                {
                    ProBeginVertexSetNew.push_back(i);
                }
            }
            ProBeginVertexSet = ProBeginVertexSetNew;
        }
    }

    cout << "** finish bottom-up refresh" << endl;

    Timer clock;
    clock.tick();
    EachNodeIncMainThread(ProBeginVertexSet);
    clock.tock();
    cout << "** time for top-down refresh: " << clock.duration().count() << " seconds" << endl;
}

void Graph::EachNodeIncMainThread(const vector<int> &sources)
{
    vector<int> uRanks;
    uRanks.reserve(sources.size());
    for (int ProBeginVertexID: sources)
    {
        uRanks.push_back(rank[ProBeginVertexID]);
    }

    Timer clock;
    clock.tick();
    while (!uRanks.empty() and uRanks.size() < threadNum)
    {
        boost::thread_group threads;
        for (const auto &u: uRanks)
        {
            threads.create_thread(boost::bind(&Graph::EachNodeIncMain, this, u));
        }
        threads.join_all();

        vector<int> newURanks;
        for (const auto &ur: uRanks)
        {
            for (int j: Tree[ur].vChildren)
            {
                newURanks.push_back(j);
            }
        }
        uRanks = newURanks;
    }
    clock.tock();
//    cout << "** time for top-down refresh 1: " << clock.duration().count() << " ms" << endl;

    clock.tick();
    boost::thread_group threads;
    for (const auto &ur: uRanks)
    {
        threads.create_thread(boost::bind(&Graph::EachNodeIncMainRec, this, ur));
    }
    threads.join_all();
    clock.tock();
//    cout << "** time for top-down refresh 2: " << clock.duration().count() << " ms" << endl;
}

void Graph::EachNodeIncMainRec(int uRank)
{
    EachNodeIncMain(uRank);
    for (int j: Tree[uRank].vChildren)
        EachNodeIncMainRec(j);
}

void Graph::EachNodeIncMain(int uRank)
{
    int u = Tree[uRank].uniqueVertex;
    int uH = Tree[uRank].height - 1;

    vector<int> line = Tree[uRank].ancIDs;

    for (int aH = 0; aH < line.size(); aH++)
    {
        int a = line[aH];
        LPFunction disAUOld = Tree[uRank].disIn[aH];

        if (Tree[uRank].cntIn[aH] > 0) continue;
        //firstly, calculate the actual distance
        int b, bH;
        LPFunction disAUMin, scBU;
        for (int j = 0; j < Tree[uRank].vertIn.size(); j++)
        {
            // a -> b -> u
            b = Tree[uRank].vertIn[j].first;
            scBU = Tree[uRank].vertIn[j].second.first;
            scBU.scSupToLbSup(u);
            bH = Tree[rank[b]].height - 1;
            bool empty = disAUMin.vX.size() < 2;
            LPFunction disAUTmp;
            CatSupRec &catSupRec = intermediateLBsIn[u][b][a];
            if (bH < aH)
            {
                LPFunction disAB = Tree[rank[a]].disOut[bH];
                if (empty)
                {
                    disAUTmp = disAB.LPFCatSupport(scBU, lBound, uBound, catSupRec, acc);
                    disAUMin = disAUTmp;
                } else if (disAB.minY + scBU.minY <= disAUMin.maxY)
                {
                    disAUTmp = disAB.LPFCatSupport(scBU, lBound, uBound, catSupRec, acc);
                }
            } else if (bH == aH)
            {
                // b == a
                if (empty)
                {
                    disAUTmp = scBU;
                    disAUMin = disAUTmp;
                } else if (scBU.minY <= disAUMin.maxY or scBU.vX.back() > disAUMin.vX.back())
                    disAUTmp = scBU;
            } else
            {
                LPFunction disAB = Tree[rank[b]].disIn[aH];
                if (empty)
                {
                    disAUTmp = disAB.LPFCatSupport(scBU, lBound, uBound, catSupRec, acc);
                    disAUMin = disAUTmp;
                } else if (disAB.minY + scBU.minY <= disAUMin.maxY)
                {
                    disAUTmp = disAB.LPFCatSupport(scBU, lBound, uBound, catSupRec, acc);
                }
            }

            if (disAUTmp.vX.size() > 1 and u == 1 and a == 1)
            {
                cout << "~~~~~~~~~~~~~~~~~~~" << endl;
                disAUMin.display();
                disAUTmp.display();
            }
            if (!empty and disAUTmp.vX.size() > 1 and !disAUMin.dominate(disAUTmp, acc))
            {
                disAUMin = disAUMin.LPFMinSupByPlaneSweep(disAUTmp, acc);
            }
            assert(disAUMin.ID1 == a and disAUMin.ID2 == u);
        }

        vector<int> changedPos;
        bool updated = disAUMin.vX.size() > 1 and !disAUOld.equal(disAUMin, changedPos);
        if (!updated)
            continue;

        Tree[uRank].disIn[aH] = disAUMin;
        Tree[uRank].cntIn[aH] = 1;

        bool newLabDec = !disAUOld.dominate(disAUMin, acc);
        // firstly, check which dis can be infected
        for (int vRank: VidtoTNid[u])
        {
            // lb(a, u) + sc(u, v)
            int &cnt = Tree[vRank].cntIn[aH];
            if (cnt == 0)
                continue;
            LPFunction &disAV = Tree[vRank].disIn[aH];
            bool changed = newLabDec or lpfIsSupportBy(disAV, disAUOld, changedPos, false);
            if (changed)
                cnt = 0;
        }

        for (int vRank: VidtoTNid[a])
        {
            // sc(v, a) + lb(a, u)
            if (Tree[vRank].height <= Tree[uRank].height or Tree[vRank].ancIDs[uH] != u)
                continue;
            int &cnt = Tree[vRank].cntOut[uH];
            if (cnt == 0)
                continue;
            LPFunction &disVU = Tree[vRank].disOut[uH];
            bool changed = newLabDec or lpfIsSupportBy(disVU, disAUOld, changedPos, false);
            if (changed)
                cnt = 0;
        }
    }

    for (int aH = 0; aH < line.size(); aH++)
    {
        int a = line[aH];
        LPFunction disUAOld = Tree[uRank].disOut[aH];

        if (Tree[uRank].cntOut[aH] > 0) continue;
        //firstly, calculate the actual distance
        int b, bH;

        LPFunction disUAMin, scUB, scUA;
        for (int j = 0; j < Tree[uRank].vertOut.size(); j++)
        {
            // a <- b <- u
            b = Tree[uRank].vertOut[j].first;
            scUB = Tree[uRank].vertOut[j].second.first;
            scUB.scSupToLbSup(u);
            bH = Tree[rank[b]].height - 1;
            bool empty = disUAMin.vX.size() < 2;
            LPFunction disUATmp;
            CatSupRec &catSupRec = intermediateLBsOut[u][b][a];
            if (bH < aH)
            {
                LPFunction disBA = Tree[rank[a]].disIn[bH];
                if (empty)
                {
                    disUATmp = scUB.LPFCatSupport(disBA, lBound, uBound, catSupRec, acc);
                    disUAMin = disUATmp;
                } else
                {
                    if (disBA.minY + scUB.minY <= disUAMin.maxY)
                    {
                        disUATmp = scUB.LPFCatSupport(disBA, lBound, uBound, catSupRec, acc);
                    }
                }

            } else if (bH == aH)
            {
                // b == a
                scUA = scUB;
                if (empty)
                {
                    disUATmp = scUB;
                    disUAMin = disUATmp;
                } else if (scUB.minY <= disUAMin.maxY or scUB.vX.back() > disUAMin.vX.back())
                {
                    disUATmp = scUB;
                }
            } else
            {
                LPFunction disBA = Tree[rank[b]].disOut[aH];
                if (empty)
                {
                    disUATmp = scUB.LPFCatSupport(disBA, lBound, uBound, catSupRec, acc);
                    disUAMin = disUATmp;
                } else if (disBA.minY + scUB.minY <= disUAMin.maxY)
                {
                    disUATmp = scUB.LPFCatSupport(disBA, lBound, uBound, catSupRec, acc);
                }
            }

            if (disUATmp.vX.size() > 1 and !disUAMin.dominate(disUATmp, acc))
            {
                disUAMin = disUAMin.LPFMinSupByPlaneSweep(disUATmp, acc);
            }
            assert(disUAMin.ID1 == u and disUAMin.ID2 == a);
        }

        vector<int> changedPos;
        bool updated = disUAMin.vX.size() > 1 and !disUAOld.equal(disUAMin, changedPos);

        if (!updated)
            continue;

        Tree[uRank].disOut[aH] = disUAMin;
        Tree[uRank].cntOut[aH] = 1;

        bool newLabDec = !disUAOld.dominate(disUAMin, acc);
        // secondly, check which dis can be infected
        for (int vRank: VidtoTNid[u])
        {
            // sc(v,u) + lb(u,a)
            int &cnt = Tree[vRank].cntOut[aH];
            if (cnt == 0)
                continue;
            LPFunction &disVA = Tree[vRank].disOut[aH];
            bool changed = newLabDec or lpfIsSupportBy(disVA, disUAOld, changedPos, false);
            if (changed)
                cnt = 0;
        }

        for (int vRank: VidtoTNid[a])
        {
            // lb(u,a) + sc(a,v)
            if (Tree[vRank].height <= Tree[uRank].height or Tree[vRank].ancIDs[uH] != u)
                continue;

            int &cnt = Tree[vRank].cntIn[uH];
            if (cnt == 0)
                continue;
            LPFunction &disUV = Tree[vRank].disIn[uH];
            bool changed = newLabDec or lpfIsSupportBy(disUV, disUAOld, changedPos, false);
            if (changed)
                cnt = 0;
        }
    }
}

void Graph::H2HExtendBatDiGraph(vector<pair<int, int>> &wBatch)
{
//    cout << "Index maintenance for " << wBatch.size() << " extension updates" << endl;
    unordered_map < pair<int, int>, LPFunction, boost::hash<pair<int, int>> > OCdisOld;

    // NodeOrderss.clear();
    NodeOrderss.assign(NodeOrder.begin(), NodeOrder.end());
    vector<set<int>> SCreIn, SCreOut; //SCre.clear();
    set<int> ss; //ss.clear();
    SCreIn.assign(nodeNum, ss), SCreOut.assign(nodeNum, ss); //{vertexID, set<int>}
    set < OrderCompp > OC;
    OC.clear();//vertexID in decreasing Node order

    int debugID1 = 1, debugID2 = 1;
    for (const auto &k: wBatch)
    {
        int a = k.first;
        int b = k.second;
        if (a == debugID1 and b == debugID2)
            cout << 1;

        if (NodeOrder[a] < NodeOrder[b])
        {
            for (auto &j: Tree[rank[a]].vertOut)
            {
                if (j.first == b and j.second.first.vSupportContains(a))
                {
//                    j.second.first.extendFunction(lBound, uBound);
                    OCdisOld[{a, b}] = j.second.first;
                    SCreOut[a].insert(b);
                    OC.insert(OrderCompp(a));
                    break;
                }
            }
        } else
        {
            for (auto &j: Tree[rank[b]].vertIn)
            {
                if (j.first == a and j.second.first.vSupportContains(b))
                {
//                    j.second.first.extendFunction(lBound, uBound);
                    OCdisOld[{a, b}] = j.second.first;
                    SCreIn[b].insert(a);
                    OC.insert(OrderCompp(b));
                    break;
                }
            }
        }
    }
//    cout << "** finish sc refresh" << endl;

    vector<int> ProBeginVertexSet;
    vector<int> ProBeginVertexSetNew;
    int ProID;
    bool influence;

    Timer clock;
    clock.tick();
    while (!OC.empty())
    {
        ProID = (*OC.begin()).x;
        if (ProID == -1)
            cout << 1;

        OC.erase(OC.begin());
        influence = false;

//        cout << "ProID: " << ProID << " " << Tree[rank[ProID]].disIn.size() << endl;
        for (const auto &Cid: SCreIn[ProID])
        {
            if (Cid == debugID1 and ProID == debugID2)
                cout << 1;
            int CidH = Tree[rank[Cid]].height - 1;
            /* check if sc(cid,pid) is really influenced, if so, update it with the new function */
            LPFunction scCPNew(Cid, ProID, lBound, uBound);
            for (const auto &i: adjEdge[Cid])
            {
                if (i.first == ProID)
                {
                    scCPNew = vEdge[i.second].lpf;
                    break;
                }
            }
            LPFunction scCPOld = OCdisOld[{Cid, ProID}];

            for (auto &vs: SCconNodesMT[Cid][ProID])
            {
                LPFunction CidToVS, vsToProId;
                for (auto &i: Tree[rank[vs]].vertIn)
                {
                    if (i.first == Cid)
                    {
//                        if (i.second.first.upperBound < uBound)
//                            i.second.first.extendFunction(lBound, uBound);
                        CidToVS = i.second.first;
                        break;
                    }
                }
                for (auto &i: Tree[rank[vs]].vertOut)
                {
                    if (i.first == ProID)
                    {
//                        if (i.second.first.upperBound < uBound)
//                            i.second.first.extendFunction(lBound, uBound);
                        vsToProId = i.second.first;
                        break;
                    }
                }

                LPFunction CPTmp;
                CatSupRec &catSupRec = intermediateSCs[Cid][vs][ProID];
                if (CidToVS.minY + vsToProId.minY <= scCPNew.maxY or scCPNew.vX.size() < 2)
                {
                    CPTmp = CidToVS.LPFCatSupportExtend(vsToProId, itvLen, catSupRec, acc);
                }

                if (CPTmp.vX.size() > 1 and !scCPNew.dominate(CPTmp, acc))
                {
                    scCPNew = scCPNew.LPFMinSupByPlaneSweep(CPTmp, acc);
                }
            }
            assert(scCPNew.vX.size() > 1);

            vector<int> changedPos;
            if (scCPOld.equal(scCPNew, changedPos))
                continue;

            for (auto &j: Tree[rank[ProID]].vertIn)
            {
                if (j.first == Cid)
                {
                    assert(scCPNew.ID1 == Cid and scCPNew.ID2 == ProID);
                    j.second = {scCPNew, 1};
                    break;
                }
            }

            // HneiOut: <accessories contract after Cid, lpf from ProId to hid
            unordered_set<int> HneiOut;
            // LneiOut: accessories contract after ProId and before Cid, lpf from ProId to lid
            vector<int> LneiOut;
            for (auto &j: Tree[rank[ProID]].vertOut)
            {
                // j.first is the accessory of ProID
                if (NodeOrder[j.first] > NodeOrder[Cid])
                    HneiOut.insert(j.first);
                else if (NodeOrder[j.first] < NodeOrder[Cid])
                    LneiOut.emplace_back(j.first);
            }

            // in theory, scCPOld should dominate scCPNew, but due to rounding, domination may not hold
            bool newSCDec = !scCPOld.dominate(scCPNew, acc); // assert(newSCDec == false);

            for (auto &i: Tree[rank[Cid]].vertOut)
            {
                int hid = i.first;
                if (HneiOut.find(hid) != HneiOut.end())
                {
                    bool changed =
                            newSCDec or i.second.second <= 0 or
                            lpfIsSupportBy(i.second.first, scCPOld, changedPos, true);
                    if (i.second.second > 0 and changed)
                    {
//                        if (i.second.first.upperBound < uBound)
//                            i.second.first.extendFunction(lBound, uBound);
                        OCdisOld[{Cid, hid}] = i.second.first;
                        SCreOut[Cid].insert(hid);
                        OC.insert(OrderCompp(Cid));
                        i.second.second = 0;
                    }
                }
            }

            for (auto &i: LneiOut)
            {
                int lid = i;
                for (auto &j: Tree[rank[lid]].vertIn)
                {
                    if (j.first == Cid)
                    {
                        bool changed =
                                newSCDec or j.second.second <= 0 or
                                lpfIsSupportBy(j.second.first, scCPOld, changedPos, true);
                        if (j.second.second > 0 and changed)
                        {
//                            if (j.second.first.upperBound < uBound)
//                                j.second.first.extendFunction(lBound, uBound);
                            OCdisOld[{Cid, lid}] = j.second.first;
                            SCreIn[lid].insert(Cid);
                            OC.insert(OrderCompp(lid));
                            j.second.second = 0;
                        }
                        break;
                    }
                }
            }

            int pAncs = (int) Tree[rank[ProID]].ancIDs.size();
            for (int i = 0; i < pAncs; i++)
            {
                // i -> cid -> pid
                LPFunction &disIP = Tree[rank[ProID]].disIn[i];
                if (newSCDec)
                {
                    if (i == CidH)
                        influence = true;
                    Tree[rank[ProID]].cntIn[i] = 0;
                } else if (lpfIsSupportBy(disIP, scCPOld, changedPos, false))
                {
                    Tree[rank[ProID]].cntIn[i] = 0;
                    if (i == CidH)
                        influence = true;
                }
            }

//            scCPNew.scSupToLbSup(ProID);
            CatSupRec catSupRec;
            catSupRec.vX1 = scCPNew.vX;
            catSupRec.vX2 = scCPNew.vX;
            catSupRec.vY = scCPNew.vY;
//            catSupRec.extendEnd(uBound + itvLen, itvLen);
            intermediateLBsIn[ProID][ProID][Cid] = catSupRec;
        }

        for (const auto &Cid: SCreOut[ProID])
        {
            if (Cid == debugID2 and ProID == debugID1)
                cout << 1;
            int CidH = Tree[rank[Cid]].height - 1;

            /* check if sc(pid,cid) is really influenced, if so, update it with the new function */
            LPFunction scPCNew(ProID, Cid, lBound, uBound);
            for (const auto &i: adjEdge[ProID])
            {
                if (i.first == Cid)
                {
                    scPCNew = vEdge[i.second].lpf;
                    break;
                }
            }

            LPFunction scPCOld = OCdisOld[{ProID, Cid}];

            for (auto &vs: SCconNodesMT[ProID][Cid])
            {
                LPFunction ProidToVS, vsToCid;
                for (auto &i: Tree[rank[vs]].vertIn)
                {
                    if (i.first == ProID)
                    {
//                        if (i.second.first.upperBound < uBound)
//                            i.second.first.extendFunction(lBound, uBound);
                        ProidToVS = i.second.first;
                        break;
                    }
                }

                for (auto &i: Tree[rank[vs]].vertOut)
                {
                    if (i.first == Cid)
                    {
//                        if (i.second.first.upperBound < uBound)
//                            i.second.first.extendFunction(lBound, uBound);
                        vsToCid = i.second.first;
                        break;
                    }
                }

                assert(ProidToVS.ID1 == ProID and vsToCid.ID2 == Cid);
                LPFunction PCTmp;
                if (ProidToVS.minY + vsToCid.minY <= scPCNew.maxY or scPCNew.vX.size() < 2)
                {
                    CatSupRec &catSupRec = intermediateSCs[ProID][vs][Cid];
                    PCTmp = ProidToVS.LPFCatSupportExtend(vsToCid, itvLen, catSupRec, acc);
//                    catSupRec.extendEnd(uBound + itvLen, itvLen);
                }

                if (PCTmp.vX.size() > 1 and !scPCNew.dominate(PCTmp, acc))
                {
                    scPCNew = scPCNew.LPFMinSupByPlaneSweep(PCTmp, acc);
                }
            }
            assert(scPCNew.vX.size() > 1);

            vector<int> changedPos;
            if (scPCOld.equal(scPCNew, changedPos))
                continue;
            for (auto &j: Tree[rank[ProID]].vertOut)
            {
                if (j.first == Cid)
                {
                    assert(scPCNew.ID1 == ProID and scPCNew.ID2 == Cid);
                    j.second = {scPCNew, 1};
                    break;
                }
            }

            unordered_set<int> HneiIn;
            vector<int> LneiIn;
            for (auto &i: Tree[rank[ProID]].vertIn)
            {
                // j.first is the accessory of ProID
                if (NodeOrder[i.first] > NodeOrder[Cid])
                    HneiIn.insert(i.first);
                else if (NodeOrder[i.first] < NodeOrder[Cid])
                    LneiIn.emplace_back(i.first);
            }

            // in theory, scCPOld should dominate scCPNew, but due to rounding, domination may not hold
            bool newSCDec = !scPCOld.dominate(scPCNew, acc); // assert(newSCDec == false);

            for (auto &i: Tree[rank[Cid]].vertIn)
            {
                int hid = i.first;
                if (HneiIn.find(hid) != HneiIn.end())
                {
                    bool changed =
                            newSCDec or i.second.second <= 0 or
                            lpfIsSupportBy(i.second.first, scPCOld, changedPos, true);
                    if (i.second.second > 0 and changed)
                    {
//                        if (i.second.first.upperBound < uBound)
//                            i.second.first.extendFunction(lBound, uBound);
                        OCdisOld[{hid, Cid}] = i.second.first;
                        SCreIn[Cid].insert(hid);
                        OC.insert(OrderCompp(Cid));
                        i.second.second = 0;
                    }
                }
            }

            for (auto &lid: LneiIn)
            {
                for (auto &i: Tree[rank[lid]].vertOut)
                {
                    if (i.first == Cid)
                    {
                        bool changed =
                                newSCDec or i.second.second <= 0 or
                                lpfIsSupportBy(i.second.first, scPCOld, changedPos, true);
                        if (i.second.second > 0 and changed)
                        {
//                            if (i.second.first.upperBound < uBound)
//                                i.second.first.extendFunction(lBound, uBound);
                            OCdisOld[{lid, Cid}] = i.second.first;
                            SCreOut[lid].insert(Cid);
                            OC.insert(OrderCompp(lid));
                            i.second.second = 0;
                        }
                        break;
                    }
                }
            }

            influence = true;
            int pAncs = (int) Tree[rank[ProID]].ancIDs.size();
            for (int i = 0; i < pAncs; i++)
            {
                LPFunction &disPI = Tree[rank[ProID]].disOut[i], disCI;
                if (newSCDec)
                {
                    if (i == CidH)
                        influence = true;
                    Tree[rank[ProID]].cntOut[i] = 0;
                } else if (lpfIsSupportBy(disPI, scPCOld, changedPos, false))
                {
                    Tree[rank[ProID]].cntOut[i] = 0;
                    if (i == CidH)
                        influence = true;
                }
            }

//            scPCNew.scSupToLbSup(ProID);
            CatSupRec catSupRec;
            catSupRec.vX1 = scPCNew.vX;
            catSupRec.vX2 = scPCNew.vX;
            catSupRec.vY = scPCNew.vY;
//            catSupRec.extendEnd(uBound + itvLen, itvLen);
            intermediateLBsOut[ProID][ProID][Cid] = catSupRec;
        }

        if (influence)
        {
            ProBeginVertexSetNew.clear();
            ProBeginVertexSetNew.reserve(ProBeginVertexSet.size() + 1);
            ProBeginVertexSetNew.push_back(ProID);
            int rnew = rank[ProID], r;
            for (int i: ProBeginVertexSet)
            {
                r = rank[i];
                if (LCAQuery(rnew, r) != rnew)
                {
                    ProBeginVertexSetNew.push_back(i);
                }
            }
            ProBeginVertexSet = ProBeginVertexSetNew;
        }
    }
    clock.tock();
//    cout << "** finish bottom-up refresh: " << clock.duration().count() << " ms" << endl;

    clock.tick();
    EachNodeExtendMainThread(ProBeginVertexSet);
    clock.tock();
//    cout << "** time for top-down refresh: " << clock.duration().count() << " ms" << endl;
//    lBound += itvLen;
//    uBound += itvLen;

//    clock.tick();
//    extendIndex();
//    clock.tock();
//    cout << "    Extending time: " << clock.duration().count() << " ms" << endl;
}

void Graph::EachNodeExtendMainThread(const vector<int> &sources)
{
    vector<int> uRanks;
    uRanks.reserve(sources.size());
    int avgHeight = 0;
    for (int ProBeginVertexID: sources)
    {
//        cout << rank[ProBeginVertexID] << " ";
        uRanks.push_back(rank[ProBeginVertexID]);
        avgHeight += Tree[rank[ProBeginVertexID]].height;
    }
    cout << "    sources size: " << sources.size() << " avg height: " << avgHeight / sources.size() << endl;

    int thisThreadNum = threadNum;
    Timer clock;
    clock.tick();
    bool cond = false;
    while (!uRanks.empty() and (uRanks.size() > thisThreadNum * 2 or uRanks.size() < thisThreadNum * 0.1))
    {
//        vector<pair<int, int>> vpairs;
//        vpairs.reserve(uRanks.size());
//        int step = uRanks.size()/thisThreadNum;
//        for (int i = 0; i < thisThreadNum - 1; i++)
//        {
//            vpairs.push_back({i * step, (i + 1) * step});
//        }
//        vpairs.push_back({(thisThreadNum - 1) * step, uRanks.size()});
//
//        boost::thread_group threads;
//        for (const auto &p:vpairs)
//        {
//            threads.create_thread(
//                    boost::bind(&Graph::EachNodeExtendMainRange,
//                                this, p.first, p.second, boost::ref(uRanks)));
//        }
//        threads.join_all();


        boost::thread_group threads;
        for (const auto &p: uRanks)
        {
            threads.create_thread(
                    boost::bind(&Graph::EachNodeExtendMain, this, p));
        }
        threads.join_all();

        vector<int> newURanks;
        for (const auto &ur: uRanks)
        {
            for (int j: Tree[ur].vChildren)
            {
                newURanks.push_back(j);
            }
        }
        uRanks = newURanks;
    }
    clock.tock();

    if (uRanks.empty())
        return;
    boost::thread_group threads;
    for (const auto &urank: uRanks)
    {
        threads.create_thread(boost::bind(&Graph::EachNodeExtendMainRec, this, urank));
    }
    threads.join_all();
}

void Graph::EachNodeExtendMainRange(int begin, int end, const vector<int> &uRank)
{
    for (int i = begin; i < end; i++)
    {
        EachNodeExtendMain(uRank[i]);
    }
}

void Graph::EachNodeExtendMainRec(int uRank)
{
    EachNodeExtendMain(uRank);
    for (int j: Tree[uRank].vChildren)
        EachNodeExtendMainRec(j);
}

void Graph::EachNodeExtendMain(int uRank)
{
    int u = Tree[uRank].uniqueVertex;
    int uH = Tree[uRank].height - 1;

    vector<int> line = Tree[uRank].ancIDs;

    for (int aH = 0; aH < line.size(); aH++)
    {
        int a = line[aH];
//        if (Tree[uRank].disIn[aH].upperBound < uBound)
//            Tree[uRank].disIn[aH].extendFunction(lBound, uBound);
        LPFunction disAUOld = Tree[uRank].disIn[aH];

        if (Tree[uRank].cntIn[aH] > 0) continue;
        //firstly, calculate the actual distance
        int b, bH;
        LPFunction disAUMin(a, u, lBound, uBound), scBU;
        for (int j = 0; j < Tree[uRank].vertIn.size(); j++)
        {
            // a -> b -> u
            b = Tree[uRank].vertIn[j].first;
//            if (Tree[uRank].vertIn[j].second.first.upperBound < uBound)
//                Tree[uRank].vertIn[j].second.first.extendFunction(lBound, uBound);
            scBU = Tree[uRank].vertIn[j].second.first;
            scBU.scSupToLbSup(u);
            bH = Tree[rank[b]].height - 1;
            LPFunction disAUTmp;
            CatSupRec &catSupRec = intermediateLBsIn[u][b][a];
            if (bH < aH)
            {
                LPFunction &disAB = Tree[rank[a]].disOut[bH];
                if (disAB.minY + scBU.minY <= disAUMin.maxY)
                {
//                    if (disAB.upperBound < uBound)
//                        disAB.extendFunction(lBound, uBound);
                    disAUTmp = disAB.LPFCatSupportExtend(scBU, itvLen, catSupRec, acc);
                }
            } else if (bH == aH)
            {
                // b == a
                if (scBU.minY <= disAUMin.maxY)
                    disAUTmp = scBU;
            } else
            {
                LPFunction disAB = Tree[rank[b]].disIn[aH];
                if (disAB.minY + scBU.minY <= disAUMin.maxY or disAUMin.vX.size() < 2)
                {
//                    if (disAB.upperBound < uBound)
//                        disAB.extendFunction(lBound, uBound);
                    disAUTmp = disAB.LPFCatSupportExtend(scBU, itvLen, catSupRec, acc);
//                    catSupRec.extendEnd(uBound + itvLen, itvLen);
                }
            }

            if (disAUTmp.vX.size() > 1 and u == 1 and a == 1)
            {
                cout << "~~~~~~~~~~~~~~~~~~~" << endl;
                disAUMin.display();
                disAUTmp.display();
            }
            if (!disAUMin.dominate(disAUTmp, acc))
            {
                disAUMin = disAUMin.LPFMinSupByPlaneSweep(disAUTmp, acc);
            }
            assert(disAUMin.ID1 == a and disAUMin.ID2 == u);
        }

        vector<int> changedPos;
        bool updated = disAUMin.vX.size() > 1 and !disAUOld.equal(disAUMin, changedPos);
        if (!updated)
            continue;

        Tree[uRank].disIn[aH] = disAUMin;
        Tree[uRank].cntIn[aH] = 1;

        bool newLabDec = !disAUOld.dominate(disAUMin, acc);
        // firstly, check which dis can be infected
        for (int vRank: VidtoTNid[u])
        {
            // lb(a, u) + sc(u, v)
            int &cnt = Tree[vRank].cntIn[aH];
            if (cnt == 0)
                continue;
            LPFunction &disAV = Tree[vRank].disIn[aH];
            bool changed = newLabDec or lpfIsSupportBy(disAV, disAUOld, changedPos, false);
            if (changed)
                cnt = 0;
        }

        for (int vRank: VidtoTNid[a])
        {
            // sc(v, a) + lb(a, u)
            if (Tree[vRank].height <= Tree[uRank].height or Tree[vRank].ancIDs[uH] != u)
                continue;
            int &cnt = Tree[vRank].cntOut[uH];
            if (cnt == 0)
                continue;
            LPFunction &disVU = Tree[vRank].disOut[uH];
            bool changed = newLabDec or lpfIsSupportBy(disVU, disAUOld, changedPos, false);
            if (changed)
                cnt = 0;
        }
    }

    for (int aH = 0; aH < line.size(); aH++)
    {
        int a = line[aH];
//        if (Tree[uRank].disOut[aH].upperBound < uBound)
//            Tree[uRank].disOut[aH].extendFunction(lBound, uBound);
        LPFunction disUAOld = Tree[uRank].disOut[aH];

        if (Tree[uRank].cntOut[aH] > 0) continue;
        //firstly, calculate the actual distance
        int b, bH;

        LPFunction disUAMin(u, a, lBound, uBound), scUB;
        for (int j = 0; j < Tree[uRank].vertOut.size(); j++)
        {
            // a <- b <- u
            b = Tree[uRank].vertOut[j].first;
//            if (Tree[uRank].vertOut[j].second.first.upperBound < uBound)
//                Tree[uRank].vertOut[j].second.first.extendFunction(lBound, uBound);
            scUB = Tree[uRank].vertOut[j].second.first;
            scUB.scSupToLbSup(u);
            bH = Tree[rank[b]].height - 1;
            bool empty = disUAMin.vX.size() < 2;
            LPFunction disUATmp;
            CatSupRec &catSupRec = intermediateLBsOut[u][b][a];
            if (bH < aH)
            {
                LPFunction disBA = Tree[rank[a]].disIn[bH];
                if (disBA.minY + scUB.minY <= disUAMin.maxY)
                {
//                    if (disBA.upperBound < uBound)
//                        disBA.extendFunction(lBound, uBound);
                    disUATmp = scUB.LPFCatSupportExtend(disBA, itvLen, catSupRec, acc);
//                    catSupRec.extendEnd(uBound + itvLen, itvLen);
                }
            } else if (bH == aH)
            {
                // b == a
                if (scUB.minY <= disUAMin.maxY)
                    disUATmp = scUB;
            } else
            {
                LPFunction disBA = Tree[rank[b]].disOut[aH];
                if (disBA.minY + scUB.minY <= disUAMin.maxY or disUAMin.vX.size() < 2)
                {
//                    if (disBA.upperBound < uBound)
//                        disBA.extendFunction(lBound, uBound);
                    disUATmp = scUB.LPFCatSupportExtend(disBA, itvLen, catSupRec, acc);
//                    catSupRec.extendEnd(uBound + itvLen, itvLen);
                }
            }

            if (!disUAMin.dominate(disUATmp, acc))
            {
                disUAMin = disUAMin.LPFMinSupByPlaneSweep(disUATmp, acc);
            }
            assert(disUAMin.ID1 == u and disUAMin.ID2 == a);
        }

        vector<int> changedPos;
        bool updated = disUAMin.vX.size() > 1 and !disUAOld.equal(disUAMin, changedPos);

        if (!updated)
            continue;

        Tree[uRank].disOut[aH] = disUAMin;
        Tree[uRank].cntOut[aH] = 1;

        bool newLabDec = !disUAOld.dominate(disUAMin, acc);
        // secondly, check which dis can be infected
        for (int vRank: VidtoTNid[u])
        {
            // sc(v,u) + lb(u,a)
            int &cnt = Tree[vRank].cntOut[aH];
            if (cnt == 0)
                continue;
            LPFunction &disVA = Tree[vRank].disOut[aH];
            bool changed = newLabDec or
                           lpfIsSupportBy(disVA, disUAOld, changedPos, false);
            if (changed)
                cnt = 0;
        }

        for (int vRank: VidtoTNid[a])
        {
            // lb(u,a) + sc(a,v)
            if (Tree[vRank].height <= Tree[uRank].height or Tree[vRank].ancIDs[uH] != u)
                continue;

            int &cnt = Tree[vRank].cntIn[uH];
            if (cnt == 0)
                continue;
            LPFunction &disUV = Tree[vRank].disIn[uH];
            bool changed = newLabDec or
                           lpfIsSupportBy(disUV, disUAOld, changedPos, false);
            if (changed)
                cnt = 0;
        }
    }
}

void Graph::timeWinSwitchSimpleExtend(int win1, int win2)
{
    vector<pair<int, int>> vpID;
    int thisThreadNum = threadNum;
    int step = nodeNum / thisThreadNum;

    vpID.reserve(thisThreadNum - 1);
    for (int i = 0; i < thisThreadNum - 1; i++)
    {
        vpID.emplace_back(i * step, (i + 1) * step);
    }
    vpID.emplace_back((thisThreadNum - 1) * step, nodeNum);

    boost::thread_group threads;
    for (int i = 0; i < thisThreadNum; i++)
    {
        threads.add_thread(
                new boost::thread(
                        &Graph::timeWinSwitchSimpleExtendRange, this, vpID[i].first, vpID[i].second
                ));
    }
    threads.join_all();
}

void Graph::timeWinSwitchSimpleExtendRange(int begin, int end)
{
    for (int i = begin; i < end; i++)
    {
//        for (auto &j: adjEdge[i])
//        {
//            vEdge[j.second].lpf.extendFunction(lBound, uBound);
//        }
//
//        for (auto &p: Tree[i].vertIn)
//        {
//            p.second.first.extendFunction(lBound, uBound);
//        }
//
//        for (auto &p: Tree[i].vertOut)
//        {
//            p.second.first.extendFunction(lBound, uBound);
//        }

//        int pid = Tree[i].uniqueVertex;

//        for (const auto &ps: SCconNodesMT[pid])
//        {
//            for (const auto &j: ps.second)
//            {
//                intermediateSCs[pid][j][ps.first].extendEnd(uBound, itvLen);
//            }
//        }

//        for (auto &j: adjEdge[pid])
//        {
//            intermediateSCs[pid][pid][j.first].extendEnd(uBound, itvLen);
//        }
//
//        for (const auto &p: Tree[i].vertIn)
//        {
//            int vid = p.first;
//            for (const auto &aid: Tree[i].ancIDs)
//            {
//                intermediateLBsIn[pid][vid][aid].extendEnd(uBound, itvLen);
//                intermediateLBsOut[pid][vid][aid].extendEnd(uBound, itvLen);
//            }
//            intermediateLBsIn[pid][pid][vid].extendEnd(uBound, itvLen);
//            intermediateLBsOut[pid][pid][vid].extendEnd(uBound, itvLen);
//        }

//        for (auto &p: Tree[i].disIn)
//        {
//            p.extendFunction(lBound, uBound);
//        }
//
//        for (auto &p: Tree[i].disOut)
//        {
//            p.extendFunction(lBound, uBound);
//        }
    }
}
