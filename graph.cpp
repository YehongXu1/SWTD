#include "graph.h"

vector<int> _DD, _DD2;

struct DegComp
{
    int x;

    DegComp(int _x)
    {
        x = _x;
    }

    bool operator<(const DegComp d) const
    {
        if (_DD[x] != _DD[d.x])
            return _DD[x] < _DD[d.x];
        if (_DD2[x] != _DD2[d.x])
            return _DD2[x] < _DD2[d.x];
        return x < d.x;
    }
};

void Graph::detNodeOrder()
{
    //Contracted Graph E
    //Neighbor, Function
    vector<set<int>> ES(nodeNum);
    vector<set<int>> ESR(nodeNum);

    for (int i = 0; i < (int) adjEdge.size(); i++)
    {
        for (int j = 0; j < (int) adjEdge[i].size(); j++)
        {
            ES[i].insert(adjEdge[i][j].first);
        }
    }

    for (int i = 0; i < (int) adjEdgeR.size(); i++)
    {
        for (int j = 0; j < (int) adjEdgeR[i].size(); j++)
        {
            // there is an edge (j, i) in the original graph
            ESR[i].insert(adjEdgeR[i][j].first);
        }
    }

    _DD.assign(nodeNum, 0);
    _DD2.assign(nodeNum, 0);
    DD.assign(nodeNum, 0);
    DD2.assign(nodeNum, 0);

    set < DegComp > Deg;
    int degree;
    for (int i = 0; i < nodeNum; i++)
    {
        degree = (int) (adjEdge[i].size() + adjEdgeR[i].size());
        if (degree != 0)
        {
            _DD[i] = degree;
            _DD2[i] = degree;
            DD[i] = degree;
            DD2[i] = degree;
            Deg.insert(DegComp(i));
        }
    }

    vector<bool> exist(nodeNum, true);
    //If Degree Changed by Contraction
    vector<bool> change(nodeNum, false);

    int count = 0;
    while (!Deg.empty())
    {
        count++;
        int x = (*Deg.begin()).x;

        while (true)
        {
            if (change[x])
            {
                Deg.erase(DegComp(x));
                _DD[x] = DD[x];
                _DD2[x] = DD2[x];
                Deg.insert(DegComp(x));
                change[x] = false;
                x = (*Deg.begin()).x;
            } else
                break;
        }

        vNodeOrder.emplace_back(x);
        Deg.erase(Deg.begin());
        exist[x] = false;

        //Tmp neighbors of the contracted graph
        vector<int> Neigh, NeighR;
        for (auto &it: ES[x])
            if (exist[it])
                Neigh.emplace_back(it);

        for (auto &it: ESR[x])
            if (exist[it])
                NeighR.emplace_back(it);

        //Maintain E
        for (auto &y: Neigh)
        {
//            deleteE(x, y);
            int u = x, v = y;
            if (ES[u].find(v) != ES[u].end())
            {
                ES[u].erase(ES[u].find(v));
                DD[u]--;
            }

            if (ESR[v].find(u) != ESR[v].end())
            {
                ESR[v].erase(ESR[v].find(u));
                DD[v]--;
            }
            change[y] = true;
        }

        for (auto &y: NeighR)
        {
            int u = y, v = x;
            if (ES[u].find(v) != ES[u].end())
            {
                ES[u].erase(ES[u].find(v));
                DD[u]--;
            }

            if (ESR[v].find(u) != ESR[v].end())
            {
                ESR[v].erase(ESR[v].find(u));
                DD[v]--;
            }
        }

        for (auto &ID1: NeighR)
        {
            for (auto &ID2: Neigh)
            {
                if (ID1 == ID2)
                    continue;

                if (ES[ID1].find(ID2) == ES[ID1].end())
                {
                    ES[ID1].insert(ID2);
                    DD[ID1]++;
                    DD2[ID1]++;
                }

                if (ESR[ID2].find(ID1) == ESR[ID2].end())
                {
                    ESR[ID2].insert(ID1);
                    DD[ID2]++;
                    DD2[ID2]++;
                }
                //Degree Changed
                change[ID1] = true;
                change[ID2] = true;
            }
        }
    }

//    assert(vNodeOrder.size() == nodeNum);
    NodeOrder.assign(nodeNum, -1);
    for (int k = 0; k < (int) vNodeOrder.size(); k++)
    {
        NodeOrder[vNodeOrder[k]] = k;
//        cout << vNodeOrder[k] << ": " << k << endl;
    }
}


int Graph::readBeijingMapDirectedNew(string filenameGraph, string filenameMap)
{
    ifstream inGraph(filenameGraph.c_str());
    if (!inGraph)
    {
        cout << "Cannot open Beijing map " << filenameGraph << endl;
        return -1;
    }

    cout << "Reading Graph Data " << filenameGraph << endl;
    inGraph >> this->nodeNum;
    cout << "Nodenum:" << this->nodeNum << endl;
    vector<pair<int, int> > vp;
    adjList.assign(this->nodeNum, vp);
    adjListR.assign(this->nodeNum, vp);
    adjEdge.assign(this->nodeNum, vp);
    adjEdgeR.assign(this->nodeNum, vp);

    int nodeID, neighborNum, neighborNodeID, neighborLength, neighborRoadID;
    double lon, lat;
    for (int i = 0; i < this->nodeNum; i++)
    {
        inGraph >> nodeID >> neighborNum;
        for (int j = 0; j < neighborNum; j++)
        {
            inGraph >> neighborNodeID >> neighborLength >> neighborRoadID;
            adjList[nodeID].push_back(make_pair(neighborNodeID, neighborLength));
            adjListR[neighborNodeID].push_back(make_pair(nodeID, neighborLength));

            adjEdge[nodeID].push_back(make_pair(neighborNodeID, neighborRoadID));
            adjEdgeR[neighborNodeID].push_back(make_pair(nodeID, neighborRoadID));
        }
    }
    inGraph.close();

    ifstream inMap(filenameMap);
    if (!inMap)
    {
        cout << "Cannot open clean map data " << filenameMap << endl;
        return -1;
    }
    cout << "Reading Node Map Data " << filenameMap << endl;
    int ID1, ID2;
    vNodeMapR.reserve(this->nodeNum);

    while (inMap >> ID1)
    {
        inMap >> ID2;
        vNodeMap.push_back(ID2);
        if (ID2 != -1)
            vNodeMapR[ID2] = ID1;
    }
    inMap.close();

    return 0;
}

int Graph::readBeijingTD(string filenameRoad, string filenameSpeed)
{
    int roadNum;
    ifstream inRoad(filenameRoad);
    if (!inRoad)
    {
        cout << "Cannot open Road file " << filenameRoad << endl;
        return -1;
    }
    inRoad >> roadNum;
    cout << "Reading Road Data " << filenameRoad << endl;
    Edge eTmp;
    vEdge.assign(roadNum, eTmp);
    for (int i = 0; i < roadNum; i++)
    {
        Edge e;
        inRoad >> e.edgeID >> e.ID1 >> e.ID2 >> e.length >> e.twinID;
        vEdge[e.edgeID] = e;
    }
    inRoad.close();

    ifstream inSpeed(filenameSpeed);
    if (!inSpeed)
    {
        cout << "Cannot open speed file " << filenameSpeed << endl;
        return -1;
    }
    cout << "Reading Speed Data " << filenameSpeed << endl;
    string s;
    stringstream ss;
    int slotNum = itvLen / 300;
    for (int i = 0; i < roadNum; i++)
    {
        inSpeed >> s;
        vector<string> vs = Tools::split(s, " ");
        vector<int> vX, vY;
        vX.reserve(288 / slotNum);
        vY.assign(288 / slotNum, -1);
        int count = 0;
        for (int j = 0; j < 288; j++)
        {
            if (j % slotNum == 0)
            {
                vX.push_back(j * 600);
                ss.clear();
                ss.str("");
                ss << vs[j];
                ss >> vY[int(j / slotNum)];
            }
        }
        vEdge[i].vXFull = vX;
        vEdge[i].vYFull = vY;

        vEdge[i].vX.assign(vX.begin(), vX.begin() + 5);
        vEdge[i].vY.assign(vY.begin(), vY.begin() + 5);
//        vector<int> vSupportFull(vX.size() - 1, vEdgeInBG[i].ID1);
        vector<vector<int>> vSupport;
        if (NodeOrder[vEdge[i].ID1] < NodeOrder[vEdge[i].ID2])
            vSupport.assign(4, {NodeOrder[vEdge[i].ID1]});
        else
            vSupport.assign(4, {NodeOrder[vEdge[i].ID2]});
        vEdge[i].lpf = LPFunction(vEdge[i].ID1, vEdge[i].ID2, lBound, uBound, vEdge[i].vX, vEdge[i].vY, acc);

    }
    inSpeed.close();

    return 0;
}

void Graph::readExtension(vector<pair<int,int>> &testdata, int newUpperbound)
{
    testdata.clear();
    for (auto &edge: vEdge)
    {
        assert(newUpperbound == edge.vX.back());
        int ID1 = edge.ID1, ID2 = edge.ID2;
        int actualWeight = edge.vYFull[newUpperbound/itvLen];
//        if (edge.lpf.upperBound < appendPoint * 300)
//            edge.lpf.extendFunction(appendPoint * 300 - 7200, appendPoint * 300);
        if (actualWeight != edge.lpf.vY.back())
        {
            testdata.push_back({ID1,ID2});
            edge.lpf.vY[edge.lpf.vY.size() - 1] = actualWeight;
            edge.vY[edge.vY.size() - 1] = actualWeight;
//            intermediateSCs[ID1][ID1][ID2].vY.end()[-1] = actualWeight;
        }

        edge.lpf = LPFunction(
                ID1, ID2, lBound, newUpperbound, edge.vX, edge.vY, acc);
        int lVid = NodeOrder[ID1] < NodeOrder[ID2] ? ID1 : ID2;
        for (int j = 0; j < edge.lpf.vX.size() - 1; j++)
        {
            edge.lpf.vSupportPre.push_back({{lVid, {j}}});
        }
    }
}

void Graph::readUndirectedGraph(string mapfile, string speedfile, int lowerB, int upperB, int slotNum)
{
    ifstream in(mapfile.c_str());
    if (!in)
    {
        cout << "Cannot open file " << mapfile << endl;
        return;
    }
    LPFunction lpf;
    lBound = lowerB, uBound = upperB, deltaT = 300 * slotNum;
    int winNum1 = lBound / deltaT, winNum2 = upperB / deltaT;


    string x;
    in >> x >> nodeNum >> edgeNum;


    adjList.assign(nodeNum + 1, vector<pair<int, int> >());
    adjListR.assign(nodeNum + 1, vector<pair<int, int> >());
    adjEdge.assign(nodeNum + 1, vector<pair<int, int> >());
    adjEdgeR.assign(nodeNum + 1, vector<pair<int, int> >());

    Edge eTmp;
    vEdge.assign(edgeNum, eTmp);

    int ID1, ID2, roadID = 0;
    while (in >> ID1 >> ID2)
    {
        adjList[ID1].push_back(make_pair(ID2, 10));
        adjListR[ID2].push_back(make_pair(ID1, 10));

        adjEdge[ID1].push_back(make_pair(ID2, roadID));
        adjEdgeR[ID2].push_back(make_pair(ID1, roadID));

        Edge e;
        e.edgeID = roadID;
        e.ID1 = ID1;
        e.ID2 = ID2;
        vEdge[e.edgeID] = e;
        roadID += 1;
    }
    in.close();

    detNodeOrder();
    ifstream inSpeed(speedfile.c_str());
    if (!inSpeed)
    {
        cout << "Cannot open speed file " << speedfile << endl;
    }

    string s;
    stringstream ss;
    map<int, int> turningPointMap;
    for (int i = 0; i < roadID; i++)
    {
        int u = vEdge[i].ID1, v = vEdge[i].ID2, revRID = -1;
        for (auto &j: adjEdge[v])
        {
            if (j.first == u)
            {
                // there is same edge with reverse direction
                revRID = j.second;
                break;
            }
        }

        int lVid = NodeOrder[u] < NodeOrder[v] ? u : v;
        revRID = edgeNum + 1; // if revRID is very large, then w(a,b) and w(b,a) would be different
//        int slotnum = itvLen / 300;
        if (revRID > i)
        {
            // i'm edge(u,v), edge(v, u) is not initialized yet
            // if vEdgeInBG[i] already has
            inSpeed >> s;
            vector<string> vs = Tools::split(s, ",");
            vector<int> vX, vY;
            vX.reserve(288); // 1 day = 288 5-minutes
            vY.assign(288, -1);
            for (int j = 0; j < 288; j++)
            {
                vX.push_back(j * 300);
                ss.clear();
                ss.str("");
                ss << vs[j];
                ss >> vY[j];
            }
            // vX = {600 * i, i = 0, 1, 2,..., 285}
            vEdge[i].vXFull = vX;
            vEdge[i].vYFull = vY;

            for (int j = winNum1; j <= winNum2; j++)
            {
                assert(j * slotNum < vEdge[i].vYFull.size());
                vEdge[i].vX.push_back(vX[j * slotNum]);
                vEdge[i].vY.push_back(vY[j * slotNum]);
            }
            vEdge[i].vY[vEdge[i].vY.size() - 1] = vEdge[i].vY[vEdge[i].vY.size() - 2];

            vEdge[i].lpf = LPFunction(
                    u, v, lBound, uBound, vEdge[i].vX, vEdge[i].vY, acc);
            if (turningPointMap.find(vEdge[i].lpf.vX.size()) == turningPointMap.end())
                turningPointMap[vEdge[i].lpf.vX.size()] = 1;
            else
                turningPointMap[vEdge[i].lpf.vX.size()] += 1;
            for (int j = 0; j < vEdge[i].lpf.vX.size() - 1; j++)
            {
                vEdge[i].lpf.vSupportPre.push_back({{lVid, {j}}});
            }
        } else
        {
            // i'm edge(u,v), edge(v, u) is initialized
            vEdge[i].vXFull = vEdge[revRID].vXFull;
            vEdge[i].vYFull = vEdge[revRID].vYFull;

            vEdge[i].vX.assign(vEdge[revRID].vXFull.begin() + winNum1, vEdge[revRID].vXFull.begin() + winNum2);
            vEdge[i].vY.assign(vEdge[revRID].vYFull.begin() + winNum1, vEdge[revRID].vYFull.begin() + winNum2);
            vEdge[i].lpf = LPFunction(u, v, lBound, uBound, vEdge[revRID].vX, vEdge[revRID].vY, acc);

            for (int j = 0; j < vEdge[i].lpf.vX.size() - 1; j++)
            {
                map<int, vector<int>> mp = {{lVid, {j}}};
                vEdge[i].lpf.vSupportPre.push_back(mp);
            }
        }
    }
    inSpeed.close();

    vector<int> turningPoint, count;
    double turnSum = 0;
    for (auto &it: turningPointMap)
    {
        turningPoint.push_back(it.first);
        count.push_back(it.second);
        turnSum += it.first * it.second;
    }

//    for (const auto &v: turningPoint)
//        cout << v << ", ";
//    cout << endl;
//    for (const auto &v: count)
//        cout << v << ", ";
//    cout << endl;
    cout << "Finish loading graph, " << "|V|: " << nodeNum << ", |E|: " << edgeNum << ", Thread num: "
         << threadNum << ", Accuracy: " << acc << ", Time domain: " << lBound << " " << uBound
         << " SlotNum: " << slotNum << ", Avg turnP num: " << 1.0 * turnSum / edgeNum << endl;
}

void Graph::readUndirectedGraph(string mapfile, string speedfile, string orderfile, int slotnum)
{
    ifstream in(mapfile.c_str());
    if (!in)
    {
        cout << "Cannot open file " << mapfile << endl;
        assert(false);
    }

    string x;
    in >> x >> nodeNum >> edgeNum;

    adjList.assign(nodeNum + 1, vector<pair<int, int> >());
    adjListR.assign(nodeNum + 1, vector<pair<int, int> >());
    adjEdge.assign(nodeNum + 1, vector<pair<int, int> >());
    adjEdgeR.assign(nodeNum + 1, vector<pair<int, int> >());

    Edge eTmp;
    vEdge.assign(edgeNum, eTmp);

    int ID1, ID2, roadID = 0;
    while (in >> ID1 >> ID2)
    {
        adjList[ID1].push_back(make_pair(ID2, 10));
        adjListR[ID2].push_back(make_pair(ID1, 10));

        adjEdge[ID1].push_back(make_pair(ID2, roadID));
        adjEdgeR[ID2].push_back(make_pair(ID1, roadID));

        Edge e;
        e.edgeID = roadID;
        e.ID1 = ID1;
        e.ID2 = ID2;
        vEdge[e.edgeID] = e;
        roadID += 1;
    }
    in.close();

    ifstream in2(orderfile.c_str());
    if (!in2)
    {
        cout << "Cannot open order file " << orderfile << endl;
        assert(false);
    }

    vNodeOrder.reserve(nodeNum);
    NodeOrder.assign(nodeNum, -1);
    int vid, k = 0;
    while (in2 >> vid)
    {
        vNodeOrder.push_back(vid);
        NodeOrder[vNodeOrder[k]] = k;
        k += 1;
    }

    ifstream inSpeed(speedfile.c_str());
    if (!inSpeed)
    {
        cout << "Cannot open speed file " << speedfile << endl;
        assert(false);
    }
    string s;
    stringstream ss;
    for (int i = 0; i < roadID; i++)
    {
        int u = vEdge[i].ID1, v = vEdge[i].ID2, revRID = -1;
        for (auto &j: adjEdge[v])
        {
            if (j.first == u)
            {
                // there is same edge with reverse direction
                revRID = j.second;
                break;
            }
        }

        int lVid = NodeOrder[u] < NodeOrder[v] ? u : v;
        revRID = edgeNum + 1; // if revRID is very large, then w(a,b) and w(b,a) would be different
        if (revRID > i)
        {
            // i'm edge(u,v), edge(v, u) is not initialized yet
            // if vEdgeInBG[i] already has
            inSpeed >> s;
            vector<string> vs = Tools::split(s, ",");
            vector<int> vX, vY;
            vX.reserve(288 / slotnum); // 1 day = 288 5-minutes
            vY.assign(288 / slotnum, -1);
            for (int j = 0; j < 288; j++)
            {
                if (j % slotnum == 0)
                {
                    vX.push_back(j * 300);
                    ss.clear();
                    ss.str("");
                    ss << vs[j];
                    ss >> vY[int(j / slotnum)];
                }
            }
            // vX = {600 * i, i = 0, 1, 2,..., 285}
            vEdge[i].vXFull = vX;
            vEdge[i].vYFull = vY;

//            vEdge[i].vX.assign(vX.begin(), vX.begin() + 5);
            vEdge[i].vX = {0, 900, 1800, 2700, 3600};
            vEdge[i].vY.assign(vY.begin(), vY.begin() + 5);
            vEdge[i].lpf = LPFunction(u, v, lBound, uBound, vEdge[i].vX, vEdge[i].vY, acc);

            for (int j = 0; j < vEdge[i].lpf.vX.size() - 1; j++)
            {
                vEdge[i].lpf.vSupportPre.push_back({{lVid, {j}}});
            }
        } else
        {
            // i'm edge(u,v), edge(v, u) is initialized
            vEdge[i].vXFull = vEdge[revRID].vXFull;
            vEdge[i].vYFull = vEdge[revRID].vYFull;

            vEdge[i].vX.assign(vEdge[revRID].vXFull.begin(), vEdge[revRID].vXFull.begin() + 5);
            vEdge[i].vY.assign(vEdge[revRID].vYFull.begin(), vEdge[revRID].vYFull.begin() + 5);
            vEdge[i].lpf = LPFunction(u, v, lBound, uBound, vEdge[revRID].vX, vEdge[revRID].vY, acc);

            for (int j = 0; j < vEdge[i].lpf.vX.size() - 1; j++)
            {
                map<int, vector<int>> mp = {{lVid, {j}}};
                vEdge[i].lpf.vSupportPre.push_back(mp);
            }
        }
    }
    inSpeed.close();
    cout << "Finish loading graph, nodeNum: " << nodeNum << endl;
}

int Graph::DijkstraPath(int ID1, int ID2, vector<int> &vPath, vector<int> &vPathEdge)
{
    benchmark::heap<2, int, int> queue(nodeNum);
    queue.update(ID1, 0);

    vector<int> vDistance(nodeNum, INF);
    vector<int> vPrevious(nodeNum, -1);
    vector<int> vPreviousEdge(nodeNum, -1);
    vector<bool> vbVisited(nodeNum, false);
    int topNodeID, neighborNodeID, neighborLength, neighborRoadID;

    vDistance[ID1] = 0;

    while (!queue.empty())
    {
        int topDistance;
        queue.extract_min(topNodeID, topDistance);
        vbVisited[topNodeID] = true;
        if (topNodeID == ID2)
            break;

        for (int i = 0; i < (int) adjList[topNodeID].size(); i++)
        {
            neighborNodeID = adjList[topNodeID][i].first;
            neighborLength = adjList[topNodeID][i].second;
            neighborRoadID = adjEdge[topNodeID][i].second;
            int d = vDistance[topNodeID] + neighborLength;
            if (!vbVisited[neighborNodeID])
            {
                if (vDistance[neighborNodeID] > d)
                {
                    vDistance[neighborNodeID] = d;
                    queue.update(neighborNodeID, d);
                    vPrevious[neighborNodeID] = topNodeID;
                    vPreviousEdge[neighborNodeID] = neighborRoadID;
                }
            }
        }
    }

    vPath.clear();
    vPathEdge.clear();
    vPath.push_back(ID2);
    int p = vPrevious[ID2];
    int e = vPreviousEdge[ID2];
    while (p != -1)
    {
        vPath.push_back(p);
        vPathEdge.push_back(e);
        e = vPreviousEdge[p];
        p = vPrevious[p];
    }

//	if(vPathEdge.size() > 0)
//		vPathEdge.erase(vPathEdge.end()-1);

    reverse(vPath.begin(), vPath.end());
    reverse(vPathEdge.begin(), vPathEdge.end());

    return vDistance[ID2];
}

int Graph::DijkstraPath2(int ID1, int ID2, unordered_set<int> &sRemovedNode, vector<int> &vPath, vector<int> &vPathEdge)
{
    benchmark::heap<2, int, int> queue(adjList.size());
    queue.update(ID1, 0);

    vector<int> vDistance(adjList.size(), INF);
    vector<int> vPrevious(adjList.size(), -1);
    vector<int> vPreviousEdge(adjList.size(), -1);
    vector<bool> vbVisited(adjList.size(), false);
    int topNodeID, neighborNodeID, neighborLength, neighborRoadID;

    vDistance[ID1] = 0;

    while (!queue.empty())
    {
        int topDistance;
        queue.extract_min(topNodeID, topDistance);
        vbVisited[topNodeID] = true;
        if (topNodeID == ID2)
            break;

        for (int i = 0; i < (int) adjList[topNodeID].size(); i++)
        {
            neighborNodeID = adjList[topNodeID][i].first;
            if (sRemovedNode.find(neighborNodeID) != sRemovedNode.end())
                continue;
            neighborLength = adjList[topNodeID][i].second;
            neighborRoadID = adjEdge[topNodeID][i].second;
            int d = vDistance[topNodeID] + neighborLength;
            if (!vbVisited[neighborNodeID])
            {
                if (vDistance[neighborNodeID] > d)
                {
                    vDistance[neighborNodeID] = d;
                    queue.update(neighborNodeID, d);
                    vPrevious[neighborNodeID] = topNodeID;
                    vPreviousEdge[neighborNodeID] = neighborRoadID;
                }
            }
        }
    }

    vPath.clear();
    vPathEdge.clear();
    vPath.push_back(ID2);
    int p = vPrevious[ID2];
    int e = vPreviousEdge[ID2];
    while (p != -1)
    {
        vPath.push_back(p);
        vPathEdge.push_back(e);
        e = vPreviousEdge[p];
        p = vPrevious[p];
    }

//	if(vPathEdge.size() > 0)
//		vPathEdge.erase(vPathEdge.end()-1);

    reverse(vPath.begin(), vPath.end());
    reverse(vPathEdge.begin(), vPathEdge.end());

    return vDistance[ID2];
}

void Graph::Dijkstra(int ID1, int ID2, int t, int &result)
{
    benchmark::heap<2, int, int> queue(nodeNum);
    queue.update(ID1, t);

    vector<int> vDistance(nodeNum, INF);
    vector<bool> vbVisited(nodeNum, false);
    int topNodeID, neighborNodeID, neighborLength;
    vector<pair<int, int> >::iterator ivp;

    vDistance[ID1] = t;
    while (!queue.empty())
    {
        int topDistance;
        queue.extract_min(topNodeID, topDistance);
        vbVisited[topNodeID] = true;
//		cout << topNodeID << "\t" << vDistance[topNodeID] << endl;
        if (topNodeID == ID2)
            break;

        for (ivp = adjEdge[topNodeID].begin(); ivp != adjEdge[topNodeID].end(); ivp++)
        {
            neighborNodeID = (*ivp).first;
            int edgeID = (*ivp).second;
            if (topDistance > uBound)
            {
                // out-of-bound

                result = -1;
                return;
            }
            neighborLength = vEdge[edgeID].lpf.getY(topDistance);
            int d = topDistance + neighborLength;
            if (!vbVisited[neighborNodeID])
            {
                if (vDistance[neighborNodeID] > d)
                {
                    vDistance[neighborNodeID] = d;
                    queue.update(neighborNodeID, d);
                }
            }
        }
    }
    result = vDistance[ID2] - t;
    return;
}

inline int Graph::EuclideanDistance(int ID1, int ID2)
{
    int lat = (int) (abs(vCoor[ID1].first - vCoor[ID2].first) * 111319);
    int lon = (int) (abs(vCoor[ID1].second - vCoor[ID2].second) * 83907);
    int min, max;
    min = (lat > lon) ? lon : lat;
    max = (lat > lon) ? lat : lon;
    int approx = max * 1007 + min * 441;
    if (max < (min << 4))
        approx -= max * 40;
    return (approx + 512) >> 10;
}

inline int Graph::EuclideanDistanceAdaptive(int ID1, int ID2, int latU, int lonU)
{
    int lat = (int) (abs(vCoor[ID1].first - vCoor[ID2].first) * latU);
    int lon = (int) (abs(vCoor[ID1].second - vCoor[ID2].second) * lonU);
    int min, max;
    min = (lat > lon) ? lon : lat;
    max = (lat > lon) ? lat : lon;
    int approx = max * 1007 + min * 441;
    if (max < (min << 4))
        approx -= max * 40;
    return (approx + 512) >> 10;
}

//Need to create 2hop only with functions
//Init 2hop to empty
//normal fastest path query
LPFunction Graph::forwardSearch(int ID1, int ID2, int &miny)
{
    benchmark::heap<2, int, int> queue(nodeNum + 1);
    vector<int> vTime(nodeNum + 1, INF);
    vector<bool> bVisited(nodeNum + 1, false);
    vector<bool> vbHeap(nodeNum + 1, false); //true if in Heap
    vector<bool> vbLPF(nodeNum + 1, false); // true if ever been reached
    LPFunction lpfInit;
    vector<LPFunction> vLPF(nodeNum + 1, lpfInit);

    int topNodeID, topNodeDist, neighborNodeID, neighborRoadID;
    int i, minCost;
    bool bUpdated;

    vTime[ID1] = lBound;
    vbHeap[ID1] = true;
    bVisited[ID1] = true;
    vbLPF[ID1] = true;
    queue.update(ID1, lBound);
    while (!queue.empty())
    {
        queue.extract_min(topNodeID, topNodeDist);
        vbHeap[topNodeID] = false;
        bVisited[topNodeID] = true;

        if (topNodeID == ID1)
        {
            LPFunction lpf;
            vLPF[topNodeID] = lpf;
        }

        if (topNodeID > nodeNum or topNodeID == ID2)
            break;

        for (i = 0; i < adjEdge[topNodeID].size(); i++)
        {
            bUpdated = false;
            neighborNodeID = adjEdge[topNodeID][i].first;
            neighborRoadID = adjEdge[topNodeID][i].second;

            LPFunction lpf;
            if (topNodeID == ID1)
            {
                lpf = vEdge[neighborRoadID].lpf; // lpf: ID1 -> x
                if (ID1 == vEdge[neighborRoadID].ID2)
                {
                    assert(false);
                    // lpf.ID1 == topNodeID == ID1 == ID2 cannot happen
                    lpf.ID1 = vEdge[neighborRoadID].ID2;
                    lpf.ID2 = vEdge[neighborRoadID].ID1;
                }
            } else
            {
                LPFunction lpftmp = vEdge[neighborRoadID].lpf;
                if (topNodeID == vEdge[neighborRoadID].ID2)
                {
                    // topNodeID == lpftmp.ID1 == lpftmp.ID1
                    assert(false);
                    lpftmp.ID1 = vEdge[neighborRoadID].ID2;
                    lpftmp.ID2 = vEdge[neighborRoadID].ID1;
                }
                lpf = vLPF[topNodeID].LPFCatSupport(lpftmp, lBound, uBound, 0); // lpf: u -> x
            }

            if (lpf.vX.size() == 1)
                continue;

            if (vbLPF[neighborNodeID])
            {
                //vLPF[neighborNodeID].display()
                //lpf.display();
                if (vLPF[neighborNodeID].vX.size() > 1 and vLPF[neighborNodeID].dominate(lpf, 5))
                    continue;
//                vLPF[neighborNodeID] = vLPF[neighborNodeID].LPFMinNew3(lpf);
                vLPF[neighborNodeID] = vLPF[neighborNodeID].LPFMinSupForDec(lpf, 0);
                bUpdated = true;
            } else
            {
                vLPF[neighborNodeID] = lpf;
                vbLPF[neighborNodeID] = true;
                bUpdated = true;
            }

            minCost = vLPF[neighborNodeID].minY;

            //Updated and not in Heap
            if ((!vbHeap[neighborNodeID] && bUpdated) || !bVisited[neighborNodeID])
            {
                vTime[neighborNodeID] = minCost;
                queue.update(neighborNodeID, minCost);
                vbHeap[neighborNodeID] = true;
                bVisited[neighborNodeID] = true;
            }
                //Updated and in Heap, key changed
            else if (bUpdated && vbHeap[neighborNodeID])
            {
                vTime[neighborNodeID] = minCost;
                queue.update(neighborNodeID, minCost);
            }
        }
    }
    miny = vLPF[ID2].minY;
    return vLPF[ID2];
}

void Graph::readODs(string filename, vector<pair<int, int>> &vOD, int nodeNum)
{
    ifstream in(filename.c_str());
    if (!in)
    {
        cout << "Cannot open file " << filename << endl;
        return;
    }

    int ID1, ID2;
    while (in >> ID1 >> ID2)
    {
        ID1 = ID1 % nodeNum;
        ID2 = ID2 % nodeNum;
        vOD.push_back(make_pair(ID1, ID2));
    }
    in.close();
}

void Graph::readUpdates(string filename, int nodeNum, vector<pair<pair<int, int>, pair<int, double>>> &testdata)
{
    int ID1, ID2, vX;
    double changeFrac;
    ifstream IF(filename);
    if (!IF)
    {
        cout << "Cannot open update file " << filename << endl;
        return;
    }
    while (IF >> ID1 >> ID2 >> vX >> changeFrac)
    {
        ID1 = ID1 % nodeNum;
        ID2 = ID2 % nodeNum;
        testdata.emplace_back(make_pair(ID1, ID2), make_pair(vX, changeFrac));
    }
    IF.close();
    cout << "Finish reading updates" << endl;
}

void Graph::readUpdates2(std::string filename, int nodeNum, int timeSlot,
                         double frac, vector<pair<pair<int, int>, pair<int, double>>> &testdata)
{
    int ID1, ID2;
    testdata.clear();
    ifstream IF(filename);
    if (!IF)
    {
        cout << "Cannot open update file " << filename << endl;
        return;
    }
    int cnt = 0;
    while (IF >> ID1 >> ID2)
    {
        ID1 = ID1 % nodeNum;
        ID2 = ID2 % nodeNum;
        testdata.emplace_back(make_pair(ID1, ID2), make_pair(timeSlot, frac));

        cnt++;
        if (cnt >= 1000)
            break;
    }
    IF.close();
    cout << "Finish reading " << testdata.size() << " updates at time " << timeSlot << endl;
}

void Graph::decUpdate(vector<pair<pair<int, int>, pair<int, double>>> &testdata, int timeWId, int timeWLen)
{
    cout << "Start dec updates of " << testdata.size() << " edge weights" << endl;
    for (const auto &k: testdata)
    {
        int a, b, vX, newW, oldW;
        a = k.first.first;
        b = k.first.second;
        vX = (k.second.first + timeWId) * timeWLen;
        double frac = k.second.second;
        for (const auto &i: adjEdge[a])
        {
            if (i.first == b)
            {
                LPFunction &lpf = vEdge[i.second].lpf;
                oldW = lpf.getY(vX);
                newW = ceil(frac * oldW);

                if (vX < lpf.vX[0] or vX > lpf.vX.back())
                    continue;

                vector<int>::const_iterator low;
                low = lower_bound(lpf.vX.begin(), lpf.vX.end(), vX);
                int pos = (int) (low - lpf.vX.begin());

                if (*low == vX)
                {
                    lpf.vY[pos] = newW;
                } else if (*low > vX)
                {
                    lpf.vX.insert(low, vX);
                    lpf.vY.insert(lpf.vY.begin() + pos, newW);
                    assert(lpf.vX[pos] == vX);
                }

                if (newW < lpf.minY)
                    lpf.minY = newW;

                vector<int> newVX = lpf.vX, newVY = lpf.vY;
                lpf.setValue(newVX, newVY, -1, 0);

                int lid = NodeOrder[a] < NodeOrder[b] ? a : b;
                lpf.vSupportPre.clear();
                for (int j = 0; j < lpf.vX.size() - 1; j++)
                {
                    lpf.vSupportPre.push_back({{lid, {j}}});
                }

                break;
            }
        }
    }
}

void Graph::incUpdate(vector<pair<pair<int, int>, pair<int, double>>> &testdata, int timeWId, int timeWLen)
{
    cout << "Start inc updates of " << testdata.size() << " edge weights" << endl;

    for (auto &k: testdata)
    {
        int a, b, vX, newW, pos;
        a = k.first.first;
        b = k.first.second;
        vX = (k.second.first + timeWId) * timeWLen;
        double frac = k.second.second;

        for (auto &i: adjEdge[a])
        {
            if (i.first == b)
            {
                LPFunction &lpf = vEdge[i.second].lpf;

                if (vX < lpf.vX[0] or vX > lpf.vX.back())
                    continue;

                int oldW = lpf.getY(vX);
                newW = oldW + int(frac * 10);

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
                    pos -= 1;
                }

                vector<int> newVX = lpf.vX, newVY = lpf.vY;
                lpf.setValue(newVX, newVY, -1, acc);

                int lid = NodeOrder[a] < NodeOrder[b] ? a : b;
                lpf.vSupportPre.clear();
                for (int j = 0; j < lpf.vX.size() - 1; j++)
                {
                    lpf.vSupportPre.push_back({{lid, {j}}});
                }
                break;
            }
        }
    }
}

void Graph::indexSize()
{
    int thisThreadNum = threadNum;
    int step = nodeNum / thisThreadNum;
    vector<pair<int, int>> intervals;
    for (int i = 0; i < thisThreadNum - 1; i++)
    {
        intervals.emplace_back(i * step, (i + 1) * step);
    }
    intervals.emplace_back((thisThreadNum - 1) * step, nodeNum);

    vector<int> scSize(thisThreadNum, 0), lbSize(thisThreadNum, 0);
    boost::thread_group threads;
    for (int i = 0; i < thisThreadNum; i++)
    {
        threads.add_thread(
                new boost::thread(&Graph::indexSizeThread, this,
                                  intervals[i].first, intervals[i].second,
                                  boost::ref(scSize[i]), boost::ref(lbSize[i])
                ));
    }
    threads.join_all();

    int scTtSize = 0, lbTtSize = 0;
    for (int i = 0; i < thisThreadNum; i++)
    {
        scTtSize += scSize[i];
        lbTtSize += lbSize[i];
    }

    cout << "Size of SC / LB index: " << scTtSize / 1000 << "k " << lbTtSize / 1000 << endl;
}

void Graph::indexSizeThread(int begin, int end, int &scSize, int &lbSize)
{
    for (int i = begin; i < end; i++)
    {
        for (auto &p: Tree[i].vertIn)
        {
            scSize += p.second.first.vX.size();
        }

        for (auto &p: Tree[i].vertOut)
        {
            scSize += p.second.first.vX.size();
        }

        for (auto &p: Tree[i].disIn)
        {
            lbSize += p.vX.size();
        }

        for (auto &p: Tree[i].disOut)
        {
            lbSize += p.vX.size();
        }
    }
}



