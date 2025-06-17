#include "DHL.h"


bool equalSupInfo(const LPFunction &f1, const LPFunction &f2)
{
    set<int> s1;
    for (const auto &p: f1.vSupportPre)
    {
        for (const auto &pp: p)
        {
            s1.insert(pp.first);
        }
    }

    set<int> s2;
    for (const auto &p: f2.vSupportPre)
    {
        for (const auto &pp: p)
        {
            s2.insert(pp.first);
        }
    }

    if (s1.size() != s2.size())
    {
        return false;
    }

    for (const auto &i: s1)
    {
        if (s2.find(i) == s2.end())
        {
            return false;
        }
    }

    for (const auto &i: s2)
    {
        if (s1.find(i) == s1.end())
        {
            return false;
        }
    }

    return true;

}

long long int totalTuningPointSum(Graph &g)
{
    long long int dinInSum = 0, disOutSum = 0, vertInSize = 0, vertOutSize = 0;
    int lastPoint = INF;
    for (int i = 0; i < g.nodeNum; i++)
    {
        if (i >= g.Tree.size())
            continue;

        for (const auto &p: g.Tree[i].disIn)
        {
            dinInSum += p.vX.size();
            if (p.vX.empty())
                lastPoint = g.lBound;
            else if (p.vX.back() < lastPoint)
                lastPoint = p.vX.back();
        }

        for (const auto &p: g.Tree[i].disOut)
        {
            disOutSum += p.vX.size();
            if (p.vX.empty())
                lastPoint = g.lBound;
            else if (p.vX.back() < lastPoint)
                lastPoint = p.vX.back();
        }

        for (const auto &p: g.Tree[i].vertIn)
        {
            vertInSize += p.second.first.vX.size();
        }

        for (const auto &p: g.Tree[i].vertOut)
        {
            vertOutSize += p.second.first.vX.size();
        }
    }

    int sizeInt = sizeof(int);
    long long int indexSize = sizeInt * (dinInSum + disOutSum + vertInSize + vertOutSize);
//    cout << "dinInSum: " << dinInSum << " disOutSum: " << disOutSum << " vertInSize: " << vertInSize << " vertOutSize: "
//         << vertOutSize << endl;
    cout << "Earliest upperbound: " << lastPoint << " Index size: " << indexSize * 0.000000954 << " mb" << endl;
    return indexSize * 0.000000954;

}

bool compareSCs(vector<pair<int, pair<LPFunction, int>>> &v1, vector<pair<int, pair<LPFunction, int>>> &v2)
{
//    cout << "**Check scs**" << endl;
    bool same = true;
    if (v1.size() != v2.size())
    {
        cout << "error 1" << endl;
        return false;
    }
//    cout << "sc num: " << v1.size() << endl;

    for (const auto &i: v1)
    {
        for (const auto &j: v2)
        {
            if (i.first == j.first)
            {
                if (!i.second.first.equal(j.second.first))
//                if (!equalSupInfo(i.second.first, j.second.first))
                {
                    cout << "**************************" << endl;
                    i.second.first.display();
                    j.second.first.display();
                    same = false;
                }
                break;
            }
        }
    }
//    cout << "sc check done" << endl;
    return same;
}

bool compareLabels(vector<LPFunction> &labs1, vector<LPFunction> &labs2)
{
    if (labs1.size() != labs2.size())
    {
        return false;
    }

    for (int i = 0; i < labs1.size(); i++)
    {
        if (!labs1[i].equal(labs2[i]))
        {
            return false;
        }
    }

    return true;
}

//int checkIndexMaintain(Graph &g1, Graph &g2, int noQueries, int accDev)
//{
//    srand(time(nullptr));
//    int o, d, nbQ = 0, failCnt = 0;
//    unordered_map<int, unordered_set<int>> ods;
//    while (nbQ < noQueries)
//    {
//        o = rand() % g1.nodeNum;
//        d = rand() % g1.nodeNum;
//        if (o == d or (ods.find(o) != ods.end() and ods[o].find(d) != ods[o].end()))
//            continue;
//
//        nbQ += 1;
//
//        if (ods.find(o) == ods.end())
//        {
//            unordered_set<int> tmp;
//            tmp.insert(d);
//            ods[o] = tmp;
//        } else
//        {
//            ods[o].insert(d);
//        }
//
//        int test;
//        long long int lcaSize = 0;
//        LPFunction lpf = g1.forwardSearch(o, d, test);
//        LPFunction lpf2 = g2.QueryH2H(o, d, lcaSize);
//
//        int def = abs(lpf.minY - lpf2.minY);
//        cout << lpf.minY << " " << lpf2.minY << " " << def << endl;
//        if (def > accDev)
//        {
//            failCnt += 1;
//        }
//    }
//    return failCnt;
//}

void findFirstIncorrectSCs(const Graph &g1, const Graph &g2)
{
    cout << "Find first incorrect SCs" << endl;
    vector<TreeNode> tr1 = g1.Tree, tr2 = g2.Tree;
    vector<int> leaves;
    stack<int> q;
    q.push(0);
    while (!q.empty())
    {
        int uR = q.top();
        q.pop();

        for (const auto &cR: tr1[uR].vChildren)
        {
            q.push(cR);
        }

        bool inSm = compareSCs(tr1[uR].vertIn, tr2[uR].vertIn);
        bool outSm = compareSCs(tr1[uR].vertOut, tr2[uR].vertOut);

        if (!inSm or !outSm)
        {
            leaves.emplace_back(uR);
        }
    }

    sort(leaves.begin(), leaves.end(), greater<int>());
    cout << leaves.size() << " incorrect nodes" << endl;

    int i = 0;
    for (int leave: leaves)
    {
        cout << g1.Tree[leave].uniqueVertex << "\t";
        i++;
        if (i > 100)
            break;
    }
    cout << endl;
}

void findFirstIncorrectLabs(Graph &g1, Graph &g2)
{
    cout << "Find first incorrect Labels" << endl;
    vector<TreeNode> tr1 = g1.Tree, tr2 = g2.Tree;
    vector<int> leaves;
    queue<int> q;
    q.push(0);
    int debugID1 = 4913, debugID2 = 1111;
    while (!q.empty())
    {
        int uR = q.front();
        q.pop();

        bool inSm = compareLabels(tr1[uR].disIn, tr2[uR].disIn);
        bool outSm = compareLabels(tr1[uR].disOut, tr2[uR].disOut);

        if (!inSm or !outSm)
        {
            leaves.emplace_back(uR);
        }

        for (const auto &cR: tr1[uR].vChildren)
        {
            q.push(cR);
        }
    }
    sort(leaves.begin(), leaves.end());

    cout << leaves.size() << " incorrect nodes" << endl;
    int printS = 100000000;
    for (const auto &vRank: leaves)
    {
        int j = 0;
        for (int i = 0; i < g1.Tree[vRank].disIn.size(); i++)
        {
            LPFunction l1 = g1.Tree[vRank].disIn[i];
            LPFunction l2 = g2.Tree[vRank].disIn[i];
            if (l1.ID1 == debugID1 or l1.ID2 == debugID2)
            {
                if (!l1.equal(l2))
                {
                    cout << "~~~~~~~~~~~~~~~~~~~~~~~~ in ~~~~~~~~~~~~~~~~~~~~~~~~~\n";
                    l1.display();
                    l2.display();
                    j += 1;
                }
            }
        }
        for (int i = 0; i < g1.Tree[vRank].disOut.size(); i++)
        {
            LPFunction l1 = g1.Tree[vRank].disOut[i];
            LPFunction l2 = g2.Tree[vRank].disOut[i];
            if (l1.ID1 == debugID1 or l1.ID2 == debugID2)
            {
                if (!l1.equal(l2))
                {
                    cout << "~~~~~~~~~~~~~~~~~~~~~~~~ out ~~~~~~~~~~~~~~~~~~~~~~~~~\n";
                    l1.display();
                    l2.display();

                    j += 1;
                }
            }
        }
//        if (j > printS)
//            break;
    }
}

void testTDGTree(Graph &g, string &DHLPath, vector<string> &queryFiles, vector<string> &updateFiles)
{
    Timer clock;
    TDGTree tdg = TDGTree(g);
    tdg.ReadGTree(DHLPath);

    clock.tick();
    tdg.TwoPMatrixBuild();
    clock.tock();
    cout << "TD-GTree construction time: " << clock.duration().count() << endl;

    long long int turningPointSum = 0;
    for (const auto &nd: tdg.HGPTree)
    {
        for (const auto &p: nd.dynMat)
        {
            for (const auto &pp: p)
            {
                turningPointSum += pp.second.vX.size();
            }
        }
    }

    long long int indexSize = sizeof(int) * turningPointSum;
    cout << "Index size: " << indexSize * 0.000000954 << endl;


    for (string &queryFile: queryFiles)
    {
        cout << "Query set: " << queryFile << endl;
        vector<pair<int, int>> ODs;
        Graph::readODs(queryFile, ODs, g.nodeNum);

        int cnt = 0;
        clock.tick();
        for (const auto &p: ODs)
        {
            tdg.TIPQuery(p.first, p.second);
            cnt += 1;
            if (cnt % 100 == 0)
                cout << cnt << endl;
        }
        clock.tock();
        cout << "Avg Flex Query time: " << 1.0 * clock.duration().count() / ODs.size() << " ms" << endl;

        cnt = 0;
        clock.tick();
        for (const auto &p: ODs)
        {
            tdg.TDSPQuery(p.first, p.second, g.lBound);
            cnt += 1;
            if (cnt % 100 == 0)
                cout << cnt << endl;
        }
        clock.tock();
        cout << "Avg Fix query time: " << 1.0 * clock.duration().count() / ODs.size() << " ms" << endl;
    }


    for (const auto &updateFile: updateFiles)
    {
        vector<float> fractions = {1, 0.75, 0.5, 0.25};
        vector<pair<int, int>> testdata;
        Graph::readODs(updateFile, testdata, g.nodeNum);
        for (float &f: fractions)
        {
            clock.tick();
            int t = (g.uBound + 600) * f - 600;
            tdg.Update(testdata, t);
            clock.tock();
            cout << "UpdSize: " << testdata.size() << " updT: " << t << " spent time: " << clock.duration().count()
                 << endl;
        }
    }

    cout << "tdg check done" << endl;
}

void testDHLTree(Graph &g, string &DHLPath, vector<string> &querFiles, vector<string> &updateFiles)
{

//    Graph::readUpdates2(updateFile, g.nodeNum, changedvx, frac, testdata);
    Timer clock;

    DHL dhl = DHL(g);
    dhl.ReadGTree(DHLPath);
    clock.tick();
    dhl.BuildDHLIndex();
    clock.tock();
    cout << "DHL construction time: " << clock.duration().count() << endl;

    long long int turningPointSum = 0;
    for (const auto &nd: dhl.HGPTree)
    {
        for (const auto &p: nd.dynMat)
        {
            for (const auto &pp: p)
            {
                turningPointSum += pp.second.vX.size();
            }
        }
    }


    turningPointSum += totalTuningPointSum(dhl.borderGraph);

//    long long int indexSize = sizeof(int) * turningPointSum;
//    cout << "Index size: " << indexSize * 0.000000954 << " mb" << endl;

//    string queryFile = destFile + "queryODs.txt";
//    vector<pair<int, int>> ODs;
//    Graph::readODs(queryFile, ODs, g.nodeNum);
//    clock.tick();
//    int cnt = 0;
//    for (const auto &p: ODs)
//    {
//        dhl.DHLQuery(p.first, p.second);
//        cnt += 1;
//        if (cnt % 10 == 0)
//            cout << cnt << endl;
//    }
//    clock.tock();
//    cout << " Query time: " << clock.duration().count() << endl;

    for (const auto &updateFile: updateFiles)
    {
        vector<pair<int, int>> testdata;
        Graph::readODs(updateFile, testdata, g.nodeNum);

        vector<float> fractions = {1, 0.75, 0.5, 0.25};
        for (float &f: fractions)
        {
            clock.tick();
            int t = (g.uBound + 600) * f - 600;
            dhl.Update(testdata, t);
            clock.tock();
            cout << "UpdSize: " << testdata.size() << " updT: " << t << "; spent time: " << clock.duration().count()
                 << endl;
        }
    }

//    for (const auto &p: testdata)
//    {
//        int u = p.first.first;
//        int v = p.first.second;
//        cout << dhl.TIPQuery(u, v).minY << endl;
//    }

    cout << "dhl check done" << endl;
}

void testH2H(Graph &g, vector<string> &queryFiles, vector<string> &updateFiles)
{
    Timer clock;

    clock.tick();
    g.loadTD();
    clock.tock();
    cout << "SWTD build time: " << clock.duration().count() << endl;
//    clock.tick();
//    g.extendIndex();
//    clock.tock();
//    cout << "SWTD extend time: " << clock.duration().count() << endl;
    totalTuningPointSum(g);

//    string queryFile = destFile + "queryODs.txt";
    for (auto &queryFile: queryFiles)
    {
        vector<pair<int, int>> testdata;
        Graph::readODs(queryFile, testdata, g.nodeNum);
        if (testdata.empty())
            continue;
        cout << "query set: " << queryFile << " " << testdata.size() << " queries" << endl;
        long long int lcaSizeFlex = 0, lcaSizeFix = 0, flexT = 0, fixT = 0;
        using namespace std::chrono;
        clock.tick();
        for (const auto &p: testdata)
        {
            int o = p.first;
            int d = p.second;
            g.QueryH2H(o, d, lcaSizeFlex);
        }
        clock.tock();
        cout << "Flex query time (Avg): " << 1.0 * clock.duration().count() / testdata.size() << " avg LCA: "
             << lcaSizeFlex / testdata.size() << endl;

        clock.tick();
        for (const auto &p: testdata)
        {
            int o = p.first;
            int d = p.second;

            g.QueryH2HFixedT(o, d, g.lBound, lcaSizeFix);
        }
        clock.tock();
        cout << "Fix query time (Avg): " << 1.0 * clock.duration().count() / testdata.size() << " avg LCA: "
             << lcaSizeFix / testdata.size() << endl;
    }
    for (const string &updateFile: updateFiles)
    {
        vector<pair<int, int>> updateOD;
        Graph::readODs(updateFile, updateOD, g.nodeNum);

        vector<double> fractions = {1};
        for (double &f: fractions)
        {
            clock.tick();
            g.H2HIncBatDiGraph(updateOD, (g.uBound + 600) * f - 600);
            clock.tock();
            cout << "Update time: " << clock.duration().count() << endl;
        }
    }
}

void testCH(Graph &g, string &map, vector<string> &queryFiles, vector<string> &updateFiles)
{
    Timer clock;

    clock.tick();
    g.CHConstruction();
    g.makeTree();
    clock.tock();
    cout << "CH build time: " << clock.duration().count() << endl;
//    clock.tick();
//    g.extendIndex();
//    clock.tock();
//    cout << "SWTD extend time: " << clock.duration().count() << endl;
    totalTuningPointSum(g);

//    string queryFile = destFile + "queryODs.txt";
    for (auto &queryFile: queryFiles)
    {
        vector<pair<int, int>> testdata;
        Graph::readODs(queryFile, testdata, g.nodeNum);
        cout << "query set: " << queryFile << " " << testdata.size() << " queries" << endl;
        long long int lcaSizeFlex = 0, lcaSizeFix = 0, flexT = 0, fixT = 0;
        using namespace std::chrono;
        clock.tick();
        for (const auto &p: testdata)
        {
            int o = p.first;
            int d = p.second;
            int x;
            g.QueryCHItvT(o, d, x);
        }
        clock.tock();
        cout << "Flex query time (Avg): " << 1.0 * clock.duration().count() / testdata.size() << " avg LCA: "
             << lcaSizeFlex / testdata.size() << endl;

        clock.tick();
        for (const auto &p: testdata)
        {
            int o = p.first;
            int d = p.second;
            g.QueryCHFixedT(o, d, g.lBound);
        }
        clock.tock();
        cout << "Fix query time (Avg): " << 1.0 * clock.duration().count() / testdata.size() << " avg LCA: "
             << lcaSizeFix / testdata.size() << endl;
    }
    for (string updateFile: updateFiles)
    {
        vector<pair<int, int>> updateOD;
        Graph::readODs(updateFile, updateOD, g.nodeNum);
        clock.tick();
        vector<double> fractions = {1, 0.75, 0.5, 0.25};
        for (double &f: fractions)
        {
            g.CHUpdate(updateOD, (g.uBound + 600) * f - 600);
            clock.tock();
            cout << "Update time: " << clock.duration().count() << endl;
        }
    }
}

vector<int> sumError(Graph &g1, Graph &g2)
{
    // summarize error in labels
    cout << "Summarize error in labels" << endl;
    vector<int> error;
    error.reserve(g1.nodeNum * g1.Tree[g1.rank.size() - 1].disIn.size() * 2);

    int maxError = 0;
    LPFunction lpf1, lpf2;
    int totalLabNo = 0;
    for (int i = 0; i < g1.nodeNum; i++)
    {
        totalLabNo += g1.Tree[i].disIn.size();
        for (int j = 0; j < g1.Tree[i].disIn.size(); j++)
        {
            LPFunction l1 = g1.Tree[i].disIn[j];
            LPFunction l2 = g2.Tree[i].disIn[j];
            int def = l1.equalValue2(l2);
            if (def > 0)
                error.emplace_back(def);
            if (def > maxError)
            {
                maxError = def;
                lpf1 = l1;
                lpf2 = l2;
            }
        }
    }

    for (int i = 0; i < g1.nodeNum; i++)
    {
        totalLabNo += g1.Tree[i].disOut.size();
        for (int j = 0; j < g1.Tree[i].disOut.size(); j++)
        {
            LPFunction l1 = g1.Tree[i].disOut[j];
            LPFunction l2 = g2.Tree[i].disOut[j];
            int def = l1.equalValue2(l2);
            if (def > 0)
                error.emplace_back(def);
            if (def > maxError)
            {
                maxError = def;
                lpf1 = l1;
                lpf2 = l2;
            }
        }
    }

    sort(error.begin(), error.end(), greater<int>());
    cout << "max dev: " << maxError << endl;
    cout << "total labels: " << totalLabNo << " prob labels: " << error.size() << endl;

    cout << endl;
    if (maxError > 0)
    {
        lpf1.display();
        lpf2.display();
    }
//    lpf1.equalValue2(lpf2);
    return error;
}

void checkLPF()
{
    vector<int> vX1 = {0, 20, 25, 100, 110};
    vector<int> vY1 = {5, 20, 80, 60, 60};

    vector<int> vX2 = {0, 15, 100, 110};
    vector<int> vY2 = {10, 75, 25, 25};

    LPFunction lpf(1, 2, 0, 100, vX1, vY1);
    for (int i = 0; i < vX1.size() - 1; i++)
        lpf.vSupportPre.push_back({{i, {i}}});
    lpf.cntOfEachInt.assign(vX1.size() - 1, 1);

    lpf.LPFTruncate(0);

    LPFunction lpf2(2, 3, 0, 100, vX2, vY2);
    for (int i = 0; i < vX2.size() - 1; i++)
        lpf2.vSupportPre.push_back({{1, {i}}});

    int test = lpf.LPFGetX(lpf2.vX.back());
//    vector<unordered_map<int, unordered_map<int, CatSupRec>>> intermediateLPF;
    CatSupRec catSupRec;
    LPFunction lpf3 = lpf.LPFCatSupport(lpf2, catSupRec);
    cout << 1;


//    vector<int> v1, v2, v3;
//    lpf.LPFMinSupForDec(lpf2);

//    lpf.equalValue(lpf2);
//    lpf.lastItvSubToChange = 4;
//    lpf2.lastItvSubToChange = 0;
//    CatSupRec catSupRec;
    catSupRec.vX1 = {0, 3600};
    catSupRec.vX2 = catSupRec.vX1;
    catSupRec.vY = {249, 249};
    unordered_map<int, unordered_map<int, CatSupRec>> intermediateLPF2 = {{1, {{2, catSupRec}}}};
    lpf3.vSupportLastPartContains(1, intermediateLPF2);
//    lpf.LPFCatSupportExtend(lpf2, 0, 3600, intermediateLPF);
//    LPFunction f = lpf.LPFMinSupByPlaneSweepNoSV(lpf2);

}

void checkCatRec(int u, int w, int v, vector<unordered_map<int, unordered_map<int, CatSupRec>>> &intermediates)
{
    cout << w << endl;
    if (intermediates[w].find(u) == intermediates[w].end())
    {
        cout << "not found1" << endl;
    }

    if (intermediates[w][u].find(v) == intermediates[w][u].end())
    {
        cout << "not found2" << endl;
    }

    intermediates[w][u][v].displayCatRec();
}

Graph extendGraph(Graph &g, int winNum1, int winNum2, int deltaWin, int acc)
{
    Timer clock;
    for (auto &p: g.vEdge)
    {
        vector<int> newX, newY;
        newX.assign(p.vXFull.begin() + winNum1, p.vXFull.begin() + winNum2 + deltaWin);
        newY.assign(p.vYFull.begin() + winNum1, p.vYFull.begin() + winNum2);
        for (int i = 0; i < deltaWin; i++)
        {
            newY.push_back(p.lpf.vY.back());
        }
        int lid = g.NodeOrder[p.lpf.ID1] < g.NodeOrder[p.lpf.ID2] ? p.lpf.ID1 : p.lpf.ID2;
        LPFunction lpf(p.lpf.ID1, p.lpf.ID2, winNum1 * g.itvLen,
                       (winNum2 - 1) * g.itvLen, newX, newY);
        lpf.vSupportPre.clear();
        for (int i = 0; i < lpf.vX.size() - 1; i++)
            lpf.vSupportPre.push_back({{lid, {i}}});
        p.lpf = lpf;
    }
    clock.tick();
    g.loadTD();
    clock.tock();
    cout << "Construction time: " << clock.duration().count() << "\n" << endl;
    return g;
}

void contSliding(string &graphFile, string &filenameSpeed, int lBound, int uBound, int deltaWinNum, int acc)
{
    vector<pair<int, int>> testdata;
    vector<int> switchTimes, updTimes, updSize;

    int iterNum = 0;
    while (uBound < 86400)
    {
        Graph g = Graph();
        g.readDirectedGraph(graphFile, filenameSpeed, lBound, uBound, deltaWinNum);

        Timer clock;
        clock.tick();
        g.loadTD();
        clock.tock();
        cout << "*** Iteration: " << iterNum << ", lBound: " << lBound << ", uBound: " << uBound << endl;
        cout << "   Construction time: " << clock.duration().count() << "\n" << endl;

        g.readExtension(testdata, uBound);
        updSize.emplace_back(testdata.size());
        cout << "    Update size: " << testdata.size() << endl;

        clock.tick();
//        g.H2HExtendBatDiGraph(testdata);
        clock.tock();
        int appendUpdT = clock.duration().count();
        cout << "    Inc update time: " << appendUpdT << " ms" << endl;
        updTimes.push_back(appendUpdT);

        lBound += 300 * deltaWinNum;
        uBound += 300 * deltaWinNum;

        iterNum++;
    }

    cout << "Update sizes: ";
    for (int i: updSize)
    {
        cout << i << ",";
    }
    cout << endl;

    cout << "Update times: ";
    for (int updTime: updTimes)
    {
        cout << updTime << ",";
    }
    cout << endl;

    cout << "Switch times: ";
    for (int switchTime: switchTimes)
    {
        cout << switchTime << ",";
    }
    cout << endl;

}

void contReconstruction(string &graphFile, string &filenameSpeed, int winNum1, int winNum2, int deltaWin, int acc)
{
    int iterNum = 0;

    vector<pair<pair<int, int>, pair<int, double>>> testdata;
    vector<int> constructT, switchTimes, updTimes, updsize;

    Timer clock;
    while (winNum1 >= 0)
    {
        cout << "*** Iteration: " << iterNum << ", lBound: " << winNum1 << ", uBound: " << winNum2 << endl;

        Graph g = Graph();
        g.readDirectedGraph(graphFile, filenameSpeed, winNum1, winNum2, deltaWin);
        clock.tick();
        g.loadTD();
        clock.tock();
        int constT = clock.duration().count();
        cout << "    Construction time: " << constT << endl;
//        constructT.push_back(clock.duration().count());
//        totalTuningPointSum(g);
//        return;
        winNum1 += 300 * deltaWin;
        winNum2 += 300 * deltaWin;

//        clock.tick();
//        g.extendIndex();
//        clock.tock();
//        cout << "    Extending idx time: " << clock.duration().count() << " ms" << endl;
//        int extTime = clock.duration().count();
//        switchTimes.push_back(extTime);

//        g.readExtension(testdata, winNum2 / 300, acc);
//        updsize.emplace_back(testdata.size());
//        cout << "    Update size: " << testdata.size() << endl;
//
//        Graph g2 = g;
//        clock.tick();
//        g.H2HExtendBatDiGraph(testdata);
//        clock.tock();
//        int appendUpdT = clock.duration().count();
//        cout << "    Inc update time: " << appendUpdT << " ms" << endl;
//        updTimes.push_back(appendUpdT);
//
//        cout << "SWTD slide eff from " << winNum1 - 300 * deltaWin << " to " << winNum2 << " " << constT << " "
//             << testdata.size() << " " << extTime << " " << appendUpdT << endl;
//        sumError(g, g2);

        iterNum++;
    }

//    cout << "Construction times: ";
//    for (int i: constructT)
//    {
//        cout << i << ",";
//    }
//    cout << endl;
//
    cout << "Update sizes: ";
    for (int i: updsize)
    {
        cout << i << ",";
    }
    cout << endl;
//
//    cout << "Update times: ";
//    for (int updTime: updTimes)
//    {
//        cout << updTime << ",";
//    }
//    cout << endl;
//
//    cout << "Switch times: ";
//    for (int switchTime: switchTimes)
//    {
//        cout << switchTime << ",";
//    }
//    cout << endl;

}

void testAcc(string graphFile, string filenameSpeed, string DHLPath, string queryFile, int winNum1, int winNum2)
{
    Graph g2 = Graph();
    g2.readDirectedGraph(graphFile, filenameSpeed, winNum1, winNum2, 1);
    for (int cc: {5})
    {
        for (int deltaT: {6, 3, 2, 1})
        {
            Graph g = Graph();
            g.readDirectedGraph(graphFile, filenameSpeed, winNum1, winNum2, deltaT);
            vector<pair<int, int>> testdata;
            Graph::readODs(queryFile, testdata, g.nodeNum);

            Timer clock;
            clock.tick();
            g.loadTD();
            clock.tock();
            cout << "~~~~~~~~~ SWTD Time: " << clock.duration().count() << "\n";

            clock.tick();
            TDGTree tdg = TDGTree(g);
            tdg.ReadGTree(DHLPath);
            tdg.TwoPMatrixBuild();
            clock.tock();
            cout << "TDG Time: " << clock.duration().count() << " ~~~~~~~~~" << endl;

            double devSumCH = 0, devSumTDG = 0, devSumH2H = 0;
            vector<int> dijkResultsFix(testdata.size(), -1);
            for (auto &p: testdata)
            {

                int id1 = p.first;
                int id2 = p.second;

                std::random_device rd; // obtain a random number from hardware
                std::mt19937 gen(rd()); // seed the generator
                std::uniform_int_distribution<> distr(0, 1800);
                int departT = distr(gen);

                int dijResult = 0;
                g2.Dijkstra(id1, id2, departT, dijResult);
                int tdgResult = tdg.TDSPQuery(id1, id2, departT);
                if (tdgResult == INF)
                    continue;
                devSumTDG += 1.0 * abs(tdgResult - dijResult) / dijResult * 100.0;

//                LPFunction lpf = g.forwardSearch(id1, id2, dijResult2);
                int chResult = g.QueryCHFixedT(id1, id2, departT);
                devSumCH += 1.0 * abs(chResult - dijResult) / dijResult * 100.0;

                long long int lca = 0;
                int h2hResult = g.QueryH2HFixedT(id1, id2, departT, lca);
                devSumH2H += 1.0 * abs(h2hResult - dijResult) / dijResult * 100.0;

//                cout << dijResult << " " << chResult << " " << tdgResult << " " << h2hResult << endl;
            }
            devSumCH = 1.0 * devSumCH / testdata.size();
            devSumTDG = 1.0 * devSumTDG / testdata.size();
            devSumH2H = 1.0 * devSumH2H / testdata.size();
            cout << cc << " " << deltaT << " " << devSumCH << " " << devSumTDG << " " << devSumH2H << endl;
//            testAcc(graphFile, filenameSpeed, DHLPath, queryFile);
        }
    }
}

int main()
{
    checkLPF();
    return 1;
    srand(time(nullptr));

    std::chrono::high_resolution_clock::time_point t1, t2;
    std::chrono::duration<double> time_span{}, total_span{};
    vector<int> vPath, vPathEdge;

    string exampleGraph = "../Data/ExampleGraph";
    string exampleSpeed = "../Data/ExampleSpeedProfile";
    Graph g;
    g.readExampleGraph(exampleGraph, exampleSpeed, 0, 100);
    g.loadTD();

    return 1;

    string destFile = "./Data/";
    string filenameSpeed = destFile + "SpeedProfile_Bangkok.txt";
//    string filenameSpeed = destFile + "SpeedProfile_xj.txt";
    string map = "Bangkok";
    string graphFile = destFile + map + "Conn.txt";
    string DHLPath = destFile + map + "_gtree64.txt";
    string queryFile = destFile + "queryODs.txt";
    vector<string> updateFiles =
            {
                    destFile + map + "UpdFr0.1.txt",
                    destFile + map + "UpdFr1.txt",
                    destFile + map + "UpdFr10.txt",
                    destFile + map + "UpdFr20.txt"
            };


//    contSliding(graphFile, filenameSpeed, 0, 3600, 2, 5);
//    contReconstruction(graphFile, filenameSpeed, 0, 6600, 2, 5);
//    return 1;

    for (int i = 0; i <= 132; i++)
    {
//        int deltaWin = 2, slotNum = deltaWin, acc = 5;
//        int winNum1 = i * 300 * deltaWin, winNum2 = 7200 + winNum1;
//        Graph g = Graph(deltaWin, acc);
//        g.readDirectedGraph(graphFile, filenameSpeed, winNum1, winNum2, slotNum);
//        Timer clock;
//        clock.tick();
//        g.CHConstruction();
//        g.makeTree();
//        clock.tock();
//        cout << "CH build time from " << winNum1 << " to " << winNum2 << ": " << clock.duration().count() << endl;
//
//        Timer clock;
//        DHL dhl = DHL(g, 5);
//        dhl.ReadGTree(DHLPath);
//        clock.tick();
//        dhl.BuildDHLIndex();
//        clock.tock();
//        cout << "DHL build time from " << winNum1 << " to " << winNum2 << ": " << clock.duration().count() << endl;

//        Timer clock;
//        TDGTree tdg = TDGTree(g, 5);
//        tdg.ReadGTree(DHLPath);
//        clock.tick();
//        tdg.TwoPMatrixBuild();
//        clock.tock();
//        cout << "TD-GTree build time from " << winNum1 << " to " << winNum2 << ": " << clock.duration().count() << endl;
    }

    for (int length: {2, 4, 8, 12})
    {
        Timer clock;
        int deltaWin = 2, slotNum = deltaWin;
        int winNum1 = 0, winNum2 = 3600 * length + winNum1 - deltaWin * 300;

        Graph g = Graph();
        g.readDirectedGraph(graphFile, filenameSpeed, winNum1, winNum2, slotNum);
//
        clock.tick();
        g.loadTD();
        clock.tock();
        cout << "SWTD build time: " << clock.duration().count() << endl;

        clock.tick();
        vector<TreeNode> tree = g.precomputeSCForA();
        g.makeIndexForA(tree);
        clock.tock();
        cout << "MakeIndexDFSForA time: " << clock.duration().count() << endl;
        continue;

        vector<string> queryFiles =
                {
//                        destFile + map + "_queries_hop10.txt",
//                        destFile + map + "_queries_hop25.txt",
//                        destFile + map + "_queries_hop50.txt",
//                        destFile + map + "_queries_hop100.txt",
//                        destFile + map + "_queries_hop200.txt"
                };

//        for (auto &qf: queryFiles)
//        {
//            cout << "Query set: " << qf << endl;
//            vector<pair<int, int>> testdata;
//            Graph::readODs(qf, testdata, g.nodeNum);
//
//            int lcaSizeFlex = 0, lcaSizeFix = 0;
//            long long int flexT = 0, fixT = 0;
//            using namespace std::chrono;
//            for (const auto &p: testdata)
//            {
//                int o = p.first;
//                int d = p.second;
//                clock.tick();
//                g.forwardSearch(o, d, lcaSizeFlex);
//                clock.tock();
//                flexT += clock.duration().count();
//            }
//            cout << "Avg Flex query time (Avg): " << 1.0 * flexT / testdata.size() << " ms" << endl;
//
//            for (const auto &p: testdata)
//            {
//                int o = p.first;
//                int d = p.second;
//                clock.tick();
//                g.Dijkstra(o, d, g.lBound, lcaSizeFix);
//                clock.tock();
//                fixT += clock.duration().count();
//            }
//            cout << "Avg Fix query time (Avg): " << 1.0 * fixT / testdata.size() << " ms" << endl;
//        }
//        testDHLTree(g, DHLPath, queryFiles, updateFiles);
//        testTDGTree(g, DHLPath, queryFiles, updateFiles);
//        testH2H(g, queryFiles, updateFiles);
//        testCH(g, map, queryFiles, updateFiles);
    }

    cout << "done" << endl;
    return 1;
}
	
