#ifndef LINEARPIECEWISEFUNCTION_H
#define LINEARPIECEWISEFUNCTION_H

#include <vector>
#include <list>
#include <stdexcept>
#include <iostream>
#include <algorithm>
#include <map>
#include <set>
#include <math.h>
#include <string>
#include <sstream>
#include <queue>
#include <unordered_map>
#include <unordered_set>
//#include <google/dense_hash_map>
#include "heap.h"
//#define INF 99999999
#include "tools.h"
#include "Semaphore.h"
#include <boost/functional/hash.hpp>

using namespace std;

struct CatSupRec
{
    vector<int> vX1, vX2, vY;

    void displayCatRec()
    {
        cout << "vX: ";
        for (const auto &x: vX1)
        {
            cout << x << " ";
        }
        cout << endl;

        cout << "vY: ";
        for (const auto &y: vY)
        {
            cout << y << " ";
        }
        cout << endl;
    }

//    void extendEnd(int ubound, int itv)
//    {
//        if (vX1.size() < 2)
//            return;
//        if (false and Tools::redundant(*(vX1.end() - 2), *(vY.end() - 2), vX1.back(), vY.back(), ubound, vY.back()))
//        {
//            vX1[vX1.size() - 1] = ubound;
//            vX2[vX2.size() - 1] = itv + ubound;
//        } else
//        {
//            vX1.push_back(ubound);
//            vX2.push_back(itv + ubound);
//            vY.push_back(vY.back());
//        }
//    }
};

class LPFunction
{
public:

    vector<int> vX, vY; // x, f(x)
    //    unordered_set<pair<int, int>, boost::hash<pair<int, int>>> vSupport; //supportive nodes

    // supVId, the interval # of underlying lpf support -> interval # of this lpf
//    unordered_map<int, unordered_map<int, int>> vSupport;
    vector<map<int, vector<int>>> vSupportPre; // supVId, the interval # of underlying lpf support
    // the number of supportive options in each interval, some options may need multiple supportive nodes
    vector<int> cntOfEachInt;
    int ID1, ID2;  //ID1->ID2
    int minY{}; // min{vY}
    int maxY{}; // max{vY}
    //True: a min cost function, in-label
    //False: a cost function, out-label, need to refine the upper bound
//    int lastItvSubToChange{};
    int lowerBound, upperBound;

    static int computeY(int x1, int x2, int y1, int y2, int x);

    //test if (x2,y2) can be removed

    int setValue(vector<int> &vX2, vector<int> &vY2, int compStartX = -1);

    int setValueNoComp(vector<int> &vX2, vector<int> &vY2);

    int setValue(vector<int> &vX2, vector<int> &vY2,
                 vector<map<int, vector<int>>> &vSup, vector<int> &cntRec, int compStartX = -1);

    int setValue(vector<int> &vX2, vector<int> &vY2, vector<map<int, vector<int>>> &vSup, vector<int> &cntRec,
                 const vector<unordered_map<int, unordered_map<int, CatSupRec>>> &intermediateLPFs,
                 int compStartX = -1);

    LPFunction(int id1, int id2, int lowerBound, int upperBound);

    LPFunction(int id1, int id2, int lowerBound, int upperBound, vector<int> &vX, vector<int> &vY);

    LPFunction(int id1, int id2, int lowerBound, int upperBound,
               vector<int> &vX, vector<int> &vY, vector<map<int, vector<int>>> &vSupport);

    LPFunction();

    void extendFunction(int lBound, int uBound);

    void dummyLastItv(int itvLen);

    bool vSupportLastPartContains(
            int x, const unordered_map<int, unordered_map<int, CatSupRec>> &intermediateLPFs) const;

    bool vSupportContains(int x) const;

    int getY(int x) const;

    vector<int> getVY2(const vector<int> &inX, const vector<int> &inY) const;

    void getXF2NoComp(const LPFunction &f2, vector<int> &vXAll1, vector<int> &vXAll2, vector<int> &vYAll);

    int getX(int x1, int y1, int x2, int y2, int f2x) const;

    bool equal(const LPFunction &f2) const;

    bool equal(LPFunction &f2, vector<int> &changedPos) const;

    int LPFGetX(int y);

    void LPFConcactBeginPos(const LPFunction &f2, int &pos1, int &pos2, int &needInc) const;

    void LPFMinBeginPos(const LPFunction &f2, int &pos1, int &pos2, int itvLen) const;

    LPFunction LPFCatSupportNoComp(const LPFunction &f2, int lBound, int uBound);

    LPFunction LPFCatSupport(LPFunction &f2, int lBound, int uBound);

    LPFunction LPFCatSupport(LPFunction &f2, CatSupRec &catRec, int startT=0);

    LPFunction LPFCatSupportExtend(LPFunction &f2, int itvLen, CatSupRec &catSupRec);

    void scSupToLbSup(int lid);

    void LPFMinSupByPlaneSweepCompPart(
            const LPFunction &lpf, int endX, int &pos1, int &pos2, vector<int> &newVX, vector<int> &newVY,
            vector<map<int, vector<int>>> &supPreInfo, vector<int> &cntRec,
            pair<int, map<int, vector<int>>> &curXY, pair<int, map<int, vector<int>>> &nextXY,
            vector<vector<int>> &orderedLPFs, vector<int> &orderOfLPFs) const;

    LPFunction LPFMinSupByPlaneSweepNoSV(const LPFunction &lpf) const;

    LPFunction LPFMinSupByPlaneSweepNoSV2(const LPFunction &lpf, int itvLen = 1800) const;

    LPFunction LPFMinSupForDec(const LPFunction &lpf) const;

    LPFunction LPFMinSupForExtend(const LPFunction &lpf, int itvLen) const;

    LPFunction LPFMinSupByPlaneSweep(const LPFunction &lpf) const;

    bool dominate(const LPFunction &f2) const;

    int equalValue(const LPFunction &f2) const;

    int equalValue2(const LPFunction &f2) const;

    void display() const; // display the function (vX, vY, minX, minY, maxX, maxY

    void leftTrim(int lowerBound);

    LPFunction LPFTruncate(int startT);
};

#endif
