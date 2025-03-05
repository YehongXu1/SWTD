#include "LPFunction.h"
#include "tools.h"

//int itvLen = 1800;

int LPFunction::setValueNoComp(vector<int> &vX2, vector<int> &vY2)
{
    assert(vX2.size() == vY2.size() and !vX2.empty());

    minY = INF;
    maxY = -1;
    vX.clear();
    vY.clear();
    cntOfEachInt.clear();
    vSupportPre.clear();
    vSupportPre.reserve(vX2.size());

    vX.reserve(vX2.size());
    vY.reserve(vY2.size());
    vector<int>::iterator ivX, ivY;
    for (ivX = vX2.begin(), ivY = vY2.begin(); ivX != vX2.end(); ivX++, ivY++)
    {
        if (ivX != vX2.begin())
        {
            if(*ivX <= *(ivX - 1))
                continue;
        }

        vX.emplace_back(*ivX);
        vY.emplace_back(*ivY);

        if (*ivY < minY)
            minY = *ivY;
        if (*ivY > maxY)
            maxY = *ivY;
    }

    assert(!vX.empty());
    cntOfEachInt.assign(vX.size() - 1, 1);

    if (vX[0] != vX2[0] or vX.back() != vX2.back() or vY[0] != vY2[0] or vY.back() != vY2.back())
    {
        for (int i = 0; i < vY2.size(); i++)
        {
            cout << "(" << vX2[i] << " " << vY2[i] << ") ";
        }
        cout << endl;
        display();
    }
    assert(vX[0] == vX2[0] and vX.back() == vX2.back());
    assert(vY[0] == vY2[0] and vY.back() == vY2.back());

    return vX.size();
}

int LPFunction::setValue(vector<int> &vX2, vector<int> &vY2, int compStartX)
{
    if (vX2.size() > 1)
        assert(vX2.front() < vX2.back());
    minY = INF;
    maxY = -1;
    vector<int>::iterator ivX, ivY;
    vX.clear();
    vY.clear();
    cntOfEachInt.clear();
    vSupportPre.clear();
    vSupportPre.reserve(vX2.size());

    vX.reserve(vX2.size());
    vY.reserve(vY2.size());
    for (ivX = vX2.begin(), ivY = vY2.begin(); ivX != vX2.end(); ivX++, ivY++)
    {
        if (ivX != vX2.begin())
        {
            if (*ivX <= *(ivX - 1))
            {
                if (ivX != vX2.end() - 1)
                    continue;
                else
                {
                    while (vX.size() > 1 and *ivX <= vX.back())
                    {
                        vX.pop_back();
                        vY.pop_back();
                    }
                }
            }
        }

        while (vX.size() > 1 and vX.back() > compStartX and Tools::redundant(
                *(vX.end() - 2), *(vY.end() - 2), *(vX.end() - 1), *(vY.end() - 1), *ivX, *ivY))
        {
            // compression
            vX.pop_back();
            vY.pop_back();
        }

        vX.emplace_back(*ivX);
        vY.emplace_back(*ivY);

        if (*ivY < minY)
            minY = *ivY;
        if (*ivY > maxY)
            maxY = *ivY;
    }

    cntOfEachInt.assign(vX.size() - 1, 1);
    assert(vX[0] == vX2[0] and vX.back() == vX2.back());
    assert(vY[0] == vY2[0] and vY.back() == vY2.back());
    return vX.size();
}


int LPFunction::setValue(
        vector<int> &vX2, vector<int> &vY2, vector<map<int, vector<int>>> &vSup, vector<int> &cntRec, int compStartX)
{
    assert(vX2.size() == vY2.size() && vY2.size() == vSup.size() + 1 && cntRec.size() == vSup.size());
    // this function is called in LPFMinSupForDec
    minY = INF;
    maxY = -1;
    vector<int>::iterator ivX, ivY, cntPtr;
    vector<map<int, vector<int>>>::iterator ivS;
    vX.clear();
    vY.clear();
    vSupportPre.clear();
    cntOfEachInt.clear();

    vX.reserve(vX2.size());
    vY.reserve(vY2.size());
    vSupportPre.reserve(vSup.size());
    cntOfEachInt.reserve(vSup.size());

    for (ivX = vX2.begin(), ivY = vY2.begin(), ivS = vSup.begin(), cntPtr = cntRec.begin();
         ivX != vX2.end(); ivX++, ivY++, ivS++, cntPtr++)
    {
        if (ivX != vX2.begin())
            assert(*ivX > *(ivX - 1));

        while (vX.size() > 1 and vX.back() > compStartX and Tools::redundant(
                *(vX.end() - 2), *(vY.end() - 2), *(vX.end() - 1), *(vY.end() - 1), *ivX, *ivY))
        {
            // compression
            vX.pop_back();
            vY.pop_back();
            assert(vSupportPre.size() > 1);
            // if there exist a vertex that supports both interval,
            // then we only record them and disregard those that only supports one interval
            map<int, vector<int>> &tmp = *(vSupportPre.end() - 2), tmp2 = *(vSupportPre.end() - 1);
            auto hint = tmp.begin();
            // add values in tmp2 to tmp
            for (const auto &p: tmp2)
            {
                if (tmp.find(p.first) != tmp.end())
                {
                    int last = tmp[p.first].back();
                    for (const auto &v: p.second)
                    {
                        if (last != v)
                            tmp[p.first].emplace_back(v);
                    }
                } else
                {
                    tmp.insert(p);
                }
            }

            vSupportPre.pop_back();

            int cnt = min(cntOfEachInt.back(), *(cntOfEachInt.end() - 2));
            cntOfEachInt.pop_back();
            cntOfEachInt[cntOfEachInt.size() - 1] = cnt;
        }

        vX.emplace_back(*ivX);
        vY.emplace_back(*ivY);

        if (ivS != vSup.end())
        {
            vSupportPre.push_back(*ivS);
            cntOfEachInt.push_back(*cntPtr);
        }

        if (*ivY < minY)
            minY = *ivY;
        if (*ivY > maxY)
            maxY = *ivY;
    }

    assert(vX.size() > 1 and vX.size() == vY.size() and
           vY.size() == vSupportPre.size() + 1 and cntOfEachInt.size() == vSupportPre.size());
    assert(vX[0] == vX2[0] and vX.back() == vX2.back());
    assert(vY[0] == vY2[0] and vY.back() == vY2.back());
    return vX.size();
}


int LPFunction::setValue(
        vector<int> &vX2, vector<int> &vY2, vector<map<int, vector<int>>> &vSup, vector<int> &cntRec,
        const vector<unordered_map<int, unordered_map<int, CatSupRec>>> &intermediateLPFs, int compStart)
{
    assert(vX2.size() == vY2.size() && vY2.size() == vSup.size() + 1 && cntRec.size() == vSup.size());
    // this function is called in LPFMinSupForDec
    minY = INF;
    maxY = -1;
    vector<int>::iterator ivX, ivY, cntPtr;
    vector<map<int, vector<int>>>::iterator ivS;
    vX.clear();
    vY.clear();
    vSupportPre.clear();
    cntOfEachInt.clear();

    vX.reserve(vX2.size());
    vY.reserve(vY2.size());
    vSupportPre.reserve(vSup.size());
    cntOfEachInt.reserve(vSup.size());

    for (ivX = vX2.begin(), ivY = vY2.begin(), ivS = vSup.begin(), cntPtr = cntRec.begin();
         ivX != vX2.end(); ivX++, ivY++, ivS++, cntPtr++)
    {
        if (ivX != vX2.begin())
        {
            assert(*ivX > *(ivX - 1));
        }

        map<int, vector<int>> cleanedSups;
        if (ivS != vSup.end())
        {
            map<int, vector<int>> sups = *ivS;
            assert(!sups.empty());
            if (sups.size() == 1)
            {
                cleanedSups = sups;
            } else
            {
                int t1 = *ivX, t2 = *(ivX + 1);
                for (const auto &p: sups)
                {
                    int x;
                    bool catFunc = true;
                    if (p.first != ID1 and p.first != ID2)
                        x = p.first;
                    else
                    {
                        x = ID1;
                        catFunc = false;
                    }

                    if (intermediateLPFs[x].find(ID1) == intermediateLPFs[x].end()
                        or intermediateLPFs[x].find(ID1)->second.find(ID2) ==
                           intermediateLPFs[x].find(ID1)->second.end())
                    {
                        cout << ID1 << " " << x << " " << ID2 << endl;
                    }
                    assert(intermediateLPFs[x].find(ID1) != intermediateLPFs[x].end()
                           and intermediateLPFs[x].find(ID1)->second.find(ID2) !=
                               intermediateLPFs[x].find(ID1)->second.end());

                    vector<int> supLPFvX;
                    if (catFunc)
                        supLPFvX = intermediateLPFs[x].find(ID1)->second.find(ID2)->second.vX1;
                    else
                        supLPFvX = intermediateLPFs[ID1].find(ID1)->second.find(ID2)->second.vX1;
                    assert(supLPFvX.size() > 1);

                    vector<int> cleanedInvs;
                    cleanedInvs.reserve(p.second.size());
                    for (const auto &intId: p.second)
                    {
                        assert(intId < supLPFvX.size() - 1);
                        int supLPFT1 = supLPFvX[intId], supLPFT2 = supLPFvX[intId + 1];
                        if (supLPFT2 < t1 or supLPFT1 >= t2) // < instead of <= because break points in the lpf
                            continue;
                        cleanedInvs.emplace_back(intId);
                    }
                    if (!cleanedInvs.empty())
                        cleanedSups.insert({p.first, cleanedInvs});
                }
            }
            if (cleanedSups.empty())
            {
                cout << ID1 << " " << ID2 << " " << *ivX << endl;
                for (int i = 0; i < vX2.size(); i++)
                {
                    cout << "(" << vX2[i] << " " << vY2[i];
                    if (i < vSup.size())
                    {
                        cout << ", " << cntRec[i] << ", ";
                        for (const auto &p: vSup[i])
                        {
                            cout << p.first << ": ";
                            for (const auto &v: p.second)
                                cout << v << " ";
                            cout << ", ";
                        }
                        cout << ")" << endl;
                    }
                }

                cout << ")" << endl;

                for (const auto &p: sups)
                {
                    CatSupRec supCatRec;
                    if (p.first != ID1 and p.first != ID2)
                        supCatRec = intermediateLPFs[p.first].find(ID1)->second.find(ID2)->second;
                    else
                        supCatRec = intermediateLPFs[ID1].find(ID1)->second.find(ID2)->second;

                    cout << "w:" << p.first << " ";
                    for (int i = 0; i < supCatRec.vX1.size(); i++)
                    {
                        cout << supCatRec.vX1[i] << " ";
                    }
                    cout << endl;
                }
                cout << "~~~~~~~~~~~" << endl;
            }
            assert(!cleanedSups.empty());
        }

        /** compression **/
        if (vX.size() > 1 and vX.back() > compStart and Tools::redundant(
                *(vX.end() - 2), *(vY.end() - 2), *(vX.end() - 1), *(vY.end() - 1), *ivX, *ivY))
        {
            vX.pop_back();
            vY.pop_back();
            assert(vSupportPre.size() > 1);
            map<int, vector<int>> &tmp = *(vSupportPre.end() - 2), tmp2 = *(vSupportPre.end() - 1);
            auto hint = tmp.begin();
            // add values in tmp2 to tmp
            for (const auto &p: tmp2)
            {
                if (tmp.find(p.first) != tmp.end())
                {
                    int last = tmp[p.first].back();
                    for (const auto &v: p.second)
                    {
                        if (last != v)
                            tmp[p.first].emplace_back(v);
                    }
                } else
                {
                    tmp.insert(p);
                }
            }

            vSupportPre.pop_back();

            int cnt = min(cntOfEachInt.back(), *(cntOfEachInt.end() - 2));
            cntOfEachInt.pop_back();
            cntOfEachInt[cntOfEachInt.size() - 1] = cnt;
        }
        vX.emplace_back(*ivX);
        vY.emplace_back(*ivY);

        if (ivS != vSup.end())
        {
            vSupportPre.push_back(cleanedSups);
            cntOfEachInt.push_back(*cntPtr);
        }

        if (*ivY < minY)
            minY = *ivY;
        if (*ivY > maxY)
            maxY = *ivY;
    }

    assert(vX.size() > 1 and vX.size() == vY.size() and vY.size() == vSupportPre.size() + 1
           and cntOfEachInt.size() == vSupportPre.size());
    assert(vX[0] == vX2[0] and vX.back() == vX2.back());
    assert(vY[0] == vY2[0] and vY.back() == vY2.back());
    return vX.size();
}

LPFunction::LPFunction(int id1, int id2, int lowerBound, int upperBound)
{
    ID1 = id1;
    ID2 = id2;
    this->lowerBound = lowerBound;
    this->upperBound = upperBound;
}

//Build from raw data, remove redundant slower
LPFunction::LPFunction(int id1, int id2, int lowerBound, int upperBound, vector<int> &vX, vector<int> &vY)
{
    ID1 = id1;
    ID2 = id2;
    this->lowerBound = lowerBound;
    this->upperBound = upperBound;
    this->setValue(vX, vY);
}

LPFunction::LPFunction(int id1, int id2, int lowerBound, int upperBound,
                       vector<int> &vX, vector<int> &vY, vector<map<int, vector<int>>> &vSupport)
{
    ID1 = id1;
    ID2 = id2;
    this->lowerBound = lowerBound;
    this->upperBound = upperBound;
    vector<int> cntRec(vX.size() - 1, 1);
    this->setValue(vX, vY, vSupport, cntRec);
}

LPFunction::LPFunction()
{
    ID1 = -1;
    ID2 = -1;
    lowerBound = -1;
    upperBound = -1;
    minY = INF;
    maxY = -1;
}

void LPFunction::scSupToLbSup(int lid)
{
    if (vSupportPre.empty())
        return;
    int size = vSupportPre.size();
    vSupportPre.clear();
    vSupportPre.reserve(size);
    for (int i = 0; i < size; i++)
    {
        vSupportPre.push_back({{lid, {i}}});// only one support
    }
}


int LPFunction::getY(int x) const
{
    if (vX.empty())
        throw runtime_error("The function is Empty!");

    if (x < vX.front())
    {
        return vY[0];
    }

    if (x > vX.back())
    {
        return vY.back();
    }

    vector<int>::const_iterator low, up;
    low = lower_bound(vX.begin(), vX.end(), x);
    if (*low == x)
        return vY[low - vX.begin()];

    low -= 1;
    up = low + 1;

    int pos = low - vX.begin();

    if (*low == x)
        return vY[pos];

    int x1 = *low;
    int x2 = *up;
    int y1 = vY[pos];
    int y2 = vY[pos + 1];

    if (x1 == x2)
        return y1;
    else
        return computeY(x1, x2, y1, y2, x);

}

vector<int> LPFunction::getVY2(const vector<int> &inX, const vector<int> &inY) const
{
    int x1, y1, x2, y2, x, y;
    vector<int> vrY;
    vector<int>::const_iterator ivinX, ivinY, ivX, ivY;
    vrY.reserve(inX.size());

    if (inX.empty())
        return vrY;

    ivinX = inX.begin();
    ivinY = inY.begin();
    while (*ivinX < *(vX.begin()) && ivinX < inX.end())
    {
        vrY.push_back(*ivinY);
        ivinX++;
        ivinY++;
    }

    for (ivX = vX.begin(), ivY = vY.begin(); ivX + 1 != vX.end();)
    {
        if (ivinX >= inX.end())
            break;
        if (*ivX <= *ivinX && *(ivX + 1) >= *ivinX)
        {
            x1 = *ivX;
            y1 = *ivY;
            x2 = *(ivX + 1);
            y2 = *(ivY + 1);
            x = *ivinX;
            y = computeY(x1, x2, y1, y2, x);
            vrY.push_back(y);

            ivinX++;
        } else if (*ivX > *ivinX)
        {
            ivinX++;
        } else if (*(ivX + 1) < *ivinX)
        {
            ivX++;
            ivY++;
        }
    }

    return vrY;
}

int LPFunction::computeY(int x1, int x2, int y1, int y2, int x)
{
    assert(x1 <= x && x <= x2);
    if (x == x1)
        return y1;
    if (x == x2)
        return y2;
    if (x1 == x2)
        return y1;
    return Tools::Round(1.0 * (y2 - y1) * (x - x1) / (x2 - x1) + y1);
}

int LPFunction::LPFGetX(int y)
{
    if (vX.empty() or vX.front() + vY.front() > y)
        return -1;
    int pos = -1;
    for (int i = 0; i < vX.size() - 1; i++)
    {
        if (vX[i] + vY[i] <= y and vX[i + 1] + vY[i + 1] >= y)
        {
            pos = i;
            break;
        }
    }
    assert(pos >= 0);
    int x = getX(vX[pos], vY[pos], vX[pos + 1], vY[pos + 1], y);
    return x;
}

int LPFunction::getX(int x1, int y1, int x2, int y2, int f2x) const
{
    assert(f2x >= x1 + y1 and f2x <= x2 + y2);
    int y = Tools::Round(((double) (x2 - x1) * (f2x - y1) + (y2 - y1) * x1) / (double) (x2 - x1 + y2 - y1));
    return y;
}

void LPFunction::LPFConcactBeginPos(const LPFunction &f2, int &pos1, int &pos2, int &needInc) const
{
    pos1 = 0, pos2 = 0, needInc = -1;
    int needInc1 = INF, needInc1y = INF, needInc2 = INF;
    if (pos2 < f2.vX.size())
        needInc2 = f2.vX[pos2];

    if (pos1 < vX.size())
    {
        needInc1 = vX[pos1];
        needInc1y = vY[pos1];

        // ensure needInc2 is included during concatenation
        // this important when needInc2 < needInc1
        while (pos1 > 0 and vX[pos1] + vY[pos1] > needInc2)
            pos1--;

        // ensure needInc1 is included during concatenation
        // this important when needInc1 < needInc2
        while (pos2 > 0 and f2.vX[pos2 - 1] >= needInc1 + needInc1y)
            pos2--;
    } else
    {
        do
        {
            pos1--;
        } while (pos1 > 0 and vX[pos1] + vY[pos1] > needInc2);
    }

    assert(pos1 >= 0 and pos1 < vX.size());

    if (needInc1 + needInc1y <= needInc2)
    {
        // needInc1 is the real bound
        assert(needInc1 == vX[pos1]);
        if (pos2 < f2.vX.size())
            assert(vX[pos1] + vY[pos1] <= f2.vX[pos2]);

        needInc = needInc1 + needInc1y;
    } else
    {
        // f2's changed part is front of f1's changed part
        assert(pos2 < f2.vX.size() and needInc2 == f2.vX[pos2]);
        if (vX[pos1] + vY[pos1] <= needInc2)
            needInc = needInc2;
        else
            needInc = -1; // pos1 == 0 or need2 <= vX[0] + vY[0]
    }

    assert(pos2 >= 0 and pos2 < f2.vX.size());
}

void LPFunction::LPFMinBeginPos(const LPFunction &f2, int &pos1, int &pos2, int itvLen) const
{
    int start1 = vX.size() - 1, start2 = f2.vX.size() - 1;
    while (start1 >= 0 and vX[start1] > upperBound - itvLen)
        start1--;
    while (start2 >= 0 and f2.vX[start2] > upperBound - itvLen)
        start2--;
    pos1 = min(start1, 0), pos2 = min(start2, (int) f2.vX.size() - 1);


    int needInc1 = vX[pos1], needInc2 = f2.vX[pos2];

    if (needInc2 < needInc1)
    {
        while (pos1 > 0 and vX[pos1] > needInc2)
            pos1 -= 1;
        assert(pos1 >= 0 and pos1 < vX.size());
    }
    pos2 = 0;
    while (pos2 < f2.vX.size() - 1 and f2.vX[pos2] < vX[pos1])
        pos2++;
    if (f2.vX[pos2] > vX[pos1])
        pos2--;
    assert(pos2 >= 0 and pos2 < f2.vX.size() and f2.vX[pos2] <= vX[pos1]);
}

void LPFunction::getXF2NoComp(const LPFunction &f2, vector<int> &vXAll1, vector<int> &vXAll2, vector<int> &vYAll)
{
    if (ID2 != f2.ID1 or lowerBound != f2.lowerBound or upperBound != f2.upperBound)
    {
        display();
        f2.display();
        assert(false);
    }
    int x1, y1, x2, y2, x;
    vector<int>::const_iterator ivX, ivY, ivX2, ivY2;

    vector<int> vXAll, vY1, vY2;
    vXAll.reserve(vX.size() + f2.vX.size());
    vY1.reserve(vX.size() + f2.vX.size());
    vY2.reserve(vX.size() + f2.vX.size());


    int pos1, pos2, needInc;
    LPFConcactBeginPos(f2, pos1, pos2, needInc);

    for (ivX = vX.begin() + pos1, ivY = vY.begin() + pos1, ivX2 = f2.vX.begin() + pos2, ivY2 = f2.vY.begin() + pos2;
         ivX + 1 != vX.end() && ivX2 != f2.vX.end();)
    {
        // might be a bug here if the function is not FIFO
        if (ivX2 != f2.vX.begin())
        {
            if (*(ivX2 - 1) >= *ivX2)
            {
                display();
                f2.display();
                assert(false);
            }
        }

        x1 = *ivX;
        y1 = *ivY;
        x2 = *(ivX + 1);
        y2 = *(ivY + 1);
        if (*ivX + *ivY <= *ivX2 and *ivX2 <= *(ivX + 1) + *(ivY + 1))
        {
            x = getX(x1, y1, x2, y2, *ivX2);

            if (x1 + y1 >= needInc and x1 + y1 < *ivX2 and ivX2 > f2.vX.begin() and x1 + y1 > *(ivX2 - 1))
            {
                vXAll.push_back(x1);
                vY1.push_back(y1);
                int y = computeY(*(ivX2 - 1), *ivX2, *(ivY2 - 1), *ivY2, x1 + y1);
                vY2.push_back(y);
            }

            if (*ivX2 >= needInc)
            {
                vXAll.push_back(x);
                vY1.push_back(*ivX2 - x);
                vY2.push_back(*ivY2);
            }

            ivX2++;
            ivY2++;
        } else if (x1 + y1 > *ivX2)
        {
            ivX2++;
            ivY2++;
        } else if (x2 + y2 < *ivX2)
        {
            if (ivX2 > f2.vX.begin() and x1 + y1 > *(ivX2 - 1))
            {
                vXAll.push_back(x1);
                vY1.push_back(y1);
                int y = computeY(*(ivX2 - 1), *ivX2, *(ivY2 - 1), *ivY2, x1 + y1);
                vY2.push_back(y);
            }

            if (ivX == vX.end() - 2 and ivX2 > f2.vX.begin() and x2 + y2 > *(ivX2 - 1))
            {
                vXAll.push_back(x2);
                vY1.push_back(y2);
                int y = computeY(*(ivX2 - 1), *ivX2, *(ivY2 - 1), *ivY2, x2 + y2);
                vY2.push_back(y);
            }

            ivX++;
            ivY++;
        }
    }

    vXAll1.clear();
    vXAll2.clear();
    vYAll.clear();

    vXAll1.reserve(vXAll.size());
    vXAll2.reserve(vXAll.size());
    vYAll.reserve(vXAll.size());
    for (int i = 0; i < vXAll.size(); i++)
    {
        if (vXAll1.empty() or vXAll[i] > vXAll1.back())
        {
            vXAll1.push_back(vXAll[i]);
            vXAll2.push_back(vXAll[i] + vY1[i]);
            vYAll.push_back(vY1[i] + vY2[i]);
        } else if (vXAll[i] == vXAll1.back())
        {
            int y11 = vYAll.back();
            int y22 = vY1[i] + vY2[i];
            if (y22 < y11)
            {
                vXAll2[vXAll2.size() - 1] = vXAll[i] + vY1[i];
                vYAll[vYAll.size() - 1] = y22;
            }
        } else if (i == vXAll.size() - 1)
        {
            vXAll1.pop_back();
            vXAll2.pop_back();
            vYAll.pop_back();
            vXAll1.push_back(vXAll[i]);
            vXAll2.push_back(vXAll[i] + vY1[i]);
            vYAll.push_back(vY1[i] + vY2[i]);
        }
    }

    int back = vXAll1.empty() ? -1 : vXAll1.back();
    for (int i = 0; i < vX.size(); i++)
    {
        if (vX[i] <= back)
            continue;
        vXAll1.push_back(vX[i]);
        vXAll2.push_back(vY[i] + vX[i]);
        vYAll.push_back(vY[i] + f2.vY.back());
    }

    assert(!vXAll1.empty() and vXAll1.size() == vXAll2.size() and vXAll2.size() == vYAll.size());
}

LPFunction LPFunction::LPFCatSupportNoComp(const LPFunction &f2, int lBound, int uBound)
{
    if (ID2 != f2.ID1)
    {
        display();
        f2.display();
        assert(false);
    }
//    assert(lBound <= vX.front() and uBound >= vX.back());
//    assert(lBound <= f2.vX.front() and uBound >= f2.vX.back());

    // suppose concatenate u -> w and w -> v

    if (vX.empty() || f2.vX.empty())
    {
        LPFunction fr(ID1, f2.ID2, lBound, uBound);
        return fr;
    }


    // get the timestamps t such that you can arrive w at a timestamp in f2.vX by departing from u at t
    vector<int> vXAll1, vXAll2, vYAll;
    getXF2NoComp(f2, vXAll1, vXAll2, vYAll);
    if (vXAll1.empty())
    {
        LPFunction lpf;
        return lpf;
    }

    LPFunction fr(ID1, f2.ID2, lBound, uBound);
    fr.setValueNoComp(vXAll1, vYAll);
    fr.vSupportPre.reserve(fr.vX.size() - 1);
    for (int i = 0; i < fr.vX.size() - 1; i++)
        fr.vSupportPre.push_back({{ID2, {i}}});

    int debugID1 = 1, debugID2 = 1, debugID3 = 1;
    if (ID1 == debugID1 && ID2 == debugID2 && f2.ID2 == debugID3)
    {
        cout << "\nLPFCatSupport" << endl;
        display();
        f2.display();
        fr.display();
    }

    return fr;
}

LPFunction LPFunction::LPFCatSupport(LPFunction &f2, int lBound, int uBound)
{
    // used in query processing, decrease updates, and ground truth
    int debugID1 = 1, debugID2 = 1, debugID3 = 1;
    if (vX.empty() || f2.vX.empty())
    {
        LPFunction fr(ID1, f2.ID2, lBound, uBound);
        return fr;
    }

    if (ID1 == debugID1 and ID2 == debugID2 and f2.ID2 == debugID3)
    {
        cout << "\nLPFCatSupportForDec" << endl;
        display();
        f2.display();
    }
    // for internal decrease updates
    LPFunction lpf = LPFCatSupportNoComp(f2, lBound, uBound);
    assert(lpf.vX.size() > 1);
    LPFunction lpf2(ID1, f2.ID2, lBound, uBound);
    lpf2.setValue(lpf.vX, lpf.vY, lpf.vSupportPre, lpf.cntOfEachInt);
    if (ID1 == debugID1 and ID2 == debugID2 and f2.ID2 == debugID3)
        lpf.display();

    return lpf2;
}

void LPFunction::dummyLastItv(int itvLen)
{
    for (int i = (int) vX.size() - 1; i >= 0; i--)
    {
        if (vX[i] <= upperBound - itvLen)
        {
            break;
        }
        vY[i] = INF;
        maxY = INF;
    }
}

void LPFunction::extendFunction(int lBound, int uBound)
{
    assert(vX.size() > 1);
    lowerBound = lBound;
    upperBound = uBound;
    vX.emplace_back(upperBound);
    vY.emplace_back(vY.back());
    vSupportPre.emplace_back(vSupportPre.back());
    cntOfEachInt.emplace_back(cntOfEachInt.back());
}


LPFunction LPFunction::LPFCatSupportExtend(LPFunction &f2, int itvLen, CatSupRec &pp)
{
    // no compression in existing part
    if (vX.empty() || f2.vX.empty())
    {
        LPFunction fr(ID1, f2.ID2, lowerBound, upperBound);
        return fr;
    }

    int debugID1 = 1, debugID2 = 1, debugID3 = 1, ubound = 1;
    if (ID1 == debugID1 and ID2 == debugID2 and f2.ID2 == debugID3 and upperBound == ubound)
    {
        cout << "\nLPFCatSupportExtend" << endl;
        display();
        f2.display();
    }
    // suppose concatenate u -> w and w -> v
    vector<int>::iterator ivX1, ivX2, ivY;
    map<int, int> mX;
    map<int, int>::iterator imX;

    vector<int> vXAll1, vXAll2, vYAll;

    getXF2NoComp(f2, vXAll1, vXAll2, vYAll);

//    if (pp.vX1.empty())
//        cout << 1;
    if (vXAll1.empty())
    {
        display();
        f2.display();
        assert(!vXAll1.empty());
    }

    vector<int> vrX, vrY, vf2X;
    vrX.reserve(vXAll1.size() + pp.vX1.size() + 1);
    vf2X.reserve(vXAll1.size() + pp.vX1.size() + 1);
    vrY.reserve(vXAll1.size() + pp.vX1.size() + 1);

    int mark = 0, minY2 = INF, maxY2 = -1;
    if (pp.vX1.empty())
    {
        pp.vX1.push_back(vX[0]);
        pp.vX2.push_back(vX[0] + vY[0]);
        pp.vY.push_back(vYAll[0]);
    }

    while (mark < pp.vX1.size() and pp.vX1[mark] < vXAll1.front())
    {
        vrX.emplace_back(pp.vX1[mark]);
        vf2X.emplace_back(pp.vX2[mark]);
        vrY.emplace_back(pp.vY[mark]);
        if (pp.vY[mark] < minY2)
            minY2 = pp.vY[mark];
        if (pp.vY[mark] > maxY2)
            maxY2 = pp.vY[mark];
        mark++;
    }

    LPFunction fr(ID1, f2.ID2, lowerBound, upperBound);
    for (ivX1 = vXAll1.begin(), ivX2 = vXAll2.begin(), ivY = vYAll.begin();
         ivX1 != vXAll1.end() && ivX2 != vXAll2.end() && ivY != vYAll.end(); ivX1++, ivX2++, ivY++)
    {
        if (vrX.size() > 1 and vrX.back() != upperBound - itvLen and Tools::redundant(
                *(vrX.end() - 2), *(vrY.end() - 2), vrX.back(), vrY.back(), *ivX1, *ivY))
        {
            vrX.pop_back();
            vf2X.pop_back();
            vrY.pop_back();
        }

        if (*ivX1 > upperBound - itvLen)
        {
            assert(!vrX.empty());
            if (vrX.back() < upperBound - itvLen)
            {
                vf2X.push_back(vf2X.back() - vrX.back() + upperBound - itvLen);
                vrX.push_back(upperBound - itvLen);
                vrY.push_back(vrY.back());
            }


        }

        vrX.push_back(*ivX1);
        vf2X.push_back(*ivX2);
        vrY.push_back(*ivY);
        if (*ivY < minY2)
            minY2 = *ivY;
        if (*ivY > maxY2)
            maxY2 = *ivY;
    }

    pp.vX1 = vrX;
    pp.vX2 = vf2X;
    pp.vY = vrY;
    fr.vX = vrX;
    fr.vY = vrY;
    fr.minY = minY2;
    fr.maxY = maxY2;
    if (pp.vX1.empty() or pp.vX1.front() != 0)
        assert(false);
    fr.vSupportPre.reserve(fr.vX.size() - 1);
    for (int k = 0; k < fr.vX.size() - 1; k++)
        fr.vSupportPre.push_back({{ID2, {k}}});

    if (ID1 == debugID1 && ID2 == debugID2 && f2.ID2 == debugID3)
        fr.display();
    return fr;
}


LPFunction LPFunction::LPFCatSupport(LPFunction &f2, CatSupRec &pp, int startT)
{
    if (vX.size() < 2 || f2.vX.size() < 2 || startT >= vX.back())
    {
        LPFunction fr(ID1, f2.ID2, lowerBound, upperBound);
        return fr;
    }

    int debugID1 = 1, debugID2 = 1, debugID3 = 1;
    if (ID1 == debugID1 and f2.ID2 == debugID3)
    {
        cout << "\nLPFCatSupport" << endl;
        display();
        f2.display();
    }

    int upperBound2 = upperBound;

    vector<int> vXX1 = vX, vYX1 = vY, vXX2 = f2.vX, vYX2 = f2.vY;

    if (startT > lowerBound)
    {
        upperBound2 += 600;
        if (vX.back() == upperBound)
        {
            vXX1.push_back(upperBound2);
            vYX1.push_back(vY.back());
        }
        if (f2.vX.back() == f2.upperBound)
        {
            vXX2.push_back(upperBound2);
            vYX2.push_back(f2.vY.back());
        }

        auto it = std::lower_bound(vXX1.begin(), vXX1.end(), startT);

        // Create new vectors starting from t
        std::vector<int> newVX;
        std::vector<int> newVY;

        if (it != vX.end() && *it == startT)
        {
            // If t is in vX, start directly from t (no need to add (t, y_t))
            newVX.insert(newVX.end(), it, vXX1.end());
            newVY.insert(newVY.end(), vYX1.begin() + (it - vXX1.begin()), vYX1.end());
        } else
        {
            // If t is not in vX, add the interpolated point (t, y_t)
            newVX.push_back(startT);
            newVY.push_back(getY(startT));

            // Add the remaining points from vX and vY
            newVX.insert(newVX.end(), it, vXX1.end());
            newVY.insert(newVY.end(), vYX1.begin() + (it - vXX1.begin()), vYX1.end());
        }

        vXX1 = newVX;
        vYX1 = newVY;
    }
    // suppose concatenate u -> w and w -> v

    // get the timestamps t such that you can arrive w at a timestamp in f2.vX by departing from u at t
    vector<int> vXAll, vXAll2, vYAll;

    int x1, y1, x2, y2, x;
    vector<int>::const_iterator ivX = vXX1.begin(), ivY, ivX2, ivY2;


    for (ivY = vYX1.begin(), ivX2 = vXX2.begin(), ivY2 = vYX2.begin();
         ivX + 1 != vXX1.end() && ivX2 != vXX2.end();)
    {
        // might be a bug here if the function is not FIFO
        if (ivX2 != vXX2.begin())
        {
            if (*(ivX2 - 1) >= *ivX2)
            {
                continue;
            }
        }

        x1 = *ivX;
        y1 = *ivY;
        x2 = *(ivX + 1);
        y2 = *(ivY + 1);
        if (x1 + y1 <= *ivX2 and *ivX2 <= *(ivX + 1) + *(ivY + 1))
        {
            x = getX(x1, y1, x2, y2, *ivX2);

            if (x1 + y1 < *ivX2 and ivX2 > vXX2.begin() and x1 + y1 > *(ivX2 - 1))
            {
                int y = computeY(*(ivX2 - 1), *ivX2, *(ivY2 - 1), *ivY2, x1 + y1);
                if (vXAll.size() > 1 and Tools::redundant(
                        *(vXAll.end() - 2), *(vYAll.end() - 2), vXAll.back(), vYAll.back(), x1, y + y1))
                {
                    vXAll.pop_back();
                    vXAll2.pop_back();
                    vYAll.pop_back();
                }
                vXAll.push_back(x1);
                vXAll2.push_back(x1 + y1);
                vYAll.push_back(y + y1);

                if (x1 + y + y1 > upperBound2)
                    break;
            }

            if (vXAll.size() > 1 and
                Tools::redundant(*(vXAll.end() - 2), *(vYAll.end() - 2), vXAll.back(), vYAll.back(),
                                 x, *ivY2 + *ivX2 - x))
            {
                vXAll.pop_back();
                vXAll2.pop_back();
                vYAll.pop_back();
            }

            vXAll.push_back(x);
            vXAll2.push_back(*ivX2);
            vYAll.push_back(*ivY2 + *ivX2 - x);

            if (*ivY2 + *ivX2 > upperBound2)
                break;

            ivX2++;
            ivY2++;
        } else if (x1 + y1 > *ivX2)
        {
            ivX2++;
            ivY2++;
        } else if (x2 + y2 < *ivX2)
        {
            if (ivX2 > vXX2.begin() and x1 + y1 > *(ivX2 - 1))
            {
                int y = computeY(*(ivX2 - 1), *ivX2, *(ivY2 - 1), *ivY2, x1 + y1);
                if (vXAll.size() > 1 and Tools::redundant(
                        *(vXAll.end() - 2), *(vYAll.end() - 2), vXAll.back(), vYAll.back(), x1, y + y1))
                {
                    vXAll.pop_back();
                    vXAll2.pop_back();
                    vYAll.pop_back();
                }

                vXAll.push_back(x1);
                vXAll2.push_back(x1 + y1);
                vYAll.push_back(y + y1);

                if (x1 + y + y1 > upperBound2)
                    break;
            }

            if (ivX == vXX1.end() - 2 and ivX2 > vXX2.begin() and x2 + y2 > *(ivX2 - 1))
            {
                int y = computeY(*(ivX2 - 1), *ivX2, *(ivY2 - 1), *ivY2, x2 + y2);
                if (vXAll.size() > 1 and Tools::redundant(
                        *(vXAll.end() - 2), *(vYAll.end() - 2), vXAll.back(), vYAll.back(), x2, y + y2))
                {
                    vXAll.pop_back();
                    vXAll2.pop_back();
                    vYAll.pop_back();
                }
                vXAll.push_back(x2);
                vXAll2.push_back(x2 + y2);

                vYAll.push_back(y + y2);

                if (x2 + y + y2 > upperBound2)
                    break;
            }

            ivX++;
            ivY++;
        }
    }

    if (!vXAll.empty() and vXAll[0] != startT)
    {
        display();
        f2.display();
        cout << "Start: " << startT << endl;
        for (int i = 0; i < vXAll.size(); i++)
            cout << vXAll[i] << " ";
        cout << endl;
        for (int i = 0; i < vXX1.size(); i++)
            cout << vXX1[i] << " ";
        assert(false);
    }

    if (vXAll.size() < 2)
    {
        LPFunction fr(ID1, f2.ID2, lowerBound, upperBound);
        return fr;
    }

    if (vXAll.back() + vYAll.back() > upperBound2)
    {
        double tangent = 1.0 * (vYAll.back() - vYAll[vYAll.size() - 2]) / (vXAll.back() - vXAll[vXAll.size() - 2]);
        vXAll.erase(vXAll.end() - 1);
        vXAll2.erase(vXAll2.end() - 1);
        vYAll.erase(vYAll.end() - 1);

        int t = ceil(upperBound2 - vYAll.back() + tangent * vXAll.back()) / (1 + tangent);
        if (t + 1 < vXAll.back())
        {
            display();
            f2.display();
            for (int i = 0; i < vXAll.size(); i++)
                cout << vXAll[i] << " ";
            cout << t << endl;
            for (int i = 0; i < vXAll2.size(); i++)
                cout << vXAll2[i] << " ";
            cout << t + getY(t) << endl;
            for (int i = 0; i < vYAll.size(); i++)
                cout << vYAll[i] << " ";
            cout << upperBound2 - t << endl;
            cout << endl;
            assert(false);
        }

        if (vXAll.size() > 1 and Tools::redundant(
                *(vXAll.end() - 2), *(vYAll.end() - 2), vXAll.back(), vYAll.back(), t, upperBound2 - t))
        {
            vXAll.pop_back();
            vXAll2.pop_back();
            vYAll.pop_back();
        }

        vXAll.push_back(t);
        vXAll2.push_back(t + getY(t));
        vYAll.push_back(upperBound2 - t);
    }

    if (upperBound2 > upperBound and vXAll.back() > upperBound)
    {
        // truncate at upperBound

        while (vXAll.size() > 1 and vXAll[vXAll.size() - 2] > upperBound)
        {
            vXAll.pop_back();
            vXAll2.pop_back();
            vYAll.pop_back();
        }

        if (vXAll.back() > upperBound)
        {
            int lastX = vXAll.back();
            int lastY = vYAll.back();
            vXAll.pop_back();
            vXAll2.pop_back();
            vYAll.pop_back();

            if (vXAll.back() < upperBound)
            {
                vXAll.push_back(upperBound);
                vXAll2.push_back(upperBound + vY.back());
                vYAll.push_back(vY.back() + f2.vY.back());
            }
        }

    }

    if (vXAll[vXAll.size() - 2] == vXAll.back())
        cout << 1;
    LPFunction fr(ID1, f2.ID2, lowerBound, upperBound);
    fr.setValueNoComp(vXAll, vYAll);
    fr.vSupportPre.reserve(fr.vX.size() - 1);
    for (int i = 0; i < fr.vX.size() - 1; i++)
        fr.vSupportPre.push_back({{ID2, {i}}});

    if (ID1 == debugID1 && f2.ID2 == debugID3)
    {
        cout << lowerBound << " " << upperBound << endl;
        fr.display();
    }

    pp.vX1 = vXAll;
    pp.vX2 = vXAll2;
    pp.vY = vYAll;

    return fr;
}

LPFunction LPFunction::LPFTruncate(int startT)
{
    if (vX.empty() or startT >= vX.back())
    {
        LPFunction fr(ID1, ID2, lowerBound, upperBound);
        return fr;
    }

    if (startT == vX[0])
        return *this;

    // Find the first index where x >= startT
    auto it = std::lower_bound(vX.begin(), vX.end(), startT);
    // Create new vectors for the truncated function
    vector<int> newVX;
    vector<int> newVY;
    vector<std::map<int, std::vector<int>>> newVSupportPre;
    vector<int> newCntOfEachInt;

    size_t startIndex = std::distance(vX.begin(), it);
    // If startT is not in vX, insert it and interpolate f(startT)
    if (it == vX.end() || *it != startT)
    {
        int yStart = getY(startT);
        newVX.push_back(startT);
        newVY.push_back(yStart);

        // Add the segment starting from the next index
        newVX.insert(newVX.end(), vX.begin() + startIndex, vX.end());
        newVY.insert(newVY.end(), vY.begin() + startIndex, vY.end());

        // Add the corresponding vSupportPre and cntOfEachInt
        newVSupportPre.insert(newVSupportPre.end(), vSupportPre.begin() + startIndex - 1, vSupportPre.end());
        newCntOfEachInt.insert(newCntOfEachInt.end(), cntOfEachInt.begin() + startIndex - 1, cntOfEachInt.end());
    } else
    {
        // If startT is in vX, start directly from startIndex
        newVX.insert(newVX.end(), vX.begin() + startIndex, vX.end());
        newVY.insert(newVY.end(), vY.begin() + startIndex, vY.end());

        // Add the corresponding vSupportPre and cntOfEachInt
        newVSupportPre.insert(newVSupportPre.end(), vSupportPre.begin() + startIndex, vSupportPre.end());
        newCntOfEachInt.insert(newCntOfEachInt.end(), cntOfEachInt.begin() + startIndex, cntOfEachInt.end());
    }

    LPFunction fr = LPFunction(ID1, ID2, lowerBound, upperBound);
    fr.vX = newVX;
    fr.vY = newVY;
    fr.vSupportPre = newVSupportPre;
    fr.cntOfEachInt = newCntOfEachInt;
    fr.minY = INF;
    fr.maxY = -1;
    for (int y: fr.vY)
    {
        if (y < fr.minY)
            fr.minY = y;
        if (y > fr.maxY)
            fr.maxY = y;
    }
    return fr;
}


//LPFunction LPFunction::LPFCatSupport(LPFunction &f2, int lBound, int uBound, CatSupRec &pp)
//{
//    if (vX.empty() || f2.vX.empty())
//    {
//        LPFunction fr(ID1, f2.ID2, lBound, uBound);
//        return fr;
//    }
//
//    int debugID1 = 1, debugID2 = 1, debugID3 = 1;
//    if (ID1 == debugID1 and ID2 == debugID2 and f2.ID2 == debugID3)
//    {
//        cout << "\nLPFCatSupport" << endl;
//        display();
//        f2.display();
//    }
//    // suppose concatenate u -> w and w -> v
//    vector<int>::iterator ivX1, ivX2, ivY;
//    map<int, int> mX;
//    map<int, int>::iterator imX;
//
//    // get the timestamps t such that you can arrive w at a timestamp in f2.vX by departing from u at t
//    vector<int> vXAll1, vXAll2, vYAll;
//    getXF2NoComp(f2, vXAll1, vXAll2, vYAll);
//
//    if (vXAll1.empty())
//    {
//        LPFunction fr(ID1, f2.ID2, lBound, uBound);
//        if (pp.vX1.empty())
//            return fr;
//        fr.setValueNoComp(pp.vX1, pp.vY);
//        fr.vSupportPre.reserve(fr.vX.size() - 1);
//        fr.lastItvSubToChange = fr.vX.size();
//        for (int i = 0; i < fr.vX.size() - 1; i++)
//            fr.vSupportPre.push_back({{ID2, {i}}});
//        return fr;
//    }
//
//    vector<int> vrX, vrY, vf2X;
//    vrX.reserve(vX.size() + f2.vX.size());
//    vrY.reserve(vX.size() + f2.vX.size());
//    vf2X.reserve(vX.size() + f2.vX.size());
//
//    for (ivX1 = vXAll1.begin(), ivX2 = vXAll2.begin(), ivY = vYAll.begin();
//         ivX1 != vXAll1.end() && ivX2 != vXAll2.end() && ivY != vYAll.end(); ivX1++, ivX2++, ivY++)
//    {
//        while (vrX.size() > 1 and Tools::redundant(
//                *(vrX.end() - 2), *(vrY.end() - 2), vrX.back(), vrY.back(), *ivX1, *ivY))
//        {
//            vrX.pop_back();
//            vf2X.pop_back();
//            vrY.pop_back();
//        }
//
//        vrX.push_back(*ivX1);
//        vf2X.push_back(*ivX2);
//        vrY.push_back(*ivY);
//    }
//
//    LPFunction fr(ID1, f2.ID2, lBound, uBound);
//    if (vrX.empty())
//    {
//        display();
//        f2.display();
//        for (int i = 0; i < vXAll1.size(); i++)
//            cout << "(" << vXAll1[i] << " " << vXAll2[i] << " " << vYAll[i] << ")" << endl;
//        assert(false);
//    }
//    fr.setValueNoComp(vrX, vrY);
//    fr.vSupportPre.reserve(fr.vX.size() - 1);
//    for (int i = 0; i < fr.vX.size() - 1; i++)
//        fr.vSupportPre.push_back({{ID2, {i}}});
//
//    if (ID1 == debugID1 && ID2 == debugID2 && f2.ID2 == debugID3)
//    {
//        cout << lBound << " " << uBound << endl;
//        fr.display();
//    }
//
//    pp.vX1 = vrX;
//    pp.vX2 = vf2X;
//    pp.vY = vrY;
//
//    return fr;
//}


void LPFunction::display() const
{
    cout << "From " << ID1 << " to " << ID2 << "|" << endl << " End Points:";
    for (int i = 0; i < vX.size(); i++)
    {
        cout << "(" << vX[i] << "," << vY[i] << ")\t";
    }
    cout << endl;
    cout << " Support V: ";
    for (const auto &interval: vSupportPre)
    {
        cout << "(";
        for (const auto &p: interval)
        {
            cout << p.first << ":";
            for (const auto &v: p.second)
                cout << v << ",";
            cout << "|";
        }
        cout << ")\t";
    }
    cout << "\n Cnt: ";
    for (const auto &cnt: cntOfEachInt)
    {
        cout << cnt << "\t";
    }
    cout << "\n Min cost:" << minY << " Max cost:" << maxY;
    cout << "\n Lowerbound:" << lowerBound << " Upperbound:" << upperBound << endl;
    if (vX.size() > 1)
    {
//        assert(vX.size() == vY.size() and vY.size() == vSupportPre.size() + 1);
    }
}

//0: no intersection
//1: p11 and xy safe
//2: p21 and xy safe
//3: on the same line
bool LPFunction::equal(LPFunction &f2, vector<int> &changedPos) const
{
    int acc = 5;
    assert(f2.vX.size() > 1);
    assert(ID1 >= 0 and ID2 >= 0 and ID1 == f2.ID1 and ID2 == f2.ID2);
    changedPos.clear();
    changedPos.reserve(vX.size());
    int pos = 0, firstDiffPos = min(f2.vSupportPre.size(), vSupportPre.size());
    while (pos < vX.size() and pos < f2.vX.size())
    {
        int x1 = vX[pos], y1 = vY[pos];
        int x2 = f2.vX[pos], y2 = f2.vY[pos];
        bool diffX = (x1 - x2 > acc || x2 - x1 > acc);
        bool diffY = (y1 - y2 > acc || y2 - y1 > acc);

        if (diffX or diffY)
        {
            if (pos > 0 and (changedPos.empty() or changedPos.back() != pos - 1))
                changedPos.push_back(pos - 1);

            if (pos <= firstDiffPos)
            {
                firstDiffPos = pos > 0 ? pos - 1 : 0;
            }
            changedPos.push_back(pos);
        }

        pos += 1;
    }


    while (pos < vX.size())
    {
        if (pos > 0 and (changedPos.empty() or changedPos.back() != pos - 1))
            changedPos.push_back(pos - 1);

        changedPos.push_back(pos);

        pos += 1;
    }

    while (!changedPos.empty() and changedPos.back() >= vX.size() - 1)
        changedPos.pop_back();

    if (changedPos.empty() and vX.size() != f2.vX.size())
    {
        // in this case, there is no change in the function itself, but new turning points are added to the function
        // so we still think that the two functions are different
        // we mark this phenomenon by
        changedPos.push_back(-1);
    }

    return changedPos.empty() and vX.size() == f2.vX.size();
}

bool LPFunction::equal(const LPFunction &f2) const
{
    if (vX.size() < 2 or f2.vX.size() < 2)
        return false;

    if (f2.ID1 < 0 or f2.ID2 < 0)
        return false;

    assert(ID1 == f2.ID1 and ID2 == f2.ID2);

    if (minY != f2.minY)
        return false;

    if (vX.size() != f2.vX.size())
        return false;

    int acceptError = 10;
    vector<int>::const_iterator ivX1, ivX2, ivY1, ivY2;
    for (ivX1 = vX.begin(), ivX2 = f2.vX.begin(), ivY1 = vY.begin(), ivY2 = f2.vY.begin();
         ivX1 != vX.end(); ivX1++, ivX2++, ivY1++, ivY2++)
    {
        if (*ivX1 - *ivX2 > acceptError || *ivX2 - *ivX1 > acceptError)
            return false;
        else if (*ivY1 - *ivY2 > acceptError || *ivY2 - *ivY1 > acceptError)
            return false;
    }

    return true;
}

int LPFunction::equalValue(const LPFunction &f2) const
{
    // return the maximum difference between the two functions
    if (ID1 == -1 || ID2 == -1)
        return false;

    if (vX.size() < 2)
        return false;
    else if (f2.vX.size() < 2)
        return true;

//    assert(ID1 == f2.ID1 and ID2 == f2.ID2);

    // return true if this.vY <= f2.vY for all vX
    vector<int>::const_iterator ivY, ivYtmp;
    vector<int> vYtmp;
    vector<int>::const_iterator ivX, ivX2, ivY2;

    vYtmp = f2.getVY2(vX, vY); // getVY2 is used to get the Y values for the given X values of the function f2

    ivX = vX.begin();
    ivX2 = f2.vX.begin();
    ivY = vY.begin();
    ivY2 = f2.vY.begin();
    ivYtmp = vYtmp.begin();

    int maxDiff = 0;
    for (; ivX != vX.end() - 1 && ivY != vY.end() - 1 && ivYtmp != vYtmp.end() - 1; ivY++, ivYtmp++, ivX++)
    {
        while (*ivX2 < *ivX)
        {
            ivX2++;
            ivY2++;
        }

        int diff;
        if (*ivY > *ivYtmp)
        {
            diff = *ivY - *ivYtmp;
            if (ivX > vX.begin() and *ivX - 1 == *(ivX - 1))
            {
                // might be a breaking point
                diff = -1;
            } else if (diff > maxDiff and ivY != vY.begin())
            {
                diff = min(diff, abs(*ivY2 - *(ivY - 1)));
            }
        } else
        {
            diff = *ivYtmp - *ivY;
            if (ivX < vX.end() - 1 and *ivX + 1 == *(ivX + 1))
            {
                // might be a breaking point
                diff = -1;
            } else if (diff > maxDiff and ivY != vY.end() - 1)
            {
                diff = min(diff, abs(*ivY2 - *(ivY + 1)));
            }
        }

        if (diff > maxDiff)
        {
            maxDiff = diff;
        }
    }

    vYtmp.clear();
    vYtmp = getVY2(f2.vX, f2.vY);

    ivY = f2.vY.begin();
    ivY2 = vY.begin();
    ivYtmp = vYtmp.begin();
    ivX = f2.vX.begin();
    ivX2 = vX.begin();

    for (; ivX != f2.vX.end() - 1 && ivY != f2.vY.end() - 1 && ivYtmp != vYtmp.end() - 1; ivY++, ivYtmp++, ivX++)
    {
        while (*ivX2 < *ivX)
        {
            ivX2++;
            ivY2++;
        }

        int diff;
        if (*ivY > *ivYtmp)
        {
            diff = *ivY - *ivYtmp;
            if (ivX > f2.vX.begin() and *ivX - 1 == *(ivX - 1))
            {
                // might be a breaking point
                diff = -1;
            } else if (diff > maxDiff and ivY != f2.vY.begin())
            {
                diff = min(diff, abs(*(ivY - 1) - *ivY2));
            }
        } else
        {
            diff = *ivYtmp - *ivY;
            if (ivX < f2.vX.end() - 1 and *ivX + 1 == *(ivX + 1))
            {
                // might be a breaking point
                diff = -1;
            } else if (diff > maxDiff and ivY < f2.vY.end() - 1)
            {
                diff = min(diff, abs(*(ivY + 1) - *ivY2));
            }
        }

        if (diff > maxDiff)
            maxDiff = diff;
    }
    return maxDiff;
}

int LPFunction::equalValue2(const LPFunction &f2) const
{
    // return the maximum difference between the two functions
    if (ID1 == -1 || ID2 == -1)
        return false;

    if (vX.size() < 2)
        return false;
    else if (f2.vX.size() < 2)
        return true;

//    assert(ID1 == f2.ID1 and ID2 == f2.ID2);

    // return true if this.vY <= f2.vY for all vX
    vector<int>::const_iterator ivY, ivYtmp;
    vector<int> vYtmp;
    vector<int>::const_iterator ivX;

    vYtmp = f2.getVY2(vX, vY); // getVY2 is used to get the Y values for the given X values of the function f2

    ivX = vX.begin();
    ivY = vY.begin();
    ivYtmp = vYtmp.begin();

    int maxDiff = 0;
    int maxY2 = max(maxY, f2.maxY);
    int bound = vX.back() - max(maxY, f2.maxY);
    for (; ivX != vX.end() && ivY != vY.end() && ivYtmp != vYtmp.end(); ivY++, ivYtmp++, ivX++)
    {
//        if (*ivX >= bound)
//            break;
        int diff = -1;
        if (*ivY > *ivYtmp)
        {
            diff = *ivY - *ivYtmp;
        } else
        {
            diff = *ivYtmp - *ivY;
        }

        if (diff > maxDiff)
        {
            maxDiff = diff;
        }
    }

    vYtmp.clear();
    vYtmp = getVY2(f2.vX, f2.vY);

    ivY = f2.vY.begin();
    ivYtmp = vYtmp.begin();
    ivX = f2.vX.begin();

    for (; ivX != f2.vX.end() && ivY != f2.vY.end() && ivYtmp != vYtmp.end() - 1; ivY++, ivYtmp++, ivX++)
    {
//        if (*ivX >= bound)
//            break;
        int diff = -1;
        if (*ivY > *ivYtmp)
        {
            diff = *ivY - *ivYtmp;
        } else
        {
            diff = *ivYtmp - *ivY;
        }

        if (diff > maxDiff)
            maxDiff = diff;
    }
    return maxDiff;
}

bool LPFunction::dominate(const LPFunction &f2) const
{
    int acc = 5;
    if (vX.size() < 2 or ID1 == -1 or ID2 == -1)
        return false;

    if (f2.vX.size() < 2 or f2.ID1 == -1 || f2.ID2 == -1)
        return true;

    if (vX.back() < f2.vX.back())
        return false;

    assert(ID1 == f2.ID1 and ID2 == f2.ID2);

    // return true if this.vY <= f2.vY for all vX
    int debugID1 = 1, debugID2 = 1;
    vector<int>::const_iterator ivY, ivYtmp;
    vector<int> vYtmp;
    vector<int>::const_iterator ivX;

    vYtmp = f2.getVY2(vX, vY); // getVY2 is used to get the Y values for the given X values of the function f2

    if (ID1 == debugID1 and ID2 == debugID2)
    {
        cout << "dominate test: \ntempY1: ";
        for (const auto &y: vYtmp)
            cout << y << " ";
        cout << endl;
    }

    ivX = vX.begin();
    ivY = vY.begin();
    ivYtmp = vYtmp.begin();

    for (; ivX != vX.end() && ivY != vY.end() && ivYtmp != vYtmp.end(); ivY++, ivYtmp++, ivX++)
    {
        if (ID1 == debugID1 and ID2 == debugID2)
            cout << *ivX << "\t" << *ivY << "\t" << *ivYtmp << endl;

        if (*ivY - acc > *ivYtmp)
        {
            return false;
        }
    }

    vYtmp.clear();
    vYtmp = getVY2(f2.vX, f2.vY);

    if (ID1 == debugID1 and ID2 == debugID2)
    {
        cout << "tempY2: ";
        for (const auto &y: vYtmp)
            cout << y << " ";
        cout << endl;
    }

    ivY = f2.vY.begin();
    ivYtmp = vYtmp.begin(); // performance of vX
    ivX = f2.vX.begin();

    for (ivY = f2.vY.begin(), ivYtmp = vYtmp.begin(), ivX = f2.vX.begin();
         ivY != f2.vY.end() && ivYtmp != vYtmp.end(); ivY++, ivYtmp++, ivX++)
    {

        if (ID1 == debugID1 and ID2 == debugID2)
            cout << *ivX << "\t" << *ivY << "\t" << *ivYtmp << endl;

        if (*ivYtmp - acc > *ivY)
            return false;
    }

    return true;
}

//void LPFunction::LPFMinSupByPlaneSweepCompPart(
//        const LPFunction &lpf, int endX, int &pos0, int &pos1, vector<int> &newVX, vector<int> &newVY,
//        vector<map<int, vector<int>>> &supPreInfo, vector<int> &cntRec,
//        pair<int, map<int, vector<int>>> &curXY, pair<int, map<int, vector<int>>> &nextXY,
//        vector<vector<int>> &orderedLPFs, vector<int> &orderOfLPFs) const
//{
//    if (lowerBound != lpf.lowerBound or upperBound != lpf.upperBound)
//    {
//        display();
//        lpf.display();
////        assert(false);
//    }
//
//    while (true)
//    {
//        int curX = curXY.first, curMinY = curXY.second.begin()->first;
//
//        map<int, vector<int>> yToInvLPFIds = curXY.second;
//
//        if (yToInvLPFIds.size() == 1 and yToInvLPFIds.begin()->second.size() == 2)
//        {
//            // two lpfs have turning points with same x value and y value
//            int x11 = vX[pos0], x12 = vX[pos0 + 1], x21 = lpf.vX[pos1], x22 = lpf.vX[pos1 + 1];
//            int y11 = vY[pos0], y12 = vY[pos0 + 1], y21 = lpf.vY[pos1], y22 = lpf.vY[pos1 + 1];
//            double y1, y2;
//            assert(x11 == x21 and y11 == y21);
//            y1 = Tools::lineGradient(x11, y11, x12, y12);
//            y2 = Tools::lineGradient(x21, y21, x22, y22);
//
//            if (y1 < y2)
//            {
//                orderedLPFs[0] = {0};
//                orderedLPFs[1] = {1};
//                orderOfLPFs[0] = 0;
//                orderOfLPFs[1] = 1;
//            } else if (y1 == y2)
//            {
//                orderedLPFs[0] = {0, 1};
//                orderedLPFs[1] = {};
//                orderOfLPFs[0] = 0;
//                orderOfLPFs[1] = 0;
//            } else
//            {
//                orderedLPFs[0] = {1};
//                orderedLPFs[1] = {0};
//                orderOfLPFs[0] = 1;
//                orderOfLPFs[1] = 0;
//            }
//
//            map<int, vector<int>> supIDs;
//            int cnt = 0;
//            for (const int &lpfId: orderedLPFs[0])
//            {
//                if (lpfId == 0)
//                {
//                    cnt += cntOfEachInt[pos0];
//                    supIDs.insert(vSupportPre[pos0].begin(), vSupportPre[pos0].end());
//                } else if (lpfId == 1)
//                {
//                    cnt += 1;
//                    supIDs.insert(lpf.vSupportPre[pos1].begin(), lpf.vSupportPre[pos1].end());
//                }
//            }
//
//            assert(!supIDs.empty());
//            if (newVX.empty() or curX > newVX.back())
//            {
//                newVX.emplace_back(curX);
//                newVY.emplace_back(curMinY);
//                supPreInfo.emplace_back(supIDs);
//                cntRec.emplace_back(cnt);
//            }
//
//            int ix = -1, iy = -1;
//            bool aInt = Tools::doIntersect(x11, y11, x12, y12, x21, y21, x22, y22, ix, iy);
//            if (aInt and ix > curX)
//            {
//                if (ix < nextXY.first)
//                {
//                    nextXY = {ix, {{iy, {-1}}}};
//                } else if (ix == nextXY.first and nextXY.second.size() == 1)
//                    nextXY = {ix, {{iy, {-2}}}};
//            }
//        } else if (yToInvLPFIds.size() == 2)
//        {
//            // two turning points with same x but different y
//            int x11 = vX[pos0], x12 = vX[pos0 + 1], x21 = lpf.vX[pos1], x22 = lpf.vX[pos1 + 1];
//            int y11 = vY[pos0], y12 = vY[pos0 + 1], y21 = lpf.vY[pos1], y22 = lpf.vY[pos1 + 1];
//
//            assert(y11 != y21);
//            if (y11 < y21)
//            {
//                orderedLPFs[0] = {0};
//                orderedLPFs[1] = {1};
//                orderOfLPFs[0] = 0;
//                orderOfLPFs[1] = 1;
//            } else
//            {
//                orderedLPFs[0] = {1};
//                orderedLPFs[1] = {0};
//                orderOfLPFs[0] = 1;
//                orderOfLPFs[1] = 0;
//            }
//
//            int lpfId = orderedLPFs[0][0];
//            map<int, vector<int>> supIDs;
//            int cnt = 0;
//            if (lpfId == 0)
//            {
//                cnt += cntOfEachInt[pos0];
//                supIDs.insert(vSupportPre[pos0].begin(), vSupportPre[pos0].end());
//            } else if (lpfId == 1)
//            {
//                cnt += 1;
//                supIDs.insert(lpf.vSupportPre[pos1].begin(), lpf.vSupportPre[pos1].end());
//            }
//
//            assert(!supIDs.empty());
//            if (newVX.empty() or curX > newVX.back())
//            {
//                newVX.emplace_back(curX);
//                newVY.emplace_back(curMinY);
//                supPreInfo.emplace_back(supIDs);
//                cntRec.emplace_back(cnt);
//            }
//
//            int ix = -1, iy = -1;
//            bool aInt = Tools::doIntersect(x11, y11, x12, y12, x21, y21, x22, y22, ix, iy);
//            if (aInt and ix > curX)
//            {
//                if (ix < nextXY.first)
//                {
//                    nextXY = {ix, {{iy, {-1}}}};
//                } else if (ix == nextXY.first and nextXY.second.size() == 1)
//                    nextXY = {ix, {{iy, {-2}}}};
//            }
//
//        } else if (yToInvLPFIds.size() == 1 and yToInvLPFIds.begin()->second.size() == 1)
//        {
//            // tuning point of one of LPF or intersection of two LPFs
//            int lpfId = yToInvLPFIds.begin()->second[0];
//            if (lpfId < 0)
//            {
//                assert(orderedLPFs[0].size() == 1 and orderedLPFs[1].size() == 1);
//                // intersection
//                vector<vector<int>> tmpOrderedLPFs = orderedLPFs;
//                orderedLPFs[0] = tmpOrderedLPFs[1];
//                orderedLPFs[1] = tmpOrderedLPFs[0];
//
//                orderOfLPFs[orderedLPFs[0][0]] = 0;
//                orderOfLPFs[orderedLPFs[1][0]] = 1;
//
//                int newSupId = orderedLPFs[0][0];
//                map<int, vector<int>> supIDs;
//                int cnt = 0;
//                if (newSupId == 0)
//                {
//                    cnt += cntOfEachInt[pos0];
//                    supIDs.insert(vSupportPre[pos0].begin(), vSupportPre[pos0].end());
//                } else if (newSupId == 1)
//                {
//                    cnt += 1;
//                    supIDs.insert(lpf.vSupportPre[pos1].begin(), lpf.vSupportPre[pos1].end());
//                }
//
//                assert(!supIDs.empty());
//                if (newVX.empty() or curX > newVX.back())
//                {
//                    newVX.emplace_back(curX);
//                    newVY.emplace_back(curMinY);
//                    supPreInfo.emplace_back(supIDs);
//                    cntRec.emplace_back(cnt);
//                }
//
//                if (lpfId == -2)
//                {
//                    int x11 = vX[pos0], x12 = vX[pos0 + 1], x21 = lpf.vX[pos1], x22 = lpf.vX[pos1 + 1];
//                    int y11 = vY[pos0], y12 = vY[pos0 + 1], y21 = lpf.vY[pos1], y22 = lpf.vY[pos1 + 1];
//                    int ix = -1, iy = -1;
//                    bool aInt = Tools::doIntersect(x11, y11, x12, y12, x21, y21, x22, y22, ix, iy);
//                    if (aInt > 0 and ix > curX)
//                    {
//                        if (ix < nextXY.first)
//                        {
//                            nextXY = {ix, {{iy, {-1}}}};
//                        } else if (ix == nextXY.first and nextXY.second.size() == 1)
//                            nextXY = {ix, {{iy, {-2}}}};
//                    }
//                }
//            } else
//            {
//                // tuning point
//                int x11 = vX[pos0], x12 = vX[pos0 + 1], x21 = lpf.vX[pos1], x22 = lpf.vX[pos1 + 1];
//                int y11 = vY[pos0], y12 = vY[pos0 + 1], y21 = lpf.vY[pos1], y22 = lpf.vY[pos1 + 1];
//
//                if (curX > x12 or curX < x11)
//                    cout << 1;
//                double y0 = Tools::getY(x11, y11, x12, y12, curX);
//                double y1 = Tools::getY(x21, y21, x22, y22, curX);
//
//                bool same = false;
//                if (y0 == y1)
//                {
//                    y0 = Tools::lineGradient(x11, y11, x12, y12);
//                    y1 = Tools::lineGradient(x21, y21, x22, y22);
//                    same = true;
//                }
//
//                if (y0 < y1)
//                {
//                    orderedLPFs[0] = {0};
//                    orderedLPFs[1] = {1};
//                    orderOfLPFs[0] = 0;
//                    orderOfLPFs[1] = 1;
//                } else if (y0 == y1)
//                {
//                    orderedLPFs[0] = {0, 1};
//                    orderedLPFs[1] = {};
//                    orderOfLPFs[0] = 0;
//                    orderOfLPFs[1] = 0;
//                } else
//                {
//                    orderedLPFs[0] = {1};
//                    orderedLPFs[1] = {0};
//                    orderOfLPFs[0] = 1;
//                    orderOfLPFs[1] = 0;
//                }
//
//                if (orderOfLPFs[lpfId] == 0 or same)
//                {
//                    map<int, vector<int>> supIDs;
//                    int cnt = 0;
//                    for (const auto &newSupId: orderedLPFs[0])
//                    {
//                        if (newSupId == 0)
//                        {
//                            cnt += cntOfEachInt[pos0];
//                            supIDs.insert(vSupportPre[pos0].begin(), vSupportPre[pos0].end());
//                        } else if (newSupId == 1)
//                        {
//                            cnt += 1;
//                            supIDs.insert(lpf.vSupportPre[pos1].begin(), lpf.vSupportPre[pos1].end());
//                        }
//                    }
//
//                    assert(!supIDs.empty());
//                    if (newVX.empty() or curX > newVX.back())
//                    {
//                        newVX.emplace_back(curX);
//                        newVY.emplace_back(curMinY);
//                        supPreInfo.emplace_back(supIDs);
//                        cntRec.emplace_back(cnt);
//                    }
//                }
//
//                int ix = -1, iy = -1;
//                bool aInt = Tools::doIntersect(x11, y11, x12, y12, x21, y21, x22, y22, ix, iy);
//                if (aInt > 0 and ix > curX and y12 != INF and y22 != INF)
//                {
//                    if (ix < nextXY.first)
//                    {
//                        nextXY = {ix, {{iy, {-1}}}};
//                    } else if (ix == nextXY.first and nextXY.second.size() == 1)
//                        nextXY = {ix, {{iy, {-2}}}};
//                }
//            }
//        } else
//        {
//            assert(false);
//        }
//
//        curXY = nextXY;
//        if (vX[pos0 + 1] <= nextXY.first)
//            pos0 += 1;
//        if (lpf.vX[pos1 + 1] <= nextXY.first)
//            pos1 += 1;
//        if (curXY.first == endX)
//        {
//            if (curXY.second.begin()->second.size() == 1 and curXY.second.begin()->second[0] < 0)
//            {
//                assert(orderedLPFs[0].size() == 1 and orderedLPFs[1].size() == 1);
//                // intersection
//                vector<vector<int>> tmpOrderedLPFs = orderedLPFs;
//                orderedLPFs[0] = tmpOrderedLPFs[1];
//                orderedLPFs[1] = tmpOrderedLPFs[0];
//
//                orderOfLPFs[orderedLPFs[0][0]] = 0;
//                orderOfLPFs[orderedLPFs[1][0]] = 1;
//            } else if (orderedLPFs[0].size() == 2)
//            {
//                double y1 = -1;
//                double y2 = -1;
//                if (vX[pos0] == endX and lpf.vX[pos1] != endX)
//                {
//                    y1 = vY[pos0];
//                    y2 = Tools::getY(lpf.vX[pos1], lpf.vY[pos1], lpf.vX[pos1 + 1], lpf.vY[pos1 + 1], endX);
//                } else if (vX[pos0] == endX and lpf.vX[pos1] == endX)
//                {
//                    y1 = vY[pos0];
//                    y2 = lpf.vY[pos1];
//                } else if (vX[pos0] != endX and lpf.vX[pos1] == endX)
//                {
//                    y1 = Tools::getY(vX[pos0], vY[pos0], vX[pos0 + 1], vY[pos0 + 1], endX);
//                    y2 = lpf.vY[pos1];
//                } else
//                {
//                    assert(false);
//                }
//
//                if (y1 < y2)
//                {
//                    orderedLPFs[0] = {0};
//                    orderedLPFs[1] = {1};
//                    orderOfLPFs[0] = 0;
//                    orderOfLPFs[1] = 1;
//                } else if (y1 == y2)
//                {
//                    orderedLPFs[0] = {0, 1};
//                    orderedLPFs[1] = {};
//                    orderOfLPFs[0] = 0;
//                    orderOfLPFs[1] = 0;
//                } else
//                {
//                    orderedLPFs[0] = {1};
//                    orderedLPFs[1] = {0};
//                    orderOfLPFs[0] = 1;
//                    orderOfLPFs[1] = 0;
//                }
//            }
//            break;
//        }
//
//        if (vX[pos0 + 1] < lpf.vX[pos1 + 1])
//            nextXY = {vX[pos0 + 1], {{vY[pos0 + 1], {0}}}};
//        else if (vX[pos0 + 1] == lpf.vX[pos1 + 1])
//        {
//            if (vY[pos0 + 1] == lpf.vY[pos1 + 1])
//                nextXY = {vX[pos0 + 1], {{vY[pos0 + 1], {0, 1}}}};
//            else
//                nextXY = {vX[pos0 + 1], {{vY[pos0 + 1], {0}}, {lpf.vY[pos1 + 1], {1}}}};
//        } else
//            nextXY = {lpf.vX[pos1 + 1], {{lpf.vY[pos1 + 1], {1}}}};
//    }
//}
void LPFunction::LPFMinSupByPlaneSweepCompPart(
        const LPFunction &lpf, int endX, int &pos0, int &pos1, vector<int> &newVX, vector<int> &newVY,
        vector<map<int, vector<int>>> &supPreInfo, vector<int> &cntRec,
        pair<int, map<int, vector<int>>> &curXY, pair<int, map<int, vector<int>>> &nextXY,
        vector<vector<int>> &orderedLPFs, vector<int> &orderOfLPFs) const
{
    if (lowerBound != lpf.lowerBound or upperBound != lpf.upperBound)
    {
        display();
        lpf.display();
//        assert(false);
    }

    while (true)
    {
        int curX = curXY.first, curMinY = curXY.second.begin()->first;

        map<int, vector<int>> yToInvLPFIds = curXY.second;

        if (yToInvLPFIds.size() == 1 and yToInvLPFIds.begin()->second.size() == 2)
        {
            // two lpfs have turning points with same x value and y value
            int x11 = vX[pos0], x12 = vX[pos0 + 1], x21 = lpf.vX[pos1], x22 = lpf.vX[pos1 + 1];
            int y11 = vY[pos0], y12 = vY[pos0 + 1], y21 = lpf.vY[pos1], y22 = lpf.vY[pos1 + 1];
            double y1, y2;
            assert(x11 == x21 and y11 == y21);
            y1 = Tools::lineGradient(x11, y11, x12, y12);
            y2 = Tools::lineGradient(x21, y21, x22, y22);

            if (y1 < y2)
            {
                orderedLPFs[0] = {0};
                orderedLPFs[1] = {1};
                orderOfLPFs[0] = 0;
                orderOfLPFs[1] = 1;
            } else if (y1 == y2)
            {
                orderedLPFs[0] = {0, 1};
                orderedLPFs[1] = {};
                orderOfLPFs[0] = 0;
                orderOfLPFs[1] = 0;
            } else
            {
                orderedLPFs[0] = {1};
                orderedLPFs[1] = {0};
                orderOfLPFs[0] = 1;
                orderOfLPFs[1] = 0;
            }

            map<int, vector<int>> supIDs;
            int cnt = 0;
            for (const int &lpfId: orderedLPFs[0])
            {
                if (lpfId == 0)
                {
                    cnt += cntOfEachInt[pos0];
                    supIDs.insert(vSupportPre[pos0].begin(), vSupportPre[pos0].end());
                } else if (lpfId == 1)
                {
                    cnt += 1;
                    supIDs.insert(lpf.vSupportPre[pos1].begin(), lpf.vSupportPre[pos1].end());
                }
            }

            assert(!supIDs.empty());
            if (newVX.empty() or curX > newVX.back())
            {
                newVX.emplace_back(curX);
                newVY.emplace_back(curMinY);
                supPreInfo.emplace_back(supIDs);
                cntRec.emplace_back(cnt);
            }

            int ix = -1, iy = -1;
            bool aInt = Tools::doIntersect(x11, y11, x12, y12, x21, y21, x22, y22, ix, iy);
            if (aInt and ix > curX)
            {
                if (ix < nextXY.first)
                {
                    nextXY = {ix, {{iy, {-1}}}};
                } else if (ix == nextXY.first and nextXY.second.size() == 1)
                    nextXY = {ix, {{iy, {-2}}}};
            }
        } else if (yToInvLPFIds.size() == 2)
        {
            // two turning points with same x but different y
            int x11 = vX[pos0], x12 = vX[pos0 + 1], x21 = lpf.vX[pos1], x22 = lpf.vX[pos1 + 1];
            int y11 = vY[pos0], y12 = vY[pos0 + 1], y21 = lpf.vY[pos1], y22 = lpf.vY[pos1 + 1];

            assert(y11 != y21);
            if (y11 < y21)
            {
                orderedLPFs[0] = {0};
                orderedLPFs[1] = {1};
                orderOfLPFs[0] = 0;
                orderOfLPFs[1] = 1;
            } else
            {
                orderedLPFs[0] = {1};
                orderedLPFs[1] = {0};
                orderOfLPFs[0] = 1;
                orderOfLPFs[1] = 0;
            }

            int lpfId = orderedLPFs[0][0];
            map<int, vector<int>> supIDs;
            int cnt = 0;
            if (lpfId == 0)
            {
                cnt += cntOfEachInt[pos0];
                supIDs.insert(vSupportPre[pos0].begin(), vSupportPre[pos0].end());
            } else if (lpfId == 1)
            {
                cnt += 1;
                supIDs.insert(lpf.vSupportPre[pos1].begin(), lpf.vSupportPre[pos1].end());
            }

            assert(!supIDs.empty());
            if (newVX.empty() or curX > newVX.back())
            {
                newVX.emplace_back(curX);
                newVY.emplace_back(curMinY);
                supPreInfo.emplace_back(supIDs);
                cntRec.emplace_back(cnt);
            }

            int ix = -1, iy = -1;
            if (y11 < 0)
            {
                cout << curX << " " << pos0 << " " << pos1 << endl;
                display();
                lpf.display();
            }

            bool aInt = Tools::doIntersect(x11, y11, x12, y12, x21, y21, x22, y22, ix, iy);
            if (aInt and ix > curX)
            {
                if (ix < nextXY.first)
                {
                    nextXY = {ix, {{iy, {-1}}}};
                } else if (ix == nextXY.first and nextXY.second.size() == 1)
                    nextXY = {ix, {{iy, {-2}}}};
            }

        } else if (yToInvLPFIds.size() == 1 and yToInvLPFIds.begin()->second.size() == 1)
        {
            // tuning point of one of LPF or intersection of two LPFs
            int lpfId = yToInvLPFIds.begin()->second[0];
            if (lpfId < 0)
            {
                assert(orderedLPFs[0].size() == 1 and orderedLPFs[1].size() == 1);
                // intersection
                vector<vector<int>> tmpOrderedLPFs = orderedLPFs;
                orderedLPFs[0] = tmpOrderedLPFs[1];
                orderedLPFs[1] = tmpOrderedLPFs[0];

                orderOfLPFs[orderedLPFs[0][0]] = 0;
                orderOfLPFs[orderedLPFs[1][0]] = 1;

                int newSupId = orderedLPFs[0][0];
                map<int, vector<int>> supIDs;
                int cnt = 0;
                if (newSupId == 0)
                {
                    cnt += cntOfEachInt[pos0];
                    supIDs.insert(vSupportPre[pos0].begin(), vSupportPre[pos0].end());
                } else if (newSupId == 1)
                {
                    cnt += 1;
                    supIDs.insert(lpf.vSupportPre[pos1].begin(), lpf.vSupportPre[pos1].end());
                }

                assert(!supIDs.empty());
                if (newVX.empty() or curX > newVX.back())
                {
                    newVX.emplace_back(curX);
                    newVY.emplace_back(curMinY);
                    supPreInfo.emplace_back(supIDs);
                    cntRec.emplace_back(cnt);
                }

                if (lpfId == -2)
                {
                    int x11 = vX[pos0], x12 = vX[pos0 + 1], x21 = lpf.vX[pos1], x22 = lpf.vX[pos1 + 1];
                    int y11 = vY[pos0], y12 = vY[pos0 + 1], y21 = lpf.vY[pos1], y22 = lpf.vY[pos1 + 1];
                    int ix = -1, iy = -1;

                    bool aInt = Tools::doIntersect(x11, y11, x12, y12, x21, y21, x22, y22, ix, iy);
                    if (aInt > 0 and ix > curX)
                    {
                        if (ix < nextXY.first)
                        {
                            nextXY = {ix, {{iy, {-1}}}};
                        } else if (ix == nextXY.first and nextXY.second.size() == 1)
                            nextXY = {ix, {{iy, {-2}}}};
                    }
                }
            } else
            {
                // tuning point
                int x11 = vX[pos0], x12 = vX[pos0 + 1], x21 = lpf.vX[pos1], x22 = lpf.vX[pos1 + 1];
                int y11 = vY[pos0], y12 = vY[pos0 + 1], y21 = lpf.vY[pos1], y22 = lpf.vY[pos1 + 1];

                if (x11 > curX or x12 < curX or x21 > curX or x22 < curX)
                {
                    cout << 2 << endl;
                    display();
                    lpf.display();
                    cout << "curX: " << curX << endl;
                    cout << "pos0: " << pos0 << endl;
                    cout << "pos1: " << pos1 << endl;
                    for (int i: newVX)
                        cout << i << " ";
                    cout << endl;
                    for (int i: newVY)
                        cout << i << " ";
                    cout << endl;
                    assert(false);
                }
                double y0 = Tools::getY(x11, y11, x12, y12, curX);
                double y1 = Tools::getY(x21, y21, x22, y22, curX);

                bool same = false;
                if (y0 == y1)
                {
                    y0 = Tools::lineGradient(x11, y11, x12, y12);
                    y1 = Tools::lineGradient(x21, y21, x22, y22);
                    same = true;
                }

                if (y0 < y1)
                {
                    orderedLPFs[0] = {0};
                    orderedLPFs[1] = {1};
                    orderOfLPFs[0] = 0;
                    orderOfLPFs[1] = 1;
                } else if (y0 == y1)
                {
                    orderedLPFs[0] = {0, 1};
                    orderedLPFs[1] = {};
                    orderOfLPFs[0] = 0;
                    orderOfLPFs[1] = 0;
                } else
                {
                    orderedLPFs[0] = {1};
                    orderedLPFs[1] = {0};
                    orderOfLPFs[0] = 1;
                    orderOfLPFs[1] = 0;
                }

                if (orderOfLPFs[lpfId] == 0 or same)
                {
                    map<int, vector<int>> supIDs;
                    int cnt = 0;
                    for (const auto &newSupId: orderedLPFs[0])
                    {
                        if (newSupId == 0)
                        {
                            cnt += cntOfEachInt[pos0];
                            supIDs.insert(vSupportPre[pos0].begin(), vSupportPre[pos0].end());
                        } else if (newSupId == 1)
                        {
                            cnt += 1;
                            supIDs.insert(lpf.vSupportPre[pos1].begin(), lpf.vSupportPre[pos1].end());
                        }
                    }

                    if (supIDs.empty())
                    {
                        display();
                        lpf.display();
                        cout << "curX: " << curX << endl;
                        assert(!supIDs.empty());
                    }

                    if (newVX.empty() or curX > newVX.back())
                    {
                        newVX.emplace_back(curX);
                        newVY.emplace_back(curMinY);
                        supPreInfo.emplace_back(supIDs);
                        cntRec.emplace_back(cnt);
                    }
                }

                int ix = -1, iy = -1;

                bool aInt = Tools::doIntersect(x11, y11, x12, y12, x21, y21, x22, y22, ix, iy);
                if (aInt > 0 and ix > curX and y12 != INF and y22 != INF)
                {
                    if (ix < nextXY.first)
                    {
                        nextXY = {ix, {{iy, {-1}}}};
                    } else if (ix == nextXY.first and nextXY.second.size() == 1)
                        nextXY = {ix, {{iy, {-2}}}};
                }
            }
        } else
        {
            assert(false);
        }

        curXY = nextXY;
        if (vX[pos0 + 1] <= nextXY.first)
            pos0 += 1;
        if (lpf.vX[pos1 + 1] <= nextXY.first)
            pos1 += 1;
        if (curXY.first == endX)
        {
            if (curXY.second.begin()->second.size() == 1 and curXY.second.begin()->second[0] < 0)
            {
                assert(orderedLPFs[0].size() == 1 and orderedLPFs[1].size() == 1);
                // intersection
                vector<vector<int>> tmpOrderedLPFs = orderedLPFs;
                orderedLPFs[0] = tmpOrderedLPFs[1];
                orderedLPFs[1] = tmpOrderedLPFs[0];

                orderOfLPFs[orderedLPFs[0][0]] = 0;
                orderOfLPFs[orderedLPFs[1][0]] = 1;
            } else if (orderedLPFs[0].size() == 2)
            {
                double y1 = -1;
                double y2 = -1;
                if (vX[pos0] == endX and lpf.vX[pos1] != endX)
                {
                    y1 = vY[pos0];
                    if (endX > lpf.vX[pos1 + 1] or endX < lpf.vX[pos1])
                    {
                        cout << 3 << endl;
                        display();
                        lpf.display();
                        cout << "curX: " << curX << endl;
                        cout << "endX: " << endX << endl;
                        cout << "pos0: " << pos0 << endl;
                        cout << "pos1: " << pos1 << endl;
                        for (int i: newVX)
                            cout << i << " ";
                        cout << endl;
                        for (int i: newVY)
                            cout << i << " ";
                        cout << endl;

                    }
                    y2 = Tools::getY(lpf.vX[pos1], lpf.vY[pos1], lpf.vX[pos1 + 1], lpf.vY[pos1 + 1], endX);
                } else if (vX[pos0] == endX and lpf.vX[pos1] == endX)
                {

                    y1 = vY[pos0];
                    y2 = lpf.vY[pos1];
                } else if (vX[pos0] != endX and lpf.vX[pos1] == endX)
                {
                    if (endX > vX[pos0 + 1] or endX < vX[pos0])
                    {
                        cout << 4 << endl;
                        display();
                        lpf.display();
                        cout << "curX: " << curX << endl;
                        cout << "endX: " << endX << endl;
                        cout << "pos0: " << pos0 << endl;
                        cout << "pos1: " << pos1 << endl;
                        for (int i: newVX)
                            cout << i << " ";
                        cout << endl;
                        for (int i: newVY)
                            cout << i << " ";
                        cout << endl;

                    }
                    y1 = Tools::getY(vX[pos0], vY[pos0], vX[pos0 + 1], vY[pos0 + 1], endX);
                    y2 = lpf.vY[pos1];
                } else
                {
                    assert(false);
                }

                if (y1 < y2)
                {
                    orderedLPFs[0] = {0};
                    orderedLPFs[1] = {1};
                    orderOfLPFs[0] = 0;
                    orderOfLPFs[1] = 1;
                } else if (y1 == y2)
                {
                    orderedLPFs[0] = {0, 1};
                    orderedLPFs[1] = {};
                    orderOfLPFs[0] = 0;
                    orderOfLPFs[1] = 0;
                } else
                {
                    orderedLPFs[0] = {1};
                    orderedLPFs[1] = {0};
                    orderOfLPFs[0] = 1;
                    orderOfLPFs[1] = 0;
                }
            }
            break;
        }

        if (vX[pos0 + 1] < lpf.vX[pos1 + 1])
            nextXY = {vX[pos0 + 1], {{vY[pos0 + 1], {0}}}};
        else if (vX[pos0 + 1] == lpf.vX[pos1 + 1])
        {
            if (vY[pos0 + 1] == lpf.vY[pos1 + 1])
                nextXY = {vX[pos0 + 1], {{vY[pos0 + 1], {0, 1}}}};
            else
                nextXY = {vX[pos0 + 1], {{vY[pos0 + 1], {0}}, {lpf.vY[pos1 + 1], {1}}}};
        } else
            nextXY = {lpf.vX[pos1 + 1], {{lpf.vY[pos1 + 1], {1}}}};
    }
}

LPFunction LPFunction::LPFMinSupByPlaneSweepNoSV(const LPFunction &lpf) const
{
    if (this->vX.size() < 2)
        return lpf;
    if (lpf.vX.size() < 2)
        return *this;

    assert(ID1 == lpf.ID1 and ID2 == lpf.ID2);

    vector<int> newVX, newVY;
    vector<map<int, vector<int>>> supPreInfo;
    supPreInfo.reserve(vX.size() * lpf.vX.size());
    vector<int> cntRec;
    cntRec.reserve(vX.size() * lpf.vX.size());

    int startX = max(vX[0], lpf.vX[0]);
    int endX = min(vX.back(), lpf.vX.back());
    assert(startX <= endX);

    int pos0 = 0, pos1 = 0;
    while (pos0 < vX.size() and vX[pos0] < startX)
    {
        newVX.emplace_back(vX[pos0]);
        newVY.emplace_back(vY[pos0]);
        if (pos0 < vSupportPre.size())
        {
            supPreInfo.emplace_back(vSupportPre[pos0]);
            cntRec.emplace_back(cntOfEachInt[pos0]);
        }
        pos0 += 1;
    }

    while (pos1 < lpf.vX.size() and lpf.vX[pos1] < startX)
    {
        newVX.emplace_back(lpf.vX[pos1]);
        newVY.emplace_back(lpf.vY[pos1]);
        if (pos1 < lpf.vSupportPre.size())
        {
            supPreInfo.emplace_back(lpf.vSupportPre[pos1]);
            cntRec.emplace_back(1);
        }
        pos1 += 1;
    }

    if (startX == endX)
    {
        if (pos0 > 0)
        {
            assert(pos1 == 0 and lpf.vX.front() == startX and vX[pos0] == startX);
            while (pos1 < lpf.vX.size())
            {
                newVX.emplace_back(lpf.vX[pos1]);
                newVY.emplace_back(lpf.vY[pos1]);
                if (pos1 < lpf.vSupportPre.size())
                {
                    supPreInfo.emplace_back(lpf.vSupportPre[pos1]);
                    cntRec.emplace_back(1);
                }
                pos1 += 1;
            }
        } else if (pos1 > 0)
        {
            assert(pos0 == 0 and vX.front() == startX and lpf.vX[pos1] == startX);
            while (pos0 < vX.size())
            {
                newVX.emplace_back(vX[pos0]);
                newVY.emplace_back(vY[pos0]);
                if (pos0 < vSupportPre.size())
                {
                    supPreInfo.emplace_back(vSupportPre[pos0]);
                    cntRec.emplace_back(cntOfEachInt[pos0]);
                }
                pos0 += 1;
            }
        }

        LPFunction fr(ID1, ID2, lpf.lowerBound, lpf.upperBound);
        fr.setValueNoComp(newVX, newVY);
        fr.vSupportPre = supPreInfo;
        return fr;
    }

    pair<int, map<int, vector<int>>> curXY, nextXY;

    if (vX.front() < lpf.vX.front())
    {
        assert(pos0 > 0 and vX[pos0] >= startX and pos1 == 0 and lpf.vX.front() == startX);
        curXY = {startX, {{lpf.vY.front(), {1}}}};
        pos0 -= 1;
    } else if (vX.front() > lpf.vX.front())
    {
        assert(pos0 == 0 and vX[0] == startX and pos1 > 0 and lpf.vX[pos1] >= startX);
        curXY = {startX, {{vY.front(), {0}}}};

        pos1 -= 1;
    } else
    {
        assert(pos0 == 0 and vX.front() == startX and pos1 == 0 and lpf.vX.front() == startX);
        if (vY.front() != lpf.vY.front())
            curXY = {startX, {{vY.front(), {0}}, {lpf.vY.front(), {1}}}};
        else
            curXY = {startX, {{vY.front(), {0, 1}}}};
    }

    vector<vector<int>> orderedLPFs(2, vector<int>());
    vector<int> orderOfLPFs(2, -1);

    assert(pos0 + 1 < vX.size() and pos1 + 1 < lpf.vX.size());
    if (vX[pos0 + 1] < lpf.vX[pos1 + 1])
        nextXY = {vX[pos0 + 1], {{vY[pos0 + 1], {0}}}};
    else if (vX[pos0 + 1] == lpf.vX[pos1 + 1])
    {
        if (vY[pos0 + 1] == lpf.vY[pos1 + 1])
            nextXY = {vX[pos0 + 1], {{vY[pos0 + 1], {0, 1}}}};
        else
            nextXY = {vX[pos0 + 1], {{vY[pos0 + 1], {0}}, {lpf.vY[pos1 + 1], {1}}}};
    } else
        nextXY = {lpf.vX[pos1 + 1], {{lpf.vY[pos1 + 1], {1}}}};

    LPFMinSupByPlaneSweepCompPart(
            lpf, endX, pos0, pos1, newVX, newVY, supPreInfo, cntRec,
            curXY, nextXY, orderedLPFs, orderOfLPFs);

    bool shouldPadding = true;
    if (vX.back() == endX and lpf.vX.back() == endX)
    {
        shouldPadding = false;
        assert(pos0 == vX.size() - 1 and pos1 == lpf.vX.size() - 1);
        newVX.push_back(endX);
        newVY.push_back(min(vY[pos0], lpf.vY[pos1]));
    } else
    {
        // lpf0 ends at endX
        if (orderedLPFs[0].size() == 1 and curXY.second.begin()->second[0] < 0)
        {
            newVX.push_back(endX);
            newVY.push_back(curXY.second.begin()->first);
            if (lpf.vX[pos1] == endX)
            {
                assert(pos0 < vSupportPre.size());
                supPreInfo.emplace_back(vSupportPre[pos0]);
                cntRec.emplace_back(cntOfEachInt[pos0]);
            } else
            {
                assert(pos1 < lpf.vSupportPre.size());
                supPreInfo.emplace_back(lpf.vSupportPre[pos1]);
                cntRec.emplace_back(lpf.cntOfEachInt[pos1]);
            }
        } else if (orderedLPFs[0].size() == 1)
        {
            if (vX[pos0] == lpf.vX[pos1] and vY[pos0] == lpf.vY[pos1])
            {
                newVX.push_back(endX);
                newVY.push_back(vY[pos0]);
                if (vX.back() == endX)
                {
                    assert(pos1 < lpf.vSupportPre.size());
                    supPreInfo.emplace_back(lpf.vSupportPre[pos1]);
                    cntRec.emplace_back(lpf.cntOfEachInt[pos1]);
                } else
                {
                    assert(pos0 < vSupportPre.size());
                    supPreInfo.emplace_back(vSupportPre[pos0]);
                    cntRec.emplace_back(cntOfEachInt[pos0]);
                }
            } else if (vX.back() == endX and orderOfLPFs[0] == 0)
            {
                shouldPadding = false;
                newVX.push_back(endX);
                newVY.push_back(vY[pos0]);
                assert(pos0 == vSupportPre.size());
                assert(pos1 < lpf.vX.size() - 1 and lpf.vX[pos1] <= endX);

            } else if (vX.back() == endX and orderOfLPFs[0] == 1 and lpf.vX[pos1] == endX)
            {
                newVX.push_back(endX);
                assert(pos1 < lpf.vSupportPre.size());
                newVY.push_back(lpf.vY[pos1]);
                supPreInfo.emplace_back(lpf.vSupportPre[pos1]);
                cntRec.emplace_back(lpf.cntOfEachInt[pos1]);

            } else if (lpf.vX.back() == endX and orderOfLPFs[1] == 0)
            {
                shouldPadding = false;
                newVX.push_back(endX);
                newVY.push_back(lpf.vY[pos1]);
                assert(pos1 == lpf.vSupportPre.size());

            } else if (lpf.vX.back() == endX and orderOfLPFs[1] == 1 and vX[pos0] == endX)
            {
                newVX.push_back(endX);
                assert(pos0 < vSupportPre.size());
                newVY.push_back(vY[pos0]);
                supPreInfo.emplace_back(vSupportPre[pos0]);
                cntRec.emplace_back(cntOfEachInt[pos0]);
            }
        }

        if (orderedLPFs[0].size() == 2)
        {
//            shouldPadding = false;
            newVX.push_back(endX);
            newVY.push_back(vY[pos0]);
            map<int, vector<int>> info;
            if (vSupportPre.size() > pos0)
                info.insert(vSupportPre[pos0].begin(), vSupportPre[pos0].end());
            if (lpf.vSupportPre.size() > pos1)
                info.insert(lpf.vSupportPre[pos1].begin(), lpf.vSupportPre[pos1].end());
            supPreInfo.push_back(info);
            cntRec.emplace_back(1);
        }

    }

    pos0 += 1;
    pos1 += 1;

    if (shouldPadding)
    {
        for (; pos0 < vX.size(); pos0++)
        {
            newVX.emplace_back(vX[pos0]);
            newVY.emplace_back(vY[pos0]);
            if (pos0 < vX.size() - 1)
            {
                supPreInfo.emplace_back(vSupportPre[pos0]);
                cntRec.emplace_back(cntOfEachInt[pos0]);
            }
        }

        for (; pos1 < lpf.vX.size(); pos1++)
        {
            newVX.emplace_back(lpf.vX[pos1]);
            newVY.emplace_back(lpf.vY[pos1]);
            if (pos1 < lpf.vX.size() - 1)
            {
                supPreInfo.emplace_back(lpf.vSupportPre[pos1]);
                cntRec.emplace_back(1);
            }
        }
    }
    assert(supPreInfo.size() == newVX.size() - 1 and cntRec.size() == supPreInfo.size());

    LPFunction fr(ID1, ID2, lpf.lowerBound, lpf.upperBound);
    fr.setValueNoComp(newVX, newVY);
    fr.vSupportPre = supPreInfo;
    fr.cntOfEachInt = cntRec;

    int debugID1 = 1, debugID2 = 1; // debugID2 -> debugID1 is wrong
    if (ID1 == debugID1 and ID2 == debugID2)
    {
        cout << "\nLPFMinSupByPlaneSweep2" << endl;
        display();
        lpf.display();
        fr.display();
    }

    return fr;
}

LPFunction LPFunction::LPFMinSupByPlaneSweepNoSV2(const LPFunction &lpf, int itvLen) const
{
    if (this->vX.size() < 2)
        return lpf;
    if (lpf.vX.size() < 2)
        return *this;

    assert(ID1 == lpf.ID1 and ID2 == lpf.ID2);

    vector<int> newVX, newVY;
    vector<map<int, vector<int>>> supPreInfo;
    vector<int> cntRec;
    newVX.reserve(vX.size() + lpf.vX.size());
    newVY.reserve(vX.size() + lpf.vX.size());
    supPreInfo.reserve(vX.size() * lpf.vX.size());
    cntRec.reserve(vX.size() * lpf.vX.size());

    int pos0 = 0, pos1 = 0;
    LPFMinBeginPos(lpf, pos0, pos1, itvLen);

    // assume this lpf is correct
    for (int i = 0; i < pos0; i++)
    {
        newVX.emplace_back(vX[i]);
        newVY.emplace_back(vY[i]);
        if (i < vSupportPre.size())
        {
            supPreInfo.emplace_back(vSupportPre[i]);
            cntRec.emplace_back(cntOfEachInt[i]);
        }
    }

    pair<int, map<int, vector<int>>> curXY, nextXY;

    int nextPos0 = pos0, nextPos1 = pos1;
    int x, y1, y2;
    if (vX[pos0] > lpf.vX[pos1 + 1] or vX[pos0] < lpf.vX[pos1])
        assert(false);
    int y = computeY(lpf.vX[pos1], lpf.vX[pos1 + 1], lpf.vY[pos1], lpf.vY[pos1 + 1], vX[pos0]);
    if (y == vY[pos0] and vX[pos0] == lpf.vX[pos1])
    {
        curXY = {vX[pos0], {{vY[pos0], {0, 1}}}};
    } else if (y >= vY[pos0])
    {
        curXY = {vX[pos0], {{vY[pos0], {0}}}};
    } else
    {
        curXY = {vX[pos0], {{y, {1}}}};
    }
    nextPos0 += 1;
    nextPos1 += 1;
//    assert(lpf.vX[pos1] >= vX[pos0]);
//    if (lpf.vX[pos1] == vX[pos0])
//        nextPos1 += 1;

    vector<vector<int>> orderedLPFs(2, vector<int>());
    vector<int> orderOfLPFs(2, -1);

    if (pos0 + 1 >= vX.size() or pos1 + 1 >= lpf.vX.size())
    {
        display();
        lpf.display();
        cout << pos0 << " " << pos1 << endl;
        assert(false);
    }
    if (vX[nextPos0] < lpf.vX[nextPos1])
        nextXY = {vX[nextPos0], {{vY[nextPos0], {0}}}};
    else if (vX[nextPos0] == lpf.vX[nextPos1])
    {
        if (vY[nextPos0] == lpf.vY[nextPos1])
            nextXY = {vX[nextPos0], {{vY[nextPos0], {0, 1}}}};
        else
            nextXY = {vX[nextPos0], {{vY[nextPos0], {0}}, {lpf.vY[nextPos1], {1}}}};
    } else
        nextXY = {lpf.vX[nextPos1], {{lpf.vY[nextPos1], {1}}}};

    LPFMinSupByPlaneSweepCompPart(
            lpf, vX.back(), pos0, pos1, newVX, newVY, supPreInfo, cntRec,
            curXY, nextXY, orderedLPFs, orderOfLPFs);

    if (pos0 != vX.size() - 1 or pos1 != lpf.vX.size() - 1)
    {
        display();
        lpf.display();
        cout << pos0 << " " << pos1 << endl;
        for (int i = 0; i < newVX.size(); i++)
        {
            cout << newVX[i] << " " << newVY[i] << endl;
        }
        assert(false);
    }
    newVX.push_back(vX.back());
    newVY.push_back(min(vY.back(), lpf.vY.back()));
    assert(supPreInfo.size() == newVX.size() - 1 and cntRec.size() == supPreInfo.size());

    if (newVX.front() != 0)
        assert(false);

    LPFunction fr(ID1, ID2, lpf.lowerBound, lpf.upperBound);
    fr.setValueNoComp(newVX, newVY);
    fr.vSupportPre = supPreInfo;
    fr.cntOfEachInt = cntRec;

    int debugID1 = 1, debugID2 = 1; // debugID2 -> debugID1 is wrong
    if (ID1 == debugID1 and ID2 == debugID2)
    {
        cout << "\nLPFMinSupByPlaneSweep2" << endl;
        display();
        lpf.display();
        fr.display();
    }

    return fr;
}

LPFunction LPFunction::LPFMinSupForDec(const LPFunction &lpf) const
{
    if (vX.size() < 2 or minY >= lpf.maxY)
        return lpf;
    if (lpf.vX.size() < 2 or maxY <= lpf.minY)
    {
        LPFunction lpf2(ID1, ID2, lpf.lowerBound, lpf.upperBound);
        lpf2.vX = vX;
        lpf2.vY = vY;
        lpf2.vSupportPre = vSupportPre;
        lpf2.cntOfEachInt = cntOfEachInt;
        return lpf2;
    }

    int debugID1 = 1, debugID2 = 1, supCheck = 1; // debugID2 -> debugID1 is wrong
    if (ID1 == debugID1 and ID2 == debugID2)
    {
        cout << "\nLPFMinSupForDec" << endl;
        display();
        lpf.display();
    }

    if (vX.back() == lpf.vX.front())
    {
        LPFunction lpf2(ID1, ID2, lpf.lowerBound, lpf.upperBound);
        lpf2.vX = vX;
        lpf2.vY = vY;
        lpf2.vSupportPre = vSupportPre;
        lpf2.cntOfEachInt = cntOfEachInt;

        for (int i = 1; i < lpf.vX.size(); i++)
        {
            lpf2.vX.emplace_back(lpf.vX[i]);
            lpf2.vY.emplace_back(lpf.vY[i]);

            lpf2.vSupportPre.emplace_back(lpf.vSupportPre[i - 1]);
            lpf2.cntOfEachInt.emplace_back(1);
        }

        return lpf2;
    } else if (vX.back() > lpf.vX.front())
    {
        LPFunction fr = LPFMinSupByPlaneSweepNoSV(lpf);
        LPFunction lpf2(ID1, ID2, lpf.lowerBound, lpf.upperBound);
        lpf2.setValue(fr.vX, fr.vY, fr.vSupportPre, fr.cntOfEachInt);
        int pos = 0;
        while (pos < vX.size() and pos < lpf2.vX.size())
        {
            if (vX[pos] == lpf2.vX[pos] and vY[pos] == lpf2.vY[pos])
                pos += 1;
            else
                break;
        }

        pos = pos > 0 ? pos - 1 : 0;
        lpf2.lowerBound = lpf.lowerBound;
        lpf2.upperBound = lpf.upperBound;

        if (ID1 == debugID1 and ID2 == debugID2)
        {
            lpf2.display();
        }
        return lpf2;
    } else
    {
        assert(false);
        return LPFunction();
    }
}

LPFunction LPFunction::LPFMinSupForExtend(const LPFunction &lpf, int itvLen) const
{
    assert(vX.size() > 1);
//    int debugID1 = 3450, debugID2 = 3315, supCheck = -1; // debugID2 -> debugID1 is wrong
    int debugID1 = 1, debugID2 = 1, supCheck = -1; // debugID2 -> debugID1 is wrong
    if (ID1 == debugID1 and ID2 == debugID2 and upperBound == -1)
    {
        cout << "\nLPFMinSupForExtend" << endl;
        display();
        lpf.display();
    }

    LPFunction lpf2 = LPFMinSupByPlaneSweepNoSV2(lpf);
    LPFunction lpf3(ID1, ID2, lpf.lowerBound, lpf.upperBound);
    lpf3.setValue(lpf2.vX, lpf2.vY, lpf2.vSupportPre, lpf2.cntOfEachInt, upperBound - itvLen);
    int pos = (int) lpf3.vX.size() - 1;
    while (pos >= 0 and lpf3.vX[pos] > upperBound - itvLen)
    {
        pos -= 1;
    }
    assert(pos >= 0);
//    if (lpf3.vX[pos] != upperBound - itvLen)
//    {
//        display();
//        lpf.display();
//        lpf3.display();
//        assert(false);
//    }
    return lpf3;
}

LPFunction LPFunction::LPFMinSupByPlaneSweep(const LPFunction &lpf) const
{
    if (vX.size() < 2)
        return lpf;
    if (lpf.vX.size() < 2)
        return *this;

    int debugID1 = 1, debugID2 = 1, supCheck = 1;
    bool test = vSupportPre[0].begin()->first == supCheck or lpf.vSupportPre[0].begin()->first == supCheck;
    if (ID1 == debugID1 and ID2 == debugID2 and lpf.vX.back() == 25200)
    {
        cout << "\nLPFMinSupByPlaneSweep3" << endl;
    }
    LPFunction fr2 = LPFMinSupByPlaneSweepNoSV(lpf);

    LPFunction fr(ID1, ID2, lpf.lowerBound, lpf.upperBound);
    if (ID1 == debugID1 and ID2 == debugID2 and lpf.vX.back() == 25200)
    {
        display();
        lpf.display();
    }

    fr.setValue(fr2.vX, fr2.vY, fr2.vSupportPre, fr2.cntOfEachInt);
    if (ID1 == debugID1 and ID2 == debugID2 and lpf.vX.back() == 25200)
        fr.display();

    return fr;
}

bool LPFunction::vSupportContains(int x) const
{
    for (const auto &p: vSupportPre)
    {
        if (p.find(x) != p.end())
            return true;
    }
    return false;
}

bool LPFunction::vSupportLastPartContains(
        int x, const unordered_map<int, unordered_map<int, CatSupRec>> &intermediateLPFs) const
{
    if (intermediateLPFs.find(ID1) == intermediateLPFs.end())
        return false;
    if (intermediateLPFs.at(ID1).find(ID2) == intermediateLPFs.at(ID1).end())
        return false;

    CatSupRec catRec = intermediateLPFs.find(ID1)->second.at(ID2);
    if (catRec.vX1.size() <= 1)
    {
//        if (ID1 == 30 and x == 40 and ID2 == 3812)
//            cout << "latPartContains " << ID1 << " " << x << " " << ID2 << endl;
        return false;
    }
    int participate = false;
    for (int k = vSupportPre.size() - 1; k >= 0; k--)
    {
        if (participate or vX[k + 1] < *(catRec.vX1.end() - 2))
            break;

        for (const auto &pp: vSupportPre[k])
        {
            if (pp.first == x)
            {
                for (const auto &itv: pp.second)
                {
                    if (itv == catRec.vX1.size() - 2 and catRec.vY.back() == vY.back())
                    {
                        participate = true;
                        break;
                    }
                }
                break;
            }
            if (participate)
                break;
        }
    }

    return participate;
}

void LPFunction::leftTrim(int lowerBound)
{
    assert(lowerBound >= vX.front());

    if (lowerBound >= vX.back())
    {
        int delNo = vX.size() - 1;
        vX.clear();
        vY.clear();
        cntOfEachInt.clear();
        vSupportPre.clear();
        minY = INF;
        maxY = -1;
    }

    vector<int>::iterator ivX, ivY;

    ivX = vX.begin();
    ivY = vY.begin();
    for (; ivX < vX.end() - 1;)
    {
        if (lowerBound <= *ivX)
            break;
        else if (lowerBound >= *(ivX + 1))
        {
            vX.erase(ivX);
            vY.erase(ivY);
        } else if (lowerBound < *(ivX + 1))
        {
            int y = computeY(*ivX, *(ivX + 1), *(ivY + 1), *ivY, lowerBound);
            vX.erase(ivX);
            vY.erase(ivY);
            vX.insert(vX.begin(), lowerBound);
            vY.insert(vY.begin(), y);
        }
    }

    assert(vX.size() > 1);

    cntOfEachInt.clear();
    cntOfEachInt.assign(vX.size() - 1, 0);
    vSupportPre.clear();
    vSupportPre.assign(vX.size() - 1, {{-1, {-1}}});

    minY = INF;
    maxY = -1;

    for (const auto y: vY)
    {
        if (minY > y)
            minY = y;
        if (maxY < y)
            maxY = y;
    }
}
