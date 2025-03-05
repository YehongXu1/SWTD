#ifndef TOOLS_H
#define TOOLS_H

#include <vector>
#include <string>
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <iostream>
#include <cmath>
#include <tgmath.h>
#include <iomanip>
#include <algorithm>
#include <iomanip>
#include <cassert>
#include <chrono>
#include <random>

//#define INF 999999999
class Tools
{
public:
    static int Round(double r);

    static double parseFloat(const std::string &input);    //string -> double

    static int parseInt(const std::string &input);    //string -> int

    static std::vector<std::string> split(const std::string &s, const std::string &seperator);    //split string into vector

    static int getdir(std::string dir, std::vector<std::string> &files);

    static double nodeRoadDistance(double nx, double ny, double rx1, double ry1, double rx2,double ry2);    //Distance from a Node to a Road

    static void GetFootOfPerpendicular(double nx, double ny, double rx1, double ry1, double rx2, double ry2, double &x, double &y);

    static bool nodeOnRoad(double nx, double ny, double rx1, double ry1, double rx2, double ry2);

    static int findCost(std::vector<int> &vCostT, std::vector<int> &vCostC, int t);

    static int hasIntersection(int x11, int y11, int x12, int y12, int x21, int y21, int x22, int y22, int &x, int &y);

    static double direction(std::pair<int, int> pi, std::pair<int, int> pj, std::pair<int, int> pk);

    static bool onSegment(std::pair<int, int> pi, std::pair<int, int> pj, std::pair<int, int> pk);

    static bool onSegment(int x1, int y1, int x2, int y2, int x3, int y3);

    static bool lineCoincide(std::pair<int, int> p1, std::pair<int, int> p2, std::pair<int, int> p3, std::pair<int, int> p4,
                      std::pair<int, int> &pr1, std::pair<int, int> &pr2);

    static double lineGradient(int x1, int y1, int x2, int y2);

    static double getY(int x1, int y1, int x2, int y2, int x);

    static double min(double x, double y);

    static double max(double x, double y);

    static int hasIntersectionDouble(double x11, double y11, double x12, double y12, double x21, double y21, double x22,
                              double y22, double &x, double &y);

    static double direction(std::pair<double, double> pi, std::pair<double, double> pj, std::pair<double, double> pk);

    static int orientation(long int x1, long int y1, long int x2, long int y2, long int x3, long int y3);

    static bool doIntersect(int x11, int y11, int x12, int y12, int x21, int y21, int x22, int y22, int &ix, int &iy);

    static bool findFile(std::string filename);

    static bool redundant(int x1, int y1, int x2, int y2, int x3, int y3);

};

template<class DT = std::chrono::milliseconds , class ClockT = std::chrono::high_resolution_clock>
class Timer
{
    using timep_t = typename ClockT::time_point;
    timep_t _start = ClockT::now(), _end = {};

public:
    void tick()
    {
        _end = timep_t{};
        _start = ClockT::now();
    }

    void tock() { _end = ClockT::now(); }

    template<class T = DT>
    auto duration() const
    {
        return std::chrono::duration_cast<T>(_end - _start);
    }
};

#endif
