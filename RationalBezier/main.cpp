#include <iomanip>
#include <iostream>
#include <exception>

#include "RationalBezierCurve.h"

typedef float Decimal;

struct DecimalPt{
    Decimal x;
    Decimal y;
    Decimal z;

    DecimalPt():x(0.0), y(0.0), z(0.0) {}
    DecimalPt(Decimal t):x(t), y(t), z(t) {}
    DecimalPt(Decimal _x, Decimal _y, Decimal _z) :x(_x), y(_y), z(_z) {}
    DecimalPt(const DecimalPt& pt) :x(pt.x), y(pt.y), z(pt.z) {}
    DecimalPt(DecimalPt&& pt) noexcept :x(pt.x),y(pt.y),z(pt.z) {}
};


int main(int argc, char** argv)
{
    int N = 50;
    try {
        if (argc == 2)
        {
            N = std::atoi(argv[1]);
        }
    }
    catch (std::exception& e)
    {
        std::cout << "wrong arg " << argv[1] << std::endl;
        std::cout << "---HELP---" << std::endl;
        std::cout << "main.exe <number of interpolations> " << std::endl;
        return 0;
    }

    std::vector<Decimal> weight;
    std::vector<DecimalPt> points;
    weight.push_back(.1);
    weight.push_back(.1);
    weight.push_back(.1);
    weight.push_back(.1);
    weight.push_back(.1);
    weight.push_back(.1);
    weight.push_back(.1);
    points.push_back(DecimalPt(1.0f, 0.0f,0.0f));
    points.push_back(DecimalPt(2.0f, 2.0f, 1.0f));
    points.push_back(DecimalPt(3.0f, 0.0f, 2.0f));
    points.push_back(DecimalPt(4.0f, -2.0f,1.0f));
    points.push_back(DecimalPt(5.0f, 0.0f, 0.0f));
    points.push_back(DecimalPt(6.0f, 2.0f, -1.0f));
    points.push_back(DecimalPt(7.0f, 0.0f, -2.0f));
    RationalBezierCurve<Decimal, DecimalPt> spline(points,weight);

    std::cout << "Rational Bezier curve original pts: \n";
    for (int i = 0; i < points.size(); i++)
    {
        std::cout << std::fixed;
        std::cout << std::setprecision(4);
        std::cout << std::setw(8) << " Index : " << i << std::setw(12) << " : pt : (" << points[i].x << "," << points[i].y <<"," << points[i].z << ")" << std::setw(12) << weight[i] << std::endl;
    }

    std::cout << "\nRational Bezier curve interpolated pts :\n";
    for (int i = 0; i <= N; i++)
    {
        auto val1 = spline[i/static_cast<Decimal>(N)];
        std::cout << std::fixed;
        std::cout << std::setprecision(4);
        std::cout<< std::setw(8) <<" U : "<<i/ static_cast<Decimal>(N) << std::setw(12)<<" : pt : (" << val1.x << "," << val1.y << "," << val1.z <<")"<< std::endl;
    }
}