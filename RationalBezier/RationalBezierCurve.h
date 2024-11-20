#pragma once
#include <iostream>
#include <mutex>
#include <vector>

#include "ICurve.h"

template<class T,class V>
class RationalBezierCurve : public ICurve<T,V>
{
public:
	RationalBezierCurve() = delete;
	RationalBezierCurve(std::vector<V> pts, std::vector<T> weights);
	~RationalBezierCurve();

	/*
	* Description:
	* This is O(N) time complexity algrithm for store NCR array per Rational bezier Curve.
	* This is require for per curve since number of point can change among curves.
	* T u : value to interpolate in curve [0,1]
	* return V a 3d coordinate.
	*/
	V operator[](T u) override;
private:
	static void updateFactorial(long N);
	void updateNCR(long N);
	T BernsteinPoly(long i, long N,const std::vector<T>& u, const std::vector<T>& uInverse);
private:
	std::vector<T> m_NCR;
	std::vector<T> m_weights;
	std::vector<V> m_pts;
	static std::vector<long> g_factorial;
	static std::mutex g_mutex;
};
template<class T, class V>
std::vector<long> RationalBezierCurve<T, V>::g_factorial;

template<class T, class V>
std::mutex RationalBezierCurve<T, V>::g_mutex;

template<class T,class V>
RationalBezierCurve<T,V>::RationalBezierCurve(std::vector<V> pts, std::vector<T> weights):ICurve<T,V>(), m_pts(pts), m_weights(weights)
{
}

template<class T, class V>
RationalBezierCurve<T, V>::~RationalBezierCurve() {}

/*
* this is O(N) algrithm for generate Factorial sequence. This is a static
* array and will be use within class and update as N goes higher.
*/
template<class T, class V>
void RationalBezierCurve<T, V>::updateFactorial(long N)
{
	std::lock_guard<std::mutex> lk(g_mutex);
	size_t currentSize = g_factorial.size();
	if (N < currentSize) return;

	g_factorial.resize(N+1);
	for (size_t i = currentSize; i <= N; i++)
	{
		// Handling Factorial(0)
		if (i == 0) {
			g_factorial[i] = 1;
			continue;
		}
		g_factorial[i] = g_factorial[i - 1] * i;
	}
}

template<class T, class V>
void RationalBezierCurve<T, V>::updateNCR(long N)
{
	if (g_factorial.size() <= N)
		updateFactorial(N);
	m_NCR.clear();
	m_NCR.resize(N+1);
	for (long i = 0; i <= N; i++)
	{
		m_NCR[i] = g_factorial[N] / (g_factorial[i] * g_factorial[N - i]);
	}
}

/*
* This is O(N) algrithm for get series of bernsteinPolynomial. 
* u^x and (1-u)^x arrays are inputs to this function
*/

template<class T,class V>
T RationalBezierCurve<T,V>::BernsteinPoly(long i,long N,const std::vector<T>& u,const std::vector<T>& uInverse)
{
	return m_NCR[i] * u[i] * uInverse[i];
}


/*
* This is O(N) algrithm for get interpolated C(u) for both [x,y] coordinates.
* This is written as [] operator of the Rational Brezier curve.
* Function only calculate u^x and (1-u)^x array each call and use existing NCR
* array for subseqence usage.
*/
template<class T, class V>
V RationalBezierCurve<T, V>::operator[](T u)
{
	if (u < static_cast<T>(0.0))
	{
		return m_pts[0];
	}

	if (u > static_cast<T>(1.0))
	{
		return m_pts[m_pts.size()-1];
	}

	if (m_pts.size() != m_weights.size())
	{
		std::cout << "Error, wieghts and point counts are different" << std::endl;
		return static_cast<T>(std::nan("NAN"));
	}
	// m_pts.size() give N+1 since N rational bspline has N+1 pts.
	long N = m_pts.size() -1;

	// O(N) to fill factorial
	updateFactorial(N);
	updateNCR(N);

	T denom = static_cast<T>(0.0);
	T numeratorX = static_cast<T>(0.0);
	T numeratorY = static_cast<T>(0.0);
	T numeratorZ = static_cast<T>(0.0);
	std::vector<T> uArray;
	std::vector<T> inverseU;
	uArray.resize(N + 1);
	inverseU.resize(N + 1);
	T uValue = static_cast<T>(1.0);
	T uInverseValue = static_cast<T>(1.0);

	uArray[0] = static_cast<T>(1.0);
	inverseU[N] = static_cast<T>(1.0);
	for (size_t i = 1; i < m_pts.size(); i++)
	{
		uValue *= u;
		uArray[i] = uValue;
		uInverseValue *= (1.0 - u);
		inverseU[N - i] = uInverseValue;
	}

	for (long i = 0; i < m_pts.size(); i++)
	{
		auto bernsteinMulWeight = BernsteinPoly(i, N, uArray, inverseU)* m_weights[i];
		denom += bernsteinMulWeight;
		numeratorX += bernsteinMulWeight * m_pts[i].x;
		numeratorY += bernsteinMulWeight * m_pts[i].y;
		numeratorZ += bernsteinMulWeight * m_pts[i].z;
	}

	return V(numeratorX/denom, numeratorY/denom, numeratorZ/denom);
}