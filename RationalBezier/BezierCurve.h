#pragma once

#include <vector>
#include <mutex>
#include <memory>
#include <iostream>

#include "ICurve.h"

template<class T, class V>
class BezierCurve : public ICurve<T, V>
{
public:
	BezierCurve(std::vector<V> pts);
	~BezierCurve();

	virtual V operator[](T u) override;

private:
	void updateFactorial(long N);
	void updateNCR(long N);
	T BernsteinPoly(long i, long N, const std::vector<T>& u, const std::vector<T>& uInverse);

private:
	std::vector<T> m_NCR;
	std::vector<V> m_pts;
	static std::vector<long> g_factorial;
	static std::mutex g_mutex;
};

template<class T, class V>
std::vector<long> BezierCurve<T, V>::g_factorial;

template<class T, class V>
std::mutex BezierCurve<T, V>::g_mutex;

template<class T, class V>
BezierCurve<T, V>::BezierCurve(std::vector<V> pts) :ICurve<T, V>(), m_pts(pts)
{
}

template<class T, class V>
BezierCurve<T, V>::~BezierCurve() {}

/*
* this is O(N) algrithm for generate Factorial sequence. This is a static
* array and will be use within class and update as N goes higher.
*/
template<class T, class V>
void BezierCurve<T, V>::updateFactorial(long N)
{
	std::lock_guard<std::mutex> lk(g_mutex);
	size_t currentSize = g_factorial.size();
	if (N < currentSize) return;

	g_factorial.resize(N + 1);
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
void BezierCurve<T, V>::updateNCR(long N)
{
	if (g_factorial.size() <= N)
		updateFactorial(N);
	m_NCR.clear();
	m_NCR.resize(N + 1);
	for (long i = 0; i <= N; i++)
	{
		m_NCR[i] = g_factorial[N] / (g_factorial[i] * g_factorial[N - i]);
	}
}

/*
* This is O(N) algrithm for get series of bernsteinPolynomial.
* u^x and (1-u)^x arrays are inputs to this function
*/

template<class T, class V>
T BezierCurve<T, V>::BernsteinPoly(long i, long N, const std::vector<T>& u, const std::vector<T>& uInverse)
{
	return m_NCR[i] * u[i] * uInverse[i];
}


template<class T, class V>
V BezierCurve<T, V>::operator[](T u)
{
	if (u < static_cast<T>(0.0))
	{
		return m_pts[0];
	}

	if (u > static_cast<T>(1.0))
	{
		return m_pts[m_pts.size() - 1];
	}

	// m_pts.size() give N+1 since N rational bspline has N+1 pts.
	long N = m_pts.size() - 1;

	// O(N) to fill factorial
	updateFactorial(N);
	updateNCR(N);

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
		auto bernstein = BernsteinPoly(i, N, uArray, inverseU);
		numeratorX += bernstein * m_pts[i].x;
		numeratorY += bernstein * m_pts[i].y;
		numeratorZ += bernstein * m_pts[i].z;
	}

	return V(numeratorX, numeratorY, numeratorZ);
}