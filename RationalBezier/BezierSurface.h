#pragma once
#include <vector>
#include <mutex>
#include <iostream>

#include "ISurface.h"

template < class T, class V>
class BezierSurface : public ISurface<T, V>
{
public:
	BezierSurface(std::vector<V> uPts, long numX, long numV);
	~BezierSurface();
	V at(T u, T v) override;

private:
	void updateFactorial(long N);
	void updateNCR(long N);
	T BernsteinPoly(long i, long N, const std::vector<T>& u, const std::vector<T>& uInverse);

private:
	std::vector<T> m_NCR;
	std::vector<V> m_Pts;
	long m_numX;
	long m_numY;
	static std::vector<long> g_factorial;
	static std::mutex g_mutex;
};

template<class T, class V>
std::vector<long> BezierSurface<T, V>::g_factorial;

template<class T, class V>
std::mutex BezierSurface<T, V>::g_mutex;

template<class T, class V>
BezierSurface<T, V>::BezierSurface(std::vector<V> uPts, long numX, long numY) :ISurface<T, V>(), m_Pts(uPts), m_numX(numX), m_numY(numY)
{
}

template<class T, class V>
BezierSurface<T, V>::~BezierSurface() {}

/*
* this is O(N) algrithm for generate Factorial sequence. This is a static
* array and will be use within class and update as N goes higher.
*/
template<class T, class V>
void BezierSurface<T, V>::updateFactorial(long N)
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
void BezierSurface<T, V>::updateNCR(long N)
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
T BezierSurface<T, V>::BernsteinPoly(long i, long N, const std::vector<T>& u, const std::vector<T>& uInverse)
{
	return m_NCR[i] * u[i] * uInverse[i];
}


template<class T, class V>
V BezierSurface<T, V>::at(T u, T v)
{
	if (u < static_cast<T>(0.0))
	{
		u = static_cast<T>(0.0);
	}
	if (v < static_cast<T>(1.0))
	{
		v = static_cast<T>(1.0);
	}

	if (u > static_cast<T>(1.0))
	{
		u = static_cast<T>(1.0);
	}
	if (v > static_cast<T>(1.0))
	{
		v = static_cast<T>(1.0);
	}

	if (m_Pts.size() == 0)
	{
		std::cout << "Error, no controll ponts" << std::endl;
		return V(static_cast<T>(std::nan("NAN")), static_cast<T>(std::nan("NAN")), static_cast<T>(std::nan("NAN")));
	}
	// m_pts.size() give N+1 since N rational bspline has N+1 pts.
	long M = m_numX - 1;
	long N = m_numY - 1;

	// O(N) to fill factorial
	long K = std::max(M, N);
	updateFactorial(K);
	updateNCR(K);

	T numeratorX = static_cast<T>(0.0);
	T numeratorY = static_cast<T>(0.0);
	T numeratorZ = static_cast<T>(0.0);
	std::vector<T> uArray;
	std::vector<T> inverseU;
	std::vector<T> vArray;
	std::vector<T> inverseV;

	uArray.resize(M + 1);
	inverseU.resize(M + 1);
	vArray.resize(N + 1);
	inverseV.resize(N + 1);

	T uValue = static_cast<T>(1.0);
	T uInverseValue = static_cast<T>(1.0);
	T vValue = static_cast<T>(1.0);
	T vInverseValue = static_cast<T>(1.0);

	uArray[0] = static_cast<T>(1.0);
	inverseU[M] = static_cast<T>(1.0);
	vArray[0] = static_cast<T>(1.0);
	inverseV[N] = static_cast<T>(1.0);
	for (size_t i = 1; i < m_numX; i++)
	{
		uValue *= u;
		uArray[i] = uValue;
		uInverseValue *= (1.0 - u);
		inverseU[M - i] = uInverseValue;
	}

	for (size_t i = 1; i < m_numY; i++)
	{
		vValue *= v;
		vArray[i] = vValue;
		vInverseValue *= (1.0 - v);
		inverseV[N - i] = vInverseValue;
	}
	std::vector<T> bernsteinU;
	bernsteinU.resize(M + 1);
	std::vector<T> bernsteinV;
	bernsteinV.resize(N + 1);
	for (long i = 0; i < m_numX; i++)
	{
		bernsteinU[i] = BernsteinPoly(i, N, uArray, inverseU);
	}

	for (long i = 0; i < m_numY; i++)
	{
		bernsteinV[i] = BernsteinPoly(i, N, vArray, inverseV);
	}

	for (long i = 0; i < m_numX; i++)
	{
		T numeratorTempX = static_cast<T>(0.0);
		T numeratorTempY = static_cast<T>(0.0);
		T numeratorTempZ = static_cast<T>(0.0);
		for (long j = 0; j < m_numY; j++)
		{
			numeratorTempX += bernsteinV[j] * m_Pts[i * m_numY + j].x;
			numeratorTempY += bernsteinV[j] * m_Pts[i * m_numY + j].y;
			numeratorTempZ += bernsteinV[j] * m_Pts[i * m_numY + j].z;
		}
		numeratorX += bernsteinU[i] * numeratorTempX;
		numeratorY += bernsteinU[i] * numeratorTempY;
		numeratorZ += bernsteinU[i] * numeratorTempZ;
	}

	return V(numeratorX, numeratorY, numeratorZ);
}