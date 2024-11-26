#pragma once

include "ICurve.h"

template<class T,class V>
class BezierCurve : public ICurve<T,V>
{
public:
	BezierCurve();
	~BezierCurve();

	virtual V operator[](T u) override;
};
template<class T, class V>
BezierCurve<T, V>::BezierCurve):ICurve<T,V>() {}

template<class T, class V>
BezierCurve<T,V>::~BezierCurve() {}


template<class T, class V>
V BezierCurve<T, V>::operator[](T u)
{
	V ret =static_cast<V>(0.0);
	return ret;
}