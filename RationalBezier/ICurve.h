#pragma once

template<class T, class V>
class ICurve
{
public:
	ICurve() {};
	virtual ~ICurve() = 0;

	virtual V operator[](T u) = 0;
};

template<class T, class V>
ICurve<T, V>::~ICurve() {}