#pragma once

template<class T, class V>
class ISurface
{
public:
	ISurface() {};
	virtual ~ISurface() = 0;

	virtual V at(T u,T v) = 0;
};

template<class T, class V>
ISurface<T, V>::~ISurface(){}
