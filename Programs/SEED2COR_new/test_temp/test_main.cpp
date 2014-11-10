#include "test.h"
#include <cstdio>
#include <iostream>

struct Base {
	virtual void Print() = 0;
	virtual void Print( float input ) = 0;
};

struct X : public Base {
	float datax;
	virtual void Print() {
		std::cerr<<"Inside X"<<std::endl;
	}
	virtual void Print(float input) {
		std::cerr<<"Inside X with input "<<input<<std::endl;
	}
};

struct Y {
	float datay;
};
struct XYZ : public X, public Y {
	float dataz;
	XYZ() {}
	XYZ( const X& x ) : X(x) {}
	XYZ( const Y& y ) : Y(y) {}

	virtual void Print() {
		std::cerr<<"Inside XYZ"<<std::endl;
	}
	virtual void Print(float input) {
		std::cerr<<"Inside XYZ with input "<<input<<std::endl;
	}
};

int main() {
	Foo foo1[10];
	Foo foo2[100];
	return 0;
}
