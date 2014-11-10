#include "test.h"
#include <iostream>

struct Foo::Fooimpl {
public:
	Fooimpl() { oid = num_obj++; std::cerr<<oid<<" "<<num_obj<<std::endl; }
private:
	static int num_obj;
	int oid;
};

Foo::Foo() : pimpl( new Fooimpl() ) {}

Foo::~Foo() {}

int Foo::Fooimpl::num_obj = 0;
