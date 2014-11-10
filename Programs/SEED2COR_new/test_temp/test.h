#ifndef TEST_H
#define TEST_H

#include <memory>

class Foo {
public:
	Foo();
	~Foo();
private:
	struct Fooimpl;
	std::unique_ptr<Fooimpl> pimpl;
};

#endif
