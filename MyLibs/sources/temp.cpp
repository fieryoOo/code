#include "SacRec.h"

int main() {
	SacRec sac;
	//sac.DumpHD();
	sac.PrintHD("dist");

	sac.ChHdr("DisT", "138.");
	sac.ChHdr("kSTNM", "J58A");

	//std::cerr<<sac.shd;
	sac.PrintHD("dist");
}
