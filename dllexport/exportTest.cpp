#include "stdafx.h"
#include "exportTest.h"


exportTest::exportTest()
{
}


exportTest::~exportTest()
{
}

void exportTest::init()
{
}

exportTest& GiveClassFactory()
{
	// TODO: �ڴ˴����� return ���
	static exportTest ans;
	return ans;
}

exportTest &exportT = GiveClassFactory();