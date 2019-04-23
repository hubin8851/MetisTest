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
	// TODO: 在此处插入 return 语句
	static exportTest ans;
	return ans;
}

exportTest &exportT = GiveClassFactory();