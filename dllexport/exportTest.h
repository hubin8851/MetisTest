#pragma once

#ifndef _DLL_API
#define _DLL_API _declspec(dllexport)
#else
#define _DLL_API _declspec(dllimport)
#endif

_DLL_API class exportTest
{
public:
	exportTest();
	~exportTest();

	_DLL_API void init();
};


extern exportTest& exportT;

_DLL_API exportTest &GiveClassFactory();