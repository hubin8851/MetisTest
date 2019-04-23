#pragma once

//hbx,2019
//InputRecord和傻逼的mesh进行数据交互接口函数
//需要metis和CUFEM的两者头文件

#include <GKlib.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <limits.h>
#include <signal.h>
#include <setjmp.h>
#include <assert.h>

#include <libCUFEM\InputRecord.h>
extern "C" {
#include <metis.h>
#include "..\programs\metisbin.h"
}


namespace HBXFEMDef
{

	//以图论划分单元的划分类
	class CGraphDepart
	{
	public:
		CGraphDepart();
		~CGraphDepart();

		mesh_t* ReadMesh(HBXFEMDef::InputRecord* _param);
	private:
		mesh_t* _MyMesh;//mesh类指针
	};


	//公共函数，从一个
	mesh_t* DynamicRecordToMesh(HBXFEMDef::InputRecord* _param);
}






