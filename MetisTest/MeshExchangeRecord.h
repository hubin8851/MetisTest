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
		CGraphDepart( char* SourceFile, int _nPart = 5 );
		~CGraphDepart();

		//初始化传入参数结构体
		void Initial( params_t *params );

		mesh_t* ReadMesh(HBXFEMDef::InputRecord* _param);
		//直接copy metis的mpmetis的函数
		mesh_t* ReadMesh(params_t *params);

		//主函数
		InputFileResult_t SetInputData();

		//输出划分文件
		InputFileResult_t MeshPartition();

	private:
		idx_t options[METIS_NOPTIONS];
		params_t* _MyParams;//划分参数结构体
		mesh_t* _MyMesh;//mesh类指针
		idx_t objval;
		idx_t* epart;//单元划分向量
		idx_t* npart;//节点划分向量

		int status; //部分函数返回状态值
	};


	//公共函数，从一个动态记录转换成mesh结构体
	mesh_t* DynamicRecordToMesh(HBXFEMDef::InputRecord* _param);
}






