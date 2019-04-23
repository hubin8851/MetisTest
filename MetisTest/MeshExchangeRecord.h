#pragma once

//hbx,2019
//InputRecord��ɵ�Ƶ�mesh�������ݽ����ӿں���
//��Ҫmetis��CUFEM������ͷ�ļ�

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

	//��ͼ�ۻ��ֵ�Ԫ�Ļ�����
	class CGraphDepart
	{
	public:
		CGraphDepart();
		~CGraphDepart();

		mesh_t* ReadMesh(HBXFEMDef::InputRecord* _param);
	private:
		mesh_t* _MyMesh;//mesh��ָ��
	};


	//������������һ��
	mesh_t* DynamicRecordToMesh(HBXFEMDef::InputRecord* _param);
}






