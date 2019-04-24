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
		CGraphDepart( char* SourceFile, int _nPart = 5 );
		~CGraphDepart();

		//��ʼ����������ṹ��
		void Initial( params_t *params );

		mesh_t* ReadMesh(HBXFEMDef::InputRecord* _param);
		//ֱ��copy metis��mpmetis�ĺ���
		mesh_t* ReadMesh(params_t *params);

		//������
		InputFileResult_t SetInputData();

		//��������ļ�
		InputFileResult_t MeshPartition();

	private:
		idx_t options[METIS_NOPTIONS];
		params_t* _MyParams;//���ֲ����ṹ��
		mesh_t* _MyMesh;//mesh��ָ��
		idx_t objval;
		idx_t* epart;//��Ԫ��������
		idx_t* npart;//�ڵ㻮������

		int status; //���ֺ�������״ֵ̬
	};


	//������������һ����̬��¼ת����mesh�ṹ��
	mesh_t* DynamicRecordToMesh(HBXFEMDef::InputRecord* _param);
}






