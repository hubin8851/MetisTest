// MetisTest.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include "pch.h"
#include <HBXPreDef.h>
#include <boost\program_options.hpp>
#include <metis.h> 
#include "..\programs\metisbin.h"
#include <libReader\InpDataReader.h>
#include <exportTest.h>
#include "MeshExchangeRecord.h"


HBXFEMDef::InputRecord* g_Record = nullptr;

int main( int argc, char const *argv[] )
{
	exportT.init();


	if (false)
	{
// 	/************************************************************************/
// 	/*     使用boost自带的cmd描述器，测试用例								*/
// 	/*    参考 https://blog.csdn.net/Windgs_YF/article/details/81201456     */
// 	/************************************************************************/
	namespace bpo = boost::program_options;
	//步骤一: 构造选项描述器和选项存储器
	// 选项描述器, 其参数为该描述器的名字
		bpo::options_description opts("all options");
		//选项存储器,继承自map容器
		bpo::variables_map vm;

		//步骤二: 为选项描述器增加选项,其参数依次为: key, value的类型，该选项的描述
		opts.add_options()
			("filename", bpo::value<std::string>(), "the file name which want to be found")
			("help", "this is a program to find a specified file");

		try
		{
			bpo::store(bpo::parse_command_line(argc, argv, opts), vm);
		}
		catch (...)
		{
			std::cout << "输入的参数中存在未定义的选项！\n";
			return 0;
		}
	}


	if (1)
	{
		using namespace HBXFEMDef;
		BaseReader* g_EltReader = InstanceElemtPropReader();
		g_EltReader->SetSourceFilePath("EltProperty.xml", "F:\\data_from_HBX_phd\\vs2015\\FEM_CUDA_Boost_V2\\");
		InpDataReader g_InpDataReader( "B31BEAMTEST.inp", "F:\\data_from_HBX_phd\\database\\B31Test\\" );
		g_InpDataReader.SetEltPropFilePath("EltProperty.xml", "F:\\data_from_HBX_phd\\vs2015\\FEM_CUDA_Boost_V2\\");
		g_InpDataReader.SetInputData();
		g_Record = g_InpDataReader.GetInputRecord();
	}

	HBXFEMDef::CGraphDepart g_Depart( "E:\\Desktop\\metis-5.1.0\\graphs\\2D4NShell.mesh", 5 );
	g_Depart.ReadMesh(g_Record);

	graph_t* graph;		//生成的无向图


	int status = 0;	//状态值，当前选用

}

// 运行程序: Ctrl + F5 或调试 >“开始执行(不调试)”菜单
// 调试程序: F5 或调试 >“开始调试”菜单

// 入门提示: 
//   1. 使用解决方案资源管理器窗口添加/管理文件
//   2. 使用团队资源管理器窗口连接到源代码管理
//   3. 使用输出窗口查看生成输出和其他消息
//   4. 使用错误列表窗口查看错误
//   5. 转到“项目”>“添加新项”以创建新的代码文件，或转到“项目”>“添加现有项”以将现有代码文件添加到项目
//   6. 将来，若要再次打开此项目，请转到“文件”>“打开”>“项目”并选择 .sln 文件
