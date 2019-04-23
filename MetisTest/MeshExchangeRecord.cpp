#include <memory>
#include "pch.h"
#include <HbxDefMacro.h>
#include <helper_hbx.h>
#include <libCUFEM\NSortMat.h>
#include "MeshExchangeRecord.h"


namespace HBXFEMDef
{


CGraphDepart::CGraphDepart()
{

}

CGraphDepart::~CGraphDepart()
{

}

mesh_t* CGraphDepart::ReadMesh(HBXFEMDef::InputRecord * _param)
{
	using namespace std;
	using namespace HBXDef;
	typedef std::map< std::string, std::shared_ptr< HBXFEMDef::NSortMat<HBXDef::UserReadPrec>> > _ElementData;
	typedef std::vector<HBXFEMDef::Node> _NodeData;

	int32_t nfields, ncon, node;
	int32_t *eptr, *eind, *ewgt;
	mesh_t* _mesh;
	void* _Element;
	_NodeData _Node;
	_mesh = CreateMesh();
	_mesh->ncon = 0;
	CheckUserDefErrors(_param->GiveField(_Element, HBXFEMDef::InputRecord::DataType_t::ALL_ELEMENT));
	_param->GiveField(_Node, "Node");
	_ElementData::iterator _iter;
	while (_iter != ((_ElementData*)_Element)->end() )
	{
		_mesh->ne = _iter->second->GetRowNum();
		eptr = _mesh->eptr = ismalloc(_mesh->ne + 1, 0, "ReadMesh: eptr");
		eind = _mesh->eind = imalloc(ntokens, "ReadMesh: eind");
		ewgt = _mesh->ewgt = ismalloc((ncon == 0 ? 1 : ncon)*mesh->ne, 1, "ReadMesh: ewgt");
	}
	


	return nullptr;
}




}