#include "pch.h"
#include <memory>
#include "metis.h"
#include <HbxDefMacro.h>
#include <helper_hbx.h>
#include <libCUFEM\NSortMat.h>
#include "MeshExchangeRecord.h"


namespace HBXFEMDef
{


CGraphDepart::CGraphDepart( char* SourceFile, int _nPart ):_MyParams(nullptr), _MyMesh(nullptr), status(0)
{
	_MyParams = (params_t*)gk_malloc(sizeof(params_t), "parse_cmdline");
	memset((void *)_MyParams, 0, sizeof(params_t));
	/* initialize the params data structure */
	_MyParams->gtype = METIS_GTYPE_DUAL;
	_MyParams->ptype = METIS_PTYPE_KWAY;
	_MyParams->objtype = METIS_OBJTYPE_CUT;
	_MyParams->ctype = METIS_CTYPE_SHEM;
	_MyParams->iptype = METIS_IPTYPE_GROW;
	_MyParams->rtype = -1;

	_MyParams->minconn = 0;
	_MyParams->contig = 0;

	_MyParams->nooutput = 1;//暂时不要输出
	_MyParams->wgtflag = 3;

	_MyParams->ncuts = 1;
	_MyParams->niter = 10;
	_MyParams->ncommon = 1;

	_MyParams->dbglvl = 0;
	_MyParams->balance = 0;
	_MyParams->seed = -1;
	_MyParams->dbglvl = 0;

	_MyParams->tpwgtsfile = NULL;

	_MyParams->filename = NULL;
	_MyParams->nparts = 1;

	_MyParams->ufactor = -1;

	gk_clearcputimer(_MyParams->iotimer);
	gk_clearcputimer(_MyParams->parttimer);
	gk_clearcputimer(_MyParams->reporttimer);

	_MyParams->filename = SourceFile;
	_MyParams->nparts = _nPart;

	if (_MyParams->nparts < 2)
		errexit("The number of partitions should be greater than 1!\n");

	/* Set the ptype-specific defaults */
	if (_MyParams->ptype == METIS_PTYPE_RB) {
		_MyParams->rtype = METIS_RTYPE_FM;
	}
	if (_MyParams->ptype == METIS_PTYPE_KWAY) {
		_MyParams->iptype = METIS_IPTYPE_METISRB;	//初始分区方案
		_MyParams->rtype = METIS_RTYPE_GREEDY; //优化方案
	}

#ifdef _DEBUG DEBUG下做检查
	/* Check for invalid parameter combination */
	if (_MyParams->ptype == METIS_PTYPE_RB) {
		if (_MyParams->contig)
			errexit("The -contig option cannot be specified with rb partitioning.\n");
		if (_MyParams->minconn)
			errexit("The -minconn option cannot be specified with rb partitioning.\n");
		if (_MyParams->objtype == METIS_OBJTYPE_VOL)
			errexit("The -objtype=vol option cannot be specified with rb partitioning.\n");
	}
#endif
}

CGraphDepart::~CGraphDepart()
{
	FreeMesh(&_MyMesh);

	gk_free((void **)&epart, &npart, LTERM);
	gk_free((void **)&_MyParams->filename, &_MyParams->tpwgtsfile, &_MyParams->tpwgts,
		&_MyParams->ubvecstr, &_MyParams->ubvec, &_MyParams, LTERM);
}

void CGraphDepart::Initial( params_t *params )
{
	FreeMesh(&_MyMesh);

	if (nullptr != params)
	{
		_MyParams = params;
	}
	return;
}

mesh_t* CGraphDepart::ReadMesh(HBXFEMDef::InputRecord * _param)
{
	gk_startcputimer(_MyParams->iotimer);

	_MyMesh = DynamicRecordToMesh(_param);

	ReadTPwgts(_MyParams, _MyMesh->ncon);
	gk_stopcputimer(_MyParams->iotimer);

	return _MyMesh;
}

mesh_t * CGraphDepart::ReadMesh(params_t * params)
{
	gk_startcputimer(_MyParams->iotimer);

	idx_t i, j, k, l, nfields, ncon, node;
	idx_t *eptr, *eind, *ewgt;
	size_t nlines, ntokens;
	char *line = NULL, *curstr, *newstr;
	size_t lnlen = 0;
	FILE *fpin;

	if (!gk_fexists(params->filename))
		errexit("File %s does not exist!\n", params->filename);

	_MyMesh = CreateMesh();

	/* get some file stats */
	gk_getfilestats(params->filename, &nlines, &ntokens, NULL, NULL);

	fpin = gk_fopen(params->filename, "r", __func__);

	/* Skip comment lines until you get to the first valid line */
	do {
		if (gk_getline(&line, &lnlen, fpin) == -1)
			errexit("Premature end of input file: file: %s\n", params->filename);
	} while (line[0] == '%');


	_MyMesh->ncon = 0;
	nfields = sscanf(line, "%d %d", &(_MyMesh->ne), &(_MyMesh->ncon));

	if (nfields < 1)
		std::cerr<<"The input file does not specify the number of elements."<<std::endl;

	if (_MyMesh->ne <= 0)
		std::cerr<<"The supplied number of elements:" << _MyMesh->ne << "must be positive," <<std::endl;

	if (_MyMesh->ne > nlines)
		std::cerr << "The file has " << nlines << "lines which smaller than the number of elements of"
			<< _MyMesh->ne << "specified in the header line."<<std::endl;

	ncon = _MyMesh->ncon;
	eptr = _MyMesh->eptr = ismalloc(_MyMesh->ne + 1, 0, "ReadMesh: eptr");
	eind = _MyMesh->eind = imalloc(ntokens, "ReadMesh: eind");
	ewgt = _MyMesh->ewgt = ismalloc((ncon == 0 ? 1 : ncon)*_MyMesh->ne, 1, "ReadMesh: ewgt");


	/*----------------------------------------------------------------------
	 * Read the mesh file
	 *---------------------------------------------------------------------*/
	for (eptr[0] = 0, k = 0, i = 0; i < _MyMesh->ne; i++) {
		do {
			if (gk_getline(&line, &lnlen, fpin) == -1)
				std::cerr << "Premature end of input file while reading element " << i + 1 << std::endl;
		} while (line[0] == '%');

		curstr = line;
		newstr = NULL;

		/* Read element weights */
		for (l = 0; l < ncon; l++) {
			ewgt[i*ncon + l] = strtoidx(curstr, &newstr, 10);
			if (newstr == curstr)
				std::cerr << "The line for vertex" <<i + 1 <<" does not have enough weights " << std::endl;
			if (ewgt[i*ncon + l] < 0)
				std::cerr << "The weight for element "<< i + 1 <<" and constraint %"<< l <<" must be >= 0" << std::endl;
			curstr = newstr;
		}

		while (1) {
			node = strtoidx(curstr, &newstr, 10);
			if (newstr == curstr)
				break; /* End of line */
			curstr = newstr;

			if (node < 1)
				std::cerr << "Node "<< node <<" for element "<< i + 1 <<" is out of bounds\n"<<std::endl;

			eind[k++] = node - 1;
		}
		eptr[i + 1] = k;
	}
	gk_fclose(fpin);

	_MyMesh->ncon = (ncon == 0 ? 1 : ncon);
	_MyMesh->nn = imax(eptr[_MyMesh->ne], eind) + 1;

//	gk_free((void *)&line, LTERM);

	ReadTPwgts(_MyParams, _MyMesh->ncon);
	gk_stopcputimer(_MyParams->iotimer);

	return _MyMesh;
}

InputFileResult_t CGraphDepart::SetInputData()
{
	MPPrintInfo(_MyParams, _MyMesh);
	
	epart = imalloc(_MyMesh->ne, "main: epart");
	npart = imalloc(_MyMesh->nn, "main: npart");

	METIS_SetDefaultOptions(options);
	options[METIS_OPTION_PTYPE] = _MyParams->ptype;
	options[METIS_OPTION_OBJTYPE] = _MyParams->objtype;
	options[METIS_OPTION_CTYPE] = _MyParams->ctype;
	options[METIS_OPTION_IPTYPE] = _MyParams->iptype;
	options[METIS_OPTION_RTYPE] = _MyParams->rtype;
	options[METIS_OPTION_DBGLVL] = _MyParams->dbglvl;
	options[METIS_OPTION_UFACTOR] = _MyParams->ufactor;
	options[METIS_OPTION_MINCONN] = _MyParams->minconn;
	options[METIS_OPTION_CONTIG] = _MyParams->contig;
	options[METIS_OPTION_SEED] = _MyParams->seed;
	options[METIS_OPTION_NITER] = _MyParams->niter;
	options[METIS_OPTION_NCUTS] = _MyParams->ncuts;

	gk_malloc_init();
	gk_startcputimer(_MyParams->parttimer);
	switch (_MyParams->gtype) {
	case METIS_GTYPE_DUAL:
		status = METIS_PartMeshDual(&_MyMesh->ne, &_MyMesh->nn, _MyMesh->eptr, _MyMesh->eind,
			_MyMesh->ewgt, NULL, &_MyParams->ncommon, &_MyParams->nparts,
			_MyParams->tpwgts, options, &objval, epart, npart);
		break;

	case METIS_GTYPE_NODAL:
		status = METIS_PartMeshNodal(&_MyMesh->ne, &_MyMesh->nn, _MyMesh->eptr, _MyMesh->eind,
			NULL, NULL, &_MyParams->nparts, _MyParams->tpwgts, options, &objval,
			epart, npart);
		break;
	}

	gk_stopcputimer(_MyParams->parttimer);

	return InputFileResult_t::IRRT_OK;
}

InputFileResult_t CGraphDepart::MeshPartition()
{
#ifdef _DEBUG
	if (gk_GetCurMemoryUsed() != 0)
	{
		printf("***It seems that Metis did not free all of its memory! Report this.\n");
		InputFileResult_t::IRRT_BAD_FORMAT;
	}
	_MyParams->maxmemory = gk_GetMaxMemoryUsed();
	gk_malloc_cleanup(0);
	if (status != METIS_OK) {
		printf("\n***Metis returned with an error.\n");
		InputFileResult_t::IRRT_BAD_FORMAT;
	}
	else {
#endif // _DEBUG
		if (!_MyParams->nooutput) {
			/* 写入数据 */
			gk_startcputimer(_MyParams->iotimer);
			WriteMeshPartition(_MyParams->filename, _MyParams->nparts, _MyMesh->ne, epart, _MyMesh->nn, npart);
			gk_stopcputimer(_MyParams->iotimer);
		}
		MPReportResults(_MyParams, _MyMesh, epart, npart, objval);
#ifdef _DEBUG
	}
#endif // _DEBUG

	return InputFileResult_t::IRRT_OK;
}




mesh_t * DynamicRecordToMesh(HBXFEMDef::InputRecord * _param)
{
	using namespace std;
	using namespace HBXDef;
	typedef std::map< std::string, std::shared_ptr< HBXFEMDef::NSortMat<HBXDef::UserReadPrec>> > _ElementData;
	typedef std::vector<HBXFEMDef::Node> _NodeData;

	int32_t k;
	int32_t nfields, ncon, node;
	int32_t *eptr, *eind, *ewgt;
	mesh_t* _mesh;
	void* _Element = nullptr;
	_NodeData _Node;
	_mesh = CreateMesh();
	_mesh->ncon = 0;
	CheckUserDefErrors(_param->GiveField(_Element, HBXFEMDef::InputRecord::DataType_t::ALL_ELEMENT));
	_param->GiveField(_Node, "Node");
	_ElementData::iterator _iter;
	while (_iter != ((_ElementData*)_Element)->end())
	{
		std::shared_ptr< HBXFEMDef::NSortMat<HBXDef::UserReadPrec>> _pMeshMat(_iter->second);
		_mesh->ne = _iter->second->GetRowNum();
		eptr = _mesh->eptr = ismalloc(_mesh->ne + 1, 0, "ReadMesh: eptr");
		eptr[0] = 0;
		//假定所有的元素就是单元个数*单元节点数+1
		eind = _mesh->eind = imalloc( (_mesh->ne)*(_iter->second->GetColNum())+1, "ReadMesh: eind");
		ewgt = _mesh->ewgt = ismalloc((_mesh->ncon == 0 ? 1 : _mesh->ncon)*_mesh->ne, 1, "ReadMesh: ewgt");

		//读取mesh数据
		for (size_t i = 0, k = 0; i < _mesh->ne; i++)
		{
			//不用getline函数读取每一行
			//不用读取单元权重
			for (size_t j = 0; j < _iter->second->GetColNum(); j++)
			{

				eind[k++] = (*_pMeshMat)(i, j);
			}
			eptr[i + 1] = k;//第i个单元在该向量第i+1位置上记录当前所有节点数的累加和
		}
		//不用关文件

		_mesh->ncon = (_mesh->ncon == 0 ? 1 : _mesh->ncon);
		_mesh->nn = imax(eptr[_mesh->ne], eind) + 1;

	}

	return _mesh;
}



}