#pragma once
#include <vector>
#include "gxcType.h"
// #include "gxcMeshfreeInfo.h"
#include "gxcEcho.h"

class gxcControl
{
public:
    bool MeshfreeSwitch {false};
    gxcId GridNum[3] {};     
    gxcData Dlt {};
    gxcData TimeEnd {};
    gxcData TimeCurrent {};
    gxcData TimeOutputPeriod {};
    gxcData TimetoOutput {};
    gxcId TimeStepID {1};
    gxcId TimeOutputStepID {1};
    gxcId BoundarySetNum {};
    short MaterialType {};
    short LagrangianType {};
    gxcData MaterialDensity {};
    gxcData ModelBound[2][3] {};
    gxcData GridSize[3] {};
    gxcData YoungMod, PoissonRatio, LamdaM, MuM;
    gxcData VelocityInitial[3] {};
    int NodalVarNumOutput {9};
    int CellVarNumOutput {14};
    char* CellVarNameOutput[14];
    char* NodalVarNameOutput[9];
    std::vector<std::vector<short>> NodeSetType;
    gxcData InternalEnergy {};
    gxcData ExternalEnergy {};
    gxcData KineticEnergy {};
    gxcData ThreadNumber {1};
    gxcId np {};
    gxcId nc {};
    gxcId BlockNum {};
    gxcId NodeSetNum {};
    gxcId SideSetNum {};
    std::string FileName {"etest.exo"};
    std::vector<std::vector<gxcId>> BlockCellId;


    // std::vector<gxcId> BoundaryConditionStep;
    // gxcId BoundaryConditionSetNum{};
    std::vector<std::vector<std::vector<std::vector<gxcId>>>> CellInGrid;
    // gxcMeshfreeInfo MeshfreeInfo;
    std::vector<std::vector<gxcData>> CellPosition;
};