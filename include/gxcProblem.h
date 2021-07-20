#pragma once
#include "gxcControl.h"
#include <vector>
#include <PostProcessing.h>
// #include "gxcMeshfreeInfo.h"
#include "gxcFace.h"
#include "gxcNodeSet.h"
#include <chrono>
#include "exodusII.h"
#include <stdlib.h>
#include "gxcBlock.h"
#include <stdio.h>
#include <fstream>
#include "gxcType.h"
#include "gxcEcho.h"
#include "HexCellFem.h"
#include "TetCellFem.h"
#include "HexCellRK.h"
#include "TetCell.h"
#include <sstream>
#include <iostream>
#include <string>
#include <cstring>
#include <iostream>
#include "gxcNode.h"
#include <cmath>
#include "Search.h"

class gxcProblem
{
public:
//data
    gxcControl SimulationParameter;
    std::vector<gxcBlock> Block;
    std::vector<gxcNode> Node;
    std::vector<gxcLoadingFaceSet> LoadingFaceSet;
    std::vector<gxcEssentialNodeSet> EssentialNodeSet;
    std::vector<gxcTargetFaceSet*> TargetFaceSet;
    std::vector<gxcContactNodeSet> ContactNodeSet;
    std::vector<std::vector<std::vector<gxcData>>> FextTemp;
    std::vector<std::vector<std::vector<gxcData>>> FintTemp;
    std::vector<std::vector<std::vector<gxcId>>> FaceTemp;
    std::vector<std::vector<gxcId>> NodeSetTemp;

    void ReadModel();
    void ReadControl();
    void CheckControl();
    void StateInitiate();
    void Assemble();
    void UpdateEnergy();
    gxcProblem ();
    void TimeIntegration();
    ~gxcProblem(){EchoStr("Simulation completed");};
};