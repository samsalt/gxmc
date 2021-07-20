#pragma once
#include <iostream>
#include <vector>
#include "gxcMath.h"
#include "gxcBlockInfo.h"
#include "gxcBlock.h"
#include "gxcNode.h"
#include "gxcCell.h"
#include "gxcControl.h"
#include "gxcEcho.h"
// #include "gxcMeshfreeInfo.h"
#include <math.h>
#include "omp.h"
#include <fstream>
// #include "gxcFace.h"

class gxcTargetFaceSet;

class gxcEssentialNodeSet
{
public:
    std::vector<gxcId> NodeId;
    gxcId NodeSetId{};
    std::vector<std::vector<gxcData>> EssentialBoundaryConditionSet;
    std::vector<bool> EBCSwitch;
    std::vector<gxcData> CurrentEssentialBoundaryCondition;
    std::vector<gxcData> CurrentVelocity;

    void UpdateCurrentBoundaryCondition(const gxcControl &SimulationParameter);
    void Predictor(std::vector<gxcNode> &Node, const double dlt);
    void Corrector(std::vector<gxcNode> &Node, const double dlt);
};
class gxcContactNode
{
public:
    gxcId NodeId;
    gxcData Position[3]{};
    gxcData ContactLength{};
    gxcData Depth{};
    gxcId ContactFaceId[2]{};
    bool ImpactFlag {false};
    gxcData CorrectForce {};
    void debugoutput()
    {
        std::ofstream debugFile("debug.dat", std::ios_base::app);
        if (debugFile)
        {
            debugFile << "NodeId " << NodeId << " Position= " << Position[0] << ", " << Position[1] << ", " << Position[2] << ", Depth= " << Depth << " ContactFaceId= " << ContactFaceId[1] << std::endl;
            debugFile.close();
        }
    }
};
class gxcContactNodeSet
{
public:
    std::vector<gxcContactNode> ContactNode;
    gxcId NodeSetId{};
    void UpdatePosition(const std::vector<gxcNode> &Node);
    void Corrector(const gxcData Dlt, std::vector<gxcTargetFaceSet*> &TargetFaceSet,  std::vector<gxcNode> &Node);
    void ForceCorrector(const gxcData Dlt, std::vector<gxcTargetFaceSet*> &TargetFaceSet,  std::vector<gxcNode> &Node);
};

