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
#include "gxcNodeSet.h"
// #include "gxcMeshfreeInfo.h"
#include <math.h>
#include <omp.h>


class gxcLoadingFace
{
public:
    int FaceId {};
    gxcData Area {};
    gxcData Normal[3] {};
    gxcId NodeList[4] {};
    gxcData ShapeTot[4] {};
};


class gxcLoadingFaceSet
{
public:
    std::vector<gxcLoadingFace> Face;
    gxcId FaceSetId {};
    gxcId VertexNum {};
    std::vector<std::vector<gxcData>> NaturalBoundaryConditionSet;
    std::vector<gxcData> CurrentNaturalBoundaryConditionSet;
    void UpdateSurface(std::vector<gxcNode> Node, std::vector<gxcBlock> &Block);
    void UpdateCurrentBoundaryCondition(const gxcControl &SimulationParameter);
    void UpdateFext(std::vector<std::vector<std::vector<gxcData>>> &FextTemp);
};
class gxcTargetFace
{
public:
    int FaceId {};
    int CellId {};
    gxcData Normal[3] {};
    gxcData Centroid[3] {};
    gxcData CurrentCentroid[3] {};
    gxcData CorrectForce {};
    gxcData SlaveMass {};
    gxcId NodeList[4] {};
    bool ImpactFlag {false};
    std::vector<std::vector<bool>> FaceToNodeDirection;
};
class gxcTargetFaceSet
{
public:
    std::vector<gxcTargetFace> Face;
    gxcId FaceSetId{};
    gxcId VertexNum {};
    virtual void UpdateSurface(std::vector<gxcNode> Node, const gxcControl SimulationParameter, std::vector<gxcBlock> &Block)=0;
    void InitImpactCondition(std::vector<gxcContactNodeSet> &ContactNodeSet);
    void CheckImpact(std::vector<gxcContactNodeSet> &ContactNodeSet);
    void AssignToNode(std::vector<gxcContactNodeSet> &ContactNodeSet, std::vector<gxcNode> &Node);
    void Corrector(std::vector<gxcContactNodeSet> &ContactNodeSet, std::vector<gxcNode> &Node);
};
class gxcTargetFaceSetHex: public gxcTargetFaceSet 
{
public:
    gxcTargetFaceSetHex(){VertexNum=4;};
    void UpdateSurface(std::vector<gxcNode> Node, const gxcControl SimulationParameter, std::vector<gxcBlock> &Block);
};
class gxcTargetFaceSetTet: public gxcTargetFaceSet 
{
public:
    gxcTargetFaceSetTet(){VertexNum=3;};
    void UpdateSurface(std::vector<gxcNode> Node, const gxcControl SimulationParameter, std::vector<gxcBlock> &Block);
};