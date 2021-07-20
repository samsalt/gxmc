#pragma once
#include <iostream>
#include <vector>
#include "gxcMath.h"
#include "gxcBlockInfo.h"
#include "gxcNode.h"
#include "gxcControl.h"
// #include "gxcMeshfreeInfo.h"
#include "gxcEcho.h"
#include "math.h"
#include "stdlib.h"
#include "omp.h"

using namespace std;

class gxcCell
{
public:
    // short GetCellType() {return CellType;}
    gxcData GetVolume() { return Volume; }
    // virtual void SetCell(const gxcId *vlist) =0;
    virtual void UpdateVolume(const std::vector<gxcNode> Node) = 0;
    void inline Constitution(std::array<std::array<gxcData, 3>, 3> &SPK, const std::array<std::array<gxcData, 3>, 3> &GreenStrain, const gxcBlockInfo &BlockInfo, gxcId QuadraturePointId);
    void UpdateFintTot(std::vector<gxcNode> &Node, const gxcBlockInfo &BlockInfo,std::vector<std::vector<std::vector<gxcData>>> &FintTemp);
    void UpdateFintUp(std::vector<gxcNode> &Node, const gxcBlockInfo &BlockInfo, const gxcData Dlt,std::vector<std::vector<std::vector<gxcData>>> &FintTemp);
    void UpdateLumpedMass(std::vector<gxcNode> &Node, const gxcBlockInfo &BlockInfo);
    void inline GetLocalDisplacement(const std::vector<gxcNode> &Node, std::vector<std::vector<gxcData>> &DisplacementLocal);
    void inline GetLocalVelocity(const std::vector<gxcNode> &Node, std::vector<std::vector<gxcData>> &VelocityLocal);
    void inline GetLocalTotalShapeDx(const gxcId NodeNum, std::vector<std::vector<gxcData>> &TotalShapeDxLocal, const gxcId QuadraturePointId);
    virtual void UpdateTotalShape(const std::vector<gxcNode> &Node,const gxcControl &SimulationParameter) = 0;
    virtual void UpdateShapeDx(const std::vector<gxcNode> &Node, std::vector<std::vector<gxcData>> &UpShapeDxLocal, const gxcId QuadraturePointId) = 0;
    void AddVertexId(gxcId i, gxcId vid) { VertexId[i] = vid; }
    gxcId GetVertexID(gxcId i) { return VertexId[i]; }
    void SaveStress(const std::array<std::array<gxcData, 3>, 3> &mat);
    void SaveStress(const std::array<std::array<gxcData, 3>, 3> &mat, const gxcId QuadraturePointId);
    void inline GetCauchyMat(std::array<std::array<gxcData, 3>, 3> &mat,const gxcId QuadraturePointId);
    void SaveStrain(const std::array<std::array<gxcData, 1>, 6> &StrainVoigt, const gxcId QuadraturePointId);
    void SetGlobalId(const gxcId IdTemp) {GlobalId=IdTemp;}
    gxcId GetGlobalId() {return GlobalId;}
    void FindGrid(gxcControl &SimulationParameter);
    virtual void FindNeighbor(const gxcControl &SimulationParameter) {};
    gxcData GetPosition(gxcId dim) {return Position[dim];}

    //Data

    // gxcData Strain[6]{};
    short   SideNum {};
    gxcId VertexId[20];
    // gxcData State[10] {};
    std::vector<std::vector<gxcData>> State;    //0, eps;
    std::vector<std::vector<gxcData>> Stress;
    std::vector<std::vector<gxcData>> Strain;

protected:
    gxcData Volume;
    std::vector<gxcData> Jacobian;
    gxcData Position[3] {};


    std::vector<std::vector<gxcData>> TotalShapeSurface;
    std::vector<std::vector<std::vector<gxcData>>> TotalShapeDx;
    std::vector<std::vector<gxcData>> TotalShape;
    
    short   NodeNum {};
    gxcId NeighborNode[40];
    gxcId GridIndex[3]  {};
    gxcId GlobalId {};
    gxcId QuadratureNum {};
    gxcData QWeight {};
    gxcData ParentVolume {};
};