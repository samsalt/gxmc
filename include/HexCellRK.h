#pragma once
#include <iostream>
#include <HexCell.h>
#include <vector>
#include "gxcControl.h"
#include "gxcEcho.h"

class HexCellRK : public HexCell
{
    public:
    void UpdateTotalShape(const std::vector<gxcNode> &Node,const gxcControl &SimulationParameter) {};
    void UpdateShapeDx(const std::vector<gxcNode> &Node, std::vector<std::vector<gxcData>> &UpShapeDxLocal, const gxcId QuadraturePointId) {};
    void FindNeighbor(const gxcControl &SimulationParameter) {};
    void HardSearch(const std::vector<gxcId> &NeighborList,const gxcControl& Simulationparameter);
};