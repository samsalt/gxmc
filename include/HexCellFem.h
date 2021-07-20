#pragma once
#include <iostream>
#include <HexCell.h>
#include <vector>
#include "gxcMath.h"
#include "gxcEcho.h"
#include "gxcControl.h"

class HexCellFem : public HexCell
{
public:
    static const gxcData ShapeHexDxi[8][3];
    static const gxcData ShapeHex8Dxi[8][8][3];
    static const gxcData ShapeHex8[8][8];
    HexCellFem();// constructer
    // gxcData ShapeHexDxi[][3];
    // void SetCell(const gxcId *vlist);
    void UpdateTotalShape(const std::vector<gxcNode> &Node,const gxcControl &SimulationParameter);
    void UpdateShapeDx(const std::vector<gxcNode> &Node, std::vector<std::vector<gxcData>> &UpShapeDxLocal, const gxcId QuadraturePointId);
};