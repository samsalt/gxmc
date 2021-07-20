#pragma once
#include <iostream>
#include <vector>
#include "gxcMath.h"
#include "gxcEcho.h"
#include "gxcControl.h"
#include "TetCell.h"

class TetCellFem : public TetCell
{
public:
    static const gxcData ShapeTetDxi[4][3];
    static const gxcData ShapeTet4Dxi[4][4][3];
    static const gxcData ShapeTet4[4][4];
    TetCellFem();// constructer
    // gxcData ShapeHexDxi[][3];
    // void SetCell(const gxcId *vlist);
    void UpdateTotalShape(const std::vector<gxcNode> &Node,const gxcControl &SimulationParameter);
    void UpdateShapeDx(const std::vector<gxcNode> &Node, std::vector<std::vector<gxcData>> &UpShapeDxLocal, const gxcId QuadraturePointId);
};