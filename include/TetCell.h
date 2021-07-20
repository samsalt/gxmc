#pragma once

#include <iostream>
#include "gxcCell.h"
#include <vector>
#include "gxcControl.h"

class TetCell: public gxcCell
{
    public:
    // void SetCell(const gxcId *vlist);
    void UpdateVolume(const std::vector<gxcNode> Node);
    virtual void UpdateTotalShape(const std::vector<gxcNode> &Node,const gxcControl &SimulationParameter)=0;
    virtual void UpdateShapeDx(const std::vector<gxcNode> &Node, std::vector<std::vector<gxcData>> &UpShapeDxLocal, const gxcId QuadraturePointId)=0;
};