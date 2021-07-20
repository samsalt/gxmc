#pragma once

#include <vector>
#include "gxcType.h"
#include <array>
class gxcBlockInfo
{
    public:
    gxcId CellType;
    gxcId VertexPerCellNum;
    gxcId AttributeNum;
    gxcId CellNum;
    gxcData YoungMod, PoissonRatio, LamdaM, MuM, Density;
    gxcData VelocityInitial[3] {};
    gxcData MaterialProperty[20] {};
    short MaterialType;
    std::array<std::array<gxcData, 6>, 6> Cmat;
    void FormCmat();
    void UpdateBasicMaterialInfo();
};