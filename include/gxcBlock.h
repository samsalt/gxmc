#ifndef gxc_gxcBLOCK_H
#define gxc_gxcBLOCK_H

#include <vector>
#include <gxcCell.h>
#include "gxcType.h"
#include <gxcBlockInfo.h>

class gxcBlock
{
    public:
    gxcBlockInfo Info;
    std::vector<gxcCell*> Cell;
};
#endif