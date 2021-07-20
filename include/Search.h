#pragma once
#include <iostream>
#include <vector>
#include "gxcMath.h"
#include "gxcBlockInfo.h"
#include "gxcBlock.h"
#include "gxcNode.h"
// #include "gxcMeshfreeInfo.h"
#include "gxcEcho.h"
#include "gxcControl.h"

void UpdateGrid(std::vector<gxcNode> &Node, std::vector<gxcBlock> &Block,gxcControl &SimulationParameter);