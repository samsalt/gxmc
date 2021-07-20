#ifndef gxc_POSTPROCESSING_H
#define gxc_POSTPROCESSING_H

#include <vector>
#include <string>
#include <iostream>
#include "gxcControl.h"
#include "gxcBlock.h"
#include "gxcNode.h"
#include <exodusII.h>
#include <stdlib.h>
#include "gxcEcho.h"
#include <iomanip>

void Output(std::vector<gxcBlock> &Block, std::vector<gxcNode> &Node, gxcControl &SimulationParameter);

// void MeshfreeOutputInit(const gxcControl &SimulationParameter);

#endif
