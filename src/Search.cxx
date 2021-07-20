#include "Search.h"
// void UpdateGrid(   std::vector<gxcNode> &Node, std::vector<gxcBlock> &Block, gxcControl &SimulationParameter)
// {
//     gxcData XMax, YMax, ZMax, XMin, YMin, ZMin;
//     XMax = Node[0].InitCoordinate[0];
//     YMax = Node[0].InitCoordinate[1];
//     ZMax = Node[0].InitCoordinate[2];
//     XMin = Node[0].InitCoordinate[0];
//     YMin = Node[0].InitCoordinate[1];
//     ZMin = Node[0].InitCoordinate[2];
//     for (gxcId i = 0; i < SimulationParameter.np; i++)
//     {
//         if (Node[i].InitCoordinate[0] < XMin)
//             XMin = Node[i].InitCoordinate[0];
//         if (Node[i].InitCoordinate[0] > XMax)
//             XMax = Node[i].InitCoordinate[0];
//         if (Node[i].InitCoordinate[1] < YMin)
//             YMin = Node[i].InitCoordinate[1];
//         if (Node[i].InitCoordinate[1] > YMax)
//             YMax = Node[i].InitCoordinate[1];
//         if (Node[i].InitCoordinate[2] < ZMin)
//             ZMin = Node[i].InitCoordinate[2];
//         if (Node[i].InitCoordinate[2] > ZMax)
//             ZMax = Node[i].InitCoordinate[2];
//     }
//     SimulationParameter.ModelBound[0][0] = XMin;
//     SimulationParameter.ModelBound[0][1] = YMin;
//     SimulationParameter.ModelBound[0][2] = ZMin;
//     SimulationParameter.ModelBound[1][0] = XMax;
//     SimulationParameter.ModelBound[1][1] = YMax;
//     SimulationParameter.ModelBound[1][2] = ZMax;
//     for (gxcId i = 0; i < 3; i++)
//     {   
//         SimulationParameter.GridSize[i] = SimulationParameter.MeshfreeInfo.WinMax[i] * 1.5;
//         SimulationParameter.GridNum[i] = floor((SimulationParameter.ModelBound[1][i] - SimulationParameter.ModelBound[0][i]) / (SimulationParameter.GridSize[i])) + 1;
//     }
//     SimulationParameter.CellInGrid.resize(SimulationParameter.GridNum[0]);
//     for (gxcId i = 0; i < SimulationParameter.GridNum[0]; i++)
//     {
//         SimulationParameter.CellInGrid[i].resize(SimulationParameter.GridNum[1]);
//         for (gxcId j = 0; j < SimulationParameter.GridNum[1]; j++)
//         {
//             SimulationParameter.CellInGrid[i][j].resize(SimulationParameter.GridNum[2]);
//         }
//     }
// }
