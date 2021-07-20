#include "gxcProblem.h"

// void UpdateWin(std::vector<gxcNode> Node, std::vector<gxcBlock> &Block, gxcControl &SimulationParameter)
// {
//     for (gxcId BlockId = 0; BlockId < SimulationParameter.BlockNum; BlockId++)
//         for (gxcId i = 0; i < Block[BlockId].Info.CellNum; i++)
//         {
//             gxcData CellDiameter[3]{};
//             gxcId BaseId, ComId;
//             gxcData DiamaterTemp{};
//             for (gxcId j = 0; j < Block[BlockId].Info.VertexPerCellNum; j++)
//             {
//                 BaseId = Block[BlockId].Cell[i]->GetVertexID(j);
//                 for (gxcId k = j + 1; k < Block[BlockId].Info.VertexPerCellNum; k++)
//                 {
//                     ComId = Block[BlockId].Cell[i]->GetVertexID(k);
//                     for (gxcId kk = 0; kk < 3; kk++)
//                     {
//                         DiamaterTemp = std::abs(Node[BaseId].InitCoordinate[kk] - Node[ComId].InitCoordinate[kk]);
//                         if (DiamaterTemp > CellDiameter[kk])
//                             CellDiameter[kk] = DiamaterTemp;
//                     }
//                 }
//             }
//             gxcId GlobalIdLocal;
//             GlobalIdLocal = Block[BlockId].Cell[i]->GetGlobalId();
//             for (gxcId j = 0; j < 3; j++)
//             {
//                 SimulationParameter.MeshfreeInfo.Win[GlobalIdLocal][j] = CellDiameter[j] * SimulationParameter.MeshfreeInfo.NormalWin;
//                 if (SimulationParameter.MeshfreeInfo.Win[GlobalIdLocal][j] > SimulationParameter.MeshfreeInfo.WinMax[j])
//                     SimulationParameter.MeshfreeInfo.WinMax[j] = SimulationParameter.MeshfreeInfo.Win[GlobalIdLocal][j];
//             }
//         }
// }

// void UpdateSingularKernel()
// {

// }

// void NodeSetSwap(const std::vector<gxcBlock> &Block)
// {
//     std::vector<std::vector<gxcId>> CellContainingNode(SimulationParameter.np);
//     for (gxcId BlockId = 0; BlockId < SimulationParameter.BlockNum; BlockId++)
//         for (gxcId i = 0; i < Block[BlockId].Info.CellNum; i++)
//             for (gxcId j = 0; j < Block[BlockId].Info.VertexPerCellNum; j++)
//                 CellContainingNode[Block[BlockId].Cell[i]->GetVertexID(j)].push_back(Block[BlockId].Cell[i]->GetGlobalId());

//     std::vector<std::vector<gxcId>> NodeSetTemp(Model.EssentialNodeSet);
//     for (gxcId i = 0; i < SimulationParameter.NodeSetNum; i++)
//     {
//         gxcId NodeInCellLocal{};
//         std::vector<bool> DuplicationCheck(SimulationParameter.nc, true);
//         for (std::vector<gxcId>::iterator NodeIt = NodeSetTemp[i].begin(); NodeIt != NodeSetTemp[i].end(); ++NodeIt)
//         {
//             for (std::vector<gxcId>::iterator CellIt = CellContainingNode[*NodeIt].begin(); CellIt != CellContainingNode[*NodeIt].end(); ++CellIt)
//             {
//                 if (DuplicationCheck[*CellIt])
//                 {
//                     DuplicationCheck[*CellIt] = false;
//                     Model.EssentialNodeSet[i][NodeInCellLocal] = *CellIt;
//                     NodeInCellLocal++;
//                 }
//             }
//         }
//         Model.EssentialNodeSet[i].erase(Model.EssentialNodeSet[i].begin() + NodeInCellLocal, Model.EssentialNodeSet[i].end());
//     }
// }

void gxcProblem::StateInitiate()
{

    // update  volume
    #pragma omp parallel for
    for (unsigned BlockId = 0; BlockId < Block.size(); BlockId++)
        for (unsigned i = 0; i < Block[BlockId].Cell.size(); i++)
            Block[BlockId].Cell[i]->UpdateVolume(Node);

    // update surface
    #pragma omp parallel for
    for (unsigned i = 0; i < LoadingFaceSet.size(); i++)
        LoadingFaceSet[i].UpdateSurface(Node, Block);

    // update surface
    #pragma omp parallel for
    for (unsigned i = 0; i < TargetFaceSet.size(); i++)
        TargetFaceSet[i]->UpdateSurface(Node, SimulationParameter, Block);

    if (SimulationParameter.MeshfreeSwitch)
    {
        // Node.resize(SimulationParameter.nc);
        // SimulationParameter.np = SimulationParameter.nc;
        // SimulationParameter.MeshfreeInfo.StateInitiate(SimulationParameter.nc);
        // UpdateWin(Node, Block, SimulationParameter);
        // UpdateGrid(Node, Block, SimulationParameter);
        // SimulationParameter.CellPosition.resize(SimulationParameter.nc, std::vector<gxcData>(3, 0));
        // // NodeSetSwap(Model, Block);
        // for (gxcId i = 0; i < SimulationParameter.nc; i++)
        // {
        //     gxcId BlockId = SimulationParameter.BlockCellId[i][0];
        //     gxcId CellId = SimulationParameter.BlockCellId[i][1];
        //     for (gxcId j = 0; j < 3; j++)
        //         SimulationParameter.CellPosition[i][j] = Block[BlockId].Cell[CellId]->GetPosition(j);
        // }
        // // find grid index
        // for (unsigned BlockId = 0; BlockId < Block.size(); BlockId++)
        //     for (unsigned i = 0; i < Block[BlockId].Cell.size(); i++)
        //         Block[BlockId].Cell[i]->FindGrid(SimulationParameter);

        // // find neighbor
        // for (unsigned BlockId = 0; BlockId < Block.size(); BlockId++)
        //     for (unsigned i = 0; i < Block[BlockId].Cell.size(); i++)
        //         Block[BlockId].Cell[i]->FindNeighbor(SimulationParameter);

        // MeshfreeOutputInit(SimulationParameter);
    }
    else
    {
        Node.resize(SimulationParameter.np);
        SimulationParameter.np = SimulationParameter.np;
    }

    // update essential node set

    // for (gxcId i = 0; i < EssentialNodeSet.size(); i++)
    // {
    //     EssentialNodeSet[i].NodeId=NodeSetTemp[EssentialNodeSet[i].NodeSetId];
    // }

    // update block initial velocity and c mat
    for (unsigned BlockId = 0; BlockId < Block.size(); BlockId++)
    {
        if ((Block[BlockId].Info.VelocityInitial[0] != 0) || (Block[BlockId].Info.VelocityInitial[1] != 0) || (Block[BlockId].Info.VelocityInitial[2] != 0))
        {
            if (SimulationParameter.MeshfreeSwitch)
            {
                for (unsigned i = 0; i < Block[BlockId].Cell.size(); i++)
                {
                    gxcId CurrentNode = Block[BlockId].Cell[i]->GetGlobalId();
                    Node[CurrentNode].Velocity[0] = Block[BlockId].Info.VelocityInitial[0];
                    Node[CurrentNode].Velocity[1] = Block[BlockId].Info.VelocityInitial[1];
                    Node[CurrentNode].Velocity[2] = Block[BlockId].Info.VelocityInitial[2];
                }
            }
            else
            {
                for (unsigned i = 0; i < Block[BlockId].Cell.size(); i++)
                {
                    for (int j = 0; j < Block[BlockId].Info.VertexPerCellNum; j++)
                    {
                        gxcId CurrentNode = Block[BlockId].Cell[i]->GetVertexID(j);
                        Node[CurrentNode].Velocity[0] = Block[BlockId].Info.VelocityInitial[0];
                        Node[CurrentNode].Velocity[1] = Block[BlockId].Info.VelocityInitial[1];
                        Node[CurrentNode].Velocity[2] = Block[BlockId].Info.VelocityInitial[2];
                    }
                }
            }
        }
        Block[BlockId].Info.UpdateBasicMaterialInfo();
        if (Block[BlockId].Info.MaterialType != 2)
            Block[BlockId].Info.FormCmat();
    }

    // update shape function
    for (unsigned BlockId = 0; BlockId < Block.size(); BlockId++)
    {
        #pragma omp parallel for
        for (auto i = 0; i < Block[BlockId].Info.CellNum; i++)
        {
            Block[BlockId].Cell[i]->UpdateTotalShape(Node, SimulationParameter);
        }
    }

    // openMP
    omp_set_num_threads(SimulationParameter.ThreadNumber);
    FextTemp.resize(SimulationParameter.np, std::vector<std::vector<gxcData>>(3, std::vector<gxcData>(SimulationParameter.ThreadNumber, 0)));
    FintTemp.resize(SimulationParameter.np, std::vector<std::vector<gxcData>>(3, std::vector<gxcData>(SimulationParameter.ThreadNumber, 0)));

    // update mass matrix

    for (unsigned BlockId = 0; BlockId < Block.size(); BlockId++)
        for (unsigned i = 0; i < Block[BlockId].Cell.size(); i++)
            Block[BlockId].Cell[i]->UpdateLumpedMass(Node, Block[BlockId].Info);

    for (unsigned i = 0; i < ContactNodeSet.size(); i++)
        ContactNodeSet[i].UpdatePosition(Node);

    for (unsigned i = 0; i < TargetFaceSet.size(); i++)
        TargetFaceSet[i]->InitImpactCondition(ContactNodeSet);
}

void gxcProblem::Assemble()
{
    #pragma omp parallel for
    for (unsigned i = 0; i < Node.size(); i++)
    {
        for (auto j = 0; j < 3; j++)
        {
            for (auto k = 0; k < SimulationParameter.ThreadNumber; k++)
            {
                Node[i].Fint[j] += FintTemp[i][j][k];
                FintTemp[i][j][k] = 0;
                Node[i].Fext[j] += FextTemp[i][j][k];
                FextTemp[i][j][k] = 0;
            }
        }
    }
}
