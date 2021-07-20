#include "gxcProblem.h"

void gxcProblem:: ReadControl()
{
    std::ifstream ControlFile("control.dat");
    std::string line;
    int EssentialNodeSetNum {};

    if (ControlFile)
    {
        while (ControlFile)
        {
            getline(ControlFile, line);
            if (line == "# time duration")
            {
                ControlFile >> SimulationParameter.TimeEnd;
            }
            else if (line == "# time step size")
            {
                ControlFile >> SimulationParameter.Dlt;
            }
            else if (line == "# time output period")
            {
                ControlFile >> SimulationParameter.TimeOutputPeriod;
            }
            else if (line == "# LagrangianType")
            {
                ControlFile >> SimulationParameter.LagrangianType;
            }
            else if (line == "# material property")
            {
                gxcId BlockId;
                ControlFile >> BlockId;
                BlockId--;
                ControlFile >> Block[BlockId].Info.MaterialType;
                for (int i = 0; i < 20; i++)
                    ControlFile >> Block[BlockId].Info.MaterialProperty[i];
            }
            else if (line == "# block set initial velocity")
            {
                gxcId BlockId;
                ControlFile >> BlockId;
                BlockId--;
                ControlFile >> Block[BlockId].Info.VelocityInitial[0];
                ControlFile >> Block[BlockId].Info.VelocityInitial[1];
                ControlFile >> Block[BlockId].Info.VelocityInitial[2];
            }
            else if (line == "# essential boundary condition")
            {
                gxcId NodeSetId;
                ControlFile >> NodeSetId;
                gxcEssentialNodeSet ENTemp;
                EssentialNodeSet.push_back(ENTemp);
                EssentialNodeSet[EssentialNodeSetNum].NodeSetId=NodeSetId-1;
                EssentialNodeSet[EssentialNodeSetNum].NodeId=NodeSetTemp[NodeSetId-1];
                std::string DisplacementName = "Displacement" + std::to_string(NodeSetId) + ".dat";
                std::ifstream DisplacementFile(DisplacementName);
                std::string lineDisplacement;
                int DisplacementStep{};
                NodeSetId--;
                EssentialNodeSet[EssentialNodeSetNum].CurrentEssentialBoundaryCondition.resize(3, 0);
                EssentialNodeSet[EssentialNodeSetNum].CurrentVelocity.resize(3, 0);
                EssentialNodeSet[EssentialNodeSetNum].EBCSwitch.resize(3, false);
                int swtemp{};
                ControlFile >> swtemp;
                if (swtemp)
                    EssentialNodeSet[EssentialNodeSetNum].EBCSwitch[0] = true;
                ControlFile >> swtemp;
                if (swtemp)
                    EssentialNodeSet[EssentialNodeSetNum].EBCSwitch[1] = true;
                ControlFile >> swtemp;
                if (swtemp)
                    EssentialNodeSet[EssentialNodeSetNum].EBCSwitch[2] = true;

                if (DisplacementFile)
                {
                    std::getline(DisplacementFile, lineDisplacement);
                    DisplacementFile >> DisplacementStep;
                    EssentialNodeSet[EssentialNodeSetNum].EssentialBoundaryConditionSet.resize(DisplacementStep, std::vector<gxcData>(4, 0));
                    std::getline(DisplacementFile, lineDisplacement);
                    std::getline(DisplacementFile, lineDisplacement);
                    for (int i = 0; i < DisplacementStep; i++)
                    {
                        DisplacementFile >> EssentialNodeSet[EssentialNodeSetNum].EssentialBoundaryConditionSet[i][0];
                        DisplacementFile >> EssentialNodeSet[EssentialNodeSetNum].EssentialBoundaryConditionSet[i][1];
                        DisplacementFile >> EssentialNodeSet[EssentialNodeSetNum].EssentialBoundaryConditionSet[i][2];
                        DisplacementFile >> EssentialNodeSet[EssentialNodeSetNum].EssentialBoundaryConditionSet[i][3];
                    }
                    if (EssentialNodeSet[EssentialNodeSetNum].EssentialBoundaryConditionSet[DisplacementStep-1][0]<SimulationParameter.TimeEnd) EchoError("The time duration in "+DisplacementName+" cannot be smaller than the total simulation time");
                    DisplacementFile.close();
                    EssentialNodeSetNum++;
                }
                else
                {
                    EchoError("Cannot find file" + DisplacementName+" for Essential boundary condition");
                }
            }
            else if (line == "# contact nodeset")
            {
                gxcId NodeSetId;
                std::vector<gxcId> contactSet;
                getline(ControlFile, line);
                std::stringstream stream(line);
                while (1)
                {
                    stream >> NodeSetId;
                    if (!stream)
                        break;
                    contactSet.push_back(NodeSetId-1);
                }
                ContactNodeSet.resize(contactSet.size());
                for (unsigned j = 0; j < contactSet.size(); j++)
                {
                    ContactNodeSet[j].NodeSetId=contactSet[j];
                    ContactNodeSet[j].ContactNode.resize(NodeSetTemp[contactSet[j]].size());
                    for (unsigned k = 0; k < ContactNodeSet[j].ContactNode.size(); k++)
                    {
                        ContactNodeSet[j].ContactNode[k].NodeId=NodeSetTemp[contactSet[j]][k];
                        ContactNodeSet[j].ContactNode[k].Position[0]=Node[ContactNodeSet[j].ContactNode[k].NodeId].InitCoordinate[0];
                        ContactNodeSet[j].ContactNode[k].Position[1]=Node[ContactNodeSet[j].ContactNode[k].NodeId].InitCoordinate[1];
                        ContactNodeSet[j].ContactNode[k].Position[2]=Node[ContactNodeSet[j].ContactNode[k].NodeId].InitCoordinate[2];
                    }
                }
            }
            else if (line == "# target faceset")
            {
                gxcId SidesetId;
                std::vector<gxcId> targetSet;
                getline(ControlFile, line);
                std::stringstream stream(line);
                while (1)
                {
                    stream >> SidesetId;
                    if (!stream)
                        break;
                    targetSet.push_back(SidesetId-1);
                }
                TargetFaceSet.resize(targetSet.size());
                for (unsigned j = 0; j < targetSet.size(); j++)
                {   
                    gxcId FaceSetIdTemp {};
                    FaceSetIdTemp=targetSet[j];
                    gxcId BlockId1 = SimulationParameter.BlockCellId[FaceTemp[FaceSetIdTemp][0][1]][0];
                    switch (Block[BlockId1].Info.VertexPerCellNum)
                    {
                    case 4:
                        TargetFaceSet[j]=new gxcTargetFaceSetTet;
                        break;
                    case 8:
                        TargetFaceSet[j]=new gxcTargetFaceSetHex;
                        break;                    
                    default:
                        EchoError("Unknown element for face defination");
                    }
                    TargetFaceSet[j]->FaceSetId=j;
                    TargetFaceSet[j]->Face.resize(FaceTemp[FaceSetIdTemp].size());
                    for (unsigned k = 0; k < FaceTemp[FaceSetIdTemp].size(); k++)
                    {
                        TargetFaceSet[j]->Face[k].FaceId=FaceTemp[FaceSetIdTemp][k][0];
                        TargetFaceSet[j]->Face[k].CellId=FaceTemp[FaceSetIdTemp][k][1];
                    }
                }
            }
            else if (line == "# natural boundary condition")
            {
                gxcId SidesetId;
                std::vector<gxcId> LoadingSet;
                getline(ControlFile, line);
                std::stringstream stream(line);
                while (1)
                {
                    stream >> SidesetId;
                    if (!stream)
                        break;
                    LoadingSet.push_back(SidesetId-1);
                }
                LoadingFaceSet.resize(LoadingSet.size());
                for (unsigned j = 0; j < LoadingSet.size(); j++)
                {
                    gxcId BlockId, CellId, SetId;
                    SetId=LoadingSet[j];
                    BlockId = SimulationParameter.BlockCellId[FaceTemp[SetId][0][1]][0];
                    CellId = SimulationParameter.BlockCellId[FaceTemp[SetId][0][1]][1];
                    LoadingFaceSet[j].VertexNum = Block[BlockId].Cell[CellId]->SideNum;
                    if (LoadingFaceSet[j].VertexNum == 6)
                    {
                        LoadingFaceSet[j].VertexNum = 4;
                    }
                    else if (LoadingFaceSet[j].VertexNum == 4)
                    {
                        LoadingFaceSet[j].VertexNum = 3;
                    }
                    gxcId FaceNodeId[6][4] = {{0, 4, 5, 1}, {1, 5, 6, 2}, {2, 6, 7, 3}, {3, 7, 4, 0}, {0, 1, 2, 3}, {4, 7, 6, 5}};
                    LoadingFaceSet[SetId].Face.resize(FaceTemp[SetId].size());
                    for (unsigned jj = 0; jj < FaceTemp[SetId].size(); jj++)
                    {
                        BlockId = SimulationParameter.BlockCellId[FaceTemp[SetId][jj][1]][0];
                        CellId = SimulationParameter.BlockCellId[FaceTemp[SetId][jj][1]][1];

                        LoadingFaceSet[SetId].Face[jj].FaceId = FaceTemp[SetId][jj][0];
                        for (unsigned k = 0; k < 4; k++)
                            LoadingFaceSet[SetId].Face[jj].NodeList[k] = Block[BlockId].Cell[CellId]->VertexId[FaceNodeId[FaceTemp[SetId][jj][0]][k]];
                    }
                    SidesetId = LoadingSet[j];
                    std::string LoadingName = "Loading" + std::to_string(SidesetId) + ".dat";
                    std::ifstream LoadingFile(LoadingName);
                    std::string lineLoading;
                    int LoadingStep{};
                    SidesetId--;
                    if (LoadingFile)
                    {
                        std::getline(LoadingFile, lineLoading);
                        std::cout << lineLoading;
                        LoadingFile >> LoadingStep;
                        LoadingFaceSet[SidesetId].NaturalBoundaryConditionSet.resize(LoadingStep, std::vector<gxcData>(4, 0));
                        LoadingFaceSet[SidesetId].CurrentNaturalBoundaryConditionSet.resize(3, 0);
                        std::getline(LoadingFile, lineLoading);
                        std::getline(LoadingFile, lineLoading);
                        std::cout << lineLoading;
                        for (int i = 0; i < LoadingStep; i++)
                        {
                            LoadingFile >> LoadingFaceSet[SidesetId].NaturalBoundaryConditionSet[i][0];
                            LoadingFile >> LoadingFaceSet[SidesetId].NaturalBoundaryConditionSet[i][1];
                            LoadingFile >> LoadingFaceSet[SidesetId].NaturalBoundaryConditionSet[i][2];
                            LoadingFile >> LoadingFaceSet[SidesetId].NaturalBoundaryConditionSet[i][3];
                        }
                        if (LoadingFaceSet[SidesetId].NaturalBoundaryConditionSet[LoadingStep-1][0]<SimulationParameter.TimeEnd) EchoError("Time duration in "+LoadingName+"cannot be smallter than the total simulation time");
                        LoadingFile.close();
                    }
                    else
                    {
                        EchoError("Cannot find file" + LoadingName+ " for natural boundary condition");
                    }
                }
            }
            else if (line == "# thread number")
            {
                ControlFile >> SimulationParameter.ThreadNumber;
            }
            else if (line == "# rk information")
            {
                // SimulationParameter.MeshfreeSwitch = true;
                // ControlFile >> SimulationParameter.MeshfreeInfo.Deg;
                // ControlFile >> SimulationParameter.MeshfreeInfo.NormalWin;
            }
        }
        // EchoStr("Finish reading control file");
    }
    else
        EchoError("Cannot find control.dat file!");

    ControlFile.close();
    CheckControl();

}

void gxcProblem::CheckControl()
{
    EchoVar("Number of essential nodeset",EssentialNodeSet.size());
    EchoVar("Number of contact nodeset",ContactNodeSet.size());
    EchoVar("Number of loading faceset",LoadingFaceSet.size());
    EchoVar("Number of target faceset",TargetFaceSet.size());
    if  ((SimulationParameter.NodeSetNum-ContactNodeSet.size()-EssentialNodeSet.size())>0) EchoVar("Redundant nodeset number",SimulationParameter.NodeSetNum-EssentialNodeSet.size());
    if  ((SimulationParameter.SideSetNum-TargetFaceSet.size()-LoadingFaceSet.size())>0) EchoVar("Redundant sideset number",SimulationParameter.NodeSetNum-EssentialNodeSet.size()); 
}