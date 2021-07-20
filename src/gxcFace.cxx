#include "gxcFace.h"

void gxcLoadingFaceSet::UpdateSurface(std::vector<gxcNode> Node, std::vector<gxcBlock> &Block)
{
    for (unsigned i = 0; i < Face.size(); i++)
    {
        gxcData xyzFace[4][3];
        for (int j = 0; j < 4; j++)
            for (int k = 0; k < 3; k++)
                xyzFace[j][k] = Node[Face[i].NodeList[j]].InitCoordinate[k];

        gxcData aArray[3], bArray[3], cross[3];
        for (int j = 0; j < 3; j++)
        {
            aArray[j] = xyzFace[3][j] - xyzFace[1][j];
            bArray[j] = xyzFace[2][j] - xyzFace[0][j];
        }
        cross[0] = aArray[1] * bArray[2] - aArray[2] * bArray[1];
        cross[1] = aArray[2] * bArray[0] - aArray[0] * bArray[2];
        cross[2] = aArray[0] * bArray[1] - aArray[1] * bArray[0];
        Face[i].Area = 0.5 * std::sqrt(cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]);
        for (int j = 0; j < 3; j++)
            Face[i].Normal[j] = cross[j] / (2 * Face[i].Area);
    }
}

void gxcLoadingFaceSet::UpdateCurrentBoundaryCondition(const gxcControl &SimulationParameter)
{
    //find the current step
    for (unsigned j = 0; j < NaturalBoundaryConditionSet.size(); j++)
    {
        if (SimulationParameter.TimeCurrent <= NaturalBoundaryConditionSet[j][0])
        {
            gxcData LastStepVal{}, LastTime{};

            for (gxcId k = 0; k < 3; k++)
            {
                if (j != 0)
                {
                    LastStepVal = NaturalBoundaryConditionSet[j - 1][k + 1];
                    LastTime = NaturalBoundaryConditionSet[j - 1][0];
                }
                CurrentNaturalBoundaryConditionSet[k] = LinearInterpolation(SimulationParameter.TimeCurrent, LastTime, NaturalBoundaryConditionSet[j][0], LastStepVal, NaturalBoundaryConditionSet[j][k + 1]);
            }
        }
    }
    // EchoVarDebug(CurrentNaturalBoundaryConditionSet[0]);
}
void gxcLoadingFaceSet::UpdateFext(std::vector<std::vector<std::vector<gxcData>>> &FextTemp)
{
    gxcData FaceShape[6][4]{{0.0625, 0.1875, 0.5625, 0.1875}, {0.1875, 0.5625, 0.1875, 0.0625}, {0.09375, 0.28125, 0.09375, 0.03125}, {0.03125, 0.09375, 0.28125, 0.09375}, {0.09375, 0.28125, 0.09375, 0.03125}, {0.1875, 0.0625, 0.1875, 0.5625}};

#pragma omp parallel for
    for (unsigned i = 0; i < Face.size(); i++)
    {
        gxcId NodeId{};
        gxcData *CurrentShape;
        gxcData LocalArea;
        int threadId{};
        threadId = omp_get_thread_num();
        gxcData Loading[3]{};
        for (int j = 0; j < 3; j++)
            Loading[j] = Face[i].Normal[j] * CurrentNaturalBoundaryConditionSet[j];

        CurrentShape = FaceShape[Face[i].FaceId];
        LocalArea = Face[i].Area;

        for (int j = 0; j < 4; j++)
        {
            NodeId = Face[i].NodeList[j];
            for (gxcId k = 0; k < 3; k++)
            {
                FextTemp[NodeId][k][threadId] += CurrentShape[j] * Loading[k] * LocalArea;
            }
        }
    }
}

void gxcTargetFaceSetHex::UpdateSurface(std::vector<gxcNode> Node, const gxcControl SimulationParameter, std::vector<gxcBlock> &Block)
{
    gxcId FaceNodeId[6][4] = {{0, 4, 5, 1}, {1, 5, 6, 2}, {2, 6, 7, 3}, {3, 7, 4, 0}, {0, 1, 2, 3}, {4, 7, 6, 5}};
#pragma omp parallel for
    for (unsigned i = 0; i < Face.size(); i++)
    {
        gxcId BlockIdTemp, CellIdTemp;
        BlockIdTemp = SimulationParameter.BlockCellId[Face[i].CellId][0];
        CellIdTemp = SimulationParameter.BlockCellId[Face[i].CellId][1];

        for (int k = 0; k < 4; k++)
            Face[i].NodeList[k] = Block[BlockIdTemp].Cell[CellIdTemp]->VertexId[FaceNodeId[Face[i].FaceId][k]];
        gxcData xyzFace[4][3];
        for (int j = 0; j < 4; j++)
            for (int k = 0; k < 3; k++)
            {
                xyzFace[j][k] = Node[Face[i].NodeList[j]].InitCoordinate[k];
                Face[i].Centroid[k] += 0.25 * xyzFace[j][k];
                Face[i].CurrentCentroid[k] = Face[i].Centroid[k];
            }

        gxcData aArray[3], bArray[3], cross[3];
        for (int j = 0; j < 3; j++)
        {
            aArray[j] = xyzFace[3][j] - xyzFace[1][j];
            bArray[j] = xyzFace[2][j] - xyzFace[0][j];
        }
        cross[0] = aArray[1] * bArray[2] - aArray[2] * bArray[1];
        cross[1] = aArray[2] * bArray[0] - aArray[0] * bArray[2];
        cross[2] = aArray[0] * bArray[1] - aArray[1] * bArray[0];
        gxcData Area{};
        Area = 0.5 * std::sqrt(cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]);
        for (int j = 0; j < 3; j++)
            Face[i].Normal[j] = cross[j] / (2 * Area);
    }
}

void gxcTargetFaceSet::InitImpactCondition(std::vector<gxcContactNodeSet> &ContactNodeSet)
{
#pragma omp parallel for
    for (unsigned i = 0; i < Face.size(); i++)
    {
        Face[i].FaceToNodeDirection.resize(ContactNodeSet.size());
        for (unsigned j = 0; j < ContactNodeSet.size(); j++)
        {
            Face[i].FaceToNodeDirection[j].resize(ContactNodeSet[j].ContactNode.size(), false);
            double FaceToNodeVec[3]{}, CrossProduct{};
            for (unsigned k = 0; k < ContactNodeSet[j].ContactNode.size(); k++)
            {
                FaceToNodeVec[0] = ContactNodeSet[j].ContactNode[k].Position[0] - Face[i].CurrentCentroid[0];
                FaceToNodeVec[1] = ContactNodeSet[j].ContactNode[k].Position[1] - Face[i].CurrentCentroid[1];
                FaceToNodeVec[2] = ContactNodeSet[j].ContactNode[k].Position[2] - Face[i].CurrentCentroid[2];
                CrossProduct = FaceToNodeVec[0] * Face[i].Normal[0] + FaceToNodeVec[1] * Face[i].Normal[1] + FaceToNodeVec[2] * Face[i].Normal[2];
                if (CrossProduct >= 0)
                    Face[i].FaceToNodeDirection[j][k] = true;
            }
        }
    }
}

void gxcTargetFaceSet::CheckImpact(std::vector<gxcContactNodeSet> &ContactNodeSet)
{
#pragma omp parallel for
    for (unsigned i = 0; i < Face.size(); i++)
    {
        for (unsigned j = 0; j < ContactNodeSet.size(); j++)
        {
            double FaceToNodeVec[3]{}, CrossProduct{}, FaceToNodeLength{};
            bool directionTemp{false};
            for (unsigned k = 0; k < ContactNodeSet[j].ContactNode.size(); k++)
            {
                directionTemp = false;
                FaceToNodeVec[0] = ContactNodeSet[j].ContactNode[k].Position[0] - Face[i].CurrentCentroid[0];
                FaceToNodeVec[1] = ContactNodeSet[j].ContactNode[k].Position[1] - Face[i].CurrentCentroid[1];
                FaceToNodeVec[2] = ContactNodeSet[j].ContactNode[k].Position[2] - Face[i].CurrentCentroid[2];
                CrossProduct = FaceToNodeVec[0] * Face[i].Normal[0] + FaceToNodeVec[1] * Face[i].Normal[1] + FaceToNodeVec[2] * Face[i].Normal[2];
                if (CrossProduct >= 0)
                    directionTemp = true;
                if (directionTemp ^ Face[i].FaceToNodeDirection[j][k])
                {
#pragma omp critical
                    {
                        FaceToNodeLength = VecNorm3(FaceToNodeVec);
                        // FaceToNodeLength=std::abs(CrossProduct);

                        if (!ContactNodeSet[j].ContactNode[k].ImpactFlag)
                        {
                            ContactNodeSet[j].ContactNode[k].ImpactFlag = true;
                            ContactNodeSet[j].ContactNode[k].ContactLength = FaceToNodeLength;
                            ContactNodeSet[j].ContactNode[k].Depth = std::abs(CrossProduct);
                            ContactNodeSet[j].ContactNode[k].ContactFaceId[0] = FaceSetId;
                            ContactNodeSet[j].ContactNode[k].ContactFaceId[1] = i;
                        }
                        else if (FaceToNodeLength < ContactNodeSet[j].ContactNode[k].ContactLength)
                        {

                            ContactNodeSet[j].ContactNode[k].ContactLength = FaceToNodeLength;
                            ContactNodeSet[j].ContactNode[k].Depth = std::abs(CrossProduct);
                            ContactNodeSet[j].ContactNode[k].ContactFaceId[0] = FaceSetId;
                            ContactNodeSet[j].ContactNode[k].ContactFaceId[1] = i;
                        }
                    }
                }
            }
        }
    }
}

void gxcTargetFaceSet::AssignToNode(std::vector<gxcContactNodeSet> &ContactNodeSet, std::vector<gxcNode> &Node)
{
    for (unsigned i = 0; i < Face.size(); i++)
    {
        if (Face[i].ImpactFlag)
        {
            gxcData CorrectForceVec[3] {};
            CorrectForceVec[0]=-Face[i].CorrectForce*Face[i].Normal[0];
            CorrectForceVec[1]=-Face[i].CorrectForce*Face[i].Normal[1];
            CorrectForceVec[2]=-Face[i].CorrectForce*Face[i].Normal[2];
            for (auto j = 0; j < VertexNum; j++)
            {
                Node[Face[i].NodeList[j]].ImpactFlag=true;
                Node[Face[i].NodeList[j]].CorrectForce[0]+=CorrectForceVec[0]/VertexNum;
                Node[Face[i].NodeList[j]].CorrectForce[1]+=CorrectForceVec[1]/VertexNum;
                Node[Face[i].NodeList[j]].CorrectForce[2]+=CorrectForceVec[2]/VertexNum;
                Node[Face[i].NodeList[j]].SlaveMass+=Face[i].SlaveMass/VertexNum;
            }
            Face[i].SlaveMass=0;
            Face[i].CorrectForce=0;
            Face[i].ImpactFlag=false;         
        }
    }
}

void gxcTargetFaceSetTet::UpdateSurface(std::vector<gxcNode> Node, const gxcControl SimulationParameter, std::vector<gxcBlock> &Block)
{
    gxcId FaceNodeId[4][3] = {{0, 3, 1}, {1,3,2}, {0,2,3},{0,1,2}};
#pragma omp parallel for
    for (unsigned i = 0; i < Face.size(); i++)
    {
        gxcId BlockIdTemp, CellIdTemp;
        BlockIdTemp = SimulationParameter.BlockCellId[Face[i].CellId][0];
        CellIdTemp = SimulationParameter.BlockCellId[Face[i].CellId][1];

        for (int k = 0; k < 3; k++)
            Face[i].NodeList[k] = Block[BlockIdTemp].Cell[CellIdTemp]->VertexId[FaceNodeId[Face[i].FaceId][k]];
        gxcData xyzFace[3][3];
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++)
            {
                xyzFace[j][k] = Node[Face[i].NodeList[j]].InitCoordinate[k];
                Face[i].Centroid[k] += xyzFace[j][k]/3;
                Face[i].CurrentCentroid[k] = Face[i].Centroid[k];
            }

        gxcData aArray[3], bArray[3], cross[3];
        for (int j = 0; j < 3; j++)
        {
            aArray[j] = xyzFace[1][j] - xyzFace[0][j];
            bArray[j] = xyzFace[2][j] - xyzFace[0][j];
        }
        cross[0] = aArray[1] * bArray[2] - aArray[2] * bArray[1];
        cross[1] = aArray[2] * bArray[0] - aArray[0] * bArray[2];
        cross[2] = aArray[0] * bArray[1] - aArray[1] * bArray[0];
        gxcData temp{};
        temp = std::sqrt(cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]);
        for (int j = 0; j < 3; j++)
            Face[i].Normal[j] = -cross[j] / temp;
    }
}