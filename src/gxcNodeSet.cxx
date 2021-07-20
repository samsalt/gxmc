#include "gxcNodeSet.h"
#include "gxcFace.h"

void gxcEssentialNodeSet::UpdateCurrentBoundaryCondition(const gxcControl &SimulationParameter)
{
    //find the current step
    for (unsigned j = 0; j < EssentialBoundaryConditionSet.size(); j++)
    {
        if (SimulationParameter.TimeCurrent <= EssentialBoundaryConditionSet[j][0])
        {
            gxcData LastStepVal{}, LastTime{};

            for (gxcId k = 0; k < 3; k++)
            {
                if (j != 0)
                {
                    LastStepVal = EssentialBoundaryConditionSet[j - 1][k + 1];
                    LastTime = EssentialBoundaryConditionSet[j - 1][0];
                }
                CurrentEssentialBoundaryCondition[k] = LinearInterpolation(SimulationParameter.TimeCurrent, LastTime, EssentialBoundaryConditionSet[j][0], LastStepVal, EssentialBoundaryConditionSet[j][k + 1]);
                CurrentVelocity[k] = (EssentialBoundaryConditionSet[j][k + 1] - LastStepVal) / (EssentialBoundaryConditionSet[j][0] - LastTime);
            }
        }
    }
    // EchoVarDebug(CurrentEssentialBoundaryCondition[0]);
}

void gxcEssentialNodeSet::Predictor(std::vector<gxcNode> &Node, const double dlt)
{
#pragma omp parallel for
    for (unsigned i = 0; i < NodeId.size(); i++)
    {
        for (gxcId j = 0; j < 3; j++)
        {
            if (EBCSwitch[j])
            {
                Node[NodeId[i]].Displacement[j] = CurrentEssentialBoundaryCondition[j];
                Node[NodeId[i]].Velocity[j] = CurrentVelocity[j];
                Node[NodeId[i]].DisplacementIncrement[j] = CurrentVelocity[j] * dlt;
            }
        }
    }
}

void gxcEssentialNodeSet::Corrector(std::vector<gxcNode> &Node, const double dlt)
{
#pragma omp parallel for
    for (unsigned i = 0; i < NodeId.size(); i++)
    {
        for (gxcId j = 0; j < 3; j++)
        {
            if (EBCSwitch[j])
            {
                Node[NodeId[i]].Displacement[j] = CurrentEssentialBoundaryCondition[j];
                Node[NodeId[i]].Velocity[j] = CurrentVelocity[j];
                Node[NodeId[i]].DisplacementIncrement[j] = CurrentVelocity[j] * dlt;
            }
        }
    }
}

void gxcContactNodeSet::UpdatePosition(   const std::vector<gxcNode> &Node)
{
#pragma omp parallel for
    for (unsigned i = 0; i < ContactNode.size(); i++)
    {
        ContactNode[i].Position[0] = Node[ContactNode[i].NodeId].InitCoordinate[0] + Node[ContactNode[i].NodeId].Displacement[0];
        ContactNode[i].Position[1] = Node[ContactNode[i].NodeId].InitCoordinate[1] + Node[ContactNode[i].NodeId].Displacement[1];
        ContactNode[i].Position[2] = Node[ContactNode[i].NodeId].InitCoordinate[2] + Node[ContactNode[i].NodeId].Displacement[2];
    }
}

void gxcContactNodeSet::ForceCorrector(const gxcData Dlt, std::vector<gxcTargetFaceSet*> &TargetFaceSet,  std::vector<gxcNode> &Node)
{
    gxcData Dlts;
    Dlts = Dlt * Dlt;
    // #pragma omp parallel for
    for (unsigned i = 0; i < ContactNode.size(); i++)
    {
        if (ContactNode[i].ImpactFlag)
        {
            ContactNode[i].CorrectForce = 2 * Node[ContactNode[i].NodeId].Mass*ContactNode[i].Depth / Dlts;
            TargetFaceSet[ContactNode[i].ContactFaceId[0]]->Face[ContactNode[i].ContactFaceId[1]].ImpactFlag=true;
            TargetFaceSet[ContactNode[i].ContactFaceId[0]]->Face[ContactNode[i].ContactFaceId[1]].CorrectForce+=ContactNode[i].CorrectForce;
            TargetFaceSet[ContactNode[i].ContactFaceId[0]]->Face[ContactNode[i].ContactFaceId[1]].SlaveMass+=Node[ContactNode[i].NodeId].Mass;
            
        }
    }
}

void gxcContactNodeSet::Corrector(const gxcData Dlt, std::vector<gxcTargetFaceSet*> &TargetFaceSet,  std::vector<gxcNode> &Node)
{
#pragma omp parallel for
    for (unsigned i = 0; i < ContactNode.size(); i++)
    {
        if (ContactNode[i].ImpactFlag)
        {
            gxcData CorrectForceVec[3] {};
            CorrectForceVec[0]=-ContactNode[i].CorrectForce*TargetFaceSet[ContactNode[i].ContactFaceId[0]]->Face[ContactNode[i].ContactFaceId[1]].Normal[0];
            CorrectForceVec[1]=-ContactNode[i].CorrectForce*TargetFaceSet[ContactNode[i].ContactFaceId[0]]->Face[ContactNode[i].ContactFaceId[1]].Normal[1];
            CorrectForceVec[2]=-ContactNode[i].CorrectForce*TargetFaceSet[ContactNode[i].ContactFaceId[0]]->Face[ContactNode[i].ContactFaceId[1]].Normal[2];
            gxcData MasterAcceleration[3] {};
            gxcId VertexNumLocal=TargetFaceSet[ContactNode[i].ContactFaceId[0]]->VertexNum;
            for (int j = 0; j < VertexNumLocal; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    MasterAcceleration[k]+=Node[TargetFaceSet[ContactNode[i].ContactFaceId[0]]->Face[ContactNode[i].ContactFaceId[1]].NodeList[j]].Acceleration[k]/VertexNumLocal;
                }
            }
            for (int k = 0; k < 3; k++)
            {
                Node[ContactNode[i].NodeId].Velocity[k] -= 0.5 * Dlt * Node[ContactNode[i].NodeId].Acceleration[k];
                Node[ContactNode[i].NodeId].Acceleration[k]=MasterAcceleration[k]-CorrectForceVec[k]/Node[ContactNode[i].NodeId].Mass;
                Node[ContactNode[i].NodeId].Velocity[k] += 0.5 * Dlt * Node[ContactNode[i].NodeId].Acceleration[k];
            }
            ContactNode[i].ImpactFlag=false;
        }
    }
}