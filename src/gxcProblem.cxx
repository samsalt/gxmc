#include "gxcProblem.h"

gxcProblem::gxcProblem()
{
    ReadModel();
    ReadControl();
    StateInitiate();
    EchoVar("Using number of threads", SimulationParameter.ThreadNumber);
    EchoStr("Problem initiated successfully");
}
void gxcProblem::TimeIntegration()
{
    EchoStr("Start time integrating");
    chrono::steady_clock sc;
    auto start = sc.now();
    while (SimulationParameter.TimeCurrent < SimulationParameter.TimeEnd)
    {

        for (unsigned i = 0; i < LoadingFaceSet.size(); i++)
            LoadingFaceSet[i].UpdateCurrentBoundaryCondition(SimulationParameter);

        for (unsigned i = 0; i < EssentialNodeSet.size(); i++)
            EssentialNodeSet[i].UpdateCurrentBoundaryCondition(SimulationParameter);

#pragma omp parallel for
        for (unsigned i = 0; i < Node.size(); i++)
            Node[i].Predictor(SimulationParameter);

        for (unsigned i = 0; i < EssentialNodeSet.size(); i++)
            EssentialNodeSet[i].Predictor(Node, SimulationParameter.Dlt);

        for (unsigned i = 0; i < ContactNodeSet.size(); i++)
            ContactNodeSet[i].UpdatePosition(Node);

        if (SimulationParameter.LagrangianType == 0)
        {
            for (unsigned BlockId = 0; BlockId < Block.size(); BlockId++)
            {
#pragma omp parallel for
                for (unsigned i = 0; i < Block[BlockId].Cell.size(); i++)
                    Block[BlockId].Cell[i]->UpdateFintTot(Node, Block[BlockId].Info, FintTemp);
            }
        }
        else
        {
            for (unsigned BlockId = 0; BlockId < Block.size(); BlockId++)
            {
#pragma omp parallel for
                for (unsigned i = 0; i < Block[BlockId].Cell.size(); i++)
                {
                    Block[BlockId].Cell[i]->UpdateFintUp(Node, Block[BlockId].Info, SimulationParameter.Dlt, FintTemp);
                }
            }
        }

        for (unsigned i = 0; i < LoadingFaceSet.size(); i++)
            LoadingFaceSet[i].UpdateFext(FextTemp);

        Assemble();

        for (unsigned i = 0; i < TargetFaceSet.size(); i++)
            TargetFaceSet[i]->CheckImpact(ContactNodeSet);

        for (unsigned i = 0; i < ContactNodeSet.size(); i++)
            ContactNodeSet[i].ForceCorrector(SimulationParameter.Dlt, TargetFaceSet, Node);

        for (unsigned i = 0; i < TargetFaceSet.size(); i++)
            TargetFaceSet[i]->AssignToNode(ContactNodeSet, Node);

#pragma omp parallel for
        for (unsigned i = 0; i < Node.size(); i++)
            Node[i].Corrector(SimulationParameter);

        for (unsigned i = 0; i < ContactNodeSet.size(); i++)
            ContactNodeSet[i].Corrector(SimulationParameter.Dlt, TargetFaceSet, Node);

        for (unsigned i = 0; i < EssentialNodeSet.size(); i++)
            EssentialNodeSet[i].Corrector(Node, SimulationParameter.Dlt);

        UpdateEnergy();

        if (SimulationParameter.TimeCurrent >= SimulationParameter.TimetoOutput)
            Output(Block, Node, SimulationParameter);

        // EchoVarDebug(SimulationParameter.InternalEnergy);

        SimulationParameter.TimeCurrent += SimulationParameter.Dlt;
    }

    auto end = sc.now();
    auto simulationTime = static_cast<chrono::duration<double>>(end - start);
    EchoVar("Simulation duration (sec): ", simulationTime.count());
}
void gxcProblem::UpdateEnergy()
{
    gxcData LocalInternalEnergy{}, LocalKineticEnergy{};
#pragma omp parallel for reduction(+: LocalInternalEnergy)
        for (unsigned i = 0; i < Node.size(); i++)
            LocalInternalEnergy+=Node[i].GetInternalEnergy();

    SimulationParameter.InternalEnergy = LocalInternalEnergy;

#pragma omp parallel for reduction(+: LocalKineticEnergy)
    for (unsigned i = 0; i < Node.size(); i++)
    {
        LocalKineticEnergy+=Node[i].GetKineticEnergy();
    }
    SimulationParameter.KineticEnergy = LocalKineticEnergy;
}