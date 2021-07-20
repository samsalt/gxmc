#include "gxcNode.h"

void gxcNode::Predictor(const gxcControl &SimulationParameter)
{

    gxcData Dlt = SimulationParameter.Dlt;

    for (int j = 0; j < 3; j++)
    {

        DisplacementIncrement[j] = Dlt * Velocity[j] + Dlt * Dlt * 0.5 * Acceleration[j];
        Displacement[j] += DisplacementIncrement[j];
        Velocity[j] += 0.5 * Dlt * Acceleration[j];
    }
}

void gxcNode::Corrector(const gxcControl &SimulationParameter)
{
    gxcData Dlt = SimulationParameter.Dlt;
    
    if (!ImpactFlag)
    for (int j = 0; j < 3; j++)
    {
        Acceleration[j] = (Fext[j] - Fint[j]) / Mass;
        Velocity[j] += 0.5 * Dlt * Acceleration[j];
    }
    else
    {
        for (int j = 0; j < 3; j++)
        {
            Acceleration[j]=CorrectForce[j]/(Mass+SlaveMass);
            Velocity[j] += 0.5 * Dlt * Acceleration[j];
            CorrectForce[j]=0;
        }
        SlaveMass=0;
        ImpactFlag=false;
    }
}
gxcData gxcNode::GetInternalEnergy()
{
    // InternalEnergy=0;
    for (int j = 0; j < 3; j++)
    {
        gxcData InternalEnergyInc {};
        InternalEnergyInc=0.5*DisplacementIncrement[j]*(Fint[j]+FintPre[j]);
        InternalEnergy+=InternalEnergyInc;
        FintPre[j]=Fint[j];
        Fint[j] = 0;
        Fext[j] = 0;
    }
        return InternalEnergy;
}
gxcData gxcNode::GetKineticEnergy()
{
    return(0.5*Mass*(Velocity[0]*Velocity[0]+Velocity[1]*Velocity[1]+Velocity[2]*Velocity[2]));
}