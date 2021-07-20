#pragma once
#include <vector>
#include "gxcType.h"
#include "gxcEcho.h"
#include "gxcControl.h"


#include <iostream>

class gxcNode
{
public:
    gxcData Displacement[3]{};
    gxcData DisplacementIncrement[3]{};
    gxcData Velocity[3]{};
    gxcData Acceleration[3]{};
    gxcData Fint[3]{};
    gxcData FintPre[3]{};
    gxcData Fext[3]{};
    gxcData Mass{};
    gxcData InternalEnergy{};
    gxcData InitCoordinate[3] {};
    bool ImpactFlag {false};
    gxcData CorrectForce[3] {};
    gxcData SlaveMass {};
    // gxcId BoundarySetId{};
    // short BoundaryCondition[3]{};
    void Predictor(const gxcControl &SimulationParameter);
    void Corrector(const gxcControl &SimulationParameter);
    gxcData GetInternalEnergy();
    gxcData GetKineticEnergy();
};