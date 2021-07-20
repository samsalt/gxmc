#include <gxcCell.h>

inline void gxcCell::GetCauchyMat(std::array<std::array<gxcData, 3>, 3> &mat, const gxcId QuadraturePointId)
{

    mat[0][0] = Stress[QuadraturePointId][0];
    mat[1][1] = Stress[QuadraturePointId][1];
    mat[2][2] = Stress[QuadraturePointId][2];
    mat[1][2] = Stress[QuadraturePointId][3];
    mat[0][2] = Stress[QuadraturePointId][4];
    mat[0][1] = Stress[QuadraturePointId][5];
    mat[2][1] = mat[1][2];
    mat[2][0] = mat[0][2];
    mat[1][0] = mat[0][1];
}

void gxcCell::SaveStress(const std::array<std::array<gxcData, 3>, 3> &mat, const gxcId QuadraturePointId)
{
    Stress[QuadraturePointId][0] = mat[0][0];
    Stress[QuadraturePointId][1] = mat[1][1];
    Stress[QuadraturePointId][2] = mat[2][2];
    Stress[QuadraturePointId][3] = mat[1][2];
    Stress[QuadraturePointId][4] = mat[0][2];
    Stress[QuadraturePointId][5] = mat[0][1];
}

inline void gxcCell::GetLocalDisplacement(const std::vector<gxcNode> &Node, std::vector<std::vector<gxcData>> &DisplacementLocal)
{
    for (unsigned i = 0; i < NodeNum; i++)
        for (unsigned j = 0; j < 3; j++)
            DisplacementLocal[j][i] = Node[NeighborNode[i]].Displacement[j];
}
inline void gxcCell::GetLocalVelocity(const std::vector<gxcNode> &Node, std::vector<std::vector<gxcData>> &VelocityLocal)
{
    for (unsigned i = 0; i < NodeNum; i++)
        for (unsigned j = 0; j < 3; j++)
            VelocityLocal[j][i] = Node[NeighborNode[i]].Velocity[j];
}
inline void gxcCell::GetLocalTotalShapeDx(const gxcId NodeNum, std::vector<std::vector<gxcData>> &TotalShapeDxLocal, const gxcId QuadraturePointId)
{
    for (int i = 0; i < NodeNum; i++)
        for (int j = 0; j < 3; j++)
            TotalShapeDxLocal[i][j] = TotalShapeDx[QuadraturePointId][i][j];
}

void gxcCell::UpdateLumpedMass(std::vector<gxcNode> &Node, const gxcBlockInfo &BlockInfo)
{
    for (int i = 0; i < NodeNum; i++)
        for (int j = 0; j < NodeNum; j++)
            for (auto QuadraturePointId = 0; QuadraturePointId < QuadratureNum; QuadraturePointId++)
            {
                Node[NeighborNode[i]].Mass += BlockInfo.Density * TotalShape[QuadraturePointId][i] * TotalShape[QuadraturePointId][j] * Jacobian[QuadraturePointId] * QWeight;
            }
}

void gxcCell::SaveStrain(const std::array<std::array<gxcData, 1>, 6> &GreenStrain, const gxcId QuadraturePointId)
{
    Strain[QuadraturePointId][0] = GreenStrain[0][0];
    Strain[QuadraturePointId][1] = GreenStrain[1][1];
    Strain[QuadraturePointId][2] = GreenStrain[2][2];
    Strain[QuadraturePointId][3] = GreenStrain[1][2];
    Strain[QuadraturePointId][4] = GreenStrain[0][2];
    Strain[QuadraturePointId][5] = GreenStrain[0][1];
}

void gxcCell::UpdateFintTot(std::vector<gxcNode> &Node, const gxcBlockInfo &BlockInfo, std::vector<std::vector<std::vector<gxcData>>> &FintTemp)
{

    // zero_matrix<gxcData> zero_m(3, 3);
    // identity_matrix<gxcData> i_m(3);
    // matrix<gxcData> H(3, 3), HTrans(3, 3), DisplacementLocal(3, NodeNum), GreenStrain(3, 3), SPK(3, 3), NominalStress(3, 3), DeformationGradient(3, 3), TotalShapeDxLocal(NodeNum, 3), CauchyStressMat(3, 3), StressVoigt(6, 1), StrainVoigt(6, 1);

    std::array<std::array<gxcData, 3>, 3> H, HTrans, GreenStrain, SPK, NominalStress, DeformationGradient, CauchyStressMat;
    std::array<std::array<gxcData, 1>, 6> StressVoigt, StrainVoigt;
    std::vector<std::vector<gxcData>> DisplacementLocal, TotalShapeDxLocal;
    DisplacementLocal.resize(3, std::vector<gxcData>(NodeNum));
    TotalShapeDxLocal.resize(NodeNum, std::vector<gxcData>(3));

    GetLocalDisplacement(Node, DisplacementLocal);
    for (auto QuadraturePointId = 0; QuadraturePointId < QuadratureNum; QuadraturePointId++)
    {
        GetLocalTotalShapeDx(NodeNum, TotalShapeDxLocal, QuadraturePointId);

        H = prodA33(DisplacementLocal, TotalShapeDxLocal);
        DeformationGradient = H;
        for (int i = 0; i < 3; i++)
            DeformationGradient[i][i] += 1.0;

        HTrans = trans(H);
        GreenStrain = scaleA(0.5, addA(H, addA(HTrans, prodA(HTrans, H))));

        // EchoVar("GreenStrain",GreenStrain);

        Constitution(SPK, GreenStrain, BlockInfo, QuadraturePointId);
        NominalStress = prodA(SPK, trans(DeformationGradient));
        SaveStress(NominalStress, QuadraturePointId);
        // local element stiffness matrix
        // EchoVar("h",H);
        // EchoVar("DisplacementLocal",DisplacementLocal);
        // EchoVar("TotalShapeDxLocal",TotalShapeDxLocal);
        // EchoVar("DeformationGradient",DeformationGradient);

        // assemble to global; this can be moved out
        int threadId{};
        threadId = omp_get_thread_num();

        for (int i = 0; i < NodeNum; i++)
            for (int j = 0; j < 3; j++)
            {
                gxcData FintLocal{};
                for (auto p = 0; p < 3; p++)
                {
                    FintLocal += TotalShapeDxLocal[i][p] * NominalStress[p][j];
                }
                // FintTemp[NeighborNode[i]][j][threadId] += 0;
                FintTemp[NeighborNode[i]][j][threadId] += FintLocal * QWeight * Jacobian[QuadraturePointId];
            }
        // Node[NeighborNode[i]].Fint[j] += FintLocal[i][j] * Volume;
    }
    // EchoVar("cell #",GlobalId);
    // EchoVar("DeformationGradient",DeformationGradient[0][0]);
}

void gxcCell::UpdateFintUp(std::vector<gxcNode> &Node, const gxcBlockInfo &BlockInfo, const gxcData Dlt, std::vector<std::vector<std::vector<gxcData>>> &FintTemp)
{
    //using updated Lagrangian for internal force

    std::array<std::array<gxcData, 3>, 3> H, HTrans, GreenStrain, DeformationGradient;
    std::array<std::array<gxcData, 1>, 6> StrainVoigt;
    std::vector<std::vector<gxcData>> DisplacementLocal, TotalShapeDxLocal;
    DisplacementLocal.resize(3, std::vector<gxcData>(NodeNum));
    TotalShapeDxLocal.resize(NodeNum, std::vector<gxcData>(3));

    GetLocalDisplacement(Node, DisplacementLocal);
    H = {};
    GetLocalTotalShapeDx(NodeNum, TotalShapeDxLocal, 0);

    H = prodA33(DisplacementLocal, TotalShapeDxLocal);
    DeformationGradient = H;
    for (int i = 0; i < 3; i++)
        DeformationGradient[i][i] += 1.0;

    HTrans = trans(H);
    GreenStrain = scaleA(0.5, addA(H, addA(HTrans, prodA(HTrans, H))));
    Mat2VoigtStrain(GreenStrain, StrainVoigt);
    SaveStrain(StrainVoigt, 0);

    std::array<std::array<gxcData, 3>, 3> L, D, W, JaumannRateMat, CauchyRate, CauchyStressMat;
    std::array<std::array<gxcData, 1>, 6> DVoigt, JaumannRateVoigt;
    std::vector<std::vector<gxcData>> VelocityLocal, UpShapeDxLocal;
    VelocityLocal.resize(3, std::vector<gxcData>(NodeNum));
    UpShapeDxLocal.resize(NodeNum, std::vector<gxcData>(3));

    GetLocalVelocity(Node, VelocityLocal);

    for (auto QuadraturePointId = 0; QuadraturePointId < QuadratureNum; QuadraturePointId++)
    {

        UpdateShapeDx(Node, UpShapeDxLocal, QuadraturePointId);

        // Update Velocity Gradiant and the symetric part
        L = prodA33(VelocityLocal, UpShapeDxLocal);

        D = scaleA(0.5, addA(L, trans(L)));
        W = scaleA(0.5, subA(L, trans(L)));

        // get Cauchy stress by using Jaumann rate
        GetCauchyMat(CauchyStressMat, QuadraturePointId);
        Mat2VoigtStrain(D, DVoigt);
        JaumannRateVoigt = prodA(BlockInfo.Cmat, DVoigt);
        Voigt2MatStress(JaumannRateVoigt, JaumannRateMat);
        CauchyRate = addA(prodA(W, CauchyStressMat), prodA(CauchyStressMat, trans(W)));
        CauchyRate = addA(JaumannRateMat, CauchyRate);
        CauchyStressMat = addA(CauchyStressMat,scaleA(Dlt,CauchyRate));
        SaveStress(CauchyStressMat, QuadraturePointId);

        if (BlockInfo.MaterialType == 3)
            Constitution(CauchyStressMat, GreenStrain, BlockInfo, QuadraturePointId);

        // EchoVar("VelocityLocal",VelocityLocal);
        // EchoVar("D",D);
        // EchoVar("JaumannRateVoigt",JaumannRateVoigt);


        // assemble to global; this can be moved out
        int threadId{};
        threadId = omp_get_thread_num();
        for (int i = 0; i < NodeNum; i++)
            for (int j = 0; j < 3; j++)
            {
                gxcData FintLocal{};
                for (auto p = 0; p < 3; p++)
                    FintLocal += UpShapeDxLocal[i][p] * CauchyStressMat[p][j];
                FintTemp[NeighborNode[i]][j][threadId] += FintLocal * Jacobian[QuadraturePointId] * QWeight;
            }
    }

    // FintTemp[threadId][NeighborNode[i]][j] += FintLocal[i][j] * Volume;
}

void gxcCell::FindGrid(gxcControl &SimulationParameter)
{
    for (gxcId i = 0; i < 3; i++)
    {
        if (Position[i] > SimulationParameter.ModelBound[1][i])
        {
            GridIndex[i] = SimulationParameter.GridNum[i] - 1;
        }
        else if (Position[i] < SimulationParameter.ModelBound[0][i])
        {
            GridIndex[i] = 0;
        }
        else
        {
            GridIndex[i] = floor((Position[i] - SimulationParameter.ModelBound[0][i]) / SimulationParameter.GridSize[i]);
        }
    }
    SimulationParameter.CellInGrid[GridIndex[0]][GridIndex[1]][GridIndex[2]].push_back(GlobalId);
}

inline void gxcCell::Constitution(std::array<std::array<gxcData, 3>, 3> &SPK, const std::array<std::array<gxcData, 3>, 3> &GreenStrain, const gxcBlockInfo &BlockInfo, gxcId QuadraturePointId)
{

    switch (BlockInfo.MaterialType)
    {
    case 1:
    {
        std::array<std::array<gxcData, 1>, 6> GreenStrainVoigt, SPKVoigt;
        Mat2VoigtStrain(GreenStrain, GreenStrainVoigt);
        SaveStrain(GreenStrainVoigt, QuadraturePointId);
        SPKVoigt = prodA(BlockInfo.Cmat, GreenStrainVoigt);
        Voigt2MatStress(SPKVoigt, SPK);
        break;
    }
    case 2:
    {
        gxcData Q1, I1, I2, I3, I1bar, I3T1, JT, qtemp, ptemp;

        std::array<std::array<gxcData, 3>, 3> GreenDeformation, InvG;

        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                GreenDeformation[i][j] = 2 * GreenStrain[i][j];

        GreenDeformation[0][0] += 1.;
        GreenDeformation[1][1] += 1.;
        GreenDeformation[2][2] += 1.;

        I1 = GreenDeformation[0][0] + GreenDeformation[1][1] + GreenDeformation[2][2];

        I2 = GreenDeformation[0][0] * GreenDeformation[1][1] + GreenDeformation[1][1] * GreenDeformation[2][2] + GreenDeformation[0][0] * GreenDeformation[2][2] - GreenDeformation[0][1] * GreenDeformation[1][0] - GreenDeformation[1][2] * GreenDeformation[2][1] - GreenDeformation[0][2] * GreenDeformation[2][0];

        I3 = determinant_3a(GreenDeformation);

        JT = sqrt(I3);

        I1bar = I1 / cbrt(I3);

        Q1 = BlockInfo.MaterialProperty[2] + 2 * BlockInfo.MaterialProperty[3] * (I1bar - 3.) + 3 * BlockInfo.MaterialProperty[4] * (I1bar - 3.) * (I1bar - 3.);

        qtemp = 2.0 * Q1 / cbrt(I3);

        ptemp = BlockInfo.MaterialProperty[1] * (JT - 1.0) * JT;

        InvMat3a(GreenDeformation, InvG);

        const gxcData ident[3][3]{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};

        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                SPK[i][j] = qtemp * (ident[i][j] - (1. / 3.) * I1 * InvG[i][j]) + ptemp * InvG[i][j];

        // EchoVar("InvG", InvG);
        break;
    }
    case 3:
    {
        gxcData MuM, Sig, B, Ce, eps;
        MuM = BlockInfo.MaterialProperty[19];
        // Sig = std::sqrt(2 / 3) * BlockInfo.MaterialProperty[5];
        Sig = BlockInfo.MaterialProperty[5];
        B = BlockInfo.MaterialProperty[6];
        Ce = BlockInfo.MaterialProperty[7];
        eps = State[QuadraturePointId][0];

        gxcData DevStress[6]{};
        gxcData SphStress{};
        gxcData DevStressNorm{};

        DeviatoricStress(Stress[QuadraturePointId], DevStress, SphStress);
        DevStressNorm = StressNorm(DevStress);

        gxcData YieldStress{}, YieldState{};

        YieldStress = Sig * std::pow((1 + B * eps), Ce);
        YieldState = DevStressNorm - YieldStress;

        if (YieldState > 0)
        {
            gxcData dGamma{}, dYieldStress{}, deltaGamma{};

            for (int i = 0; i < 20; i++)
            {

                YieldStress = Sig * std::pow((1 + B * eps), Ce);
                YieldState = std::sqrt(3 / 2) * DevStressNorm - 3 * MuM * deltaGamma - YieldStress;
                if (std::abs(YieldState) > Sig * 0.00001)
                {
                    if (i == 19)
                    {
                        // EchoVar("YieldStress", YieldStress);
                        EchoError("Plasticity iteration cannot converge");
                    }
                    // if (i == 1) EchoStr("second iteration!");
                    dYieldStress = B * Ce * Sig * std::pow((1 + B * eps), (Ce - 1));
                    dGamma = YieldState / (3 * MuM + dYieldStress);
                    eps += dGamma;
                    deltaGamma += dGamma;
                }
                else
                {
                    State[QuadraturePointId][0] = eps;
                    for (int j = 0; j < 6; j++)
                    {
                        DevStress[j] -= 3 * MuM * deltaGamma * DevStress[j] / DevStressNorm;
                    }

                    break;
                }
            }
        }

        Stress[QuadraturePointId][0] = DevStress[0] + SphStress;
        Stress[QuadraturePointId][1] = DevStress[1] + SphStress;
        Stress[QuadraturePointId][2] = DevStress[2] + SphStress;
        Stress[QuadraturePointId][3] = DevStress[3];
        Stress[QuadraturePointId][4] = DevStress[4];
        Stress[QuadraturePointId][5] = DevStress[5];

        break;
    }
    default:
        EchoStr("Material type not found!");
    }
}