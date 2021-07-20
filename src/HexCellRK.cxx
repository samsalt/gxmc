#include "HexCellRK.h"

// void GetHx(const gxcData (&Xs)[3], const short Deg, bool (&HSign)[3], boost::numeric::ublas::matrix<gxcData> &HxLocal, std::vector<boost::numeric::ublas::matrix<gxcData>> &DHxLocal)
// {
//     HSign[0] = signbit(Xs[0]);
//     HSign[1] = signbit(Xs[1]);
//     HSign[2] = signbit(Xs[2]);
//     switch (Deg)
//     {
//     case 1:
//         HxLocal(0, 0) = 1;
//         HxLocal(1, 0) = Xs[0];
//         HxLocal(2, 0) = Xs[1];
//         HxLocal(3, 0) = Xs[2];
//         DHxLocal[0](1, 0) = 1;
//         DHxLocal[1](2, 0) = 1;
//         DHxLocal[2](3, 0) = 1;
//         break;
//     default:
//         EchoVar("This RK degree is not supported:", Deg);
//         break;
//     }
// }

// void GetPhi(std::vector<gxcData> &Phi, std::vector<std::vector<gxcData>> &DPhi, const std::vector<gxcData> &WinLocal, const bool (&HSign)[3], const gxcData (&Zl)[3],gxcId i)
// {
//     gxcData PhiLocal[3]{};
//     gxcData DPhiLocal[3]{};

//     for (int i = 0; i < 3; i++)
//     {
//         if (Zl[i] <= 0.5)
//         {
//             PhiLocal[i] = 2.0 / 3.0 - 4 * Zl[i] * Zl[i] + 4 * Zl[i] * Zl[i] * Zl[i];
//             DPhiLocal[i] = -8 * Zl[i] + 12 * Zl[i] * Zl[i];
//             if (HSign[i])
//                 DPhiLocal[i] = -DPhiLocal[i];
//         }
//         else if ((Zl[i] >= 0.5) && (Zl[i] <= 1))
//         {
//             PhiLocal[i] = 4.0 / 3.0 - 4 * Zl[i] + 4 * Zl[i] * Zl[i] - 4.0 / 3.0 * Zl[i] * Zl[i] * Zl[i];
//             DPhiLocal[i] = -4 + 8 * Zl[i] - 4 * Zl[i] * Zl[i];
//             if (HSign[i])
//                 DPhiLocal[i] = -DPhiLocal[i];
//         }
//         else
//         {
//             PhiLocal[i] = 0;
//             DPhiLocal[0] = 0;
//             DPhiLocal[1] = 0;
//             DPhiLocal[2] = 0;
//         }
//     }
//     Phi[i] = PhiLocal[0] * PhiLocal[1] * PhiLocal[2];
//     DPhi[i][0] = PhiLocal[1] * PhiLocal[2] * DPhiLocal[0] / WinLocal[0];
//     DPhi[i][1] = PhiLocal[0] * PhiLocal[2] * DPhiLocal[1] / WinLocal[1];
//     DPhi[i][2] = PhiLocal[0] * PhiLocal[1] * DPhiLocal[2] / WinLocal[2];
// }

// void HexCellRK::UpdateTotalShape(const std::vector<gxcNode> &Node,const gxcControl &SimulationParameter)
// {
//     using namespace boost::numeric::ublas;

//     // zero_matrix<double> v30 (3, 1);
//     short MDimensionLcoal=SimulationParameter.MeshfreeInfo.MDimension;
//     zero_matrix<gxcData> VM0(MDimensionLcoal, 1);
//     zero_matrix<gxcData> M0(MDimensionLcoal, MDimensionLcoal);
//     matrix<gxcData> M(MDimensionLcoal, MDimensionLcoal),DM1(MDimensionLcoal, MDimensionLcoal),DM2(MDimensionLcoal, MDimensionLcoal),DM3(MDimensionLcoal, MDimensionLcoal);
//     matrix<gxcData> InvM(MDimensionLcoal, MDimensionLcoal),InvDM1(MDimensionLcoal, MDimensionLcoal),InvDM2(MDimensionLcoal, MDimensionLcoal),InvDM3(MDimensionLcoal, MDimensionLcoal), Hx0(MDimensionLcoal,1),B0(1,MDimensionLcoal),B1(1,MDimensionLcoal),B2(1,MDimensionLcoal),B3(1,MDimensionLcoal);
//     M=M0;
//     DM1=M0;
//     DM2=M0;
//     DM3=M0;
//     bool HSign[3];
//     std::vector<matrix<gxcData>> Hx;
//     Hx.resize(NodeNum, VM0);
//     std::vector<std::vector<matrix<gxcData>>> DHx;
//     DHx.resize(NodeNum, std::vector<matrix<gxcData>>(3, VM0));
//     std::vector<gxcData> Phi(NodeNum);
//     std::vector<std::vector<gxcData>> DPhi;
//     DPhi.resize(NodeNum,std::vector<gxcData>(3, 0));

//     gxcData Xs[3]{};
//     gxcData Zl[3]{};
//     gxcId NodeId{};
//     for (int i = 0; i < NodeNum; i++)
//     {
//         NodeId = NeighborNode[i];
//         Xs[0] = Position[0] - SimulationParameter.CellPosition[NodeId][0];
//         Xs[1] = Position[1] - SimulationParameter.CellPosition[NodeId][1];
//         Xs[2] = Position[2] - SimulationParameter.CellPosition[NodeId][2];
//         Zl[0]=std::abs(Xs[0]/SimulationParameter.MeshfreeInfo.Win[NodeId][0]);
//         Zl[1]=std::abs(Xs[1]/SimulationParameter.MeshfreeInfo.Win[NodeId][1]);
//         Zl[2]=std::abs(Xs[2]/SimulationParameter.MeshfreeInfo.Win[NodeId][2]);
//         GetHx(Xs,SimulationParameter.MeshfreeInfo.Deg,HSign,Hx[i],DHx[i]);
//         GetPhi(Phi, DPhi,SimulationParameter.MeshfreeInfo.Win[NodeId],HSign,Zl,i);

//         M=M+Phi[i]*prod(Hx[i],trans(Hx[i]));
//         DM1=DM1+Phi[i]*prod(DHx[i][0],trans(Hx[i]))+Phi[i]*prod(Hx[i],trans(DHx[i][0]))+DPhi[i][0]*prod(Hx[i],trans(Hx[i]));
//         DM2=DM2+Phi[i]*prod(DHx[i][1],trans(Hx[i]))+Phi[i]*prod(Hx[i],trans(DHx[i][1]))+DPhi[i][1]*prod(Hx[i],trans(Hx[i]));
//         DM3=DM3+Phi[i]*prod(DHx[i][2],trans(Hx[i]))+Phi[i]*prod(Hx[i],trans(DHx[i][2]))+DPhi[i][2]*prod(Hx[i],trans(Hx[i]));
//     }

//     InvertMatrix(M,InvM);
//     InvDM1=prod(-InvM,DM1);
//     InvDM1=prod(InvDM1,InvM);
//     InvDM2=prod(-InvM,DM2);
//     InvDM2=prod(InvDM2,InvM);
//     InvDM3=prod(-InvM,DM3);
//     InvDM3=prod(InvDM3,InvM);

//     // EchoVar("invM",InvM);

//     B0=prod(SimulationParameter.MeshfreeInfo.H0,InvM);
//     B1=prod(SimulationParameter.MeshfreeInfo.H0,InvDM1);
//     B2=prod(SimulationParameter.MeshfreeInfo.H0,InvDM2);
//     B3=prod(SimulationParameter.MeshfreeInfo.H0,InvDM3);

//     // gxcData sumshape {};
//     // gxcData sumshape1 {};
//     // gxcData sumshape2 {};
//     // gxcData sumshape3 {};

//     for (int i = 0; i < NodeNum; i++)
//     {
//         TotalShape[i]=Phi[i]*prod(B0,Hx[i])(0,0);
//         TotalShapeDx[i][0]=Phi[i]*prod(B1,Hx[i])(0,0)+Phi[i]*prod(B0,DHx[i][0])(0,0)+DPhi[i][0]*prod(B0,Hx[i])(0,0);
//         TotalShapeDx[i][1]=Phi[i]*prod(B2,Hx[i])(0,0)+Phi[i]*prod(B0,DHx[i][1])(0,0)+DPhi[i][1]*prod(B0,Hx[i])(0,0);
//         TotalShapeDx[i][2]=Phi[i]*prod(B3,Hx[i])(0,0)+Phi[i]*prod(B0,DHx[i][2])(0,0)+DPhi[i][2]*prod(B0,Hx[i])(0,0);
//         // sumshape+=TotalShape[i];
//         // sumshape1+=TotalShapeDx[i][0];
//         // sumshape2+=TotalShapeDx[i][1];
//         // sumshape3+=TotalShapeDx[i][2];
//     }

//     // EchoVar("sumshape",sumshape);
//     // EchoVar("TotalShape",TotalShape[0]);
//     // EchoVar("DM3",DM3);
//     // EchoVar("InvDM3",InvDM3);
//     // EchoVar("sumshape3",sumshape3);
//     // EchoVar("B3",B3);
//     // EchoVar("DM1",DM1);
//     // // EchoVar("M",M);
//     // // EchoVar("H0",SimulationParameter.MeshfreeInfo.H0);
//     // // EchoVar("B0",B0);
//     // EchoVar("invM",InvM);
    
// }

// void HexCellRK::FindNeighbor(const gxcControl &SimulationParameter)
// {
//     gxcId XId[3] {};
//     gxcId YId[3] {};
//     gxcId ZId[3] {};
//     std::vector<gxcId> NeighborList;
//     for (gxcId i = 0; i<3 ; i++)
//     {
//         XId[i]=GridIndex[0]-1+i;
//         YId[i]=GridIndex[1]-1+i;
//         ZId[i]=GridIndex[2]-1+i;
//     }

//     for (gxcId i = 0; i < 3; i++)
//     {
//         if ((XId[i] >= 0) && (XId[i] < SimulationParameter.GridNum[0]))
//         {
//             for (gxcId j = 0; j < 3; j++)
//             {
//                 if ((YId[j] >= 0) && (YId[j] < SimulationParameter.GridNum[1]))
//                 {
//                     for (gxcId k = 0; k < 3; k++)
//                     {
//                         if ((ZId[k] >= 0) && (ZId[k] < SimulationParameter.GridNum[2]))
//                         {
//                            HardSearch(SimulationParameter.CellInGrid[XId[i]][YId[j]][ZId[k]],SimulationParameter);
//                         }
//                     }
//                 }
//             }
//         }
//     }
// }


// void HexCellRK::HardSearch(const std::vector<gxcId> &NeighborList,const gxcControl& Simulationparameter)
// {
//     for (auto nid:NeighborList)
//     {
//         gxcData NeiCoor[3] {};
//         NeiCoor[0]=Simulationparameter.CellPosition[nid][0];
//         NeiCoor[1]=Simulationparameter.CellPosition[nid][1];
//         NeiCoor[2]=Simulationparameter.CellPosition[nid][2];
//         if (std::abs(NeiCoor[0]-Position[0])<Simulationparameter.MeshfreeInfo.Win[nid][0])
//             if (std::abs(NeiCoor[1]-Position[1])<Simulationparameter.MeshfreeInfo.Win[nid][1])
//                 if (std::abs(NeiCoor[2]-Position[2])<Simulationparameter.MeshfreeInfo.Win[nid][2])
//                 {
//                     NeighborNode[NodeNum]=nid;
//                     NodeNum++;
//                 }
//     }
// }