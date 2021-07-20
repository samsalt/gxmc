#include "TetCellFem.h"

void TetCellFem::UpdateTotalShape(const std::vector<gxcNode> &Node, const gxcControl &SimulationParameter)
{
    for (int i = 0; i < 4; i++)
        NeighborNode[i] = VertexId[i];

    std::array<std::array<gxcData, 3>, 3> XDxi, InvXDxi;
    std::array<std::array<gxcData, 4>, 3> xyzel;

    // get vertex
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 4; j++)
            xyzel[i][j] = Node[VertexId[j]].InitCoordinate[i];

    for (auto QuadraturePointId = 0; QuadraturePointId < QuadratureNum; QuadraturePointId++)
    {

        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
                XDxi[i][j] = 0.0;
                for (int k = 0; k < 4; k++)
                    XDxi[i][j] += ShapeTet4Dxi[QuadraturePointId][k][j] * xyzel[i][k];
            }

        Jacobian[QuadraturePointId] = determinant_3a(XDxi);
        if (Jacobian[QuadraturePointId] < 0)
        {

            EchoStr("Jacobian is negative!");
            EchoVar("QuadraturePointId", QuadraturePointId);
            EchoVar("GlobalId", GlobalId);
            exit(1);
        }
        InvMat3a(XDxi, InvXDxi);

        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 3; j++)
            {
                TotalShapeDx[QuadraturePointId][i][j] = 0;
                for (int k = 0; k < 3; k++)
                    TotalShapeDx[QuadraturePointId][i][j] += InvXDxi[k][j] * ShapeTet4Dxi[QuadraturePointId][i][k];
            }

        for (int i = 0; i < 4; i++)
            TotalShape[QuadraturePointId][i] = ShapeTet4[QuadraturePointId][i];
    }

    // for debug
    // for (int i=0; i<3; i++)
    // {
    //     double sum=0;
    //     for (int j=0; j<8; j++)
    //     {
    //         sum+=TotalShapeDx(j,i);
    //     }

    //     std::cout<<sum<<"\n";
    // }
    // std::cout<<TotalShapeDx<<"\n";
}

void TetCellFem::UpdateShapeDx(const std::vector<gxcNode> &Node, std::vector<std::vector<gxcData>> &UpShapeDxLocal, const gxcId QuadraturePointId)
{
    std::array<std::array<gxcData, 3>, 3> xDxi, InvxDxi;
    std::array<std::array<gxcData, 4>, 3> xyzel;

    // get vertex
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 4; j++)
            xyzel[i][j] = Node[VertexId[j]].InitCoordinate[i] + Node[VertexId[j]].Displacement[i];

    // for (auto QuadraturePointId = 0; QuadraturePointId < QuadratureNum; QuadraturePointId++)
    // {
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
        {
            xDxi[i][j] = 0.0;
            for (int k = 0; k < 4; k++)
                xDxi[i][j] += ShapeTet4Dxi[QuadraturePointId][k][j] * xyzel[i][k];
        }

    Jacobian[QuadraturePointId] = determinant_3a(xDxi);
    InvMat3a(xDxi, InvxDxi);

    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 3; j++)
        {
            UpShapeDxLocal[i][j] = 0;
            for (int k = 0; k < 3; k++)
                UpShapeDxLocal[i][j] += InvxDxi[k][j] * ShapeTet4Dxi[QuadraturePointId][i][k];
        }
    // }
}
TetCellFem::TetCellFem()
{
    NodeNum = 4;
    SideNum = 4;

    QuadratureNum = 4;
    QWeight = 0.25;
    ParentVolume = 1;
    Jacobian.resize(QuadratureNum, 0);
    TotalShapeDx.resize(QuadratureNum, std::vector<std::vector<gxcData>>(4, std::vector<gxcData>(3, 0)));
    TotalShape.resize(QuadratureNum, std::vector<gxcData>(4, 0));
    Stress.resize(QuadratureNum, std::vector<gxcData>(6, 0));
    Strain.resize(QuadratureNum, std::vector<gxcData>(6, 0));
    State.resize(QuadratureNum, std::vector<gxcData>(5, 0));
}