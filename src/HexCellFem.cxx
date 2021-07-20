#include "HexCellFem.h"

void HexCellFem::UpdateTotalShape(const std::vector<gxcNode> &Node, const gxcControl &SimulationParameter)
{
    for (int i = 0; i < 8; i++)
        NeighborNode[i] = VertexId[i];


    std::array<std::array<gxcData, 3>, 3> XDxi, InvXDxi;
    std::array<std::array<gxcData, 8>, 3> xyzel;

    // get vertex
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 8; j++)
            xyzel[i][j] = Node[VertexId[j]].InitCoordinate[i];

    for (auto QuadraturePointId = 0; QuadraturePointId < QuadratureNum; QuadraturePointId++)
    {

        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
                XDxi[i][j] = 0.0;
                for (int k = 0; k < 8; k++)
                    XDxi[i][j] += ShapeHex8Dxi[QuadraturePointId][k][j] * xyzel[i][k];
            }

        Jacobian[QuadraturePointId] = determinant_3a(XDxi);
        InvMat3a(XDxi, InvXDxi);

        for (int i = 0; i < 8; i++)
            for (int j = 0; j < 3; j++)
            {
                TotalShapeDx[QuadraturePointId][i][j] = 0;
                for (int k = 0; k < 3; k++)
                    TotalShapeDx[QuadraturePointId][i][j] += InvXDxi[k][j] * ShapeHex8Dxi[QuadraturePointId][i][k];
            }

        for (int i = 0; i < 8; i++)
            TotalShape[QuadraturePointId][i] = ShapeHex8[QuadraturePointId][i];
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

void HexCellFem::UpdateShapeDx(const std::vector<gxcNode> &Node, std::vector<std::vector<gxcData>> &UpShapeDxLocal,    const gxcId QuadraturePointId)
{
    std::array<std::array<gxcData, 3>, 3> xDxi, InvxDxi;
    std::array<std::array<gxcData, 8>, 3> xyzel;

    // get vertex
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 8; j++)
            xyzel[i][j] = Node[VertexId[j]].InitCoordinate[i] + Node[VertexId[j]].Displacement[i];

    // for (auto QuadraturePointId = 0; QuadraturePointId < QuadratureNum; QuadraturePointId++)
    // {
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
                xDxi[i][j] = 0.0;
                for (int k = 0; k < 8; k++)
                    xDxi[i][j] += ShapeHex8Dxi[QuadraturePointId][k][j] * xyzel[i][k];
            }

        Jacobian[QuadraturePointId] = determinant_3a(xDxi);
        InvMat3a(xDxi, InvxDxi);

        for (int i = 0; i < 8; i++)
            for (int j = 0; j < 3; j++)
            {
                UpShapeDxLocal[i][j] = 0;
                for (int k = 0; k < 3; k++)
                    UpShapeDxLocal[i][j] += InvxDxi[k][j] * ShapeHex8Dxi[QuadraturePointId][i][k];
            }
    // }
}
HexCellFem::HexCellFem()
{
    NodeNum = 8;
    SideNum = 6;

    QuadratureNum = 8;
    QWeight=1;
    ParentVolume=8;
    Jacobian.resize(QuadratureNum, 0);
    TotalShapeDx.resize(QuadratureNum, std::vector<std::vector<gxcData>>(8, std::vector<gxcData>(3, 0)));
    TotalShape.resize(QuadratureNum, std::vector<gxcData>(8, 0));
    Stress.resize(QuadratureNum, std::vector<gxcData>(6, 0));
    Strain.resize(QuadratureNum, std::vector<gxcData>(6, 0));
    State.resize(QuadratureNum, std::vector<gxcData>(5, 0));
}