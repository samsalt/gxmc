#include <iostream>
#include <HexCell.h>
#include "gxcEcho.h"
#include "math.h"

void HexCell::UpdateVolume(const std::vector<gxcNode> Node)
{
    // gxcData volume
    std::vector<std::vector<gxcData>> at(3, std::vector<gxcData>(3, 0));
    std::vector<std::vector<gxcData>> bt(3, std::vector<gxcData>(3, 0));
    std::vector<std::vector<gxcData>> ct(3, std::vector<gxcData>(3, 0));

    std::vector<std::vector<gxcData>> xyzel(8, std::vector<gxcData>(3, 0));
    // get vertex
    for (int i = 0; i < 8; i++)
        for (int j = 0; j < 3; j++)
            xyzel[i][j] = Node[VertexId[i]].InitCoordinate[j];

    for (int i = 0; i < 3; i++)
    {
        at[0][i] = xyzel[6][i] - xyzel[0][i];
        bt[0][i] = xyzel[6][i] - xyzel[0][i];
        ct[0][i] = xyzel[6][i] - xyzel[0][i];
    }
    for (int i = 0; i < 3; i++)
    {
        at[1][i] = xyzel[1][i] - xyzel[0][i];
        bt[1][i] = xyzel[4][i] - xyzel[0][i];
        ct[1][i] = xyzel[3][i] - xyzel[0][i];
    }
    for (int i = 0; i < 3; i++)
    {
        at[2][i] = xyzel[2][i] - xyzel[5][i];
        bt[2][i] = xyzel[5][i] - xyzel[7][i];
        ct[2][i] = xyzel[7][i] - xyzel[2][i];
    }

    Volume = (determinant_3(at) + determinant_3(bt) + determinant_3(ct)) / 6;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 8; j++)
        {
            Position[i] += xyzel[j][i];
        }
        Position[i] = Position[i] * 0.125;
    }
}