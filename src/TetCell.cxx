#include <iostream>
#include <TetCell.h>

using namespace std;

void TetCell::UpdateVolume(const std::vector<gxcNode> Node)
{
    std::array<std::array<int, 3>, 4> xyzel;
    // gxcData volume

    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 3; j++)
            xyzel[i][j] = Node[VertexId[i]].InitCoordinate[j];

    Volume = (xyzel[3][0] - xyzel[0][0]) * ((xyzel[1][1] - xyzel[0][1]) * (xyzel[2][2] - xyzel[0][2]) - (xyzel[1][2] - xyzel[0][2]) * (xyzel[2][1] - xyzel[0][1])) + (xyzel[3][1] - xyzel[0][1]) * ((xyzel[1][2] - xyzel[0][2]) * (xyzel[2][0] - xyzel[0][0]) - (xyzel[1][0] - xyzel[0][0]) * (xyzel[2][2] - xyzel[0][2])) + (xyzel[3][2] - xyzel[0][2]) * ((xyzel[1][0] - xyzel[0][0]) * (xyzel[2][1] - xyzel[0][1]) - (xyzel[1][1] - xyzel[0][1]) * (xyzel[2][0] - xyzel[0][0]));
    Volume=Volume/6;
}
