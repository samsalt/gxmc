// #pragma once
// #include <vector>
// #include "gxcType.h"

// class gxcMeshfreeInfo
// {
//     public:
//     std::vector<bool> SingularNode;
//     std::vector<std::vector<gxcData>> Win;
//     boost::numeric::ublas::matrix<gxcData> H0;
//     short Deg {1};
//     short MDimension {4};
//     gxcData NormalWin {1.2};
//     gxcData WinMax[3] {};
//     void StateInitiate(gxcId np)
//     {
//         SingularNode.resize(np,false);
//         Win.resize(np,std::vector<gxcData>(3,0));
//         H0.resize(1,MDimension);
//         H0(0,0)=1.0;
//     };
// };