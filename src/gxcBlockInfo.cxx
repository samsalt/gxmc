#include <gxcBlockInfo.h>

void gxcBlockInfo::FormCmat()
{

    Cmat = {};

    YoungMod = MaterialProperty[1];
    PoissonRatio = MaterialProperty[2];
    LamdaM = YoungMod * PoissonRatio / ((1 + PoissonRatio) * (1 - 2 * PoissonRatio));
    MuM = 0.5 * YoungMod / (1 + PoissonRatio);
    gxcData LamdaPlus2Mu = 2 * MuM + LamdaM;
    MaterialProperty[18]=LamdaM;
    MaterialProperty[19]=MuM;


    Cmat[0][0] = LamdaPlus2Mu;
    Cmat[1][1] = LamdaPlus2Mu;
    Cmat[2][2] = LamdaPlus2Mu;
    Cmat[3][3] = MuM;
    Cmat[4][4] = MuM;
    Cmat[5][5] = MuM;
    Cmat[0][1] = LamdaM;
    Cmat[1][0] = LamdaM;
    Cmat[0][2] = LamdaM;
    Cmat[2][0] = LamdaM;
    Cmat[1][2] = LamdaM;
    Cmat[2][1] = LamdaM;
}

void gxcBlockInfo::UpdateBasicMaterialInfo()
{
    Density=MaterialProperty[0];   

}