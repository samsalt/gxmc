#include "PostProcessing.h"

void Output(std::vector<gxcBlock> &Block, std::vector<gxcNode> &Node, gxcControl &SimulationParameter)
{
    int CPU_word_size, IO_word_size, exoid, error;
    float version;
    CPU_word_size = sizeof(gxcData);
    IO_word_size = sizeof(gxcData);

    exoid = ex_open("gxcOut.exo", EX_WRITE, &CPU_word_size, &IO_word_size, &version);
    
    error = ex_put_time(exoid, SimulationParameter.TimeOutputStepID, &SimulationParameter.TimeCurrent);

    // write nodal values
    gxcData *dsp_x;
    gxcData *vel_x;
    gxcData *acl_x;
    gxcData *dsp_y;
    gxcData *vel_y;
    gxcData *acl_y;
    gxcData *dsp_z;
    gxcData *vel_z;
    gxcData *acl_z;
    dsp_x = (gxcData *)calloc(SimulationParameter.np, sizeof(gxcData));
    vel_x = (gxcData *)calloc(SimulationParameter.np, sizeof(gxcData));
    acl_x = (gxcData *)calloc(SimulationParameter.np, sizeof(gxcData));
    dsp_y = (gxcData *)calloc(SimulationParameter.np, sizeof(gxcData));
    vel_y = (gxcData *)calloc(SimulationParameter.np, sizeof(gxcData));
    acl_y = (gxcData *)calloc(SimulationParameter.np, sizeof(gxcData));
    dsp_z = (gxcData *)calloc(SimulationParameter.np, sizeof(gxcData));
    vel_z = (gxcData *)calloc(SimulationParameter.np, sizeof(gxcData));
    acl_z = (gxcData *)calloc(SimulationParameter.np, sizeof(gxcData));
    for (int i = 0; i < SimulationParameter.np; i++)
    {
        dsp_x[i] = Node[i].Displacement[0];
        vel_x[i] = Node[i].Velocity[0];
        acl_x[i] = Node[i].Acceleration[0];
        dsp_y[i] = Node[i].Displacement[1];
        vel_y[i] = Node[i].Velocity[1];
        acl_y[i] = Node[i].Acceleration[1];
        dsp_z[i] = Node[i].Displacement[2];
        vel_z[i] = Node[i].Velocity[2];
        acl_z[i] = Node[i].Acceleration[2];
    }

    error = ex_put_nodal_var(exoid, SimulationParameter.TimeOutputStepID, 1, SimulationParameter.np, dsp_x);
    error = ex_put_nodal_var(exoid, SimulationParameter.TimeOutputStepID, 2, SimulationParameter.np, dsp_y);
    error = ex_put_nodal_var(exoid, SimulationParameter.TimeOutputStepID, 3, SimulationParameter.np, dsp_z);
    error = ex_put_nodal_var(exoid, SimulationParameter.TimeOutputStepID, 4, SimulationParameter.np, vel_x);
    error = ex_put_nodal_var(exoid, SimulationParameter.TimeOutputStepID, 5, SimulationParameter.np, vel_y);
    error = ex_put_nodal_var(exoid, SimulationParameter.TimeOutputStepID, 6, SimulationParameter.np, vel_z);
    error = ex_put_nodal_var(exoid, SimulationParameter.TimeOutputStepID, 7, SimulationParameter.np, acl_x);
    error = ex_put_nodal_var(exoid, SimulationParameter.TimeOutputStepID, 8, SimulationParameter.np, acl_y);
    error = ex_put_nodal_var(exoid, SimulationParameter.TimeOutputStepID, 9, SimulationParameter.np, acl_z);

    free(dsp_x);
    free(vel_x);
    free(acl_x);
    free(dsp_y);
    free(vel_y);
    free(acl_y);
    free(dsp_z);
    free(vel_z);
    free(acl_z);

    // write cell values

    // output cell var

    if (!SimulationParameter.MeshfreeSwitch)
    {

        // gxcId CurrentCell{};
        for (int BlockId = 0; BlockId < SimulationParameter.BlockNum; BlockId++)
        {
            gxcData *strain_xx;
            gxcData *strain_yy;
            gxcData *strain_zz;
            gxcData *strain_yz;
            gxcData *strain_xz;
            gxcData *strain_xy;
            strain_xx = (gxcData *)calloc(Block[BlockId].Info.CellNum, sizeof(gxcData));
            strain_yy = (gxcData *)calloc(Block[BlockId].Info.CellNum, sizeof(gxcData));
            strain_zz = (gxcData *)calloc(Block[BlockId].Info.CellNum, sizeof(gxcData));
            strain_yz = (gxcData *)calloc(Block[BlockId].Info.CellNum, sizeof(gxcData));
            strain_xz = (gxcData *)calloc(Block[BlockId].Info.CellNum, sizeof(gxcData));
            strain_xy = (gxcData *)calloc(Block[BlockId].Info.CellNum, sizeof(gxcData));
            gxcData *stress_xx;
            gxcData *stress_yy;
            gxcData *stress_zz;
            gxcData *stress_yz;
            gxcData *stress_xz;
            gxcData *stress_xy;
            stress_xx = (gxcData *)calloc(Block[BlockId].Info.CellNum, sizeof(gxcData));
            stress_yy = (gxcData *)calloc(Block[BlockId].Info.CellNum, sizeof(gxcData));
            stress_zz = (gxcData *)calloc(Block[BlockId].Info.CellNum, sizeof(gxcData));
            stress_yz = (gxcData *)calloc(Block[BlockId].Info.CellNum, sizeof(gxcData));
            stress_xz = (gxcData *)calloc(Block[BlockId].Info.CellNum, sizeof(gxcData));
            stress_xy = (gxcData *)calloc(Block[BlockId].Info.CellNum, sizeof(gxcData));
            gxcData *eps;
            gxcData *damage;
            eps = (gxcData *)calloc(Block[BlockId].Info.CellNum, sizeof(gxcData));
            damage = (gxcData *)calloc(Block[BlockId].Info.CellNum, sizeof(gxcData));
            for (int i = 0; i < Block[BlockId].Info.CellNum; i++)
            {
                strain_xx[i] = Block[BlockId].Cell[i]->Strain[0][0];
                strain_yy[i] = Block[BlockId].Cell[i]->Strain[0][1];
                strain_zz[i] = Block[BlockId].Cell[i]->Strain[0][2];
                strain_yz[i] = Block[BlockId].Cell[i]->Strain[0][3];
                strain_xz[i] = Block[BlockId].Cell[i]->Strain[0][4];
                strain_xy[i] = Block[BlockId].Cell[i]->Strain[0][5];
                stress_xx[i] = Block[BlockId].Cell[i]->Stress[0][0];
                stress_yy[i] = Block[BlockId].Cell[i]->Stress[0][1];
                stress_zz[i] = Block[BlockId].Cell[i]->Stress[0][2];
                stress_yz[i] = Block[BlockId].Cell[i]->Stress[0][3];
                stress_xz[i] = Block[BlockId].Cell[i]->Stress[0][4];
                stress_xy[i] = Block[BlockId].Cell[i]->Stress[0][5];
                eps[i] = Block[BlockId].Cell[i]->State[0][0];
                damage[i] = Block[BlockId].Cell[i]->State[0][1];
                // CurrentCell++;
            }
            error = ex_put_elem_var(exoid, SimulationParameter.TimeOutputStepID, 1, BlockId + 1, Block[BlockId].Info.CellNum, strain_xx);
            error = ex_put_elem_var(exoid, SimulationParameter.TimeOutputStepID, 2, BlockId + 1, Block[BlockId].Info.CellNum, strain_yy);
            error = ex_put_elem_var(exoid, SimulationParameter.TimeOutputStepID, 3, BlockId + 1, Block[BlockId].Info.CellNum, strain_zz);
            error = ex_put_elem_var(exoid, SimulationParameter.TimeOutputStepID, 4, BlockId + 1, Block[BlockId].Info.CellNum, strain_yz);
            error = ex_put_elem_var(exoid, SimulationParameter.TimeOutputStepID, 5, BlockId + 1, Block[BlockId].Info.CellNum, strain_xz);
            error = ex_put_elem_var(exoid, SimulationParameter.TimeOutputStepID, 6, BlockId + 1, Block[BlockId].Info.CellNum, strain_xy);

            error = ex_put_elem_var(exoid, SimulationParameter.TimeOutputStepID, 7, BlockId + 1, Block[BlockId].Info.CellNum, stress_xx);
            error = ex_put_elem_var(exoid, SimulationParameter.TimeOutputStepID, 8, BlockId + 1, Block[BlockId].Info.CellNum, stress_yy);
            error = ex_put_elem_var(exoid, SimulationParameter.TimeOutputStepID, 9, BlockId + 1, Block[BlockId].Info.CellNum, stress_zz);
            error = ex_put_elem_var(exoid, SimulationParameter.TimeOutputStepID, 10, BlockId + 1, Block[BlockId].Info.CellNum, stress_yz);
            error = ex_put_elem_var(exoid, SimulationParameter.TimeOutputStepID, 11, BlockId + 1, Block[BlockId].Info.CellNum, stress_xz);
            error = ex_put_elem_var(exoid, SimulationParameter.TimeOutputStepID, 12, BlockId + 1, Block[BlockId].Info.CellNum, stress_xy);
            error = ex_put_elem_var(exoid, SimulationParameter.TimeOutputStepID, 13, BlockId + 1, Block[BlockId].Info.CellNum, eps);
            error = ex_put_elem_var(exoid, SimulationParameter.TimeOutputStepID, 14, BlockId + 1, Block[BlockId].Info.CellNum, damage);
            free(strain_xx);
            free(strain_yy);
            free(strain_zz);
            free(strain_yz);
            free(strain_xz);
            free(strain_xy);
            free(stress_xx);
            free(stress_yy);
            free(stress_zz);
            free(stress_yz);
            free(stress_xz);
            free(stress_xy);
            free(eps);
            free(damage);
        }
    }

    error = ex_close(exoid);

    // update time information
    SimulationParameter.TimetoOutput += SimulationParameter.TimeOutputPeriod;
    SimulationParameter.TimeOutputStepID++;
    std::cout << std::fixed;
    std::cout<<"time="<<std::setw(6) <<SimulationParameter.TimeCurrent<<"   Internal energy="<<std::setprecision(6)<<SimulationParameter.InternalEnergy<<"   Kinetic energy="<<std::setprecision(6)<<SimulationParameter.KineticEnergy<<"   Total energy="<<std::setprecision(6)<<SimulationParameter.KineticEnergy+SimulationParameter.InternalEnergy<<std::endl;
    std::fstream logFile("gxmc.log", std::fstream::app);
    logFile << std::fixed;
    logFile << "time="<<std::setw(6) <<SimulationParameter.TimeCurrent<<"   Internal energy="<<std::setprecision(6)<<SimulationParameter.InternalEnergy<<"   Kinetic energy="<<std::setprecision(6)<<SimulationParameter.KineticEnergy<<"   Total energy="<<std::setprecision(6)<<SimulationParameter.KineticEnergy+SimulationParameter.InternalEnergy<<std::endl;
    logFile.close();
    // EchoVar("Time",SimulationParameter.TimeCurrent);
}

// void MeshfreeOutputInit(   const gxcControl &SimulationParameter)
// {
//     int CPU_word_size, IO_word_size, error, OutputExoid, QANum{1};
//     IO_word_size = sizeof(gxcData);
//     CPU_word_size = sizeof(gxcData);

//     CPU_word_size = sizeof(gxcData);
//     IO_word_size = sizeof(gxcData);
//     OutputExoid = ex_create("gxcOut.exo", EX_CLOBBER, &CPU_word_size, &IO_word_size);
//     error = ex_put_init(OutputExoid, "gxcDatabase", 3, SimulationParameter.np, SimulationParameter.nc, 1, 0, 0);

//     gxcData *xcoo, *ycoo, *zcoo;
//     xcoo = (gxcData *)calloc(SimulationParameter.np, sizeof(gxcData));
//     ycoo = (gxcData *)calloc(SimulationParameter.np, sizeof(gxcData));
//     zcoo = (gxcData *)calloc(SimulationParameter.np, sizeof(gxcData));
//     for (gxcId i = 0; i < SimulationParameter.np; i++)
//     {
//         xcoo[i] = SimulationParameter.CellPosition[i][0];
//         ycoo[i] = SimulationParameter.CellPosition[i][1];
//         zcoo[i] = SimulationParameter.CellPosition[i][2];
//     }
//     error = ex_put_coord(OutputExoid, xcoo, ycoo, zcoo);
//     free(xcoo);
//     free(ycoo);
//     free(zcoo);

//     char *coord_names[2];
//     coord_names[0] = "X";
//     coord_names[1] = "Y";
//     coord_names[2] = "Z";
//     error = ex_put_coord_names(OutputExoid, coord_names);
//     error = ex_put_elem_block(OutputExoid, 1, "SPHERE", SimulationParameter.np, 1, 1);

//     int *connect;
//     connect = (int *)calloc(SimulationParameter.np, sizeof(int));
//     for (gxcId i = 0; i < SimulationParameter.np; i++)
//     {
//         connect[i] = i + 1;
//     }
//     error = ex_put_elem_conn(OutputExoid, 1, connect);
//     free(connect);
//     error = ex_close(OutputExoid);
// }