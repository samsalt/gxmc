#include "gxcProblem.h"

void SetModelFileName(std::string &ModelFileName)
{
  std::ifstream ControlFile("control.dat");
  std::string line;
  bool HasFileName = 0;
  if (ControlFile)
  {
    while (ControlFile)
    {
      getline(ControlFile, line);
      if (line == "# model file name")
        getline(ControlFile, ModelFileName);
      HasFileName = 1;
    }
    if (!HasFileName)
      EchoStr("File name not included in control.dat");
    else
    EchoStr("Using model file: "+ ModelFileName);
  }
  else
    EchoStr("control.dat file not found");
  ControlFile.close();
}

void gxcProblem::ReadModel()
{
  SetModelFileName(SimulationParameter.FileName);
  bool MeshfreeSwitchLocal{false};
  int sn = SimulationParameter.FileName.length();
  char CharName[sn + 1];
  strcpy(CharName, SimulationParameter.FileName.c_str());

  int *ids;

  // open and create input and output files
  int CPU_word_size, IO_word_size, error, InputExoid, OutputExoid;
  IO_word_size = sizeof(gxcData);
  CPU_word_size = sizeof(gxcData);

  float version;
  InputExoid = ex_open(CharName, EX_READ, &CPU_word_size, &IO_word_size, &version);

#ifdef DEBUG
  EchoVar("IO_word_size", IO_word_size);
  EchoVar("CPU_word_size", CPU_word_size);
#endif

  CPU_word_size = sizeof(gxcData);
  IO_word_size = sizeof(gxcData);
  OutputExoid = ex_create("gxcOut.exo", EX_CLOBBER, &CPU_word_size, &IO_word_size);

  // read and write model information
  char title[MAX_LINE_LENGTH + 1];
  int DimensionNum;
  error = ex_get_init(InputExoid, title, &DimensionNum, &SimulationParameter.np, &SimulationParameter.nc, &SimulationParameter.BlockNum, &SimulationParameter.NodeSetNum, &SimulationParameter.SideSetNum);
  error = ex_put_init(OutputExoid, "gxcDatabase", 3, SimulationParameter.np, SimulationParameter.nc, SimulationParameter.BlockNum, 0, 0);
  // error = ex_put_qa(OutputExoid, QANum, QARecords);
  Block.resize(SimulationParameter.BlockNum);
  Node.resize(SimulationParameter.np);
  SimulationParameter.BlockCellId.resize(SimulationParameter.nc, std::vector<gxcId>(2, 0));

  EchoVar("Number of Blocks", SimulationParameter.BlockNum);
  EchoVar("Number of Points", SimulationParameter.np);
  EchoVar("Number of Cells", SimulationParameter.nc);

  // read and write node coordinates
  {
    gxcData *xcoo, *ycoo, *zcoo;
    xcoo = (gxcData *)calloc(SimulationParameter.np, sizeof(gxcData));
    ycoo = (gxcData *)calloc(SimulationParameter.np, sizeof(gxcData));
    zcoo = (gxcData *)calloc(SimulationParameter.np, sizeof(gxcData));
    error = ex_get_coord(InputExoid, xcoo, ycoo, zcoo);
    error = ex_put_coord(OutputExoid, xcoo, ycoo, zcoo);
    for (int i = 0; i < SimulationParameter.np; i++)
    {
      Node[i].InitCoordinate[0] = xcoo[i];
      Node[i].InitCoordinate[1] = ycoo[i];
      Node[i].InitCoordinate[2] = zcoo[i];
    }
    free(xcoo);
    free(ycoo);
    free(zcoo);
  }

  char **coord_names;
  coord_names = (char **) malloc(3);
  coord_names[0] = (char *) malloc(2);
  coord_names[1] = (char *) malloc(2);
  coord_names[2] = (char *) malloc(2);
  std::string tempStr {"X"};
  strcpy(coord_names[0], tempStr.c_str());
  tempStr="Y";
  strcpy(coord_names[1], tempStr.c_str());
  tempStr="Z";
  strcpy(coord_names[2], tempStr.c_str());
  error = ex_put_coord_names(OutputExoid, coord_names);

  // read all element blocks information and new cells in each block
  for (gxcId BlockId = 0; BlockId < SimulationParameter.BlockNum; BlockId++)
  {
    char elem_type[MAX_STR_LENGTH + 1];
    error = ex_get_elem_block(InputExoid, BlockId + 1, elem_type, &Block[BlockId].Info.CellNum, &Block[BlockId].Info.VertexPerCellNum, &Block[BlockId].Info.AttributeNum);
    error = ex_put_elem_block(OutputExoid, BlockId + 1, elem_type, Block[BlockId].Info.CellNum, Block[BlockId].Info.VertexPerCellNum, Block[BlockId].Info.AttributeNum);

    Block[BlockId].Cell.resize(Block[BlockId].Info.CellNum);

    int *connect;
    connect = (int *)calloc(Block[BlockId].Info.VertexPerCellNum * Block[BlockId].Info.CellNum, sizeof(int));
    error = ex_get_elem_conn(InputExoid, BlockId + 1, connect);
    error = ex_put_elem_conn(OutputExoid, BlockId + 1, connect);

    std::string elem_type_string;
    elem_type_string = elem_type;
    elem_type_string.erase(elem_type_string.find_last_not_of(" ") + 1);
    if (elem_type_string == "HEX8")
    {
      Block[BlockId].Info.CellType = 12;
      if (MeshfreeSwitchLocal)
      {
        for (int i = 0; i < Block[BlockId].Info.CellNum; i++)
        {
          Block[BlockId].Cell[i] = new HexCellRK;
          for (int j = 0; j < 8; j++)
            Block[BlockId].Cell[i]->AddVertexId(j, connect[i * 8 + j] - 1);
        }
      }
      else
      {
        for (int i = 0; i < Block[BlockId].Info.CellNum; i++)
        {
          Block[BlockId].Cell[i] = new HexCellFem();
          for (int j = 0; j < 8; j++)
            Block[BlockId].Cell[i]->AddVertexId(j, connect[i * 8 + j] - 1);
        }
      }
    }
    else if (elem_type_string == "TETRA")
    {
      Block[BlockId].Info.CellType = 10;
      for (int i = 0; i < Block[BlockId].Info.CellNum; i++)
        {
          Block[BlockId].Cell[i] = new TetCellFem();
          for (int j = 0; j < 4; j++)
            Block[BlockId].Cell[i]->AddVertexId(j, connect[i * 4 + j] - 1);
        }
    }
    else{
      EchoStr("Not implemented element type");
      EchoStr(elem_type_string);
    }
    free(connect);
  }

  // EchoStr("~~~~ Finish reading  geometry ~~~~");


  // get the global id of each cell
  gxcId GlobalIdLocal{};
  for (gxcId BlockId = 0; BlockId < SimulationParameter.BlockNum; BlockId++)
    for (gxcId i = 0; i < Block[BlockId].Info.CellNum; i++)
    {
      Block[BlockId].Cell[i]->SetGlobalId(GlobalIdLocal);
      SimulationParameter.BlockCellId[GlobalIdLocal][0] = BlockId;
      SimulationParameter.BlockCellId[GlobalIdLocal][1] = i;
      GlobalIdLocal++;
    }
  // EchoStr("Get global id.");

  NodeSetTemp.resize(SimulationParameter.NodeSetNum);
  for (int NodeSetId = 0; NodeSetId < SimulationParameter.NodeSetNum; NodeSetId++)
  {
    gxcId NodesNuminSet;
    int num_df_in_set; // i don't know what is this
    error = ex_get_node_set_param(InputExoid, NodeSetId + 1, &NodesNuminSet, &num_df_in_set);
    NodeSetTemp[NodeSetId].resize(NodesNuminSet);
    gxcId *node_list;
    node_list = (gxcId *)calloc(NodesNuminSet, sizeof(gxcId));
    error = ex_get_node_set(InputExoid, NodeSetId + 1, node_list);

    for (int i = 0; i < NodesNuminSet; i++)
      NodeSetTemp[NodeSetId][i] = node_list[i] - 1;
    free(node_list);
  }
  // EchoStr("~~~~ Finish reading node set ~~~~");
  EchoVar("Number of nodeset",SimulationParameter.NodeSetNum);

  // read side set
  ids = (int *)calloc(SimulationParameter.SideSetNum, sizeof(int));
  error = ex_get_side_set_ids(InputExoid, ids);
  FaceTemp.resize(SimulationParameter.SideSetNum);

  for (int i = 0; i < SimulationParameter.SideSetNum; i++)
  {
    int num_sides_in_set, num_df_in_set, num_elem_in_set;
    int *elem_list, *side_list;
    error = ex_get_side_set_param(InputExoid, ids[i], &num_sides_in_set,
                                  &num_df_in_set);
    /* Note: The # of elements is same as # of sides! there could be repeated elements*/
    num_elem_in_set = num_sides_in_set;
    FaceTemp[i].resize(num_sides_in_set,std::vector<gxcId>(2,0));
    // FaceTemp[i].FaceSetId=i;
    elem_list = (int *)calloc(num_elem_in_set, sizeof(int));
    side_list = (int *)calloc(num_sides_in_set, sizeof(int));
    error = ex_get_side_set(InputExoid, ids[i], elem_list, side_list);

    for (int j = 0; j < num_elem_in_set; j++)
    {
      FaceTemp[i][j][0]=side_list[j] - 1;
      FaceTemp[i][j][1]=elem_list[j] - 1;
    }    
    free(elem_list);
    free(side_list);
  }

  EchoVar("Number of Sideset",SimulationParameter.SideSetNum);

  error = ex_close(InputExoid);

  std::vector<std::string> NameList{"Displacement_X", "Displacement_Y", "Displacement_Z", "Velocity_X", "Velocity_Y", "Velocity_Z", "Acceleration_X", "Acceleration_Y", "Acceleration_Z"};
  for (int i = 0; i < 9; i++)
  {
    SimulationParameter.NodalVarNameOutput[i] = new char[NameList[i].length() + 1];
    strcpy(SimulationParameter.NodalVarNameOutput[i], NameList[i].c_str());
  }
  std::vector<std::string> NameList2{"Strain_XX", "Strain_YY", "Strain_ZZ", "Strain_YZ", "Strain_XZ", "Strain_XY", "Stress_XX", "Stress_YY", "Stress_ZZ", "Stress_YZ", "Stress_XZ", "Stress_XY","Effective Plastic Strain","Damage"};
  for (int i = 0; i < SimulationParameter.CellVarNumOutput; i++)
  {
    SimulationParameter.CellVarNameOutput[i] = new char[NameList2[i].length() + 1];
    strcpy(SimulationParameter.CellVarNameOutput[i], NameList2[i].c_str());
  }
  error = ex_put_var_param(OutputExoid, "n", SimulationParameter.NodalVarNumOutput);
  error = ex_put_var_names(OutputExoid, "n", SimulationParameter.NodalVarNumOutput, SimulationParameter.NodalVarNameOutput);
  error = ex_put_var_param(OutputExoid, "e", SimulationParameter.CellVarNumOutput);
  error = ex_put_var_names(OutputExoid, "e", SimulationParameter.CellVarNumOutput, SimulationParameter.CellVarNameOutput);
  error = ex_close(OutputExoid);
}