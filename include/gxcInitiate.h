#pragma once
#include <iostream>
#include <string>
#include <vector>
#include "HexCellFem.h"
#include "TetCellFem.h"
#include "gxcType.h"


const gxcData HexCellFem::ShapeHexDxi[8][3] = {
    {-0.125, -0.125, 0.125},
    {0.125, -0.125, 0.125},
    {0.125, -0.125, -0.125},
    {-0.125, -0.125, -0.125},
    {-0.125, 0.125, 0.125},
    {0.125, 0.125, 0.125},
    {0.125, 0.125, -0.125},
    {-0.125, 0.125, -0.125}};

const gxcData HexCellFem::ShapeHex8Dxi[8][8][3] = {
{{-0.08333333333333331,-0.08333333333333331,0.3110042339640731},
{0.08333333333333331,-0.022329099369260218,0.08333333333333331},
{0.3110042339640731,-0.08333333333333331,-0.08333333333333331},
{-0.3110042339640731,-0.3110042339640731,-0.3110042339640731},
{-0.022329099369260218,0.08333333333333331,0.08333333333333331},
{0.022329099369260218,0.022329099369260218,0.022329099369260218},
{0.08333333333333331,0.08333333333333331,-0.022329099369260218},
{-0.08333333333333331,0.3110042339640731,-0.08333333333333331}},
{{-0.3110042339640731,-0.3110042339640731,0.3110042339640731},
{0.3110042339640731,-0.08333333333333331,0.08333333333333331},
{0.08333333333333331,-0.022329099369260218,-0.08333333333333331},
{-0.08333333333333331,-0.08333333333333331,-0.3110042339640731},
{-0.08333333333333331,0.3110042339640731,0.08333333333333331},
{0.08333333333333331,0.08333333333333331,0.022329099369260218},
{0.022329099369260218,0.022329099369260218,-0.022329099369260218},
{-0.022329099369260218,0.08333333333333331,-0.08333333333333331}},
{{-0.022329099369260218,-0.08333333333333331,0.08333333333333331},
{0.022329099369260218,-0.022329099369260218,0.022329099369260218},
{0.08333333333333331,-0.08333333333333331,-0.022329099369260218},
{-0.08333333333333331,-0.3110042339640731,-0.08333333333333331},
{-0.08333333333333331,0.08333333333333331,0.3110042339640731},
{0.08333333333333331,0.022329099369260218,0.08333333333333331},
{0.3110042339640731,0.08333333333333331,-0.08333333333333331},
{-0.3110042339640731,0.3110042339640731,-0.3110042339640731}},
{{-0.08333333333333331,-0.3110042339640731,0.08333333333333331},
{0.08333333333333331,-0.08333333333333331,0.022329099369260218},
{0.022329099369260218,-0.022329099369260218,-0.022329099369260218},
{-0.022329099369260218,-0.08333333333333331,-0.08333333333333331},
{-0.3110042339640731,0.3110042339640731,0.3110042339640731},
{0.3110042339640731,0.08333333333333331,0.08333333333333331},
{0.08333333333333331,0.022329099369260218,-0.08333333333333331},
{-0.08333333333333331,0.08333333333333331,-0.3110042339640731}},
{{-0.08333333333333331,-0.022329099369260218,0.08333333333333331},
{0.08333333333333331,-0.08333333333333331,0.3110042339640731},
{0.3110042339640731,-0.3110042339640731,-0.3110042339640731},
{-0.3110042339640731,-0.08333333333333331,-0.08333333333333331},
{-0.022329099369260218,0.022329099369260218,0.022329099369260218},
{0.022329099369260218,0.08333333333333331,0.08333333333333331},
{0.08333333333333331,0.3110042339640731,-0.08333333333333331},
{-0.08333333333333331,0.08333333333333331,-0.022329099369260218}},
{{-0.3110042339640731,-0.08333333333333331,0.08333333333333331},
{0.3110042339640731,-0.3110042339640731,0.3110042339640731},
{0.08333333333333331,-0.08333333333333331,-0.3110042339640731},
{-0.08333333333333331,-0.022329099369260218,-0.08333333333333331},
{-0.08333333333333331,0.08333333333333331,0.022329099369260218},
{0.08333333333333331,0.3110042339640731,0.08333333333333331},
{0.022329099369260218,0.08333333333333331,-0.08333333333333331},
{-0.022329099369260218,0.022329099369260218,-0.022329099369260218}},
{{-0.022329099369260218,-0.022329099369260218,0.022329099369260218},
{0.022329099369260218,-0.08333333333333331,0.08333333333333331},
{0.08333333333333331,-0.3110042339640731,-0.08333333333333331},
{-0.08333333333333331,-0.08333333333333331,-0.022329099369260218},
{-0.08333333333333331,0.022329099369260218,0.08333333333333331},
{0.08333333333333331,0.08333333333333331,0.3110042339640731},
{0.3110042339640731,0.3110042339640731,-0.3110042339640731},
{-0.3110042339640731,0.08333333333333331,-0.08333333333333331}},
{{-0.08333333333333331,-0.08333333333333331,0.022329099369260218},
{0.08333333333333331,-0.3110042339640731,0.08333333333333331},
{0.022329099369260218,-0.08333333333333331,-0.08333333333333331},
{-0.022329099369260218,-0.022329099369260218,-0.022329099369260218},
{-0.3110042339640731,0.08333333333333331,0.08333333333333331},
{0.3110042339640731,0.3110042339640731,0.3110042339640731},
{0.08333333333333331,0.08333333333333331,-0.3110042339640731},
{-0.08333333333333331,0.022329099369260218,-0.08333333333333331}}
};
const gxcData HexCellFem::ShapeHex8[8][8]={{0.13144585576580212,0.035220810900864506,0.13144585576580212,0.4905626121623441,0.035220810900864506,0.009437387837655926,0.035220810900864506,0.13144585576580212},
{0.4905626121623441,0.13144585576580212,0.035220810900864506,0.13144585576580212,0.13144585576580212,0.035220810900864506,0.009437387837655926,0.035220810900864506},
{0.035220810900864506,0.009437387837655926,0.035220810900864506,0.13144585576580212,0.13144585576580212,0.035220810900864506,0.13144585576580212,0.4905626121623441},
{0.13144585576580212,0.035220810900864506,0.009437387837655926,0.035220810900864506,0.4905626121623441,0.13144585576580212,0.035220810900864506,0.13144585576580212},
{0.035220810900864506,0.13144585576580212,0.4905626121623441,0.13144585576580212,0.009437387837655926,0.035220810900864506,0.13144585576580212,0.035220810900864506},
{0.13144585576580212,0.4905626121623441,0.13144585576580212,0.035220810900864506,0.035220810900864506,0.13144585576580212,0.035220810900864506,0.009437387837655926},
{0.009437387837655926,0.035220810900864506,0.13144585576580212,0.035220810900864506,0.035220810900864506,0.13144585576580212,0.4905626121623441,0.13144585576580212},
{0.035220810900864506,0.13144585576580212,0.035220810900864506,0.009437387837655926,0.13144585576580212,0.4905626121623441,0.13144585576580212,0.035220810900864506}};

const gxcData TetCellFem::ShapeTet4Dxi[4][4][3] ={{ {1,0,0},{0,1,0},{-1,-1,-1},{0,0,1}},
    {{1,0,0},{0,1,0},{-1,-1,-1},{0,0,1}},
    {{1,0,0},{0,1,0},{-1,-1,-1},{0,0,1}},
    {{1,0,0},{0,1,0},{-1,-1,-1},{0,0,1}}};
const gxcData TetCellFem::ShapeTet4[4][4]={{0.5854102,0.1381966,0.1381966,0.1381966},{0.1381966,0.5854102,0.1381966,0.1381966},{0.1381966,0.1381966,0.5854102,0.1381966},{0.1381966,0.1381966,0.1381966,0.5854102}};
// static const gxcData HexCellFem::WeightHex8[8]={1,1,1,1,1,1,1,1};