`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// 
// Create Date: 05/14/2016 01:20:37 AM
// Design Name: 
// Project Name: Fast Fourier Transform for Fast Multiplication
// Description: 
// 
// CMPE 297 Project
//////////////////////////////////////////////////////////////////////////////////

module topModule(input [3:0]A[7:0], [3:0]B[7:0], output reg [63:0]Fx);

wire [63:0]Px[15:0][1:0];
wire [63:0]Py[15:0][1:0];
wire [63:0]Multout[15:0][1:0];
wire [63:0]Intout[15:0][1:0];

//zero extend because product is 64-bit
DFT f1({4'h0,4'h0,4'h0,4'h0,4'h0,4'h0,4'h0,4'h0,A},Px);//DFT(A)
DFT f2({4'h0,4'h0,4'h0,4'h0,4'h0,4'h0,4'h0,4'h0,B},Py);//DFT(B)
Mult m1 (Px, Py, Multout); //Px*Py=Cx
Interpol i1 (Multout, Intout);//Intepol(Cx)
Flatten flat (Intout, Fx);//Flatten(Cx)

endmodule

/* returns value representation of polynomial 
p  := polynomial
Px := value representation of polynomial 
*/
module DFT(input reg [3:0]p[15:0],output reg [63:0]Px[15:0][1:0]);
       
real M[15:0][31:0]; //16wide;32deep -- halve multiply accumulates... DFT matrix (used Excel to compute nxn matrix; code to print it in verilog)
real temp;

integer i, j, k;

assign M[0][0] = 1.000000;     assign M[0][1] = 0.000000;     assign M[0][2] = 1.000000;     assign M[0][3] = 0.000000;     assign M[0][4] = 1.000000;     assign M[0][5] = 0.000000;     assign M[0][6] = 1.000000;     assign M[0][7] = 0.000000;     
assign M[0][8] = 1.000000;     assign M[0][9] = 0.000000;     assign M[0][10] = 1.000000;     assign M[0][11] = 0.000000;     assign M[0][12] = 1.000000;     assign M[0][13] = 0.000000;     assign M[0][14] = 1.000000;     assign M[0][15] = 0.000000;     
assign M[0][16] = 1.000000;     assign M[0][17] = 0.000000;     assign M[0][18] = 1.000000;     assign M[0][19] = 0.000000;     assign M[0][20] = 1.000000;     assign M[0][21] = 0.000000;     assign M[0][22] = 1.000000;     assign M[0][23] = 0.000000;     
assign M[0][24] = 1.000000;     assign M[0][25] = 0.000000;     assign M[0][26] = 1.000000;     assign M[0][27] = 0.000000;     assign M[0][28] = 1.000000;     assign M[0][29] = 0.000000;     assign M[0][30] = 1.000000;     assign M[0][31] = 0.000000;     
assign M[1][0] = 1.000000;     assign M[1][1] = 0.000000;     assign M[1][2] = 0.923880;     assign M[1][3] = -0.382683;     assign M[1][4] = 0.707107;     assign M[1][5] = -0.707107;     assign M[1][6] = 0.382683;     assign M[1][7] = -0.923880;     
assign M[1][8] = 0.000000;     assign M[1][9] = -1.000000;     assign M[1][10] = -0.382683;     assign M[1][11] = -0.923880;     assign M[1][12] = -0.707107;     assign M[1][13] = -0.707107;     assign M[1][14] = -0.923880;     assign M[1][15] = -0.382683;     
assign M[1][16] = -1.000000;     assign M[1][17] = 0.000000;     assign M[1][18] = -0.923880;     assign M[1][19] = 0.382683;     assign M[1][20] = -0.707107;     assign M[1][21] = 0.707107;     assign M[1][22] = -0.382683;     assign M[1][23] = 0.923880;     
assign M[1][24] = 0.000000;     assign M[1][25] = 1.000000;     assign M[1][26] = 0.382683;     assign M[1][27] = 0.923880;     assign M[1][28] = 0.707107;     assign M[1][29] = 0.707107;     assign M[1][30] = 0.923880;     assign M[1][31] = 0.382683;     
assign M[2][0] = 1.000000;     assign M[2][1] = 0.000000;     assign M[2][2] = 0.707107;     assign M[2][3] = -0.707107;     assign M[2][4] = 0.000000;     assign M[2][5] = -1.000000;     assign M[2][6] = -0.707107;     assign M[2][7] = -0.707107;     
assign M[2][8] = -1.000000;     assign M[2][9] = 0.000000;     assign M[2][10] = -0.707107;     assign M[2][11] = 0.707107;     assign M[2][12] = 0.000000;     assign M[2][13] = 1.000000;     assign M[2][14] = 0.707107;     assign M[2][15] = 0.707107;     
assign M[2][16] = 1.000000;     assign M[2][17] = 0.000000;     assign M[2][18] = 0.707107;     assign M[2][19] = -0.707107;     assign M[2][20] = 0.000000;     assign M[2][21] = -1.000000;     assign M[2][22] = -0.707107;     assign M[2][23] = -0.707107;     
assign M[2][24] = -1.000000;     assign M[2][25] = 0.000000;     assign M[2][26] = -0.707107;     assign M[2][27] = 0.707107;     assign M[2][28] = 0.000000;     assign M[2][29] = 1.000000;     assign M[2][30] = 0.707107;     assign M[2][31] = 0.707107;     
assign M[3][0] = 1.000000;     assign M[3][1] = 0.000000;     assign M[3][2] = 0.382683;     assign M[3][3] = -0.923880;     assign M[3][4] = -0.707107;     assign M[3][5] = -0.707107;     assign M[3][6] = -0.923880;     assign M[3][7] = 0.382683;     
assign M[3][8] = 0.000000;     assign M[3][9] = 1.000000;     assign M[3][10] = 0.923880;     assign M[3][11] = 0.382683;     assign M[3][12] = 0.707107;     assign M[3][13] = -0.707107;     assign M[3][14] = -0.382683;     assign M[3][15] = -0.923880;     
assign M[3][16] = -1.000000;     assign M[3][17] = 0.000000;     assign M[3][18] = -0.382683;     assign M[3][19] = 0.923880;     assign M[3][20] = 0.707107;     assign M[3][21] = 0.707107;     assign M[3][22] = 0.923880;     assign M[3][23] = -0.382683;     
assign M[3][24] = 0.000000;     assign M[3][25] = -1.000000;     assign M[3][26] = -0.923880;     assign M[3][27] = -0.382683;     assign M[3][28] = -0.707107;     assign M[3][29] = 0.707107;     assign M[3][30] = 0.382683;     assign M[3][31] = 0.923880;     
assign M[4][0] = 1.000000;     assign M[4][1] = 0.000000;     assign M[4][2] = 0.000000;     assign M[4][3] = -1.000000;     assign M[4][4] = -1.000000;     assign M[4][5] = 0.000000;     assign M[4][6] = 0.000000;     assign M[4][7] = 1.000000;     
assign M[4][8] = 1.000000;     assign M[4][9] = 0.000000;     assign M[4][10] = 0.000000;     assign M[4][11] = -1.000000;     assign M[4][12] = -1.000000;     assign M[4][13] = 0.000000;     assign M[4][14] = 0.000000;     assign M[4][15] = 1.000000;     
assign M[4][16] = 1.000000;     assign M[4][17] = 0.000000;     assign M[4][18] = 0.000000;     assign M[4][19] = -1.000000;     assign M[4][20] = -1.000000;     assign M[4][21] = 0.000000;     assign M[4][22] = 0.000000;     assign M[4][23] = 1.000000;     
assign M[4][24] = 1.000000;     assign M[4][25] = 0.000000;     assign M[4][26] = 0.000000;     assign M[4][27] = -1.000000;     assign M[4][28] = -1.000000;     assign M[4][29] = 0.000000;     assign M[4][30] = 0.000000;     assign M[4][31] = 1.000000;     
assign M[5][0] = 1.000000;     assign M[5][1] = 0.000000;     assign M[5][2] = -0.382683;     assign M[5][3] = -0.923880;     assign M[5][4] = -0.707107;     assign M[5][5] = 0.707107;     assign M[5][6] = 0.923880;     assign M[5][7] = 0.382683;     
assign M[5][8] = 0.000000;     assign M[5][9] = -1.000000;     assign M[5][10] = -0.923880;     assign M[5][11] = 0.382683;     assign M[5][12] = 0.707107;     assign M[5][13] = 0.707107;     assign M[5][14] = 0.382683;     assign M[5][15] = -0.923880;     
assign M[5][16] = -1.000000;     assign M[5][17] = 0.000000;     assign M[5][18] = 0.382683;     assign M[5][19] = 0.923880;     assign M[5][20] = 0.707107;     assign M[5][21] = -0.707107;     assign M[5][22] = -0.923880;     assign M[5][23] = -0.382683;     
assign M[5][24] = 0.000000;     assign M[5][25] = 1.000000;     assign M[5][26] = 0.923880;     assign M[5][27] = -0.382683;     assign M[5][28] = -0.707107;     assign M[5][29] = -0.707107;     assign M[5][30] = -0.382683;     assign M[5][31] = 0.923880;     
assign M[6][0] = 1.000000;     assign M[6][1] = 0.000000;     assign M[6][2] = -0.707107;     assign M[6][3] = -0.707107;     assign M[6][4] = 0.000000;     assign M[6][5] = 1.000000;     assign M[6][6] = 0.707107;     assign M[6][7] = -0.707107;     
assign M[6][8] = -1.000000;     assign M[6][9] = 0.000000;     assign M[6][10] = 0.707107;     assign M[6][11] = 0.707107;     assign M[6][12] = 0.000000;     assign M[6][13] = -1.000000;     assign M[6][14] = -0.707107;     assign M[6][15] = 0.707107;     
assign M[6][16] = 1.000000;     assign M[6][17] = 0.000000;     assign M[6][18] = -0.707107;     assign M[6][19] = -0.707107;     assign M[6][20] = 0.000000;     assign M[6][21] = 1.000000;     assign M[6][22] = 0.707107;     assign M[6][23] = -0.707107;     
assign M[6][24] = -1.000000;     assign M[6][25] = 0.000000;     assign M[6][26] = 0.707107;     assign M[6][27] = 0.707107;     assign M[6][28] = 0.000000;     assign M[6][29] = -1.000000;     assign M[6][30] = -0.707107;     assign M[6][31] = 0.707107;     
assign M[7][0] = 1.000000;     assign M[7][1] = 0.000000;     assign M[7][2] = -0.923880;     assign M[7][3] = -0.382683;     assign M[7][4] = 0.707107;     assign M[7][5] = 0.707107;     assign M[7][6] = -0.382683;     assign M[7][7] = -0.923880;     
assign M[7][8] = 0.000000;     assign M[7][9] = 1.000000;     assign M[7][10] = 0.382683;     assign M[7][11] = -0.923880;     assign M[7][12] = -0.707107;     assign M[7][13] = 0.707107;     assign M[7][14] = 0.923880;     assign M[7][15] = -0.382683;     
assign M[7][16] = -1.000000;     assign M[7][17] = 0.000000;     assign M[7][18] = 0.923880;     assign M[7][19] = 0.382683;     assign M[7][20] = -0.707107;     assign M[7][21] = -0.707107;     assign M[7][22] = 0.382683;     assign M[7][23] = 0.923880;     
assign M[7][24] = 0.000000;     assign M[7][25] = -1.000000;     assign M[7][26] = -0.382683;     assign M[7][27] = 0.923880;     assign M[7][28] = 0.707107;     assign M[7][29] = -0.707107;     assign M[7][30] = -0.923880;     assign M[7][31] = 0.382683;     
assign M[8][0] = 1.000000;     assign M[8][1] = 0.000000;     assign M[8][2] = -1.000000;     assign M[8][3] = 0.000000;     assign M[8][4] = 1.000000;     assign M[8][5] = 0.000000;     assign M[8][6] = -1.000000;     assign M[8][7] = 0.000000;     
assign M[8][8] = 1.000000;     assign M[8][9] = 0.000000;     assign M[8][10] = -1.000000;     assign M[8][11] = 0.000000;     assign M[8][12] = 1.000000;     assign M[8][13] = 0.000000;     assign M[8][14] = -1.000000;     assign M[8][15] = 0.000000;     
assign M[8][16] = 1.000000;     assign M[8][17] = 0.000000;     assign M[8][18] = -1.000000;     assign M[8][19] = 0.000000;     assign M[8][20] = 1.000000;     assign M[8][21] = 0.000000;     assign M[8][22] = -1.000000;     assign M[8][23] = 0.000000;     
assign M[8][24] = 1.000000;     assign M[8][25] = 0.000000;     assign M[8][26] = -1.000000;     assign M[8][27] = 0.000000;     assign M[8][28] = 1.000000;     assign M[8][29] = 0.000000;     assign M[8][30] = -1.000000;     assign M[8][31] = 0.000000;     
assign M[9][0] = 1.000000;     assign M[9][1] = 0.000000;     assign M[9][2] = -0.923880;     assign M[9][3] = 0.382683;     assign M[9][4] = 0.707107;     assign M[9][5] = -0.707107;     assign M[9][6] = -0.382683;     assign M[9][7] = 0.923880;     
assign M[9][8] = 0.000000;     assign M[9][9] = -1.000000;     assign M[9][10] = 0.382683;     assign M[9][11] = 0.923880;     assign M[9][12] = -0.707107;     assign M[9][13] = -0.707107;     assign M[9][14] = 0.923880;     assign M[9][15] = 0.382683;     
assign M[9][16] = -1.000000;     assign M[9][17] = 0.000000;     assign M[9][18] = 0.923880;     assign M[9][19] = -0.382683;     assign M[9][20] = -0.707107;     assign M[9][21] = 0.707107;     assign M[9][22] = 0.382683;     assign M[9][23] = -0.923880;     
assign M[9][24] = 0.000000;     assign M[9][25] = 1.000000;     assign M[9][26] = -0.382683;     assign M[9][27] = -0.923880;     assign M[9][28] = 0.707107;     assign M[9][29] = 0.707107;     assign M[9][30] = -0.923880;     assign M[9][31] = -0.382683;     
assign M[10][0] = 1.000000;     assign M[10][1] = 0.000000;     assign M[10][2] = -0.707107;     assign M[10][3] = 0.707107;     assign M[10][4] = 0.000000;     assign M[10][5] = -1.000000;     assign M[10][6] = 0.707107;     assign M[10][7] = 0.707107;     
assign M[10][8] = -1.000000;     assign M[10][9] = 0.000000;     assign M[10][10] = 0.707107;     assign M[10][11] = -0.707107;     assign M[10][12] = 0.000000;     assign M[10][13] = 1.000000;     assign M[10][14] = -0.707107;     assign M[10][15] = -0.707107;     
assign M[10][16] = 1.000000;     assign M[10][17] = 0.000000;     assign M[10][18] = -0.707107;     assign M[10][19] = 0.707107;     assign M[10][20] = 0.000000;     assign M[10][21] = -1.000000;     assign M[10][22] = 0.707107;     assign M[10][23] = 0.707107;     
assign M[10][24] = -1.000000;     assign M[10][25] = 0.000000;     assign M[10][26] = 0.707107;     assign M[10][27] = -0.707107;     assign M[10][28] = 0.000000;     assign M[10][29] = 1.000000;     assign M[10][30] = -0.707107;     assign M[10][31] = -0.707107;     
assign M[11][0] = 1.000000;     assign M[11][1] = 0.000000;     assign M[11][2] = -0.382683;     assign M[11][3] = 0.923880;     assign M[11][4] = -0.707107;     assign M[11][5] = -0.707107;     assign M[11][6] = 0.923880;     assign M[11][7] = -0.382683;     
assign M[11][8] = 0.000000;     assign M[11][9] = 1.000000;     assign M[11][10] = -0.923880;     assign M[11][11] = -0.382683;     assign M[11][12] = 0.707107;     assign M[11][13] = -0.707107;     assign M[11][14] = 0.382683;     assign M[11][15] = 0.923880;     
assign M[11][16] = -1.000000;     assign M[11][17] = 0.000000;     assign M[11][18] = 0.382683;     assign M[11][19] = -0.923880;     assign M[11][20] = 0.707107;     assign M[11][21] = 0.707107;     assign M[11][22] = -0.923880;     assign M[11][23] = 0.382683;     
assign M[11][24] = 0.000000;     assign M[11][25] = -1.000000;     assign M[11][26] = 0.923880;     assign M[11][27] = 0.382683;     assign M[11][28] = -0.707107;     assign M[11][29] = 0.707107;     assign M[11][30] = -0.382683;     assign M[11][31] = -0.923880;     
assign M[12][0] = 1.000000;     assign M[12][1] = 0.000000;     assign M[12][2] = 0.000000;     assign M[12][3] = 1.000000;     assign M[12][4] = -1.000000;     assign M[12][5] = 0.000000;     assign M[12][6] = 0.000000;     assign M[12][7] = -1.000000;     
assign M[12][8] = 1.000000;     assign M[12][9] = 0.000000;     assign M[12][10] = 0.000000;     assign M[12][11] = 1.000000;     assign M[12][12] = -1.000000;     assign M[12][13] = 0.000000;     assign M[12][14] = 0.000000;     assign M[12][15] = -1.000000;     
assign M[12][16] = 1.000000;     assign M[12][17] = 0.000000;     assign M[12][18] = 0.000000;     assign M[12][19] = 1.000000;     assign M[12][20] = -1.000000;     assign M[12][21] = 0.000000;     assign M[12][22] = 0.000000;     assign M[12][23] = -1.000000;     
assign M[12][24] = 1.000000;     assign M[12][25] = 0.000000;     assign M[12][26] = 0.000000;     assign M[12][27] = 1.000000;     assign M[12][28] = -1.000000;     assign M[12][29] = 0.000000;     assign M[12][30] = 0.000000;     assign M[12][31] = -1.000000;     
assign M[13][0] = 1.000000;     assign M[13][1] = 0.000000;     assign M[13][2] = 0.382683;     assign M[13][3] = 0.923880;     assign M[13][4] = -0.707107;     assign M[13][5] = 0.707107;     assign M[13][6] = -0.923880;     assign M[13][7] = -0.382683;     
assign M[13][8] = 0.000000;     assign M[13][9] = -1.000000;     assign M[13][10] = 0.923880;     assign M[13][11] = -0.382683;     assign M[13][12] = 0.707107;     assign M[13][13] = 0.707107;     assign M[13][14] = -0.382683;     assign M[13][15] = 0.923880;     
assign M[13][16] = -1.000000;     assign M[13][17] = 0.000000;     assign M[13][18] = -0.382683;     assign M[13][19] = -0.923880;     assign M[13][20] = 0.707107;     assign M[13][21] = -0.707107;     assign M[13][22] = 0.923880;     assign M[13][23] = 0.382683;     
assign M[13][24] = 0.000000;     assign M[13][25] = 1.000000;     assign M[13][26] = -0.923880;     assign M[13][27] = 0.382683;     assign M[13][28] = -0.707107;     assign M[13][29] = -0.707107;     assign M[13][30] = 0.382683;     assign M[13][31] = -0.923880;     
assign M[14][0] = 1.000000;     assign M[14][1] = 0.000000;     assign M[14][2] = 0.707107;     assign M[14][3] = 0.707107;     assign M[14][4] = 0.000000;     assign M[14][5] = 1.000000;     assign M[14][6] = -0.707107;     assign M[14][7] = 0.707107;     
assign M[14][8] = -1.000000;     assign M[14][9] = 0.000000;     assign M[14][10] = -0.707107;     assign M[14][11] = -0.707107;     assign M[14][12] = 0.000000;     assign M[14][13] = -1.000000;     assign M[14][14] = 0.707107;     assign M[14][15] = -0.707107;     
assign M[14][16] = 1.000000;     assign M[14][17] = 0.000000;     assign M[14][18] = 0.707107;     assign M[14][19] = 0.707107;     assign M[14][20] = 0.000000;     assign M[14][21] = 1.000000;     assign M[14][22] = -0.707107;     assign M[14][23] = 0.707107;     
assign M[14][24] = -1.000000;     assign M[14][25] = 0.000000;     assign M[14][26] = -0.707107;     assign M[14][27] = -0.707107;     assign M[14][28] = 0.000000;     assign M[14][29] = -1.000000;     assign M[14][30] = 0.707107;     assign M[14][31] = -0.707107;     
assign M[15][0] = 1.000000;     assign M[15][1] = 0.000000;     assign M[15][2] = 0.923880;     assign M[15][3] = 0.382683;     assign M[15][4] = 0.707107;     assign M[15][5] = 0.707107;     assign M[15][6] = 0.382683;     assign M[15][7] = 0.923880;     
assign M[15][8] = 0.000000;     assign M[15][9] = 1.000000;     assign M[15][10] = -0.382683;     assign M[15][11] = 0.923880;     assign M[15][12] = -0.707107;     assign M[15][13] = 0.707107;     assign M[15][14] = -0.923880;     assign M[15][15] = 0.382683;     
assign M[15][16] = -1.000000;     assign M[15][17] = 0.000000;     assign M[15][18] = -0.923880;     assign M[15][19] = -0.382683;     assign M[15][20] = -0.707107;     assign M[15][21] = -0.707107;     assign M[15][22] = -0.382683;     assign M[15][23] = -0.923880;     
assign M[15][24] = 0.000000;     assign M[15][25] = -1.000000;     assign M[15][26] = 0.382683;     assign M[15][27] = -0.923880;     assign M[15][28] = 0.707107;     assign M[15][29] = -0.707107;     assign M[15][30] = 0.923880;     assign M[15][31] = -0.382683;     

//  Px = p * M

always@(*)
begin
for(i = 0 ; i < 16; i = i + 1)//for all cols of Cx
    begin
        for(j = 0; j < 2; j = j + 1)//for all rows of Cx
        begin
           for(k = 0 ; k < 16; k = k + 1 )
           begin
                  temp = temp + (p[k]*M[i][2*k+j]);  
           end
           Px[i][j] = $realtobits(temp);
    //       $display("p[k] = %f", temp, p[k]);
           temp = 0;
        end
    end
end

endmodule

module Mult(input [63:0]Ax[15:0][1:0], input [63:0]Bx[15:0][1:0], output reg [63:0]Cx[15:0][1:0]);
integer i;
real Realtemp, Complextemp, AR, AC, BR, BC;

//normal multiplication component by component
always@(*)
begin
    for(i = 0; i < 16; i = i + 1)
        begin
            AR = $bitstoreal(Ax[i][0]);
            AC = $bitstoreal(Ax[i][1]);
            BR = $bitstoreal(Bx[i][0]);
            BC = $bitstoreal(Bx[i][1]);
            // (AR + ACi) * (BR + BCi) = AR*BR + AR*BCi + AC*BRi - AC*BC
            // i^2 = -1
            Realtemp = (AR * BR) - (AC * BC);
            Complextemp = (AR * BC) + (AC * BR);
            Cx[i][0] = $realtobits(Realtemp);
            Cx[i][1] = $realtobits(Complextemp);
        end
end

endmodule

//The International Criminal Police Organization - lol 
module Interpol(input reg [63:0]Cx[15:0][1:0],output reg [63:0]Px[15:0][1:0]);
       
real MI[15:0][31:0];
real Realtemp, Complextemp;
real AR, AC, BR, BC;

integer i, j;

assign MI[0][0] = 0.062500;     assign MI[0][1] = 0.000000;     assign MI[0][2] = 0.062500;     assign MI[0][3] = 0.000000;     assign MI[0][4] = 0.062500;     assign MI[0][5] = 0.000000;     assign MI[0][6] = 0.062500;     assign MI[0][7] = 0.000000;     
assign MI[0][8] = 0.062500;     assign MI[0][9] = 0.000000;     assign MI[0][10] = 0.062500;     assign MI[0][11] = 0.000000;     assign MI[0][12] = 0.062500;     assign MI[0][13] = 0.000000;     assign MI[0][14] = 0.062500;     assign MI[0][15] = 0.000000;     
assign MI[0][16] = 0.062500;     assign MI[0][17] = 0.000000;     assign MI[0][18] = 0.062500;     assign MI[0][19] = 0.000000;     assign MI[0][20] = 0.062500;     assign MI[0][21] = 0.000000;     assign MI[0][22] = 0.062500;     assign MI[0][23] = 0.000000;     
assign MI[0][24] = 0.062500;     assign MI[0][25] = 0.000000;     assign MI[0][26] = 0.062500;     assign MI[0][27] = 0.000000;     assign MI[0][28] = 0.062500;     assign MI[0][29] = 0.000000;     assign MI[0][30] = 0.062500;     assign MI[0][31] = 0.000000;     
assign MI[1][0] = 0.062500;     assign MI[1][1] = 0.000000;     assign MI[1][2] = 0.057742;     assign MI[1][3] = 0.023918;     assign MI[1][4] = 0.044194;     assign MI[1][5] = 0.044194;     assign MI[1][6] = 0.023918;     assign MI[1][7] = 0.057742;     
assign MI[1][8] = 0.000000;     assign MI[1][9] = 0.062500;     assign MI[1][10] = -0.023918;     assign MI[1][11] = 0.057742;     assign MI[1][12] = -0.044194;     assign MI[1][13] = 0.044194;     assign MI[1][14] = -0.057742;     assign MI[1][15] = 0.023918;     
assign MI[1][16] = -0.062500;     assign MI[1][17] = 0.000000;     assign MI[1][18] = -0.057742;     assign MI[1][19] = -0.023918;     assign MI[1][20] = -0.044194;     assign MI[1][21] = -0.044194;     assign MI[1][22] = -0.023918;     assign MI[1][23] = 0.023918;     
assign MI[1][24] = 0.000000;     assign MI[1][25] = -0.062500;     assign MI[1][26] = 0.023918;     assign MI[1][27] = -0.057742;     assign MI[1][28] = 0.044194;     assign MI[1][29] = -0.044194;     assign MI[1][30] = 0.057742;     assign MI[1][31] = -0.023918;     
assign MI[2][0] = 0.062500;     assign MI[2][1] = 0.000000;     assign MI[2][2] = 0.044194;     assign MI[2][3] = 0.044194;     assign MI[2][4] = 0.000000;     assign MI[2][5] = 0.062500;     assign MI[2][6] = -0.044194;     assign MI[2][7] = 0.044194;     
assign MI[2][8] = -0.062500;     assign MI[2][9] = 0.000000;     assign MI[2][10] = -0.044194;     assign MI[2][11] = -0.044194;     assign MI[2][12] = 0.000000;     assign MI[2][13] = -0.062500;     assign MI[2][14] = 0.044194;     assign MI[2][15] = -0.044194;     
assign MI[2][16] = 0.062500;     assign MI[2][17] = 0.000000;     assign MI[2][18] = 0.044194;     assign MI[2][19] = 0.044194;     assign MI[2][20] = 0.000000;     assign MI[2][21] = 0.062500;     assign MI[2][22] = -0.044194;     assign MI[2][23] = 0.044194;     
assign MI[2][24] = -0.062500;     assign MI[2][25] = 0.000000;     assign MI[2][26] = -0.044194;     assign MI[2][27] = -0.044194;     assign MI[2][28] = 0.000000;     assign MI[2][29] = -0.062500;     assign MI[2][30] = 0.044194;     assign MI[2][31] = -0.044194;     
assign MI[3][0] = 0.062500;     assign MI[3][1] = 0.000000;     assign MI[3][2] = 0.023918;     assign MI[3][3] = 0.057742;     assign MI[3][4] = -0.044194;     assign MI[3][5] = 0.044194;     assign MI[3][6] = -0.057742;     assign MI[3][7] = -0.023918;     
assign MI[3][8] = 0.000000;     assign MI[3][9] = -0.062500;     assign MI[3][10] = 0.057742;     assign MI[3][11] = -0.023918;     assign MI[3][12] = 0.044194;     assign MI[3][13] = 0.044194;     assign MI[3][14] = -0.023918;     assign MI[3][15] = 0.057742;     
assign MI[3][16] = -0.062500;     assign MI[3][17] = 0.000000;     assign MI[3][18] = -0.023918;     assign MI[3][19] = -0.057742;     assign MI[3][20] = 0.044194;     assign MI[3][21] = -0.044194;     assign MI[3][22] = 0.057742;     assign MI[3][23] = 0.057742;     
assign MI[3][24] = 0.000000;     assign MI[3][25] = 0.062500;     assign MI[3][26] = -0.057742;     assign MI[3][27] = 0.023918;     assign MI[3][28] = -0.044194;     assign MI[3][29] = -0.044194;     assign MI[3][30] = 0.023918;     assign MI[3][31] = -0.057742;     
assign MI[4][0] = 0.062500;     assign MI[4][1] = 0.000000;     assign MI[4][2] = 0.000000;     assign MI[4][3] = 0.062500;     assign MI[4][4] = -0.062500;     assign MI[4][5] = 0.000000;     assign MI[4][6] = 0.000000;     assign MI[4][7] = -0.062500;     
assign MI[4][8] = 0.062500;     assign MI[4][9] = 0.000000;     assign MI[4][10] = 0.000000;     assign MI[4][11] = 0.062500;     assign MI[4][12] = -0.062500;     assign MI[4][13] = 0.000000;     assign MI[4][14] = 0.000000;     assign MI[4][15] = -0.062500;     
assign MI[4][16] = 0.062500;     assign MI[4][17] = 0.000000;     assign MI[4][18] = 0.000000;     assign MI[4][19] = 0.062500;     assign MI[4][20] = -0.062500;     assign MI[4][21] = 0.000000;     assign MI[4][22] = 0.000000;     assign MI[4][23] = 0.062500;     
assign MI[4][24] = 0.062500;     assign MI[4][25] = 0.000000;     assign MI[4][26] = 0.000000;     assign MI[4][27] = 0.062500;     assign MI[4][28] = -0.062500;     assign MI[4][29] = 0.000000;     assign MI[4][30] = 0.000000;     assign MI[4][31] = -0.062500;     
assign MI[5][0] = 0.062500;     assign MI[5][1] = 0.000000;     assign MI[5][2] = -0.023918;     assign MI[5][3] = 0.057742;     assign MI[5][4] = -0.044194;     assign MI[5][5] = -0.044194;     assign MI[5][6] = 0.057742;     assign MI[5][7] = -0.023918;     
assign MI[5][8] = 0.000000;     assign MI[5][9] = 0.062500;     assign MI[5][10] = -0.057742;     assign MI[5][11] = -0.023918;     assign MI[5][12] = 0.044194;     assign MI[5][13] = -0.044194;     assign MI[5][14] = 0.023918;     assign MI[5][15] = 0.057742;     
assign MI[5][16] = -0.062500;     assign MI[5][17] = 0.000000;     assign MI[5][18] = 0.023918;     assign MI[5][19] = -0.057742;     assign MI[5][20] = 0.044194;     assign MI[5][21] = 0.044194;     assign MI[5][22] = -0.057742;     assign MI[5][23] = 0.057742;     
assign MI[5][24] = 0.000000;     assign MI[5][25] = -0.062500;     assign MI[5][26] = 0.057742;     assign MI[5][27] = 0.023918;     assign MI[5][28] = -0.044194;     assign MI[5][29] = 0.044194;     assign MI[5][30] = -0.023918;     assign MI[5][31] = -0.057742;     
assign MI[6][0] = 0.062500;     assign MI[6][1] = 0.000000;     assign MI[6][2] = -0.044194;     assign MI[6][3] = 0.044194;     assign MI[6][4] = 0.000000;     assign MI[6][5] = -0.062500;     assign MI[6][6] = 0.044194;     assign MI[6][7] = 0.044194;     
assign MI[6][8] = -0.062500;     assign MI[6][9] = 0.000000;     assign MI[6][10] = 0.044194;     assign MI[6][11] = -0.044194;     assign MI[6][12] = 0.000000;     assign MI[6][13] = 0.062500;     assign MI[6][14] = -0.044194;     assign MI[6][15] = -0.044194;     
assign MI[6][16] = 0.062500;     assign MI[6][17] = 0.000000;     assign MI[6][18] = -0.044194;     assign MI[6][19] = 0.044194;     assign MI[6][20] = 0.000000;     assign MI[6][21] = -0.062500;     assign MI[6][22] = 0.044194;     assign MI[6][23] = 0.044194;     
assign MI[6][24] = -0.062500;     assign MI[6][25] = 0.000000;     assign MI[6][26] = 0.044194;     assign MI[6][27] = -0.044194;     assign MI[6][28] = 0.000000;     assign MI[6][29] = 0.062500;     assign MI[6][30] = -0.044194;     assign MI[6][31] = -0.044194;     
assign MI[7][0] = 0.062500;     assign MI[7][1] = 0.000000;     assign MI[7][2] = -0.057742;     assign MI[7][3] = 0.023918;     assign MI[7][4] = 0.044194;     assign MI[7][5] = -0.044194;     assign MI[7][6] = -0.023918;     assign MI[7][7] = 0.057742;     
assign MI[7][8] = 0.000000;     assign MI[7][9] = -0.062500;     assign MI[7][10] = 0.023918;     assign MI[7][11] = 0.057742;     assign MI[7][12] = -0.044194;     assign MI[7][13] = -0.044194;     assign MI[7][14] = 0.057742;     assign MI[7][15] = 0.023918;     
assign MI[7][16] = -0.062500;     assign MI[7][17] = 0.000000;     assign MI[7][18] = 0.057742;     assign MI[7][19] = -0.023918;     assign MI[7][20] = -0.044194;     assign MI[7][21] = 0.044194;     assign MI[7][22] = 0.023918;     assign MI[7][23] = 0.023918;     
assign MI[7][24] = 0.000000;     assign MI[7][25] = 0.062500;     assign MI[7][26] = -0.023918;     assign MI[7][27] = -0.057742;     assign MI[7][28] = 0.044194;     assign MI[7][29] = 0.044194;     assign MI[7][30] = -0.057742;     assign MI[7][31] = -0.023918;     
assign MI[8][0] = 0.062500;     assign MI[8][1] = 0.000000;     assign MI[8][2] = -0.062500;     assign MI[8][3] = 0.000000;     assign MI[8][4] = 0.062500;     assign MI[8][5] = 0.000000;     assign MI[8][6] = -0.062500;     assign MI[8][7] = 0.000000;     
assign MI[8][8] = 0.062500;     assign MI[8][9] = 0.000000;     assign MI[8][10] = -0.062500;     assign MI[8][11] = 0.000000;     assign MI[8][12] = 0.062500;     assign MI[8][13] = 0.000000;     assign MI[8][14] = -0.062500;     assign MI[8][15] = 0.000000;     
assign MI[8][16] = 0.062500;     assign MI[8][17] = 0.000000;     assign MI[8][18] = -0.062500;     assign MI[8][19] = 0.000000;     assign MI[8][20] = 0.062500;     assign MI[8][21] = 0.000000;     assign MI[8][22] = -0.062500;     assign MI[8][23] = 0.000000;     
assign MI[8][24] = 0.062500;     assign MI[8][25] = 0.000000;     assign MI[8][26] = -0.062500;     assign MI[8][27] = 0.000000;     assign MI[8][28] = 0.062500;     assign MI[8][29] = 0.000000;     assign MI[8][30] = -0.062500;     assign MI[8][31] = 0.000000;     
assign MI[9][0] = 0.062500;     assign MI[9][1] = 0.000000;     assign MI[9][2] = -0.057742;     assign MI[9][3] = -0.023918;     assign MI[9][4] = 0.044194;     assign MI[9][5] = 0.044194;     assign MI[9][6] = -0.023918;     assign MI[9][7] = -0.057742;     
assign MI[9][8] = 0.000000;     assign MI[9][9] = 0.062500;     assign MI[9][10] = 0.023918;     assign MI[9][11] = -0.057742;     assign MI[9][12] = -0.044194;     assign MI[9][13] = 0.044194;     assign MI[9][14] = 0.057742;     assign MI[9][15] = -0.023918;     
assign MI[9][16] = -0.062500;     assign MI[9][17] = 0.000000;     assign MI[9][18] = 0.057742;     assign MI[9][19] = 0.023918;     assign MI[9][20] = -0.044194;     assign MI[9][21] = -0.044194;     assign MI[9][22] = 0.023918;     assign MI[9][23] = -0.023918;     
assign MI[9][24] = 0.000000;     assign MI[9][25] = -0.062500;     assign MI[9][26] = -0.023918;     assign MI[9][27] = 0.057742;     assign MI[9][28] = 0.044194;     assign MI[9][29] = -0.044194;     assign MI[9][30] = -0.057742;     assign MI[9][31] = 0.023918;     
assign MI[10][0] = 0.062500;     assign MI[10][1] = 0.000000;     assign MI[10][2] = -0.044194;     assign MI[10][3] = -0.044194;     assign MI[10][4] = 0.000000;     assign MI[10][5] = 0.062500;     assign MI[10][6] = 0.044194;     assign MI[10][7] = -0.044194;     
assign MI[10][8] = -0.062500;     assign MI[10][9] = 0.000000;     assign MI[10][10] = 0.044194;     assign MI[10][11] = 0.044194;     assign MI[10][12] = 0.000000;     assign MI[10][13] = -0.062500;     assign MI[10][14] = -0.044194;     assign MI[10][15] = 0.044194;     
assign MI[10][16] = 0.062500;     assign MI[10][17] = 0.000000;     assign MI[10][18] = -0.044194;     assign MI[10][19] = -0.044194;     assign MI[10][20] = 0.000000;     assign MI[10][21] = 0.062500;     assign MI[10][22] = 0.044194;     assign MI[10][23] = -0.044194;     
assign MI[10][24] = -0.062500;     assign MI[10][25] = 0.000000;     assign MI[10][26] = 0.044194;     assign MI[10][27] = 0.044194;     assign MI[10][28] = 0.000000;     assign MI[10][29] = -0.062500;     assign MI[10][30] = -0.044194;     assign MI[10][31] = 0.044194;     
assign MI[11][0] = 0.062500;     assign MI[11][1] = 0.000000;     assign MI[11][2] = -0.023918;     assign MI[11][3] = -0.057742;     assign MI[11][4] = -0.044194;     assign MI[11][5] = 0.044194;     assign MI[11][6] = 0.057742;     assign MI[11][7] = 0.023918;     
assign MI[11][8] = 0.000000;     assign MI[11][9] = -0.062500;     assign MI[11][10] = -0.057742;     assign MI[11][11] = 0.023918;     assign MI[11][12] = 0.044194;     assign MI[11][13] = 0.044194;     assign MI[11][14] = 0.023918;     assign MI[11][15] = -0.057742;     
assign MI[11][16] = -0.062500;     assign MI[11][17] = 0.000000;     assign MI[11][18] = 0.023918;     assign MI[11][19] = 0.057742;     assign MI[11][20] = 0.044194;     assign MI[11][21] = -0.044194;     assign MI[11][22] = -0.057742;     assign MI[11][23] = -0.057742;     
assign MI[11][24] = 0.000000;     assign MI[11][25] = 0.062500;     assign MI[11][26] = 0.057742;     assign MI[11][27] = -0.023918;     assign MI[11][28] = -0.044194;     assign MI[11][29] = -0.044194;     assign MI[11][30] = -0.023918;     assign MI[11][31] = 0.057742;     
assign MI[12][0] = 0.062500;     assign MI[12][1] = 0.000000;     assign MI[12][2] = 0.000000;     assign MI[12][3] = -0.062500;     assign MI[12][4] = -0.062500;     assign MI[12][5] = 0.000000;     assign MI[12][6] = 0.000000;     assign MI[12][7] = 0.062500;     
assign MI[12][8] = 0.062500;     assign MI[12][9] = 0.000000;     assign MI[12][10] = 0.000000;     assign MI[12][11] = -0.062500;     assign MI[12][12] = -0.062500;     assign MI[12][13] = 0.000000;     assign MI[12][14] = 0.000000;     assign MI[12][15] = 0.062500;     
assign MI[12][16] = 0.062500;     assign MI[12][17] = 0.000000;     assign MI[12][18] = 0.000000;     assign MI[12][19] = -0.062500;     assign MI[12][20] = -0.062500;     assign MI[12][21] = 0.000000;     assign MI[12][22] = 0.000000;     assign MI[12][23] = -0.062500;     
assign MI[12][24] = 0.062500;     assign MI[12][25] = 0.000000;     assign MI[12][26] = 0.000000;     assign MI[12][27] = -0.062500;     assign MI[12][28] = -0.062500;     assign MI[12][29] = 0.000000;     assign MI[12][30] = 0.000000;     assign MI[12][31] = 0.062500;     
assign MI[13][0] = 0.062500;     assign MI[13][1] = 0.000000;     assign MI[13][2] = 0.023918;     assign MI[13][3] = -0.057742;     assign MI[13][4] = -0.044194;     assign MI[13][5] = -0.044194;     assign MI[13][6] = -0.057742;     assign MI[13][7] = 0.023918;     
assign MI[13][8] = 0.000000;     assign MI[13][9] = 0.062500;     assign MI[13][10] = 0.057742;     assign MI[13][11] = 0.023918;     assign MI[13][12] = 0.044194;     assign MI[13][13] = -0.044194;     assign MI[13][14] = -0.023918;     assign MI[13][15] = -0.057742;     
assign MI[13][16] = -0.062500;     assign MI[13][17] = 0.000000;     assign MI[13][18] = -0.023918;     assign MI[13][19] = 0.057742;     assign MI[13][20] = 0.044194;     assign MI[13][21] = 0.044194;     assign MI[13][22] = 0.057742;     assign MI[13][23] = -0.057742;     
assign MI[13][24] = 0.000000;     assign MI[13][25] = -0.062500;     assign MI[13][26] = -0.057742;     assign MI[13][27] = -0.023918;     assign MI[13][28] = -0.044194;     assign MI[13][29] = 0.044194;     assign MI[13][30] = 0.023918;     assign MI[13][31] = 0.057742;     
assign MI[14][0] = 0.062500;     assign MI[14][1] = 0.000000;     assign MI[14][2] = 0.044194;     assign MI[14][3] = -0.044194;     assign MI[14][4] = 0.000000;     assign MI[14][5] = -0.062500;     assign MI[14][6] = -0.044194;     assign MI[14][7] = -0.044194;     
assign MI[14][8] = -0.062500;     assign MI[14][9] = 0.000000;     assign MI[14][10] = -0.044194;     assign MI[14][11] = 0.044194;     assign MI[14][12] = 0.000000;     assign MI[14][13] = 0.062500;     assign MI[14][14] = 0.044194;     assign MI[14][15] = 0.044194;     
assign MI[14][16] = 0.062500;     assign MI[14][17] = 0.000000;     assign MI[14][18] = 0.044194;     assign MI[14][19] = -0.044194;     assign MI[14][20] = 0.000000;     assign MI[14][21] = -0.062500;     assign MI[14][22] = -0.044194;     assign MI[14][23] = -0.044194;     
assign MI[14][24] = -0.062500;     assign MI[14][25] = 0.000000;     assign MI[14][26] = -0.044194;     assign MI[14][27] = 0.044194;     assign MI[14][28] = 0.000000;     assign MI[14][29] = 0.062500;     assign MI[14][30] = 0.044194;     assign MI[14][31] = 0.044194;     
assign MI[15][0] = 0.062500;     assign MI[15][1] = 0.000000;     assign MI[15][2] = 0.057742;     assign MI[15][3] = -0.023918;     assign MI[15][4] = 0.044194;     assign MI[15][5] = -0.044194;     assign MI[15][6] = 0.023918;     assign MI[15][7] = -0.057742;     
assign MI[15][8] = 0.000000;     assign MI[15][9] = -0.062500;     assign MI[15][10] = -0.023918;     assign MI[15][11] = -0.057742;     assign MI[15][12] = -0.044194;     assign MI[15][13] = -0.044194;     assign MI[15][14] = -0.057742;     assign MI[15][15] = -0.023918;     
assign MI[15][16] = -0.062500;     assign MI[15][17] = 0.000000;     assign MI[15][18] = -0.057742;     assign MI[15][19] = 0.023918;     assign MI[15][20] = -0.044194;     assign MI[15][21] = 0.044194;     assign MI[15][22] = -0.023918;     assign MI[15][23] = -0.023918;     
assign MI[15][24] = 0.000000;     assign MI[15][25] = 0.062500;     assign MI[15][26] = 0.023918;     assign MI[15][27] = 0.057742;     assign MI[15][28] = 0.044194;     assign MI[15][29] = 0.044194;     assign MI[15][30] = 0.057742;     assign MI[15][31] = 0.023918;
//  Px = p * M

always@(*)
begin
for(i = 0 ; i < 16; i = i + 1)//for all cols of Cx
    begin
    for(j = 0 ; j < 16; j = j + 1 )
        begin
            AR = $bitstoreal(Cx[j][0]);
            AC = $bitstoreal(Cx[j][1]);
            BR = $bitstoreal(MI[i][2*j]);
            BC = $bitstoreal(MI[i][2*j+1]);
            // (AR + ACi) * (BR + BCi) = AR*BR + AR*BCi + AC*BRi - AC*BC
            // i^2 = -1
            Realtemp = Realtemp + (AR * BR) - (AC * BC);
            Complextemp = Complextemp + (AR * BC) + (AC * BR);
        end
        Px[i][0] = $realtobits(Realtemp);
        Px[i][1] = $realtobits(Complextemp); 
        Realtemp = 0;
        Complextemp = 0;           
//        $display("# %d - In Complex Form: Real: %f, Imaginary: %f i", i, Realtemp, Complextemp);
    end
end

endmodule

module Flatten(input [63:0]Ax[15:0][1:0], output reg [63:0]Fx);

integer i;
real Realtemp, Powertemp, total;

//normal multiplication component by component

always@(*)
begin
    for(i = 0; i < 16; i = i + 1)
        begin
            Realtemp = $bitstoreal(Ax[i][0]);
            if(Realtemp > 1 || Realtemp < -1)
                begin
                    Realtemp = $bitstoreal(Ax[i][0]);
                    Powertemp = Realtemp * (10**i);
                    total = total + Powertemp;
                    $display("# %d - Total is %f, Realtemp is %f, Powertemp is %f", i, total, Realtemp, Powertemp);
                end
            else
                begin
                    $display("# %d - Total is %f, Realtemp is %f, SKIPPED", i, total, Realtemp);
                end

        end
    Fx = $realtobits(total);
    total = 0;
end

endmodule
