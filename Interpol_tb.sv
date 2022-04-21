`timescale 1ns / 1ps

module Interpol_tb;

real C[15:0];
reg [63:0]Cxin[15:0][1:0];
reg [63:0]Cxout[15:0][1:0];
real Cxoutb[15:0][1:0];

integer i, j, k;

Interpol i2 (Cxin, Cxout);

initial
    begin
    #5;
        for(i = 0; i < 8; i = i + 1)
            begin
                C[i] = i + 1;
            end
        C[8] = 0;
        C[9] = 0;
        C[10] = 0;
        C[11] = 0;
        C[12] = 0;
        C[13] = 0;
        C[14] = 0;
        C[15] = 0;
        for(j = 0 ; j < 16; j = j + 1)
            begin
                Cxin[j][0] = $realtobits(C[j]);
                Cxin[j][1] = $realtobits(C[j]);
            end
        #10;
        for(k = 0 ; k < 16; k = k + 1)
            begin
                Cxoutb[k][0] = $bitstoreal(Cxout[k][0][63:0]);
                Cxoutb[k][1] = $bitstoreal(Cxout[k][1][63:0]);
                $display("# %d - In Complex Form: Real: %f, Imaginary: %f i", k, Cxoutb[k][0], Cxoutb[k][1]);
            end
    end

endmodule
