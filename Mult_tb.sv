`timescale 1ns / 1ps

module Mult_tb;

reg [63:0]Ax[15:0][1:0];
reg [63:0]Bx[15:0][1:0];
reg [63:0]Px[15:0][1:0];
real Pxb[15:0][1:0];
real Atemp[15:0];
real Btemp[15:0];

integer i, r;

Mult DUT2(Ax, Bx, Px);

initial
    begin
        #10;
        for(r = 0; r < 16; r = r + 1)
            begin
                #10;
                Atemp[r] = 2 * r + 0;
                Btemp[r] = 2 * r + 1;
                Ax[r][0] = $realtobits(Atemp[r]);       
                Bx[r][0] = $realtobits(Btemp[r]); 
                Ax[r][1] = 0;      
                Bx[r][1] = 0;   
            end
            
        for(i = 0 ; i < 16; i = i + 1)//for all cols of Cx
            begin
                Pxb[i][0] = $bitstoreal(Px[i][0][63:0]);
                Pxb[i][1] = $bitstoreal(Px[i][1][63:0]);
                #10;
                $display("# %d - In Complex Form: Real: %f, Imaginary: %f i", i, Pxb[i][0], Pxb[i][1]);
            end
        #10;
    end
endmodule
