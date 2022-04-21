`timescale 1ns / 1ps

module DFT_tb;

reg [3:0]x[15:0];
reg [63:0]Px[15:0][1:0];
real Pxr[15:0][1:0];

integer i;

DFT f1(x,Px);

initial
    begin
    #5;
        for(i = 0 ; i < 16; i = i + 1)
            begin
                Pxr[i][0] = $bitstoreal(Px[i][0][63:0]);
                Pxr[i][1] = $bitstoreal(Px[i][1][63:0]);
                $display("# %d - In Complex Form: Real: %f, Imaginary: %f i", i, Pxr[i][0], Pxr[i][1]);
            end
    end

assign x = {4'h0,4'h0,4'h0,4'h0,4'h0,4'h0,4'h0,4'h0,4'h8,4'h7,4'h6,4'h5,4'h4,4'h3,4'h2,4'h1};

endmodule
