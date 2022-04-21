`timescale 1ns / 1ps

module topModule_tb;

reg [3:0]A[7:0];
reg [3:0]B[7:0];
reg [63:0]C;
real Cr;
//real Cr[15:0][1:0];

integer i;

topModule top(A, B, C);

initial
    begin
        #20;
        Cr = $bitstoreal(C);
        #10;
        $display("Final Result is %f", Cr);
    //        for(i = 0 ; i < 16; i = i + 1)
    //            begin
    //                Cr[i][0] = $bitstoreal(C[i][0][63:0]);
    //                Cr[i][1] = $bitstoreal(C[i][1][63:0]);
    //                $display("# %d - In Complex Form: Real: %f, Imaginary: %f i", i, Cr[i][0], Cr[i][1]);
    //            end
    end

assign A = {4'h0,4'h0,4'h0,4'h0,4'h0,4'h2,4'h3,4'h7};
assign B = {4'h0,4'h0,4'h0,4'h0,4'h0,4'h4,4'h5,4'h1};

endmodule
