# DFT-in-SystemVerilog

Discrete Fourier transform for Fast Multiplication in SystemVerilog - May 2016

Submitted as a special project assignment as part of CMPE 297 Special Topics -- Cryptographic Application-Specific Hardware 

CMPE 297 Special Topics Course on: Application-Specific Hardware Design for computer arithmetic as applied in cryptography, coding, and digital signal processing. Micro-code driven state-machine and its application in cryptography.


Team: Saar Sagir and Samira C. Oliva Madrigal

Lab computes DFT on 32-bit integers using matrix form

Easily convert to FFT by applying butterfly technique which halves the number of multiply accumulates

Basis for NTT (Numer Theoretic Transform): DFT applied to rings and fields, becomes fast if transform length is power of two or may opt for Wang/Zhu algorithms.

Steps

1 Convert integers to polynomial form and zero-extend to length of product--> A`, B`

2 Covert DFT( A`), --> Ax DFT(B`) --> Bx

3 Multiply in DFT domain (Ax, Bx) --> Cx 

4 Interpolation (Cx) computes inverser DFT --> C` in coefficient representation (real numbers)

5 Flatten (C`) to decimal representation as a single number (normal component by component multiplication)

We use a lookup up table with precomputed values Vandermonde matrix, for the Nth roots of unity, from Excel sheet.

May use Excel sheet or reference calculator. 

For our example, result had  3.45 error, (expected - actual) / expected, based on FFT/DFT calculator. 

References: 

[1] https://inst.eecs.berkeley.edu/~cs170/fa07/handouts/cs170FFT.pdf

[2] Algorithms by Dasgupta, Papadimitriou & Vazirani

[3] https://www.youtube.com/watch?v=D5ueRUyCP58

[4] http://www.pi314.net/eng/multipli_fft.php

[5] http://calculator.vhex.net/calculator/fast-fourier-transform-calculator-fft/1d-discrete-fourier-transform

[6] http://scistatcalc.blogspot.com/2013/12/fft-calculator.html
