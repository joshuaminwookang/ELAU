
function automatic integer log2floor;
    input integer n;
    integer m;
    integer p;
    begin
        m = -1;
        p = 1;
        while (p <= n) begin
            m = m + 1;
            p = p * 2;
        end
        log2floor = m;
    end
endfunction


// Copyright 2024 ETH Zurich and University of Bologna.
// Solderpad Hardware License, Version 0.51, see LICENSE for details.
// SPDX-License-Identifier: SHL-0.51
//
// Based on the work by Reto Zimmermann 1998 - ETH Zürich
// Originally written in VHDL, available under: 
// https://iis-people.ee.ethz.ch/~zimmi/arith_lib.html#library
//
// Authors:
// - Thomas Benz <tbenz@iis.ee.ethz.ch>
// - Philippe Sauter <phsauter@iis.ee.ethz.ch>
// - Paul Scheffler <paulsc@iis.ee.ethz.ch>
//
// Description :
// Multiplier for unsigned numbers with one input operand and the
// result in carry-save number representation (Brown). First adds two
// numbers, then multiplies the result with the multiplicand without
// performing final addition. Result is only valid if sum of
// carry-save input operands does not overflow.
// The exact bit-pattern of S and C change depending on the structure 
// of the implemented tree because operands are commutative and can be added in
// different parts of the tree, contributing either to S or C. 
// The only guarantee is that S+C is the same.

module MulCsvUns #(
	parameter int widthX = 8,     // word width of XS, XC (<= widthY)
	parameter int widthY = 8,     // word width of Y
	parameter int speed = 2  // performance parameter
) (
	input logic [widthX-1:0] XS,  // multiplier
	input logic [widthX-1:0] XC,  // multiplier
	input logic [widthY-1:0] Y,  // multiplicand
	output logic [widthX+widthY-1:0] PS,  // sum
	output logic [widthX+widthY-1:0] PC  // carry
);

	logic [widthX*(widthX+widthY)-1:0] PP;  // partial products
	logic [widthX+widthY-1:0] ST, CT;  // intermediate sum/carry bits

	// generation of partial products
	AddMulPPGenUns #(widthX, widthY) ppGen (
		.XS(XS),
		.XC(XC),
		.Y (Y),
		.PP(PP)
	);

	// carry-save addition of partial products
	AddMopCsv #(widthX+widthY, widthX, speed) csvAdd (
		.A (PP),
		.S (PS),
		.C (PC)
	);

endmodule


module Cpr #(
	parameter int              depth = 4,             // number of input bits
	parameter int speed = 2  // performance parameter
) (
	input  logic [depth-1:0] A,   // input bits
	input  logic [depth-4:0] CI,  // intermediate carries in
	output logic             S,   // sum out
	output logic             C,   // carry out
	output logic [depth-4:0] CO   // intermediate carries out
);

	logic [depth+2*(depth-2)-1:0] F;  // FIFO vector of internal signals
	logic [            depth-3:0] CIT, COT;  // temp. int. carries

	// put input bits to beginning of FIFO vector
	assign F[depth-1:0] = A;

	// temporary intermediate carries in
	assign CIT[depth-3:0] = {1'b0, CI[depth-4:0]};

	// compressor with linear structure
	if (speed == 0) begin : slowCpr
		// first full-adder
		FullAdder fa0 (
			.A (F[0]),
			.B (F[1]),
			.CI(F[2]),
			.S (F[depth]),
			.CO(COT[0])
		);

		// linear arrangement of full-adders
		for (genvar i = 1; i < depth-2; i++) begin : linear
			FullAdder fa (
				.A (F[i+2]),
				.B (CIT[i-1]),
				.CI(F[depth+(i-1)*2]),
				.S (F[depth+i*2]),
				.CO(COT[i])
			);
		end
	end  // compressor with tree structure
	else begin : fastCpr
		// tree arrangement of full-adders

		for (genvar i = 0; i < depth-2; i++) begin : tree
			// take inputs from beginning of FIFO vector
			// attach sum output to end of FIFO vector
			// put carry output to intermediate carry-out
			FullAdder fa (
				.A (F[i*3]),
				.B (F[i*3+1]),
				.CI(F[i*3+2]),
				.S (F[depth+i*2]),
				.CO(COT[i])
			);
			// attach intermediate carry-in to end of FIFO vector
			assign F[depth+i*2+1] = CIT[i];
		end
	end

	// intermediate carries out
	assign CO = COT[depth-4:0];

	// sum and carry out
	assign S  = F[3*depth-6];
	assign C  = COT[depth-3];

endmodule

module AddMopCsv #(
	parameter int              width = 8,             // word width
	parameter int              depth = 4,             // number of operands
	parameter int speed = 2  // performance parameter
) (
	input  logic [(depth*width)-1:0] A,  // operands
	output logic [        width-1:0] S,  // sum
	output logic [        width-1:0] C   // carry vector
);

	logic [      (depth*width)-1:0] AT;  // re-arranged inputs
	logic [(depth-3)*(width+1)-1:0] CI;  // intermediate carries
	logic [              width-1:0] CT;  // unshifted output carries

	// re-arrange input bits: group bits of same magnitude
	always_comb begin : swizzle
		for (int k = 0; k < depth; k++) begin
			for (int i = 0; i < width; i++) begin
				AT[i*depth+k] = A[k*width+i];
			end
		end
	end

	// set intermediate carries into first slice
	assign CI[depth-4:0] = 1'b0;

	// carry-save addition using (m,2) compressor bit-slices
	for (genvar i = 0; i < width; i++) begin : bits
		Cpr #(
			.depth(depth),
			.speed(speed)
		) slice (
			.A (AT[(i+1)*depth -1 : i*depth]),
			.CI(CI[(i+1)*(depth-3) -1 : i*(depth-3)]),
			.S (S[i]),
			.C (CT[i]),
			.CO(CI[(i+2)*(depth-3) -1 : (i+1)*(depth-3)])
		);
	end

	// shift left output carries by one position
	assign C = {CT[width-2:0], 1'b0};

endmodule

module AddMulPPGenUns #(
	parameter int widthX = 8,  // word width of XS, XC
	parameter int widthY = 8,  // word width of Y
	localparam int widthP = widthX+widthY
) (
	input logic [widthX-1:0] XS,  // multipliers
	input logic [widthX-1:0] XC,
	input logic [widthY-1:0] Y,  // multiplicand
	output logic [widthX*(widthP)-1:0] PP  // partial products
);

	logic [widthX-1:0] M1, M2;  // recoded multiplier
	logic [widthY+1:0] YT;  // expanded Y

	// recode multiplier
	assign M1 = XS ^ XC;
	assign M2 = XS & XC;

	// expand Y (used for term 2y)
	assign YT = {1'b0, Y, 1'b0};

	logic [widthX*widthP-1:0] ppt;

	// partial product generation
	always_comb begin
		ppt = 0;
		for (int i = 0; i < widthX; i++) begin
			for (int k = 0; k <= widthY; k++) begin
				ppt[i*widthP+i+k] = (M1[i] & YT[k+1]) | (M2[i] & YT[k]);
			end
		end
		PP = ppt;
	end

endmodule
