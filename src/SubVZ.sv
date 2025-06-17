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
// Binary subtractor using parallel-prefix carry-lookahead logic with:
//   - carry-in (CI), subtracted
//   - 2's complement overflow flag (V)
//   - zero flag (only valid for CI = 0)
// S = A -B -CI
// V = (S>MAX) | (S<MIN)
// Z = (S==0)

module SubVZ #(
	parameter width = 8,  // word width
	parameter lau_pkg::speed_e speed = lau_pkg::FAST  // performance parameter
) (
	input logic [width-1:0] A,  // operands
	input logic [width-1:0] B,
	input logic CI,  // carry in (subtracted)
	output logic [width-1:0] S,  // sum
	output logic V, // overflow flag
	output logic Z  // zero flag (only valid for CI = 0)
);

	logic [width-1:0] BI;  // B inverted
	logic CII;  // CI inverted
	logic [width-1:0] GI, PI;  // prefix gen./prop. in
	logic [width-1:0] GO, PO;  // prefix gen./prop. out
	logic [width-1:0] PT;  // adder propagate temp

	// invert B and CI for subtraction
	assign BI = ~B;
	assign CII = ~CI;

	// calculate prefix input generate/propagate signal (0)
	assign GI[0] = (A[0] & BI[0]) | (A[0] & CII) | (BI[0] & CII);
	assign PI[0] = A[0] ^ BI[0];
	// calculate adder propagate signal (0) (PT = A xor B)
	assign PT[0] = PI[0];
	// calculate prefix input generate/propagate signals (1 to width-1)
	generate
		for (genvar i = 1; i < width; i++) begin : preproc
			assign GI[i] = A[i] & BI[i];
			assign PI[i] = A[i] ^ BI[i];
			// calculate adder propagate signal (1 to width-1) (PT = A xor B)
			assign PT[i] = PI[i];
		end
	endgenerate

	// calculate prefix output generate/propagate signals
	PrefixAndOr #(
		.width(width),
		.speed(speed)
	) prefix_inst (
		.GI(GI),
		.PI(PI),
		.GO(GO),
		.PO(PO)
	);

	// calculate sum and overflow bits
	assign S = PT ^ {GO[width-2:0], CII};
	assign V = GO[width-1] ^ GO[width-2];
	assign Z = PO[width-1];

endmodule



// module behavioural_SubVZ #(
// 	parameter width = 8,  // word width
// 	parameter lau_pkg::speed_e speed = lau_pkg::FAST  // performance parameter
// ) (
// 	input logic [width-1:0] A,  // operands
// 	input logic [width-1:0] B,
// 	input logic CI,  // carry in (subtracted)
// 	output logic [width-1:0] S,  // sum
// 	output logic V, // overflow flag
// 	output logic Z  // zero flag (only valid for CI = 0)
// );
// 	localparam MAX = 2**(width-1) -1;
// 	localparam MIN = - 2**(width-1);
// 	logic signed [width+1:0] Sext, Aext, Bext, CIext;
// 	assign Aext = signed'(A);
// 	assign Bext = signed'(B);
// 	assign CIext = {{width-1{1'b0}},CI};

// 	assign Sext = (Aext - Bext -CIext);
// 	assign S = Sext[width-1:0];
// 	assign V = (Sext>MAX) | (Sext<MIN);
// 	assign Z = CI ? 'x : (Sext == '0);
// endmodule
