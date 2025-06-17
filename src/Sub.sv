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
// Binary subtractor using parallel-prefix carry-lookahead logic.
// S = A-B

module Sub #(
	parameter width = 8,  // word width
	parameter lau_pkg::speed_e speed = lau_pkg::FAST  // performance parameter
) (
	input  logic [width-1:0] A, // operands
	input  logic [width-1:0] B,
	output logic [width-1:0] S  // sum
);

	logic [width-1:0] BI;  // B inverted
	logic [width-1:0] GI, PI;  // prefix gen./prop. in
	logic [width-1:0] GO, PO;  // prefix gen./prop. out
	logic [width-1:0] PT;  // adder propagate temp

	// invert B for subtraction
	assign BI = ~B;

	// calculate prefix input generate/propagate signal (0)
	assign GI[0] = A[0] | BI[0];
	assign PI[0] = 1'b0;
	// calculate adder propagate signal (0) (PT = A xor B)
	assign PT[0] = A[0] ^ BI[0];
	// calculate prefix input generate/propagate signals (1 to width-1)
	for (genvar i = 1; i < width; i++) begin : preproc
		assign GI[i] = A[i] & BI[i];
		assign PI[i] = A[i] | BI[i];
		// calculate adder propagate signal (1 to width-1) (PT = A xor B)
		assign PT[i] = ~GI[i] & PI[i];
	end

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

	// calculate sum bits
	assign S = PT ^ {GO[width-2:0], 1'b1};

endmodule



// module behavioural_Sub #(
// 	parameter width = 8,  // word width
// 	parameter lau_pkg::speed_e speed = lau_pkg::FAST  // performance parameter
// ) (
// 	input  logic [width-1:0] A, // operands
// 	input  logic [width-1:0] B,
// 	output logic [width-1:0] S  // sum
// );
// 	assign S = A - B;
// endmodule