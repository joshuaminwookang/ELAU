
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
// Three-operand carry-save adder using full-adders.

module AddCsv #(
	parameter int width = 8  // word width
) (
	input  logic [width-1:0] A1,  // operands
	input  logic [width-1:0] A2,
	input  logic [width-1:0] A3,
	output logic [width-1:0] S,   // sum / carry vector
	output logic [width-1:0] C
);

	logic [width-1:0] CT;  // unshifted output carries

	// carry-save addition using full-adders
	for (genvar i = 0; i < width; i++) begin : bits
		FullAdder fa (
		.A  (A1[i]),
		.B  (A2[i]),
		.CI (A3[i]),
		.S  (S[i] ),
		.CO (CT[i])
		);
	end

	// rotate output carries by one position
	assign C = {CT[width-2:0], 1'b0};

endmodule

module FullAdder (
	input  logic A,
	input  logic B,
	input  logic CI,  // operands
	output logic S,
	output logic CO  // sum and carry out
);

	logic [1:0] Auns, Buns, CIuns, Suns;  // unsigned temp

	// type conversion: std_logic -> 2-bit unsigned
	assign Auns  = {1'b0, A};
	assign Buns  = {1'b0, B};
	assign CIuns = {1'b0, CI};

	// should force the compiler to use a full-adder cell
	assign Suns = Auns + Buns + CIuns;

	// type conversion: 2-bit unsigned -> std_logic
	assign {CO, S} = Suns;

endmodule

// module behavioural_AddCsv #(
// 	parameter int width = 8  // word width
// ) (
// 	input  logic [width-1:0] A1,  // operands
// 	input  logic [width-1:0] A2,
// 	input  logic [width-1:0] A3,
// 	output logic [width-1:0] S,   // sum / carry vector
// 	output logic [width-1:0] C
// );
// 	assign S = A1 ^ A2 ^ A3; // carry-less sum
// 	assign C = ((A1 & A2) | (A1 & A3) | (A2 & A3)) << 1;
// endmodule
