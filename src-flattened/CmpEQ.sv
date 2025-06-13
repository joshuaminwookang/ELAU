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

package lau_pkg;


	// Computes floor(log2(n))
	// Input n should be greater than 0 for meaningful results
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

	typedef enum logic [1:0] {
		SLOW = 2'b00,
		MEDIUM = 2'b01,
		FAST = 2'b10
	} speed_e;

endpackage


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
// Equality comparison of two numbers.
// EQ = (A == B)

module CmpEQ #(
	parameter int width = 8  // word width
) (
	input  logic [width-1:0] A,  // operands
	input  logic [width-1:0] B,
	output logic             EQ  // equal flag
);

	logic [width-1:0] EQT;  // temp.

	// comparison of individual bits
	for (genvar i = 0; i < width; i++) begin : cmpBit
		assign EQT[i] = ~(A[i] ^ B[i]);
	end

	// AND all comparison bits
	assign EQ = &EQT;

endmodule



module behavioural_CmpEQ #(
	parameter int width = 8  // word width
) (
	input  logic [width-1:0] A,  // operands
	input  logic [width-1:0] B,
	output logic             EQ  // equal flag
);
	assign EQ = (A == B);
endmodule
