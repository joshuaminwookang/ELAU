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
// Detection of the all-zeroes vector.
// Z = (A == '0)

module AllZeroDet #(
	parameter int width = 8  // word width
) (
	input  logic [width-1:0] A,  // operand
	output logic             Z   // all-zeroes flag
);

	logic ZT;  // temp.

	// all-zeroes detection
	RedOr #(width) zeroFlag (
		.A(A),
		.Z(ZT)
	);

	// result
	assign Z = ~ZT;

endmodule



module behavioural_AllZeroDet #(
	parameter int width = 8  // word width
) (
	input  logic [width-1:0] A,  // operand
	output logic             Z   // all-zeroes flag
);
	assign Z = (A == '0);
endmodule

module RedOr #(
	parameter int width = 8  // word width
) (
	input logic [width-1:0] A,  // input vector
	output logic Z  // output bit
);

	logic zv;
	
	// OR all bits
	// behavioral description used (well handled by all synthesizers)
	always_comb begin
		zv = A[0];
		for (int i = 1; i < width; i++) begin
			zv |= A[i];
		end
		Z = zv;
	end

endmodule
