
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
// Detection of the all-ones vector.
// Z = (A == {width{1'b1}})

module AllOneDet #(
	parameter int width = 8  // word width
) (
	input  logic [width-1:0] A,  // operand
	output logic             Z   // all-ones flag
);
	// all-ones detection
	RedAnd #(width) zeroFlag (
		.A(A),
		.Z(Z)
	);

endmodule



module behavioural_AllOneDet #(
	parameter int width = 8  // word width
) (
	input  logic [width-1:0] A,  // operand
	output logic             Z   // all-ones flag
);
	assign Z = (A == {width{1'b1}});
endmodule

module RedAnd #(
	parameter int width = 8  // word width
) (
	input  logic [width-1:0] A,  // input vector
	output logic             Z   // output bit
);

	logic zv;

	// AND all bits
	// Using a procedural block for behavioral description
	always_comb begin
		zv = A[0];
		for (int i = 1; i < width; i++) begin
			zv &= A[i];
		end
		Z = zv;
	end

endmodule