
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
// (m,k)-counter adds m bits. Result is sum vector of k bits.
// Composed of (m,k)-counter slices. Condition: depth > 1.
// S = A[0]+A[1]+A[2]...+A[depth-1]

module Cnt #(
	parameter int              depth = 18,            // number of input bits
	parameter int speed = 2  // performance parameter
) (
	input logic [depth-1:0] A,  // input bits
	output logic [log2floor(depth):0] S  // sum output
);

	localparam m = depth;  // number of input bits
	localparam n = log2floor(depth) +1;  // number of output bits
	logic [m*n:0] CT;  // intermediate carries

	// input bits are first intermediate carries
	assign CT[m-1:0] = A;

	// linear arrangement of (m,k)-counter slices
	for (genvar i = 0; i < n-2; i++) begin : bits
		CntSlice #(m/(2**i), speed) slice (
			.A (CT[i*m + m/(2**i) -1 : i*m]),
			.S (S[i]),
			.CO(CT[(i+1)*m + m/(2**(i+1)) -1 : (i+1)*m])
		);
	end

	// add third carry if only two exist
	if (m/ 2**(n-2) == 2) begin : even
		assign CT[(n-2)*m+2] = 1'b0;
	end

	// full-adder for adding the last three carries
	FullAdder fa0 (
		.A (CT[(n-2)*m]),
		.B (CT[(n-2)*m+1]),
		.CI(CT[(n-2)*m+2]),
		.S (S[n-2]),
		.CO(S[n-1])
	);

endmodule



// module behavioural_Cnt #(
// 	parameter int              depth = 18,            // number of input bits
// 	parameter int speed = 2  // performance parameter
// ) (
// 	input logic [depth-1:0] A,  // input bits
// 	output logic [log2floor(depth):0] S  // sum output
// );
// 	always_comb begin
// 		S = '0;
// 		for(int i = 0; i < depth; i++) begin
// 			S += A[i];
// 		end
// 	end
// endmodule

module CntSlice #(
	parameter int              depth = 4,             // number of input bits
	parameter int speed = 2  // performance parameter
) (
	input  logic [  depth-1:0] A,  // input bits
	output logic               S,  // sum out
	output logic [depth/2-1:0] CO  // carries out
);

	localparam noFA = depth / 2;  // number of used full-adders
	localparam depthOdd = depth + ((depth + 1) % 2);  // next higher odd
	logic [3*noFA:0] F;  // FIFO vector of int. signals

	// put input bits to beginning of FIFO vector
	assign F[depth-1:0] = A;
	// add a zero if even number of input bits
	if (depth < depthOdd) begin
		assign F[depthOdd-1] = 1'b0;
	end

	// counter with linear structure
	if (speed == 0) begin : slowCnt
		// first full-adder
		FullAdder fa0 (
			.A (F[0]),
			.B (F[1]),
			.CI(F[2]),
			.S (F[depthOdd]),
			.CO(CO[0])
		);

		// linear arrangement of full-adders
		for (genvar i = 1; i < noFA; i++) begin : linear
			FullAdder fa (
				.A (F[i*2+1]),
				.B (F[i*2+2]),
				.CI(F[depthOdd+i-1]),
				.S (F[depthOdd+i]),
				.CO(CO[i])
			);
		end
	end  // counter with tree structure
	else begin : fastCnt
		// tree arrangement of full-adders
		for (genvar i = 0; i < noFA; i++) begin : tree
			FullAdder fa (
				.A (F[i*3]),
				.B (F[i*3+1]),
				.CI(F[i*3+2]),
				.S (F[depthOdd+i]),
				.CO(CO[i])
			);
		end
	end

	// sum out
	assign S = F[3*noFA];

endmodule