
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
// Computes integer logarithm to base 2. 
// Z = floor(log A))
// Example: A = "00010110" -> Z = "100".

module Log2 #(
    parameter int width = 8,     // word width
    parameter int speed = 1  // performance parameter
) (
    input logic [width-1:0] A,  // operand
    output logic [$clog2(width)-1:0] Z  // result
);

    logic [width-1:0] ZT;  // temp.

    // leading zero detection (i.e. most significant '1')
    LeadZeroDet #(width, speed) loz (
        .A(A),
        .Z(ZT)
    );

    // binary encode
    Encode #(width) enc (
        .A(ZT),
        .Z(Z)
    );

endmodule



module behavioural_Log2 #(
    parameter int width = 8,     // word width
    parameter int speed = 1  // performance parameter
) (
    input logic [width-1:0] A,  // operand
    output logic [$clog2(width)-1:0] Z  // result
);
    always_comb begin
		Z = '0;
		for (int i = 0; i < width ; i++ ) begin
			if(A[i] == 1'b1) begin
				Z = i;
			end
		end
	end
endmodule

module Encode #(
	parameter int width = 8  // word width
) (
	input logic [width-1:0] A,  // input vector
	output logic [$clog2(width)-1:0] Z  // enc. output
);

	localparam int n = width;
	localparam int m = $clog2(width);

	logic zv;

	// example: n = 8, m = 3
	//   Z[0] = A[7] || A[5] || A[3] || A[1]
	//   Z[1] = A[7] || A[6] || A[3] || A[2]
	//   Z[2] = A[7] || A[6] || A[5] || A[4]
	// indices correspond to position of black nodes in Sklansky parallel-prefix
	// algorithm
	always_comb begin
		for (int l = 1; l <= m; l++) begin : outbit
			zv = 1'b0;
			for (int k = 0; k < 2**(m-l); k++) begin
				for (int i = 0; i < 2**(l-1); i++) begin
					if (k * 2**l + 2**(l-1) + i < n) begin
						zv |= A[k * 2**l + 2**(l-1) +i];
					end
				end
			end
			Z[l-1] = zv;
		end
	end

endmodule

module PrefixAnd #(
	parameter int width = 8,  // word width
	parameter int speed = 1  // performance parameter
) (
	input  logic [width-1:0] PI,  // propagate in
	output logic [width-1:0] PO   // propagate out
);

	localparam int n = width;  // prefix structure width
	localparam int m = $clog2(width);  // prefix structure depth

	// Sklansky parallel-prefix propagate-lookahead structure
	if (speed == 2) begin : fastPrefix
		logic [(m+1)*n-1:0] PT;

		assign PT[n-1:0] = PI;
		for (genvar l = 1; l <= m; l++) begin : levels
			for (genvar k = 0; k < 2**(m-l); k++) begin : groups
				for (genvar i = 0; i < 2**(l-1); i++) begin : bits
					if ((k * 2**l +i) < n) begin : white
						assign PT[l*n + k*2**l +i] = PT[(l-1)*n + k* 2**l +i];
					end
					if ((k * 2**l + 2**(l-1) +i) < n) begin : black
						assign PT[l*n + k*2**l + 2**(l-1) + i] =
								PT[(l-1)*n + k*2**l + 2**(l-1) + i] 
							  & PT[(l-1)*n + k*2**l + 2**(l-1) - 1];
					end
				end : bits
			end : groups
		end : levels
		assign PO = PT[(m+1)*n-1 : m*n];
	end // fastPrefix

	// Brent-Kung parallel-prefix propagate-lookahead structure
	if (speed == 1) begin : mediumPrefix
		logic [(2*m)*width-1:0] PT;

		assign PT[n-1:0] = PI;
		for (genvar l = 1; l <= m; l++) begin : levels1
			for (genvar k = 0; k < 2 ** (m - l); k++) begin : groups
				for (genvar i = 0; i < 2 ** l - 1; i++) begin : bits
					if ((k * 2 ** l + i) < n) begin : white
						assign PT[l*width+k*2**l+i] = PT[(l-1)*width+k*2**l+i];
					end
				end : bits
				if ((k * 2 ** l + 2 ** l - 1) < n) begin : black
					assign PT[l*width + k*2**l + 2**l - 1] =
		PT[(l-1)*width + k*2**l + 2**l - 1] & PT[(l-1)*width + k*2**l + 2**(l-1) - 1];
				end
			end : groups
		end : levels1
		for (genvar l = m + 1; l <= 2 * m - 1; l++) begin : levels2
			for (genvar i = 0; i < 2 ** (2 * m - l); i++) begin : bits
				if (i < n) begin : white
					assign PT[l*width+i] = PT[(l-1)*width+i];
				end
			end : bits
			for (genvar k = 1; k < 2 ** (l - m); k++) begin : groups
				if (l < 2 * m - 1) begin : empty
					for (genvar i = 0; i < 2 ** (2 * m - l - 1) - 1; i++) begin : bits
						if ((k * 2 ** (2 * m - l) + i) < n) begin : white
							assign PT[l*width+k*2**(2*m-l)+i] = PT[(l-1)*width+k*2**(2*m-l)+i];
						end
					end : bits
				end
				if ((k * 2 ** (2 * m - l) + 2 ** (2 * m - l - 1) - 1) < n) begin : black
					assign PT[l*width + k*2**(2*m-l) + 2**(2*m-l-1) - 1] =
									PT[(l-1)*width + k*2**(2*m-l) + 2**(2*m-l-1) - 1] 
								  & PT[(l-1)*width + k*2**(2*m-l) - 1];
				end
				for (genvar i = 2 ** (2 * m - l - 1); i < 2 ** (2 * m - l); i++) begin : bits
					if ((k * 2 ** (2 * m - l) + i) < n) begin : white
						assign PT[l*width+k*2**(2*m-l)+i] = PT[(l-1)*width+k*2**(2*m-l)+i];
					end
				end : bits
			end : groups
		end : levels2
		assign PO = PT[2*m*width-1 : (2*m-1)*width];
	end // mediumPrefix

	// Serial-prefix propagate-lookahead structure
	if (speed == 0) begin : slowPrefix
		logic [width-1:0] PT;

		assign PT[0] = PI[0];
		for (genvar i = 1; i < n; i++) begin : bits
			assign PT[i] = PI[i] & PT[i-1];
		end
		assign PO = PT;
	end // slowPrefix  

endmodule

module LeadZeroDet #(
	parameter int width = 8,     // word width
	parameter int speed = 1  // performance parameter
) (
	input  logic [width-1:0] A,  // operand
	output logic [width-1:0] Z   // LZD output
);

	logic [width-1:0] PI;  // prefix prop. in
	logic [width-1:0] PO;  // prefix prop. out
	logic [width-1:0] PIT;  // temp.
	logic [width-1:0] POT;  // temp.

	// calculate prefix propagate in signals
	assign PI = ~A;

	// reverse bit order of PI
	for (genvar i = width - 1; i >= 0; i--) begin
		assign PIT[i] = PI[width-i-1];
	end

	// solve reverse prefix problem for leading-zeroes detection
	// (example: "000101" -> "111100")
	PrefixAnd #(width, speed) prefix_and (
		.PI(PIT),
		.PO(POT)
	);

	// reverse bit order of PO
	for (genvar i = width - 1; i >= 0; i--) begin
		assign PO[i] = POT[width-i-1];
	end

	// code output: only bit indicating position of first '1' is '1'
	// (example: "000101" -> "000100")
	assign Z[width-1] = A[width-1];
	for (genvar i = width - 2; i >= 0; i--) begin
		assign Z[i] = PO[i+1] & A[i];
	end

endmodule
