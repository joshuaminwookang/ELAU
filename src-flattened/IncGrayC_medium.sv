
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
// Incrementer for Gray numbers using parallel-prefix propagate-lookahead
// logic with:
//   - carry-in (CI)
// Bases on the following algorithm:
//   P = A(n-1) xor A(n-2) xor ... xor A(0)
//   Z(0) = A(0) xnor P
//   Z(i) = A(i) xor (A(i-1) and not A(i-2) and not A(i-3) ...
//                           and not A(0) and P)              ; i = 1, ..., n-2
//   Z(n-1) = A(n-1) xor (not A(n-3) and not A(n-4) ... and not A(0) and P)

module IncGrayC #(
	parameter int width = 16,  // word width
	parameter int speed = 1  // performance parameter
) (
	input  logic [width-1:0] A,   // operand
	input  logic             CI,  // carry in
	output logic [width-1:0] Z    // result
);

	logic [width-2:0] PI;  // prefix prop. in
	logic [width-2:0] PO;  // prefix prop. out
	logic P, PT;  // parity bit
	logic [width-1:2] T1;  // temp.
	logic [width-1:0] T2;  // temp.

	// calculate parity bit P
	RedXor #(width) parity (
		.A(A),
		.Z(P)
	);

	// calculate prefix input propagate signal (PI = not A(i))
	assign PI[width-2:1] = ~A[width-3:0];

	// feed slow P signal into prefix circuit for slow architecture
	// consider CI here for bits 2, ..., n-1
	if (speed == 0) begin
		assign PI[0] = P & CI;
	end else begin
		assign PI[0] = CI;
	end

	// calculate prefix output prop. signal (PO = not A(i-2) ... and not A(0))
	PrefixAnd #(width-1, speed) prefix_and (
		.PI(PI),
		.PO(PO)
	);

	// calculate (A and PO)
	assign T1[width-2:2] = A[width-3:1] & PO[width-3:1];
	// special case
	assign T1[width-1] = PO[width-2];

	// calculate (T1 and P) in fast architecture
	for (genvar i = width-1; i >= 2; i--) begin
		if (speed == 0) begin
			assign T2[i] = T1[i];
		end else begin
			assign T2[i] = T1[i] & P;
		end
	end

	// special cases, consider CI hier for bits 0 and 1
	assign T2[1] = A[0] & P & CI;
	assign T2[0] = ~(P | ~CI);

	// calculate result bits
	assign Z = A ^ T2;

endmodule



module behavioural_IncGrayC #(
	parameter int width = 16,  // word width
	parameter int speed = 1  // performance parameter
) (
	input  logic [width-1:0] A,   // operand
	input  logic             CI,  // carry in
	output logic [width-1:0] Z    // result
);
	logic [width-1:0] Abin, Zbin;
	behavioural_Gray2Bin #(width, speed) i_gray2bin (
		.G(A),
		.B(Abin)
	);
	
	assign Zbin = Abin + CI;

	behavioural_Bin2Gray #(width) i_bin2gray (
		.B(Zbin),
		.G(Z)
	);

endmodule


module RedXor #(
	parameter int width = 8  // word width
) (
	input logic [width-1:0] A,  // input vector
	output logic Z  // output bit
);

	logic zv;

	// XOR all bits
	// behavioral description used (well handled by all synthesizers)
	always_comb begin
		zv = A[0];
		for (int i = 1; i < width; i++) begin
			zv ^= A[i];
		end
		Z = zv;
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

module behavioural_Bin2Gray #(
	parameter int width = 8  // word width
) (
	input  logic [width-1:0] B,  // binary input
	output logic [width-1:0] G   // Gray output
);

	assign G = B ^ (B >> 1);

endmodule

module behavioural_Gray2Bin #(
	parameter int width = 8,  // word width
	parameter int speed = 1  // performance parameter
) (
	input  logic [width-1:0] G,  // Gray input
	output logic [width-1:0] B   // binary output
);

    for (genvar i = 0; i < width; i++)
        assign B[i] = ^G[width-1:i];

endmodule