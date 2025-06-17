
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
// Binary subtractor using parallel-prefix carry-lookahead logic with:
//   - carry-in (CI), subtracted
//   - carry-out (CO), '1' if subtraction result is negative
// {CO,S} = A - B -CI

module SubC #(
	parameter width = 8,  // word width
	parameter int speed = 0  // performance parameter
) (
	input logic [width-1:0] A,  // operands
	input logic [width-1:0] B,
	input logic CI,  // carry in (subtracted)
	output logic [width-1:0] S,  // sum
	output logic CO // carry out ('1' if S negative)
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
	assign PI[0] = 1'b0;
	// calculate adder propagate signal (0) (PT = A xor B)
	assign PT[0] = A[0] ^ BI[0];
	// calculate prefix input generate/propagate signals (1 to width-1)
	generate
		for (genvar i = 1; i < width; i++) begin : preproc
			assign GI[i] = A[i] & BI[i];
			assign PI[i] = A[i] | BI[i];
			// calculate adder propagate signal (1 to width-1) (PT = A xor B)
			assign PT[i] = ~GI[i] & PI[i];
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

	// calculate sum and carry-out bits
	assign S  = PT ^ {GO[width-2:0], CII};
	assign CO = ~GO[width-1];

endmodule



// module behavioural_SubC #(
// 	parameter width = 8,  // word width
// 	parameter int speed = 0  // performance parameter
// ) (
// 	input logic [width-1:0] A,  // operands
// 	input logic [width-1:0] B,
// 	input logic CI,  // carry in (subtracted)
// 	output logic [width-1:0] S,  // sum
// 	output logic CO // carry out ('1' if S negative)
// );
// 	assign {CO,S} = A - B -CI;
// endmodule

module PrefixAndOr #(
	parameter int              width = 8,             // word width
	parameter int speed = 0  // performance parameter
) (
	input  logic [width-1:0] GI,  // gen./prop. in
	input  logic [width-1:0] PI,  // gen./prop. in
	output logic [width-1:0] GO,  // gen./prop. out
	output logic [width-1:0] PO   // gen./prop. out
);

	// Constants
	localparam int n = width;  // prefix structure width
	localparam int m = $clog2(width);  // prefix structure depth

	// Sklansky parallel-prefix carry-lookahead structure
	if (speed == 2) begin : fastPrefix
		logic [(m+1)*n-1:0] GT, PT;  // gen./prop. temp
			assign GT[n-1:0] = GI;
			assign PT[n-1:0] = PI;
			for (genvar l = 1; l <= m; l++) begin : levels
				for (genvar k = 0; k < 2 ** (m - l); k++) begin : groups
					for (genvar i = 0; i < 2 ** (l - 1); i++) begin : bits
						// pass prop and gen to following nodes
						if ((k * 2 ** l + i) < n) begin : white
							assign GT[l*n+k*2**l+i] = GT[(l-1)*n+k*2**l+i];
							assign PT[l*n+k*2**l+i] = PT[(l-1)*n+k*2**l+i];
						end
						// calculate new propagate and generate
						if ((k * 2 ** l + 2 ** (l - 1) + i) < n) begin : black
							assign GT[l*n + k*2**l + 2**(l-1) + i] = 
												GT[(l-1)*n + k*2**l + 2**(l-1) + i]
											  | (  PT[(l-1)*n + k*2**l + 2**(l-1) + i]
												 & GT[(l-1)*n + k*2**l + 2**(l-1) - 1] );
							assign PT[l*n + k*2**l + 2**(l-1) + i] = 
													PT[(l-1)*n + k*2**l + 2**(l-1) + i]
												  & PT[(l-1)*n + k*2**l + 2**(l-1) - 1];
						end
					end
				end
			end
			assign GO = GT[(m+1)*n-1 : m*n];
			assign PO = PT[(m+1)*n-1 : m*n];
	end

	// Brent-Kung parallel-prefix carry-lookahead structure
	if (speed == 1) begin : mediumPrefix
		logic [(2*m)*n -1:0] GT, PT;  // gen./prop. temp
		assign GT[n-1:0] = GI;
		assign PT[n-1:0] = PI;

		for (genvar l = 1; l <= m; l++) begin : levels1
			for (genvar k = 0; k < 2**(m-l); k++) begin : groups
				for (genvar i = 0; i < 2**l -1; i++) begin : bits
					if ((k* 2**l +i) < n) begin : white
						assign GT[l*n + k* 2**l +i] = GT[(l-1)*n + k* 2**l +i];
						assign PT[l*n + k* 2**l +i] = PT[(l-1)*n + k* 2**l +i];
					end // white
				end // bits
				if ((k* 2**l + 2**l -1) < n) begin : black
					assign GT[l*n + k* 2**l + 2**l -1] =
								GT[(l-1)*n + k* 2**l + 2**l - 1] |
							  | (  PT[(l-1)*n + k* 2**l + 2**l     -1] 
							     & GT[(l-1)*n + k* 2**l + 2**(l-1) -1]);
					assign PT[l*n + k* 2**l + 2**l -1] =
								PT[(l-1)*n + k*2**l + 2**l     -1] 
							  & PT[(l-1)*n + k*2**l + 2**(l-1) -1];
				end // black
			end
		end // level1
		for (genvar l = m +1; l < 2*m; l++) begin : levels2
			for (genvar i = 0; i < 2**(2*m -l); i++) begin : bits
				if (i < n) begin : white
					assign GT[l*n +i] = GT[(l-1)*n +i];
					assign PT[l*n +i] = PT[(l-1)*n +i];
				end // white
			end // bits
			for (genvar k = 1; k < 2**(l-m); k++) begin : groups
				if (l < 2*m -1) begin : empty
					for (genvar i = 0; i < 2**(2*m -l -1) -1; i++) begin : bits
						if ((k* 2**(2*m -l) +i) < n) begin : white
							assign GT[l*n + k* 2**(2*m -l) +i] = GT[(l-1)*n + k* 2**(2*m -l) +i];
							assign PT[l*n + k* 2**(2*m -l) +i] = PT[(l-1)*n + k* 2**(2*m -l) +i];
						end // white
					end
				end // empty
				if ((k* 2**(2*m -l) + 2**(2*m -l -1) -1) < n) begin : black
					assign GT[l*n + k* 2**(2*m -l) + 2**(2*m -l-1) -1] = 
								GT[(l-1)*n + k* 2**(2*m-l) + 2**(2*m -l-1) -1]
							  | (  PT[(l-1)*n + k* 2**(2*m-l) + 2**(2*m-l-1) -1] 
								 & GT[(l-1)*n + k* 2**(2*m-l) -1] );
					assign PT[l*n + k* 2**(2*m -l) + 2**(2*m -l-1) -1] = 
								PT[(l-1)*n + k* 2**(2*m -l) + 2**(2*m -l-1) -1]
							  & PT[(l-1)*n + k* 2**(2*m -l) -1];
				end // black
				for (genvar i = 2**(2*m -l-1); i < 2**(2*m -l); i++) begin : bits
					if ((k* 2**(2*m -l) +i) < n) begin : white
						assign GT[l*n + k* 2**(2*m -l) +i] = GT[(l-1)*n + k* 2**(2*m -l) +i];
						assign PT[l*n + k* 2**(2*m -l) +i] = PT[(l-1)*n + k* 2**(2*m -l) +i];
					end // white
				end
			end
		end // level2
		assign GO = GT[2*m*n -1 : (2*m -1) * n];
		assign PO = PT[2*m*n -1 : (2*m -1) * n];
	end  // Serial-prefix carry-lookahead structure
	else if (speed == 0) begin : slowPrefix
		logic [n-1:0] GT, PT;  // gen./prop. temp
		assign GT[0] = GI[0];
		assign PT[0] = PI[0];
		
		for (genvar i = 1; i < n; i++) begin : bits
			assign GT[i] = GI[i] | (PI[i] & GT[i-1]);
			assign PT[i] = PI[i] & PT[i-1];
		end
		assign GO = GT;
		assign PO = PT;
	end

endmodule