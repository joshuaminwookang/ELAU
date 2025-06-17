
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
// Multiplier-adder for unsigned numbers (Brown). First multiplies two numbers,
// then adds an additional operand to the result.
// P = (X*Y) +A

module MulAddUns #(
	parameter int              widthX = 8,             // word width of XS, XC (<= widthY)
	parameter int              widthY = 8,             // word width of Y
	parameter int              widthA = 20,            // word width of A (>= widthX+widthY)
	parameter int speed  = 1  // performance parameter
) (
	input  logic [widthX-1:0] X,  // multiplier
	input  logic [widthY-1:0] Y,  // multiplicand
	input  logic [widthA-1:0] A,  // augend
	output logic [widthA-1:0] P   // product
);

	logic [widthX*(widthX+widthY)-1:0] PP;  // partial products
	logic [widthA-1:0] ST1, CT1, ST2, CT2;  // intermediate sum/carry bits

	// generation of partial products
	MulPPGenUns #(widthX, widthY) ppGen (
		.X (X),
		.Y (Y),
		.PP(PP)
	);

	// carry-save addition of partial products
	AddMopCsv #(widthX + widthY, widthX, speed) csvAdd1 (
		.A (PP),
		.S (ST1[widthX+widthY-1:0]),
		.C (CT1[widthX+widthY-1:0])
	);

	// extend
	for (genvar i = widthA - 1; i >= widthX + widthY; i--) begin
		assign ST1[i] = 1'b0;
		assign CT1[i] = 1'b0;
	end

	// carry-save addition of augend
	AddCsv #(widthA) csvAdd2 (
		.A1(A),
		.A2(CT1),
		.A3(ST1),
		.S (ST2),
		.C (CT2)
	);

	// final carry-propagate addition
	Add #(widthA, speed) cpAdd (
		.A(ST2),
		.B(CT2),
		.S(P)
	);

endmodule



module behavioural_MulAddUns #(
	parameter int              widthX = 8,             // word width of XS, XC (<= widthY)
	parameter int              widthY = 8,             // word width of Y
	parameter int              widthA = 20,            // word width of A (>= widthX+widthY)
	parameter int speed  = 1  // performance parameter
) (
	input  logic [widthX-1:0] X,  // multiplier
	input  logic [widthY-1:0] Y,  // multiplicand
	input  logic [widthA-1:0] A,  // augend
	output logic [widthA-1:0] P   // product
);
	assign P = (X * Y) + A;
endmodule

module PrefixAndOr #(
	parameter int              width = 8,             // word width
	parameter int speed = 1  // performance parameter
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

module Add #(
	parameter int              width = 8,             // word width
	parameter int speed = 1  // performance parameter
) (
	input  logic [width-1:0] A,  // operands
	input  logic [width-1:0] B,
	output logic [width-1:0] S   // sum
);

	// Function: Binary adder using parallel-prefix carry-lookahead logic.

	logic [width-1:0] GI, PI;  // prefix gen./prop. in
	logic [width-1:0] GO, PO;  // prefix gen./prop. out
	logic [width-1:0] PT;  // adder propagate temp

	// Internal signals for unsigned operands
	logic [width-1:0] Auns, Buns, Suns;

	// default ripple-carry adder as slow implementation
	if (speed == 0) begin
		// type conversion: std_logic_vector -> unsigned
		assign Auns = A;
		assign Buns = B;

		// addition
		assign Suns = Auns + Buns;

		// type conversion: unsigned -> std_logic_vector
		assign S = Suns;
	end else begin
		// parallel-prefix adders as medium and fast implementations

		// calculate prefix input generate/propagate signals
		assign GI = A & B;
		assign PI = A | B;
		// calculate adder propagate signals (PT = A xor B)
		assign PT = ~GI & PI;

		// calculate prefix output generate/propagate signals
		PrefixAndOr #(
			.width(width),
			.speed(speed)
		) prefix (
			.GI(GI),
			.PI(PI),
			.GO(GO),
			.PO(PO)
		);

		// calculate sum bits
		assign S = PT ^ {GO[width-2:0], 1'b0};
	end
endmodule

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

module Cpr #(
	parameter int              depth = 4,             // number of input bits
	parameter int speed = 1  // performance parameter
) (
	input  logic [depth-1:0] A,   // input bits
	input  logic [depth-4:0] CI,  // intermediate carries in
	output logic             S,   // sum out
	output logic             C,   // carry out
	output logic [depth-4:0] CO   // intermediate carries out
);

	logic [depth+2*(depth-2)-1:0] F;  // FIFO vector of internal signals
	logic [            depth-3:0] CIT, COT;  // temp. int. carries

	// put input bits to beginning of FIFO vector
	assign F[depth-1:0] = A;

	// temporary intermediate carries in
	assign CIT[depth-3:0] = {1'b0, CI[depth-4:0]};

	// compressor with linear structure
	if (speed == 0) begin : slowCpr
		// first full-adder
		FullAdder fa0 (
			.A (F[0]),
			.B (F[1]),
			.CI(F[2]),
			.S (F[depth]),
			.CO(COT[0])
		);

		// linear arrangement of full-adders
		for (genvar i = 1; i < depth-2; i++) begin : linear
			FullAdder fa (
				.A (F[i+2]),
				.B (CIT[i-1]),
				.CI(F[depth+(i-1)*2]),
				.S (F[depth+i*2]),
				.CO(COT[i])
			);
		end
	end  // compressor with tree structure
	else begin : fastCpr
		// tree arrangement of full-adders

		for (genvar i = 0; i < depth-2; i++) begin : tree
			// take inputs from beginning of FIFO vector
			// attach sum output to end of FIFO vector
			// put carry output to intermediate carry-out
			FullAdder fa (
				.A (F[i*3]),
				.B (F[i*3+1]),
				.CI(F[i*3+2]),
				.S (F[depth+i*2]),
				.CO(COT[i])
			);
			// attach intermediate carry-in to end of FIFO vector
			assign F[depth+i*2+1] = CIT[i];
		end
	end

	// intermediate carries out
	assign CO = COT[depth-4:0];

	// sum and carry out
	assign S  = F[3*depth-6];
	assign C  = COT[depth-3];

endmodule

module AddMopCsv #(
	parameter int              width = 8,             // word width
	parameter int              depth = 4,             // number of operands
	parameter int speed = 1  // performance parameter
) (
	input  logic [(depth*width)-1:0] A,  // operands
	output logic [        width-1:0] S,  // sum
	output logic [        width-1:0] C   // carry vector
);

	logic [      (depth*width)-1:0] AT;  // re-arranged inputs
	logic [(depth-3)*(width+1)-1:0] CI;  // intermediate carries
	logic [              width-1:0] CT;  // unshifted output carries

	// re-arrange input bits: group bits of same magnitude
	always_comb begin : swizzle
		for (int k = 0; k < depth; k++) begin
			for (int i = 0; i < width; i++) begin
				AT[i*depth+k] = A[k*width+i];
			end
		end
	end

	// set intermediate carries into first slice
	assign CI[depth-4:0] = 1'b0;

	// carry-save addition using (m,2) compressor bit-slices
	for (genvar i = 0; i < width; i++) begin : bits
		Cpr #(
			.depth(depth),
			.speed(speed)
		) slice (
			.A (AT[(i+1)*depth -1 : i*depth]),
			.CI(CI[(i+1)*(depth-3) -1 : i*(depth-3)]),
			.S (S[i]),
			.C (CT[i]),
			.CO(CI[(i+2)*(depth-3) -1 : (i+1)*(depth-3)])
		);
	end

	// shift left output carries by one position
	assign C = {CT[width-2:0], 1'b0};

endmodule

module MulPPGenUns #(
	parameter widthX = 8,  // word width of X
	parameter widthY = 8
)  // word width of Y
(
	input logic [widthX-1:0] X,  // multiplier
	input logic [widthY-1:0] Y,  // multiplicand
	output logic [widthX*(widthX+widthY)-1:0] PP
);

	// width of single part. prod.
	localparam widthP = widthX + widthY;

	logic [widthX*widthP-1:0] ppt;  // partial products

	always_comb begin : ppGen
		ppt = '0;  // partial products

		// partial products x(i)y(k)
		for (int i = 0; i < widthX; i++) begin
			for (int k = 0; k < widthY; k++) begin
				ppt[i*widthP+i+k] = X[i] & Y[k];
			end
		end

		PP = ppt;
	end

endmodule
