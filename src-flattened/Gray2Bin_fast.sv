
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
// Converts a number from Gray to binary representation. Corresponds to a
// prefix problem in reversed bit order (compared to addition)

module Gray2Bin #(
	parameter int width = 8,  // word width
	parameter int speed = 2  // performance parameter
) (
	input  logic [width-1:0] G,  // Gray input
	output logic [width-1:0] B   // binary output
);

	logic [width-1:0] BT, GT;  // temp.

	// reverse bit order of G
	for (genvar i = 0; i < width; i++) begin : revG
		assign GT[i] = G[width-i-1];
	end


	// convert Gray to binary using prefix XOR computation
	PrefixXor #(
		.width(width),
		.speed(speed)
	) prefix_xor (
		.PI(GT),
		.PO(BT)
	);

	// reverse bit order of B

	for (genvar i = 0; i < width; i++) begin : revB
		assign B[i] = BT[width-1-i];
	end


endmodule



// module behavioural_Gray2Bin #(
// 	parameter int width = 8,  // word width
// 	parameter int speed = 2  // performance parameter
// ) (
// 	input  logic [width-1:0] G,  // Gray input
// 	output logic [width-1:0] B   // binary output
// );

//     for (genvar i = 0; i < width; i++)
//         assign B[i] = ^G[width-1:i];

// endmodule


module PrefixXor #(
	parameter int width = 8,  // word width
	parameter int speed = 2  // performance parameter
)  // performance parameter
(
	input  logic [width-1:0] PI,  // propagate in
	output logic [width-1:0] PO
);  // propagate out

	// prefix structure width and depth
	localparam n = width;
	localparam m = $clog2(width);

	// Sklansky parallel-prefix propagate-lookahead structure
	if (speed == 2) begin : fastPrefix
		logic [(m+1)*n-1:0] PT;

		// Initialize PT
		assign PT[n-1:0] = PI;

		// Compute PT levels
		for (genvar l = 1; l <= m; l++) begin : levels
			for (genvar k = 0; k < 2 ** (m - l); k++) begin : groups
				for (genvar i = 0; i < 2 ** (l - 1); i++) begin : bits_white
					if ((k * 2 ** l + i) < n) begin
						assign PT[l*n+k*2**l+i] = PT[(l-1)*n+k*2**l+i];
					end
				end

				for (genvar i = 0; i < 2 ** (l - 1); i++) begin : bits_black
					if ((k * 2 ** l + 2 ** (l - 1) + i) < n) begin
						assign PT[l*n + k*2**l + 2**(l-1) + i] =
				PT[(l-1)*n + k*2**l + 2**(l-1) + i] ^ PT[(l-1)*n + k*2**l + 2**(l-1) - 1];
					end
				end
			end
		end

		// Assign output
		assign PO = PT[(m+1)*n-1 : m*n];

	end  // Brent-Kung parallel-prefix propagate-lookahead structure
	else if (speed == 1) begin : mediumPrefix
		logic [2*m*n -1:0] PT;

		// Initialize PT
		assign PT[n-1:0] = PI;

		// Compute levels1
		for (genvar l = 1; l <= m; l++) begin : levels1
			for (genvar k = 0; k < 2**(m-l); k++) begin : groups
				for (genvar i = 0; i < 2**l -1; i++) begin : bits
					if ((k* 2**l +i) < n) begin : white
						assign PT[l*n + k* 2**l +i] = PT[(l-1)*n + k* 2**l +i];
					end // white
				end // bits
				if ((k * 2**l + 2**l -1) < n) begin : black
					assign PT[l*n + k* 2**l + 2**l - 1] =
								PT[(l-1)*n + k* 2**l + 2**l     -1] 
							  ^ PT[(l-1)*n + k* 2**l + 2**(l-1) -1];
				end // black
			end
		end // level1

		// Compute levels2
		for (genvar l = m + 1; l < 2*m; l++) begin : levels2
			for (genvar i = 0; i < 2**(2*m -l); i++) begin : bits
				if (i < n) begin : white
					assign PT[l*n +i] = PT[(l-1)*n +i];
				end // white
			end // bits

			for (genvar k = 1; k < 2**(l-m); k++) begin : groups
				if (l < 2*m -1) begin : empty
					for (genvar i = 0; i < 2**(2*m -l-1) -1; i++) begin : bits
						if ((k* 2**(2*m -l) +i) < n) begin : white
							assign PT[l*n + k* 2**(2*m -l) +i] = PT[(l-1)*n + k* 2**(2*m -l) +i];
						end // white
					end // bits
				end // empty

				if ((k* 2**(2*m -l) + 2**(2*m -l-1) -1) < n) begin : black
					assign PT[l*n + k* 2**(2*m -l) + 2**(2*m -l-1) -1] =
								PT[(l-1)*n + k* 2**(2*m-l) + 2**(2*m -l-1) -1] 
							  ^ PT[(l-1)*n + k* 2**(2*m-l) -1];
				end // black

				for (genvar i = 2**(2*m -l-1); i < 2**(2*m -l); i++) begin : bits
					if ((k* 2**(2*m -l) +i) < n) begin : white
						assign PT[l*n + k* 2**(2*m -l) +i] = PT[(l-1)*n + k* 2**(2*m -l) +i];
					end // white
				end // bits
			end
		end // levels2

		// Assign output
		assign PO = PT[2*m*n -1 : (2*m -1) * n];

	end  // Serial-prefix propagate-lookahead structure
	else begin : slowPrefix
		logic [n-1:0] PT;

		// Initialize PT
		assign PT[0] = PI[0];
		for (genvar i = 1; i < n; i++) begin
			assign PT[i] = PI[i] ^ PT[i-1];
		end

		// Assign output
		assign PO = PT;
	end

endmodule