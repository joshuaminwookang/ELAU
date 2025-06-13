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
// Prefix structures of different depth (i.e. speed) for propagate calculation
// in different arithmetic units. Compute in m levels new propagate signal
// pairs for always larger groups of bits. Basic logic operation: OR for
// propagate signals.

module PrefixOr #(
	parameter int width = 8,  // word width
	parameter lau_pkg::speed_e speed = lau_pkg::FAST  // performance parameter
) (
	input  logic [width-1:0] PI,  // propagate in
	output logic [width-1:0] PO   // propagate out
);

	localparam int n = width;  // prefix structure width
	localparam int m = $clog2(width);  // prefix structure depth

	// Sklansky parallel-prefix propagate-lookahead structure
	if (speed == lau_pkg::FAST) begin : fastPrefix
		logic [(m+1)*n-1:0] PT;  // gen./prop. temp
		assign PT[n-1:0] = PI;

		for (genvar l = 1; l <= m; l++) begin : levels
			for (genvar k = 0; k < 2 ** (m - l); k++) begin : groups
				for (genvar i = 0; i < 2 ** (l - 1); i++) begin : bits
					if ((k * 2 ** l + i) < n) begin : white
						assign PT[l*n+k*2**l+i] = PT[(l-1)*n+k*2**l+i];
					end
					if ((k * 2 ** l + 2 ** (l - 1) + i) < n) begin : black
						assign PT[l*n + k*2**l + 2**(l-1) + i] =
			  PT[(l-1)*n + k*2**l + 2**(l-1) + i] |
			  PT[(l-1)*n + k*2**l + 2**(l-1) - 1];
					end
				end
			end
		end
		assign PO = PT[(m+1)*n-1 : m*n];
	end  // Brent-Kung parallel-prefix propagate-lookahead structure
	else if (speed == lau_pkg::MEDIUM) begin : mediumPrefix
		logic [(2*m)*n -1:0] PT;  // gen./prop. temp
		assign PT[n-1:0] = PI;

		for (genvar l = 1; l <= m; l++) begin : levels1
			for (genvar k = 0; k < 2**(m-l); k++) begin : groups
				for (genvar i = 0; i < 2**l -1; i++) begin : bits
					if ((k * 2**l +i) < n) begin : white
						assign PT[l*n + k* 2**l +i] = PT[(l-1)*n + k* 2**l +i];
					end
				end // bits
				if ((k*2**l + 2**l -1) < n) begin : black
					assign PT[l*n + k*2**l + 2**l -1] =
									PT[(l-1)*n + k*2**l + 2**l     -1]
									| PT[(l-1)*n + k*2**l + 2**(l-1) -1];
				end
			end
		end // level1

		for (genvar l = m +1; l < 2*m; l++) begin : levels2
			for (genvar i = 0; i < 2**(2*m -l); i++) begin : bits
				if (i < n) begin : white
					assign PT[l*n +i] = PT[(l-1)*n +i];
				end
			end // bits
			for (genvar k = 1; k < 2**(l-m); k++) begin : groups
				if (l < 2*m -1) begin : empty
					for (genvar i = 0; i < 2**(2*m -l -1) -1; i++) begin : bits
						if ((k * 2**(2*m -l) + i) < n) begin : white
							assign PT[l*n + k* 2**(2*m-l) +i] = PT[(l-1)*n + k* 2**(2*m-l) +i];
						end
					end
				end // empty
				if ((k* 2**(2*m -l) + 2**(2*m -l -1) -1) < n) begin : black
					assign PT[l*n + k* 2**(2*m -l) + 2**(2*m -l-1) -1] =
									PT[(l-1)*n + k* 2**(2*m -l) + 2**(2*m -l-1) -1]
								  | PT[(l-1)*n + k* 2**(2*m -l) - 1];
				end // black
				for (genvar i = 2**(2*m -l -1); i < 2**(2*m -l); i++) begin : bits
					if ((k* 2**(2*m -l) +i) < n) begin : white
						assign PT[l*n + k* 2**(2*m -l) +i] = PT[(l-1)*n + k* 2**(2*m -l) +i];
					end
				end // bits
			end
		end // level2
		assign PO = PT[2*m*n -1 : (2*m -1) * n];
	end  // Serial-prefix propagate-lookahead structure
	else if (speed == lau_pkg::SLOW) begin : slowPrefix
		logic [n-1:0] PT;  // gen./prop. temp
		assign PT[0] = PI[0];

		for (genvar i = 1; i < n; i++) begin : bits
			assign PT[i] = PI[i] | PT[i-1];
		end
		assign PO = PT;
	end

endmodule
