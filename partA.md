# DIGITAL INTEGRATED CIRCUITS DESIGN LAB PC459EC

## PART-A Experiments: Write the code using Verilog HDL and simulate the following:

### 1.	Write structural and dataflow Verilog HDL models for

#### a) 4-bit ripple carry adder.

```
//Design of 4-bit ripple carry adder using 4, 1-bit Full adders using Structural modelling

module full_adder(a,b,cin,sum,co);
input a,b,cin;
output sum,co;
assign sum=a^b^cin;
assign co=(a&b)|(a&cin)|(cin&b);
endmodule

module rca_4bit (a,b,sum,cout);
output [3:0]sum;
output cout;
input [3:0]a;
input [3:0]b;
wire [2:0]c;
  full_adder u0 (a[0],b[0],1'b0,sum[0],c[0]);
  full_adder u1 (a[1],b[1],c[0],sum[1],c[1]);
  full_adder u2 (a[2],b[2],c[1],sum[2],c[2]);
  full_adder u3 (a[3],b[3],c[2],sum[3],cout);
endmodule

```
```
// Design of 4-bit ripple carry adder using Dataflow modelling

module rca_4bit (a,b,sum,cout);
output [3:0]sum;
output cout;
input [3:0]a;
input [3:0]b;
input cin;
assign {cout, sum}=a+b;
endmodule
```
```
// Testbench code of 4-bit ripple carry adder

`timescale 1ns / 1ps
module tb_rca_4bit;
// Inputs
reg [3:0] a;
reg [3:0] b;

// Outputs
wire [3:0] sum;
wire cout;

// Instantiate the Unit Under Test (UUT)
rca_4bit dut(a,b,sum,cout);
initial begin
  $monitor("a=%0d b=%0d sum=%0d cout=%b",a,b,sum,cout);
// Initialize Inputs
a = 0;
b = 0;
// Wait 100 ns for global reset to finish
#100;
a = 5;
b = 6;
// Wait 100 ns for global reset to finish
#100;
end
endmodule

```
#### b) 4-bit Adder–cum-Subtractor.

```
// Design of 4-bit Adder–cum-Subtractor using Structural modelling

module fa(a,b,cin,sum,co);

input a,b,cin;
output sum,co;

assign sum=a^b^cin;
assign co=(a&b)|(a&cin)|(cin&b);

endmodule

module addsub4bit ( a ,b ,add0sub1, sumdiff ,carryborrow );

output [3:0] sumdiff ;
output carryborrow ;

input [3:0] a ;
input [3:0] b ;

input add0sub1;

wire [3:0]x;

xor(x[0],b[0],add0sub1);
xor(x[1],b[1],add0sub1);
xor(x[2],b[2],add0sub1);
xor(x[3],b[3],add0sub1);

fa u0 (a[0],x[0],add0sub1,sumdiff[0],c1);
fa u1 (a[1],x[1],c1,sumdiff[1],c2);
fa u2 (a[2],x[2],c2,sumdiff[2],c3);
fa u3 (a[3],x[3],c3,sumdiff[3],carryborrow);

endmodule
```
```
// Design of 4-bit Adder–cum-Subtractor using Dataflow modelling

module addsub4bit ( a ,b ,add0sub1, sumdiff ,carryborrow );

output [3:0] sumdiff ;
output carryborrow ;

input [3:0] a ;
input [3:0] b ;

input add0sub1;

assign {carryborrow, sumdiff} = add0sub1 ? a-b : a+b;

endmodule

```
```
// Testbench of 4-bit Adder–cum-Subtractor

`timescale 1ns/1ns
module tb_addsub4bit;

wire [3:0] sumdiff ;
wire carryborrow ;

reg [3:0] a ;
reg [3:0] b ;

reg add0sub1;

integer i;

addsub4bit ins1( a ,b ,add0sub1, sumdiff ,carryborrow );

initial begin 
    a=4; b=4; add0sub1=0;
#10 a=10; b=4; add0sub1=0;
#10 a=6; b=4; add0sub1=0;
#10 a=3; b=7; add0sub1=0;
#10 a=4; b=4; add0sub1=1;
#10 a=10; b=4; add0sub1=1;
#10 a=6; b=4; add0sub1=1;
#10 a=3; b=7; add0sub1=1;
#10;
$finish;
end

initial
$monitor("simtime=%g, a=%0d, b=%0d, add0sub1=%b, sumdiff=%0d, carryborrow=%b", $time,a,b,add0sub1,sumdiff,carryborrow);

endmodule
```

#### c) 2-digit BCD adder / subtractor.

```
// Design of 2-digit BCD adder/subtractor using Structural modelling 

module fulladder1bit(a, b, c_in,sum, c_out);
  input a, b, c_in;
  output sum, c_out;
  wire s1, c1, c2;
  xor (s1,a, b);
  and ( c1,a, b);
  xor (sum,s1, c_in);
  and (c2,s1, c_in);
  or (c_out,c2, c1);
endmodule

module fulladder4bit(a, b,sum, c_out);
  //i/o port declaration
  input [3:0] a, b;
  output [3:0] sum;
  output c_out;
  //internal net
  wire c1, c2, c3;
  fulladder1bit fa0(a[0], b[0], 1'b0, sum[0], c1);
  fulladder1bit fa1(a[1], b[1], c1, sum[1], c2);
  fulladder1bit fa2(a[2], b[2], c2, sum[2], c3);
  fulladder1bit fa3(a[3], b[3], c3, sum[3], c_out);
endmodule

module onedigit_bcd_adder( A, B, CIN, F, COUT);
  input [3:0] A, B;
  input CIN;
  output [3:0] F;
  output COUT;
  wire [3:0] Z,S;
  wire k,w1,w2,w3;
  fulladder4bit add0(A, B, CIN, Z, k);
  and (w1,Z[3],Z[2]);
  and (w2,Z[3],Z[1]);
  or (COUT,k,w1,w2);
  assign S = {1'b0,COUT,COUT,1'b0};
  fulladder4bit add1(Z, S, 0,F,w3);
endmodule

// 9's complement_generator

module complement_generator(B, M, x);
  input [3:0]B;
  input M;
  output [3:0]x;
  wire w1,w2,w3,w4,w5,w6,w7,w8,w9;
  xor (x[0],B[0],M);
  assign x[1]=B[1];
  xor (w5,B[1],B[2]);
  and(w9,w5,M);
  not (w1,M);
  and (w6,B[2],w1);
  or (x[2],w9,w6);
  not (w2,B[1]);
  not (w3,B[2]);
  not (w4,B[3]);
  and (w8,M,w2,w3,w4);
  and (w7,B[3],w1);
  or (x[3],w8,w7);
endmodule

// 1-digit bcd Adder-cum-subtractor

module onedigit_bcd_add_sub(A, B, CIN, M, F, COUT);
  input [3:0]A,B;
  input CIN,M;
  output [3:0]F;
  output COUT;
  wire [3:0]W;
  complement_generator comgen0(B,M,W);
  onedigit_bcd_adder add0(A, W, CIN, F, COUT);
endmodule

// 2-digit bcd Adder-cum-subtractor

module twodigit_bcd_add_sub(A2,A1,B2,B1,CIN,M,F2,F1,COUT);
  input [3:0]A1,B1;
  input [7:4]A2,B2;
  input CIN,M;
  output [3:0]F1;
  output [7:4]F2;
  output COUT;
  wire COUT1;
  wire [3:0]W1;
  wire [7:4]W2;
  complement_generator comgen0(B1,M,W1);
  complement_generator comgen1(B2,M,W2);
  onedigit_bcd_adder add0(A1,W1,CIN, F1, COUT1);
  onedigit_bcd_adder add1(A2,W2,COUT1,F2, COUT);
endmodule
```
```
// Design of 2-digit BCD adder/subtractor using Dataflow modelling

module twodigit_bcd_add_sub(A2,A1,B2,B1,CIN,M,F2,F1,COUT);
  input [3:0]A1,B1;
  input [3:0]A2,B2;
  input CIN,M;
  output [3:0]F1;
  output [3:0]F2;
  output COUT;
  wire COUT1;
  wire [3:0]W1;
  wire [3:0]W2;
  
  wire [3:0]B1_9scomplement;
  wire [3:0]B2_9scomplement;
 
  assign B1_9scomplement= M ? ( 4'b1001 - B1) : B1;
  
  // First digit adder
  
	assign {C1,W1}= A1 + B1_9scomplement + CIN;

	assign and1a=W1[3]&W1[2];
	assign and1b=W1[3]&W1[1];
	assign COUT1=C1|and1a|and1b;

	assign F1 = W1 + {1'b0,COUT1,COUT1,1'b0};

    assign B2_9scomplement= M ? ( 4'b1001 - B2) : B2;
 
// Second digit adder

  	assign {C2, W2}= A2 + B2_9scomplement + COUT1;

	assign and2a=W2[3]&W2[2];
	assign and2b=W2[3]&W2[1];
	assign COUT=C2|and2a|and2b;

	assign F2=W2+{1'b0,COUT,COUT,1'b0};

endmodule
```
```
// Testbench of 2-digit BCD adder/subtractor

`timescale 1ns/1ns

module tb;

  reg [3:0]A1,B1;
  reg [3:0]A2,B2;
  reg CIN,M;
  wire [3:0]F1;
  wire [3:0]F2;
  wire COUT;
  
 twodigit_bcd_add_sub dut(A2,A1,B2,B1,CIN,M,F2,F1,COUT);
 
  initial begin
  
       CIN=0;
  
       A2=5; A1=5; B2=9; B1=9; M=0;
  
  #10  A2=5; A1=5; B2=3; B1=3;  M=0;
  
  #10  A2=5; A1=5; B2=1; B1=1;  M=0;
  
  #10  A2=6; A1=6; B2=9; B1=9;  M=0;
  
  #10  A2=5; A1=5; B2=3; B1=3; M=1; 
  
  #10  A2=5; A1=5; B2=1; B1=1; M=1; 
  
  #10  A2=9; A1=8; B2=2; B1=3; M=1; 
  #10  A2=9; A1=8; B2=2; B1=3; M=0;     
   
  #10;
  
  $finish;
  
  end
  
  initial
  $monitor($time, "ns A2=%0d A1=%0d  B2=%0d  B1=%0d  CIN=%b M=%b COUT=%b  F2=%0d F1=%0d", A2,A1,B2,B1,CIN,M,COUT,F2,F1);
  
  endmodule
```

#### d) 4-bit carry look-ahead adder
```
// Design of 4-bit carry look-ahead adder using Structural modelling

module cla_4bit(
    input [3:0] a,
    input [3:0] b,
    input cin,
    output [3:0] sum,
    output cout
);

wire [3:0] G, P, C;
wire Cout;

Generate G1(a[0], b[0], G[0]);
Generate G2(a[1], b[1], G[1]);
Generate G3(a[2], b[2], G[2]);
Generate G4(a[3], b[3], G[3]);

Propagate P1(a[0], b[0], P[0]);
Propagate P2(a[1], b[1], P[1]);
Propagate P3(a[2], b[2], P[2]);
Propagate P4(a[3], b[3], P[3]);

CarryLookAhead CLA(G, P, cin, C);
  
xorgate ins1(P[0],cin, sum[0]);
xorgate ins2(P[1],C[0], sum[1]);
xorgate ins3(P[2],C[1], sum[2]);
xorgate ins4(P[3],C[2], sum[3]);
    
assign cout=C[3];
endmodule

// Generate carry signals
module Generate (
    input [3:0] A,
    input [3:0] B,
    output [3:0] G
);
    
assign G = A & B;

endmodule

// Propagate carry signals
module Propagate (
    input [3:0] A,
    input [3:0] B,
    output [3:0] P
);
    
assign P = A ^ B;

endmodule

// Carry Look-Ahead Logic
module CarryLookAhead (
    input [3:0] G,
    input [3:0] P,
    input Cin,
    output [3:0] C
);
    
assign C[0] = G[0] | (P[0] & Cin);
assign C[1] = G[1] | (P[1] & (G[0] | (P[0] & Cin)));
assign C[2] = G[2] | (P[2] & (G[1] | (P[1] & (G[0] | (P[0] & Cin)))));
assign C[3] = G[3] | (P[3] & (G[2] | (P[2] & (G[1] | (P[1] & (G[0] | (P[0] & Cin)))))));

endmodule

module xorgate(a,b,y);
  
  input a,b;
  output y;
  
  assign y=a^b;
  
endmodule

```
```
// Design of 4-bit carry look-ahead adder using Dataflow modelling

module cla_4bit(a,b,cin,sum,cout);
input[3:0] a,b;
input cin;
output [3:0] sum;
output cout;
wire p0,p1,p2,p3,g0,g1,g2,g3,c1,c2,c3,c4;

assign p0=(a[0]^b[0]);
assign p1=(a[1]^b[1]);
assign p2=(a[2]^b[2]);
assign p3=(a[3]^b[3]);

assign g0=(a[0]&b[0]);
assign g1=(a[1]&b[1]);
assign g2=(a[2]&b[2]);
assign g3=(a[3]&b[3]);

assign c0=cin;

assign c1=g0|(p0&cin);
assign c2=g1|(p1&g0)|(p1&p0&cin);
assign c3=g2|(p2&g1)|(p2&p1&g0)|(p1&p1&p0&cin);
assign c4=g3|(p3&g2)|(p3&p2&g1)|(p3&p2&p1&g0)|(p3&p2&p1&p0&cin);

assign sum[0]=p0^c0;
assign sum[1]=p1^c1;
assign sum[2]=p2^c2;
assign sum[3]=p3^c3;

assign cout=c4;

endmodule
```
```
// Testbench of 4-bit carry look-ahead adder

`timescale 1ns / 1ps

module tb_cla_4bit;
// Inputs
reg [3:0] a;
reg [3:0] b;
reg cin;

// Outputs
wire [3:0] sum;
wire cout;

// Instantiate the Unit Under Test (UUT)
cla_4bit uut (
.a(a),
.b(b),
.cin(cin),
.sum(sum),
.cout(cout)
);
initial begin
// Initialize Inputs
a = 0;
b = 0;
cin = 0;
// Wait 100 ns for global reset to finish
#100;
a = 5;
b = 6;
cin = 1;
// Wait 100 ns for global reset to finish
#100;
end
endmodule
```

#### e) 4-bit comparator
```
// Design of 4-bit comparator using Structural modelling

// comparator 1-bit

module comp1bit(a,b,gt,lt,eq,agtb,altb,aeqb);
input a,b,gt,lt,eq;
output agtb,altb,aeqb;

assign agtb=(a&(~b))|(a&gt)|(gt&(~b));
assign altb=(~a&(b))|(b&lt)|(lt&(~a));
assign aeqb=(~a&(~b)&eq)|((a)&(b)&eq);
endmodule

//Comparator 4-bit

module comp4bit(a,b,agtb,altb,aeqb);
input [3:0]a,b;
output agtb,altb,aeqb;

wire agtb1,altb1,aeqb1,agtb2,altb2,aeqb2,agtb3,altb3,aeqb3;

 comp1bit ins1(a[0],b[0],1'b0,1'b0,1'b1,agtb1,altb1,aeqb1);
 comp1bit ins2(a[1],b[1],agtb1,altb1,aeqb1,agtb2,altb2,aeqb2);
 comp1bit ins3(a[2],b[2],agtb2,altb2,aeqb2,agtb3,altb3,aeqb3);
 comp1bit ins4(a[3],b[3],agtb3,altb3,aeqb3,agtb,altb,aeqb);

endmodule
```
```
// Design of 4-bit comparator using Dataflow modelling

module comp4bit(a,b,agtb,altb,aeqb);
input [3:0]a,b;
output agtb,altb,aeqb;

assign agtb = a > b;
assign altb = a < b;
assign aeqb = a == b;

endmodule
```
```
// Testbench of 4-bit comparator

`timescale 1ns / 1ps

module tb_comp4bit;

	// Inputs
	reg [3:0] a;
	reg [3:0] b;

	// Outputs
	wire agtb;
	wire altb;
	wire aeqb;

// Instantiate the Unit Under Test (UUT)
	comp4bit uut (
		.a(a), 
		.b(b), 
		.agtb(agtb), 
		.altb(altb), 
		.aeqb(aeqb)
	);

	initial 
			begin
			$monitor($time, "a=%0d b=%0d agtb=%b aeqb=%b altb=%b",a,b,agtb,aeqb,altb);
				// Initialize Inputs
					   a = 0; b = 0;
					#5 a = 5; b = 3;
					#5 a = 5; b = 7;
					#5 a = 5; b = 5;	
			end
      
endmodule
```

## 2. Write a Verilog HDL program in behavioral model for

#### a) 8:1 multiplexer
```
// Design of Mux 8x1 using behavioral modelling

module mux8x1(i,s,y);
  
  input [7:0]i;
  input [2:0]s;
  output reg y;
  
  always@*
    case(s)
      3'b000: y=i[0];
      3'b001: y=i[1];
      3'b010: y=i[2];
      3'b011: y=i[3];
      3'b100: y=i[4];
      3'b101: y=i[5];
      3'b110: y=i[6];
      3'b111: y=i[7];
    endcase
  
endmodule
```
```
// Testbench of Mux 8x1

`timescale 1ns/1ps

module tb_mux8x1;
  
  reg [7:0]i;
  reg [2:0]s;
  wire y;
  
  integer j;
  mux8x1 dut(i,s,y);
  
  initial 
    
    begin
      $monitor($time, "ns i=%b s=%0d y=%b", i,s,y);
      i=8'b1010_1010;
      for(j=0; j<8; j=j+1)
        begin
          s=j;
          #10;
        end
    end
endmodule
```

#### b) 3:8 decoder
```
// Design of Decoder 3:8 using behavioral modelling

module decoder_3to8 (input [2:0] in, output reg [7:0] out);

always @* begin
    case (in)
        3'b000: out = 8'b00000001;
        3'b001: out = 8'b00000010;
        3'b010: out = 8'b00000100;
        3'b011: out = 8'b00001000;
        3'b100: out = 8'b00010000;
        3'b101: out = 8'b00100000;
        3'b110: out = 8'b01000000;
        3'b111: out = 8'b10000000;
        default: out = 8'b00000000;
    endcase
end

endmodule
```
```
// Testbench of Decoder 3-to-8

`timescale 1ns/1ps
module tb_decoder_3to8;

    // Define the input and output signals
    reg [2:0] in;
    wire [7:0] out;
    integer i;
    
     // Instantiate the module under test
    decoder_3to8 dut (.in(in), .out(out));
     
    // Initialize the inputs
    initial begin
      for(i=0; i<8; i=i+1) 
        begin
        in = i;
        #10;
        $display("Input: %b, Output: %b", in, out);
       
    end
    end
endmodule
```
#### c) 8:3 encoder
```
// Design of encoder 8-to-3 using behavioral modelling

module encoder_8to3 (input [7:0] in, output reg [2:0] out);

always @* begin
    case (in)
        8'b00000001: out = 3'b000;
        8'b00000010: out = 3'b001;
        8'b00000100: out = 3'b010;
        8'b00001000: out = 3'b011;
        8'b00010000: out = 3'b100;
        8'b00100000: out = 3'b101;
        8'b01000000: out = 3'b110;
        8'b10000000: out = 3'b111;
        default: out = 3'b000;
    endcase
end

endmodule
```
```
// Testbench of encoder 8-to-3

`timescale 1ns/1ps
module tb_encoder_8to3;

   // Define the input and output signals
    reg [7:0] in;
    wire [2:0] out;
    
     // Instantiate the module under test
    encoder_8to3 dut (.in(in),.out(out));
    
    // Initialize the inputs
    initial 
      begin
        in = 8'b00000001;
        #10;
        $display("Input: %b, Output: %b", in, out);
      
	repeat(8)
		begin
	    in = in << 1;
        #10;
        $display("Input: %b, Output: %b", in, out);
       end		
      end
endmodule
```
#### d) 8-bit parity generator and checker
```
// Design of 8-bit parity generator and checker using behavioural modeling

module parity_generator_checker (
  input [7:0] data_in,
  input odd_parity_control,
  input even_parity_control,
  input odd_parity_check_control,
  input even_parity_check_control,
  output [8:0] data_out_with_parity_bit,
  output parity_out,
  output reg parity_error
);

  reg parity_bit;

  assign parity_out = parity_bit;
  
  always @* begin  // Parity generator
    if (odd_parity_control) begin
      parity_bit = ~^data_in;
    end else if (even_parity_control) begin
      	  parity_bit = ^data_in;
    end else if (odd_parity_check_control) 
	begin

           if (~^{data_in,parity_bit})
	       begin
	        $display($time,"ns For ODD parity: Let's say from DUT the parity_bit is %b\n",parity_bit);
	        parity_error = 1'b1;
	       end
           else
               begin
	         $display($time,"ns For ODD parity: Let's say from DUT the parity_bit is %b\n",parity_bit);
	         parity_error = 1'b0;
	       end
        end	  
	else if (even_parity_check_control) 
	begin
                 if (^{data_in,parity_bit})
	           begin
	             $display($time,"ns For EVEN parity: Let's say from DUT the parity_bit is %b\n",parity_bit);
	             parity_error = 1'b1;
                   end
		else
               	   begin
	            $display($time,"ns For EVEN parity: Let's say from DUT the parity_bit is %b\n",parity_bit);
	            parity_error = 1'b0;
	           end
        end
	else 
	begin
        parity_error = 1'b0;
	parity_bit=1'b0;
        end
    end

  // Output data with parity bit
  assign data_out_with_parity_bit = {data_in, parity_bit};

endmodule
```
```
// Testbench for 8-bit parity generator and checker

module parity_generator_checker_tb;

  reg [7:0] data_in;
  reg odd_parity_control;
  reg even_parity_control;
  reg odd_parity_check_control;
  reg even_parity_check_control;
  wire [8:0] data_out_with_parity_bit;
  wire parity_out;
  wire parity_error;

//reg parity_bit;
parity_generator_checker dut(data_in,odd_parity_control,even_parity_control,odd_parity_check_control,
                             even_parity_check_control,data_out_with_parity_bit,parity_out,parity_error);
  
initial 
begin

// Test for Parity generator 
  
// Initialize odd_parity_check_control,even_parity_check_control to ZERO
  
// Start test for even parity
  
$display($time, "ns -------------Start of Even parity Generation--------------\n");
#2;
$display($time, "ns ==Let's give EVEN number of 1's in Data and see at the end bit which is parity==\n");
data_in=8'b1010_1010; odd_parity_control=1'b0; even_parity_control=1'b1;
 
// Wait for some time
#2;
  
$display($time, "ns Data given is data_in=%b  ", data_in,
  "even_parity_control is high %b and odd_parity_control is low %b ", even_parity_control,odd_parity_control,
  "and the generated data_out_with_parity_bit is %b\n", data_out_with_parity_bit);
			  
// Wait for some time
#2;
$display($time, "ns ==Let's give ODD number of 1's in Data and see at the end bit which is parity==\n");
data_in=8'b0111_1010; odd_parity_control=1'b0; even_parity_control=1'b1;
 
// Wait for some time
#2;
$display($time,"ns Data given is data_in=%b  ", data_in,
"even_parity_control is high %b and odd_parity_control is low %b ", even_parity_control,odd_parity_control,
  "and the generated data_out_with_parity_bit is %b\n", data_out_with_parity_bit);
			  
// Wait for some time
 #2;
$display($time,"ns -------------END of Even parity Generation--------------\n");
#2;

// Start test for odd parity
$display($time,"ns ------------Start of ODD parity Generation------------\n");
#2;
$display($time, "ns ==Let's give EVEN number of 1's in Data and see at the end bit which is parity==\n");
data_in=8'b1010_1010; odd_parity_control=1'b1; even_parity_control=1'b0;
 
// Wait for some time
#2;
$display($time,"ns Data given is data_in=%b  ", data_in,
"even_parity_control is low %b and odd_parity_control is high %b ", even_parity_control,odd_parity_control,
  "and the generated data_out_with_parity_bit is %b\n", data_out_with_parity_bit);
			  
// Wait for some time
#2;
$display($time, "ns ==Let's give ODD number of 1's in Data and see at the end bit which is Parity==\n");
   data_in=8'b111_1010; odd_parity_control=1'b1; even_parity_control=1'b0;
 
// Wait for some time
#2;
$display($time,"ns Data given is data_in=%b  ", data_in,
  "even_parity_control is low %b and odd_parity_control is high %b ", even_parity_control,odd_parity_control,
  "and the generated data_out_with_parity_bit is %b\n", data_out_with_parity_bit);
			  
// Wait for some time
#2;
$display($time,"ns ------------END of ODD parity Generation------------\n");
#2;
  
$display($time,"ns ------------Start of Even Parity check------------\n");
  
data_in=8'b0111_1010; dut.parity_bit=1'b0;odd_parity_control=1'b0; even_parity_control=1'b0;
even_parity_check_control=1'b1; odd_parity_check_control=1'b0;
#2;
$display($time, "ns data_in=%b ",data_in,
			"parity_bit=%b ", dut.parity_bit,
			"odd_parity_control=%b ", odd_parity_control,
			"even_parity_control=%b ", even_parity_control,
			"even_parity_check_control=%b ",even_parity_check_control,
			"odd_parity_check_control=%b ",odd_parity_check_control,
			"parity_error=%b\n ",parity_error
			
			);
#2;
data_in=8'b0111_1010; dut.parity_bit=1'b1;odd_parity_control=1'b0; even_parity_control=1'b0;
even_parity_check_control=1'b1; odd_parity_check_control=1'b0;
#2;
$display($time, "ns data_in=%b ",data_in,
			"parity_bit=%b ", dut.parity_bit,
			"odd_parity_control=%b ", odd_parity_control,
			"even_parity_control=%b ", even_parity_control,
			"even_parity_check_control=%b ",even_parity_check_control,
			"odd_parity_check_control=%b ",odd_parity_check_control,
			"parity_error=%b\n ",parity_error
			
			);
#2;
$display($time, "ns -----------END of Even Parity check------------\n");
#2;
$display($time,"ns ----------Start of ODD Parity check-------------\n");
#2;
data_in=8'b0111_1010; dut.parity_bit=1'b1; odd_parity_control=1'b0; even_parity_control=1'b0;
even_parity_check_control=1'b0; odd_parity_check_control=1'b1;
#2;
$display($time, "ns data_in=%b ",data_in,
			"parity_bit=%b ", dut.parity_bit,
			"odd_parity_control=%b ", odd_parity_control,
			"even_parity_control=%b ", even_parity_control,
			"even_parity_check_control=%b ",even_parity_check_control,
			"odd_parity_check_control=%b ",odd_parity_check_control,
			"parity_error=%b\n",parity_error
			
			);
#2;
data_in=8'b0111_1010; dut.parity_bit=1'b0; odd_parity_control=1'b0; even_parity_control=1'b0;
even_parity_check_control=1'b0; odd_parity_check_control=1'b1;
  
#2;
$display($time, "ns data_in=%b ",data_in,
			"parity_bit=%b ", dut.parity_bit,
			"odd_parity_control=%b ", odd_parity_control,
			"even_parity_control=%b ", even_parity_control,
			"even_parity_check_control=%b ",even_parity_check_control,
			"odd_parity_check_control=%b ",odd_parity_check_control,
			"parity_error=%b\n",parity_error
			
			);
#2;
    
$display($time,"ns -----------END of ODD Parity check------------\n");
#2;
$display($time, "ns ==When no control signal is active,The default parity bit and error both are set to ZERO==\n");
#2;
data_in=8'b0111_1010; odd_parity_control=1'b0; even_parity_control=1'b0;
even_parity_check_control=1'b0; odd_parity_check_control=1'b0;
#2;
$display($time,"ns data_in=%b  ",data_in,
		   "odd_parity_control=%b  ", odd_parity_control,
			"even_parity_control=%b  ", even_parity_control,
			"even_parity_check_control=%b  ",even_parity_check_control,
			"odd_parity_check_control=%b  ",odd_parity_check_control,
			"data_out_with_parity_bit_0=%b   ", data_out_with_parity_bit,
			"parity_out=%b  ",parity_out,
			"parity_error=%b\n",parity_error
			
			);
  
$display($time, "ns ------------All tests ends----------");
end
endmodule
```
### 3. Write a Verilog HDL program in a Hierarchical structural model for
#### a) 16:1 multiplexer realization using 4:1 multiplexer
```
// MUX 4:1
module mux4x1(input [3:0]i,
              input [1:0]select,
              output reg y);
  
  
  always@*
    case(select)
      
      2'b00: y=i[0];
      2'b01: y=i[1];
      2'b10: y=i[2];
      2'b11: y=i[3];
      
    endcase
endmodule

// MUX 16:1 using MUX 4:1
module mux16x1( input [15:0]I, input [3:0]select, output y);
    
  // level 1
  mux4x1 dut1(i[3:0],select[1:0],y1);
  mux4x1 dut2(i[7:4],select[1:0],y2);
  mux4x1 dut3(i[11:8],select[1:0],y3);
  mux4x1 dut4(i[15:12],select[1:0],y4);
  
  //level 2
  
  mux4x1 dut5({y4,y3,y2,y1},select[3:2],y);
  
endmodule
```
```
// Testbench for MUX 16:1 
module test_mux16x1;

   // Inputs
   reg [15:0] data_in;
   reg [3:0] select;

   // Outputs
   wire out;

   // Instantiate the DUT
   mux16x1 dut (.i(data_in),.select(select),.y(out));

   // Initialize inputs
   initial begin
      data_in = 16'b1010_1010_1010_1010;
      select = 4'b0000;
   end

   // Stimulus
   always #5 select = select + 1'b1;
  
  initial
   #50 $finish;

   // Display output
 initial
   $monitor($time, "data_in=%b select=%0d out = %b", data_in,select,out);

endmodule
```
#### b) 3:8 decoder realization through 2:4 decoder
```
module decoder2x4_enable(input [1:0] sel, input enable, output reg [3:0] out);

    always @(sel, enable) begin
        if (enable) begin
            case(sel)
                2'b00: out = 4'b0001;
                2'b01: out = 4'b0010;
                2'b10: out = 4'b0100;
                2'b11: out = 4'b1000;
            endcase
        end else begin
            out = 4'b0000;
        end
    end

endmodule

module decoder3x8(input [2:0]sel, output [7:0]out);
  
  
  decoder2x4_enable uut1(sel[1:0],~sel[2], out[3:0]);
  decoder2x4_enable uut2(sel[1:0],sel[2], out[7:4]);
  
endmodule
```
```
// Testbench  for Decoder 3-to-8
module testbench();
  
   // Declare the inputs and outputs for the testbench
    reg [2:0] sel;
    wire [7:0] out;

    // Instantiate the DUT (device under test)
    decoder3x8 dut(.sel(sel), .out(out));

    // Initialize the inputs
    initial begin
        sel = 0;
        #10;
        sel = 1;
        #10;
        sel = 2;
        #10;
        sel = 3;
        #10;
        sel = 4;
        #10;
        sel = 5;
        #10;
        sel = 6;
        #10;
        sel = 7;
        #10;
        $finish;
    end

    // Monitor the outputs and print them to the console
    always @(out) begin
        $display("sel = %d, out = %b", sel, out);
    end

endmodule
```
#### c) 8-bit comparator using 4-bit comparators and additional logic
```
// Comparator 1-bit

module comp1bit(a,b,gt,lt,eq,agtb,altb,aeqb);

input a,b,gt,lt,eq;
output agtb,altb,aeqb;

assign agtb=(a&(~b))|(a&gt)|(gt&(~b));
assign altb=(~a&(b))|(b&lt)|(lt&(~a));
assign aeqb=(~a&(~b)&eq)|((a)&(b)&eq);

endmodule

// Comparator 4-bit

module comp4bit(a,b,gt,lt,eq,agtb,altb,aeqb);
input [3:0]a,b;
input gt,lt,eq;
output agtb,altb,aeqb;

 comp1bit ins1(a[0],b[0],gt,lt,eq,agtb1,altb1,aeqb1);
 comp1bit ins2(a[1],b[1],agtb1,altb1,aeqb1,agtb2,altb2,aeqb2);
 comp1bit ins3(a[2],b[2],agtb2,altb2,aeqb2,agtb3,altb3,aeqb3);
 comp1bit ins4(a[3],b[3],agtb3,altb3,aeqb3,agtb,altb,aeqb);

endmodule

// Comparator 8-bit

module comp8bit(a,b,gt,lt,eq,agtb,altb,aeqb);
  input [7:0]a,b;
  input gt,lt,eq;
  output agtb,altb,aeqb;

   

  comp4bit  dut1(a[3:0],b[3:0],gt,lt,eq,agtb1,altb1,aeqb1);
  comp4bit dut2(a[7:4],b[7:4],agtb1,altb1,aeqb1,agtb,altb,aeqb);
  
endmodule
```
```
// Testbench for Comparator 8-bit
`timescale 1ns/1ns
module tb;
  
  reg [7:0]a,b;
  reg gt,lt,eq;
  wire agtb,altb,aeqb;
  
  integer i;
  
 comp8bit dut(a,b,gt,lt,eq,agtb,altb,aeqb);
  
  initial begin
    gt=0;lt=0;eq=1;
    $monitor("a %0d and b %0d agtb=%b altb=%b aeqb=%b", a,b,agtb,altb, aeqb);
    for(i=0; i<8; i=i+1) begin
      a=$random; b=$random;
      #10;
      
    end
    a=5; b=5;
    #10;
    $finish;
  end
endmodule
```

## 4. Write a Verilog HDL program in behavioral model for D, T, and JK flip flops, shift registers, and counters.
### Flip flops
#### 4a) D Flip flop
```
module dff(clk,rst,d,q);
  
  input clk,rst,d;
  output reg q;
  
  always@(posedge clk)
    begin
      if(rst)   // Active high reset
        q<=0;
      else
        q<=d;      
    end
  
endmodule
```
```
`timescale 1ns/1ps
module tb;
  
  reg clk,rst,d;
  wire  q;
  
  dff dut(clk,rst,d,q);
  
  initial
    begin
      clk<=0;
    end
  
  always #5 clk=~clk;
  
  initial
    begin
       
       rst<=1; d<=1;
       #1;
       @(posedge clk);
       #1;
       
       rst<=0; d<=0;
       #1;
       @(posedge clk);
       #1; 
       
       d<=1;
       #1;
       @(posedge clk);
       #1;
      
      d<=0;
      #1;
      @(posedge clk);
      #1;
      
      #50;
      $stop;
    end
  
  initial 
    $monitor("clk=%b rst=%b d=%b q=%b",clk,rst,d,q);
   
endmodule
```
#### 4b) T Flip flop
```
module tff(clk,rst,t,q);
  
  input clk,rst,t;
  output reg q;
  
  always@(posedge clk)
    begin
      if(rst)
        q<=0;
      else if(t)
        q<=~q;
      else
        q<=q;
    end  
  
endmodule
```
```
module tb_tff;
  
   reg clk,rst,t;
   wire q;
  
   tff dut(clk,rst,t,q);
  
  initial
    begin
      clk<=0;
    end
  
  always #5 clk=~clk;
  
  initial
    begin
       
       rst<=1; t<=1;
       #1;
       @(posedge clk);
       #1;
       
       rst<=0; t<=0;
       #1;
       @(posedge clk);
       #1; 
       
       t<=1;
       #1;
       @(posedge clk);
       #1;
      
       t<=1;
       #1;
       @(posedge clk);
       #1;
      
       t<=1;
       #1;
       @(posedge clk);
       #1;
      
      t<=0;
      #1;
      @(posedge clk);
      #1;
      
       t<=1;
       #1;
       @(posedge clk);
       #1;
      
       t<=1;
       #1;
       @(posedge clk);
       #1;
      
      #50;
      $stop;
    end
  
  initial 
    $monitor("clk=%b rst=%b t=%b q=%b",clk,rst,t,q);
   
endmodule
```
#### 4c) JK Flip flop
```
module jkff(j,k,clk,rst,q);
  
  input j,k,clk,rst;
  output reg q;
  
  always@(posedge clk)
    begin
      if(rst)
        q<=0;
      else
        case({j,k})
          2'b00: q<=q;
          2'b01: q<=0;
          2'b10: q<=1;
          2'b11: q<=~q;
        endcase
    end
endmodule
```
```
module tb_tff;
  
  reg j,k,clk,rst;
  wire q;
  
   jkff dut(j,k,clk,rst,q);
  
  initial
    begin
      clk<=0;
    end
  
  always #5 clk=~clk;
  
  initial
    begin
       
       rst<=1; j<=1; k<=0;
       #1;
       @(posedge clk);
       #1;
       
       rst<=0; j<=1; k<=0;
       #1;
       @(posedge clk);
       #1; 
       
      j<=1; k<=1;
       #1;
       @(posedge clk);
       #1;
      
      j<=0; k<=1;
       #1;
       @(posedge clk);
       #1;
      
       j<=1; k<=0;
       #1;
       @(posedge clk);
       #1;
      
      j<=0; k<=1;
      #1;
      @(posedge clk);
      #1;
      
      j<=1; k<=1;
       #1;
       @(posedge clk);
       #1;
      
      j<=1; k<=0;
       #1;
       @(posedge clk);
       #1;
      
      #50;
      $stop;
    end
  
  initial 
    $monitor("clk=%b rst=%b j=%b; k=%b;=%b q=%b",clk,rst,j,k,q);
   
endmodule
```
### Shift registers
#### 4a) Left shift register
```
//design file

module leftshift_4bit(clk,rst,data_in,data_out);
  input clk,rst;
  input data_in;
  output data_out;
  
  reg [3:0]temp;
  
  always@(posedge clk)
    begin
      if(rst)
        temp<=4'b0000;
      else
        temp<={temp[2:0],data_in};
    end
  
  assign data_out=temp[3];
  
endmodule
```
```
//testbench file

`timescale 1ns/1ps

module tb;

reg clk,rst;
reg data_in;
wire data_out;

leftshift_4bit uut(clk,rst,data_in,data_out);

initial
begin
clk<=0;
end

always #5 clk<=~clk;

initial
begin

rst<=0;
#1;
@(posedge clk);
#1;

rst<=1;
#1;
@(posedge clk);
#1;

rst<=0;data_in<=1;
#1;
@(posedge clk);
#1;

data_in<=0;
#1;
@(posedge clk);
#1;

data_in<=1;
#1;
@(posedge clk);
#1;

data_in<=0;
#1;
@(posedge clk);
#1;

data_in<=1;
#1;
@(posedge clk);
#1;

data_in<=0;
#1;
@(posedge clk);
#1;

data_in<=1;
#1;
@(posedge clk);
#1;

data_in<=0;
#1;
@(posedge clk);
#1;

end

initial
$monitor("clk %b,rst %b,data_in %b,data_out %b",clk,rst,data_in,data_out);

endmodule
```
#### 4b) Right shift register
```
//design file

module rightshift_4bit(clk,rst,data_in,data_out);
  input clk,rst;
  input data_in;
  output data_out;
  
  reg [3:0]temp;
  
  always@(posedge clk)
    begin
      if(rst)
        temp<=4'b0000;
      else
        temp<={data_in,temp[3:1]};
    end
  
  assign data_out=temp[0];
  
endmodule
```
```
//testbench file

`timescale 1ns/1ps

module tb;

reg clk,rst;
reg data_in;
wire data_out;

rightshift_4bit uut(clk,rst,data_in,data_out);

initial
begin
clk<=0;
end

always #5 clk<=~clk;

initial
begin

rst<=0;
#1;
@(posedge clk);
#1;

rst<=1;
#1;
@(posedge clk);
#1;

rst<=0;data_in<=1;
#1;
@(posedge clk);
#1;

data_in<=0;
#1;
@(posedge clk);
#1;

data_in<=1;
#1;
@(posedge clk);
#1;

data_in<=0;
#1;
@(posedge clk);
#1;

data_in<=1;
#1;
@(posedge clk);
#1;

data_in<=0;
#1;
@(posedge clk);
#1;

data_in<=1;
#1;
@(posedge clk);
#1;

data_in<=0;
#1;
@(posedge clk);
#1;

end

initial
$monitor("clk %b,rst %b,data_in %b,data_out %b",clk,rst,data_in,data_out);

endmodule
```
#### 4c) Universal shift register
```
// Universal Shift Register

module usr(input clk,rst,
           input [1:0]reg_op,
		   input sr_input, sl_input,
           input [3:0]pl_input,
		   output sr_output, sl_output,
		   output reg [3:0]q);
		   
		   
		always@(posedge clk)
		   
		   begin
			if(rst)
			begin
		    q<=4'b0000;
			end
			else
			case(reg_op)
			2'b00: q<=q;                    // No changes
			2'b01: q<={sr_input,q[3:1]};    // Shift right
			2'b10: q<={q[2:0], sl_input};   // Shift left
			2'b11: q<=pl_input;             // Parallel load
			endcase
		  end
		   
		assign sr_output=q[0];     // right shift
		assign sl_output=q[3];     // left shift
		   
endmodule
```
```
// Testbench for Universal Shift Register

`timescale 1ns/1ps
module tb;

reg clk,rst;
reg [1:0]reg_op;
reg sr_input, sl_input;
reg [3:0]pl_input;
wire sr_output, sl_output;
wire [3:0]q;

usr dut(clk,rst,reg_op,sr_input,sl_input,pl_input,sr_output,sl_output,q);

initial 
begin
clk=0;
rst=0;
end

always #5 clk=~clk;

initial
begin
$monitor($time,"clk=%b rst=%b reg_op=%b sr_input=%b sl_input=%b pl_input=%b sr_output=%b sl_output=%b q=%b",
clk,rst,reg_op,sr_input,sl_input,pl_input,sr_output,sl_output,q);
  
#1;
rst<=1;
@(posedge clk);
#1;
  
rst<=0; reg_op<=2'b00;      // No change
@(posedge clk);
#1;
  
reg_op<=2'b01; sr_input<=1;  // Right shift
@(posedge clk);
#1;
  
reg_op<=2'b01; sr_input<=0;   // Right shift
@(posedge clk);
#1;
  
reg_op<=2'b01; sr_input<=0;   // Right shift
@(posedge clk);
#1;
  
reg_op<=2'b01; sr_input<=0;   // Right shift
@(posedge clk);
#1;
  
reg_op<=2'b10; sl_input<=0; // Left shift
@(posedge clk);
#1;
  
reg_op<=2'b10; sl_input<=1; // Left shift
@(posedge clk);
  
#1;
reg_op<=2'b10; sl_input<=1; // Left shift
@(posedge clk);
#1;
  
reg_op<=2'b10; sl_input<=1; // Left shift
@(posedge clk);
#1;
  
reg_op<=2'b11; pl_input<=4'b1110;  //Parallel load
@(posedge clk);
#1;
  
$finish;

end
endmodule
```
### Counters
#### 4a) Up Counter
```
module counter_up4bit(clk,rst,count);
  
  input clk,rst;
  output reg [3:0]count;
  
  always@(posedge clk)
    begin
      if(rst)
        count<=0;
      else
        count<=count+1'b1;
    end
  
  endmodule
```
```
`timescale 1ns/1ps
module tb_counter_up4bit;
  
  reg clk,rst;
  wire [3:0]count;
  
  
  counter_up4bit dut(clk,rst,count);
  
  initial begin
    clk<=0; rst<=0;
  end
  
  always #5 clk<=~clk;
  
  initial begin
    @(posedge clk); rst<=1;
    @(posedge clk); rst<=0;
  end
  
  initial
  #200 $stop;
  
  initial
    $monitor("clk=%b rst=%b count=%0d", clk,rst,count);
endmodule
```
#### 4b) Down Counter
```
module counter_down4bit(clk,rst,count);
  
  input clk,rst;
  output reg [3:0]count;
  
  always@(posedge clk)
    begin
      if(rst)
        count<=0;
      else
        count<=count-1'b1;
    end
  
  endmodule
```
```
`timescale 1ns/1ps
module tb_counter_up4bit;
  
  reg clk,rst;
  wire [3:0]count;
  
  
  counter_down4bit dut(clk,rst,count);
  
  initial begin
    clk<=0; rst<=0;
  end
  
  always #5 clk<=~clk;
  
  initial begin
    @(posedge clk); rst<=1;
    @(posedge clk); rst<=0;
  end
  
  initial
  #200 $stop;
  
  initial
    $monitor("clk=%b rst=%b count=%0d", clk,rst,count);
endmodule
```
#### 4c) Up-down counter
```
module counter_updown4bit(clk,rst,ud,count);
  
  input clk,rst;
  input ud; // ud=1 for up count; ud=0 for down count
  output reg [3:0]count;
  
  always@(posedge clk)
    begin
      if(rst)
        count<=0;
      else if(ud)
        count<=count+1'b1;
      else
        count<=count-1'b1;
      
    end
  
  endmodule
```
```
`timescale 1ns/1ps
module tb;
  
  reg clk,rst;
  reg ud; // ud=1 for up count; ud=0 for down count
  wire [3:0]count;
 
  counter_updown4bit dut(clk,rst,ud,count);
   
  initial
    begin
      clk<=0;
      rst<=0;
    end
  
  always #5 clk<=~clk;
  
  initial
    begin
      $monitor("clk %b,rst %b,ud %b,count %0d",clk,rst,ud,count);
      @(posedge clk); rst<=1;
      @(posedge clk); rst<=0;
      repeat(20)
        begin
          @(posedge clk); ud<=1;
        end
       repeat(20)
        begin
          @(posedge clk); ud<=0;
        end
      $stop;
    end
endmodule
```

## 5. Write a Verilog HDL program in structural and behavioral models for
### a) 8-bit asynchronous up-down counter
#### Structural modeling
#### Behavioural modeling
```
module counter(clk,reset,up_down,count);
  //define input and output ports
  input clk,reset,up_down;
  output reg [7:0] count;
  
  always@(posedge clk, posedge reset) 
  begin
    if(reset)    //Set Counter to Zero
      count <= 0;
    else if(up_down)        //count up
      count <= count + 1;
    else            //count down
      count <= count - 1;
  end
endmodule
```
```
// Testbench
module counter_tb;
  reg clk,reset,ud;
  wire [7:0] count;
  // instance counter design
   counter dut(clk,reset,ud,count);
  //clock generator
  initial
	begin
	clk = 1'b0;
	repeat(30)
	#3 clk= ~clk;
	end
  //insert all the input signal
  initial
	begin
	reset=1'b1;
	#7 reset=1'b0;
	#35 reset=1'b1;
	end
  
  initial
	begin
	#5 ud=1'b1;
	#24 ud=1'b0;
	end
  
  //monitor all the input and output ports at times when any inputs changes its state
  initial begin 
    $monitor("time=%0d,reset=%b,ud=%b,count=%d", $time,reset,ud,count);
  end
endmodule
```
### b) 8-bit synchronous up-down counter
#### Structural modeling
#### Behavioural modeling
```
module counter(clk,reset,up_down,count);
  //define input and output ports
  input clk,reset,up_down;
  output reg [7:0] count;

  always@(posedge clk) 
  begin
    if(reset)    //Set Counter to Zero
      count <= 0;
    else if(up_down)        //count up
      count <= count + 1;
    else            //count down
      count <= count - 1;
  end
endmodule
```
```
// Code your testbench here

module counter_tb;
  reg clk,reset,ud;
  wire [7:0] count;
  // instance counter design
   counter dut(clk,reset,ud,count);
  //clock generator
  initial
	begin
	clk = 1'b0;
	repeat(30)
	#3 clk= ~clk;end
  //insert all the input signal
  initial
	begin
	reset=1'b1;
	#7 reset=1'b0;
	#35 reset=1'b1;
	end
  
  initial
	begin
	#5 ud=1'b1;
	#24 ud=1'b0;
	end
  
  //monitor all the input and output ports at times when any inputs changes its state
  initial begin 
    $monitor("time=%0d,reset=%b,ud=%b,count=%d", $time,reset,ud,count);
  end
endmodule
```
### 6. Write a Verilog HDL program for a 4-bit sequence detector through Moore state machines
### 7. Write a Verilog HDL program for 4-bit sequence detector through Mealy state machines
