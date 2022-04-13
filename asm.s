	GBLA ORDER
ORDER  SETA    2

	MACRO 
	mult_exp $opA , $opB , $res , $logtable, $exptable , $tmp
;	LDRB $tmp , [ $logtable , $opA ]
;	LDRB $res , [ $logtable , $opB ]
;	ADD  $tmp, $res
;	LDRB $res , [ $exptable , $tmp ]
;	RSB  $tmp , $opA , #0
;	AND  $tmp , $opB , $tmp , ASR #32
;	RSB  $tmp , #0
;	AND  $res , $res ,$tmp , ASR #32

;	mult_ht $opA , $opB , $res , $logtable, $exptable , $tmp

	mult_ft $opA , $opB , $res , $logtable , $tmp
	MEND
	
	MACRO
	mult_ht $opA , $opB , $res , $table1, $table2 , $tmp
	EOR $tmp , $opB , $opA , LSL #8
	LSR $tmp,#4
	LDRB $res , [$table1, $tmp]
	EOR $tmp , $opA , $opB , LSL #(28)
	ROR $tmp,#28
	LDRB $tmp , [$table2, $tmp]
	EOR $res , $tmp
	MEND

;	MACRO
;	mult_ht $opA , $opB , $res , $table1, $table2 , $tmp
;	EOR $tmp , $opB , $opA , LSL #8
;	;LSR $tmp,#4
;	LDRB r1 , [r0, r1, LSR#1]
;	EOR $tmp , $opA , $opB , LSL #28
;;	ROR $tmp,#28
;	LDRB r1 , [r0, r1, ROR#2]
;	EOR $res , $tmp
;	MEND
	

;	MACRO
;	mult_ht $opA , $opB , $res , $table1, $table2 , $tmp
;;	EOR $tmp , $opB , $opA , LSL # 8
;	LSR $tmp,#4
;	LDRB $res , [$table1];, $tmp]
;;	EOR $tmp , $opA , $opB , LSL #(28)
;	ROR $tmp,#28
;	LDRB $tmp , [$table2];, $tmp]
;	EOR $res , $tmp
;	MEND
	
	MACRO
	mult_ft $opA , $opB , $res , $pttab , $tmp0
	EOR $tmp0 , $opB , $opA , LSL # 8
	LDRB $res , [ $pttab, $tmp0 ]
	MEND


	MACRO 
	innerproduct_m $c_addr, $a_addr, $tmp1, $tmp2, $tmp3, $tmp4, $res, $logtable, $exptable, $order, $interval
	MOV  $res, #0
	SUB $a_addr,#1
	GBLA count
	GBLA counti
count   SETA    1
counti  SETA    0
	WHILE   count <= ORDER
;;	.irp ci, $orders //,0,1,2,3
		LDRB $tmp1 , [$c_addr, #(counti)]
		LDRB $tmp2 , [$a_addr,#1]!
		mult_exp $tmp1 , $tmp2 , $tmp3 , $logtable, $exptable , $tmp4
		EOR $res, $tmp3
count   SETA    count+1
counti  SETA    counti+$interval
	WEND
	SUB $a_addr,#(ORDER-1)
	MEND


	MACRO 
	tensorproduct_m $a_addr, $b_addr, $tmp1, $tmp2, $tmp3, $tmp4, $res_addr, $logtable, $exptable, $order
	;; aTb
	GBLA counti
	GBLA countj
counti  SETA   0
	WHILE   counti < $order
countj  SETA   0
		WHILE   countj < $order
			LDRB $tmp1 , [$a_addr, #counti]
			LDRB $tmp2 , [$b_addr, #countj]
			mult_exp $tmp1 , $tmp2 , $tmp3 , $logtable, $exptable , $tmp4
			STRB $tmp3, [$res_addr, #(counti*$order+countj)]
countj   SETA    countj+1
		WEND
counti  SETA    counti+1
	WEND
	MEND
	
	MACRO 
	onlinemul
	;;mul; input:r1,r2, sbox_t1t2r3t4:r10, sbox_alpha:r11, sbox_zs:r12, ouput: r0, logtable:r6, exptable:r7
	;;unsaved regs: r3,r4,r5,r8
	;;unused reg:r9
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;	
	LDRB r4,[r10],#1
	EOR  r1, r4
	LDRB r4,[r10],#1
	EOR  r2, r4
;;	mult_exp opA , opB , res , logtable, exptable , tmp
	mult_exp r1,    r2,  r0,    r6,       r7,        r8
	GBLA countii
countii  SETA   0
	WHILE   countii < ORDER
		LDRB r5, [r10],#1
		LDRB r3, [r12,#countii+ORDER]
	;;	mult_exp opA , opB , res , logtable, exptable , tmp
		mult_exp r1,    r3,  r4,    r6,       r7,        r8
		EOR  r5, r4
		LDRB r3, [r12,#countii]
	;;	mult_exp opA , opB , res , logtable, exptable , tmp
		mult_exp r2,    r3,  r4,    r6,       r7,        r8
		EOR  r5, r4
		LDRB r3, [r11,#countii]
	;;	mult_exp opA , opB , res , logtable, exptable , tmp
		mult_exp r5,    r3,  r4,    r6,       r7,        r8
		EOR  r0, r4
countii  SETA    countii+1
	WEND
	LDRB r5, [r10],#1
	EOR  r0,r5
	MEND
	
	MACRO
	onlinesbox
	;;sbox; input:r1, sbox_t1t2r3t4:r10, sbox_alpha:r11, sbox_zs:r12, ouput: r0, logtable:r6, exptable:r7
	;;unsaved regs: r3,r4,r5,r8
	;;unused reg:r9
	LDR  r3, =power2table
	LTORG
	LDRB r2,[r3,r1]
	PUSH {r2}
	;;mul; input:r1,r2, sbox_t1t2r3t4:r10, sbox_alpha:r11, sbox_zs:r12, ouput: r0, logtable:r6, exptable:r7
	;;unsaved regs: r3,r4,r5,r8; unused reg:r9
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;	
	onlinemul
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;	
	MOV  r1, r0
	LDR  r3, =power4table
	LTORG
	LDRB r2,[r3,r1]
	PUSH {r2}
	ADD r11, #ORDER
	;;mul; input:r1,r2, sbox_t1t2r3t4:r10, sbox_alpha:r11, sbox_zs:r12, ouput: r0, logtable:r6, exptable:r7
	;;unsaved regs: r3,r4,r5,r8; unused reg:r9
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;	
	onlinemul
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	MOV  r1, r0
	LDR  r3, =power16table
	LTORG
	LDRB r1,[r3,r1]
	POP  {r2}
	ADD r11, #ORDER
	;;mul; input:r1,r2, sbox_t1t2r3t4:r10, sbox_alpha:r11, sbox_zs:r12, ouput: r0, logtable:r6, exptable:r7
	;;unsaved regs: r3,r4,r5,r8; unused reg:r9
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;	
	onlinemul
	LTORG
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	MOV  r1, r0
	POP  {r2}
	ADD r11, #ORDER
	;;mul; input:r1,r2, sbox_t1t2r3t4:r10, sbox_alpha:r11, sbox_zs:r12, ouput: r0, logtable:r6, exptable:r7
	;;unsaved regs: r3,r4,r5,r8; unused reg:r9
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;	
	onlinemul
	LTORG
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	ADD r11, #ORDER
	MEND
	
	

	AREA ASM, code, readonly	
	THUMB
	EXPORT gfmul
	EXPORT innerproduct
	EXPORT innerproduct_inv
	EXPORT tensorproduct
	EXPORT matproduct
	EXPORT matproduct_inv
	EXPORT matproduct_oxrsum
	EXPORT summat
	EXPORT summat_inv
	EXPORT dotproduct
	EXPORT dotproduct2
	EXPORT spproduct
	EXPORT sboxac
	EXPORT matproduct_invalpha
	EXPORT matproductalpha
	IMPORT logtable
	IMPORT exptable
	IMPORT sbox_in
	IMPORT sbox_out
	IMPORT sbox_t1t2r3t4
	IMPORT sbox_alpha
	IMPORT sbox_zs
	IMPORT power2table
	IMPORT power4table
	IMPORT power16table
	IMPORT afftable
gfmul
	PUSH {r1,r2,r3,r4,r5,r6,r7,r8}
	MOV r2,r0
	LDR r7 , =logtable
	LDR r8 , =exptable
	;;mult_exp opA , opB , res , logtable, exptable , tmp
	mult_exp   r1  , r2  , r4  , r7,       r8       , r3
	MOV  r0,r4
	POP {r1,r2,r3,r4,r5,r6,r7,r8}
	bx lr
test
	PUSH {r1,r2,r3,r4,r5,r6,r7,r8}
	MOV r2,r0
	LDR r7 , =logtable
	LDR r8 , =exptable
	LDRB r3 , [ r7 , r1 ]
	LDRB r4 , [ r7 , r2 ]
	ADD  r3, r4
	LDRB r4 , [ r8 , r3 ]
	RSB  r3 , r1 , #0
	AND  r3 , r2 , r3 , ASR #32
	RSB  r3 , #0
	AND  r4 , r4 , r3,  ASR #32
	MOV  r0,r4
	POP {r1,r2,r3,r4,r5,r6,r7,r8}
	bx lr
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
innerproduct
	;; r0: input vector address, and output
	;; r1: alpha vector address
	PUSH {r3,r4,r5,r6,r7,r8}
	LTORG
	LDR r7 , =logtable
	LTORG
	LDR r8 , =exptable
;;	innerproduct_m c_addr, a_addr, tmp1, tmp2, tmp3, tmp4, res, logtable, exptable, order, interval
	innerproduct_m r0,     r1,     r3,   r4,   r5,   r6,   r2 ,  r7,      r8,       ORDER,     1
	MOV r0, r2
	POP {r3,r4,r5,r6,r7,r8}
	bx lr
innerproduct_inv
	;; r0: input vector address,and output
	;; r1: alpha vector address
	PUSH {r3,r4,r5,r6,r7,r8}
	MOV r9,r2
	LTORG
	LDR r7 , =logtable
	LTORG
	LDR r8 , =exptable
;;	innerproduct_m c_addr, a_addr, tmp1, tmp2, tmp3, tmp4, res, logtable, exptable, order, interval
	innerproduct_m r0,     r1,     r3,   r4,   r5,   r6,   r2,  r7,       r8,       ORDER,     ORDER
	MOV r0, r2
	POP {r3,r4,r5,r6,r7,r8}
	bx lr
tensorproduct
	;; r0: a vector address
	;; r1: b vector address
	;; r2: result address
	PUSH {r3,r4,r5,r6,r7,r8}
	LTORG
	LDR r7 , =logtable
	LTORG
	LDR r8 , =exptable
;;	tensorproduct_m a_addr, b_addr, tmp1, tmp2, tmp3, tmp4, res_addr, logtable, exptable, order
	tensorproduct_m r0,     r1,     r3,   r4,   r5,   r6,   r2 ,      r7,       r8,       ORDER
	POP {r3,r4,r5,r6,r7,r8}
	bx lr
matproduct
	;; r0: input mat address
	;; r1: alpha vector address
	;; r2: output res vector address
	PUSH {r0,r2,r3,r4,r5,r6,r7,r8,r9}
	MOV r9,r2
	LTORG
	LDR r7 , =logtable
	LTORG
	LDR r8 , =exptable
	;;.irp ci,0,1,2,3
	GBLA countii
countii  SETA   0
	WHILE   countii < ORDER
	;;	innerproduct_m c_addr, a_addr, tmp1, tmp2, tmp3, tmp4, res, logtable, exptable, order, interval
		innerproduct_m r0,     r1,     r3,   r4,   r5,   r6,   r2 ,  r7,      r8,       ORDER,    ORDER
		STRB r2, [r9]
		ADD r0,#1
		ADD r9,#1
countii  SETA    countii+1
	WEND
	POP {r0,r2,r3,r4,r5,r6,r7,r8,r9}
	bx lr
matproductalpha
	;; r0: input mat address
	;; r1: alpha vector address
	;; r2: output res vector address
	PUSH {r0,r2,r3,r4,r5,r6,r7,r8,r9}
	MOV r9,r2
	LTORG
	LDR r7 , =logtable
	LTORG
	LDR r8 , =exptable
	;;.irp ci,0,1,2,3
	GBLA countii
countii  SETA   0
	WHILE   countii < ORDER
	;;	innerproduct_m c_addr, a_addr, tmp1, tmp2, tmp3, tmp4, res, logtable, exptable, order, interval
		innerproduct_m r0,     r1,     r3,   r4,   r5,   r6,   r2 ,  r7,      r8,       ORDER,     ORDER
		LDRB r3, [r1,#countii]
		;;mult_exp opA , opB , res , logtable, exptable , tmp
		mult_exp   r3  , r2  , r4  , r7,       r8       , r5
		STRB r4, [r9]
		ADD r0,#1
		ADD r9,#1
countii  SETA    countii+1
	WEND
	POP {r0,r2,r3,r4,r5,r6,r7,r8,r9}
	bx lr	

matproduct_inv
	;; r0: input mat address,and output
	;; r1: alpha vector address
	; r2: output res vector address
	PUSH {r0,r2,r3,r4,r5,r6,r7,r8,r9}
	MOV r9,r2
	LTORG
	LDR r7 , =logtable
	LTORG
	LDR r8 , =exptable
	;;.irp ci,0,1,2,3
	GBLA countii
countii  SETA   0
	WHILE   countii < ORDER
	;;	innerproduct_m c_addr, a_addr, tmp1, tmp2, tmp3, tmp4, res, logtable, exptable, order, interval
		innerproduct_m r0,     r1,     r3,   r4,   r5,   r6,   r2,  r7,       r8,       ORDER    , 1
		STRB r2, [r9]
		ADD r0,#ORDER
		ADD r9,#1
countii  SETA    countii+1
	WEND
	POP {r0,r2,r3,r4,r5,r6,r7,r8,r9}
	bx lr
matproduct_invalpha
	;; r0: input mat address,and output
	;; r1: alpha vector address
	; r2: output res vector address
	PUSH {r0,r2,r3,r4,r5,r6,r7,r8,r9}
	MOV r9,r2
	LTORG
	LDR r7 , =logtable
	LTORG
	LDR r8 , =exptable
	;;.irp ci,0,1,2,3
	GBLA countii
countii  SETA   0
	WHILE   countii < ORDER
	;;	innerproduct_m c_addr, a_addr, tmp1, tmp2, tmp3, tmp4, res, logtable, exptable, order, interval
		innerproduct_m r0,     r1,     r3,   r4,   r5,   r6,   r2,  r7,       r8,       ORDER    , 1
		LDRB r3, [r1,#countii]
		;;mult_exp opA , opB , res , logtable, exptable , tmp
		mult_exp   r3  , r2  , r4  , r7,       r8       , r5
		STRB r4, [r9]
		ADD r0,#ORDER
		ADD r9,#1
countii  SETA    countii+1
	WEND
	POP {r0,r2,r3,r4,r5,r6,r7,r8,r9}
	bx lr	
	
matproduct_oxrsum
	;; r0: input mat address
	;; r1: alpha vector address
	;; r2: oxr vector address
	;; r3: res address
	PUSH {r0,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11}
	MOV r9,r2
	MOV r10,r3
	MOV r11,#0
	LTORG
	LDR r7 , =logtable
	LTORG
	LDR r8 , =exptable
	;;.irp ci,0,1,2,3
	GBLA countii
countii  SETA   0
	WHILE   countii < ORDER
	;;	innerproduct_m c_addr, a_addr, tmp1, tmp2, tmp3, tmp4, res, logtable, exptable, order, interval
		innerproduct_m r0,     r1,     r3,   r4,   r5,   r6,   r2 ,  r7,      r8,      ORDER,      ORDER
		LDRB r3, [r9,#countii]
		EOR r2 , r3
		EOR r11, r2
		ADD r0,#1
countii  SETA    countii+1
	WEND
	STRB r11,[r10]
	POP {r0,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11}
	bx lr
summat
	;; r0: input mat address
	;; r1: output res vector address
	PUSH {r2,r3}
	GBLA countii
	GBLA countjj	
countjj  SETA   0
	WHILE   countjj < ORDER
;;	.irp cj,0,1,2,3
		MOV r3,#0
countii  SETA   0
		WHILE   countii < ORDER
	;;	.irp ci,0,1,2,3
			LDRB r2,[r0,#(countii*ORDER)+countjj]
			EOR r3,r2
countii  SETA    countii+1
		WEND
		STRB r3, [r1,#countjj]
countjj  SETA    countjj+1
	WEND
	POP {r2,r3}
	bx lr
summat_inv
	;; r0: input mat address
	;; r1: output res vector address
	PUSH {r2,r3}
	;;.irp ci,0,1,2,3
countii  SETA   0
	WHILE   countii < ORDER
		MOV r3,#0
countjj  SETA   0
		WHILE   countjj < ORDER
		;;.irp cj,0,1,2,3
			LDRB r2,[r0,#(countii*ORDER)+countjj]
			EOR r3,r2
countjj  SETA    countjj+1
		WEND
		STRB r3, [r1,#countii]
countii  SETA    countii+1
	WEND
	POP {r2,r3}
	bx lr
dotproduct
	;; r0: input and output vector address
	;; r1: alpha vector address
	PUSH {r2,r3,r4,r5,r6,r7}
	LTORG
	LDR r6 , =logtable
	LTORG
	LDR r7 , =exptable
	;;.irp ci,0,1,2,3
countii  SETA   0
	WHILE   countii < ORDER
		LDRB r2, [r0,#countii]
		LDRB r3, [r1,#countii]
	;;	mult_exp opA , opB , res , logtable, exptable , tmp
		mult_exp r2,    r3,  r4,    r6,       r7,        r5
		STRB r4, [r0,#countii]
countii  SETA    countii+1
	WEND
	POP {r2,r3,r4,r5,r6,r7}
	bx lr
dotproduct2
	;; r0: input vector address
	;; r1: alpha vector address
	;; r2: output vector address
	PUSH {r2,r3,r4,r5,r6,r7,r8}
	MOV r8,r2
	LTORG
	LDR r6 , =logtable
	LTORG
	LDR r7 , =exptable
	;;.irp ci,0,1,2,3
countii  SETA   0
	WHILE   countii < ORDER
		LDRB r2, [r0,#countii]
		LDRB r3, [r1,#countii]
	;;	mult_exp opA , opB , res , logtable, exptable , tmp
		mult_exp r2,    r3,  r4,    r6,       r7,        r5
		STRB r4, [r8,#countii]
countii  SETA    countii+1
	WEND
	POP {r2,r3,r4,r5,r6,r7,r8}
	bx lr
spproduct
	;; r0: input1 address (one)
	;; r1: input2 address
	;; r2: output address to xor
	PUSH {r2,r3,r4,r5,r6,r7}
	LTORG
	LDR r6 , =logtable
	LTORG
	LDR r7 , =exptable
	LDRB r3, [r0]
	MOV r0, r3
;;	.irp ci,0,1,2,3
countii  SETA   0
	WHILE   countii < ORDER
		LDRB r3, [r1,#countii]
	;;	mult_exp opA , opB , res , logtable, exptable , tmp
		mult_exp r0,    r3,  r4,    r6,       r7,        r5
		LDRB r5, [r2,#countii]
		EOR  r4, r5
		STRB r4, [r2,#countii]
countii  SETA    countii+1
	WEND
	POP {r2,r3,r4,r5,r6,r7}
	bx lr
sboxac
	;;*sbox_in,*sbox_out,*sbox_t1t2r3t4,*sbox_alpha,*sbox_zs;
	PUSH {r4-r12}
	;;PUSH {r2}
	LDR r6 , =logtable
	LDR r7 , =exptable
	LDR r10, =sbox_t1t2r3t4
	LDR r10, [r10]
	LDR r11, =sbox_alpha
	LDR r11, [r11]
	LDR r12, =sbox_zs
	LDR r12, [r12]
	LDR  r3, =sbox_in
	LDR r3, [r3]
	LDRB r1,[r3]
	;;sbox; input:r1, ouput: r0, sbox_t1t2r3t4:r10, sbox_alpha:r11, sbox_zs:r12 logtable:r6, exptable:r7
	;;unsaved regs: r3,r4,r5,r8
	;;unused reg:r9
	onlinesbox
	LDR r3, =afftable
	;LDR r3, [r3]
	LDRB r0, [r3,r0]
	EOR r0, #99
	LDR r3, =sbox_out
	LDR r3, [r3]
	STRB r0,[r3]
	LTORG
	
	LDR  r3, =sbox_in
	LDR r3, [r3]
	LDRB r1,[r3,#1]
	;;sbox; input:r1, ouput: r0, sbox_t1t2r3t4:r10, sbox_alpha:r11, sbox_zs:r12 logtable:r6, exptable:r7
	;;unsaved regs: r3,r4,r5,r8
	;;unused reg:r9
	onlinesbox
	LDR r3, =afftable
	;LDR r3, [r3]
	LDRB r0, [r3,r0]
	EOR r0, #99
	LDR r3, =sbox_out
	LDR r3, [r3]
	STRB r0,[r3,#1]
	LTORG
	
	LDR  r3, =sbox_in
	LDR r3, [r3]
	LDRB r1,[r3,#2]
	;;sbox; input:r1, ouput: r0, sbox_t1t2r3t4:r10, sbox_alpha:r11, sbox_zs:r12 logtable:r6, exptable:r7
	;;unsaved regs: r3,r4,r5,r8
	;;unused reg:r9
	onlinesbox
	LDR r3, =afftable
	;LDR r3, [r3]
	LDRB r0, [r3,r0]
	EOR r0, #99
	LDR r3, =sbox_out
	LDR r3, [r3]
	STRB r0,[r3,#2]
	LTORG
	
	LDR  r3, =sbox_in
	LDR r3, [r3]
	LDRB r1,[r3,#3]
	;;sbox; input:r1, ouput: r0, sbox_t1t2r3t4:r10, sbox_alpha:r11, sbox_zs:r12 logtable:r6, exptable:r7
	;;unsaved regs: r3,r4,r5,r8
	;;unused reg:r9
	onlinesbox
	LDR r3, =afftable
	;LDR r3, [r3]
	LDRB r0, [r3,r0]
	EOR r0, #99
	LDR r3, =sbox_out
	LDR r3, [r3]
	STRB r0,[r3,#3]
	POP {r4-r12}
	bx lr
	LTORG
	END
	