#include "aes.h"

#define Order 2
#define LenT 5

//#define Order 4
//#define LenT 7

//#define Order 8
//#define LenT 11

//#define Order 16
//#define LenT 19

//#define Order 1
//#define LenT 4

#define LOOP_16(X) X X X X X X X X X X X X X X X X 
#if Order == 4
#define LOOP_Order(X) X X X X
#endif
#if Order == 1
#define LOOP_Order(X) X
#endif
#if Order == 8
#define LOOP_Order(X) X X X X X X X X
#endif
#if Order == 16
#define LOOP_Order(X) X X X X X X X X X X X X X X X X
#endif
#if Order == 2
#define LOOP_Order(X) X X
#endif


extern char gfmul(char a,char b);
extern char innerproduct(char *a,char *b);
extern void tensorproduct(char *a,char *b,char *mat);
extern void matproduct(char *mat,char *alpha,char *ress);
extern void matproduct_inv(char *mat,char *alpha, char *ress);
extern void matproduct_oxrsum(char *mat,char *alpha, char *xor,char *res);
extern void summat(char *mat,char *ress);
extern void summat_inv(char *mat,char *ress);
extern void dotproduct(char *inputoutput,char *alpha);
extern void dotproduct2(char *input,char *alpha,char *res);
extern void spproduct(char *value,char *input,char *ressoxr);
extern char sspproduct();
extern char sboxac();
extern char matproductalpha(char *,char *,char *);
extern char matproduct_invalpha(char *,char *,char *);

extern char power2table[256];
extern char power4table[256];
extern char power16table[256];
extern char afftable[256];
extern char time2table[256];
extern char time3table[256];

//char alpha[128][Order]={0};


//preprecomput
char hzd1[Order],hzd2[Order],hzd3[Order],hzd4[Order];
char addhzd1[Order],addhzd2[Order],addhzd3[Order],addhzd4[Order];
char Smat[Order][Order];
char SRmat[Order][Order];

char onlinepl[16]={0};
char onlinekey[16]={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
char onlinesin[16]={0};
char onlinesout[16]={0};
char onlinemixout[16]={0};

char out[16]={0};


char alpha_add[16*Order]={0};
char alpha_sbox[(16*4)*Order]={0};
char mixalpha[4*3*4*Order]={0};

char addt1t2rt4[16*3]={0};
char t1t2rt4[16*4*LenT]={0};
char mixt1t2rt4[4*4*3*3]={0};

char zs[Order*2]={0};
extern char *sbox_in,*sbox_out,*sbox_t1t2r3t4,*sbox_alpha,*sbox_zs;
char *sbox_in,*sbox_out,*sbox_t1t2r3t4,*sbox_alpha,*sbox_zs;

char vsin[16][Order];
char vsout[16][Order];

char vkey[16][Order];
char vpl[Order]={0};

#if Order==4
char RandR1[Order][Order] = {
{1,33,3,22},
{7,8,2,4},
{52,69,41,54},
{123,45,88,66}
};
char RandR2[Order][Order] = {
{44,53,6,22},
{22,33,34,4},
{52,132,56,135},
{123,232,210,231}
};
char RandR3[Order][Order] = {
{144,153,3,2},
{123,23,13,2},
{23,32,67,78},
{34,121,111,222}
};
char RandR4[Order][Order]= {
{21,32,22,33},
{11,43,221,2},
{23,23,54,78},
{22,121,33,21}
};

char RandRAdd1[Order][Order] = {
{1,33,3,22},
{7,8,2,4},
{52,69,41,54},
{123,45,88,66}
};
char RandRAdd2[Order][Order] = {
{44,53,6,22},
{22,33,34,4},
{52,132,56,135},
{123,232,210,231}
};
char RandRAdd3[Order][Order] = {
{144,153,3,2},
{123,23,13,2},
{23,32,67,78},
{34,121,111,222}
};
char RandRAdd4[Order][Order]= {
{21,32,22,33},
{11,43,221,2},
{23,23,54,78},
{22,121,33,21}
};
#endif
#if Order==1
char RandR1[Order][Order] = {
{2}
};
char RandR2[Order][Order] = {
{3}
};
char RandR3[Order][Order] = {
{4}
};
char RandR4[Order][Order]= {
{5}
};

char RandRAdd1[Order][Order] = {
{2}
};
char RandRAdd2[Order][Order] = {
{3}
};
char RandRAdd3[Order][Order] = {
{4}
};
char RandRAdd4[Order][Order]= {
{5}
};
#endif
#if Order==8
char RandR1[Order][Order];
char RandR2[Order][Order];
char RandR3[Order][Order];
char RandR4[Order][Order];
char RandRAdd1[Order][Order];
char RandRAdd2[Order][Order];
char RandRAdd3[Order][Order];
char RandRAdd4[Order][Order];
#endif
#if Order==16
char RandR1[Order][Order];
char RandR2[Order][Order];
char RandR3[Order][Order];
char RandR4[Order][Order];
char RandRAdd1[Order][Order];
char RandRAdd2[Order][Order];
char RandRAdd3[Order][Order];
char RandRAdd4[Order][Order];
#endif
#if Order==2
char RandR1[Order][Order];
char RandR2[Order][Order];
char RandR3[Order][Order];
char RandR4[Order][Order];
char RandRAdd1[Order][Order];
char RandRAdd2[Order][Order];
char RandRAdd3[Order][Order];
char RandRAdd4[Order][Order];
#endif


char hzd1_t[Order];
char hzd2_t[Order];
char hzd3_t[Order];
char hzd4_t[Order];
char *alpha_t=mixalpha;
char *mixt1t2rt4_t=mixt1t2rt4;
char *t1t2rt4_t;

int preprecomput()
{
	summat_inv(&(RandR1[0][0]),zs);
	summat_inv(&(RandR2[0][0]),zs+Order);
	summat_inv(&(RandR4[0][0]),hzd4);
	
	summat_inv(&(RandRAdd3[0][0]),addhzd3);
	
	summat_inv(&(RandRAdd1[0][0]),addhzd1);
	summat_inv(&(RandRAdd2[0][0]),addhzd2);
	summat_inv(&(RandRAdd4[0][0]),addhzd4);

	tensorproduct(zs,zs+Order,&(Smat[0][0]));
	int i=0;
	int j=0;
	LOOP_Order(i=0;LOOP_Order(SRmat[i][j]=Smat[i][j]^RandR3[i][j];j++;)i++;)
	i=0;LOOP_16(dotproduct2(addhzd3,alpha_add+Order*i,vkey[i]);i++;)
	return 1;
}

int precomput_add()
{
	alpha_t=alpha_add;
	t1t2rt4_t = addt1t2rt4;
	int j,i;
	for (j=0;j<16;j++)
	{
		matproduct_oxrsum(&(RandRAdd1[0][0]),alpha_t,vkey[j],t1t2rt4_t+0);
		matproduct_oxrsum(&(RandRAdd2[0][0]),alpha_t,vpl,t1t2rt4_t+1);
		dotproduct2(addhzd1,alpha_t,hzd1_t);
		dotproduct2(addhzd2,alpha_t,hzd2_t);
		i=0;LOOP_Order(hzd3_t[i]=hzd1_t[i]^hzd2_t[i];i++;);
		matproduct_oxrsum(&(RandRAdd4[0][0]),alpha_t,hzd3_t,t1t2rt4_t+2);
		dotproduct2(addhzd4,alpha_t,vsin[j]);
		t1t2rt4_t = t1t2rt4_t+3;
		alpha_t = alpha_t+Order;
	}
	return 1;
}

int online_add()
{
	int j;
	t1t2rt4_t = addt1t2rt4;
	for (j=0;j<16;j++)
	{
		onlinesin[j] = (onlinepl[j]^(*t1t2rt4_t))^(onlinekey[j]^(*(t1t2rt4_t+1)))^(*(t1t2rt4_t+2));
		t1t2rt4_t = t1t2rt4_t+3;
	}
	return 1;
}


int precomput_sbox(int index)
{
	char *vm1in1;
	char vm1in2[Order];
	char vm1out[Order];
	char *vm2in1;
	char vm2in2[Order];
	char vm2out[Order];
	char *vm3in1;
	char *vm3in2;
	char vm3out[Order];
	char *vm4in1;
	char *vm4in2;
	char vm4out[Order];
	
	int i=0;
	int j=0;
	alpha_t = alpha_sbox+index*Order*4;
	t1t2rt4_t = t1t2rt4+index*LenT*4;
	
	vm1in1 = vsin[index];
	i=0;LOOP_Order(vm1in2[i]=power2table[(unsigned char)vm1in1[i]];i++;)
	
	matproduct_oxrsum(&(RandR1[0][0]),alpha_t,vm1in1,t1t2rt4_t+0);
	matproduct_oxrsum(&(RandR2[0][0]),alpha_t,vm1in2,t1t2rt4_t+1);
	matproduct_invalpha(&(RandR3[0][0]),alpha_t,hzd3);
	matproduct(SRmat[0],alpha_t,t1t2rt4_t+2);
	matproduct_oxrsum(&(RandR4[0][0]),alpha_t,hzd3,t1t2rt4_t+Order+2);
	dotproduct2(hzd4,alpha_t,vm1out);
	alpha_t = alpha_t+Order;
	t1t2rt4_t = t1t2rt4_t+LenT;
	
	vm2in1=vm1out;
	i=0;LOOP_Order(vm2in2[i]=power4table[(unsigned char)vm2in1[i]];i++;);
	
	matproduct_oxrsum(&(RandR1[0][0]),alpha_t,vm2in1,t1t2rt4_t+0);
	matproduct_oxrsum(&(RandR2[0][0]),alpha_t,vm2in2,t1t2rt4_t+1);
	matproduct_invalpha(&(RandR3[0][0]),alpha_t,hzd3);
	matproduct(SRmat[0],alpha_t,t1t2rt4_t+2);
	matproduct_oxrsum(&(RandR4[0][0]),alpha_t,hzd3,t1t2rt4_t+Order+2);
	dotproduct2(hzd4,alpha_t,vm2out);
	alpha_t = alpha_t+Order;
	t1t2rt4_t = t1t2rt4_t+LenT;
	
	vm3in1=vm2out;
	i=0;LOOP_Order(vm3in1[i]=power16table[(unsigned char)vm3in1[i]];i++;);
	vm3in2=vm2in2;
	
	matproduct_oxrsum(&(RandR1[0][0]),alpha_t,vm3in1,t1t2rt4_t+0);
	matproduct_oxrsum(&(RandR2[0][0]),alpha_t,vm3in2,t1t2rt4_t+1);
	matproduct_invalpha(&(RandR3[0][0]),alpha_t,hzd3);
	matproduct(SRmat[0],alpha_t,t1t2rt4_t+2);
	matproduct_oxrsum(&(RandR4[0][0]),alpha_t,hzd3,t1t2rt4_t+2+Order);
	dotproduct2(hzd4,alpha_t,vm3out);
	alpha_t = alpha_t+Order;
	t1t2rt4_t = t1t2rt4_t+LenT;
	
	vm4in1=vm3out;
	vm4in2=vm1in2;
	
	matproduct_oxrsum(&(RandR1[0][0]),alpha_t,vm4in1,t1t2rt4_t+0);
	matproduct_oxrsum(&(RandR2[0][0]),alpha_t,vm4in2,t1t2rt4_t+1);
	matproduct_invalpha(&(RandR3[0][0]),alpha_t,hzd3);
	matproduct(SRmat[0],alpha_t,t1t2rt4_t+2);
	matproduct_oxrsum(&(RandR4[0][0]),alpha_t,hzd3,t1t2rt4_t+2+Order);
	dotproduct2(hzd4,alpha_t,vm4out);
	alpha_t = alpha_t+Order;
	t1t2rt4_t = t1t2rt4_t+LenT;
	i=0;LOOP_Order(vsout[index][i]=afftable[(unsigned char)vm4out[i]];i++;);
	return 1;
}

char vmixin1[4][Order];
char vmixin2[4][Order];
char vmixin3[4][Order];
char vmixin4[4][Order];

int precomput_mix(int index)
{
	int i;
	t1t2rt4_t = mixt1t2rt4+4*3*3*index;
	alpha_t=mixalpha+3*4*Order*index;
	
	i=0;LOOP_Order(vmixin1[0][i]=time2table[(unsigned char)vsout[0+index*4][i]];i++;);
	i=0;LOOP_Order(vmixin1[1][i]=time3table[(unsigned char)vsout[1+index*4][i]];i++;);
	i=0;LOOP_Order(vmixin1[2][i]=vsout[2+index*4][i];i++;);
	i=0;LOOP_Order(vmixin1[3][i]=vsout[3+index*4][i];i++;);
	
	i=0;LOOP_Order(vmixin2[0][i]=vsout[0+index*4][i];i++;);
	i=0;LOOP_Order(vmixin2[1][i]=time2table[(unsigned char)vsout[1+index*4][i]];i++;);
	i=0;LOOP_Order(vmixin2[2][i]=time3table[(unsigned char)vsout[2+index*4][i]];i++;);
	i=0;LOOP_Order(vmixin2[3][i]=vsout[3+index*4][i];i++;);
	
	i=0;LOOP_Order(vmixin3[0][i]=vsout[0+index*4][i];i++;);
	i=0;LOOP_Order(vmixin3[1][i]=vsout[1+index*4][i];i++;);
	i=0;LOOP_Order(vmixin3[2][i]=time2table[(unsigned char)vsout[2+index*4][i]];i++;);
	i=0;LOOP_Order(vmixin3[3][i]=time3table[(unsigned char)vsout[3+index*4][i]];i++;);
	
	i=0;LOOP_Order(vmixin4[0][i]=time3table[(unsigned char)vsout[0+index*4][i]];i++;);
	i=0;LOOP_Order(vmixin4[1][i]=vsout[1+index*4][i];i++;);
	i=0;LOOP_Order(vmixin4[2][i]=vsout[2+index*4][i];i++;);
	i=0;LOOP_Order(vmixin4[3][i]=time2table[(unsigned char)vsout[3+index*4][i]];i++;);
	
	matproduct_oxrsum(&(RandR1[0][0]),alpha_t,vmixin1[0],t1t2rt4_t+0);
	matproduct_oxrsum(&(RandR2[0][0]),alpha_t,vmixin1[1],t1t2rt4_t+1);
	dotproduct2(zs,alpha_t,hzd1_t);
	dotproduct2(zs+Order,alpha_t,hzd2_t);
	i=0;LOOP_Order(hzd3_t[i]=hzd1_t[i]^hzd2_t[i];i++;);
	matproduct_oxrsum(&(RandR4[0][0]),alpha_t,hzd3_t,t1t2rt4_t+2);
	dotproduct2(hzd4,alpha_t,hzd4_t);
	t1t2rt4_t = t1t2rt4_t+3;
	alpha_t = alpha_t+Order;
	
	matproduct_oxrsum(&(RandR1[0][0]),alpha_t,hzd4_t,t1t2rt4_t+0);
	matproduct_oxrsum(&(RandR2[0][0]),alpha_t,vmixin1[2],t1t2rt4_t+1);
	dotproduct2(zs,alpha_t,hzd1_t);
	dotproduct2(zs+Order,alpha_t,hzd2_t);
	i=0;LOOP_Order(hzd3_t[i]=hzd1_t[i]^hzd2_t[i];i++;);
	matproduct_oxrsum(&(RandR4[0][0]),alpha_t,hzd3_t,t1t2rt4_t+2);
	dotproduct2(hzd4,alpha_t,hzd4_t);
	t1t2rt4_t = t1t2rt4_t+3;
	alpha_t = alpha_t+Order;
	
	matproduct_oxrsum(&(RandR1[0][0]),alpha_t,hzd4_t,t1t2rt4_t+0);
	matproduct_oxrsum(&(RandR2[0][0]),alpha_t,vmixin1[3],t1t2rt4_t+1);
	dotproduct2(zs,alpha_t,hzd1_t);
	dotproduct2(zs+Order,alpha_t,hzd2_t);
	i=0;LOOP_Order(hzd3_t[i]=hzd1_t[i]^hzd2_t[i];i++;);
	matproduct_oxrsum(&(RandR4[0][0]),alpha_t,hzd3_t,t1t2rt4_t+2);
	dotproduct2(hzd4,alpha_t,hzd4_t);
	t1t2rt4_t = t1t2rt4_t+3;
	alpha_t = alpha_t+Order;
//////////
	matproduct_oxrsum(&(RandR1[0][0]),alpha_t,vmixin2[0],t1t2rt4_t+0);
	matproduct_oxrsum(&(RandR2[0][0]),alpha_t,vmixin2[1],t1t2rt4_t+1);
	dotproduct2(zs,alpha_t,hzd1_t);
	dotproduct2(zs+Order,alpha_t,hzd2_t);
	i=0;LOOP_Order(hzd3_t[i]=hzd1_t[i]^hzd2_t[i];i++;);
	matproduct_oxrsum(&(RandR4[0][0]),alpha_t,hzd3_t,t1t2rt4_t+2);
	dotproduct2(hzd4,alpha_t,hzd4_t);
	t1t2rt4_t = t1t2rt4_t+3;
	alpha_t = alpha_t+Order;
	
	matproduct_oxrsum(&(RandR1[0][0]),alpha_t,hzd4_t,t1t2rt4_t+0);
	matproduct_oxrsum(&(RandR2[0][0]),alpha_t,vmixin2[2],t1t2rt4_t+1);
	dotproduct2(zs,alpha_t,hzd1_t);
	dotproduct2(zs+Order,alpha_t,hzd2_t);
	i=0;LOOP_Order(hzd3_t[i]=hzd1_t[i]^hzd2_t[i];i++;);
	matproduct_oxrsum(&(RandR4[0][0]),alpha_t,hzd3_t,t1t2rt4_t+2);
	dotproduct2(hzd4,alpha_t,hzd4_t);
	t1t2rt4_t = t1t2rt4_t+3;
	alpha_t = alpha_t+Order;
	
	matproduct_oxrsum(&(RandR1[0][0]),alpha_t,hzd4_t,t1t2rt4_t+0);
	matproduct_oxrsum(&(RandR2[0][0]),alpha_t,vmixin2[3],t1t2rt4_t+1);
	dotproduct2(zs,alpha_t,hzd1_t);
	dotproduct2(zs+Order,alpha_t,hzd2_t);
	i=0;LOOP_Order(hzd3_t[i]=hzd1_t[i]^hzd2_t[i];i++;);
	matproduct_oxrsum(&(RandR4[0][0]),alpha_t,hzd3_t,t1t2rt4_t+2);
	dotproduct2(hzd4,alpha_t,hzd4_t);
	t1t2rt4_t = t1t2rt4_t+3;
	alpha_t = alpha_t+Order;
//////////
	matproduct_oxrsum(&(RandR1[0][0]),alpha_t,vmixin3[0],t1t2rt4_t+0);
	matproduct_oxrsum(&(RandR2[0][0]),alpha_t,vmixin3[1],t1t2rt4_t+1);
	dotproduct2(zs,alpha_t,hzd1_t);
	dotproduct2(zs+Order,alpha_t,hzd2_t);
	i=0;LOOP_Order(hzd3_t[i]=hzd1_t[i]^hzd2_t[i];i++;);
	matproduct_oxrsum(&(RandR4[0][0]),alpha_t,hzd3_t,t1t2rt4_t+2);
	dotproduct2(hzd4,alpha_t,hzd4_t);
	t1t2rt4_t = t1t2rt4_t+3;
	alpha_t = alpha_t+Order;
	
	matproduct_oxrsum(&(RandR1[0][0]),alpha_t,hzd4_t,t1t2rt4_t+0);
	matproduct_oxrsum(&(RandR2[0][0]),alpha_t,vmixin3[2],t1t2rt4_t+1);
	dotproduct2(zs,alpha_t,hzd1_t);
	dotproduct2(zs+Order,alpha_t,hzd2_t);
	i=0;LOOP_Order(hzd3_t[i]=hzd1_t[i]^hzd2_t[i];i++;);
	matproduct_oxrsum(&(RandR4[0][0]),alpha_t,hzd3_t,t1t2rt4_t+2);
	dotproduct2(hzd4,alpha_t,hzd4_t);
	t1t2rt4_t = t1t2rt4_t+3;
	alpha_t = alpha_t+Order;
	
	matproduct_oxrsum(&(RandR1[0][0]),alpha_t,hzd4_t,t1t2rt4_t+0);
	matproduct_oxrsum(&(RandR2[0][0]),alpha_t,vmixin3[3],t1t2rt4_t+1);
	dotproduct2(zs,alpha_t,hzd1_t);
	dotproduct2(zs+Order,alpha_t,hzd2_t);
	i=0;LOOP_Order(hzd3_t[i]=hzd1_t[i]^hzd2_t[i];i++;);
	matproduct_oxrsum(&(RandR4[0][0]),alpha_t,hzd3_t,t1t2rt4_t+2);
	dotproduct2(hzd4,alpha_t,hzd4_t);
	t1t2rt4_t = t1t2rt4_t+3;
	alpha_t = alpha_t+Order;
//////////
	matproduct_oxrsum(&(RandR1[0][0]),alpha_t,vmixin4[0],t1t2rt4_t+0);
	matproduct_oxrsum(&(RandR2[0][0]),alpha_t,vmixin4[1],t1t2rt4_t+1);
	dotproduct2(zs,alpha_t,hzd1_t);
	dotproduct2(zs+Order,alpha_t,hzd2_t);
	i=0;LOOP_Order(hzd3_t[i]=hzd1_t[i]^hzd2_t[i];i++;);
	matproduct_oxrsum(&(RandR4[0][0]),alpha_t,hzd3_t,t1t2rt4_t+2);
	dotproduct2(hzd4,alpha_t,hzd4_t);
	t1t2rt4_t = t1t2rt4_t+3;
	alpha_t = alpha_t+Order;
	
	matproduct_oxrsum(&(RandR1[0][0]),alpha_t,hzd4_t,t1t2rt4_t+0);
	matproduct_oxrsum(&(RandR2[0][0]),alpha_t,vmixin4[2],t1t2rt4_t+1);
	dotproduct2(zs,alpha_t,hzd1_t);
	dotproduct2(zs+Order,alpha_t,hzd2_t);
	i=0;LOOP_Order(hzd3_t[i]=hzd1_t[i]^hzd2_t[i];i++;);
	matproduct_oxrsum(&(RandR4[0][0]),alpha_t,hzd3_t,t1t2rt4_t+2);
	dotproduct2(hzd4,alpha_t,hzd4_t);
	t1t2rt4_t = t1t2rt4_t+3;
	alpha_t = alpha_t+Order;
	
	matproduct_oxrsum(&(RandR1[0][0]),alpha_t,hzd4_t,t1t2rt4_t+0);
	matproduct_oxrsum(&(RandR2[0][0]),alpha_t,vmixin4[3],t1t2rt4_t+1);
	dotproduct2(zs,alpha_t,hzd1_t);
	dotproduct2(zs+Order,alpha_t,hzd2_t);
	i=0;LOOP_Order(hzd3_t[i]=hzd1_t[i]^hzd2_t[i];i++;);
	matproduct_oxrsum(&(RandR4[0][0]),alpha_t,hzd3_t,t1t2rt4_t+2);
	dotproduct2(hzd4,alpha_t,hzd4_t);
	t1t2rt4_t = t1t2rt4_t+3;
	alpha_t = alpha_t+Order;
	
	return 1;
}


void onlinecomput_sbox(int index)
{
	sbox_in=onlinesin+4*index;
	sbox_out=onlinesout+4*index;
	sbox_t1t2r3t4=t1t2rt4+4*4*index*LenT;
	sbox_alpha=alpha_sbox+4*4*Order*index;
	sbox_zs=zs;
//	
	sboxac();
}



int onine_mix(int index)
{
	char mixin[4][4];
	char temp;
//	char mixout[4];
//	char *mt=mixt1t2rt4;
	
	t1t2rt4_t=mixt1t2rt4+index*36;
	
	mixin[0][0]=time2table[(unsigned char)onlinesout[0+index*4]];
	mixin[0][1]=time3table[(unsigned char)onlinesout[1+index*4]];
	mixin[0][2]=onlinesout[2+index*4];
	mixin[0][3]=onlinesout[3+index*4];
	
	mixin[1][0]=onlinesout[0+index*4];
	mixin[1][1]=time2table[(unsigned char)onlinesout[1+index*4]];
	mixin[1][2]=time3table[(unsigned char)onlinesout[2+index*4]];
	mixin[1][3]=onlinesout[3+index*4];
	
	mixin[2][0]=onlinesout[0+index*4];
	mixin[2][1]=onlinesout[1+index*4];
	mixin[2][2]=time2table[(unsigned char)(unsigned char)onlinesout[2+index*4]];
	mixin[2][3]=time3table[(unsigned char)onlinesout[3+index*4]];
	
	mixin[3][0]=time3table[(unsigned char)onlinesout[0+index*4]];
	mixin[3][1]=onlinesout[1+index*4];
	mixin[3][2]=onlinesout[2+index*4];
	mixin[3][3]=time2table[(unsigned char)onlinesout[3+index*4]];
	
	
	temp = (mixin[0][0]^(*(t1t2rt4_t)))^(mixin[0][1]^(*(t1t2rt4_t+1)))^(*(t1t2rt4_t+2));
	t1t2rt4_t=t1t2rt4_t+3;
	temp = (temp^(*(t1t2rt4_t)))^(mixin[0][2]^(*(t1t2rt4_t+1)))^(*(t1t2rt4_t+2));
	t1t2rt4_t=t1t2rt4_t+3;
	onlinemixout[0+index*4] = (temp^(*(t1t2rt4_t)))^(mixin[0][3]^(*(t1t2rt4_t+1)))^(*(t1t2rt4_t+2));
	t1t2rt4_t=t1t2rt4_t+3;

	temp = (mixin[1][0]^(*(t1t2rt4_t)))^(mixin[1][1]^(*(t1t2rt4_t+1)))^(*(t1t2rt4_t+2));
	t1t2rt4_t+=3;
	temp = (temp^(*(t1t2rt4_t)))^(mixin[1][2]^(*(t1t2rt4_t+1)))^(*(t1t2rt4_t+2));
	t1t2rt4_t+=3;
	onlinemixout[1+index*4] = (temp^(*(t1t2rt4_t)))^(mixin[1][3]^(*(t1t2rt4_t+1)))^(*(t1t2rt4_t+2));
	t1t2rt4_t+=3;
	
	temp = (mixin[2][0]^(*(t1t2rt4_t)))^(mixin[2][1]^(*(t1t2rt4_t+1)))^(*(t1t2rt4_t+2));
	t1t2rt4_t+=3;
	temp = (temp^(*(t1t2rt4_t)))^(mixin[2][2]^(*(t1t2rt4_t+1)))^(*(t1t2rt4_t+2));
	t1t2rt4_t+=3;
	onlinemixout[2+index*4] = (temp^(*(t1t2rt4_t)))^(mixin[2][3]^(*(t1t2rt4_t+1)))^(*(t1t2rt4_t+2));
	t1t2rt4_t+=3;
	
	temp = (mixin[3][0]^(*(t1t2rt4_t)))^(mixin[3][1]^(*(t1t2rt4_t+1)))^(*(t1t2rt4_t+2));
	t1t2rt4_t+=3;
	temp = (temp^(*(t1t2rt4_t)))^(mixin[3][2]^(*(t1t2rt4_t+1)))^(*(t1t2rt4_t+2));
	t1t2rt4_t+=3;
	onlinemixout[3+index*4] = (temp^(*(t1t2rt4_t)))^(mixin[3][3]^(*(t1t2rt4_t+1)))^(*(t1t2rt4_t+2));
	t1t2rt4_t+=3;
	
	return 1;
}

int prealpha()
{
	int i,j;
	for (i=0;i<16*Order;i++)
		alpha_add[i] = i%255+1;
	for (i=0;i<(16*4)*Order;i++)
		alpha_sbox[i] = i%255+1;
	j=0;
	for (i=4*3*4*Order-1;i>=0;i--)
	{
		mixalpha[j]=i%255+1;
		j++;
	}
#if Order==1
	alpha_sbox[0] = 51;
	alpha_sbox[1] = 177;
	alpha_sbox[2] = 15;
	alpha_sbox[3] = 84;
	alpha_sbox[4] = 178;
	alpha_sbox[5] = 202;
	alpha_sbox[6] = 143;
	alpha_sbox[7] = 166;
	alpha_sbox[8] = 79;
	alpha_sbox[9] = 110;
	alpha_sbox[10] = 201;
	alpha_sbox[11] = 224;
	alpha_sbox[12] = 44;
	alpha_sbox[13] = 76;
	alpha_sbox[14] = 94;
	alpha_sbox[15] = 57;
	alpha_sbox[16] = 149;
	alpha_sbox[17] = 168;
	alpha_sbox[18] = 244;
	alpha_sbox[19] = 78;
	alpha_sbox[20] = 131;
	alpha_sbox[21] = 199;
	alpha_sbox[22] = 198;
	alpha_sbox[23] = 138;
	alpha_sbox[24] = 248;
	alpha_sbox[25] = 220;
	alpha_sbox[26] = 217;
	alpha_sbox[27] = 167;
	alpha_sbox[28] = 31;
	alpha_sbox[29] = 54;
	alpha_sbox[30] = 225;
	alpha_sbox[31] = 41;
	alpha_sbox[32] = 232;
	alpha_sbox[33] = 52;
	alpha_sbox[34] = 59;
	alpha_sbox[35] = 70;
	alpha_sbox[36] = 133;
	alpha_sbox[37] = 97;
	alpha_sbox[38] = 135;
	alpha_sbox[39] = 173;
	alpha_sbox[40] = 35;
	alpha_sbox[41] = 161;
	alpha_sbox[42] = 154;
	alpha_sbox[43] = 181;
	alpha_sbox[44] = 23;
	alpha_sbox[45] = 77;
	alpha_sbox[46] = 112;
	alpha_sbox[47] = 69;
	alpha_sbox[48] = 165;
	alpha_sbox[49] = 176;
	alpha_sbox[50] = 21;
	alpha_sbox[51] = 87;
	alpha_sbox[52] = 22;
	alpha_sbox[53] = 88;
	alpha_sbox[54] = 7;
	alpha_sbox[55] = 38;
	alpha_sbox[56] = 71;
	alpha_sbox[57] = 120;
	alpha_sbox[58] = 121;
	alpha_sbox[59] = 194;
	alpha_sbox[60] = 210;
	alpha_sbox[61] = 164;
	alpha_sbox[62] = 234;
	alpha_sbox[63] = 150;
	alpha_sbox[64] = 13;
	alpha_sbox[65] = 107;
	alpha_sbox[66] = 37;
	alpha_sbox[67] = 42;
	alpha_sbox[68] = 28;
	alpha_sbox[69] = 46;
	alpha_sbox[70] = 139;
	alpha_sbox[71] = 141;
	alpha_sbox[72] = 195;
	alpha_sbox[73] = 27;
	alpha_sbox[74] = 30;
	alpha_sbox[75] = 11;
	alpha_sbox[76] = 103;
	alpha_sbox[77] = 43;
	alpha_sbox[78] = 29;
	alpha_sbox[79] = 25;
	alpha_sbox[80] = 205;
	alpha_sbox[81] = 208;
	alpha_sbox[82] = 73;
	alpha_sbox[83] = 49;
	alpha_sbox[84] = 50;
	alpha_sbox[85] = 171;
	alpha_sbox[86] = 124;
	alpha_sbox[87] = 186;
	alpha_sbox[88] = 229;
	alpha_sbox[89] = 81;
	alpha_sbox[90] = 62;
	alpha_sbox[91] = 146;
	alpha_sbox[92] = 98;
	alpha_sbox[93] = 85;
	alpha_sbox[94] = 151;
	alpha_sbox[95] = 227;
	alpha_sbox[96] = 241;
	alpha_sbox[97] = 75;
	alpha_sbox[98] = 203;
	alpha_sbox[99] = 104;
	alpha_sbox[100] = 108;
	alpha_sbox[101] = 156;
	alpha_sbox[102] = 200;
	alpha_sbox[103] = 140;
	alpha_sbox[104] = 212;
	alpha_sbox[105] = 114;
	alpha_sbox[106] = 89;
	alpha_sbox[107] = 158;
	alpha_sbox[108] = 74;
	alpha_sbox[109] = 82;
	alpha_sbox[110] = 14;
	alpha_sbox[111] = 47;
	alpha_sbox[112] = 90;
	alpha_sbox[113] = 39;
	alpha_sbox[114] = 196;
	alpha_sbox[115] = 204;
	alpha_sbox[116] = 213;
	alpha_sbox[117] = 182;
	alpha_sbox[118] = 60;
	alpha_sbox[119] = 26;
	alpha_sbox[120] = 45;
	alpha_sbox[121] = 216;
	alpha_sbox[122] = 67;
	alpha_sbox[123] = 157;
	alpha_sbox[124] = 179;
	alpha_sbox[125] = 155;
	alpha_sbox[126] = 163;
	alpha_sbox[127] = 206;
	mixalpha[0] = 169;
	mixalpha[1] = 58;
	mixalpha[2] = 170;
	mixalpha[3] = 134;
	mixalpha[4] = 197;
	mixalpha[5] = 148;
	mixalpha[6] = 101;
	mixalpha[7] = 152;
	mixalpha[8] = 86;
	mixalpha[9] = 193;
	mixalpha[10] = 118;
	mixalpha[11] = 109;
	mixalpha[12] = 122;
	mixalpha[13] = 218;
	mixalpha[14] = 105;
	mixalpha[15] = 102;
	mixalpha[16] = 162;
	mixalpha[17] = 180;
	mixalpha[18] = 116;
	mixalpha[19] = 92;
	mixalpha[20] = 100;
	mixalpha[21] = 242;
	mixalpha[22] = 106;
	mixalpha[23] = 236;
	mixalpha[24] = 117;
	mixalpha[25] = 91;
	mixalpha[26] = 233;
	mixalpha[27] = 185;
	mixalpha[28] = 61;
	mixalpha[29] = 174;
	mixalpha[30] = 55;
	mixalpha[31] = 211;
	mixalpha[32] = 153;
	mixalpha[33] = 137;
	mixalpha[34] = 56;
	mixalpha[35] = 115;
	mixalpha[36] = 93;
	mixalpha[37] = 188;
	mixalpha[38] = 209;
	mixalpha[39] = 113;
	mixalpha[40] = 228;
	mixalpha[41] = 240;
	mixalpha[42] = 19;
	mixalpha[43] = 99;
	mixalpha[44] = 184;
	mixalpha[45] = 147;
	mixalpha[46] = 214;
	mixalpha[47] = 226;
	mixalpha[48] = 53;
	mixalpha[49] = 172;
	mixalpha[50] = 142;
	mixalpha[51] = 145;
	mixalpha[52] = 83;
	mixalpha[53] = 230;
	mixalpha[54] = 212;
	mixalpha[55] = 188;
	mixalpha[56] = 74;
	mixalpha[57] = 67;
	mixalpha[58] = 143;
	mixalpha[59] = 81;
	mixalpha[60] = 39;
	mixalpha[61] = 164;
	mixalpha[62] = 225;
	mixalpha[63] = 205;
	mixalpha[64] = 214;
	mixalpha[65] = 83;
	mixalpha[66] = 91;
	mixalpha[67] = 163;
	mixalpha[68] = 77;
	mixalpha[69] = 51;
	mixalpha[70] = 98;
	mixalpha[71] = 35;
	mixalpha[72] = 177;
	mixalpha[73] = 202;
	mixalpha[74] = 142;
	mixalpha[75] = 185;
	mixalpha[76] = 121;
	mixalpha[77] = 148;
	mixalpha[78] = 178;
	mixalpha[79] = 106;
	mixalpha[80] = 15;
	mixalpha[81] = 195;
	mixalpha[82] = 172;
	mixalpha[83] = 174;
	mixalpha[84] = 213;
	mixalpha[85] = 186;
	mixalpha[86] = 113;
	mixalpha[87] = 139;
	mixalpha[88] = 152;
	mixalpha[89] = 197;
	mixalpha[90] = 244;
	mixalpha[91] = 181;
	mixalpha[92] = 43;
	mixalpha[93] = 199;
	mixalpha[94] = 133;
	mixalpha[95] = 62;
#endif
#if Order==2
	alpha_sbox[0] = 51;
	alpha_sbox[1] = 177;
	alpha_sbox[2] = 15;
	alpha_sbox[3] = 84;
	alpha_sbox[4] = 178;
	alpha_sbox[5] = 202;
	alpha_sbox[6] = 143;
	alpha_sbox[7] = 166;
	alpha_sbox[8] = 79;
	alpha_sbox[9] = 110;
	alpha_sbox[10] = 201;
	alpha_sbox[11] = 224;
	alpha_sbox[12] = 44;
	alpha_sbox[13] = 76;
	alpha_sbox[14] = 94;
	alpha_sbox[15] = 57;
	alpha_sbox[16] = 149;
	alpha_sbox[17] = 168;
	alpha_sbox[18] = 244;
	alpha_sbox[19] = 78;
	alpha_sbox[20] = 131;
	alpha_sbox[21] = 199;
	alpha_sbox[22] = 198;
	alpha_sbox[23] = 138;
	alpha_sbox[24] = 248;
	alpha_sbox[25] = 220;
	alpha_sbox[26] = 217;
	alpha_sbox[27] = 167;
	alpha_sbox[28] = 31;
	alpha_sbox[29] = 54;
	alpha_sbox[30] = 225;
	alpha_sbox[31] = 41;
	alpha_sbox[32] = 232;
	alpha_sbox[33] = 52;
	alpha_sbox[34] = 59;
	alpha_sbox[35] = 70;
	alpha_sbox[36] = 133;
	alpha_sbox[37] = 97;
	alpha_sbox[38] = 135;
	alpha_sbox[39] = 173;
	alpha_sbox[40] = 35;
	alpha_sbox[41] = 161;
	alpha_sbox[42] = 154;
	alpha_sbox[43] = 181;
	alpha_sbox[44] = 23;
	alpha_sbox[45] = 77;
	alpha_sbox[46] = 112;
	alpha_sbox[47] = 69;
	alpha_sbox[48] = 165;
	alpha_sbox[49] = 176;
	alpha_sbox[50] = 21;
	alpha_sbox[51] = 87;
	alpha_sbox[52] = 22;
	alpha_sbox[53] = 88;
	alpha_sbox[54] = 7;
	alpha_sbox[55] = 38;
	alpha_sbox[56] = 71;
	alpha_sbox[57] = 120;
	alpha_sbox[58] = 121;
	alpha_sbox[59] = 194;
	alpha_sbox[60] = 210;
	alpha_sbox[61] = 164;
	alpha_sbox[62] = 234;
	alpha_sbox[63] = 150;
	alpha_sbox[64] = 13;
	alpha_sbox[65] = 107;
	alpha_sbox[66] = 37;
	alpha_sbox[67] = 42;
	alpha_sbox[68] = 28;
	alpha_sbox[69] = 46;
	alpha_sbox[70] = 139;
	alpha_sbox[71] = 141;
	alpha_sbox[72] = 195;
	alpha_sbox[73] = 27;
	alpha_sbox[74] = 30;
	alpha_sbox[75] = 11;
	alpha_sbox[76] = 103;
	alpha_sbox[77] = 43;
	alpha_sbox[78] = 29;
	alpha_sbox[79] = 25;
	alpha_sbox[80] = 205;
	alpha_sbox[81] = 208;
	alpha_sbox[82] = 73;
	alpha_sbox[83] = 49;
	alpha_sbox[84] = 50;
	alpha_sbox[85] = 171;
	alpha_sbox[86] = 124;
	alpha_sbox[87] = 186;
	alpha_sbox[88] = 229;
	alpha_sbox[89] = 81;
	alpha_sbox[90] = 62;
	alpha_sbox[91] = 146;
	alpha_sbox[92] = 98;
	alpha_sbox[93] = 85;
	alpha_sbox[94] = 151;
	alpha_sbox[95] = 227;
	alpha_sbox[96] = 241;
	alpha_sbox[97] = 75;
	alpha_sbox[98] = 203;
	alpha_sbox[99] = 104;
	alpha_sbox[100] = 108;
	alpha_sbox[101] = 156;
	alpha_sbox[102] = 200;
	alpha_sbox[103] = 140;
	alpha_sbox[104] = 212;
	alpha_sbox[105] = 114;
	alpha_sbox[106] = 89;
	alpha_sbox[107] = 158;
	alpha_sbox[108] = 74;
	alpha_sbox[109] = 82;
	alpha_sbox[110] = 14;
	alpha_sbox[111] = 47;
	alpha_sbox[112] = 90;
	alpha_sbox[113] = 39;
	alpha_sbox[114] = 196;
	alpha_sbox[115] = 204;
	alpha_sbox[116] = 213;
	alpha_sbox[117] = 182;
	alpha_sbox[118] = 60;
	alpha_sbox[119] = 26;
	alpha_sbox[120] = 45;
	alpha_sbox[121] = 216;
	alpha_sbox[122] = 67;
	alpha_sbox[123] = 157;
	alpha_sbox[124] = 179;
	alpha_sbox[125] = 155;
	alpha_sbox[126] = 163;
	alpha_sbox[127] = 206;
	mixalpha[0] = 169;
	mixalpha[1] = 58;
	mixalpha[2] = 170;
	mixalpha[3] = 134;
	mixalpha[4] = 197;
	mixalpha[5] = 148;
	mixalpha[6] = 101;
	mixalpha[7] = 152;
	mixalpha[8] = 86;
	mixalpha[9] = 193;
	mixalpha[10] = 118;
	mixalpha[11] = 109;
	mixalpha[12] = 122;
	mixalpha[13] = 218;
	mixalpha[14] = 105;
	mixalpha[15] = 102;
	mixalpha[16] = 162;
	mixalpha[17] = 180;
	mixalpha[18] = 116;
	mixalpha[19] = 92;
	mixalpha[20] = 100;
	mixalpha[21] = 242;
	mixalpha[22] = 106;
	mixalpha[23] = 236;
	mixalpha[24] = 117;
	mixalpha[25] = 91;
	mixalpha[26] = 233;
	mixalpha[27] = 185;
	mixalpha[28] = 61;
	mixalpha[29] = 174;
	mixalpha[30] = 55;
	mixalpha[31] = 211;
	mixalpha[32] = 153;
	mixalpha[33] = 137;
	mixalpha[34] = 56;
	mixalpha[35] = 115;
	mixalpha[36] = 93;
	mixalpha[37] = 188;
	mixalpha[38] = 209;
	mixalpha[39] = 113;
	mixalpha[40] = 228;
	mixalpha[41] = 240;
	mixalpha[42] = 19;
	mixalpha[43] = 99;
	mixalpha[44] = 184;
	mixalpha[45] = 147;
	mixalpha[46] = 214;
	mixalpha[47] = 226;
	mixalpha[48] = 53;
	mixalpha[49] = 172;
	mixalpha[50] = 142;
	mixalpha[51] = 145;
	mixalpha[52] = 83;
	mixalpha[53] = 230;
	mixalpha[54] = 212;
	mixalpha[55] = 188;
	mixalpha[56] = 74;
	mixalpha[57] = 67;
	mixalpha[58] = 143;
	mixalpha[59] = 81;
	mixalpha[60] = 39;
	mixalpha[61] = 164;
	mixalpha[62] = 225;
	mixalpha[63] = 205;
	mixalpha[64] = 214;
	mixalpha[65] = 83;
	mixalpha[66] = 91;
	mixalpha[67] = 163;
	mixalpha[68] = 77;
	mixalpha[69] = 51;
	mixalpha[70] = 98;
	mixalpha[71] = 35;
	mixalpha[72] = 177;
	mixalpha[73] = 202;
	mixalpha[74] = 142;
	mixalpha[75] = 185;
	mixalpha[76] = 121;
	mixalpha[77] = 148;
	mixalpha[78] = 178;
	mixalpha[79] = 106;
	mixalpha[80] = 15;
	mixalpha[81] = 195;
	mixalpha[82] = 172;
	mixalpha[83] = 174;
	mixalpha[84] = 213;
	mixalpha[85] = 186;
	mixalpha[86] = 113;
	mixalpha[87] = 139;
	mixalpha[88] = 152;
	mixalpha[89] = 197;
	mixalpha[90] = 244;
	mixalpha[91] = 181;
	mixalpha[92] = 43;
	mixalpha[93] = 199;
	mixalpha[94] = 133;
	mixalpha[95] = 62;
#endif
	
	
	return 1;
}

struct AES_ctx AES_ctx1;
unsigned char aesbuf[16];

int main(void)
{
	// test for unprotected AES
	AES_ECB_encrypt(&AES_ctx1, aesbuf);
	
	prealpha();
	
	// test for protected AES
	preprecomput();

	precomput_add();
	
	precomput_sbox(0);
	precomput_sbox(1);
	precomput_sbox(2);
	precomput_sbox(3);
	precomput_sbox(4);
	precomput_sbox(5);
	precomput_sbox(6);
	precomput_sbox(7);
	precomput_sbox(8);
	precomput_sbox(9);
	precomput_sbox(10);
	precomput_sbox(11);
	precomput_sbox(12);
	precomput_sbox(13);
	precomput_sbox(14);
	precomput_sbox(15);
	
	precomput_mix(0);
	precomput_mix(1);
	precomput_mix(2);
	precomput_mix(3);
	
	online_add();
	
	onlinecomput_sbox(0);
	onlinecomput_sbox(1);
	onlinecomput_sbox(2);
	onlinecomput_sbox(3);
	
	onine_mix(0);
	onine_mix(1);
	onine_mix(2);
	onine_mix(3);

	return 1;
}





















