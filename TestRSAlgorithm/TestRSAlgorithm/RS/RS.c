//
//  RS.c
//  TestCProjectOnXCode
//
//  Created by Mark on 12-11-19.
//  Copyright (c) 2012年 Mark Chen. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>

// 此处引入string.h头文件，是因为VC无法识别memset函数(在Xcode中此处标注警告)
#include <string.h>
#include "RS.h"

// 为满足VC编译器要求，需要对使用到的函数进行声明
int read_p(RSConfig *);
int generate_gf(RSConfig *);

/**
 *  Init RS Config.
 *
 *  Return: 0 is okay, otherwise is error.
 */
int initRS(RSConfig *rsConfig, int aM, int aLength, int aK, int aR){
    int status = 0;
    if(rsConfig == NULL){
        printf("Error : initRS() - rsConfig is NULL");
        status = -1;
    }else{
        rsConfig->m = aM;
        rsConfig->length = aLength;
        rsConfig->k = aK;
        rsConfig->r = aR;
        rsConfig->n = (1<<aM)-1;
        
        if(rsConfig->alpha_to == NULL){
            rsConfig->alpha_to = (int *)malloc(sizeof(int)*(rsConfig->n + 1));
            if (rsConfig->alpha_to == NULL) {
                status = -2;
            }else{
                memset(rsConfig->alpha_to, 0, sizeof(int) * (rsConfig->n + 1));
            }
        }
        
        if (rsConfig->index_of == NULL) {
            rsConfig->index_of = (int *)malloc(sizeof(int) * (rsConfig->n + 1));
            if (rsConfig->index_of == NULL) {
                status = -3;
            }else{
                memset((rsConfig->index_of), 0, sizeof(int) * (rsConfig->n + 1));
            }
        }
        
        if (rsConfig->g == NULL) {
            rsConfig->g = (int *)malloc(sizeof(int) * (rsConfig->r + 1));
            if (rsConfig->g == NULL) {
                status = -4;
            }else{
                memset((rsConfig->g), 0, sizeof(int) * (rsConfig->r + 1));
            }
        }
        
        if (rsConfig->recd == NULL) {
            rsConfig->recd = (int *)malloc(sizeof(int) * (rsConfig->length));
            if (rsConfig->recd == NULL) {
                status = -5;
            }else{
                memset((rsConfig->recd), 0, sizeof(int) * (rsConfig->length));
            }
        }
        
        if(rsConfig->data == NULL){
            rsConfig->data = (int *)malloc(sizeof(int) * (rsConfig->k));
            if (rsConfig->data == NULL) {
                status = -6;
            }else{
                memset((rsConfig->data), 0, sizeof(int) * (rsConfig->k));
            }
        }
        
        if(rsConfig->rs_b == NULL){
            rsConfig->rs_b = (int *)malloc(sizeof(int) * (rsConfig->r));
            if (rsConfig->rs_b == NULL) {
                status = -7;
            }else{
                memset((rsConfig->rs_b), 0, sizeof(int) * rsConfig->r);
            }
        }
        
        rsConfig->numera = 0;
        
        if(rsConfig->era == NULL){
            rsConfig->era = (int *)malloc(sizeof(int) * (rsConfig->r));
            if (rsConfig->era == NULL) {
                status = -8;
            }else{
                memset(rsConfig->era, 0, sizeof(int) * rsConfig->r);
            }
        }
        
        
    }
    return status;
}

/**
 *  RS编码，通过RSConfig编码并添加响应长度纠错码，纠错码长为r.
 */
int encodeWithRS(RSConfig *rsConfig){
    if (rsConfig == NULL) {
        LOGD("\nError : RS - encodeWithRS - rsConfig = NULL\n");
        return -1;
    }
    
    int status = 0;     // return Value.
    int i;
    
    read_p(rsConfig);        // Read m */
    generate_gf(rsConfig);
    
    printf("\n--> This is an RS(%d,%d,%d) code over GF(2^%d), t = %d", rsConfig->length, rsConfig->k, rsConfig->r/*red+1*/, rsConfig->m, rsConfig->r/2);
    
    printf("    with zeros: ");
    for (i=0;i<rsConfig->r;i++){
        printf("%d ", init_zero+i);
    }
    
    printf("\n\n");
    
    gen_poly(rsConfig);
    encode_rs(rsConfig);
    
    
    // 将校验数据头和encode_rs的输出数据拼接装入recd数组中输出，这才是最终编码的结果；
    // Codeword(编码) is c(X) = data(X)*X**(length-k)+ rs_b(X)
    for (i = 0; i < rsConfig->length - rsConfig->k; i++){       // Add the header data to recd array.
        rsConfig->recd[i] = rsConfig->rs_b[i];
    }
    for (i = 0; i < rsConfig->k; i++){                // Add the message data to recd array.
        rsConfig->recd[i + rsConfig->length - rsConfig->k] = rsConfig->data[i];
    }
    
    //    // 这里为了测试制造一些麻烦, recd就是编码输出。
    //    rsConfig->recd[3] = 'E';
    
    return status;
}


int prepareDecode(RSConfig *rsConfig){
    if (rsConfig == NULL) {
        LOGD("\nError : RS - prepareDecode - rsConfig = NULL\n");
        return -1;
    }
    
    int status = 0;
    int i = 0;
    
    
    // Added
    status = read_p(rsConfig);        // Read m */
    if (status < 0) {
        LOGD("\nprepareDecode - read_p - FAIL! status = %d\n", status);
        return status;
    }
    
    status = generate_gf(rsConfig);
    if (status < 0) {
        LOGD("\nprepareDecode - generate_gf - FAIL! status = %d\n", status);
        return status;
    }
    
    // End, Added.
    printf("\n--> This is an RS(%d,%d,%d) code over GF(2^%d), t = %d", rsConfig->length, rsConfig->k, rsConfig->r/*red+1*/, rsConfig->m, rsConfig->r/2);
    
    printf("    with zeros: ");
    for (i=0;i<rsConfig->r;i++){
        printf("%d ", init_zero+i);
    }
    
    printf("\n\n");
    
    status = gen_poly(rsConfig);
    if (status < 0) {
        LOGD("\nprepareDecode - gen_poly - FAIL! status = %d\n", status);
        return status;
    }
    
    return status;
}

/**
 *  RS解码，通过RSConfig解码，并进行纠错
 *
 *  Return value:
 *  0   - Solve Success.
 *  <0  - Cannot Solve.
 */
int decodeWithRS(struct RSConfig *rsConfig){
    if (rsConfig == NULL) {
        LOGD("\nError : RS - decodeWithRS - rsConfig = NULL\n");
        return -1;
    }
    int status = 0;     // return Value.
    
    // Decode the RS input code.
    status = decode_rs(rsConfig);
    
    return status;
}


/**
 *  Clean up the memory.
 */
int cleanup(RSConfig *rsConfig){
    if (rsConfig == NULL) {
        LOGD("\nError : RS - generate_gf - rsConfig = NULL\n");
        return -1;
    }
    
    int status = 0;
    
    if(rsConfig->alpha_to != NULL){
        free(rsConfig->alpha_to);
        rsConfig->alpha_to = NULL;
    }
    
    if(rsConfig->index_of != NULL){
        free(rsConfig->index_of);
        rsConfig->index_of = NULL;
    }
    
    if(rsConfig->g != NULL){
        free(rsConfig->g);
        rsConfig->g = NULL;
    }
    
    if (rsConfig->recd != NULL) {
        free(rsConfig->recd);
        rsConfig->recd = NULL;
    }
    
    if(rsConfig->data != NULL){
        free(rsConfig->data);
        rsConfig->data = NULL;
    }
    
    if(rsConfig->rs_b != NULL){
        free(rsConfig->rs_b);
        rsConfig->rs_b = NULL;
    }
    
    if(rsConfig->era != NULL){
        free(rsConfig->era);
        rsConfig->era = NULL;
    }
    
    return status;
}

/**
 *  Read m, the degree of a primitive polynomial p(x) used to
 *  compute the Galois field GF(2**m).
 *  Get precomputed coefficients p[] of p(x). Read  the code length.
 */
int read_p(RSConfig *rsConfig) {
    int status = 0;
    if(rsConfig == NULL){
        LOGD("\nError : RS - read_p - rsConfig = NULL\n");
        return -1;
    }
    
    int m = rsConfig->m;
    int i; //, ninf;
    
    for (i=1; i<m; i++){
        rsConfig->p[i] = 0;
    }
    rsConfig->p[0] = rsConfig->p[m] = 1;      //LSB First
    if (m == 2){
        rsConfig->p[1] = 1;  //111
    }else if (m == 3){
        rsConfig->p[1] = 1;  //1101
    }else if (m == 4){
        rsConfig->p[3] = 1;  //11001
    }// else if (m == 4)        p[1] = 1;  // Commented out to match example p. 68
    else if (m == 5){
        rsConfig->p[2] = 1;  //101001
    }else if (m == 6){
        rsConfig->p[1] = 1;  //1100001
    }else if (m == 7){
        rsConfig->p[1] = 1;  //10010001
    }else if (m == 8){
        rsConfig->p[2] = rsConfig->p[3] = rsConfig->p[4] = 1; //101110001
    }else if (m == 9){
        rsConfig->p[4] = 1;  //1000100001
    } else if (m == 10){
        rsConfig->p[3] = 1;  //10010000001
    }else if (m == 11){
        rsConfig->p[2] = 1;
    }else if (m == 12){
        rsConfig->p[3] = rsConfig->p[4] = rsConfig->p[7] = 1;
    }else if (m == 13){
        rsConfig->p[1] = rsConfig->p[3] = rsConfig->p[4] = 1;
    }else if (m == 14){
        rsConfig->p[1] = rsConfig->p[11] = rsConfig->p[12] = 1;
    }else if (m == 15){
        rsConfig->p[1] = 1;
    }else if (m == 16){
        rsConfig->p[2] = rsConfig->p[3] = rsConfig->p[5] = 1;
    }else if (m == 17){
        rsConfig->p[3] = 1;
    }else if (m == 18){
        rsConfig->p[7] = 1;
    }else if (m == 19){
        rsConfig->p[1] = rsConfig->p[5] = rsConfig->p[6] = 1;
    }else if (m == 20){
        rsConfig->p[3] = 1;
    }
    
    // 原始多项式
    printf("\nPrimitive polynomial of GF(2^%d), (LSB first)   p(x) = ",m);
    
    for (i = 0; i <= m; i++){
        //n *= 2;
        printf("%1d", rsConfig->p[i]);
    }
    
    return status;
}


/*
 * 生成GF域
 * 从不可约的多项式p(X)中生成GF(2^m)域，
 */
// generate GF(2^m) from the irreducible polynomial p(X) in p[0]..p[m]
//
// lookup tables:  log->vector form           alpha_to[] contains j=alpha**i;
//                 vector form -> log form    index_of[j=alpha**i] = i
// alpha=2 is the primitive element of GF(2^m)
int generate_gf(RSConfig *rsConfig){
    if (rsConfig == NULL) {
        LOGD("\nError : RS - generate_gf - rsConfig = NULL\n");
        return -1;
    }
    
    int status = 0;
    int m = rsConfig->m;
    int n = rsConfig->n;
    register int i, mask;
    mask = 1;
    rsConfig->alpha_to[m] = 0;
    for (i=0; i<m; i++) {
        rsConfig->alpha_to[i] = mask;
        rsConfig->index_of[rsConfig->alpha_to[i]] = i;
        if (rsConfig->p[i]!=0)
            rsConfig->alpha_to[m] ^= mask;
        mask <<= 1;
    }
    rsConfig->index_of[rsConfig->alpha_to[m]] = m;
    mask >>= 1;
    for (i=m+1; i<n; i++){
        if (rsConfig->alpha_to[i-1] >= mask)
            rsConfig->alpha_to[i] = rsConfig->alpha_to[m] ^ ((rsConfig->alpha_to[i-1]^mask)<<1);
        else
            rsConfig->alpha_to[i] = rsConfig->alpha_to[i-1]<<1;
        
        rsConfig->index_of[rsConfig->alpha_to[i]] = i;
    }
    rsConfig->index_of[0] = -1;
    
    
#ifdef PRINT_GF
    printf("\nTable of GF(%d):\n",n);
    printf("----------------------\n");
    printf("   i\tvector \tlog\n");
    printf("----------------------\n");
    for (i=0; i<=n; i++){
        printf("%4d\t%4d\t%4d\n", i, rsConfig->alpha_to[i], rsConfig->index_of[i]);
    }
#endif
    
    return status;
}


// Compute the generator polynomial of the t-error correcting, length
// n=(2^m -1) Reed-Solomon code from the product of (X+alpha^i), for
// i = init_zero, init_zero + 1, ..., init_zero+length-k-1
int gen_poly(RSConfig *rsConfig)
{
    int status = 0;
    if(rsConfig == NULL){
        LOGD("\nError : RS - gen_poly() Fail! rsConfig = NULL\n");
        return -1;
    }

    int n = rsConfig->n;
    int length = rsConfig->length;
    int k = rsConfig->k;
    
    register int i,j;
    rsConfig->g[0] = rsConfig->alpha_to[init_zero];  //  <--- vector form of alpha^init_zero
    rsConfig->g[1] = 1;     // g(x) = (X+alpha^init_zero)
    for (i=2; i<=length-k; i++) {
        rsConfig->g[i] = 1;
        for (j=i-1; j>0; j--)
            if (rsConfig->g[j] != 0)
                rsConfig->g[j] = rsConfig->g[j-1]^ rsConfig->alpha_to[(rsConfig->index_of[rsConfig->g[j]]+i+init_zero-1)%n];
            else
                rsConfig->g[j] = rsConfig->g[j-1];
        rsConfig->g[0] = rsConfig->alpha_to[(rsConfig->index_of[rsConfig->g[0]]+i+init_zero-1)%n];
    }
    // convert g[] to log form for quicker encoding
    for (i=0; i<=length-k; i++)
        rsConfig->g[i] = rsConfig->index_of[rsConfig->g[i]];
    
#ifdef PRINT_POLY
    printf("Generator polynomial (independent term first):\ng(x) = ");
    for (i=0; i<=length-k; i++)
        printf("%5d", rsConfig->g[i]);
    printf("\n");
#endif
    
    return status;
    
}



// Compute the 2t parity(等价) symbols in rs_b[0]..rs_b[2*t-1]
// data[] is input and rs_b[] is output in polynomial form.
// Encoding is done by using a feedback shift register with connections
// specified by the elements of g[].
int encode_rs(RSConfig *rsConfig){
    int status = 0;
    if(rsConfig == NULL){
        LOGD("\nError : RS - encode_rs() Fail! rsConfig = NULL\n");
        return -1;
    }
    
    int n = rsConfig->n;
    int length = rsConfig->length;
    int k = rsConfig->k;
    
    register int i,j;
    int feedback;
    for (i=0; i<length-k; i++)
        rsConfig->rs_b[i] = 0;
    for (i=k-1; i>=0; i--)
    {
        
        feedback = rsConfig->index_of[rsConfig->data[i]^rsConfig->rs_b[length-k-1]];
        if (feedback != -1)
        {
            for (j=length-k-1; j>0; j--)
                if (rsConfig->g[j] != -1)
                    rsConfig->rs_b[j] = rsConfig->rs_b[j-1]^rsConfig->alpha_to[(rsConfig->g[j]+feedback)%n];
                else
                    rsConfig->rs_b[j] = rsConfig->rs_b[j-1];
            rsConfig->rs_b[0] = rsConfig->alpha_to[(rsConfig->g[0]+feedback)%n];
        }
        else
        {
            for (j=length-k-1; j>0; j--)
                rsConfig->rs_b[j] = rsConfig->rs_b[j-1];
            rsConfig->rs_b[0] = 0;
        }
    }
    
    return status;
}




///* Decode the RS Code. */
int decode_rs(RSConfig *rsConfig) {
    if(rsConfig == NULL){
        LOGD("\nError : RS - decode_rs() Fail! rsConfig = NULL\n");
        return -1;
    }
    
    int status = 4;
    
    int n = rsConfig->n;
    int length = rsConfig->length;
    int k = rsConfig->k;
    int t = rsConfig->r/2;
    int t2 = rsConfig->r;
    
//    int *recd = rsConfig->recd;
//    int *alpha_to = rsConfig->alpha_to;
//    int *index_of = rsConfig->index_of;
    
//    int numera = 0;
//    int *era = rsConfig->era;
    
	// Add By wei.liu3 2013/2/25
	// 使用动态分配内存的方法，因为VC不支持数组的非常量大小定义
    register int i,j,u,q;
    // int elp[length-k+2][length-k], d[length-k+2], l[length-k+2], u_lu[length-k+2], s[n-k+1], forney[n];// length-k+2
    int count=0, syn_error=0;//, tau[t], root[t], loc[t];
    // int err[length], reg[t+1], aux[t+1], omega[length-k+1], phi[length-k+1], phiprime[length-k+1];
    int degphi, ell;
	// 规范化定义(VS支持)
	// 1
	int** elp = (int **)malloc(sizeof(int*) * (length - k + 2));
	for (int i = 0; i < length - k + 2; ++i) {
		elp[i] = (int *)malloc(sizeof(int) * (length - k));
		memset(elp[i], 0, sizeof(int) * (length - k));
	}
	int *d = (int *)malloc(sizeof(int) * (length - k + 2));
	memset(d, 0, sizeof(int) * (length - k + 2));
	int *l = (int *)malloc(sizeof(int) * (length - k + 2));
	memset(l, 0, sizeof(int) * (length - k + 2));
	int *u_lu = (int *)malloc(sizeof(int) * (length - k + 2));
	memset(u_lu, 0, sizeof(int) * (length - k + 2));
	int *s = (int *)malloc(sizeof(int) * (n - k + 1));
	memset(s, 0, sizeof(int) * (n - k + 1));
	int *forney = (int *)malloc(sizeof(int) * n);
	memset(forney, 0, sizeof(int) * n);
	// 2
	int *tau = (int *)malloc(sizeof(int) * t);
	memset(tau, 0, sizeof(int) * t);
	int *root = (int *)malloc(sizeof(int) * t);
	memset(root, 0, sizeof(int) * t);
	int *loc = (int *)malloc(sizeof(int) * t);
	memset(loc, 0, sizeof(int) * t);
	// 3
	int *err = (int *)malloc(sizeof(int) * length);
	memset(err, 0, sizeof(int) * length);
	int *reg = (int *)malloc(sizeof(int) * (t + 1));
	memset(reg, 0, sizeof(int) * (t + 1));
	int *aux = (int *)malloc(sizeof(int) * (t + 1));
	memset(aux, 0, sizeof(int) * (t + 1));
	int *omega = (int *)malloc(sizeof(int) * (length - k + 1));
	memset(omega, 0, sizeof(int) * (length - k + 1));
	int *phi = (int *)malloc(sizeof(int) * (length - k + 1));
	memset(phi, 0, sizeof(int) * (length - k + 1));
	int *phiprime = (int *)malloc(sizeof(int) * (length - k + 1));
	memset(phiprime, 0, sizeof(int) * (length - k + 1));
	// ~Add By wei.liu3 2013/2/25
    
    // Compute the syndromes
    
#ifdef PRINT_SYNDROME
    printf("\ns = 0 ");
#endif
    for (i=1; i<=t2; i++)
    {
        s[i] = 0;
        for (j=0; j<length; j++)
            if (rsConfig->recd[j]!=0){
                s[i] ^= rsConfig->alpha_to[(rsConfig->index_of[rsConfig->recd[j]]+(i+init_zero-1)*j)%n];
            }
        
        // convert syndrome from vector form to log form  */
        if (s[i]!=0)
            syn_error=1;         // set flag if non-zero syndrome => error
        //
        // Note:    If the code is used only for ERROR DETECTION, then
        //          exit program here indicating the presence of errors.
        //
        s[i] = rsConfig->index_of[s[i]];
#ifdef PRINT_SYNDROME
        printf("%4d ", s[i]);
#endif
    }
    if (syn_error)       // if syndromes are nonzero then try to correct
    {
        s[0] = 0;
        // S(x) = 1 + s_1x + ...
        // TO HANDLE ERASURES
        
        // if erasures are present, compute the erasure locator polynomial, tau(x)
        if (rsConfig->numera) {
            for (i=0; i<=t2; i++) {
                tau[i] = 0;
                aux[i] = 0;
            }
            
            aux[1] = rsConfig->alpha_to[rsConfig->era[0]];
            aux[0] = 1;       // (X + era[0])
            if (rsConfig->numera>1){
                for (i=1; i<rsConfig->numera; i++)
                {
                    rsConfig->p[1] = rsConfig->era[i];
                    rsConfig->p[0] = 0;
                    for (j=0; j<2; j++){
                        for (ell=0; ell<=i; ell++){
                            // Line below added 8/17/2003
                            if ((rsConfig->p[j] !=-1) && (aux[ell]!=0)){
                                tau[j+ell] ^= rsConfig->alpha_to[(rsConfig->p[j]+rsConfig->index_of[aux[ell]])%n];
                            }
                        }
                    }
                    
                    if (i != (rsConfig->numera-1)){
                        for (ell=0; ell<=(i+1); ell++)
                        {
                            aux[ell] = tau[ell];
                            tau[ell] = 0;
                        }
                    }
                }
            }else {
                tau[0] = aux[0]; tau[1] = aux[1];
            }
            
            // Put in index (log) form
            for (i=0; i<=rsConfig->numera; i++){
                tau[i] = rsConfig->index_of[tau[i]]; /* tau in log form */
            }
            
#ifdef PRINT_SYNDROME
            printf("\ntau =    ");
            for (i=0; i<=rsConfig->numera; i++){
                printf("%4d ", tau[i]);
            }
            printf("\nforney = ");
#endif
            // Compute FORNEY modified syndrome:
            //            forney(x) = [ s(x) tau(x) + 1 ] mod x^{t2}
            
            for (i=0; i<=n-k; i++){
                forney[i] = 0;
            }
            
            for (i=0; i<=n-k; i++){
                for (j=0; j<=rsConfig->numera; j++){
                    if (i+j <= (n-k)){ // mod x^{n-k+1}
                        if ((s[i]!=-1)&&(tau[j]!=-1)){
                            forney[i+j] ^= rsConfig->alpha_to[(s[i]+tau[j])%n];
                        }
                    }
                }
            }
            
            
            forney[0] ^= 1;
            for (i=0; i<=n-k; i++){
                forney[i] = rsConfig->index_of[forney[i]];
            }
            
#ifdef PRINT_SYNDROME
            for (i=0; i<=n-k; i++){
                printf("%4d ", forney[i]);
            }
#endif
            
        }
        else // No erasures
        {
            tau[0]=0;
			// BEGIN, RS FIX BUG, ADDED BY MARK.
            for (i=1; i<=n-k; i++){
			// END, RS FIX BUG, ADDED BY MARK.
                forney[i] = s[i];
            }
        }
#ifdef PRINT_SYNDROME
        printf("\n");
#endif
        // --------------------------------------------------------------
        //    THE BERLEKAMP-MASSEY ALGORITHM FOR ERRORS AND ERASURES
        // --------------------------------------------------------------
        
        
        // initialize table entries
        d[0] = 0;                // log form
        d[1] = forney[rsConfig->numera+1]; // log form
        elp[0][0] = 0;           // log form
        elp[1][0] = 1;           // vector form
        for (i=1; i<t2; i++)
        {
            elp[0][i] = -1;   // log form
            elp[1][i] = 0;    // vector form
        }
        l[0] = 0;
        l[1] = 0;
        u_lu[0] = -1;
        u_lu[1] = 0;
        u = 0;
        if (rsConfig->numera<t2) {  // If errors can be corrected
            do
            {
                u++;
                if (d[u]==-1)
                {
#ifdef PRINT_SYNDROME
                    printf("d[%d] is zero\n",u);
#endif
                    l[u+1] = l[u];
                    for (i=0; i<=l[u]; i++)
                    {
                        elp[u+1][i] = elp[u][i];
                        elp[u][i] = rsConfig->index_of[elp[u][i]];
                    }
                }
                else           // search for words with greatest u_lu[q] for which d[q]!=0
                {
                    q = u-1;
                    while ((d[q]==-1) && (q>0))
                        q--;
                    
                    // have found first non-zero d[q]
                    if (q>0)
                    {
                        j=q;
                        do
                        {
                            j--;
                            if ((d[j]!=-1) && (u_lu[q]<u_lu[j]))
                                q = j;
                        } while (j>0);
                    }
                    
#ifdef PRINT_SYNDROME
                    printf("u = %4d, q = %4d, d[q] = %4d d[u] = %4d\n", u, q, d[q],d[u]);
#endif
                    
                    // have now found q such that d[u]!=0 and u_lu[q] is maximum
                    // store degree of new elp polynomial
                    if (l[u]>l[q]+u-q)
                        l[u+1] = l[u];
                    else
                        l[u+1] = l[q]+u-q;
                    
#ifdef PRINT_SYNDROME
                    printf("l[q] = %4d, l[u] = %4d\n", l[q], l[u]);
#endif
                    // compute new elp(x)
                    // for (i=0; i<t2-numera; i++)    elp[u+1][i] = 0;
                    
                    for (i=0; i<t2; i++)    elp[u+1][i] = 0;
                    for (i=0; i<=l[q]; i++)
                        if (elp[q][i]!=-1)
                            elp[u+1][i+u-q] = rsConfig->alpha_to[(d[u]+n-d[q]+elp[q][i])%n];
                    for (i=0; i<=l[u]; i++)
                    {
                        elp[u+1][i] ^= elp[u][i];
                        elp[u][i] = rsConfig->index_of[elp[u][i]];
                    }
                    
#ifdef PRINT_SYNDROME
                    printf("l[u+1] = %4d, elp[u+1] = ", l[u+1]);
                    for (i=0;  i<=l[u+1]; i++)
                        printf("%4d ",rsConfig->index_of[elp[u+1][i]]); printf("\n");
#endif
                    
                }
                u_lu[u+1] = u-l[u+1];
                // compute (u+1)th discrepancy
                // if (u<t2)    // no discrepancy computed on last iteration
                if (u<(t2-rsConfig->numera)) // no discrepancy computed on last iteration
                    // --- if ( u < (l[u+1]+t-1-(numera/2)) )
                {
                    // if (s[u+1]!=-1)
                    if (forney[rsConfig->numera+u+1]!=-1)
                        d[u+1] = rsConfig->alpha_to[forney[rsConfig->numera+u+1]];
                    else
                        d[u+1] = 0;
                    
#ifdef PRINT_SYNDROME
                    printf("discrepancy for u = %d: d[u+1] = %4d\n", u, rsConfig->index_of[d[u+1]]);
#endif
                    
                    for (i=1; i<=l[u+1]; i++)
                        // if ((s[u+1-i]!=-1) && (elp[u+1][i]!=0))
                        //   d[u+1] ^= alpha_to[(s[numera+u+1-i]
                        if ((forney[rsConfig->numera+u+1-i]!=-1) && (elp[u+1][i]!=0)) {
                            d[u+1] ^= rsConfig->alpha_to[(forney[rsConfig->numera+u+1-i] + rsConfig->index_of[elp[u+1][i]])%n];
                            
#ifdef PRINT_SYNDROME
                            printf("i=%d, forney[%d] = %4d, d[u+1] = %4d\n",i,(rsConfig->numera)+u+1-i,                    forney[(rsConfig->numera)+u+1-i],rsConfig->index_of[d[u+1]]);
#endif
                        }
                    
                    d[u+1] = rsConfig->index_of[d[u+1]];     // put d[u+1] into index form
                    
#ifdef PRINT_SYNDROME
                    printf("d[u+1] = %4d\n", d[u+1]);
#endif
                    
                }
                // } while ((u<t2) && (l[u+1]<=t));
            } while ((u<(t2-rsConfig->numera)) && (l[u+1]<=((t2-rsConfig->numera)/2)));
        }
        
                
        // else
        // case of 2t erasures
        // {
        //   elp[1][0] = 0;
        //   count = 0;
        // }
        
        u++;
        if (l[u]<=t-rsConfig->numera/2)         // can correct errors
        {
            // put elp into index form
            for (i=0; i<=l[u]; i++)
                elp[u][i] = rsConfig->index_of[elp[u][i]];
            printf("\nBM algorithm, after %d iterations:\nsigma = ", (u-1));
            for (i=0; i<=l[u]; i++)
                printf("%4d ", elp[u][i]); printf("\n");
            // find roots of the error location polynomial
            for (i=1; i<=l[u]; i++)
                reg[i] = elp[u][i];
            count = 0;
            for (i=1; i<=n && count<t; i++)
            {
                q = 1;
                for (j=1; j<=l[u]; j++)
                    if (reg[j]!=-1) {
                        reg[j] = (reg[j]+j)%n;
                        q ^= rsConfig->alpha_to[reg[j]];
                    }
                if (!q)        // store root and error location number indices
                {
//                    LOGD("\n求错误位置：q = %d, n-i = %d, length = %d\n", q, n-i, length);
                    if (n-i < length) {
                        root[count] = i;
                        loc[count] = n-i;
//                        LOGD("\n错误位置：loc[%d] = %d, q = %d, \n", count, loc[count], q);
                        
#ifdef PRINT_SYNDROME
                    printf("loc[%4d] = %4d\n", count, loc[count]);
#endif
                        
                        count++;
                    }
                    else{
                        LOGE("\n非法的错误位置！错误位置为：n-i = %d\n", n-i);
                    }
                }
            }
            if (count==l[u])    // no. roots = degree of elp hence <= t errors
            {
                
                // Compute the errata evaluator polynomial, omega(x)
                forney[0] = 0;  // as a log, to construct 1+T(x)
                
                for (i=0; i<=t2; i++)
                    omega[i] = 0;
                for (i=0; i<=t2; i++)
                {
                    for (j=0; j<=l[u]; j++)
                    {
                        if (i+j <= t2) // mod x^{t2}
                            if ((forney[i]!=-1) && (elp[u][j]!=-1))
                                omega[i+j] ^= rsConfig->alpha_to[(forney[i]+elp[u][j])%n];
                    }
                }
                for (i=0; i<=t2; i++)
                    omega[i] = rsConfig->index_of[omega[i]];
                
#ifdef PRINT_SYNDROME
                printf("\nomega =    ");
                for (i=0; i<=t2; i++)
                    printf("%4d ", omega[i]);
                printf("\n");
#endif
                
                
                
                // Compute the errata locator polynomial, phi(x)
                degphi = rsConfig->numera+l[u];
                
                for (i=0; i<=degphi; i++)
                    phi[i] = 0;
                
                for (i=0; i<=rsConfig->numera; i++)
                    for (j=0; j<=l[u]; j++)
                        if ((tau[i]!=-1)&&(elp[u][j]!=-1))
                            phi[i+j] ^= rsConfig->alpha_to[(tau[i]+elp[u][j])%n];
                
                for (i=0; i<=degphi; i++)
                    phi[i] = rsConfig->index_of[phi[i]];
                
#ifdef PRINT_SYNDROME
                printf("phi =      ");
                for (i=0; i<=degphi; i++)
                    printf("%4d ", phi[i]);
                printf("\n");
#endif
                
                // Compute the "derivative" of phi(x): phiprime
                for (i=0; i<=degphi; i++)
                    phiprime[i] = -1; // as a log
                
                for (i=0; i<=degphi; i++)
                    if (i%2)  // Odd powers of phi(x) give terms in phiprime(x)
                        phiprime[i-1] = phi[i];
                
#ifdef PRINT_SYNDROME
                printf("phiprime = ");
                for (i=0; i<=degphi; i++)
                    printf("%4d ", phiprime[i]);
                printf("\n\n");
#endif
                
                LOGD("**** count = %d\n", count);
                if (rsConfig->numera)             // Add erasure positions to error locations
                    for (i=0; i<rsConfig->numera; i++) {
                        loc[count+i] = rsConfig->era[i];
                        root[count+i] = (n-rsConfig->era[i])%n;
                    }
                
                
                // evaluate errors at locations given by errata locations, loc[i]
                //                             for (i=0; i<l[u]; i++)
                for (i=0; i<degphi; i++)
                {
                    // compute numerator of error term
                    err[loc[i]] = 0;
                    for (j=0; j<=t2; j++)
                        if ((omega[j]!=-1)&&(root[i]!=-1))
                            err[loc[i]] ^= rsConfig->alpha_to[(omega[j]+j*root[i])%n];
                    // -------  The term loc[i]^{2-init_zero}
                    if ((err[loc[i]]!=0)&&(loc[i]!=-1))
                        err[loc[i]] = rsConfig->alpha_to[(rsConfig->index_of[err[loc[i]]]                                               +loc[i]*(2-init_zero+n))%n];
                    if (err[loc[i]]!=0)
                    {
                        err[loc[i]] = rsConfig->index_of[err[loc[i]]];
                        // compute denominator of error term
                        q = 0;
                        for (j=0; j<=degphi; j++)
                            if ((phiprime[j]!=-1)&&(root[i]!=-1))
                                q ^= rsConfig->alpha_to[(phiprime[j]+j*root[i])%n];
                        
                        // Division by q
                        err[loc[i]] = rsConfig->alpha_to[(err[loc[i]]-rsConfig->index_of[q]+n)%n];
                        
#ifdef PRINT_SYNDROME
                        printf("errata[%4d] = %4d (%4d) \n",loc[i],rsConfig->index_of[err[loc[i]]],err[loc[i]]);
#endif
                        
                        rsConfig->recd[loc[i]] ^= err[loc[i]];
                    }
                }
                
                printf("\nCor =");
                for (i=0; i<length; i++) {
                    printf("%4d ", rsConfig->index_of[rsConfig->recd[i]]);
                }
                printf("\n     ");
                for (i=0; i<length; i++) {
                    printf("%4d ", rsConfig->recd[i]);
                }
                printf("\n");
            }
            else{    // no. roots != degree of elp => >t errors and cannot solve
               status = -1;
            }
        }
        else{         // elp has degree has degree >t hence cannot solve
           status = -2;
        }
    }
    else       // no non-zero syndromes => no errors: output received codeword
        ;

    // BEGIN, solve the memory leak, ADDED BY MARK.
    if(elp != NULL)
    {
        for (int i = 0; i < length - k + 2; ++i) {
            if(elp[i] != NULL) {
                free(elp[i]);
                elp[i] = NULL;
            }
        }
        free(elp);
        elp = NULL;
    }

    if (d != NULL)
    {
        free(d);
        d = NULL;
    }

    if (l != NULL)
    {
        free(l);
        l = NULL;
    }

    if (u_lu != NULL)
    {
        free(u_lu);
        u_lu = NULL;
    }

    if (s != NULL)
    {
        free(s);
        s = NULL;
    }

    if (forney != NULL)
    {
        free(forney);
        forney = NULL;
    }

    if (tau != NULL)
    {
        free(tau);
        tau = NULL;
    }

    if (root != NULL)
    {
        free(root);
        root = NULL;
    }

    if (loc != NULL)
    {
        free(loc);
        loc = NULL;
    }

    // 3.
    if (reg != NULL)
    {
        free(reg);
        reg = NULL;
    }

    if (err != NULL)
    {
        free(err);
        err = NULL;
    }

    if (aux != NULL)
    {
        free(aux);
        aux = NULL;
    }

    if (omega != NULL)
    {
        free(omega);
        omega = NULL;
    }

    if (phi != NULL)
    {
        free(phi);
        phi = NULL;
    }

    if (phiprime != NULL)
    {
        free(phiprime);
        phiprime = NULL;
    }

    // END, solve the memory leak, ADDED BY MARK.
    
    return status;
}
