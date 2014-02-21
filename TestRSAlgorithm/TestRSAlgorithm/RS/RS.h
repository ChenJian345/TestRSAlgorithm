//
//  RS.h
//  TestCProjectOnXCode
//
//  Created by Mark on 12-11-19.
//  Copyright (c) 2012年 Mark Chen. All rights reserved.
//

#include "Log.h"


#ifndef TestCProjectOnXCode_RS_h
#define TestCProjectOnXCode_RS_h

//#define PRINT_GF
//#define PRINT_POLY
//#define PRINT_SYNDROME

#define init_zero   1

// 对控制帧进行RS编码所使用的RS规范
// RS(6, 10, 6, 4)
#define RS_CONTROL_M        6
#define RS_CONTROL_LENGTH   10
#define RS_CONTROL_K        6
#define RS_CONTROL_R        4

// 对数据帧进行RS编码所用的RS规范
// RS(8, 10, 6, 4)
#define RS_DATA_M        8
#define RS_DATA_LENGTH   10
#define RS_DATA_K        6
#define RS_DATA_R        4

typedef unsigned char byte;


//static int p[13];     // array with 10 elements.

/**
 *  Define the RS Configuration structure.
 */
typedef struct RSConfig{
    int n;      // The max range of RS Code.
    int length;      // The length of code after encode.
    int k;      // The message length need to encode.
    int r;      // The count of error can checked.
    int m;      // The GF element count is 2^m
    
    int *alpha_to;
    int *index_of;
    int *g;
    int *recd;
    
    int *rs_b;
    
    int numera;
    int *era;
    
    int *data;  // 可有可无，后期传入byte流后删除
    
    int p[13];
} RSConfig;


int read_p();
int generate_gf(RSConfig *rsConfig);
int gen_poly(RSConfig *rsConfig);

/*
 * Compute the 2t parity(等价) symbols in rs_b[0]..rs_b[2*t-1]
 * data[] is input and rs_b[] is output in polynomial form.
 * Encoding is done by using a feedback shift register with connections
 * specified by the elements of g[].
 */
int encode_rs(RSConfig *rsConfig);

/* Decode the RS Code. */
int decode_rs(RSConfig *rsConfig);



/**
 *  Init RS Config.
 */
int initRS(RSConfig *rsConfig, int aM, int aLength, int aK, int aR);


/**
 *  RS编码，通过RSConfig编码并添加响应长度纠错码，纠错码长为r.
 */
int encodeWithRS(struct RSConfig *rsConfig);

/*
 *  Prepare for the decode.
 */
int prepareDecode(RSConfig *rsConfig);

/**
 *  RS解码，通过RSConfig解码，并进行纠错
 */
int decodeWithRS(struct RSConfig *rsConfig);


/**
 *  Clean up the memory allocated.
 */
int cleanup(RSConfig *rsConfig);



#endif
