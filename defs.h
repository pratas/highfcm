#ifndef DEFS_H_INCLUDED
#define DEFS_H_INCLUDED

#include <stdint.h>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#define DEBUG
#define DEBUG_1

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

typedef uint16_t ACCounter;              //Size of context counters for arrays
typedef uint16_t HCCounter;         //Size of context counters for hash tables 
typedef uint64_t ULL;             
typedef uint8_t  UChar;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

uint8_t       verbose;                             // WARNING: global variable

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

typedef struct
  {
  uint64_t    iPos;                                           // Init Position
  uint64_t    ePos;                                            // End Position
  }
Positions;

typedef struct
  {
  uint64_t    size;
  uint64_t    bits;
  uint64_t    mem;
  }
Stats;

typedef struct
  {
  uint8_t     lowOrder;
  uint8_t     highOrder;
  uint8_t     usedOrder;
  uint32_t    maxCLow;
  uint32_t    maxCUsed;
  uint32_t    alphaDen;
  uint32_t    alphaDEv;
  uint16_t    blockSize;
  uint8_t     mode;
  uint8_t     ir;
  }
Level;

typedef struct
  {
  char        fn[512];
  uint64_t    size;
  uint8_t     nThreads;
  uint8_t     id;
  Level       L;
  Stats       S;
  Positions   P;
  }
Data;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#define DEFAULT_BS               100
#define MAX_ALPHA_DEN            50
#define D_CTX_LOW_ORDER          4
#define HASH_SIZE                19999999    //19999999
#define MAX_OUT_LS               4096
#define MAX_SFIELD               512
#define SMODEL_GUARD             16
#define L_GUARD_BUF              32
#define BUF_SIZE_CAT             1000000            // CAT DECOMPRESSED FILES
#define BUF_SIZE                 65536                     
#define CBUF_SIZE                65536                           //used in ac
#define BREAK_CHAR               10
#define MAX_FILESIZE             512
#define MAXSIZEA                 255
#define INVALID_S                255
#define NSIGNF                   256
#define MAX_ARG_SIZE             1024
         //ALERT: ASSERT THAT 2^Code_value_bits > nSymbols * MAXCNT * alphaDen
#define MAXCNT			 (((uint64_t) 1 << (sizeof(ACCounter)*8)) - 2)
#define MAX_ASIZE                100
#define SCHAR			 sizeof(uint8_t)
#define MAX_SIZE_BITS		 UINT32_MAX
#define NBITS_UINT32             sizeof(uint32_t) * 8
#define NBITS_UINT64             sizeof(uint64_t) * 8
#define MAX(a,b)                 ((a) > (b) ? a : b)
#define MIN(a,b)                 ((a) < (b) ? a : b)

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif

