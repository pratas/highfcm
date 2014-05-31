#ifndef AC_HEADER
#define AC_HEADER

#include <stdio.h>
#include "defs.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

typedef struct
  {
  FILE      *fp;
  uint64_t  low;
  uint64_t  high;
  uint64_t  fbits;
  uint64_t  total_bits;
  uint32_t  buffer; 
  uint32_t  bits_to_go;
  uint32_t  index;
  uint8_t   cBuf[CBUF_SIZE];
  }
ac_encoder;

typedef struct
  {
  FILE      *fp;
  uint64_t  value;
  uint64_t  low;
  uint64_t  high;
  uint32_t  buffer;
  uint32_t  indexBuf;
  uint32_t  maxBuf;
  uint8_t   buf[CBUF_SIZE];                               
  uint32_t  bits_to_go;
  }
ac_decoder;

typedef struct
  {
  uint32_t  cfreq[256];
  uint32_t  nsym;
  }
ac_model;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

uint64_t readNBits          (uint32_t, ac_decoder * , ac_model *);
void     writeNBits         (uint64_t, uint32_t , ac_encoder *, ac_model *);
void     ac_encoder_init    (ac_encoder *  , const char * );
void     ac_decoder_init    (ac_decoder *  , const char * );
void     ac_model_init      (ac_model   *  , uint32_t   );
void     shotgunEncode      (ac_encoder *  , uint32_t [], uint32_t );
void     ac_encode_symbol   (ac_encoder *  , ac_model * , uint32_t );
void     acEncodeBinary     (ac_encoder *  , ac_model * , uint32_t );
void     acEncodeBin0       (ac_encoder *  , ac_model * );
void     acEncodeBin1       (ac_encoder *  , ac_model * );
void     acEncode0          (ac_encoder *  , ac_model * );
void     acEncode1          (ac_encoder *  , ac_model * );
void     acEncode2          (ac_encoder *  , ac_model * );
void     acEncode3          (ac_encoder *  , ac_model * );
uint32_t acDecSymHighSizeVar(ac_decoder *  , ac_model * );
uint32_t acDecSymLowSizeVar (ac_decoder *  , ac_model * );
uint32_t acDecodeBinary     (ac_decoder *  , ac_model * );
uint32_t acDecode4Symbols   (ac_decoder *  , ac_model * );
void     ac_encoder_done    (ac_encoder *  );
void     ac_decoder_done    (ac_decoder *  );
uint64_t ac_encoder_bits    (ac_encoder *  );
void     ac_model_done      (ac_model   *  );

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif
