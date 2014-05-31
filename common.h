#ifndef COMMON_H_INCLUDED
#define COMMON_H_INCLUDED

#include "context.h"
#include "defs.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

FILE      *Fopen       (const char *, const char *);
void      FillLogTable (uint32_t, uint32_t, uint32_t);
double    SearchLog    (uint32_t);
uint32_t  fLog2        (uint64_t);
uint8_t   GetComp      (uint8_t);
uint8_t   GetComp2     (uint8_t);
uint64_t  NBytesInFile (const char *);
uint8_t   *reverseStr  (uint8_t *, uint32_t);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif
