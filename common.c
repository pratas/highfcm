#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "defs.h"
#include "mem.h"
#include "common.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

float *logTable;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

uint32_t fLog2(uint64_t i)
  {
  uint32_t n, m, k = 32, o = (i & (i - 1)) ? 1 : 0;
  static const uint64_t sizes[6] = 
  { 0x0000000000000002ull, 0x000000000000000Cull, 0x00000000000000F0ull, 
    0x000000000000FF00ull, 0x00000000FFFF0000ull, 0xFFFFFFFF00000000ull };

  for(n = 6 ; n-- ; )
    {
    o += (m = (i & *(sizes+n)) ? k : 0);
    i >>= m;
    k >>= 1;
    }

  return o;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

uint64_t NBytesInFile(const char *fn)
  {
  uint64_t size = 0;
  FILE *file = Fopen(fn, "r");

  fseek(file, 0, SEEK_END);
  if((size = ftell(file)) < 1000)
    {
    fprintf(stderr, "Error: input file is very small!\n");
    exit(1);
    }
    
  fclose(file);

  return size;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void FillLogTable(uint32_t nSym, uint32_t alphaDen, uint32_t maxCHigh)
  {
  uint32_t n, maxSize = nSym * maxCHigh * alphaDen;
 
  logTable = (float *) Malloc(maxSize * sizeof(float));
  for(n = 1 ; n != maxSize ; ++n)
    logTable[n] = fLog2(n);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

double SearchLog(uint32_t idx)
  {
  return logTable[idx];
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

uint8_t GetComp(uint8_t symbol)
  {
  switch(symbol)
    {
    case 'A': return 'T';
    case 'C': return 'G';
    case 'G': return 'C';
    case 'T': return 'A';
    default:  return symbol;
    }
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

uint8_t GetComp2(uint8_t symbol)
  {
  switch(symbol)
    {
    case 0: return 3;
    case 1: return 2;
    case 2: return 1;
    case 3: return 0;
    default:  return symbol;
    }
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

FILE *Fopen(const char *path, const char *mode)
  {
  FILE *file = fopen(path, mode);

  if(file == NULL)
    {
    fprintf(stderr, "Error opening: %s (mode %s).\n", path, mode);
    exit(1);
    }

  return file;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

uint8_t *reverseStr(uint8_t *str, uint32_t end)
  {
  uint32_t start = 0;

  while(start < end)
    {
    str[start] ^= str[end];
    str[end]   ^= str[start];
    str[start] ^= str[end];
    ++start;
    --end;
    }

  return str;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
