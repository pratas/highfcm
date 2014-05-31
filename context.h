// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
//    Copyright 2014 University of Aveiro, Portugal, All Rights Reserved.    -
//    This file is part of HighFCM, contacts: pratas@ua.pt | ap@ua.pt        -
//    Description: Data structures for handling finite-context models        -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#ifndef CONTEXT_H_INCLUDED
#define CONTEXT_H_INCLUDED

#include "defs.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

typedef struct
  {
  uint64_t   key;                               //The key stored in this entry 
  HCCounter  *counters;                                 //The context counters 
  }
Entry;

typedef struct
  {
  uint16_t   *entrySize;                        //Number of keys in this entry
  Entry      **entries;                    //The heads of the hash table lists
  HCCounter  **zeroCounters;
  uint32_t   nUsedEntries;
  uint32_t   nUsedKeys;
  }
HashTable;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

typedef struct
  {
  ACCounter  *counters;	                                       // FCM Counters
  uint64_t   *multipliers;	                        // Calculated products
  uint64_t   nPModels;	                //Maximum number of probability models
  uint32_t   ctxSize;	                  // Current depth of context template
  uint32_t   nSymbols;	                           // Number of coding symbols
  uint64_t   pModelIdx;                        // Index of probabilistic model
  uint64_t   idxIr;           // Inverted repeats index of probabilistic model
  uint32_t   deltaDen;       // Delta denominator (for probability estimation)
  uint32_t   maxC;                     // Maxmimum probabilistic model counter
  HashTable  hTable;
  }
CModel;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void         FreeCModelCont    (CModel *);
void         FreeCModelHash    (CModel *);
double       InfoSym           (uint32_t, uint32_t);
CModel       *CreateCModel     (uint32_t, uint32_t, uint32_t, uint8_t, 
                               uint32_t);
void         UpdateHashCounter (CModel *, uint64_t, uint8_t);
HCCounter    *GetHCCounters    (HashTable *, uint64_t);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif
