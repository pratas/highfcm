// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// -                                                                         -
// -           U n i v e r s i t y   o f   A v e i r o   2 0 1 4             -
// -                                                                         -
// -       HighFCM is under GPL,  http://www.gnu.org/licenses/gpl.txt        -
// -                                                                         -
// -       HighFCM: a fast lossless asymmetric genomic compressor            - 
// -       that explores high Markov context orders                          -
// -                                                                         -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <limits.h>
#include <math.h>
#include <float.h>
#include <ctype.h>
#include <time.h>
#include <unistd.h>
#include "ac.h"
#include "mem.h"
#include "defs.h"
#include "common.h"
#include "context.h"
#ifdef THREADS
#include <pthread.h>
#endif

//#define BLOCKPRINT 1
//#define BLOCKTOTAL 1

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - C O M P R E S S O R - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void *Compress(void *data)
  {
  Data        *d = (Data *) data;
  FILE        *IN = NULL, *INR = NULL;
  #ifdef BLOCKTOTAL
  FILE        *BLKT = NULL;
  #endif
  #ifdef BLOCKPRINT
  FILE        *BLKLR = NULL, *BLKRL = NULL;
  #endif
  char        name[MAX_FILESIZE], nameRev[MAX_FILESIZE];
  ac_encoder  aced;
  ac_model    acmHeader, acmdh, acmd;
  ACCounter   *aCountersLow, *aCountersHigh;
  HCCounter   *hCountersHigh;
  CModel      *cModelHigh, *cModelLow;
  uint64_t    size, ePosP1, idxSymbol, multiplierHigh, multiplierLow;
  uint32_t    ctxHigh, ctxLow, dP, n, idxPos, 
              nSymL1, *tmpFreqsHigh, *tmpFreqsLow, alphaDen, alphaDEv,
              blockSize, k, idxBlockSymbol, nBlock, ctxUsed, maxCLow, 
              maxCHigh, maxCUsed, nLowComplex, totalBlocks;
  int32_t     idx, percLowBlock;
  uint8_t     *mask, id, norAlp[NSIGNF], alp[NSIGNF], nSym, *buf, *dBuf, sym, 
              symC, *tmp, *block, *seq, mode;
  double      cumHigh, cumLow;

  nLowComplex = 0;
  dP          = L_GUARD_BUF;
  buf         = (uint8_t *) Calloc(BUF_SIZE,      sizeof(uint8_t)); 
  dBuf        = (uint8_t *) Calloc(BUF_SIZE + dP, sizeof(uint8_t));
  mask        = (uint8_t *) Calloc(NSIGNF,        sizeof(uint8_t));
  id          = d->id;
  ctxHigh     = d->L.highOrder;
  ctxLow      = d->L.lowOrder;
  ctxUsed     = d->L.usedOrder; 
  blockSize   = d->L.blockSize;
  maxCLow     = d->L.maxCLow;
  maxCHigh    = MAXCNT;
  maxCUsed    = d->L.maxCUsed; 
  nSym        = 0;
  idx         = 0;
  size        = 0;
  #ifdef BLOCKTOTAL
  char tmpT[1024];
  sprintf(tmpT, "%s.%d.Total.profile", d->fn, id);
  BLKT       = Fopen(tmpT, "w");
  #endif
  #ifdef BLOCKPRINT
  char tmpLR[1024], tmpRL[1024];
  sprintf(tmpLR, "%s.%d.LR.profile", d->fn, id);
  sprintf(tmpRL, "%s.%d.RL.profile", d->fn, id);
  BLKLR       = Fopen(tmpLR, "w");
  BLKRL       = Fopen(tmpRL, "w");
  #endif
  IN          = Fopen(d->fn, "r");
  sprintf(nameRev, "%s%s%d", d->fn, ".rev", id);
  INR         = Fopen(nameRev, "w");
  seq         = (uint8_t  *) Malloc((d->P.ePos - d->P.iPos)+2);
  idxSymbol   = d->P.iPos;
  ePosP1      = d->P.ePos + 2;
  fseek(IN, idxSymbol, SEEK_SET);
  while((k = fread(buf, 1, BUF_SIZE, IN)))
    for(idxPos = 0 ; idxPos != k ; ++idxPos)
      if(++idxSymbol < ePosP1)
        {
        mask[*(buf+idxPos)] = 1;
        seq[size++] = buf[idxPos];
        }
  idxSymbol = size;
  while(idxSymbol)
    fprintf(INR, "%c", seq[--idxSymbol]);
  Free(seq);
  fclose(INR);

  for(n = MAXSIZEA ; --n ; )                                 // Build alphabet
    if(mask[n] == 0x1)
      {
      norAlp[nSym] = n;
      alp[n] = nSym++;                       // Don't change to pre-increment!
      }
    else
      alp[n] = INVALID_S; 

  block     = (uint8_t *) Calloc(size/blockSize+blockSize, sizeof(uint8_t));
  dBuf     += dP;
  nSymL1    = nSym - 1;
  alphaDen  = d->L.alphaDen;
  alphaDEv  = d->L.alphaDEv;

  if(alphaDEv == 0)
    {
    alphaDEv = 10;
    FillLogTable(nSym, MAX_ALPHA_DEN, 65536);
    }
  else
    FillLogTable(nSym, alphaDEv, 65536);

//  // Use hash table in small sequences
//  if(size < 15000000) // 10M
//    mode = (d->L.mode == 2 ? 1 : d->L.mode);
//  else
//    mode = (d->L.mode == 2 ? 0 : d->L.mode);

  mode = (d->L.mode == 2 ? 0 : d->L.mode);

  ////////////////////////////////////////////////////////////////////////////
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // E V A L U A T I O N
  //
  // GET A BINARY STREAM, BLOCK BY BLOCK, WITH THE LOW COMPLEXITY REGIONS
  // HAVING 0 AND HIGH COMPLEXITY REGIONS HAVING 1.
  //

  ////////////////////////////////////////////////////////////////////////////
  // 1. LEFT TO RIGHT
  //
  if(verbose == 1)
    printf("Running [L]eft to [R]ight on thread %u\n", id + 1);
  cModelHigh     = CreateCModel(ctxHigh, nSym, alphaDEv, mode, maxCHigh);
  cModelLow      = CreateCModel(ctxLow,  nSym,        1,    0, maxCLow);
  tmpFreqsHigh   = (uint32_t *) Calloc(nSym + 1, sizeof(uint32_t));
  tmpFreqsLow    = (uint32_t *) Calloc(nSym + 1, sizeof(uint32_t));
  multiplierHigh = cModelHigh->multipliers[cModelHigh->ctxSize - 1];
  multiplierLow  = cModelLow->multipliers[cModelLow->ctxSize   - 1];
  nBlock         = 0;
  idxBlockSymbol = 0;
  idxSymbol      = d->P.iPos;
  cumHigh        = 0;
  cumLow         = 0;
  fseek(IN, idxSymbol, SEEK_SET);
  while((k = fread(buf, 1, BUF_SIZE, IN)))
    for(idxPos = 0 ; idxPos != k ; ++idxPos)
      if(++idxSymbol < ePosP1)
        {
        sym = *(buf+idxPos);
        dBuf[idx] = sym = alp[sym];                          // Get DNA symbol
        tmp = &dBuf[idx-1];

        // HIGH MODEL  - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        cModelHigh->pModelIdx = ((cModelHigh->pModelIdx - *(tmp-ctxHigh) * 
        multiplierHigh) * nSym) + *tmp;

        if(mode == 0)
          {
          aCountersHigh = &cModelHigh->counters[nSym * cModelHigh->pModelIdx];
          tmpFreqsHigh[nSymL1] = 1 + alphaDEv * aCountersHigh[nSymL1];
          for(n = nSymL1 ; n-- ; )
            tmpFreqsHigh[n] = tmpFreqsHigh[n+1]+1+alphaDEv * aCountersHigh[n];
          cumHigh += (SearchLog(tmpFreqsHigh[0]) - SearchLog(tmpFreqsHigh[sym] 
          - tmpFreqsHigh[sym+1]));
          if(++aCountersHigh[sym] == maxCHigh)  // Check counter overflow&upd.
            for(n = nSym ; n-- ; )                  // For all stored counters
              aCountersHigh[n] >>= 1;                  // Divide counters by 2
          if(d->L.ir == 1)                                 // Inverted repeats
            {
            aCountersHigh = &cModelHigh->counters[nSym * (cModelHigh->idxIr = 
            cModelHigh->idxIr/nSym + GetComp2(*(tmp+1)) * multiplierHigh)];
            symC = GetComp2(*(tmp+1-ctxHigh));
            if(++aCountersHigh[symC] == maxCHigh) //Check counter overflow&upd
              for(n = nSym ; n-- ; )                // For all stored counters
                aCountersHigh[n] >>= 1;                // Divide counters by 2
            }
          }
        else
          {
          if(!(hCountersHigh = GetHCCounters(&cModelHigh->hTable, 
          cModelHigh->pModelIdx)))
            hCountersHigh = (HCCounter *) cModelHigh->hTable.zeroCounters;

          tmpFreqsHigh[nSymL1] = 1 + alphaDEv * hCountersHigh[nSymL1];
          for(n = nSymL1 ; n-- ; )
            tmpFreqsHigh[n] = tmpFreqsHigh[n+1]+1+alphaDEv * hCountersHigh[n];
          cumHigh += (SearchLog(tmpFreqsHigh[0]) - SearchLog(tmpFreqsHigh[sym] 
          - tmpFreqsHigh[sym+1]));
          UpdateHashCounter(cModelHigh, cModelHigh->pModelIdx, sym); // Update 
          if(d->L.ir == 1)                                 // Inverted repeats
            {
            cModelHigh->idxIr = cModelHigh->idxIr/nSym + GetComp2(*(tmp+1)) 
            * multiplierHigh;
            symC = GetComp2(*(tmp+1-ctxHigh));
            UpdateHashCounter(cModelHigh, cModelHigh->idxIr, symC);  // Update
            }
          }
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

        // LOWER MODEL - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        aCountersLow = &cModelLow->counters[nSym * (cModelLow->pModelIdx =
        ((cModelLow->pModelIdx - *(tmp-ctxLow) * multiplierLow) * nSym) + 
        *tmp)];
        tmpFreqsLow[nSymL1] = aCountersLow[nSymL1] + 1;
        for(n = nSymL1 ; n-- ; )
          tmpFreqsLow[n] = tmpFreqsLow[n+1] + aCountersLow[n] + 1;
        cumLow += (SearchLog(tmpFreqsLow[0]) - SearchLog(tmpFreqsLow[sym] -
        tmpFreqsLow[sym+1]));
        if(++aCountersLow[sym] == maxCLow)  // Check counter overflow & update
          for(n = nSym ; n-- ; )                    // For all stored counters
            aCountersLow[n] >>= 1;                     // Divide counters by 2
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        if(++idx == BUF_SIZE)
          {
          memcpy(dBuf-dP, dBuf+idx-dP, dP);
          idx = 0;
          }

        if(++idxBlockSymbol == blockSize)
          {
          if(cumHigh < cumLow)
            {
            #ifdef BLOCKPRINT
            fprintf(BLKLR, "%d\n", 1);
            #endif
            block[nBlock] = 1;
              ++nLowComplex;
            }
          #ifdef BLOCKPRINT
          else
            fprintf(BLKLR, "%d\n", 0);
          #endif

          ++nBlock;
          cumHigh         = 0;
          cumLow          = 0;
          idxBlockSymbol  = 0;
          }            
        }

  // Encode de rest if: sequence size % block size != 0
  if(idxBlockSymbol != 0 && cumHigh < cumLow)
    block[nBlock] = 1;
  totalBlocks = nBlock;

  // Flush Models & calc memory
  FreeCModelCont(cModelLow);
  mode == 1 ? FreeCModelHash(cModelHigh) : FreeCModelCont(cModelHigh);

  ////////////////////////////////////////////////////////////////////////////
  // 2. RIGHT TO LEFT
  //
  if(verbose == 1) 
    printf("Running [R]ight to [L]eft on thread %u\n", id + 1);
  //CREATE MODELS & RENEW/FLUSH VARIABLES 
  buf            = (uint8_t *) Calloc(BUF_SIZE,      sizeof(uint8_t));
  dBuf           = (uint8_t *) Calloc(BUF_SIZE + dP, sizeof(uint8_t));
  dBuf          += dP;
  cModelHigh     = CreateCModel(ctxHigh, nSym, alphaDEv, mode, maxCHigh);
  cModelLow      = CreateCModel(ctxLow,  nSym,        1,    0, maxCLow);
  INR            = Fopen(nameRev, "r");
  idxBlockSymbol = 0;
  cumHigh        = 0;
  cumLow         = 0;
  idx            = 0;
  fseek(INR, size % blockSize, SEEK_SET);
  --nBlock;
  while((k = fread(buf, 1, BUF_SIZE, INR)))
    for(idxPos = 0 ; idxPos < k ; ++idxPos)
      {
      sym = *(buf+idxPos);
      dBuf[idx] = sym = alp[sym];                            // Get DNA symbol
      tmp = &dBuf[idx-1];

      cModelHigh->pModelIdx = ((cModelHigh->pModelIdx - *(tmp-ctxHigh) *
      multiplierHigh) * nSym) + *tmp;

      if(mode == 0)
        {
        aCountersHigh = &cModelHigh->counters[nSym * cModelHigh->pModelIdx];
        tmpFreqsHigh[nSymL1] = 1 + alphaDEv * aCountersHigh[nSymL1];
        for(n = nSymL1 ; n-- ; )
          tmpFreqsHigh[n] = tmpFreqsHigh[n+1]+1+alphaDEv * aCountersHigh[n];
        cumHigh += (SearchLog(tmpFreqsHigh[0]) - SearchLog(tmpFreqsHigh[sym] 
        - tmpFreqsHigh[sym+1]));
        if(++aCountersHigh[sym] == maxCHigh)    // Check counter overflow&upd.
          for(n = nSym ; n-- ; )                    // For all stored counters
            aCountersHigh[n] >>= 1;                    // Divide counters by 2
        if(d->L.ir == 1)                                   // Inverted repeats
          {
          aCountersHigh = &cModelHigh->counters[nSym * (cModelHigh->idxIr = 
          cModelHigh->idxIr/nSym + GetComp2(*(tmp+1)) * multiplierHigh)];
          symC = GetComp2(*(tmp+1-ctxHigh));
          if(++aCountersHigh[symC] == maxCHigh)   //Check counter overflow&upd.
            for(n = nSym ; n-- ; )                  // For all stored counters
              aCountersHigh[n] >>= 1;                  // Divide counters by 2
          }
        }
      else
        {
        if(!(hCountersHigh = GetHCCounters(&cModelHigh->hTable,
        cModelHigh->pModelIdx)))
          hCountersHigh = (HCCounter *) cModelHigh->hTable.zeroCounters;
        tmpFreqsHigh[nSymL1] = 1 + alphaDEv * hCountersHigh[nSymL1];
        for(n = nSymL1 ; n-- ; )
          tmpFreqsHigh[n] = tmpFreqsHigh[n+1]+1+alphaDEv * hCountersHigh[n];
        cumHigh += (SearchLog(tmpFreqsHigh[0]) - SearchLog(tmpFreqsHigh[sym]
        - tmpFreqsHigh[sym+1]));
        UpdateHashCounter(cModelHigh, cModelHigh->pModelIdx, sym); // Update 
        if(d->L.ir == 1)                                 // Inverted repeats
          {
          cModelHigh->idxIr = cModelHigh->idxIr/nSym + GetComp2(*(tmp+1)) * 
          multiplierHigh;
          symC = GetComp2(*(tmp+1-ctxHigh));
          UpdateHashCounter(cModelHigh, cModelHigh->idxIr, symC);  // Update
          }
        }
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      // LOWER MODEL - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      aCountersLow = &cModelLow->counters[nSym * (cModelLow->pModelIdx = 
      ((cModelLow->pModelIdx - *(tmp-ctxLow) * multiplierLow) *nSym) + *tmp)];
      tmpFreqsLow[nSymL1] = aCountersLow[nSymL1] + 1;
      for(n = nSymL1 ; n-- ; )
        tmpFreqsLow[n] = tmpFreqsLow[n+1] + aCountersLow[n] + 1;
      cumLow += (SearchLog(tmpFreqsLow[0]) - SearchLog(tmpFreqsLow[sym] -
      tmpFreqsLow[sym+1]));
      if(++aCountersLow[sym] == maxCLow) 
        for(n = nSym ; n-- ; )                 
          aCountersLow[n] >>= 1;     
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      if(++idx == BUF_SIZE)
        {
        memcpy(dBuf-dP, dBuf+idx-dP, dP);
        idx = 0;
        }

      if(++idxBlockSymbol == blockSize)
        {
        #ifdef BLOCKTOTAL
        if(cumHigh < cumLow || block[nBlock] == 1)
          fprintf(BLKT, "%d\n", 1);
        else
          fprintf(BLKT, "%d\n", 0);
        #endif
        if(cumHigh < cumLow /* && block[nBlock] == 0*/)
          {
	  #ifdef BLOCKPRINT
          fprintf(BLKRL, "%d\n", 1);
          #endif
          block[nBlock] = 1;
          ++nLowComplex;
          }
        #ifdef BLOCKPRINT
        else
          fprintf(BLKRL, "%d\n", 0);
        #endif
        --nBlock;
        cumHigh         = 0;
        cumLow          = 0;
        idxBlockSymbol  = 0;
        }
      }

  // Flush Models & files
  FreeCModelCont(cModelLow);
  mode == 1 ? FreeCModelHash(cModelHigh) : FreeCModelCont(cModelHigh);
  fclose(INR);
  #ifdef BLOCKTOTAL
  fclose(BLKT);
  #endif
  #ifdef BLOCKPRINT
  fclose(BLKLR);
  fclose(BLKRL);
  #endif

  unlink(nameRev);

  ////////////////////////////////////////////////////////////////////////////
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // C O M P R E S S I O N
  //

  // Tune: alphaDen, mode, etc
  percLowBlock = (int)((double) nLowComplex / totalBlocks * 100);

/*
  //TODO: CALCULATE "BEST" PARAMETERS WITH ANALYSIS INFO.

  if(alphaDen == 0)
    {
    if(percLowBlock >= 95)
      {
      alphaDen = 50;
      if(mode == 1)
        ctxUsed  = ctxHigh + 4;
      else
        ctxUsed  = ctxHigh + 3;
      }
    else if(percLowBlock >= 90)
      {
      alphaDen = 30; 
      mode     = 1;
      ctxUsed  = ctxHigh + 2;
      }
    else if(percLowBlock >= 80)
      alphaDen = 20;
    else if(percLowBlock >= 70)
      alphaDen = 15;
    else
      alphaDen = 10;
    }
*/
  
  // HEADER
  sprintf(name,  "%s.LF%d", d->fn, id + 1);
  ac_encoder_init(&aced, name);                       // Start DNAbase Encoder

  ac_model_init(&acmHeader, 2);
  writeNBits(d->P.iPos,     46, &aced, &acmHeader);
  writeNBits(d->P.ePos,     46, &aced, &acmHeader);
  writeNBits(alphaDen,      16, &aced, &acmHeader);
  writeNBits(ctxLow,         8, &aced, &acmHeader);
  writeNBits(ctxUsed,        8, &aced, &acmHeader);
  writeNBits(maxCLow,       32, &aced, &acmHeader);
  writeNBits(maxCUsed,      32, &aced, &acmHeader);
  writeNBits(nSym,           8, &aced, &acmHeader);
  for(n = 0 ; n != nSym ; ++n)
    writeNBits(norAlp[n],    8, &aced, &acmHeader);
  ac_model_done(&acmHeader);

  if(verbose == 1)
    {
    printf("Running compression on thread %u\nUsing %u (%d %%) low entropy " 
    "blocks of %u\nAlpha denominator .............. %d\n", id+1, nLowComplex, 
    percLowBlock, totalBlocks, alphaDen);
    }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // SEQUENCE
  ac_model_init(&acmdh, 2);                   // Start DNAbase header ac model
  ac_model_init(&acmd,  nSym);                       // Start DNAbase ac model

  //CREATE MODELS & RESET/FLUSH VARIABLES 
  buf            = (uint8_t *) Calloc(BUF_SIZE,      sizeof(uint8_t));
  dBuf           = (uint8_t *) Calloc(BUF_SIZE + dP, sizeof(uint8_t));
  dBuf          += dP;
  cModelHigh     = CreateCModel(ctxUsed, nSym, alphaDen, mode, maxCUsed);
  cModelLow      = CreateCModel(ctxLow,  nSym,        1,    0, maxCLow);

  multiplierHigh = cModelHigh->multipliers[ctxUsed - 1];
  idxBlockSymbol = 0;
  idx            = 0;
  nBlock         = 0;
  idxSymbol      = d->P.iPos;
  fseek(IN, d->P.iPos, SEEK_SET);

  if(block[nBlock] == 1)
    acEncodeBin0(&aced, &acmdh);
  else
    {
    acEncodeBin1(&aced, &acmdh);
    ++acmdh.cfreq[1];
    }
  if(++acmdh.cfreq[0] == 100000000) //TODO: THIS NEEDS TO BE CHANGED IF IT
    {      // GIVES ERROR (INCREASING VALUE), PROBLEM: UPDATE MUST BE 
    acmdh.cfreq[1] >>= 1; // CORRECTED... MAYBE WE CAN ALSO USE A FCM FOR THIS
    acmdh.cfreq[0] >>= 1; // DONT FORGET DECODER!
    }

  while((k = fread(buf, 1, BUF_SIZE, IN)))
    for(idxPos = 0 ; idxPos != k ; ++idxPos)
      if(++idxSymbol < ePosP1)
        {
        sym = *(buf+idxPos);
        dBuf[idx] = sym = alp[sym];                              // Get symbol
        tmp = &dBuf[idx-1];

        cModelHigh->pModelIdx = ((cModelHigh->pModelIdx - *(tmp-ctxUsed) *
        multiplierHigh) * nSym) + *tmp;

        aCountersLow = &cModelLow->counters[nSym * (cModelLow->pModelIdx =
        ((cModelLow->pModelIdx - *(tmp-ctxLow) * multiplierLow) * nSym) +
        *tmp)];        

        if(block[nBlock] == 0)
          {
          acmd.cfreq[nSymL1] = aCountersLow[nSymL1] + 1;
          for(n = nSymL1 ; n-- ; )
            acmd.cfreq[n] = acmd.cfreq[n+1] + aCountersLow[n] + 1;
          ac_encode_symbol(&aced, &acmd, sym);                // Encode symbol
          }
        else
          {
          if(mode == 0)
            {
            aCountersHigh = &cModelHigh->counters[nSym*cModelHigh->pModelIdx];
            acmd.cfreq[nSymL1] = 1 + alphaDen * aCountersHigh[nSymL1];
            for(n = nSymL1 ; n-- ; )
              acmd.cfreq[n] = acmd.cfreq[n+1] +1+ alphaDen * aCountersHigh[n];
            ac_encode_symbol(&aced, &acmd, sym);              // Encode symbol
            if(++aCountersHigh[sym] == maxCHigh)  // Check counter overflow&up
              for(n = nSym ; n-- ; )                // For all stored counters
                aCountersHigh[n] >>= 1;                // Divide counters by 2
            }
          else
            {
            if(!(hCountersHigh = GetHCCounters(&cModelHigh->hTable,
            cModelHigh->pModelIdx)))
              hCountersHigh = (HCCounter *) cModelHigh->hTable.zeroCounters;
            acmd.cfreq[nSymL1] = 1 + alphaDen * hCountersHigh[nSymL1];
            for(n = nSymL1 ; n-- ; )
              acmd.cfreq[n] = acmd.cfreq[n+1]+1 + alphaDen * hCountersHigh[n];
            ac_encode_symbol(&aced, &acmd, sym);              // Encode symbol
            UpdateHashCounter(cModelHigh, cModelHigh->pModelIdx, sym);   // Up
            }
          }

         if(d->L.ir == 1)                                  // Inverted repeats
          {
          cModelHigh->idxIr = cModelHigh->idxIr/nSym + GetComp2(*(tmp+1)) * 
          multiplierHigh;
          if(mode == 0)
            {
            if(block[nBlock] != 0)
              {
              aCountersHigh = &cModelHigh->counters[nSym * cModelHigh->idxIr];
              symC = GetComp2(*(tmp+1-ctxUsed));
              if(++aCountersHigh[symC] == maxCHigh)    //Check counter overflow
                for(n = nSym ; n-- ; )              // For all stored counters
                  aCountersHigh[n] >>= 1;              // Divide counters by 2
              }
            }
          else
            {
            if(block[nBlock] != 0)
              {
              symC = GetComp2(*(tmp+1-ctxUsed));
              UpdateHashCounter(cModelHigh, cModelHigh->idxIr, symC);    // Up
              }
            }
          }

        if(++aCountersLow[sym] == maxCLow)  // Check counter overflow & update
          for(n = nSym ; n-- ; )                    // For all stored counters
            aCountersLow[n] >>= 1;                     // Divide counters by 2

        if(++idx == BUF_SIZE)
          {
          memcpy(dBuf-dP, dBuf+idx-dP, dP);
          idx = 0;
          }

        if(++idxBlockSymbol == blockSize)
          {
          if(block[++nBlock] == 1)
            acEncodeBin0(&aced, &acmdh);
          else
            {
            acEncodeBin1(&aced, &acmdh);
            ++acmdh.cfreq[1];
            }
          if(++acmdh.cfreq[0] == 100000000) //FIXME: SEE UP...
            {
            acmdh.cfreq[1] >>= 1;
            acmdh.cfreq[0] >>= 1;
            }
          idxBlockSymbol = 0;
          }
        }

  // Flush Models & files
  FreeCModelCont(cModelLow);
  mode == 1 ? FreeCModelHash(cModelHigh) : FreeCModelCont(cModelHigh);
  Free(block);
  fclose(IN);
  ac_model_done(&acmdh);
  ac_model_done(&acmd);
  ac_encoder_done(&aced);

  d->S.bits += ac_encoder_bits(&aced);
  d->S.size  = size;

  pthread_exit(NULL);
  }

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - D E C O M P R E S S O R - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void *Decompress(void *data)
  {
  Data        *d = (Data *) data;
  FILE        *OUT = NULL;
  char        name[MAX_FILESIZE], nameOut[MAX_FILESIZE];
  ac_model    acmHeader, acmdh, acmd;
  ac_decoder  acdd;
  HCCounter   *hCountersHigh;
  ACCounter   *aCountersHigh, *aCountersLow;
  CModel      *cModelHigh, *cModelLow;
  uint64_t    iPos, ePos, ePosP1, idxSymbol, multiplierHigh, multiplierLow;
  uint32_t    n, dP, nSym, nSymL1, ctxHigh, 
              ctxLow, nThreads, blockSize, tmpModel, alphaDen, maxCLow, 
              maxCHigh;
  int32_t     dIdx, idxBlockSymbol;
  uint8_t     norAlp[NSIGNF], id, *bufOD, *dBuf, *tmp, sym, symC, mode;

  id          = d->id;
  blockSize   = d->L.blockSize;
  nThreads    = d->nThreads;

  sprintf(nameOut, "%s.LF.tmp%d", d->fn, id + 1);
  OUT = Fopen(nameOut, "w");                          // OUTPUT TEMPORARY file
  sprintf(name, "%s.LF%d", d->fn, id + 1);
  ac_decoder_init(&acdd, name);

  ac_model_init(&acmHeader, 2);
  iPos           = readNBits(46, &acdd, &acmHeader);
  ePos           = readNBits(46, &acdd, &acmHeader);
  alphaDen       = readNBits(16, &acdd, &acmHeader);
  ctxLow         = readNBits( 8, &acdd, &acmHeader);
  ctxHigh        = readNBits( 8, &acdd, &acmHeader);
  maxCLow        = readNBits(32, &acdd, &acmHeader);
  maxCHigh       = readNBits(32, &acdd, &acmHeader);
  nSym           = readNBits( 8, &acdd, &acmHeader);
  for(n = 0 ; n != nSym ; ++n)
    norAlp[n]    = readNBits( 8, &acdd, &acmHeader);
  ac_model_done(&acmHeader);

//  // small sequence ->  force mode hash
//  if(ePos - iPos < 15000000) // 15M
//    mode = (d->L.mode == 2 ? 1 : d->L.mode);
//  else
//    mode = (d->L.mode == 2 ? 0 : d->L.mode);
  
  mode = (d->L.mode == 2 ? 0 : d->L.mode);

  tmp            = 0;
  dP             = L_GUARD_BUF;
  bufOD          = (uint8_t *) Calloc(MAX_OUT_LS,    sizeof(uint8_t));
  dBuf           = (uint8_t *) Calloc(BUF_SIZE + dP, sizeof(uint8_t));
  dIdx           = 0;
  nSymL1         = nSym - 1;
  dBuf          += dP;
  sym            = 0;
  cModelHigh     = CreateCModel(ctxHigh, nSym, alphaDen, mode, maxCHigh);
  cModelLow      = CreateCModel(ctxLow,  nSym,        1,    0, maxCLow);
  multiplierHigh = cModelHigh->multipliers[ctxHigh - 1]; 
  multiplierLow  = cModelLow->multipliers[ctxLow   - 1];

  ac_model_init(&acmdh, 2);                      // Start DNAbase Header model
  ac_model_init(&acmd,  nSym);                       // Start DNAbase ac model

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
  idxBlockSymbol = 0;
  idxSymbol      = iPos;
  ePosP1         = ePos + 1;

  tmpModel = acDecodeBinary(&acdd, &acmdh);
  if(tmpModel == 1)
    ++acmdh.cfreq[1];
  if(++acmdh.cfreq[0] == 100000000) //FIXME: SEE UP...
    {
    acmdh.cfreq[1] >>= 1;
    acmdh.cfreq[0] >>= 1;
    }

  while(idxSymbol < ePosP1)
    {
    tmp = &dBuf[dIdx-1];

    cModelHigh->pModelIdx = ((cModelHigh->pModelIdx - *(tmp-ctxHigh) *
    multiplierHigh) * nSym) + *tmp;

    aCountersLow = &cModelLow->counters[nSym * (cModelLow->pModelIdx =
    ((cModelLow->pModelIdx - *(tmp-ctxLow) * multiplierLow) * nSym) + *tmp)];

    if(tmpModel != 0)
      {
      acmd.cfreq[nSymL1] = aCountersLow[nSymL1] + 1;
      for(n = nSymL1 ; n-- ; )
        acmd.cfreq[n] = acmd.cfreq[n+1] + aCountersLow[n] + 1;
      dBuf[dIdx] = sym = acDecSymHighSizeVar(&acdd, &acmd);   // Decode symbol
      }
    else
      {
      if(mode == 0)
        {
        aCountersHigh = &cModelHigh->counters[nSym*cModelHigh->pModelIdx];
        acmd.cfreq[nSymL1] = 1 + alphaDen * aCountersHigh[nSymL1];
        for(n = nSymL1 ; n-- ; )
          acmd.cfreq[n] = acmd.cfreq[n+1] +1+ alphaDen * aCountersHigh[n];
        dBuf[dIdx] = sym = acDecSymHighSizeVar(&acdd, &acmd);        // Decode 
        if(++aCountersHigh[sym] == maxCHigh)      // Check counter overflow&up
          for(n = nSym ; n-- ; )                    // For all stored counters
            aCountersHigh[n] >>= 1;                    // Divide counters by 2
        }
      else
        {
        if(!(hCountersHigh = GetHCCounters(&cModelHigh->hTable,
        cModelHigh->pModelIdx)))
          hCountersHigh = (HCCounter *) cModelHigh->hTable.zeroCounters;
        acmd.cfreq[nSymL1] = 1 + alphaDen * hCountersHigh[nSymL1];
        for(n = nSymL1 ; n-- ; )
          acmd.cfreq[n] = acmd.cfreq[n+1]+1 + alphaDen * hCountersHigh[n];
        dBuf[dIdx] = sym = acDecSymHighSizeVar(&acdd, &acmd);        // Decode 
        UpdateHashCounter(cModelHigh, cModelHigh->pModelIdx, sym);    // Up HT
        }
      }

    if(d->L.ir == 1)
      {
      cModelHigh->idxIr = (cModelHigh->idxIr/nSym) + GetComp2(*(tmp+1)) * 
      multiplierHigh;
      if(tmpModel == 0)
        {
        if(mode == 0)
          {
          aCountersHigh = &cModelHigh->counters[nSym * cModelHigh->idxIr];
          symC = GetComp2(*(tmp+1-ctxHigh));
          if(++aCountersHigh[symC] == maxCHigh)    //Check counter overflow
            for(n = nSym ; n-- ; )              // For all stored counters
              aCountersHigh[n] >>= 1;              // Divide counters by 2
          }
        else
          {
          symC = GetComp2(*(tmp+1-ctxHigh)); // alp[Comp(norAlp[*(tmp+1-ctx)])];
          UpdateHashCounter(cModelHigh, cModelHigh->idxIr, symC);      // Update 
          }
        }
      }

    if(++aCountersLow[sym] == maxCLow)      // Check counter overflow & update
      for(n = nSym ; n-- ; )                        // For all stored counters
        aCountersLow[n] >>= 1;                         // Divide counters by 2
      
    if(++dIdx == BUF_SIZE)
      {
      memcpy(dBuf-dP, dBuf+dIdx-dP, dP);
      dIdx = 0;
      }

    bufOD[idxBlockSymbol] = norAlp[sym];

    if(++idxBlockSymbol == blockSize)
      {
      fwrite(bufOD, 1, idxBlockSymbol, OUT);
      idxBlockSymbol = 0;
      tmpModel = acDecodeBinary(&acdd, &acmdh);
      if(tmpModel == 1)
        ++acmdh.cfreq[1];
      if(++acmdh.cfreq[0] == 100000000) //FIXME: SEE UP...
        {
        acmdh.cfreq[1] >>= 1;
        acmdh.cfreq[0] >>= 1;
        }
      }
    ++idxSymbol;
    }

  if(idxBlockSymbol != 1)
    {
    if(id == nThreads - 1)
      fwrite(bufOD, 1, idxBlockSymbol - 1, OUT);
    else
      fwrite(bufOD, 1, idxBlockSymbol, OUT);
    }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  FreeCModelCont(cModelLow);
  mode == 1 ? FreeCModelHash(cModelHigh) : FreeCModelCont(cModelHigh);
  fclose(OUT);

  ac_model_done(&acmdh);
  ac_model_done(&acmd);
  ac_decoder_done(&acdd);
  
  pthread_exit(NULL);
  }

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - - -  P A R T S  - - - - - - - - - - - - - - - -

void LoadPosForThreads(Data *data, uint32_t nThreads)
  {
  uint64_t  k;
  uint8_t   n;

  k = UINT64_MAX;
  for(n = 0 ; n != nThreads ; ++n)
    {
    data[n].P.iPos = k + 1;
    data[n].P.ePos = k += (data[0].size / nThreads + 1);
    }
  data[n-1].P.ePos = data[0].size;
  }

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - P R I N T   S T A T S - - - - - - - - - - - - -

void PrintStats(Data *data, uint64_t headerBits)
  {
  uint8_t   n;
  uint64_t  bits = 0, size = 0;

  for(n = 0 ; n != data->nThreads ; ++n)
    {
    bits += data[n].S.bits;
    size += data[n].S.size;
    }

  printf("------------------------------------------------------\n");
  printf("| TOTAL   | %-13"PRIu64" | %-13"PRIu64" | %-8.4lf |\n", size, bits /
  8, (double) bits / 8 / size);
  printf("------------------------------------------------------\n");
  }

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - P R I N T   P O S I T I O N S - - - - - - - - - - -

void printPositions(Data *data)
  {
  uint8_t n;

  printf("------------------------------------------------------\n");
  printf("| THREAD  |  INITIAL POSITION  |   FINAL POSITION    |\n");
  printf("------------------------------------------------------\n");
  for(n = 0 ; n != data->nThreads ; ++n)
    printf("| %-7d | %-18"PRIu64" | %-19"PRIu64" |\n", n + 1, data[n].P.iPos, 
    data[n].P.ePos);
  printf("------------------------------------------------------\n");
  }
      

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - M A I N - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

int main(int argc, char *argv[])
  {	
  char        name[MAX_FILESIZE], file[MAX_FILESIZE];
  uint32_t    n, k, action, help, nThreads, ctxLow, ctxEval, ctxUsed, bSize, 
              ir, alphaDen, alphaDEv, maxCountLow, maxCountUsed, mode, rm;
  Data        *data;
  ac_model    acm;

  help         = 0;
  verbose      = 0;
  ctxLow       = 3;
  ctxEval      = 12;
  ctxUsed      = 12;
  alphaDen     = 10;                             // 0 = predicted by evaluation
  alphaDEv     = 10;
  maxCountLow  = MAXCNT;
  maxCountUsed = MAXCNT;
  bSize        = DEFAULT_BS;
  mode         = 2;
  nThreads     = 1;
  ir           = 0;
  action       = 0;
  rm           = 0;
    
  for(n = 1 ; n != argc ; ++n)
    if(!strcmp("-h", argv[n]))
      { help = 1; break; }

  if(argc < 2 || help == 1)
    {
    fprintf(stderr, "Usage: HighFCM [OPTION]... [FILE]                     \n\n");
    fprintf(stderr, " -h               give this help                        \n");
    fprintf(stderr, " -v               verbose mode                          \n");
    fprintf(stderr, " -cl  <ctxLow>    low context order used in compression \n");
    fprintf(stderr, " -ml  <maxCnt>    low order maximum counter             \n");
    fprintf(stderr, " -ce  <ctxEval>   high context order on evaluation      \n");
    fprintf(stderr, " -cu  <ctxUsed>   high context order on compression     \n");
    fprintf(stderr, " -mu  <maxCnt>    used order maximum counter            \n");
    fprintf(stderr, " -au  <alpha>     alpha estimator denominator for cu    \n");
    fprintf(stderr, " -ae  <alpha>     alpha estimator denominator for ce    \n");
    fprintf(stderr, " -b   <blockSize> block size (default: %d)\n", DEFAULT_BS  );
    fprintf(stderr, " -ir              use inverted repeats                  \n");
    fprintf(stderr, " -tm  <tableMode> table mode: 0|1 (0=array, 1=hash)     \n");
    fprintf(stderr, " -t   <nThreads>  number of threads / parts             \n");
    fprintf(stderr, " -d   <outFile>   decompression output file             \n");
    fprintf(stderr, " -rm              remove comp file after decomp         \n");
    fprintf(stderr, " <File>           input file to compress              \n\n");
    return EXIT_SUCCESS;
    }

  for(n = 1 ; n != argc ; ++n)
    if(!strcmp("-v", argv[n]))
      { verbose  = 1; break; }

  for(n = 1 ; n != argc ; ++n)
    if(!strcmp("-cl", argv[n]))
      { ctxLow   = atoi(argv[n+1]); break; }

  for(n = 1 ; n != argc ; ++n)
    if(!strcmp("-ml", argv[n]))
      { maxCountLow = atoi(argv[n+1]); break; }

  for(n = 1 ; n != argc ; ++n)
    if(!strcmp("-ce", argv[n]))
      { ctxEval  = atoi(argv[n+1]); break; }

  for(n = 1 ; n != argc ; ++n)
    if(!strcmp("-cu", argv[n]))
      { ctxUsed  = atoi(argv[n+1]); break; }

  for(n = 1 ; n != argc ; ++n)
    if(!strcmp("-mu", argv[n]))
      { maxCountUsed = atoi(argv[n+1]); break; }

  for(n = 1 ; n != argc ; ++n)
    if(!strcmp("-au", argv[n]))
      { alphaDen = atoi(argv[n+1]); break; }

  for(n = 1 ; n != argc ; ++n)
    if(!strcmp("-ae", argv[n]))
      { alphaDEv = atoi(argv[n+1]); break; }

  for(n = 1 ; n != argc ; ++n)
    if(!strcmp("-b", argv[n]))
      { bSize    = atoi(argv[n+1]); break; }

  for(n = 1 ; n != argc ; ++n)
    if(!strcmp("-ir", argv[n]))
      { ir       = 1; break; }

  for(n = 1 ; n != argc ; ++n)
    if(!strcmp("-tm", argv[n]))
      { mode     = atoi(argv[n+1]); break; }

  for(n = 1 ; n != argc ; ++n)
    if(!strcmp("-t", argv[n]))
      { nThreads = atoi(argv[n+1]); break; }

  for(n = 1 ; n != argc ; ++n)
    if(!strcmp("-d", argv[n]))
      { action   = 1; strcpy(file, argv[n+1]); break; }

  for(n = 1 ; n != argc ; ++n)
    if(!strcmp("-rm", argv[n]))
      { rm      = 1; break; }


  if(!action)
    {
    #ifdef THREADS
    pthread_t   t[nThreads];
    ac_encoder  ace;
    uint64_t    headerBits = 0, size;

    data        = (Data *) Calloc(nThreads, sizeof(Data));
    size        = NBytesInFile(argv[argc-1]);

    for(n = 0 ; n != nThreads ; ++n)
      {
      strcpy(data[n].fn, argv[argc-1]);
      data[n].nThreads    = nThreads;
      data[n].L.lowOrder  = ctxLow;
      data[n].L.maxCLow   = maxCountLow;
      data[n].L.highOrder = ctxEval;      
      data[n].L.usedOrder = ctxUsed;
      data[n].L.maxCUsed  = maxCountUsed; 
      data[n].L.alphaDen  = alphaDen;
      data[n].L.alphaDEv  = alphaDEv;
      data[n].L.blockSize = bSize;
      data[n].L.mode      = mode;
      data[n].L.ir        = ir;
      data[n].size        = size;
      data[n].id          = n;
      }

    LoadPosForThreads(data, nThreads);
    if(verbose)
      {
      printPositions(data);
      printf("Number of threads .............. %d\n", nThreads);
      printf("Low context order .............. %d\n", ctxLow);
      printf("Low order maxCount ............. %d\n", maxCountLow);
      printf("Eval high context order ........ %d\n", ctxEval);
      printf("Used high context order ........ %d\n", ctxUsed);
      printf("Used order maxCount ............ %d\n", maxCountUsed);
      printf("Block size ..................... %d\n", bSize);
      printf("Inverted repeats ............... %s\n", ir == 1 ? "yes" : "no");
      printf("Number of bytes ................ %"PRIu64"\n", size);
      }

    sprintf(name, "%s.LF0", data->fn); // Header - - - - - - - - - - - - - - - 
    ac_encoder_init(&ace, name); 
    ac_model_init(&acm, 2);

    writeNBits(nThreads,     12, &ace, &acm);
    writeNBits(bSize,        30, &ace, &acm);
    writeNBits(ir,            1, &ace, &acm);

    ac_model_done(&acm);
    ac_encoder_done(&ace);
    headerBits += ac_encoder_bits(&ace);  // - - - - - - - - - - - - - - - - -

    for(n = 0 ; n != nThreads ; ++n)  // solo 4: mk shore it doesn't interfere
      pthread_create(&(t[n+1]), NULL, Compress, (void *) &(data[n]));

    for(n = 0 ; n != nThreads ; ++n)                        // JOIN AND OUTPUT
      pthread_join(t[n+1], NULL);

    #else
    fprintf(stderr, "ERROR: not available!\n");             // TODO: implement
    exit(1);
    #endif	
    printf("Compression done!\n");
    PrintStats(data, headerBits);
    }
  else
    {
    #ifdef THREADS  // - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    
    FILE        *OUT = NULL;
    ac_decoder  acd;
    uint8_t     *buf;

    // Read header msg - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    sprintf(name, "%s.LF0", argv[argc-1]);
    ac_decoder_init(&acd, name);
    ac_model_init(&acm, 2);

    nThreads     = readNBits(12, &acd, &acm);
    bSize        = readNBits(30, &acd, &acm);
    ir           = readNBits( 1, &acd, &acm);
    
    data = (Data *) Calloc(nThreads, sizeof(Data));
    for(n = 0 ; n != nThreads ; ++n)
      {
      data[n].nThreads    = nThreads;
      data[n].id          = n;
      strcpy(data[n].fn, argv[argc-1]);
      data[n].L.mode      = mode;
      data[n].L.blockSize = bSize;
      data[n].L.ir        = ir; 
      }

    ac_model_done(&acm);
    ac_decoder_done(&acd);  // - - - - - - - - - - - - - - - - - - - - - - - -

    pthread_t  t[nThreads];

    for(n = 0 ; n != nThreads ; ++n)  // solo 4: mk shore it doesn't interfere
      pthread_create(&(t[n+1]), NULL, Decompress, (void *) &(data[n]));
        
    OUT = Fopen(file, "w");             // CONCATENATION OF DECOMPRESSED FILES
    buf = (uint8_t *) Malloc(BUF_SIZE_CAT * sizeof(uint8_t));
    for(n = 0 ; n != nThreads ; ++n)                        // JOIN AND OUTPUT
      {
      pthread_join(t[n+1], NULL);
      char tmp[MAX_FILESIZE], tmp2[MAX_FILESIZE];
      sprintf(tmp,  "%s.LF.tmp%d", data[n].fn, n + 1);
      sprintf(tmp2, "%s.LF%d",     data[n].fn, n);
      FILE *IN = Fopen(tmp, "r");
      while((k = fread(buf, 1, BUF_SIZE_CAT, IN)))
        fwrite(buf, 1, k, OUT); 	
      fclose(IN);
      if(rm == 1)
        {
        unlink(tmp);
        unlink(tmp2);
        }
      }
    char tmp2[MAX_FILESIZE];
    sprintf(tmp2, "%s.LF%d", data[n-1].fn, n);
    if(rm == 1)
      unlink(tmp2);

    #else
    fprintf(stderr, "ERROR: not available!\n");             // TODO: implement
    exit(1);
    #endif
    printf("Decompression done!\n");
    }

  return EXIT_SUCCESS;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

