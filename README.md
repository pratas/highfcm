<h1>HighFCM</h1>
<p align="center"><img src="/logo.png" 
alt="HighFCM" width="260" height="260" border="0" /></p>
HighFCM is a compression algorithm that relies on a pre-analysis of the data, namely an algorithmic entropy filter, before compression, with the aim of filtering regions of low complexity. This strategy enables us to use deeper context models, supported by hash-tables, without requiring huge amounts of memory. As an example, context depths as large as 32 are attainable for alphabets of four symbols, as is the case of genomic sequences. These deeper context models show very high compression capabilities in very repetitive genomic sequences, yielding improvements over previous algorithms. Furthermore, this method is universal, in the sense that it can be used in any type of textual data.


<h2>INSTALLATION</h2>

In the following instructions we show the procedure to install/compile HighFCM.

<h3>Linux</h3>
<pre>
git clone https://github.com/pratas/highfcm.git
cd highfcm
make
</pre>

<h3>Windows</h3>

In windows use cygwin (https://www.cygwin.com/) and make sure that it is included in the installation: make, unzip, wget (and any dependencies). If you install the complete cygwin packet then all these will be installed. After, all steps will be the same as in Linux.

<h2>EXECUTION</h2>

Example on running HighFCM:

<pre>
./HighFCM -v -cl 4 -ce 14 -cu 16 File.seq
</pre>

<h2>PARAMETERS</h2>

To see the possible options type
<pre>
./HighFCM
</pre>
or
<pre>
./HighFCM -h
</pre>

These will print the following options:
<pre>
<p>
Usage: HighFCM [OPTION]... [FILE]                     

 -h               give this help                        
 -v               verbose mode                          
 -cl  &#60ctxLow&#62    low context order used in compression 
 -ml  &#60maxCnt&#62    low order maximum counter             
 -ce  &#60ctxEval&#62   high context order on evaluation      
 -cu  &#60ctxUsed&#62   high context order on compression     
 -mu  &#60maxCnt&#62    used order maximum counter            
 -au  &#60alpha&#62     alpha estimator denominator for cu    
 -ae  &#60alpha&#62     alpha estimator denominator for ce    
 -b   &#60blockSize&#62 block size (default: 100)
 -ir              use inverted repeats                  
 -tm  &#60tableMode&#62 table mode: 0|1 (0=array, 1=hash)     
 -t   &#60nThreads&#62  number of threads / parts             
 -d   &#60outFile&#62   decompression output file             
 -rm              remove comp file after decomp         
 &#60File&#62           input file to compress
</p>
</pre>

<h2>CITATION</h2>

On using this software/method please cite:

Pratas, D.; Pinho, A.J., "Exploring deep Markov models in genomic data compression using sequence pre-analysis", Signal Processing Conference (EUSIPCO), 2014 Proceedings of the 22nd European, pp.2395-2399, 1-5 Sept. 2014.

<h2>ISSUES</h2>

For any issue let us know at [issues link](https://github.com/pratas/highfcm/issues).

<h2>LICENSE</h2>

GPL v2.

For more information:
<pre>http://www.gnu.org/licenses/gpl-2.0.html</pre>

                                                    

