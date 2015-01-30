////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License (version 2) as
// published by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program;if not, write to the
//
// Free Software Foundation, Inc.,
// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//
////////////////////////////////////////////////////////////////

// This is the "squish it into memory" version of ProbeLists.
// example/example-ProbeListFactory.cpp has a examples of usage.

/// @todo Add a free list for freeing single ProbeLists?
/// @todo Use 3B/2B sizes for ProbeList counts?
/// @todo allow ProbeList to be malloced.

#ifndef _PROBELISTFACTORY_H_
#define _PROBELISTFACTORY_H_

// forward declare these
class ProbeList;
class ProbeListPacked;
//
class ProbeListHeap;
class ProbeListFactory;

#include "chipstream/ProbeList.h"
#include "chipstream/ProbeListStl.h"
#include "chipstream/ProbeSet.h"
//
#include "file/TsvFile/SpfFile.h"
//
#include <cstring>
#include <ctime>
#include <string>
#include <vector>
//



// Limits of what can be stored.
#define PROBELIST_BLOCK_CNT_MAX  UINT_MAX
#define PROBELIST_PROBE_CNT_MAX  UINT_MAX
#define PROBELIST_NAME_LEN_MAX   USHRT_MAX

#ifdef __GNUC__
#define PROBELIST_PACKED __attribute__ ((packed))
#else
#define PROBELIST_PACKED
#endif

// Set needed features for sparcs.  They cant do unaligned loads.
#ifdef __sparc__
#define PROBELIST_NEEDS_UNALIGNED
#endif

// For testing you might want to uncomment this even on a
// machine which does support unaligned reads.
// #define PROBELIST_NEEDS_UNALIGNED

// These macros are for abstracting out where unaligned accesses are happening.
// Be sure to match up the sizes of access with the sizes of the elements.
#ifdef PROBELIST_NEEDS_UNALIGNED

// This arch cant to unaligned stores and loads, so
// we need some helper functions to build the value from bytes.
// NOTE: These are all big-endian
unsigned int probelist_load32b(unsigned int* p);
unsigned short probelist_load16b(unsigned short* p);
void probelist_store32b(unsigned int* p,unsigned int v);
void probelist_store16b(unsigned short* p,unsigned short v);

// Now use them when we have possibly unaligned accesses.
#define PROBELIST_LOAD32B(p)    probelist_load32b((unsigned int*)p)
#define PROBELIST_STORE32B(p,v) probelist_store32b((unsigned int*)(p),v)
//
#define PROBELIST_LOAD16B(p)    probelist_load16b((unsigned short*)p)
#define PROBELIST_STORE16B(p,v) probelist_store16b((unsigned short*)(p),v)

#else

// This arch can do unaligned loads, so just deref the pointer
#define PROBELIST_LOAD32B(p) (*(p))
#define PROBELIST_LOAD16B(p) (*(p))
//
#define PROBELIST_STORE32B(p,v) *(p)=(v)
#define PROBELIST_STORE16B(p,v) *(p)=(v)

#endif

//
#define PROBELIST_LOAD8B(p)    (*(p))
#define PROBELIST_STORE8B(p,v) *(p)=(v)

// This has to be allocated from a ProbeListHeap as it does
// not free memory on allocation.
// But we could search the regions to see if we could find our
// pointer in the Factory and if not, free our memory.

/** \page file-format-spf File Format: SPF (NON-OFFICIAL-RELEASE)

Simple Probe Format (SPF) is a text representation of the probes for a particular
probeset as stored in a ProbeList class. The intention is that spf will describe
the structure of a probeset and that an additional information about probes
(i.e. gc content, genomic position, pcr fragment size, etc.) will be provided in
another file. The additional information is not to be included in the spf file
as 1) it would make the file much larger to have every possible additional field
and 2) given the innovation in both assay and algorithm development it is
impossible to predict all the fields that may be necessary so providing them in
an additional file allows them to by dynamic and flexible. The program
cdf-to-simple will convert cdf files and pgf/clf files into the spf format.

The file starts off with the obligatory header lines of chip_type,
num-probesets, num-cols and num-rows with a preceding '#%%' to indicate that the
line is a header line. Other header lines containing meta information are
permitted, but not required. After the header lines there is a single, tab
delimited, line that supplies the names of each column (more below on what each
column contains).  For example the HG-U133 array header would look like:

<pre>
#%spf-format=1
#%%chip_type=HG-U133_Plus_2
#%%num-probesets=54675
#%%num-cols=1164
#%%num-rows=1164
name	type	num_blocks	block_sizes	block_annotations	num_match	num_probes	probes
</pre>


Following this the bulk of the data is in the form of 8 tab delimited columns:

 - name - The name (sometimes called id) for a particular probeset. Examples
   include '1552314_a_at' from the U133 designs or 'SNP_A-4282894' from the 500K
   Sty array. [string]

 - type - The type of probeset as defined by ProbeSet::Type
   (Unknown=0, Expression=1, Genotype=2, Copynumber=5) [int]

 - num_blocks - How many blocks (also known as groups) are in the
   probeset. Expression and Copynumber probesets have a single block while
   Genotyping probesets have 2 (on per allele) or 4 (one per strand and allele
   combination) [int]

 - block_sizes - Comma separated list of the number of perfect match (PM) probes
   in each block. [ints]

 - block_annotations - Comma separated list with one annotation per
   block. Annotations are defined in the ProbeSet::BlockAnnotation
   enumeration. [ints]

 - num_match - How many matching probes are there for each PM probe in the
   probeset? If PM only probeset this will be 1 (the PM itself). If a probeset
   has mismatch (MM) probes then num_match will be 2. In theory for resequencing
   arrays this would be 4, but resequencing probesets are not supported. [int]

 - num_probes - Total number of probes present in the probeset. This value
   should be equal to the sum of the block_sizes multiplied by the num_match
   field. [int]

 - probes - Comma separated list of probe indexes into the array. Indexes are
   the same 1 based indexes as used in the Bioconductor affy package. This is
   equivalent to the FusionCELData::XYToIndex() function result with 1 added to
   each index. The XYToIndex() function calculates the index as (y*c) + x where
   y is the y coordinate of a probe, x is the x coordinate of the probe on the
   array and c is the number of columns on the array. For example to get the
   correct index for a SNP6 chip (2572 rows and 2680 cols) the index for probe
   at position x=1000, y=500 is 1 + (500 * 2680 + 1000) = 1341001. The value
   '-1' is used to represent the 'null probe' which indicates no probe for that
   place. Examples of the null probe usage include a probeset with only MM probes
   or where only certain PM probes have matching MM probes. [ints]

Lets look at a couple of examples to see what to expect. The simplest case is a
PM-only expression probeset:
<pre>
3394979	1	1	4	0	1	4	797737,3775011,1875296,4778657
</pre>

The probeset's name is: '3394979' and it has a single block with 4 PM probes
797737,3775011,1875296,4778657 that are used to generate a summary estimate.

Things get slightly more interesting when we look at an example of a HG-U133 probeset like '200090_at':
<pre>
200090_at	1	1	11	0	2	22	1146869,424878,458233,325661,1187125,387873,765617,614181,879885,358159,48602,1148033,426042,459397,326825,1188289,389037,766781,615345,881049,359323,49766
</pre>

Here again there is just one block, but now the probeset contains mismatches so
the num_match field is '2'. Conceptually the first half of the probe list is PM
probes and the second half is the MM probes in the matching order. This is
illustrated in Figure 1. So probes [1146869,424878,..,358159,48602] are the PM
probes and probes [1148033,426042,..,359323,49766] are the MM probes for this
probeset.  Since they are in the same order the PM probe 1146869 has the
corresponding MM probe 1148033. In fact the MM probe is always just one row away
which means that the abs(PM index - MM index) = 1164 for this design as the PM
probes were designed to be near the matching MM probes.

@image html pm-mm-matching.png "Figure 1: Layout of PM and MM probes in the ProbeList."

Genotyping probesets make extensive use of the blocks (aka groups), with a different block
for each allele and sometimes for each strand. Looking at the GenomeWideSNP_6 probeset 'SNP_A-1999704' we see:
<pre>
SNP_A-1999704	2	2	3,3	1,2	1	6	756318,104506,6862022,756317,104505,6862021
</pre>

The probeset only has PM probes, but now there are two blocks. One block is for
the A allele and the other is for the B allele. In this case each block has three
probe indexes with the A allele being the first three and the B allele probes
being the second three. The 500K probesets often have 4 blocks with a separate block
for the A allele forward strand and A allele reverse strand. Looking at 'SNP_A-4246078' as
an example:
<pre>
SNP_A-4246078	2	4	2,2,4,4	3,4,5,6	2	24	4595060,1345852,5001712,1549918,1737600,613594,849516,1535644,4217526,3887116,852076,5702732,4595061,1345853,5001713,1549919,1737601,613595,849517,1535645,4217527,3887117,852077,5702733
</pre>

This probeset has 4 blocks with the first two blocks interrogating the forward
strand and the next two blocks interrogating the reverse strand. This probeset
also contains mismatches as the num_match field is 2. Note that the block sizes
and annotations only described for the PM portion of the probelist. The
association of MM probes to the blocks is determined by which PM probe they are
mapped to. Figure 2 has a graphical representation of this structure.

@image html probelist-blocks.png "Figure 2: Blocks within the ProbeList"

<h3> Example files: </h3>

The ubiquitous U133 layout would look something like:
<pre>
#%%chip_type=HG-U133_Plus_2
#%%num-probesets=54675
#%%num-cols=1164
#%%num-rows=1164
name	type	num_blocks	block_sizes	block_annotations	num_match	num_probes	probes
AFFX-BioB-5_at	1	1	20	0	2	40	794423,794424,794425,794426,794427,794428,794429,794430,794431,794432,794433,794434,794435,794436,794437,794438,794439,794440,796751,796752,795587,795588,795589,795590,795591,795592,795593,795594,795595,795596,795597,795598,795599,795600,795601,795602,795603,795604,797915,797916
AFFX-BioB-M_at	1	1	20	0	2	40	796753,796754,796755,796756,796757,796758,796759,796760,796761,796762,796763,796764,796765,796766,796767,796768,799079,799080,799081,799082,797917,797918,797919,797920,797921,797922,797923,797924,797925,797926,797927,797928,797929,797930,797931,797932,800243,800244,800245,800246
AFFX-BioB-3_at	1	1	20	0	2	40	799083,799084,799085,799086,799087,799088,799089,799090,799091,799092,799093,799094,799095,799096,801407,801408,801409,801410,801411,801412,800247,800248,800249,800250,800251,800252,800253,800254,800255,800256,800257,800258,800259,800260,802571,802572,802573,802574,802575,802576
AFFX-BioC-5_at	1	1	20	0	2	40	801413,801414,801415,801416,801417,801418,801419,801420,801421,801422,801423,801424,803735,803736,803737,803738,803739,803740,803741,803742,802577,802578,802579,802580,802581,802582,802583,802584,802585,802586,802587,802588,804899,804900,804901,804902,804903,804904,804905,804906
....
200605_s_at	1	1	11	0	2	22	837879,840082,1293653,16148,1312065,73758,715973,530083,276984,939364,1068385,839043,841246,1294817,17312,1313229,74922,717137,531247,278148,940528,1069549
200606_at	1	1	11	0	2	22	857936,425744,977312,481153,1179648,18573,1250096,1046820,40870,201374,701899,859100,426908,978476,482317,1180812,19737,1251260,1047984,42034,202538,703063
200607_s_at	1	1	11	0	2	22	695607,1083772,625475,1107666,1047598,765756,623013,1295961,165217,558447,1000847,696771,1084936,626639,1108830,1048762,766920,624177,1297125,166381,559611,1002011
....
1570650_at	1	1	11	0	2	22	4499,1270617,639560,639611,564882,644325,478862,623171,544217,928472,302553,5663,1271781,640724,640775,566046,645489,480026,624335,545381,929636,303717
1570651_at	1	1	11	0	2	22	189741,749316,141043,1114388,234402,448338,1028876,1226035,1228039,29472,1216749,190905,750480,142207,1115552,235566,449502,1030040,1227199,1229203,30636,1217913
1570653_at	1	1	11	0	2	22	272359,732815,242066,60473,59602,1256817,1197933,935422,206562,1322096,590788,273523,733979,243230,61637,60766,1257981,1199097,936586,207726,1323260,591952
</pre>

The Human Exon Array would look like:
<pre>
#%spf-format=1
#%%chip_type=HuEx-1_0-st-v2
#%%chip_type=HuEx-1_0-st-v1
#%%num-probesets=1432154
#%%num-cols=2560
#%%num-rows=2560
name	type	num_blocks	block_sizes	block_annotations	num_match	num_probes	probes
2590411	1	1	4	0	1	4	5402769,4684894,3869021,3774604
2609210	1	1	4	0	1	4	4163051,2486218,5470815,4829084
2616957	1	1	4	0	1	4	5947753,6282330,6070258,6003449
2391254	1	1	4	0	1	4	1170562,5776597,1007224,2676375
3499280	1	1	2	0	1	2	3609563,1320879
...
3494658	1	1	2	0	1	2	2549720,2170788
3022880	1	1	4	0	1	4	4970078,3773806,2406807,5900952
3961097	1	1	4	0	1	4	642050,647044,3391849,2835665
3352569	1	1	1	0	1	1	4385208
2711465	1	1	4	0	1	4	5378569,2578444,5547545,4583128
...
2441106	1	1	4	0	1	4	6333936,1089208,2166302,6001459
3121716	1	1	4	0	1	4	3415213,3553653,2490910,1682703
3478390	1	1	4	0	1	4	5974553,4128822,3849466,2950916
</pre>

The GenomeWideSNP_6 would look like:

<pre>
#%spf-format=1
#%%chip_type=GenomeWideSNP_6
#%%num-probesets=1856069
#%%num-cols=2680
#%%num-rows=2572
name	type	num_blocks	block_sizes	block_annotations	num_match	num_probes	probes
AFFX-5Q-123	1	1	30	0	2	60	3466655,3463975,3461295,3458615,3455935,3453255,3450575,3447895,3445215,3442535,3439855,3437175,3434495,3431815,3429135,3466657,3463977,3461297,3458617,3455937,3453257,3450577,3447897,3445217,3442537,3439857,3437177,3434497,3431817,3429137,3466656,3463976,3461296,3458616,3455936,3453256,3450576,3447896,3445216,3442536,3439856,3437176,3434496,3431816,3429136,3466658,3463978,3461298,3458618,3455938,3453258,3450578,3447898,3445218,3442538,3439858,3437178,3434498,3431818,3429138
...
SNP_A-2114462	2	2	3,3	1,2	1	6	2761136,3022142,6499884,2761135,3022141,6499883
SNP_A-4224354	2	2	4,4	1,2	1	8	5246520,1846798,2807616,4012714,5246519,1846797,2807615,4012713
SNP_A-2135352	2	2	3,3	1,2	1	6	4989424,145436,950490,4989423,145435,950489
SNP_A-4226297	2	2	3,3	1,2	1	6	4787390,4906380,2552122,4787389,4906379,2552121
SNP_A-4285774	2	2	4,4	1,2	1	8	672608,6354830,5244470,2034834,672607,6354829,5244469,2034833
...
CN_943507	5	1	1	0	1	1	656599
CN_943508	5	1	1	0	1	1	2308895
CN_943512	5	1	1	0	1	1	3208252
CN_954736	5	1	1	0	1	1	1818309
</pre>

The Mapping250K_Sty array would look like:
<pre>
#%spf-format=1
#%%chip_type=Mapping250K_Sty
#%%num-probesets=238378
#%%num-cols=2560
#%%num-rows=2560
name	type	num_blocks	block_sizes	block_annotations	num_match	num_probes	probes
AFFX-5Q-123	2	1	30	0	2	60	3295979,3293419,3290859,3288299,3285739,3283179,3280619,3278059,3275499,3272939,3270379,3267819,3265259,3262699,3260139,3295981,3293421,3290861,3288301,3285741,3283181,3280621,3278061,3275501,3272941,3270381,3267821,3265261,3262701,3260141,3295980,3293420,3290860,3288300,3285740,3283180,3280620,3278060,3275500,3272940,3270380,3267820,3265260,3262700,3260140,3295982,3293422,3290862,3288302,3285742,3283182,3280622,3278062,3275502,3272942,3270382,3267822,3265262,3262702,3260142
...
SNP_A-1871949	2	4	2,2,4,4	3,4,5,6	2	24	3076108,5844012,4310028,5841452,1055318,525766,2307904,2487962,1057878,528326,2310464,2490522,3076109,5844013,4310029,5841453,1055319,525767,2307905,2487963,1057879,528327,2310465,2490523
SNP_A-4276677	2	4	6,6,4,4	3,4,5,6	2	40	5803256,3092454,1461742,3673504,5672724,1804524,5800696,3089894,1464302,3676064,5670164,1807084,4234238,2618838,3110246,3335472,4231678,2621398,3107686,3338032,5803257,3092455,1461743,3673505,5672725,1804525,5800697,3089895,1464303,3676065,5670165,1807085,4234239,2618839,3110247,3335473,4231679,2621399,3107687,3338033
SNP_A-4273756	2	4	5,5,5,5	3,4,5,6	2	40	5329184,4571594,3219756,1227684,4012936,5331744,4569034,1976094,1225124,4015496,5664532,1543118,2151874,2533228,4068804,5661972,4584806,2149314,2535788,4066244,5329185,4571595,3219757,1227685,4012937,5331745,4569035,1976095,1225125,4015497,5664533,1543119,2151875,2533229,4068805,5661973,4584807,2149315,2535789,4066245
SNP_A-1808347	2	4	3,3,3,3	3,4,5,6	2	24	3183678,2069336,2156434,3186238,2066776,2158994,3107090,3964680,2404970,3104530,3967240,2402410,3183679,2069337,2156435,3186239,2066777,2158995,3107091,3964681,2404971,3104531,3967241,2402411
...
SNP_A-2291495	2	4	2,2,4,4	3,4,5,6	2	24	1197860,3652502,4696888,2210522,1598696,271122,4115920,2267510,492096,35098,1849638,1537196,1197861,3652503,4696889,2210523,1598697,271123,4115921,2267511,492097,35099,1849639,1537197
SNP_A-4301986	2	4	7,7,3,3	3,4,5,6	2	40	2600802,3425234,2175774,2342096,3450400,5869312,2783820,2603362,3427794,2178334,2344656,3447840,5866752,2786380,2979548,6330858,2864478,2982108,6333418,2867038,2600803,3425235,2175775,2342097,3450401,5869313,2783821,2603363,3427795,2178335,2344657,3447841,5866753,2786381,2979549,6330859,2864479,2982109,6333419,2867039
</pre>

There are three versions of the spf file format.  We refer
to them as "v1", "v2" and "v3".  The spf-v1 format was
described above.  The v2 and v3 formats are extensions of
spf to add allele and context codes to SPF.

While the formats can be identified by the column names, we
suggest putting in a meta-data header "#%spf-format=X" so
the readers can quickly check the format of the file.  If
there isnt an "spf-format" header the format is assumed to
be "1".

The v2 format has the two columns added to the v1 format:
"block_alleles" and "block_contexts".  Both of these columns
are comma seperated list of values of the same length as
num_blocks.  These values are used to populate the
allele_code and context code of each probelist block.  The
advantage of this format is it is compatable with the v1
format, the extra fields can be ignored if the sw doesnt
need them.

The v3 format restructures the file do away with the comma
seperated list of values.  By placing the blocks and probes
onto seperate lines grouped by indentation, the SPF file is
eaiser to modify.  Blocks and probeids can be commented out
without having to count the commas in a long list.

<pre>
#%spf-format=3
#%header0=name type
#%header1=allele allele_code [allele_bases]
#%header2=context context_code annotation
#%header3=probes
</pre>

The "allele" and "context" columns always have the value
"allele" and "context".  Their purpose is to act as a guide
for humans who are looking at the file.  The column
"allele_bases" are the bases that allele represents, it is
only for checking the conversion. The other fields are the
same as their comma seperated versions.

To look at some example data from the DMET SPF file:

<pre>
#%spf-format=3
#%header0=name type
#%header1=allele allele_code allele_bases
#%header2=context context_code annotation
#%header3=probes
DM3A10001	1
	allele	0	A
		context	0	0
			1145888
			173152
			277766
      [...more...]
		context	1	0
			1257548
			1259122
			312266
      [...more...]
		context	2	0
			1680942
			1682516
			176300
      [...more...]
	allele	1	C
		context	0	0
			1145887
			173151
			277765
      [...more...]
		context	1	0
			1257547
			1259121
			312265
      [...more...]
		context	2	0
			1680941
			1682515
			176299
      [...more...]
</pre>

The v3 format is oriented to DMET data, which has an allele
and context focus.  The v4 format of the SpfFile has a block
focus.  The fields related to each block appear after the
filler text "block".  (Which, again, is just filler text for
humans to read.)  This makes block-level hand edits eaiser
than the v3 format.

These columns have the same meaning as they do in the other
spf formats.  (The data is the same, it is just the
arrangement which has changed.) One item which is missing is
the probe count for each block; The count is implict in the
number of probes which follow.

The "spf-format","chip_type", "num-cols", "num-rows" headers
are required.  The "num-channels" is not required but
defaults to "1" if not present.

Here is an example v4 file:

<pre>
#%spf-format=4
#%chip_type=Axiom_GW_Hu_SNP
#%num-cols=1190
#%num-rows=1190
#%num-probesets=636352
#%num-channels=2
#%header0=name  type    num_match
#%header1=      block   annotation      allele_code     context_code    channel_code    rep_type
#%header2=              probe_id
AFFX-LCP-11715001       9       1
        block   0       0       0       0       3
                98058
                98454
                333283
                333876
                567911
                568504
AFFX-LCP-11715002       9       1
        block   0       0       0       1       3
                98060
                98456
                333285
                333878
                567913
                568506
</pre>

*/

/**
 * ProbeList class contains information about a group of probes that are going
 * to be analyzed together. Should have a 1-1 relationship with the ubiquitous
 * notion of a probe set in Affymetrix nomenclature. Class is essentially a list
 * of probe ids that should be processed together and the type of ProbeList
 * determines how that list should be viewed. There is a notion of a "blocks"
 * (aka groups from cdf file) and block annotations that should be analyzed
 * together (for example different alleles for genotyping are represented by
 * different blocks). If there are mismatches for perfect match probes the
 * numMatch field is set to '2' and the list is conceptually divided in half
 * with the PM probes in the first half and MM in the second half with the PM/MM
 * probes matches by position.  For example the first PM probe in the list
 * corresponds to the first MM probe found just after the halfway point in the
 * list. There is a notion of a 'null_probe' placeholder in case you have a
 * ProbeList of all MM probes (yes it happens) or a probeset where some PM
 * probes have matching MM and others don't. An additional class for probe level
 * data (i.e. probe info class) is planned for cases where specific probe level
 * information is needed for a particular application. This avoids cluttering up
 * ProbeList with information that will only be used in specialized
 * applications.
 *
 * A little history: Most of Affy's chips come with a cdf or pgf file that
 * document which probes should be analyzed together. To have a common datatype
 * to use across arrays the ProbeSet/Atom/Probe classes were born. ProbeSet and
 * Probe are pretty obvious and Atom is conceptually a PM/MM pair of probes that
 * should be analyzed together (or even more for a resequencing array, but those
 * were never supported). That worked fine for a while, but the nested nature of
 * the ProbeSet data type and the fact that newer designs don't have MM probes
 * makes it very inefficient in terms of RAM usage.
 *
 * This brings us to the ProbeList class which is very compact in memory
 * (especially with a custom heap for contiguous memory allocation). Here a
 * block of memory is allocated and there data is just offsets into the
 * array. This avoids any vectors and their wasted space. The design is also
 * amenable to memory mapping as everything is contiguous in memory rather than
 * a nested data structure.
 */

//////////

//////////

//////////

// ProbeList Block elements.
struct ProbeList_E_Block {
  int     m_blockSize;
  short   m_blockAnn;
  char    m_blockAllele;
  char    m_blockContext;
  char    m_blockChannel;
  char    m_blockRepType;
} PROBELIST_PACKED;

// ProbeList Probe elements
struct ProbeList_E_Probe {
  int    m_probeId;
  char   m_probeGc;
} PROBELIST_PACKED;

// The head of a ProbeList
// the sizes here allow us to compute the offsets of the Block and Probe elements.
struct ProbeList_Head {
  unsigned int m_type;
  unsigned int m_numMatch;
  // int   m_block_size; // we wont allow later additions. size==cnt
  // int   m_probe_size; // ditto for probes
  unsigned int   m_block_cnt;
  unsigned int   m_probe_cnt;
  /// the AnalysisProbeId (apid) of the 0th probe of this ProbeList
  /// Apids are assigned as the probes and blocks are added to the ProbeListFactory.
  /// Thus, the apids dont equal the probeIds; Even in the case of a single-channel design.
  unsigned int   m_apid_start;
  // char* name is right at the end of the elements
  // length of the name buffer (including the 0)
  unsigned short m_name_len;
} PROBELIST_PACKED;

// If ProbeListPacked were descended from ProbeList we would have
// :public ProbeList
// here.  It isnt as that would make it one pointer bigger.
// Its current size is one pointer.
class ProbeListPacked {
public:
  ProbeListPacked();
  ProbeListPacked(ProbeList_Head* ptr);

  bool isNull() const;

  // a pointer into a ProbeListFactoryRegion
  ProbeList_Head* m_headptr;

  int get_type() const;
  void set_type(int type);

  //
  int getApidStart() const;
/* inline int */
/* getApidStart() const */
/* { */
/*   assert(m_headptr!=NULL); */
/*   return PROBELIST_LOAD32B(&m_headptr->m_apid_start); */
/* } */
  void setApidStart(int start);
  //
  int findProbeApid(int qPid,int qAllele,int qContext,int qChannel) const;

  // the size of this ProbeList in bytes
  int byte_size() const;
  // compute the byte size needed.
  static int byte_size(int b_cnt,int p_cnt,int n_len);
  // accessors
  int block_cnt() const;
  int probe_cnt() const;
  int max_name_len() const;
  //
  int get_numMatch() const;
  void set_numMatch(int numMatch);
  //
  /// @todo add set_XXX(std::vector<int>);
  int get_blockSize(int bidx) const;
  void set_blockSize(int bidx,int val);
  //
  short get_blockAnn(int bidx) const;
  void set_blockAnn(int bidx,short val);
  //
  int get_blockAllele(int bidx) const;
  void set_blockAllele(int bidx,int val);
  //
  int get_blockContext(int bidx) const;
  void set_blockContext(int bidx,int val);
  //
  int get_blockChannel(int bidx) const;
  void set_blockChannel(int bidx,int val);
  //
  int get_blockRepType(int bidx) const;
  void set_blockRepType(int bidx,int val);

  // These are just like the ProbeListStl methods.
  int findBlockMatch(int qAllele,int qContext,int qChannel) const;
  bool hasDuplicateBlocks() const;

  //
  int get_probeIdxForBlock(int blockIdx) const;

  // the address of an element
  // const ProbeList_E_Block* const get_addr_E_block(unsigned int i) const;
  // const ProbeList_E_Probe* const get_addr_E_probe(unsigned int i) const;
  ProbeList_E_Block* get_addr_E_block(unsigned int i) const;
  ProbeList_E_Probe* get_addr_E_probe(unsigned int i) const;

/* inline ProbeList_E_Probe* */
/* get_addr_E_probe(unsigned int idx) const */
/* { */
/*   APT_ERR_ASSERT(m_headptr!=NULL, */
/*                  "Head pointer is NULL."); */
/*   APT_ERR_ASSERT(idx<probe_cnt(), */
/*                  "Probeset: '" + get_name_string() + "' - wrong number of probes." */
/*                  "(Max: " + ToStr(probe_cnt() - 1) + " Got: " + ToStr(idx)+")"); */

/*   char* ptr=(char*)m_headptr+sizeof(ProbeList_Head)+ */
/*     (block_cnt()*sizeof(ProbeList_E_Block))+ */
/*     (idx*sizeof(ProbeList_E_Probe)); */
/*   return (ProbeList_E_Probe*)ptr; */
/* } */

  //
  int get_probeId(int pidx) const;

/* inline int */
/* get_probeId(int i) const */
/* { */
/*   return PROBELIST_LOAD32B(&get_addr_E_probe(i)->m_probeId); */
/* } */
  void set_probeId(int pidx,int val);

  // Fill the vector with all the probeids at once.
  void get_probeIds(std::vector<int>& probeIds) const;

  // Fill the vector with all the probeids of a block at once.
  void get_probeIdsForBlock(int blockIdx,std::vector<int>& probeIds) const;

  //
  int pidxToBidx(int pidx) const;
  bool isPm(int pidx) const;

  /// Gets the Apid of this probe.
  /// See the note next to m_apid_start about probeId!=probeApid.
  // there is no "set_probeApid" because this value is "block.m_apid_start + i".
  int get_probeApid(int pidx) const;
/* inline int */
/* get_probeApid(int i) const */
/* { */
/*   if ((i>=0)&&(i<m_headptr->m_probe_cnt)) { */
/*     return getApidStart()+i; */
/*   } */
/*   return -1; */
/* } */

  //
  int8_t get_probeGc(int pidx) const;
/* inline int8_t */
/* get_probeGc(int i) const */
/* { */
/*   return PROBELIST_LOAD8B(&get_addr_E_probe(i)->m_probeGc); */
/* } */
  void set_probeGc(int pidx,int8_t val);

  // the address of the start of the name
  char* get_addr_E_name();
  const char* get_addr_E_name() const;
  const char* get_name_cstr() const;
  std::string get_name_string() const;
  void set_name(const std::string &name);

  //
  void dump() const;

  static double timecheck10;
  static double timecheck11;


};

// A region of allocation for ProbeListPacked.
// No need to pack the fields of a region as there arent that many of them,
// as compared to ProbeListPacked.
struct ProbeListFactory_Region {
  char* m_start_ptr;  ///< start of malloced region
  char* m_fill_ptr;   ///< where the next alloc will happen
  char* m_end_ptr;    ///< end of the memory region. (Dont go beyond here.)
};

/// A Factory for allocating and retrieving ProbeListPacked objects.
/// It really should be called ProbeListPackedFactory or PackedProbeListFactory.
class ProbeListFactory {
public:
  /// regions which this factory has available to fill.
  std::vector<ProbeListFactory_Region> m_region;
  /// current region in which allocations will occur.
  unsigned int m_ridx;
  /// how big should automaticly allocated regions be?
  int m_default_region_size;
  /// keep track of the number of allocations for usage info.
  unsigned int m_allocs;

  /// Pointers to the start of each Probelist in our regions.
  std::vector<ProbeListPacked> m_probelist_vec;

  // from the spf file
  std::vector<std::string> m_chipTypes;

  // data we remember to read from and put into an spf file.
  int m_numCols; // X
  int m_numRows; // Y
  int m_numChannels;

  int getXYsize() { return m_numCols * m_numRows; };

  /// these are mutable as we lazily sort when needed.
  /// A vector of indexes sorted by name.
  mutable std::vector<int> m_name2idx_vec;
  /// Is the vector currenly sorted? (So we dont have to sort each time.)
  mutable bool m_name2idx_issorted;
  mutable int m_name2idx_sortcnt;

  ///
  mutable std::vector<int>* m_probeid2apid_vecptr;

  /// the total number of probes in this factory.
  /// this is used to assign AnalysisProbeIds as probes are added.
  unsigned int m_probecnt_total;

  //
  ProbeListFactory();
  ~ProbeListFactory();

  //
  int getProbeSetCount() const;
  int getProbeSetIndexByName(const std::string& name) const;

  //
  void reserve(int size);
  void getDimensions(int cols,int rows) const;
  void setDimensions(int cols,int rows);
  //
  int numChannels() const;
  void setNumChannels(int num_channels);

private:
  void add_probelist_finish(ProbeListPacked& pl);

public:

  //
  void incApid(int cnt);
  unsigned int getApidMax() const;

  //
  // int getApidForProbeId(int probeid) const;

  // Cant overload with a string and an int, as  (0 might be a NULL)
  int findProbeApidByName(const std::string& qName,int qPid,int qAllele,int qContext,int qChannel) const;
  int findProbeApidByIdx (int qPlIdx              ,int qPid,int qAllele,int qContext,int qChannel) const;

  /// This only works when the probeids are unique. (Ie: single channel)
  void buildProbeId2ApidIndex() const ;
  int findProbeApidByProbeId(int qPid) const ;

  //
  void set_default_region_size(int size);
  // add a new region to the factory.
  void new_region();
  void new_region(int region_size);
  // free all our memory.
  void delete_regions();
  // clear the ProbeLists, but keep our memory for reuse.
  void clear();
  // return a pointer to a block of memory of this size
  char* alloc_ProbeList_bytes(int size);
  // return a ProbeList of the right size
  ProbeListPacked add_ProbeList(int block_cnt,int probe_cnt,int name_len);
  ProbeListPacked add_ProbeList(int block_cnt,int probe_cnt, const std::string &name);
  //
  ProbeListPacked add_ProbeList(const ProbeListStl& pl);

  // Does this probeset have MM or other non-PM type probes?
  static bool hasNonPm(const ProbeSet &ps);
  // How many PM probes are present in this probeset?
  static unsigned int countPmProbes(const ProbeSet *ps);
  // Insert a probe into an atom if it isn't the null probe.
  static bool maybeInsertProbe(Atom *a, int position, int id, unsigned char gcCount, unsigned char type,int apid);
  // Convert a ProbeSet into a ProbeList this memory is managed by the ProbeListFactory
  // @todo change back.
  //ProbeListPacked pListFromPSet(const ProbeSet& ps);
  void pListFromPSet(const ProbeSet& ps);
  // Convert a ProbeList into a ProbeSet.
  // The resulting ProbeSet needs to be deleted as it has been allocated by new()
  static ProbeSet* asProbeSet(const ProbeListPacked& pl);
  // silly function
  static double safe_div(double a,double b);

  // reading
  void readSpfFile(const std::string& spfName);
  void readSpfFile(affx::SpfFile& spfFile);

  // writing
  void writeSpfFile(const std::string& spfName,int spfFormat) const;
  void writeSpfFile(const std::string& spfName,int spfFormat,int with_allelecontext,int with_channel) const;
  //
  void writeToSpfFile(affx::SpfFile& spfFile) const;
  void writeToSpfFile_v2(affx::SpfFile& spfFile) const;
  void writeToSpfFile_v3(affx::SpfFile& spfFile) const;
  static void writeProbeListToSpfFile_v4(affx::SpfFile &spfFile, ProbeListPacked pl);
  static void writeProbeListToSpfFile_v3(affx::SpfFile& spfFile, ProbeListPacked pl);
  void writeToSpfFile_v4(affx::SpfFile& spfFile) const;

  //
  ProbeListPacked getProbeListAtIndex(int i) const;
  ProbeListPacked getProbeListByName(const std::string& name) const;
  //
  int getProbeListIndexByName(const std::string& name) const;

  //
  std::string getProbeListName(int i) const;
  const char* getProbeListNameCstr(int i) const;

private:
  ///
  void name2idx_dirty();
  void sortName2Idx() const;

public:
  /// check that the probelists are indeed sorted by name.
  /// shouldnt need to call this unless you think they arent sorted.
  void assertName2IdxSorted(int verbose) const;
  /// ensures that the map is sorted.
  void ensureSortedName2Idx() const;
  /// put the probelists in order by their name.
  /// NB: this will change the index order; Thus the results of getProbeListAtIndex will change.
  void sortFactoryByName();

  /// for debugging... display the contents of the Factory.
  void dump() const;
  /// for debugging...
  void dump_regions() const;
  void dump_probelist(int i) const;
  void dump_probelists() const;

  static double timecheck1;
  static double timecheck2;
  static double timecheck3;
  static double timecheck4;
  static double timecheck5;

};

#endif /* _PROBELISTFACTORY_H_ */
