# Uniprot-ported-gnomADs
These are datasets of gnomAD missense/nonsenses ported to Uniprot numbering.

## TL;DR

**AIM**: I want to copy-paste a snippet and get the data without downloading/pulling the repository.

```python
import requests
from collections import ChainMap
from typing import Dict, List
from typing_extensions import TypedDict  # 3.7 or lower, else just use typing

class GnomadDict(TypedDict):
    id: str
    consequence: str # 'missense_variant' etc.
    type: str # same as the above but missense or nonsense
    description: str
    from_residue: str  # single letter. May be more if weirdo fs
    residue_index: int  # human/fortran index
    to_residue: str
    impact: str # gnomAD verdict
    homozygous: int # N hmz
    frequency: float
    N: int # allele count
    # this is complicated:
    frequencies: Dict[str, float]
    # Re frequencies, https://gnomad.broadinstitute.org/help/ancestry has:
    # African/African American (afr),Amish (ami), American Admixed/Latino (amr), 
    # Ashkenazi Jewish (asj), East Asian (eas), Finnish (fin), Non-Finnish European (nfe) and South Asian (sas)
    # this is important because the human population is not global homogeneous and different ancestries
    # have different frequencies in different countries and admix at different rates not only between but also within.
    # i.e. consagnuity and community segregation mean that you cannot extrapolate from ``frequency``

def get_chunk(i) -> dict:
    url = f'https://github.com/matteoferla/Uniprot-ported-gnomADs/raw/main/uniprot_gnomADs/uniprot_gnomADs.{i}.json'
    response : requests.Response = requests.get(url)
    response.raise_for_status()
    return response.json()

# download the whole thing
uniprot_gnomads: Dict[str, List[GnomadDict]] = dict( **ChainMap(*map(get_chunk, range(0, 21))) )
```
NB: The data herein is a subset: it is only for protein coding regions and does not contain silent mutations.

They contain the uniprot numbered variants with details of frequency,
see the typed dictionary above for more details.

If a specific protein is sough and you cannot be bothered finding its uniprot id,
then you can do:

```python
import requests, json
from typing import Dict, List


def retrieve_data(filename: str) -> dict:
    url = f'https://github.com/matteoferla/Uniprot-ported-gnomADs/raw/main/{filename}'
    response : requests.Response = requests.get(url)
    response.raise_for_status()
    return response.json()

gene_name = '👾👾👾'
names2uniprot: Dict[str, str] = retrieve_data('uniprot_gnomADs/taxid9606-names2uniprot.json')
uniprot_id: str = names2uniprot[gene_name]
index: Dict[str, str] = retrieve_data('index.json')
filename: str = index[uniprot_id]
data: GnomadDict = retrieve_data(f'uniprot_gnomADs/{filename}')[uniprot_id]
```

## Mapping

> See [mapping](mapping.md) for Python code

I was told the match up between Uniprot and Ensembl is pretty good. But I never checked...

1,315,425 missense variants in gnomAD match the Uniprot sequence, 380,005 do not.

{'same': 108,519, 'good_change': 143696, 'bad_change': 11974, 'error': 30933}

A major problem I spotted is that some entries have cross-database links to canonical transcripts,
whereas in reality it's to other isoforms, so I mistrusted these.
* 7% of gnomAD canonical VEP annotations do not actually match the current Ensembl canonical.
* 81% match up between Uniprot and Ensembl —82% between Uniprot and gnomAD.
* 1% have some missense between gnomAD and Uniprot

3% had just some missenses between the sequences.
3% had differences in the sequences.
And less than 1% (44) of the Ensembl canonicals did not have a valid corresponding SwissProt Uniprot entry 
—SwissProt is curated, Trembl has all the "novel protein".
The gnomAD canonicals have 371 (see Table 1).

Using the gnomAD canonical, 1968 protein can be corrected based on the gnomAD canonical
(108,519 missenses with unaltered index, 143,690 with a concordant starting residue),
739 has missenses whose positition was the same,
while 140 were problematic and could not corrected (Table 2).

## Rant
I am not a geneticist. I am a biochemist. I do not really care about transcripts etc.
but a big problem is that Ensembl and Uniprot are not the same in many cases when it comes to what is canonical,
and, to boot, geneticists, will report not the canonical, but the worst consequence in VEP, which is generally
the first in the list, i.e. sorted by ENST, which is meaningless.
In Ensembl, the canonical transcript is the longest, in Uniprot, the most common transcriptomically.
GTex can be used to check this.

SWIFT and AlphaFold run off Uniprot.
Uniprot is curated and awesome for protein.
gnomAD is really cool as it tells you want protein missenses are present in the human population,
which for the control subset is generally healthy (or healthy carriers). It runs off Ensembl. 

## Table 1 (unmatched)
There are better ways to match ENSP to Uniprot, but for now by symbol is fine,
as only 371 are lost.

| ESNP            |   ESNP_len | symbol          |
|:----------------|-----------:|:----------------|
| ENSP00000432016 |       2617 | ANKHD1-EIF4EBP3 |
| ENSP00000016946 |       1765 | RGPD5           |
| ENSP00000415941 |       1744 | C4B             |
| ENSP00000491742 |       1640 | BIVM-ERCC5      |
| ENSP00000355890 |       1512 | EPRS1           |
| ENSP00000364815 |       1264 | VARS1           |
| ENSP00000364794 |       1262 | IARS1           |
| ENSP00000264229 |       1233 | CRACD           |
| ENSP00000378236 |       1182 | STON1-GTF2A1L   |
| ENSP00000377954 |       1176 | LARS1           |
| ENSP00000363654 |       1103 | PALM2AKAP2      |
| ENSP00000299505 |       1070 | GARRE1          |
| ENSP00000413445 |       1029 | ELAPOR2         |
| ENSP00000491516 |       1027 | SCART1          |
| ENSP00000358955 |       1013 | ELAPOR1         |
| ENSP00000371886 |       1012 | JMJD7-PLA2G4B   |
| ENSP00000502797 |        968 | AARS1           |
| ENSP00000380996 |        962 | CRACDL          |
| ENSP00000412130 |        912 | IRAG1           |
| ENSP00000274140 |        910 | MARCHF6         |
| ENSP00000262027 |        900 | MARS1           |
| ENSP00000463080 |        846 | MARCHF10        |
| ENSP00000346916 |        842 | TRIM6-TRIM34    |
| ENSP00000242786 |        835 | ADGRE5          |
| ENSP00000369897 |        831 | CARS1           |
| ENSP00000338093 |        802 | TARS3           |
| ENSP00000307567 |        775 | QARS1           |
| ENSP00000419005 |        756 | UTP25           |
| ENSP00000383996 |        752 | CERT1           |
| ENSP00000400168 |        749 | ATP5MF-PTCD1    |
| ENSP00000441269 |        748 | SLCO1B3-SLCO1B7 |
| ENSP00000490573 |        712 | TCAF2C          |
| ENSP00000386830 |        704 | MARCHF7         |
| ENSP00000231572 |        660 | RARS1           |
| ENSP00000313500 |        647 | HROB            |
| ENSP00000491279 |        638 | KDM4F           |
| ENSP00000349458 |        632 | CILK1           |
| ENSP00000380427 |        625 | ARPC4-TTLL3     |
| ENSP00000500220 |        614 | ARSL            |
| ENSP00000433415 |        602 | NT5C1B-RDH14    |
| ENSP00000411848 |        573 | MARCHF8         |
| ENSP00000354732 |        564 | PRR5-ARHGAP8    |
| ENSP00000295966 |        563 | CFAP20DC        |
| ENSP00000489685 |        560 | IQANK1          |
| ENSP00000358939 |        536 | SARS1           |
| ENSP00000480571 |        535 | CYP3A7-CYP3A51P |
| ENSP00000362576 |        528 | YARS1           |
| ENSP00000353355 |        525 | PTGES3L-AARSD1  |
| ENSP00000477920 |        511 | GIMAP1-GIMAP5   |
| ENSP00000330484 |        511 | AMY1B           |
| ENSP00000359100 |        511 | AMY1A           |
| ENSP00000425634 |        509 | HARS1           |
| ENSP00000421868 |        502 | CCDC169-SOHLH2  |
| ENSP00000264161 |        501 | DARS1           |
| ENSP00000346442 |        499 | IRAG2           |
| ENSP00000489658 |        492 | CNTNAP3C        |
| ENSP00000485258 |        478 | PRAMEF9         |
| ENSP00000347495 |        471 | WARS1           |
| ENSP00000447121 |        471 | CHURC1-FNTB     |
| ENSP00000473391 |        469 | COMMD3-BMI1     |
| ENSP00000317579 |        468 | RUSF1           |
| ENSP00000466174 |        463 | ATF7-NPFF       |
| ENSP00000474850 |        452 | TRIM49D1        |
| ENSP00000299957 |        436 | GOLM2           |
| ENSP00000359444 |        423 | HSFX1           |
| ENSP00000323036 |        418 | ACP3            |
| ENSP00000413165 |        417 | CKMT1A          |
| ENSP00000273067 |        410 | MARCHF4         |
| ENSP00000265068 |        409 | CFAP92          |
| ENSP00000488257 |        407 | MAGEB6B         |
| ENSP00000282382 |        404 | TMED7-TICAM2    |
| ENSP00000363984 |        403 | PGAP4           |
| ENSP00000333181 |        402 | MARCHF11        |
| ENSP00000498431 |        402 | SPDYE21         |
| ENSP00000355812 |        399 | CEP43           |
| ENSP00000422907 |        395 | GTF2H2C         |
| ENSP00000293662 |        395 | TAMALIN         |
| ENSP00000498596 |        379 | STING1          |
| ENSP00000267103 |        376 | MYG1            |
| ENSP00000389275 |        369 | NPIPA9          |
| ENSP00000362273 |        363 | ADPRS           |
| ENSP00000262644 |        359 | BPNT2           |
| ENSP00000422067 |        351 | FAM47E-STBD1    |
| ENSP00000466379 |        348 | GET3            |
| ENSP00000319799 |        346 | H1-8            |
| ENSP00000266643 |        346 | MARCHF9         |
| ENSP00000355877 |        337 | MTARC1          |
| ENSP00000355880 |        335 | MTARC2          |
| ENSP00000411822 |        331 | ISY1-RAB43      |
| ENSP00000293826 |        330 | TNFSF12-TNFSF13 |
| ENSP00000429150 |        324 | CYRIB           |
| ENSP00000370724 |        323 | CYRIA           |
| ENSP00000493820 |        323 | OR5D3P          |
| ENSP00000445323 |        321 | ZNF559-ZNF177   |
| ENSP00000357542 |        318 | ZNF511-PRAP1    |
| ENSP00000492889 |        314 | OR1R1P          |
| ENSP00000357731 |        312 | LORICRIN        |
| ENSP00000497734 |        312 | OR51C1P         |
| ENSP00000494254 |        311 | OR5BS1P         |
| ENSP00000493454 |        310 | OR2A1           |
| ENSP00000493452 |        308 | OR9H1P          |
| ENSP00000396755 |        307 | CCNP            |
| ENSP00000443411 |        304 | CIBAR2          |
| ENSP00000370083 |        294 | SMN1            |
| ENSP00000427223 |        289 | MARCHF1         |
| ENSP00000429367 |        289 | CIBAR1          |
| ENSP00000499970 |        289 | NFILZ           |
| ENSP00000351813 |        278 | MARCHF5         |
| ENSP00000251303 |        275 | SLX1A           |
| ENSP00000485259 |        274 | NOTCH2NLR       |
| ENSP00000497977 |        273 | GET1-SH3BGR     |
| ENSP00000498945 |        265 | SPDYE9          |
| ENSP00000334805 |        255 | H1-7            |
| ENSP00000335808 |        254 | PSME3IP1        |
| ENSP00000309141 |        253 | MARCHF3         |
| ENSP00000471536 |        246 | MARCHF2         |
| ENSP00000385560 |        243 | RBAK-RBAKDN     |
| ENSP00000411198 |        235 | KRTAP4-16       |
| ENSP00000492766 |        233 | FOXL3           |
| ENSP00000375656 |        233 | ZNF816-ZNF321P  |
| ENSP00000450085 |        233 | RNASEK-C17orf49 |
| ENSP00000463379 |        228 | RPL17-C18orf32  |
| ENSP00000403937 |        221 | IFTAP           |
| ENSP00000458075 |        217 | BOLA2-SMG1P6    |
| ENSP00000441930 |        214 | DNAAF6          |
| ENSP00000329662 |        213 | H1-10           |
| ENSP00000485552 |        208 | SAA2-SAA4       |
| ENSP00000341214 |        207 | H1-6            |
| ENSP00000326110 |        206 | MACIR           |
| ENSP00000475025 |        201 | TLCD4-RWDD3     |
| ENSP00000362618 |        200 | PABPC1L2A       |
| ENSP00000491873 |        198 | TAF11L10        |
| ENSP00000492874 |        198 | TAF11L5         |
| ENSP00000492752 |        198 | LINC02218       |
| ENSP00000492047 |        198 | TAF11L2         |
| ENSP00000492463 |        198 | TAF11L9         |
| ENSP00000491174 |        198 | TAF11L8         |
| ENSP00000492787 |        198 | TAF11L3         |
| ENSP00000491846 |        198 | TAF11L11        |
| ENSP00000461830 |        198 | CEP20           |
| ENSP00000491332 |        198 | TAF11L13        |
| ENSP00000491494 |        198 | TAF11L4         |
| ENSP00000491758 |        197 | TAF11L12        |
| ENSP00000491472 |        197 | TAF11L14        |
| ENSP00000491663 |        194 | SSU72P4         |
| ENSP00000492512 |        194 | SSU72P8         |
| ENSP00000491894 |        194 | SSU72P7         |
| ENSP00000491768 |        194 | SSU72P3         |
| ENSP00000491949 |        194 | SSU72P5         |
| ENSP00000491669 |        194 | SSU72P2         |
| ENSP00000494697 |        192 | RPL9            |
| ENSP00000471224 |        192 | DMRTC1          |
| ENSP00000276927 |        189 | IFNA1           |
| ENSP00000469011 |        188 | SSX4            |
| ENSP00000479168 |        187 | PRH1            |
| ENSP00000424176 |        179 | EPPIN-WFDC6     |
| ENSP00000496813 |        174 | GET1            |
| ENSP00000301408 |        165 | CGB5            |
| ENSP00000349954 |        165 | CGB3            |
| ENSP00000429865 |        164 | TVP23C-CDRT4    |
| ENSP00000359434 |        158 | EOLA2           |
| ENSP00000454021 |        155 | TGIF2-RAB5IF    |
| ENSP00000501406 |        153 | TUG1            |
| ENSP00000497053 |        147 | H3Y2            |
| ENSP00000354723 |        147 | H2BW1           |
| ENSP00000322421 |        142 | HBA1            |
| ENSP00000333277 |        136 | H3C13           |
| ENSP00000367891 |        136 | LGALS7          |
| ENSP00000496014 |        136 | H3Y1            |
| ENSP00000484638 |        136 | H3C8            |
| ENSP00000355780 |        136 | H3-3A           |
| ENSP00000366999 |        136 | H3C4            |
| ENSP00000480826 |        136 | H3C1            |
| ENSP00000484841 |        136 | H3C2            |
| ENSP00000484658 |        136 | H3C3            |
| ENSP00000484095 |        136 | H3C7            |
| ENSP00000483283 |        136 | H3C11           |
| ENSP00000358160 |        136 | H3C10           |
| ENSP00000489282 |        136 | H3C6            |
| ENSP00000355657 |        136 | H3-4            |
| ENSP00000339835 |        135 | H3-5            |
| ENSP00000445831 |        134 | H2BC18          |
| ENSP00000297012 |        131 | H2AC1           |
| ENSP00000351589 |        130 | H2AC13          |
| ENSP00000483842 |        130 | H2AC4           |
| ENSP00000352119 |        130 | H2AC11          |
| ENSP00000482538 |        130 | H2AC16          |
| ENSP00000367022 |        130 | H2AC6           |
| ENSP00000482431 |        130 | H2AC15          |
| ENSP00000355656 |        130 | H2AW            |
| ENSP00000332790 |        130 | H2AC21          |
| ENSP00000368747 |        130 | TRBV20OR9-2     |
| ENSP00000438553 |        129 | H2AJ            |
| ENSP00000332194 |        129 | H2AC20          |
| ENSP00000308405 |        128 | H2AZ2           |
| ENSP00000366679 |        128 | H2AC12          |
| ENSP00000491733 |        127 | CSAG2           |
| ENSP00000289316 |        126 | H2BC5           |
| ENSP00000479284 |        126 | H2BU1           |
| ENSP00000489317 |        126 | H2BC6           |
| ENSP00000445633 |        126 | H2BC8           |
| ENSP00000477527 |        126 | H2BC17          |
| ENSP00000321744 |        126 | H2BC4           |
| ENSP00000479169 |        126 | H2BC9           |
| ENSP00000348924 |        126 | H2BC7           |
| ENSP00000358151 |        126 | H2BC21          |
| ENSP00000482674 |        126 | H2BC3           |
| ENSP00000455596 |        124 | TP53TG3D        |
| ENSP00000496189 |        120 | TINCR           |
| ENSP00000474509 |        119 | IGHV3OR15-7     |
| ENSP00000474384 |        119 | IGHV2OR16-5     |
| ENSP00000386655 |        118 | RPL36A-HNRNPH2  |
| ENSP00000498060 |        118 | RAMACL          |
| ENSP00000474814 |        117 | IGHV3OR16-13    |
| ENSP00000474013 |        117 | IGHV3OR16-12    |
| ENSP00000474639 |        117 | IGHV1OR15-9     |
| ENSP00000481738 |        117 | IGHV1OR21-1     |
| ENSP00000374871 |        117 | TRGV1           |
| ENSP00000474395 |        117 | IGKV1OR2-108    |
| ENSP00000374864 |        117 | TRGV10          |
| ENSP00000339511 |        117 | H2AP            |
| ENSP00000474297 |        116 | IGKV3OR2-268    |
| ENSP00000489180 |        116 | ERVK3-1         |
| ENSP00000474941 |        116 | IGHV3OR16-10    |
| ENSP00000474964 |        116 | IGHV3OR16-8     |
| ENSP00000448600 |        115 | TRBV7-1         |
| ENSP00000374919 |        115 | TRBV23-1        |
| ENSP00000374884 |        115 | TRBV7-3         |
| ENSP00000374901 |        114 | TRBV5-7         |
| ENSP00000483468 |        114 | TRBV17          |
| ENSP00000374885 |        114 | TRBV5-3         |
| ENSP00000374896 |        114 | TRBV6-7         |
| ENSP00000482333 |        114 | TRBV15          |
| ENSP00000450448 |        113 | TRAV8-7         |
| ENSP00000473871 |        110 | URGCP-MRPS24    |
| ENSP00000346892 |        110 | SERF1A          |
| ENSP00000359693 |        108 | FKBP1C          |
| ENSP00000367034 |        103 | H4C3            |
| ENSP00000366956 |        103 | H4C8            |
| ENSP00000347168 |        103 | H4C11           |
| ENSP00000480960 |        103 | H4C13           |
| ENSP00000479794 |        103 | H4C12           |
| ENSP00000244537 |        103 | H4C6            |
| ENSP00000479461 |        103 | H4C4            |
| ENSP00000374863 |        103 | TRGV11          |
| ENSP00000479106 |        103 | H4C1            |
| ENSP00000366974 |        103 | H4C2            |
| ENSP00000484789 |        103 | H4C5            |
| ENSP00000481486 |        103 | H4C9            |
| ENSP00000477870 |         98 | H4C7            |
| ENSP00000474363 |         98 | IGHV3OR16-9     |
| ENSP00000347933 |         96 | PLGLB1          |
| ENSP00000490566 |         79 | TMEM238L        |
| ENSP00000334330 |         78 | DEFB105A        |
| ENSP00000492284 |         76 | KANTR           |
| ENSP00000320813 |         72 | DEFB104A        |
| ENSP00000334681 |         70 | DEFB107A        |
| ENSP00000335307 |         65 | DEFB106A        |
| ENSP00000303532 |         64 | DEFB4A          |
| ENSP00000419781 |         50 | IGLJ3           |
| ENSP00000489624 |         48 | TUNAR           |
| ENSP00000418274 |         47 | IGLJ2           |
| ENSP00000499205 |         33 | BLACAT1         |
| ENSP00000490930 |         25 | LINC00672       |
| ENSP00000418818 |         24 | IGLJ5           |
| ENSP00000419267 |         23 | IGLJ6           |
| ENSP00000451143 |         23 | TRAJ52          |
| ENSP00000451312 |         22 | TRAJ45          |
| ENSP00000451561 |         22 | TRAJ28          |
| ENSP00000450593 |         22 | TRAJ2           |
| ENSP00000451299 |         22 | TRAJ53          |
| ENSP00000451503 |         22 | TRAJ32          |
| ENSP00000483197 |         22 | TRAJ37          |
| ENSP00000451639 |         22 | TRAJ18          |
| ENSP00000452134 |         21 | TRAJ24          |
| ENSP00000450455 |         21 | TRAJ23          |
| ENSP00000452510 |         21 | TRAJ39          |
| ENSP00000451066 |         21 | TRAJ10          |
| ENSP00000451965 |         21 | TRAJ6           |
| ENSP00000452097 |         21 | TRAJ4           |
| ENSP00000452149 |         21 | TRAJ41          |
| ENSP00000452565 |         21 | TRAJ1           |
| ENSP00000451248 |         21 | TRAJ38          |
| ENSP00000452484 |         21 | TRAJ13          |
| ENSP00000451250 |         21 | TRAJ17          |
| ENSP00000451843 |         21 | TRAJ22          |
| ENSP00000451870 |         21 | TRAJ56          |
| ENSP00000450777 |         21 | TRAJ48          |
| ENSP00000452150 |         21 | TRAJ58          |
| ENSP00000452248 |         21 | TRAJ57          |
| ENSP00000451191 |         21 | TRAJ46          |
| ENSP00000451632 |         21 | TRAJ44          |
| ENSP00000418697 |         20 | TRGJP2          |
| ENSP00000450953 |         20 | TRAJ26          |
| ENSP00000451766 |         20 | TRAJ9           |
| ENSP00000419223 |         20 | IGHJ6           |
| ENSP00000450961 |         20 | TRAJ5           |
| ENSP00000452114 |         20 | TRAJ40          |
| ENSP00000452464 |         20 | TRAJ61          |
| ENSP00000451806 |         20 | TRAJ11          |
| ENSP00000450522 |         20 | TRAJ12          |
| ENSP00000450622 |         20 | TRAJ16          |
| ENSP00000451589 |         20 | TRAJ50          |
| ENSP00000451728 |         20 | TRAJ19          |
| ENSP00000451140 |         20 | TRDJ3           |
| ENSP00000452115 |         20 | TRAJ25          |
| ENSP00000452564 |         20 | TRAJ54          |
| ENSP00000451221 |         20 | TRAJ27          |
| ENSP00000452518 |         20 | TRAJ29          |
| ENSP00000418089 |         20 | TRGJP           |
| ENSP00000478270 |         20 | TRAJ36          |
| ENSP00000451687 |         20 | TRAJ35          |
| ENSP00000419604 |         20 | TRGJP1          |
| ENSP00000452600 |         20 | TRAJ7           |
| ENSP00000451847 |         19 | TRAJ20          |
| ENSP00000451671 |         19 | TRAJ30          |
| ENSP00000452141 |         19 | TRAJ31          |
| ENSP00000451657 |         19 | TRAJ49          |
| ENSP00000450793 |         19 | TRAJ34          |
| ENSP00000451175 |         19 | TRAJ33          |
| ENSP00000451749 |         19 | TRAJ47          |
| ENSP00000450760 |         18 | TRAJ21          |
| ENSP00000450690 |         18 | TRAJ43          |
| ENSP00000450875 |         18 | TRDJ2           |
| ENSP00000451065 |         18 | TRAJ59          |
| ENSP00000417149 |         17 | TRGJ2           |
| ENSP00000419748 |         17 | IGHJ2           |
| ENSP00000450438 |         17 | TRAJ14          |
| ENSP00000420517 |         16 | IGHJ3           |
| ENSP00000452040 |         16 | TRDJ4           |
| ENSP00000418164 |         16 | IGHJ5           |
| ENSP00000419233 |         15 | IGLJ7           |
| ENSP00000420142 |         15 | TRBJ2-2P        |
| ENSP00000419074 |         15 | IGHJ4           |
| ENSP00000418768 |         13 | IGKJ2           |
| ENSP00000420166 |         13 | IGKJ5           |
| ENSP00000419518 |         13 | IGKJ3           |
| ENSP00000428366 |         12 | IGHD3-16        |
| ENSP00000417097 |         12 | IGKJ4           |
| ENSP00000419997 |         11 | IGLJ4           |
| ENSP00000420442 |         10 | IGHD3-3         |
| ENSP00000474065 |         10 | IGHD2OR15-2A    |
| ENSP00000429952 |         10 | IGHD3-22        |
| ENSP00000427969 |         10 | IGHD2-15        |
| ENSP00000428616 |         10 | IGHD2-8         |
| ENSP00000474133 |         10 | IGHD3OR15-3A    |
| ENSP00000430788 |         10 | IGHD2-2         |
| ENSP00000419583 |          9 | IGHD3-9         |
| ENSP00000429324 |          9 | IGHD2-21        |
| ENSP00000419773 |          9 | IGHD3-10        |
| ENSP00000419564 |          7 | IGHD6-13        |
| ENSP00000419283 |          7 | IGHD5-12        |
| ENSP00000473700 |          7 | IGHD5OR15-5B    |
| ENSP00000418010 |          7 | IGHD6-19        |
| ENSP00000473849 |          7 | IGHD5OR15-5A    |
| ENSP00000417751 |          6 | IGHD6-25        |
| ENSP00000419139 |          6 | IGHD5-24        |
| ENSP00000417892 |          6 | IGHD5-5         |
| ENSP00000430248 |          6 | IGHD4-23        |
| ENSP00000417555 |          6 | IGHD5-18        |
| ENSP00000418151 |          6 | IGHD6-6         |
| ENSP00000420556 |          5 | IGHD1-20        |
| ENSP00000431089 |          5 | IGHD4-17        |
| ENSP00000418765 |          5 | IGHD1-14        |
| ENSP00000474222 |          5 | IGHD1OR15-1B    |
| ENSP00000430034 |          5 | IGHD4-11        |
| ENSP00000420794 |          5 | IGHD1-7         |
| ENSP00000428393 |          5 | IGHD4-4         |
| ENSP00000452494 |          4 | TRDD3           |
| ENSP00000418639 |          3 | IGHD7-27        |
| ENSP00000451515 |          3 | TRDD2           |

## Table 2
Sequences whose alignment prevented a conversion.

|            | symbol    | gnomAD_ENSP     |   bad_change |   good_change |   same |
|:-----------|:----------|:----------------|-------------:|--------------:|-------:|
| Q6V0I7     | FAT4      | ENSP00000501473 |          471 |            31 |      0 |
| Q96JP2     | MYO15B    | ENSP00000495242 |          332 |            18 |      0 |
| Q5SNV9     | C1orf167  | ENSP00000414909 |          309 |            29 |      0 |
| Q9Y2K3     | MYH15     | ENSP00000273353 |          283 |            24 |      0 |
| P52746     | ZNF142    | ENSP00000398798 |          262 |            13 |     17 |
| Q96AX9     | MIB2      | ENSP00000426103 |          246 |            20 |      0 |
| Q07157     | TJP1      | ENSP00000348416 |          231 |            22 |      0 |
| Q14005     | IL16      | ENSP00000302935 |          222 |            20 |      0 |
| Q9UKZ4     | TENM1     | ENSP00000403954 |          185 |             9 |      0 |
| Q6ZSZ5     | ARHGEF18  | ENSP00000482647 |          180 |            11 |      0 |
| Q5JV73     | FRMPD3    | ENSP00000276185 |          179 |            13 |      0 |
| Q92538     | GBF1      | ENSP00000359000 |          179 |            18 |     34 |
| Q12769     | NUP160    | ENSP00000367721 |          177 |            17 |      0 |
| Q5T6F2     | UBAP2     | ENSP00000354039 |          155 |            20 |     31 |
| Q9UPX8     | SHANK2    | ENSP00000345193 |          151 |            18 |      0 |
| Q08379     | GOLGA2    | ENSP00000478799 |          151 |             5 |     10 |
| P23109     | AMPD1     | ENSP00000430075 |          150 |             7 |      0 |
| O95785     | WIZ       | ENSP00000263381 |          143 |           132 |     12 |
| Q9H8H2     | DDX31     | ENSP00000361232 |          143 |             2 |      0 |
| Q6ZP01     | RBM44     | ENSP00000321179 |          137 |            10 |      0 |
| Q9ULE0     | WWC3      | ENSP00000370242 |          136 |             9 |      0 |
| Q96T25     | ZIC5      | ENSP00000267294 |          121 |            15 |      0 |
| Q8NHH1     | TTLL11    | ENSP00000321346 |          121 |             7 |      0 |
| P46379     | BAG6      | ENSP00000365131 |          115 |             5 |     20 |
| H3BPF8     | GOLGA8S   | ENSP00000455298 |          112 |            10 |      0 |
| Q8TE49     | OTUD7A    | ENSP00000305926 |          110 |            19 |     31 |
| Q6ZW05     | PTCHD4    | ENSP00000341914 |          109 |            11 |      0 |
| O14772     | FPGT      | ENSP00000359935 |          107 |             8 |      0 |
| Q5VTQ0     | TTC39B    | ENSP00000422496 |          106 |            10 |      0 |
| P55157     | MTTP      | ENSP00000427679 |          104 |             7 |      0 |
| P29973     | CNGA1     | ENSP00000426862 |          104 |            16 |      0 |
| Q5T0N1     | CFAP70    | ENSP00000310829 |          101 |             3 |     11 |
| Q03518     | TAP1      | ENSP00000346206 |          101 |            15 |      0 |
| O76011     | KRT34     | ENSP00000377570 |          101 |             6 |      0 |
| A8MY62     | LACTBL1   | ENSP00000402297 |           97 |             7 |      0 |
| P22105     | TNXB      | ENSP00000407685 |           95 |            89 |      0 |
| Q9BXY5     | CAPS2     | ENSP00000386959 |           94 |             8 |      0 |
| Q8TDS5     | OXER1     | ENSP00000367930 |           93 |             3 |      0 |
| Q12894     | IFRD2     | ENSP00000402849 |           90 |            10 |      0 |
| Q6NUS6     | TCTN3     | ENSP00000483364 |           90 |             2 |      0 |
| P01019     | AGT       | ENSP00000355627 |           89 |             8 |      0 |
| Q5SQN1     | SNAP47    | ENSP00000314157 |           87 |             3 |      0 |
| P50747     | HLCS      | ENSP00000502087 |           87 |             8 |      0 |
| P59826     | BPIFB3    | ENSP00000364643 |           86 |             6 |      0 |
| Q96FT7     | ASIC4     | ENSP00000350786 |           86 |             6 |      0 |
| Q9Y5X5     | NPFFR2    | ENSP00000307822 |           85 |             2 |      0 |
| P49761     | CLK3      | ENSP00000378505 |           83 |             2 |      0 |
| Q9UKU6     | TRHDE     | ENSP00000261180 |           83 |            12 |      0 |
| Q8NH41     | OR4K15    | ENSP00000304077 |           82 |             3 |      0 |
| Q14241     | ELOA      | ENSP00000395574 |           80 |             6 |      0 |
| Q99705     | MCHR1     | ENSP00000249016 |           79 |            10 |      0 |
| Q9UF12     | PRODH2    | ENSP00000301175 |           78 |             8 |      0 |
| Q01433     | AMPD2     | ENSP00000499465 |           77 |             6 |      0 |
| Q7Z2D5     | PLPPR4    | ENSP00000359204 |           77 |             5 |      0 |
| P52736     | ZNF133    | ENSP00000439427 |           75 |             8 |      0 |
| Q6ZUS6     | CCDC149   | ENSP00000488929 |           70 |             4 |      0 |
| Q15270     | NKX1-1    | ENSP00000407978 |           70 |             6 |      0 |
| Q6IFG1     | OR52E8    | ENSP00000444054 |           69 |            10 |      0 |
| Q6ZQW0     | IDO2      | ENSP00000443432 |           68 |             6 |      0 |
| P03999     | OPN1SW    | ENSP00000249389 |           68 |             2 |      0 |
| Q6NUI2     | GPAT2     | ENSP00000389395 |           67 |             2 |     33 |
| Q0D2K0     | NIPAL4    | ENSP00000311687 |           67 |             7 |      0 |
| Q8TE04     | PANK1     | ENSP00000302108 |           66 |             3 |      0 |
| Q14678     | KANK1     | ENSP00000477725 |           65 |            10 |    297 |
| Q8NGJ2     | OR52H1    | ENSP00000493308 |           65 |             4 |      0 |
| Q9HBY8     | SGK2      | ENSP00000340608 |           65 |             2 |      0 |
| Q8NGV7     | OR5H2     | ENSP00000347418 |           65 |             6 |      0 |
| Q6P587     | FAHD1     | ENSP00000372112 |           61 |             2 |      0 |
| Q66PJ3     | ARL6IP4   | ENSP00000313422 |           61 |             3 |      0 |
| Q96KE9     | BTBD6     | ENSP00000376337 |           61 |             8 |      0 |
| Q8WUG5     | SLC22A17  | ENSP00000380437 |           60 |             2 |      0 |
| Q86VX9     | MON1A     | ENSP00000296473 |           57 |             8 |      0 |
| Q8WYJ6     | SEPTIN1   | ENSP00000324511 |           56 |             2 |      0 |
| Q12983     | BNIP3     | ENSP00000357625 |           55 |             2 |      0 |
| Q14643     | ITPR1     | ENSP00000306253 |           55 |             4 |     63 |
| O75319     | DUSP11    | ENSP00000272444 |           55 |             7 |      0 |
| Q9H2S5     | RNF39     | ENSP00000244360 |           54 |             6 |      0 |
| P11362     | FGFR1     | ENSP00000393312 |           54 |             0 |     33 |
| Q9UJX3     | ANAPC7    | ENSP00000394394 |           53 |             2 |      0 |
| Q2TAM9     | TUSC1     | ENSP00000350716 |           52 |             9 |      0 |
| Q9NZV7     | ZIM2      | ENSP00000486502 |           50 |             4 |      1 |
| P0DME0     | SETSIP    | ENSP00000480946 |           50 |             6 |      0 |
| Q8NFV4     | ABHD11    | ENSP00000222800 |           50 |             6 |      0 |
| Q9HD15     | SRA1      | ENSP00000337513 |           50 |             9 |      0 |
| Q53RT3     | ASPRV1    | ENSP00000315383 |           49 |             4 |      0 |
| Q5MJ68     | SPDYC     | ENSP00000366390 |           48 |             8 |      0 |
| Q99871     | HAUS7     | ENSP00000359230 |           48 |             2 |      0 |
| Q9GZU5     | NYX       | ENSP00000340328 |           48 |             6 |      0 |
| Q8TE02     | ELP5      | ENSP00000379869 |           47 |             1 |      0 |
| Q9HBW0     | LPAR2     | ENSP00000443256 |           46 |             5 |      0 |
| Q00532     | CDKL1     | ENSP00000379176 |           45 |             3 |      0 |
| Q12904     | AIMP1     | ENSP00000378191 |           43 |             2 |      0 |
| Q8N5Z5     | KCTD17    | ENSP00000385096 |           43 |             4 |      0 |
| A0A0U1RRI6 | CENPVL3   | ENSP00000489547 |           42 |             0 |      0 |
| Q96DA0     | ZG16B     | ENSP00000371715 |           42 |             0 |      0 |
| A0PJE2     | DHRS12    | ENSP00000411565 |           42 |             3 |      0 |
| Q6URK8     | TEPP      | ENSP00000290871 |           41 |             6 |      0 |
| Q6SJ96     | TBPL2     | ENSP00000247219 |           40 |             4 |      0 |
| P0C1H6     | H2BW2     | ENSP00000347119 |           40 |             5 |      0 |
| Q9Y233     | PDE10A    | ENSP00000355847 |           38 |             2 |      0 |
| P49842     | STK19     | ENSP00000364482 |           38 |             4 |      0 |
| O95755     | RAB36     | ENSP00000263116 |           38 |             1 |      0 |
| O00142     | TK2       | ENSP00000299697 |           37 |             3 |      0 |
| Q9HDB5     | NRXN3     | ENSP00000489551 |           36 |             3 |      0 |
| A6NGB0     | TMEM191C  | ENSP00000489706 |           36 |             1 |      0 |
| Q6PFW1     | PPIP5K1   | ENSP00000400887 |           36 |             3 |      5 |
| O14893     | GEMIN2    | ENSP00000308533 |           35 |             2 |      0 |
| Q9H7D7     | WDR26     | ENSP00000498603 |           35 |             3 |      0 |
| B4E2M5     | ANKRD66   | ENSP00000454770 |           35 |             2 |      0 |
| Q6P4F2     | FDX2      | ENSP00000377311 |           33 |             4 |      0 |
| P0DO97     | CCDC192   | ENSP00000490579 |           33 |             2 |      0 |
| Q99766     | DMAC2L    | ENSP00000451583 |           33 |             3 |      0 |
| Q96HE8     | TMEM80    | ENSP00000380646 |           33 |             2 |      0 |
| Q14203     | DCTN1     | ENSP00000354791 |           33 |             5 |    148 |
| O43581     | SYT7      | ENSP00000444201 |           31 |            20 |     11 |
| Q96GE6     | CALML4    | ENSP00000419081 |           30 |             2 |      0 |
| Q8N8F6     | YIPF7     | ENSP00000332772 |           27 |             2 |      0 |
| Q9UNG2     | TNFSF18   | ENSP00000385470 |           27 |             2 |      0 |
| Q8IX19     | MCEMP1    | ENSP00000329920 |           26 |             2 |      0 |
| P55822     | SH3BGR    | ENSP00000332513 |           25 |             2 |      0 |
| Q8WWG9     | KCNE4     | ENSP00000281830 |           25 |             3 |      0 |
| Q14626     | IL11RA    | ENSP00000450565 |           24 |             5 |     51 |
| P00387     | CYB5R3    | ENSP00000354468 |           23 |             0 |     30 |
| Q92782     | DPF1      | ENSP00000347716 |           23 |             2 |      0 |
| P35030     | PRSS3     | ENSP00000354280 |           23 |             2 |      0 |
| Q13620     | CUL4B     | ENSP00000384109 |           20 |             1 |      0 |
| A6NJ46     | NKX6-3    | ENSP00000429553 |           20 |            16 |      0 |
| Q92802     | N4BP2L2   | ENSP00000382328 |           19 |            10 |      0 |
| A6NNV3     | SPDYE16   | ENSP00000487772 |           16 |             3 |      0 |
| A1L168     | C20orf202 | ENSP00000383474 |           16 |             2 |      0 |
| Q5H9F3     | BCORL1    | ENSP00000437775 |           15 |             0 |    133 |
| Q07654     | TFF3      | ENSP00000430690 |           15 |             1 |      0 |
| Q8WUH1     | CHURC1    | ENSP00000475473 |           14 |             0 |      0 |
| P19021     | PAM       | ENSP00000306100 |           12 |             1 |    129 |
| Q8WU43     | C2orf15   | ENSP00000497775 |           11 |             1 |      0 |
| Q969V6     | MRTFA     | ENSP00000498671 |            5 |             0 |      0 |
| Q16181     | SEPTIN7   | ENSP00000413507 |            3 |             0 |      0 |
| Q5JWR5     | DOP1A     | ENSP00000237163 |            1 |             0 |      0 |
| Q13347     | EIF3I     | ENSP00000362688 |            1 |             1 |     22 |
| P01861     | IGHG4     | ENSP00000493388 |            1 |             0 |    134 |



