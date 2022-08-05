# Uniprot-ported-gnomADs
These are datasets of gnomAD missense/nonsenses ported to Uniprot numbering.

> :construction: WIP Currently going through the gnomAD data to get the ENSP ids used as 1,793 protein do not match

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

## Mapping

> See [mapping](mapping.md) for Python code

I was told the match up between Uniprot and Ensembl is pretty good. But I never checked...
A major problem I spotted is that some entries have cross-database links to canonical transcripts,
whereas in reality it's to other isoforms, so I mistrusted these.
Turns out 81% match up between Uniprot and Ensembl. 7% do not actually match the current Ensembl canonical.
3% had just some missenses between the sequences.
3% had differences in the sequences.
And less than 1% (44) did not have a valid corresponding SwissProt Uniprot entry 
—SwissProt is curated, Trembl has all the "novel protein".

