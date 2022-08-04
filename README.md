# Uniprot-ported-gnomADs
These are datasets of gnomAD missense/nonsenses ported to Uniprot numbering.

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
Turns out 81% match up between Uniprot and Ensembl.
15% had differences in the sequences. 3% I could not match the sequences.
And less than 1% did not have a corresponding SwissProt Uniprot entry 
—SwissProt is curated, Trembl has all the "novel protein".

