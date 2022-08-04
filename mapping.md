## Code

First and foremost, I am not a geneticist. I use Python and have never used R.
I may be missing a beat, e.g. some amazing table. I simply googled API and was diffident of everything.

For Ensembl, there does not seem to be a list of 
canonical proteins (ENSP), so a double jump is needed.

First let's get the list of canonical ENST (transcripts):

```python
import pandas as pd
import operator
!wget http://ftp.ensembl.org/pub/release-107/tsv/homo_sapiens/Homo_sapiens.GRCh38.107.canonical.tsv

cano = pd.read_table('Homo_sapiens.GRCh38.107.canonical.tsv', header=None, names=['ENSG', 'ENST', 'Type'])
cano_enst = cano.loc[cano.Type == 'Ensembl Canonical'].ENST.str.split('.').apply(operator.itemgetter(0)).to_list()
```

And let's also get the sequences of the ENSP:

```python
!wget http://ftp.ensembl.org/pub/release-107/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz

import gzip, re
from Bio import SeqIO
import pandas as pd
import operator
 
data = []
    
with gzip.open('Homo_sapiens.GRCh38.pep.all.fa.gz','rt') as fh:
    for record in SeqIO.parse(fh, "fasta"):
        # transcript:ENST00000634605.1 gene_biotype:TR_V_gene transcript_biotype:TR_V_gene gene_symbol:TRBV7-2
        rex = re.search(r'transcript:(?P<trans>ENST\d+\.\d+) .* gene_symbol:(?P<symbol>\S+)', record.description)
        if not rex:
            print(record.description)
            continue
        data.append(dict(ENSP=record.name,
                         sequence=str(record.seq), 
                         ENST=rex.group('trans'), 
                         symbol=rex.group('symbol')))

seqs = pd.DataFrame(data)
seqs['unversioned_ENST']= seqs.ENST.str.split('.').apply(operator.itemgetter(0))

# cross:
seqs['is_canonical'] = seqs.unversioned_ENST.isin(cano_enst)  # noqa

seqs
```

| name              | sequence      | enst               | symbol   | unversioned_ENST   | is_canonical   |
|:------------------|:--------------|:-------------------|:---------|:-------------------|:---------------|
| ENSP00000340083.3 | MKHSK...KTYDS | ENST00000347055.4  | KRCC1    | ENST00000347055    | True           |
| ENSP00000343223.6 | MASEA...SVRLG | ENST00000340722.8  | TCL1B    | ENST00000340722    | True           |
| ENSP00000398594.2 | MLLLG...LLLLP | ENST00000448314.6  | CYP21A2  | ENST00000448314    | True           |
| ENSP00000273951.8 | MKRVL...VLLLA | ENST00000273951.13 | GC       | ENST00000273951    | True           |
| ENSP00000347486.3 | MPVSK...CPKKS | ENST00000355329.7  | MORN3    | ENST00000355329    | True           |
| ENSP00000287766.4 | MATNG...SKVAD | ENST00000287766.10 | SLC6A1   | ENST00000287766    | True           |
| ENSP00000268485.3 | MEDPS...GAREP | ENST00000268485.8  | MARVELD3 | ENST00000268485    | True           |
| ENSP00000436691.1 | MADQL...LRKKR | ENST00000530950.2  | CARD18   | ENST00000530950    | True           |
| ENSP00000455079.2 | MRASR...SPPSP | ENST00000567970.2  | C16orf95 | ENST00000567970    | True           |
| ENSP00000007510.6 | MVARS...TDSLD | ENST00000007510.9  | ARHGAP33 | ENST00000007510    | True           |


For Uniprot, I have a preparsed set of files for Michelanglo w/ gnomAD 3 missenses loaded already.
But with the new Uniprot format using `requests` to get a JSON is really easy to get the former.
Anyway, gene symbol to Uniprot:

```python
import os
import json
from michelanglo_protein import ProteinAnalyser, Structure, Mutation, global_settings
global_settings.startup(os.path.expandvars('$HOME/michelanglo/protein-data'))

human = json.load(open(os.path.join(global_settings.dictionary_folder, 'taxid9606-names2uniprot.json')))
```

Now for the mapping, I have a handy function, which uses `Bio.pairwise2.align.globalxx` to align the sequences,
and return a dictionary where the key is the index (starting from one) in the first sequence,
and the value is the index in the second.
It's not rocket science but it is very handy!

```python
from Bio import pairwise2
from typing import Dict, Union 

def get_mapping(seqA:str, seqB:str) -> Dict[int, Union[int, float]]:
    """
    Return a one-indicised mapping of positions in seqA to one-indicised position in seqB
    Missing values (gaps) are NaNs
    
    Alignment scoring is 1 for a match, -4 for a mismatch, -5 for a new gap and -0.02 for extensions
    """
    alignment = pairwise2.align.globalms(seqA, seqB, 1, -4, -5, -0.02)[0]
    assert alignment.score > 0, f'Sequences {seqA} and {seqB} are very different'
    map_pos2al = lambda seq: dict(enumerate([i for i, l in enumerate(seq) if l != '-']))
    flip_dict = lambda d: dict(zip(d.values(), d.keys()))
    B_from_al = flip_dict(map_pos2al(alignment.seqB))
    # A2B is one indicised
    return {a+1: B_from_al.get(al, float('nan'))+1 for a, al in map_pos2al(alignment.seqA).items()}
```


```python
import os, json
import pandas as pd
from typing import *
from michelanglo_protein import ProteinAnalyser

seqs: pd.DataFrame  # from above --> ensembl
human: Dict[str, str]  # from above --> symbol to uniprot
get_mapping: Callable[[str, str], Dict[int, Union[int, float]]]  # from above

os.mkdir('ensp2uniprot_maps')
from typing import *

details: List[Dict[str, Any]] = []
ensp2uniprots: Dict[Tuple[str, str], Dict[int, Union[int, float]]] = {}

for i, row in seqs.loc[seqs.is_canonical].iterrows():
    detail = dict(symbol=row.symbol, 
                  ESNP=row['name'],
                  ESNP_len=len(row.sequence),
                  indentity='No match')
    if row.symbol not in human:
        details.append(detail)
        continue
    detail['uniprot'] = human[row.symbol]
    protein = ProteinAnalyser(uniprot=detail['uniprot'], taxid=9606).load()
    detail['uniprot_len']=len(protein.sequence)
    if protein.sequence == row.sequence:
        detail['identity'] = 'Identity'
        details.append(detail)
        continue
    if detail['ESNP_len'] == detail['uniprot_len']:
        detail['identity'] == 'Missenses'
    else:
        detail['identity'] = 'Differing'
    details.append(detail)
    try:
        ensp2uniprot = get_mapping(row.sequence, protein.sequence)
        with open(os.path.join('ensp2uniprot_maps', f'{detail["ESNP"]}-{detail["uniprot"]}.json'), 'w') as fh:
            json.dump(ensp2uniprot, fh)
        ensp2uniprots[detail['ESNP'], detail['uniprot']] = ensp2uniprot
    except Exception as error:
        detail['identity'] = 'Unmappable'

deets = pd.DataFrame(details)
deets.to_csv('ENSP-uniprot_check.csv', index=False)
deets
```

| symbol   | ESNP              |   ESNP_len | identity   | uniprot   |   uniprot_len |
|:---------|:------------------|-----------:|:-----------|:----------|--------------:|
| OR6X1    | ENSP00000333724.2 |        312 | Identity   | Q8NH79    |           323 |
| TSBP1    | ENSP00000480415.1 |        565 | Differing  | Q5SRN2    |           184 |
| STK17B   | ENSP00000263955.4 |        372 | Identity   | O94768    |           351 |
| WASHC2C  | ENSP00000485513.1 |       1341 | Differing  | Q9Y4E1    |           125 |
| TRAJ10   | ENSP00000451066.1 |         21 | No match   | nan       |           nan |
| SEPTIN4  | ENSP00000500383.1 |        996 | Differing  | O43236    |           526 |
| CFAP73   | ENSP00000333915.6 |        308 | Identity   | A6NFT4    |           195 |
| CNTFR    | ENSP00000368265.3 |        372 | Identity   | P26992    |            76 |
| MPPE1    | ENSP00000465894.1 |        396 | Identity   | Q53F39    |           328 |
| IDH3A    | ENSP00000299518.2 |        366 | Identity   | P50213    |           163 |

Crude plotly waffle chart:
```python
import plotly.graph_objects as go
import numpy as np
import pandas as pd

def waffle(series: pd.Series) -> go.Figure:
    """This snippet is very crude"""
    series = series.sort_values()
    d = dict(enumerate(series.unique()))
    l = len(series)
    w = round((3*l/2)**0.5)
    h = l//w + 1
    z = np.ones((h, w)) *np.nan  # row, column
    cd = np.ones((h, w)).astype(str)  # row, column
    percentages = ((series.value_counts() * 100) // len(series)).to_dict()
    dl = series.unique().tolist()
    for i, v in tuple(enumerate(series)):
        z[i//w, i%w] = dl.index(v)
        cd[i//w, i%w] = f'{v} ({percentages[v]})%' 

    fig = go.Figure(go.Heatmap(z=z,xgap=0, ygap=0,showscale=False,customdata=cd, hovertemplate="%{customdata}"))
    fig.update_layout(width=600, height=400)
    return fig

waffle(deets.indentity)
```
Fix gnomADs in my datasets —ignore this!
Briefly, a Variant is a namedtuple, which cannot be changed and `._replace()` returns a copy.
```python
changes = {'same': 0, 'good_change': 0, 'bad_change': 0, 'error': 0}
for ensp, uniprot in ensp2uniprots: # tuple index
    protein = ProteinAnalyser(uniprot=uniprot, taxid=9606).load()
    new_gnomads = []
    ensp2uniprot = ensp2uniprots[ensp, uniprot]
    for var in protein.gnomAD:
        new_n = ensp2uniprot.get(var.x, float('nan'))
        if str(new_n) == 'nan':
            changes['error'] += 1
            continue
        if new_n == var.x:
            new_gnomads.append(var)
            changes['same'] += 1
            continue
        new_d = var.description.replace(f'{var.from_residue}{var.residue_index}{var.to_residue}', f'{var.from_residue}{new_n}{var.to_residue}')
        new_gnomads.append(var._replace(x=new_n, y=new_n,
                             description=new_d,
                             residue_index=new_n,
                             precomputed_ddG=None)
                          )
        if protein.sequence[new_n - 1] == var.from_residue:
            changes['good_change'] += 1
        else:
            changes['bad_change'] += 1
print(changes)
```