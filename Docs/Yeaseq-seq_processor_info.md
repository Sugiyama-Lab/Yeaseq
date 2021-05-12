~~没有的值显示 None~~

只显示有值的信息


## 包含在每个 Transcript 的 'Description' Key 中

| Unified key | Biomart key | Offline GFF key | Description |
| ------ | ------ | ------ | ------ |
| EnsemblGeneID | ensembl_gene_id |  | Ensembl gene id |
| GeneSymbol |  | GeneSymbol | Gene symbol |
| EnsemblTranscriptID | ensembl_transcript_id |  | Ensembl transcript id |
| GenbankTranscriptID |  | RNAGenbankID | Genbank transcript id |
| EntrezGeneID |  | EntrezID | Entrez gene id |
| ExGeneName | external_gene_name |  | External gene name |
| GeneDescription | description |  | Gene description |
| ChromosomeName | chromosome_name | ChromosomeName | Chromosome name |
| GeneBiotype | gene_biotype | GeneBiotype | Gene biotype |
| TranscriptType | transcript_biotype | RNAType | Transcript type |
| GeneType |  | GeneType | Gene type |
| Strand | strand | Strand | [1, -1] |

## Seq 信息共6种

注意 UTR 区的 intron

+ upstream
```html
[('upstream', site, site)]
```

+ 5'UTR
```html
[('UTR_5', site, site), ('intron', site, site), ('UTR_5', site, site), ...]
```

+ Exon
```html
[('Exon', site, site), ('Exon', site, site), ('Exon', site, site), ...]
```

+ Intron
```html
[('Intron', site, site), ('Intron', site, site), ('Intron', site, site), ...]
```

+ 3'UTR
```html
[('UTR_3', site, site), ('intron', site, site), ('UTR_3', site, site), ...]
```

+ downstream
```html
[('downstream', site, site)]
```

## seq process result

+ max_upstream

+ max_downstream

+ Description ()

+ seq

1. upstream
2. downstream
3. gene

+ SeqSite

1. Exon
2. Intron
3. UTR_5
4. UTR_3

+ SeqFuture

保存有哪些 future，用于 ui 显示可选 future

['Exon', 'Intron', 'UTR_5', 'UTR_3']
