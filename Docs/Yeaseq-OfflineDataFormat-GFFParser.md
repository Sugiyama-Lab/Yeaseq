## 以不同 chromosome 为名称的 json 文件

* 从 NCBI 下载的 GFF 和 Fna 文件得到
* 分为两部分
* chromosome 信息

| Key | content |
| ------ | ------ |
| ChromosomeName | e.g. 'NW_011627891.1' |
| ChromosomeSite | e.g. [1, 2395] |
| Seq | e.g. 'TTATGCG...AACAGTC' |

+ Gene 信息

| Gene ID | Key | Key | Key | Key | e.g. |
| ------ | ------ | ------ | ------ | ------ | ------ |
| gene-SJAG_00941 | Gene | GeneSite |  |  | [1, 5662] |
|  |  | ID                         |              |              | 'gene-SPAC212.11'                    |
|  |  | Symbol | | | 'SPAC212.11' |
|  |  | GeneType |  |  | 'gene' |
|  |  | Strand |  |  | -1 |
|  |  | EntrezID (GeneID)          |  |  | '2541932' |
| | | Name                       |  |  | 'tlh1' |
| | | Gene (gene)                |  |  | 'tlh1' |
| | | Synonym (gene_synonym)     |  |  | ['meu1', 'meu1-2', 'meu2'] |
|  |  | GeneBiotype (gene_biotype) |  |  | 'protein_coding' |
|  | Transcript | 'rna-XM_002175672.2' | RNA | ID | 'rna-XM_002175672.2'                 |
| |  |  |  | ParentID | 'gene-SPAC212.11' |
| |  |  |  | RNASite | [491037, 496801] |
|  |  |  |  | RNAType | 'mRNA' |
|  |  |  |       | Genbank      | 'XM_002175672.2' |
| | | |       | Product      | RecQ type DNA helicase |
|  |  |  | Exon       |          | [[491037, 492168], [492308, 496801]] |
| | |  |  |  |  |
|  |  |  | CDS        |            | [[491043, 492168], [492308, 493845]] |
| | |  |     |          |                                      |
|             | |  |  |          |  |
|  | TranscriptNumber |  |  |  | 1 |

* Gene 可能有多个产物，因此 Transcript 为包含所有产物的列表，每个产物的格式和表中相同，数量储存在 TranscriptNumber 中
* 上表为一般情况
* 有以下几种其他可能：
  1. GeneType 为 pseudogene，一般 product 为 ncRNA
  2. RNA 为 tRNA 或 rRNA，此时
     1. GeneBiotype 为 tRNA 或 rRNA
     2. RNAType 相应改变
     3. RNAGenbankID 为空，即 ''
     4. 无 CDS 一项
* 其中 Strand 为 1 代表正链，-1 为负链

## 以 id_to_chromosome.json 为名称的文件

* 通过 ID 查询 chromosome name

```html
{'SJAG_01096': 'NW_011627860.1',
 'SJAG_01095': 'NW_011627860.1',
 'SJAG_01093': 'NW_011627860.1',
 'SJAG_01091': 'NW_011627860.1',
...
}
```

## 以 symbol_to_id.json 为名称的文件

* 通过 symbol 查询 ID

## 以 entrez_to_id.json 为名称的文件

* 通过 entrez gene id 查询 ID

```html
{'7047373': 'SJAG_01096',
 '7047374': 'SJAG_01095',
 '7048336': 'SJAG_01093',
 '7052043': 'SJAG_01091',
...
}
```

## 以 product_to_id.json 为名称的文件

* 通过 product 的 Genbank ID 查询 ID
* 包括 transcript 和 protein 两个 level

```html
{'XM_002172325.2': 'SJAG_01096',
 'XM_002172324.2': 'SJAG_01095',
 'XM_002172323.2': 'SJAG_01093',
 'XM_002172322.2': 'SJAG_01091',
...
}
```

## 以 alias_to_id.json 为名称的文件

* 将一个 gene 的所有可得到信息映射到 ID 上
* 包括 ID -> ID 自身的映射

