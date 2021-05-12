# GFF3 格式

Generic Feature Format Version 3

## Columns

  * Column 1: "seqid"
  * Column 2: "source"
  * Column 3: "type"
  * Columns 4 & 5: "start" and "end"
  * Column 6: "score"
  * Column 7: "strand"
  * Column 8: "phase"
  * Column 9: "attributes"

## Attributes column
* Predefined tags in col 9 (attributes)
  * ID
  * Name
  * Alias
  * Parent
  * Target
  * Gap
  * Derives_from
  * Note
  * Dbxref
  * Ontology_term
  * Is_circular

* Characters with reserved meanings in column 9
  * ; semicolon (%3B)
  * = equals (%3D)
  * & ampersand (%26)
  * , comma (%2C)

* Attributes which can have multiple values (split by comma - ,)
  * Parent, Alias, Note, Dbxref ,and Ontology_term

* Attribute names are case sensitive.
  * Parent != parent

# GFF for yeasts from NCBI

## NCBI stored GFF data
* 可以通过 ftp 和 https 两种协议得到相同的目录
* 以下两个链接等价（Pombe 的数据）
  * https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/945/GCF_000002945.1_ASM294v2/
  * ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/945/GCF_000002945.1_ASM294v2/
* 不同注释来源
  * 在 genomes 目录下有三个目录内的文件路径和名称相同，但是内容不同
  * /genomes/all
  * /genomes/genbank
  * /genomes/refseq
  * genbank 和另外两个从数据上本身就会有很大差别
  * 但 all 和 refseq 中也有一些区别，虽然 source 列均为 RefSeq（只针对 4 中酵母的 GFF 文件）

## all 和 refseq 目录下 GFF 的差别

* 4 个 yeast 的 GFF 中，Pombe 因为注释较全，在 feature 列比其他三个多出约一倍的类型

  * 其中有一个 feature 在 all 中为 sequence_difference，并在 attributes 中有两项为 Note=conflict 和 gbkey=misc_difference

  * 而在 refseq 相同的路径的文件中为 sequence_conflict，并在 attributes 没有 Note 一项，gbkey 为 gbkey=conflict

* 以下 feature 的 ID 在 all 中为基因组位置 ID=id-NC_003421.2:1338711..1338730，

  * origin_of_replication
  * repeat_region
  * long_terminal_repeat
  * gap
  * promoter

* 以下 feature 的 ID 在 all 中为基因组的位置，如 ID=id-NC_003421.2:1338711..1338730，而在 refseq 中为 id+对应基因ID，如ID=id-SPCC132.05c

  * five_prime_UTR
  * three_prime_UTR
  * 上面两个因为在 refseq 中使用了 id-基因ID，因此占用了一个 id- 位
    * 会出现在 sequence_feature 这样一个 feature 下，ID 从 all 中的 ID=id-SPCC132.05c
    * 变为 refseq 中的 ID=id-SPCC132.05c**-2**（refseq 中也使用了 id-SPCC132.05c 而不是从 1 计数的 id-SPCC132.05c-1）

* 以上是 pobme 中特有的 feature，并因此产生的不同，下面是对另外三类 yeast 的差别

* 对于 pombe 在 all 或 refseq 目录，或其他三种 yeast 的 all 目录下的 GFF 文件

  * Attributes 的第一个 key 为 ID=...
  * 并且后面的值除了上面的 pombe 特有部分外
    * region 为 ID=NW_011627860.1:1..4291692
    * gene 为 ID=gene-SJAG_01096
    * RNA 为 ID=rna-XM_002172325.2（任意种类 RNA，包括 transcript 和 pseudogenic_tRNA）
    * exon 为 ID=exon-XM_002172325.2-1
    * CDS 为 ID=cds-XP_002172361.1

* 而对于另三种 yeast 在 refseq 目录下的 GFF 文件（和上面不同的是，即使有 gene，rna，cds 这样的关键词，没有短横线连接）

  * gene 为 ID=gene0
  * RNA 为 ID=rna0
  * CDS 为 ID=cds0
  * region 和 exon 则为 ID=id0，ID=id1

* 对于 Attributes 中的 Name=...

  * 对于 pombe，Name 即可能和 ID 完全相同，如 Name=SPAC212.10，也可能为基因名，如 Name=tlh1
  * 对于另三种 yeast，Name 均为系统命名，如 SJAG_01096（因此在 all 目录下的 GFF 中 ID 和 Name 是相同的）

# 需要从 GFF3 中得到的信息

* 以下只针对 GFFParser

## 四种 yeast 的 GFF 中不同 feature 的信息

* 对于所有四种 yeast，第三列的 "type" 有以下 7 种值
  * `'region', 'gene', 'mRNA', 'rRNA', 'tRNA', 'exon', 'CDS'`
* pombe 和 japonicus 有额外一种 pseudogene，和 gene 是等价的
* pombe 中额外还有一些
  * 和 RNA 等价的 `'transcript', 'pseudogenic_tRNA'`
  * 作为某一 gene 附属信息的`'sequence_feature', 'intron', 'promoter', 'sequence_conflict', 'three_prime_UTR', 'five_prime_UTR'`（ID 中含有对应基因的 ID）
  * 作为基因组信息的`'polyA_site', 'polyA_signal_sequence', 'gap', 'repeat_region', 'long_terminal_repeat', 'origin_of_replication'`
* 其中
  * region 只会出现在每个 chrom 的第一行，标注 chrom 的长度
  * gene 和 pseudogene 等价，具有完整的 gene 信息，ID 为 gene-xxxxx
  * `'mRNA', 'rRNA', 'tRNA', 'transcript', 'pseudogenic_tRNA'` 五者等价，ID 为 rna-xxxx，Parent 为对应的 gene ID (gene-xxxx)
  * exon 和 CDS 接于 mRNA 和 transcript 之后，ID 为 exon-xxxx, cds-xxxx，Parent 为对应的 RNA ID (rna-xxxx)
* 关于同一 parent 下有多个相同 type 连续排列
  * 相同 gene 的两个 RNA 的 ID 不会相同，一般会保持连续，如 rna-NM_001356130.1 和 rna-NM_001356131.1
  * exon 无论是单个还是多个出现，均会在 RNA 名字后面加标号，如 exon-NM_001356127.1-1 和 exon-NM_001356127.1-2
  * CDS 的 ID 为蛋白 ID，因此除了 start 和 end 外其他行均不变，如两个 CDS 均为 cds-NP_594307.1
* 关于标识符
  * 目前使用 ID 和 parent 来做 parsing，且只保留共有的 8 种和 pombe 的额外 2 种 RNA 等价的 type
  * 在 attributes 中，Dbxref 这一项会出现不一致的情况，不考虑使用
* 同一项 attribute 具有多个值
  * 虽然 Parent 中也允许以逗号分隔多个值，但在四种酵母数据中均未出现，不考虑这种情况

## 需要的 feature

* 只保留 10 种 feature，7 个通用，1 个 在 pombe 和 japonicus 存在，2 个 pombe 特有

  * region
  * gene
  * pseudogene
  * mRNA
  * transcript
  * rRNA
  * tRNA
  * pseudogenic_tRNA
  * exon
  * CDS

* 结构顺序为

  |            | parent 为 gene   | parent 为 RNA |
  | ---------- | ---------------- | ------------- |
  | gene       | mRNA             | exon          |
  | pseudogene | rRNA             | CDS           |
  |            | tRNA             |               |
  |            | transcript       |               |
  |            | pseudogenic_tRNA |               |

## 需要提取的信息

1. 每个 feature 为 region 的行
   * 记录 chromosome 名称 NC_003424.3
   * 记录长度，1, 5579133（按 start 和 end site 取，不使用 ID）
2. 每个 feature 为 gene 的行
   * 第一列所属的 chromosome
   * 第 3 列 gene type: gene 或 pseudogene
   * 第 4 和 5 列 start 和 end site
   * 第 7 列正链或反链
   * Attributes 列（required 为所有行均需要这一步，没有找到则报错）
     * (required) 得到 ID 中的命名段，如 ID=gene-SPAC212.11 中的 SPAC212.11
     * (optional) 得到 Name，因为存在 pombe 中 Name 有已命名的基因情况，可能得到 ID 中的名称，也可能为基因名 Name=tlh1，可以用来以 gene name 查询系统命名
     * (optional) GeneID，在 Dbxref 字段下，如 Dbxref=GeneID:2541932（如果是 RNA 还可能有 GenBank ID）
     * (optional) 得到 gene，如 gene=tlh1
     * (optional) 得到 gene_biotype，如 gene_biotype=protein_coding
     * (optional) gene_synonym，如 gene_synonym=meu1,meu1-2,meu2,SPAC1556.06,SPAC1556.06b，以逗号分隔的同名基因名，例中的实际保留 gene name 为 gene=meu1-1
3. 每个 feature 为 RNA 类型的行
   * 第 3 列 RNA type
   * 第 4 和 5 列 start 和 end site
   * Attributes 列
     1. ID 中的命名段，如 ID=rna-NM_001018168.1 中的 NM_001018168.1
     2. Parent 中的 gene ID 段，如 Parent=gene-SPAC212.11 中的 SPAC212.11
     3. (optional) Genbank，在 Dbxref 字段下，一般和 GeneID 在一项中，如 Dbxref=GeneID:2541932,Genbank:NM_001018168.1 中的 NM_001018168.1
     4. (optional) Transcript ID，如 transcript_id=NM_001018168.1
     5. (optional) product，如 product=RecQ type DNA helicase
4. 每个 feature 为 exon 的行
   * 第 4 和 5 列 start 和 end site
   * Attributes 列
     1. ID 中的命名段，如 ID=exon-NM_001018168.1-1 中的 NM_001018168.1-1
     2. Parent 中的 gene ID 段，如 Parent=rna-NM_001018168.1 中的 NM_001018168.1
5. 每个 feature 为 CDS 的行
   * 第 4 和 5 列 start 和 end site
   * Attributes 列
     1. ID 中的命名段，如 ID=cds-NP_595040.1 中的 NP_595040.1
     2. Parent 中的 gene ID 段，如 Parent=rna-NM_001018168.1 中的 NM_001018168.1
     3. protein_id，如 protein_id=NP_595040.1

* 最终的格式在 Yeaseq-OfflineDataFormat-GFFParser.md 中

# 已知的 GFF 文件问题

## 没有匹配 RNA 的 exon

1. Pombe 的一个 exon ID 为 chrom 的一段位置
   * 在 pombe 的 GFF 文件中 (https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/945/GCF_000002945.1_ASM294v2/GCF_000002945.1_ASM294v2_genomic.gff.gz)
   * 第 36780 行
   * feature 为 exon
   * ID=id-NC_001326.1:1..5289-1
   * 这一行接在 region 的行之后
2. Pombe 的两个 gene 没有 RNA，但有 exon
   * 在 pombe 的 GFF 文件中 (https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/945/GCF_000002945.1_ASM294v2/GCF_000002945.1_ASM294v2_genomic.gff.gz)
   * 第 36805, 36806 行，exon ID 为 ID=id-ScpofMp01-2 和 ID=id-ScpofMp01-3，parent 为 gene-ScpofMp01
   * 第 36822, 36823 行，exon ID 为 ID=id-ScpofMp04-1 和 ID=id-ScpofMp04-2，parent 为 gene-ScpofMp04
3. Pombe 的 10 个 gene 没有 RNA，但有 CDS
   * 在 pombe 的 GFF 文件中 (https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/945/GCF_000002945.1_ASM294v2/GCF_000002945.1_ASM294v2_genomic.gff.gz)
   * parent 为 gene-ScpofMp01 到 gene-ScpofMp10



