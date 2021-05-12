from urllib.request import urlopen
from requests import get as requests_get
from zlib import decompress
from os import path
from os import makedirs
from re import match as re_match
from re import sub as re_sub
from json import dump
from json import load

from .yeaseq_config import WorkDir
from .yeaseq_config import NCBIFTP
from .yeaseq_config import Yeast


class GFFParser(object):
    """
    Doc: Yeaseq-GFFParser.md

    这里不区分 yeast 种类，全部按照 10 种 feature 解析
    'region', 'gene', 'pseudogene', 'mRNA', 'transcript', 'exon', 'CDS', 'rRNA', 'tRNA', 'pseudogenic_tRNA'

    这里所有 re 均使用 match 方法，不使用 compile，且单次匹配每个 item
        re.match('.*GeneID:(\\d+)', attr)           971 ns ± 6.11 ns per loop (mean ± std. dev. of 7 runs, 1000000 loops each)
        re.match(r, attr)                           1.23 µs ± 5.47 ns per loop (mean ± std. dev. of 7 runs, 1000000 loops each)
        re.match('ID=(.+);', attr)                  866 ns ± 18.7 ns per loop (mean ± std. dev. of 7 runs, 1000000 loops each)
        re.match('ID=(.+);.*GeneID:(\\d+)', attr)   3.5 µs ± 194 ns per loop (mean ± std. dev. of 7 runs, 100000 loops each)
        re.match('.*GeneID:(\\d+)', attr).groups()  1.07 µs ± 7.85 ns per loop (mean ± std. dev. of 7 runs, 1000000 loops each)
        re.findall('.*GeneID:(\\d+)', attr)         16.4 µs ± 62.8 ns per loop (mean ± std. dev. of 7 runs, 100000 loops each)
    """
    def __init__(self, gff_data):
        self.raw_gff = gff_data
        self.all_features = [
            'region',
            'gene', 'pseudogene',
            'mRNA', 'transcript', 'exon', 'CDS',
            'rRNA', 'tRNA', 'pseudogenic_tRNA'
        ]
        self.feature_level_gff = {k: [] for k in self.all_features}

        self.GFFTitle = {
            'Reference': 0,
            'Source': 1,
            'Feature': 2,
            'Start': 3,
            'End': 4,
            'Score': 5,
            'Strand': 6,
            'Frame': 7,
            'Attributes': 8,
        }
        self.GFFTitleName = list(self.GFFTitle.keys())

        self.chrom_info = dict()
        self.id_to_chromosome = dict()
        self.symbol_to_id = dict()
        self.entrez_to_id = dict()
        self.product_to_id = dict()
        self.alias_to_id = dict()

        self.parse_gff()
        self._cleanup_all_data()

    def parse_gff(self):
        print('Divide GFF')
        self._divide_feature()
        print('Parse region')
        self._parse_region_feature()
        print('Parse gene')
        self._parse_gene_feature()
        print('Parse rna')
        self._parse_rna_feature()
        print('Parse exon')
        self._parse_exon_feature()
        print('Parse cds')
        self._parse_cds_feature()
        print('GFF parse done')

    def _divide_feature(self):
        if isinstance(self.raw_gff, str):
            gff_list = self.raw_gff.split('\n')
        elif isinstance(self.raw_gff, list):
            gff_list = self.raw_gff
        else:
            raise ValueError(f'The input of GFFParser should be str or list, now is {type(self.raw_gff)}')
        for idx, row in enumerate(gff_list):
            if row.startswith('#') or row == '' or row == '\n':
                continue
            else:
                split_row = row.split('\t')
                row_feature = split_row[self.GFFTitle['Feature']]
                if row_feature in self.all_features:
                    self.feature_level_gff[row_feature].append(split_row)
                else:
                    pass
                    # raise KeyError(f'The feature type "{row_feature}" of line {idx} '
                    #                f'in current GFF data is not in valid feature list {self.all_features}')

    def _parse_region_feature(self):
        for row in self.feature_level_gff['region']:
            ch_name = row[self.GFFTitle['Reference']]
            init_ch_info = dict()
            init_ch_info['ChromosomeName'] = ch_name
            init_ch_info['ChromosomeSite'] = (row[self.GFFTitle['Start']], row[self.GFFTitle['End']])
            self.chrom_info[ch_name] = init_ch_info

    @staticmethod
    def _re_match(partten, string):
        m = re_match(partten, string)
        if m is None:
            return ''
        else:
            return m.groups()[0]

    def _parse_gene_feature(self):
        for gene_type in ['gene', 'pseudogene']:
            for row in self.feature_level_gff[gene_type]:
                ch_name = row[self.GFFTitle['Reference']]

                init_gene_info = dict()
                attr = row[self.GFFTitle['Attributes']]

                gene_id = self._re_match('.*ID=(.+?);', attr)
                if gene_id == '':
                    print(row)
                    continue
                symbol = gene_id.replace('gene-', '')

                init_gene_info['ID'] = gene_id
                init_gene_info['Symbol'] = symbol
                init_gene_info['GeneSite'] = (row[self.GFFTitle['Start']], row[self.GFFTitle['End']])
                init_gene_info['GeneType'] = gene_type
                init_gene_info['Strand'] = row[self.GFFTitle['Strand']]
                init_gene_info['EntrezID'] = self._re_match('.*GeneID:(\\d+)', attr)
                init_gene_info['Name'] = self._re_match('.*Name=(.+?);', attr)
                init_gene_info['Gene'] = self._re_match('.*gene=(.+?);', attr)
                init_gene_info['Synonym'] = self._re_match('.*gene_synonym=(.+?);', attr)
                init_gene_info['GeneBiotype'] = self._re_match('.*gene_biotype=(.+?);', attr)

                self.chrom_info[ch_name][gene_id] = dict()
                self.chrom_info[ch_name][gene_id]['Gene'] = init_gene_info
                self.chrom_info[ch_name][gene_id]['Transcript'] = dict()
                self.chrom_info[ch_name][gene_id]['TranscriptNumber'] = 0
                self.id_to_chromosome[gene_id] = ch_name

                self.symbol_to_id[symbol] = gene_id
                if init_gene_info['EntrezID'] != '':
                    self.entrez_to_id[init_gene_info['EntrezID']] = gene_id

                for k in ['ID', 'EntrezID', 'Name', 'Gene']:
                    self.alias_to_id[init_gene_info[k]] = gene_id
                if init_gene_info['Synonym'] != '':
                    for each_synonym in init_gene_info['Synonym'].split(','):
                        self.alias_to_id[each_synonym] = gene_id
                self.alias_to_id[symbol] = gene_id
                self._cleanup_data(self.alias_to_id)

    def _parse_rna_feature(self):
        for rna_type in ['mRNA', 'transcript', 'rRNA', 'tRNA', 'pseudogenic_tRNA']:
            for row in self.feature_level_gff[rna_type]:
                ch_name = row[self.GFFTitle['Reference']]
                attr = row[self.GFFTitle['Attributes']]
                rna_id = self._re_match('.*ID=(.+?);', attr)
                parent_id = self._re_match('.*Parent=(.+?);', attr)

                init_rna_info = dict()
                init_rna_info['ID'] = rna_id
                init_rna_info['ParentID'] = parent_id
                init_rna_info['RNASite'] = (row[self.GFFTitle['Start']], row[self.GFFTitle['End']])
                init_rna_info['RNAType'] = rna_type
                init_rna_info['Genbank'] = self._re_match('.*Genbank:([a-zA-Z]+?_[\\d.]+);', attr)
                init_rna_info['Product'] = self._re_match('.*;product=(.+?);', attr)

                self.chrom_info[ch_name][parent_id]['TranscriptNumber'] += 1
                self.chrom_info[ch_name][parent_id]['Transcript'][rna_id] = dict()
                self.chrom_info[ch_name][parent_id]['Transcript'][rna_id]['RNA'] = init_rna_info

                for k in ['ID', 'Genbank', 'Product']:
                    self.product_to_id[init_rna_info[k]] = parent_id
                self.product_to_id[rna_id.replace('rna-', '')] = parent_id
                self._cleanup_data(self.product_to_id)

    def _parse_exon_feature(self):
        for row in self.feature_level_gff['exon']:
            ch_name = row[self.GFFTitle['Reference']]
            attr = row[self.GFFTitle['Attributes']]
            exon_id = self._re_match('.*ID=(.+?);', attr)
            parent_id = self._re_match('.*Parent=(.+?);', attr)
            gene_id = self.product_to_id.get(parent_id)
            if gene_id is None:
                print(row)
                continue
            sites = (row[self.GFFTitle['Start']], row[self.GFFTitle['End']])

            if 'Exon' not in self.chrom_info[ch_name][gene_id]['Transcript'][parent_id]:
                self.chrom_info[ch_name][gene_id]['Transcript'][parent_id]['Exon'] = [sites]
            else:
                self.chrom_info[ch_name][gene_id]['Transcript'][parent_id]['Exon'].append(sites)

            self.product_to_id[exon_id] = gene_id
            self.product_to_id[exon_id.replace('exon-', '')] = gene_id

            self._cleanup_data(self.product_to_id)

    def _parse_cds_feature(self):
        for row in self.feature_level_gff['CDS']:
            ch_name = row[self.GFFTitle['Reference']]
            attr = row[self.GFFTitle['Attributes']]
            cds_id = self._re_match('.*ID=(.+?);', attr)
            parent_id = self._re_match('.*Parent=(.+?);', attr)
            gene_id = self.product_to_id.get(parent_id)
            if gene_id is None:
                print(row)
                continue
            sites = (row[self.GFFTitle['Start']], row[self.GFFTitle['End']])

            if 'CDS' not in self.chrom_info[ch_name][gene_id]['Transcript'][parent_id]:
                self.chrom_info[ch_name][gene_id]['Transcript'][parent_id]['CDS'] = [sites]
            else:
                self.chrom_info[ch_name][gene_id]['Transcript'][parent_id]['CDS'].append(sites)

            self.product_to_id[cds_id] = gene_id
            self.product_to_id[cds_id.replace('cds-', '')] = gene_id

            self._cleanup_data(self.product_to_id)

    @staticmethod
    def _cleanup_data(d):
        try:
            del d['']
        except KeyError:
            pass

    def _cleanup_all_data(self):
        try:
            del self.chrom_info['']
        except KeyError:
            pass
        try:
            del self.id_to_chromosome['']
        except KeyError:
            pass
        try:
            del self.symbol_to_id['']
        except KeyError:
            pass
        try:
            del self.entrez_to_id['']
        except KeyError:
            pass
        try:
            del self.product_to_id['']
        except KeyError:
            pass
        try:
            del self.alias_to_id['']
        except KeyError:
            pass


class EasyGFFFileParser(object):
    """
    [NOTICE] This is a easy GFF parser used before for only japonicus, cryophilus, and octosporus

    GFF 文件以出现连续两次 ## 开头的行开始一个 chromosome
    第三列为 feature type，共八种 ['region', 'gene', 'pseudogene', 'mRNA', 'rRNA', 'tRNA', 'exon', 'CDS']
        region 为每个 chromosome 开始的第一行，即描述 chromosome 的信息，从这里获得 chromosome 的信息并储存
        gene 和 pseudogene 为一个 gene 的开始
        mRNA、rRNA、tRNA 为一个 transcript 的开始
        exon 为编码 RNA 的描述
        CDS 只在 mRNA 中出现
    Description 部分需要
        gene info
            GeneID:22831291 为 entrez gene id
            Name=SJAG_16127 为 gene name
            gene_biotype=tRNA 为 gene type
        mrna info
            Genbank:XM_002172325.2 为 mrna genbank id

    通过处理 GFF 文件可以得到四个输出
        1. chromosome 的信息，即每个 chromosome 中所有 gene 的信息，注意一个 gene 可能会有多个产物，详细结构在 Yeaseq-offline_storage_format.md
        2. gene name 到 chromosome 的信息，用于查询 gene 时选择对应的 chromosome 取出信息
        3. entrez gene id 到 gene name 的信息，当查询 entrez id 时获得 gene name
        4. XM... 和 XP... 到 gene name 的信息
    """
    def __init__(self, gff_content):
        self.gff_data_list = gff_content.split('\n')

        self._symbol_to_chromosome = dict()  # 2
        self._entrez_to_symbol = dict()  # 3
        self._product_to_symbol = dict()  # 4

        # 1
        self._one_ch = dict()
        self._one_gene = dict()
        self._one_gene_info = dict()
        self._one_transcript = dict()
        self._rna_num = 0

        self.GFFTitle = {'Reference': 0,
                         'Source': 1,
                         'Feature': 2,
                         'Start': 3,
                         'End': 4,
                         'Score': 5,
                         'Strand': 6,
                         'Frame': 7,
                         'Attributes': 8,
                         }
        self.GFFTitleName = list(self.GFFTitle.keys())

    def _gff_row_to_dict(self, row_data):
        row_data = row_data.split('\t')
        row_data_dict = dict(zip(self.GFFTitleName, row_data))
        row_data_dict['Start'] = int(row_data_dict['Start'])
        row_data_dict['End'] = int(row_data_dict['End'])
        row_data_dict['Strand'] = 1 if row_data_dict['Strand'] == '+' else -1
        return row_data_dict

    def generate_gene_info_per_chromosome(self):
        start = False
        for i, each_row in enumerate(self.gff_data_list):
            if each_row.startswith('#'):
                start = True
                if self._one_gene_info:
                    self._save_transcript()
                    self._reset_transcript()
                    self._save_gene()
                    self._reset_gene()
                if self._one_ch:
                    yield self._one_ch
                    self._reset_ch()
                if each_row.startswith('###'):
                    break
            else:
                row_data_dict = self._gff_row_to_dict(each_row)
                row_type = row_data_dict['Feature']

                if start:
                    assert row_type == 'region'
                    self._one_ch = dict()
                    self._one_ch['ChromosomeName'] = row_data_dict['Reference']
                    self._one_ch['ChromosomeSite'] = (row_data_dict['Start'], row_data_dict['End'])
                    start = False

                elif row_type in ['gene', 'pseudogene']:
                    self._save_transcript()
                    self._reset_transcript()
                    self._save_gene()
                    self._reset_gene()

                    self._fill_one_gene_info(row_data_dict)
                    self._save_gene_info()
                    self._add_symbol_to_ch()
                    self._add_entrez_to_symbol()

                elif row_type in ['mRNA', 'rRNA', 'tRNA']:
                    self._save_transcript()
                    self._reset_transcript()

                    self._fill_one_transcript_rna(row_data_dict)
                    self._add_product_to_symbol()

                elif row_type == 'exon':
                    self._fill_one_transcript_exon(row_data_dict)

                elif row_type == 'CDS':
                    self._fill_one_transcript_cds(row_data_dict)

                else:
                    print(f'Undefined type: {row_type} in line {i}')
                    pass

    def get_conversion_symbol2ch(self):
        return self._symbol_to_chromosome

    def get_conversion_entrez2symbol(self):
        return self._entrez_to_symbol

    def get_conversion_product2symbol(self):
        try:
            del self._product_to_symbol['']
        except KeyError:
            pass
        return self._product_to_symbol

    def _fill_one_gene_info(self, row_data):
        self._one_gene_info['GeneSite'] = (row_data['Start'], row_data['End'])
        self._one_gene_info['GeneType'] = row_data['Feature']
        self._one_gene_info['Strand'] = row_data['Strand']
        description = row_data['Attributes']
        _re_match_entrezid = re_match('.*GeneID:(\\d+)', description)
        _re_match_symbol = re_match('.*Name=(.+?);', description)
        _re_match_biotype = re_match('.*gene_biotype=(.+?);', description)
        try:
            self._one_gene_info['EntrezID'] = _re_match_entrezid.group(1)
        except AttributeError:
            self._one_gene_info['EntrezID'] = ''
        try:
            self._one_gene_info['GeneSymbol'] = _re_match_symbol.group(1)
        except AttributeError:
            self._one_gene_info['GeneSymbol'] = ''
        try:
            self._one_gene_info['GeneBiotype'] = _re_match_biotype.group(1)
        except AttributeError:
            self._one_gene_info['GeneBiotype'] = ''

    def _fill_one_transcript_rna(self, row_data):
        self._rna_num += 1
        self._one_transcript['RNA'] = dict()
        self._one_transcript['RNA']['RNASite'] = (row_data['Start'], row_data['End'])
        self._one_transcript['RNA']['RNAType'] = row_data['Feature']
        _re_match_genbank_id = re_match('.*Genbank:([XN]M_[\\d.]+);', row_data['Attributes'])
        try:
            self._one_transcript['RNA']['RNAGenbankID'] = _re_match_genbank_id.group(1)
        except AttributeError:
            self._one_transcript['RNA']['RNAGenbankID'] = ''

    def _fill_one_transcript_exon(self, row_data):
        if 'Exon' not in self._one_transcript:
            self._one_transcript['Exon'] = dict([('ExonSite', [])])
        self._one_transcript['Exon']['ExonSite'].append((row_data['Start'], row_data['End']))

    def _fill_one_transcript_cds(self, row_data):
        if 'CDS' not in self._one_transcript:
            self._one_transcript['CDS'] = dict([('CDSSite', [])])
        self._one_transcript['CDS']['CDSSite'].append((row_data['Start'], row_data['End']))

    def _save_gene(self):
        if self._one_gene:
            self._one_gene['TranscriptNumber'] = self._rna_num
            self._one_ch[self._one_gene_info['GeneSymbol']] = self._one_gene

    def _save_gene_info(self):
        if self._one_gene_info:
            self._one_gene['Gene'] = self._one_gene_info

    def _save_transcript(self):
        if self._one_transcript:
            if 'Transcript' not in self._one_gene:
                self._one_gene['Transcript'] = [self._one_transcript, ]
            elif isinstance(self._one_gene['Transcript'], list):
                self._one_gene['Transcript'].append(self._one_transcript)
            else:
                print('Wrong when saves transcript info')

    def _add_symbol_to_ch(self):
        self._symbol_to_chromosome[self._one_gene_info['GeneSymbol']] = self._one_ch['ChromosomeName']

    def _add_entrez_to_symbol(self):
        self._entrez_to_symbol[self._one_gene_info['EntrezID']] = self._one_gene_info['GeneSymbol']

    def _add_product_to_symbol(self):
        self._product_to_symbol[self._one_transcript['RNA']['RNAGenbankID']] = self._one_gene_info['GeneSymbol']

    def _reset_ch(self):
        self._one_ch = dict()

    def _reset_gene(self):
        self._one_gene = dict()
        self._one_gene_info = dict()
        self._rna_num = 0

    def _reset_transcript(self):
        self._one_transcript = dict()


class FnaParser(object):
    def __init__(self, fna_content):
        self.fna_content = fna_content
        self.seq_dict = dict()

        self._process_fna()

    def _process_fna(self):
        chromosome_seq_list = re_sub(' .+?\n', ' ', self.fna_content).replace('\n', '').split('>')[1:]
        for _ in chromosome_seq_list:
            ch_name, seq = _.split(' ')
            self.seq_dict[ch_name] = seq.upper()

    def get_seq_dict(self):
        return self.seq_dict


class OfflineDataManager(object):
    """用于搜索页面的offline数据管理，在搜索frame的init中初始化"""

    '''        4. version file
    
    启动的时候check对应文件夹是否存在，搜索的时候如果不存在提示，搜索不到提示，文件不存在提示
'''
    def __init__(self):
        _initial_version = ['']
        self.VersionDict = {'pombe': 'GCF_000002945.1_ASM294v2',
                            'japonicus': 'GCF_000149845.2_SJ5',
                            'cryophilus': 'GCF_000004155.1_SCY4',
                            'octosporus': 'GCF_000150505.1_SO6'}

        self.DataExistence = dict(zip(Yeast.YeastList, self.check_data_existence()))
        self.raw_gff_dict = dict()
        self.raw_fna_dict = dict()

    @staticmethod
    def check_data_existence():
        """检查本地数据是否存在，返回True或False的list"""
        data_existence_list = [True if path.isdir(path.join(WorkDir.get_data_dir(n))) else False for n in Yeast.YeastList]
        return data_existence_list

    @staticmethod
    def _make_data_dir(n):
        _dir = WorkDir.get_data_dir(n)
        try:
            makedirs(_dir)
        except FileExistsError:
            pass
        return _dir

    def check_lastest_version(self, n):
        ftp_url = NCBIFTP.get_latest_ftp_url(n)
        version_list = urlopen(ftp_url).readlines()
        for each_version in version_list[::-1]:
            re_matched_gcf_version = re_match('.*?(GCF_.+?) ->.*', str(each_version))
            yeast_name = Yeast.get_yeast_name(n)
            try:
                self.VersionDict[yeast_name] = re_matched_gcf_version.group(1)
            except AttributeError:
                pass

    @staticmethod
    def _download_unzip_gz(ftp_url):
        """给定一个gz文件的ftp地址，返回经过解压的gz文件内容"""
        byte_urldata = urlopen(ftp_url).read()
        required_content = decompress(byte_urldata).decode()
        return required_content

    @staticmethod
    def _download_unzip_gz_request(https_url):
        # TODO 分离 download 和 unzip，更新下载状态 panel
        # TODO 找不到输入 gene 时，检查文件完整性
        """给定一个gz文件的https地址，返回经过解压的gz文件内容"""
        byte_data = requests_get(https_url).content  # TODO Check md5 after download
        try:
            required_content = decompress(byte_data, 31).decode()
        except:
            required_content = decompress(byte_data).decode()
        return required_content

    def download_data(self, n):
        # TODO version check
        # try:
        #     self.check_lastest_version(n)
        # except:
        #     pass
        short_name = Yeast.get_yeast_name(n)
        https_url = NCBIFTP.get_genomes_all_https_url(short_name)
        _version = self.VersionDict[short_name]

        https_gff_gz = f'{https_url}{_version}_genomic.gff.gz'
        https_fna_gz = f'{https_url}{_version}_genomic.fna.gz'
        # ftp_gff_gz = f'{ftp_url}{_version}/{_version}_genomic.gff.gz'
        # ftp_fna_gz = f'{ftp_url}{_version}/{_version}_genomic.fna.gz'

        gff_content = self._download_unzip_gz_request(https_gff_gz)
        fna_content = self._download_unzip_gz_request(https_fna_gz)
        # gff_content = self._download_unzip_gz(https_gff_gz)
        # fna_content = self._download_unzip_gz(https_fna_gz)

        self.raw_gff_dict[short_name] = gff_content
        self.raw_fna_dict[short_name] = fna_content

    def save_data(self, n):
        short_name = Yeast.get_yeast_name(n)

        fna_file_parser = FnaParser(self.raw_fna_dict[short_name])
        fna_seq_dict = fna_file_parser.get_seq_dict()

        gff_parser = GFFParser(self.raw_gff_dict[short_name])

        _dir = self._make_data_dir(short_name)
        print(f'Saving to {_dir}')
        for ch_name, each_ch_data in gff_parser.chrom_info.items():
            each_ch_data['Seq'] = fna_seq_dict[ch_name]
            with open(path.join(_dir, f'{ch_name}.json'), 'w') as f:
                dump(each_ch_data, f)

        with open(path.join(_dir, 'id_to_chromosome.json'), 'w') as f:
            dump(gff_parser.id_to_chromosome, f)

        with open(path.join(_dir, 'symbol_to_id.json'), 'w') as f:
            dump(gff_parser.symbol_to_id, f)

        with open(path.join(_dir, 'entrez_to_id.json'), 'w') as f:
            dump(gff_parser.entrez_to_id, f)

        with open(path.join(_dir, 'product_to_id.json'), 'w') as f:
            dump(gff_parser.product_to_id, f)

        with open(path.join(_dir, 'alias_to_id.json'), 'w') as f:
            dump(gff_parser.alias_to_id, f)


class OfflineDataLoader(object):
    def __init__(self, query_content, n='sj', content_type='Ensembl'):
        self.short_name = Yeast.get_yeast_name(n)
        self._data_dir = WorkDir.get_data_dir(self.short_name)
        self.query_content = query_content
        self._content_type = content_type

        self._ch = None
        self._symbol_to_ch = None
        self._entrez_to_id = None
        self._product_to_symbol = None

    def get_ch_data(self):
        if self._content_type == 'Ensembl':
            self._load_symbol_to_id()
            query_id = self._symbol_to_id.get(self.query_content)
        elif self._content_type == 'NCBI':
            self._load_entrez_to_id()
            query_id = self._entrez_to_id.get(self.query_content)
        else:
            return None
        if query_id is None:
            self._load_product_to_id()
            query_id = self._product_to_id.get(self.query_content)
        if query_id is None:
            self._load_alias_to_id()
            query_id = self._alias_to_id.get(self.query_content)

        self._load_id_to_chromosome()
        _ch_name = self._id_to_chromosome.get(query_id)
        if _ch_name is None:
            return None, query_id
        self._load_target_ch(_ch_name)
        return self._ch, query_id

    def _load_symbol_to_id(self):
        with open(path.join(self._data_dir, 'symbol_to_id.json'), 'r') as f:
            self._symbol_to_id = load(f)

    def _load_entrez_to_id(self):
        with open(path.join(self._data_dir, 'entrez_to_id.json'), 'r') as f:
            self._entrez_to_id = load(f)

    def _load_product_to_id(self):
        with open(path.join(self._data_dir, 'product_to_id.json'), 'r') as f:
            self._product_to_id = load(f)

    def _load_alias_to_id(self):
        with open(path.join(self._data_dir, 'alias_to_id.json'), 'r') as f:
            self._alias_to_id = load(f)

    def _load_id_to_chromosome(self):
        with open(path.join(self._data_dir, 'id_to_chromosome.json'), 'r') as f:
            self._id_to_chromosome = load(f)

    def _load_target_ch(self, ch):
        with open(path.join(self._data_dir, ch + '.json'), 'r') as f:
            self._ch = load(f)
