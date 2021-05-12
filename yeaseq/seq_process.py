from .yeaseq_config import BaseTrans


class SeqProcessor(object):
    """
    Manage one transcrpt in different transcript panel
    """
    def __init__(self, transcript_data):
        self._trans_data = transcript_data
        self._seq_dict = self._trans_data['seq']
        self._site_dict = self._trans_data['SeqSite']

        self._set_upstream = 0
        self._set_downstream = 0
        self._set_site_info = None

        self._initial_seq_future = {'upstream': 0,
                                    'UTR_5': False,
                                    'Exon': False,
                                    'Intron': False,
                                    'UTR_3': False,
                                    'downstream': 0
                                    }
        self._descrition_order = ['EnsemblGeneID', 'GeneSymbol', 'EnsemblTranscriptID',
                                  'GenbankTranscriptID', 'EntrezGeneID', 'ExGeneName', 'Name', 'GeneSynonym',
                                  'GeneDescription', 'ChromosomeName', 'GeneBiotype',
                                  'TranscriptType', 'GeneType', 'Strand']
        self._seq_future_order = ['UTR_5', 'Exon', 'Intron', 'UTR_3']

    def set_feature(self, seq_future: dict):  # TODO 注意传入的参数
        """
        {'upstream': 0,
        'UTR_5': False,
        'Exon': True,
        'Intron': False,
        'UTR_3': False,
        'downstream': 0
        }
        传入 future dict，其中部分 future 的 key 可能不存在

        得到排过序的 site [('exon', site, site), ('intron', site, site), ...]
        site 信息为不包含 upstream or downstream 且传入为 True 的四种信息
        """
        self._set_upstream = int(seq_future['upstream'])
        self._set_downstream = int(seq_future['downstream'])
        _future_on_list = [k for k, v in seq_future.items() if k not in ('upstream', 'downstream') and v]
        _show_sites = []
        for each_future in _future_on_list:
            _show_sites.extend(self._site_dict[each_future])
        if 'Intron' not in _future_on_list:
            _show_sites = [_ for _ in _show_sites if _[0] != 'Intron']
        sorted_sites = sorted(_show_sites, key=lambda x: x[1])
        self._set_site_info = sorted_sites

    def get_seq_info(self):
        """
        在 function set_future 中传入需要的 future 之后
        :return:
        """
        refined_site_info = []
        curr_length = self._set_upstream

        # upstream
        if self._set_upstream != 0:
            seq = self._seq_dict['upstream'][-self._set_upstream:]
            refined_site_info.append(('upstream', 0, self._set_upstream))
        else:
            seq = ''

        # Four or less other futures
        for each_sites in self._set_site_info:
            _type = each_sites[0]
            _start = each_sites[1]
            _end = each_sites[2]
            _length = _end - _start + 1
            seq += self._seq_dict['gene'][_start - 1: _end]
            refined_site_info.append((_type, curr_length, curr_length + _length))
            curr_length += _length

        # downstream
        seq += self._seq_dict['downstream'][:self._set_downstream]
        refined_site_info.append(('downstream', curr_length, curr_length + self._set_downstream))

        return seq, refined_site_info

    def get_protein_seq(self):
        exon_seq = ''
        for each_exon in self._site_dict['Exon']:
            exon_seq += self._seq_dict['gene'][each_exon[1] - 1: each_exon[2]]
        if exon_seq.startswith(BaseTrans.BaseStart):
            aa_seq = ''
            exon_length = len(exon_seq)
            aa_length, _ = divmod(exon_length, 3)
            if _ != 0:
                return 'Exon length is not triple multiple'
            for _i in range(aa_length):
                _each_3base = exon_seq[_i * 3: _i * 3 + 3]
                if _each_3base in BaseTrans.BaseEnd:
                    return aa_seq
                aa_seq += BaseTrans.BaseToAA[_each_3base]
            return 'No end code'
        else:
            return 'Exon doesnt start with start code ATG'

    def get_stream_num(self):
        return self._trans_data['max_upstream'], self._trans_data['max_downstream']

    def get_description(self):
        _description = self._trans_data['Description']
        description_keys = [_ for _ in self._descrition_order if _ in _description]
        for each_desc_key in description_keys:
            each_desc = _description[each_desc_key]
            if each_desc:
                yield each_desc_key, each_desc

    def get_seq_future(self):
        seq_futures = self._trans_data['SeqFuture']
        sorted_futures = [_ for _ in self._seq_future_order if _ in seq_futures]
        return sorted_futures

    def get_initial_future_dict(self):
        return self._initial_seq_future


class SeqDataProcessor(object):
    def __init__(self, max_upstream=1000, max_downstream=1000):
        self.max_upstream = max_upstream
        self.max_downstream = max_downstream

        self._transcript_list = []  # 总的 transcript name
        self._curr_transcript_name = None  # 当前显示的 transcript name
        self._transcript_num = 0  # transcript 的数量
        self._transcript_data_dict = dict()  # 预处理得到以 transcript name 为 key 的 data，将 online 和 offline 转换为统一的格式
        self._curr_raw_data = None
        '''
        以transcript为一个base进行site的操作
        先拿到序列，然后正链以 transcript_start 起始，负链以 transcript_end 起始，注意 ensembl 的 seq是已经序列反转过的（碱基和序列方向）
        '''
        self._curr_base_site = None
        self._strand = None
        self._strand_dict = {1: 'Forward',
                             -1: 'Minus'}
        self._base_reverse_dict = BaseTrans.BaseReverse

    def reverse_complement(self, seq):
        _reversed_seq = seq[::-1]
        _reverse_complement_seq = ''.join([self._base_reverse_dict[_] for _ in _reversed_seq])
        return _reverse_complement_seq

    def _normarlize_site(self, site):
        site = int(site)
        if self._strand == 1:
            norm_site = site - self._curr_base_site + 1
        else:
            norm_site = self._curr_base_site - site + 1
        return norm_site

    @staticmethod
    def _get_intron(site_list):
        if len(site_list[0]) == 3:
            expanded_sites = sum([list(_)[1:] for _ in site_list], [])
        else:
            expanded_sites = sum([list(_) for _ in site_list], [])
        intron_sites = [(expanded_sites[_] + 1, expanded_sites[_ + 1] - 1) for _ in range(1, len(expanded_sites) - 1, 2)]
        return intron_sites

    def get_transcript_data(self):
        return self._transcript_list, self._transcript_data_dict


class OnlineDataProcessor(SeqDataProcessor):
    def __init__(self, raw_data, **kwargs):
        super(OnlineDataProcessor, self).__init__(**kwargs)
        self._raw_data = raw_data
        self._seq_first_process()

    def _seq_first_process(self):
        for each_raw_data in self._raw_data:
            self._curr_raw_data = each_raw_data
            self._curr_transcript_name = self._curr_raw_data['ensembl_transcript_id']
            self._transcript_list.append(self._curr_transcript_name)
            self._transcript_num += 1
            self._transcript_data_dict[self._curr_transcript_name] = dict()

            self._fix_max_stream()

            self._save_normal_data()

            self._save_seq()

            self._process_raw_site_info()

    def _fix_max_stream(self):
        """
        接收查询时检查上下游长度是否符合要求
        :return:
        """
        self._strand = int(self._curr_raw_data['strand'])
        _trans_start = int(self._curr_raw_data['transcript_start'])
        _trans_end = int(self._curr_raw_data['transcript_end'])

        if self._strand == 1:
            if _trans_start <= self.max_upstream:
                self._transcript_data_dict[self._curr_transcript_name]['max_upstream'] = _trans_start - 1
            else:
                self._transcript_data_dict[self._curr_transcript_name]['max_upstream'] = self.max_upstream
            down_length = len(self._curr_raw_data['seq']) - (_trans_end - _trans_start + 1) - self._transcript_data_dict[self._curr_transcript_name]['max_upstream']
            if down_length < self.max_downstream:
                self._transcript_data_dict[self._curr_transcript_name]['max_downstream'] = down_length
            else:
                self._transcript_data_dict[self._curr_transcript_name]['max_downstream'] = self.max_downstream
            self._curr_base_site = _trans_start

        elif self._strand == -1:
            if _trans_start <= self.max_downstream:
                self._transcript_data_dict[self._curr_transcript_name]['max_downstream'] = _trans_start - 1
            else:
                self._transcript_data_dict[self._curr_transcript_name]['max_downstream'] = self.max_downstream
            up_length = len(self._curr_raw_data['seq']) - (_trans_end - _trans_start + 1) - self._transcript_data_dict[self._curr_transcript_name]['max_downstream']
            if up_length < self.max_upstream:
                self._transcript_data_dict[self._curr_transcript_name]['max_upstream'] = up_length
            else:
                self._transcript_data_dict[self._curr_transcript_name]['max_upstream'] = self.max_upstream
            self._curr_base_site = _trans_end

    def _save_normal_data(self):
        _description = dict()
        _description['EnsemblGeneID'] = self._curr_raw_data['ensembl_gene_id']
        _description['EnsemblTranscriptID'] = self._curr_raw_data['ensembl_transcript_id']
        _description['ExGeneName'] = self._curr_raw_data['external_gene_name']
        _description['GeneDescription'] = self._curr_raw_data['description']
        _description['ChromosomeName'] = self._curr_raw_data['chromosome_name']
        _description['GeneBiotype'] = self._curr_raw_data['gene_biotype']
        _description['TranscriptType'] = self._curr_raw_data['transcript_biotype']
        _description['Strand'] = self._strand_dict[self._strand]

        self._transcript_data_dict[self._curr_transcript_name]['Description'] = _description

    def _save_seq(self):
        _seq = dict()
        _max_upstream = self._transcript_data_dict[self._curr_transcript_name]['max_upstream']
        _max_downstream = self._transcript_data_dict[self._curr_transcript_name]['max_downstream']
        _seq['upstream'] = self._curr_raw_data['seq'][:_max_upstream]
        _seq['downstream'] = self._curr_raw_data['seq'][-_max_downstream:]
        if _max_downstream == 0:
            _seq['gene'] = self._curr_raw_data['seq'][_max_upstream:]
        else:
            _seq['gene'] = self._curr_raw_data['seq'][_max_upstream: -_max_downstream]
        self._transcript_data_dict[self._curr_transcript_name]['seq'] = _seq

    def _process_raw_site_info(self):
        site_dict = dict()
        seq_futures = ['Exon', ]

        def sort_merge_sites(start_site_content, end_site_content):
            start_site = [self._normarlize_site(_) for _ in start_site_content.split(';')]
            end_site = [self._normarlize_site(_) for _ in end_site_content.split(';')]
            merged_sites = [sorted(_) for _ in zip(start_site, end_site)]
            sorted_merged_sites = sorted(merged_sites, key=lambda x: x[0])
            return sorted_merged_sites

        # Exons
        sorted_exon_sites = sort_merge_sites(self._curr_raw_data['exon_chrom_start'], self._curr_raw_data['exon_chrom_end'])

        # Introns (may include UTR)
        intron_total_sites = self._get_intron(sorted_exon_sites)
        if intron_total_sites:
            seq_futures.append('Intron')

        # Have CDS site info
        if self._curr_raw_data['genomic_coding_start']:
            utr_5_sites = []
            utr_3_sites = []

            sorted_cds_sites = sort_merge_sites(self._curr_raw_data['genomic_coding_start'], self._curr_raw_data['genomic_coding_end'])

            # CDS
            site_dict['Exon'] = [('Exon', *_) for _ in sorted_cds_sites]

            # Intron in CDS
            cds_intron_sites = self._get_intron(sorted_cds_sites)
            if cds_intron_sites:
                site_dict['Intron'] = [('Intron', *_) for _ in cds_intron_sites]

            cds_start = sorted_cds_sites[0][0]
            cds_end = sorted_cds_sites[-1][-1]
            for each_exon_site in sorted_exon_sites:
                each_exon_start = each_exon_site[0]
                each_exon_end = each_exon_site[1]
                if cds_start > each_exon_start:
                    if cds_start > each_exon_end:
                        utr_5_sites.append(each_exon_site)
                    else:
                        utr_5_sites.append((each_exon_start, cds_start - 1))
                if cds_end < each_exon_end:
                    if cds_end > each_exon_start:
                        utr_3_sites.append((cds_end + 1, each_exon_end))
                    else:
                        utr_3_sites.append(each_exon_site)
            # 5'UTR
            if utr_5_sites:
                seq_futures.append('UTR_5')
                site_dict['UTR_5'] = [('UTR_5', *_) for _ in utr_5_sites]
                u5_intron_sites = self._get_intron(utr_5_sites)
                if u5_intron_sites:
                    site_dict['UTR_5'] += [('Intron', *_) for _ in u5_intron_sites]
            # 3'UTR
            if utr_3_sites:
                seq_futures.append('UTR_3')
                site_dict['UTR_3'] = [('UTR_3', *_) for _ in utr_3_sites]
                u3_intron_sites = self._get_intron(utr_3_sites)
                if u3_intron_sites:
                    site_dict['UTR_3'] += [('Intron', *_) for _ in u3_intron_sites]
        # No CDS info (maybe ncRNA)
        else:
            site_dict['Exon'] = [('Exon', *_) for _ in sorted_exon_sites]
            if intron_total_sites:
                site_dict['Intron'] = [('Intron', *_) for _ in intron_total_sites]

        self._transcript_data_dict[self._curr_transcript_name]['SeqSite'] = site_dict
        self._transcript_data_dict[self._curr_transcript_name]['SeqFuture'] = seq_futures


class OfflineDataProcessor(SeqDataProcessor):
    def __init__(self, ch_data, gene_id, **kwargs):
        super(OfflineDataProcessor, self).__init__(**kwargs)
        self._ch_data = ch_data
        self._gene_id = gene_id
        self._raw_data = ch_data[self._gene_id]
        self._seq_first_process()

    def get_gene_id(self):
        return self._gene_id

    def _seq_first_process(self):
        self._transcript_num = self._raw_data['TranscriptNumber']

        for trans_id, each_raw_data in self._raw_data['Transcript'].items():
            self._curr_raw_data = each_raw_data
            self._curr_transcript_name = self._curr_raw_data['RNA']['Genbank']
            self._transcript_list.append(self._curr_transcript_name)
            self._transcript_data_dict[self._curr_transcript_name] = dict()

            self._fix_max_stream()

            self._save_normal_data()

            self._save_seq()

            self._process_raw_site_info()

    def _fix_max_stream(self):
        _strand = self._raw_data['Gene']['Strand']
        if _strand in ('+', '-'):
            _strand = {'+': 1, '-': -1}[_strand]
        self._strand = int(_strand)
        _trans_start = int(self._curr_raw_data['RNA']['RNASite'][0])
        _trans_end = int(self._curr_raw_data['RNA']['RNASite'][1])
        _ch_end = int(self._ch_data['ChromosomeSite'][1])

        ch_end_length = _ch_end - _trans_end

        if self._strand == 1:
            if _trans_start <= self.max_upstream:
                self._transcript_data_dict[self._curr_transcript_name]['max_upstream'] = _trans_start - 1
            else:
                self._transcript_data_dict[self._curr_transcript_name]['max_upstream'] = self.max_upstream
            if ch_end_length < self.max_downstream:
                self._transcript_data_dict[self._curr_transcript_name]['max_downstream'] = ch_end_length
            else:
                self._transcript_data_dict[self._curr_transcript_name]['max_downstream'] = self.max_downstream
            self._curr_base_site = _trans_start

        elif self._strand == -1:
            if _trans_start <= self.max_downstream:
                self._transcript_data_dict[self._curr_transcript_name]['max_downstream'] = _trans_start - 1
            else:
                self._transcript_data_dict[self._curr_transcript_name]['max_downstream'] = self.max_downstream
            if ch_end_length < self.max_upstream:
                self._transcript_data_dict[self._curr_transcript_name]['max_upstream'] = ch_end_length
            else:
                self._transcript_data_dict[self._curr_transcript_name]['max_upstream'] = self.max_upstream
            self._curr_base_site = _trans_end

        else:
            return None

    def _save_normal_data(self):
        _description = dict()
        _description['ChromosomeName'] = self._ch_data['ChromosomeName']

        _description['GeneType'] = self._raw_data['Gene']['GeneType']
        _description['EntrezID'] = self._raw_data['Gene']['EntrezID']
        _description['GeneSymbol'] = self._raw_data['Gene']['Symbol']
        _description['GeneBiotype'] = self._raw_data['Gene']['GeneBiotype']
        _description['Name'] = self._raw_data['Gene']['Name']
        _description['GeneSynonym'] = self._raw_data['Gene']['Synonym']
        _description['Strand'] = self._strand_dict[self._strand]

        _description['TranscriptType'] = self._curr_raw_data['RNA']['RNAType']
        _description['GenbankTranscriptID'] = self._curr_raw_data['RNA']['Genbank']

        self._transcript_data_dict[self._curr_transcript_name]['Description'] = _description

    def _save_seq(self):
        _seq = dict()
        _max_upstream = self._transcript_data_dict[self._curr_transcript_name]['max_upstream']
        _max_downstream = self._transcript_data_dict[self._curr_transcript_name]['max_downstream']

        _trans_start = int(self._curr_raw_data['RNA']['RNASite'][0])
        _trans_end = int(self._curr_raw_data['RNA']['RNASite'][1])
        if self._strand == 1:
            _start_site = _trans_start - _max_upstream
            _end_site = _trans_end + _max_downstream
        else:
            _start_site = _trans_start - _max_downstream
            _end_site = _trans_end + _max_upstream

        _raw_seq = self._ch_data['Seq'][_start_site - 1: _end_site]
        if self._strand == -1:
            _raw_seq = self.reverse_complement(_raw_seq)
        _seq['upstream'] = _raw_seq[:_max_upstream]
        _seq['downstream'] = _raw_seq[-_max_downstream:]
        if _max_downstream == 0:
            _seq['gene'] = _raw_seq[_max_upstream:]
        else:
            _seq['gene'] = _raw_seq[_max_upstream: -_max_downstream]
        self._transcript_data_dict[self._curr_transcript_name]['seq'] = _seq

    def _process_raw_site_info(self):
        site_dict = dict()
        seq_futures = ['Exon', ]

        def sort_sites(sites):
            norm_sites = [sorted([self._normarlize_site(__) for __ in _]) for _ in sites]
            sorted_sites = sorted(norm_sites, key=lambda x: x[0])
            return sorted_sites
        # Exons
        sorted_exon_sites = sort_sites(self._curr_raw_data['Exon'])

        # Introns (may include UTR)
        intron_total_sites = self._get_intron(sorted_exon_sites)
        if intron_total_sites:
            seq_futures.append('Intron')

        # Have CDS site info
        if 'CDS' in self._curr_raw_data:
            utr_5_sites = []
            utr_3_sites = []

            sorted_cds_sites = sort_sites(self._curr_raw_data['CDS'])

            # CDS
            site_dict['Exon'] = [('Exon', *_) for _ in sorted_cds_sites]

            # Intron in CDS
            cds_intron_sites = self._get_intron(sorted_cds_sites)
            if cds_intron_sites:
                site_dict['Intron'] = [('Intron', *_) for _ in cds_intron_sites]

            cds_start = sorted_cds_sites[0][0]
            cds_end = sorted_cds_sites[-1][-1]
            for each_exon_site in sorted_exon_sites:
                each_exon_start = each_exon_site[0]
                each_exon_end = each_exon_site[1]
                if cds_start > each_exon_start:
                    if cds_start > each_exon_end:
                        utr_5_sites.append(each_exon_site)
                    else:
                        utr_5_sites.append((each_exon_start, cds_start - 1))
                if cds_end < each_exon_end:
                    if cds_end > each_exon_start:
                        utr_3_sites.append((cds_end + 1, each_exon_end))
                    else:
                        utr_3_sites.append(each_exon_site)
            # 5'UTR
            if utr_5_sites:
                seq_futures.append('UTR_5')
                site_dict['UTR_5'] = [('UTR_5', *_) for _ in utr_5_sites]
                u5_intron_sites = self._get_intron(utr_5_sites)
                if u5_intron_sites:
                    site_dict['UTR_5'] += [('Intron', *_) for _ in u5_intron_sites]
            # 3'UTR
            if utr_3_sites:
                seq_futures.append('UTR_3')
                site_dict['UTR_3'] = [('UTR_3', *_) for _ in utr_3_sites]
                u3_intron_sites = self._get_intron(utr_3_sites)
                if u3_intron_sites:
                    site_dict['UTR_3'] += [('Intron', *_) for _ in u3_intron_sites]
        # No CDS info (maybe ncRNA)
        else:
            site_dict['Exon'] = [('Exon', *_) for _ in sorted_exon_sites]
            if intron_total_sites:
                site_dict['Intron'] = [('Intron', *_) for _ in intron_total_sites]

        self._transcript_data_dict[self._curr_transcript_name]['SeqSite'] = site_dict
        self._transcript_data_dict[self._curr_transcript_name]['SeqFuture'] = seq_futures
