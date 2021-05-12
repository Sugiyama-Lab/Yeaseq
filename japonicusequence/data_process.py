from Bio import Entrez
from Bio import SeqIO
from Bio.Entrez import Parser
import re
import os
import pickle


Entrez.email = 'lourh@shanghaitech.edu.cn'
address_work = os.path.dirname(__file__)
address_refseq_data = os.path.join(address_work, 'refseq_data')
for each_assembly_version in os.listdir(address_refseq_data):
    if os.path.isdir(os.path.join(address_refseq_data, each_assembly_version)):
        str_dir_assembly = each_assembly_version
address_assembly = os.path.join(address_refseq_data, str_dir_assembly)


class GetSeqInfo(object):
    #   每个函数加入没有搜索到条件返回None，wx调用每个搜索时，应先接收一次函数返回值
    #   for these values below:
    """
    self.str_description

    self.str_chromesome_accession
    self.str_accession_rna
    self.int_gene_id
    self.str_accession_protein

    self.seq_chromesome
    self.seq_gene
    self.seq_rna
    self.seq_protein

    self.int_start_site # start site of this gene in chromesome
    self.int_end_site # end site of this gene in chromesome
    self.int_cds_start
    self.int_cds_end
    """
    def __init__(self, search_type, search_content, search_organism='Schizosaccharomyces japonicus'):
        self.search_type = search_type
        self.search_content = search_content
        self.search_organism = search_organism

        if self.search_type == 'genename':
            self.int_gene_id = self.search_db_gene_name()  # search 没有结果，值为None
            self.str_chromesome_accession, self.int_start_site, self.int_end_site, self.str_accession_rna, self.str_rna_direction, self.seq_chromesome, self.seq_gene = self.search_db_gene_id()
            self.str_description, self.seq_rna, self.int_gene_id, self.str_accession_protein, self.seq_protein, self.int_cds_start, self.int_cds_end = self.search_db_nucleotide_accession()

        elif self.search_type == 'geneid':
            self.int_gene_id = self.search_content
            temp_list_csecg = self.search_db_gene_id()
            if temp_list_csecg:
                self.str_chromesome_accession, self.int_start_site, self.int_end_site, self.str_accession_rna, self.str_rna_direction, self.seq_chromesome, self.seq_gene = temp_list_csecg
            temp_list_drgascc = self.search_db_nucleotide_accession()
            if temp_list_drgascc:
                self.str_description, self.seq_rna, self.int_gene_id, self.str_accession_protein, self.seq_protein, self.int_cds_start, self.int_cds_end = temp_list_drgascc

        elif self.search_type == 'accession':
            self.str_accession_rna = self.search_content
            temp_list_drgascc = self.search_db_nucleotide_accession()
            if temp_list_drgascc:
                self.str_description, self.seq_rna, self.int_gene_id, self.str_accession_protein, self.seq_protein, self.int_cds_start, self.int_cds_end = temp_list_drgascc
            temp_list_csecg = self.search_db_gene_id()
            if temp_list_csecg:
                self.str_chromesome_accession, self.int_start_site, self.int_end_site, self.str_accession_rna, self.str_rna_direction, self.seq_chromesome, self.seq_gene = temp_list_csecg

        #   rna翻转为正链，基因组翻转
        if self.str_rna_direction == 'minus':
            self.seq_rna = self.seq_rna.reverse_complement()
        self.int_len_gene = len(self.seq_gene)
        self.int_len_rna = len(self.seq_rna)
        self.int_len_cds = self.int_cds_end - self.int_cds_start

        self.dict_raw_gene_info = self.read_sequenced_file()
        self.dict_raw_gene_info['intron'] = self.find_intron_site()
        self.dict_raw_gene_info['utr5'], self.dict_raw_gene_info['utr3'] = self.find_utr_site()

    def search_db_gene_name(self):
        search_string = '{}'.format(self.search_content) + '[Gene] AND ' + '{}'.format(self.search_organism) + '[Organism]'
        handle_gene_name = Entrez.esearch(db='gene', term=search_string, rettype='xml')
        record_handle_gene_name = Entrez.read(handle_gene_name)
        try:
            gene_id = int(record_handle_gene_name['IdList'][0])
        except Exception:
            raise
        return gene_id

    def search_db_gene_id(self):
        handle_gene_id = Entrez.efetch(db='gene', id=self.int_gene_id, rettype='gb', retmode='xml')
        try:
            record_handle_gene_id = Entrez.read(handle_gene_id)
        except Entrez.Parser.ValidationError:
            raise 'Not find the gene id {}'.format(self.int_gene_id)
        str_chromesome_accession = record_handle_gene_id[0]['Entrezgene_gene-source']['Gene-source']['Gene-source_src-str1']
        int_start_site = int(record_handle_gene_id[0]['Entrezgene_locus'][0]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval']['Seq-interval_from'])
        int_end_site = int(record_handle_gene_id[0]['Entrezgene_locus'][0]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval']['Seq-interval_to'])
        str_accession_rna = record_handle_gene_id[0]['Entrezgene_locus'][0]['Gene-commentary_products'][0]['Gene-commentary_accession']
        str_rna_direction = record_handle_gene_id[0]['Entrezgene_locus'][0]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval']['Seq-interval_strand']['Na-strand'].attributes['value']

        handle_chromesome_accession = Entrez.efetch(db='nucleotide', id=str_chromesome_accession, rettype='fasta', retmode='text')
        record_handle_chromesome_accession = SeqIO.read(handle_chromesome_accession, 'fasta')
        seq_chromesome = record_handle_chromesome_accession.seq
        seq_gene = seq_chromesome[int(int_start_site):int(int_end_site)+1]

        return [str_chromesome_accession, int_start_site, int_end_site, str_accession_rna, str_rna_direction, seq_chromesome, seq_gene]

    def search_db_nucleotide_accession(self):
        handle_gene_accession = Entrez.efetch(db='nucleotide', id=self.str_accession_rna, rettype='gb', retmode='text')
        try:
            record_handle_rna_accession = SeqIO.read(handle_gene_accession, 'genbank')
        except Entrez.Parser.ValidationError:
            raise 'Not find the rna accession {}'.format(self.str_accession_rna)
        str_description = record_handle_rna_accession.description
        seq_rna = record_handle_rna_accession.seq
        for each_feature in record_handle_rna_accession.features:
            if each_feature.type == 'gene':
                feature_db_xref = each_feature.qualifiers['db_xref']  # # ['FLYBASE:FBgn0039044', 'GeneID:2768677']
                for each_db_xref in feature_db_xref:
                    re_matched_geneid = re.match('GeneID:(\d+)', each_db_xref, re.I)
                    if re_matched_geneid:
                        int_gene_id = int(re_matched_geneid.group(1))
            if each_feature.type == 'CDS':
                cds_from_to = each_feature.location
                str_accession_protein = each_feature.qualifiers['protein_id']
                seq_protein = each_feature.qualifiers['translation']
                re_matched_cds = re.match('\[(\d+):(\d+)]', str(cds_from_to), re.I)
                int_cds_start = int(re_matched_cds.group(1))
                int_cds_end = int(re_matched_cds.group(2))
        return [str_description, seq_rna, int_gene_id, str_accession_protein, seq_protein, int_cds_start, int_cds_end]

    def read_sequenced_file(self):
        for each_data_file in os.listdir(address_assembly):
            if self.str_chromesome_accession in each_data_file:
                str_filename_ref_data = each_data_file
        with open(os.path.join(address_assembly, str_filename_ref_data), 'rb') as handle_ref_data:
            dict_ref_data = pickle.load(handle_ref_data)
        dict_gene_info = dict_ref_data[str(self.int_gene_id)]
        return dict_gene_info

    def find_intron_site(self):
        if len(self.dict_raw_gene_info['exon']) > 1:
            list_flatten_spliced_seq_site = sum(self.dict_raw_gene_info['exon'], ())[1:-1]
            list_intron_seq_site = list(zip([list_flatten_spliced_seq_site[m] for m in range(0, len(list_flatten_spliced_seq_site), 2)], [list_flatten_spliced_seq_site[n] for n in range(1, len(list_flatten_spliced_seq_site), 2)]))
        else:
            list_intron_seq_site = []
        return list_intron_seq_site

    def find_utr_site(self):
        list_utr5_site = []
        list_utr3_site = []
        if len(self.dict_raw_gene_info['exon']) == 1:
            if self.dict_raw_gene_info['exon'][0][0] != self.dict_raw_gene_info['CDS'][0][0]:
                list_utr5_site.append((self.dict_raw_gene_info['exon'][0][0], str(int(self.dict_raw_gene_info['CDS'][0][0])-1)))
            else:
                list_utr5_site = []
        elif len(self.dict_raw_gene_info['exon']) > 1:
            for each_exon in self.dict_raw_gene_info['exon']:
                if int(each_exon[1]) > int(self.dict_raw_gene_info['CDS'][0][0]):
                    list_utr5_site.append((each_exon[0], str(int(self.dict_raw_gene_info['CDS'][0][0])-1)))
                    break
                else:
                    list_utr5_site.append((each_exon[0], str(int(each_exon[1])+1)))

        if len(self.dict_raw_gene_info['exon']) == 1:
            if self.dict_raw_gene_info['exon'][0][1] != self.dict_raw_gene_info['CDS'][0][1]:
                list_utr3_site.append((self.dict_raw_gene_info['CDS'][0][1], str(int(self.dict_raw_gene_info['exon'][0][0])-1)))
            else:
                list_utr3_site = []
        elif len(self.dict_raw_gene_info['exon']) > 1:
            for each_exon in self.dict_raw_gene_info['exon']:
                if int(each_exon[1]) > int(self.dict_raw_gene_info['CDS'][-1][1]):
                    if not list_utr3_site:
                        list_utr3_site.append((str(int(self.dict_raw_gene_info['CDS'][-1][1])+1), each_exon[1]))
                    else:
                        list_utr3_site.append((each_exon[0], str(int(each_exon[1])+1)))

        return list_utr5_site, list_utr3_site

    def return_data(self):
        dict_data_return = dict()
        dict_data_return['geneid'] = self.int_gene_id
        dict_data_return['genename'] = self.dict_raw_gene_info['genename']
        dict_data_return['rna_accession'] = self.str_accession_rna
        dict_data_return['genome_accession'] = self.dict_raw_gene_info['genome_accession']
        dict_data_return['protein_accession'] = self.str_accession_protein

        dict_data_return['description'] = self.str_description
        dict_data_return['rna_direction'] = self.str_rna_direction

        dict_data_return['exon'] = self.dict_raw_gene_info['exon']
        dict_data_return['intron'] = self.dict_raw_gene_info['intron']
        dict_data_return['cds'] = self.dict_raw_gene_info['CDS']

        dict_data_return['seq_chromesome'] = self.seq_chromesome
        dict_data_return['seq_protein'] = self.seq_protein

        return dict_data_return
