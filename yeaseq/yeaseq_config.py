from os import path
from os import makedirs
from os import environ
import sys
from platform import system


class GUIConfig:
    SearchSizeLow = (640, 400)
    SearchSizeHigh = (640, 520)
    SeqSize = (1350, 800)


class WorkDir:

    BaseDir = path.abspath('.')
    # BaseDir = path.dirname(path.abspath(__file__))
    work_dir_name = 'yeaseq_refseq_data'

    UserSystem = system()
    if UserSystem == 'Windows' or UserSystem.lower() == 'windows':
        USERPROFILE = environ.get('USERPROFILE')
        if USERPROFILE is not None:
            UserDir = USERPROFILE
        else:
            UserDir = BaseDir
    elif UserSystem == 'Darwin' or UserSystem.lower() == 'darwin':
        HOMEPATH = environ.get('HOME')
        if HOMEPATH is not None:
            UserDir = HOMEPATH
        else:
            UserDir = BaseDir
    else:
        UserDir = environ.get('HOME')
        if UserDir is None:
            UserDir = BaseDir

    WorkDir = path.join(UserDir, work_dir_name)

    @staticmethod
    def get_data_dir(n: str):
        return path.join(WorkDir.WorkDir, Yeast.get_yeast_name(n))

    @staticmethod
    def get_log_path():
        log_path = path.join(WorkDir.WorkDir, 'log.txt')
        if not path.isdir(log_path):
            makedirs(log_path)
        return path.join(WorkDir.WorkDir, 'log.txt')

    @staticmethod
    def get_run_source_dir():
        if hasattr(sys, '_MEIPASS'):
            run_source_dir = path.abspath(sys._MEIPASS)
        else:
            run_source_dir = None
        return run_source_dir

    @staticmethod
    def get_external_data_dir():
        run_source_dir = WorkDir.get_run_source_dir()
        if run_source_dir is not None:
            ext_data_dir = path.join(run_source_dir, 'external_data')
        else:
            ext_data_dir = path.join('.', 'external_data')
        return ext_data_dir

    @staticmethod
    def get_external_data_path(file_name: str):
        ext_data_dir = WorkDir.get_external_data_dir()
        ext_data_path = path.join(ext_data_dir, file_name)
        if path.exists(ext_data_path):
            return ext_data_path
        else:
            return None


class QueryParam:
    MaxUpstream = 1000
    MaxDownstream = 1000


class Yeast:
    pombe = ['sp', 'pombe', 'spombe', 'schizosaccharomyces pombe', 'schizosaccharomyces_pombe']
    japonicus = ['sj', 'japonicus', 'sjaponicus', 'schizosaccharomyces japonicus', 'schizosaccharomyces_japonicus']
    cryophilus = ['sc', 'cryophilus', 'scryophilus', 'schizosaccharomyces cryophilus', 'schizosaccharomyces_cryophilus']
    octosporus = ['so', 'octosporus', 'soctosporus', 'schizosaccharomyces octosporus', 'schizosaccharomyces_octosporus']
    _yeast_name_list = [japonicus, cryophilus, octosporus, pombe]

    YeastList = ['japonicus', 'cryophilus', 'octosporus', 'pombe']
    YeastNCBIName = ['Schizosaccharomyces_japonicus', 'Schizosaccharomyces_cryophilus', 'Schizosaccharomyces_octosporus', 'Schizosaccharomyces_pombe']
    YeastFullName = ['Schizosaccharomyces japonicus', 'Schizosaccharomyces cryophilus', 'Schizosaccharomyces octosporus', 'Schizosaccharomyces pombe']

    YeastNameDict = dict(zip(YeastList, _yeast_name_list))
    YeastShort2Long = dict(zip(YeastList, YeastNCBIName))

    @staticmethod
    def get_yeast_name(n: str):
        n = n.lower()
        for each_yeast in Yeast.YeastList:
            if n in Yeast.YeastNameDict[each_yeast]:
                return each_yeast
        return None

    @staticmethod
    def get_yeast_longname(n: str):
        short_name = Yeast.get_yeast_name(n)
        if short_name:
            return Yeast.YeastShort2Long[short_name]
        else:
            return None


class NCBIFTP:
    YeastFTP = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi'
    latest_cat = 'latest_assembly_versions'
    all_cat = 'all_assembly_versions'

    @staticmethod
    def get_latest_ftp_url(n):
        yeast_name = Yeast.get_yeast_longname(n)
        return f'{NCBIFTP.YeastFTP}/{yeast_name}/{NCBIFTP.latest_cat}/'

    @staticmethod
    def get_all_ftp_url(n):
        yeast_name = Yeast.get_yeast_longname(n)
        return f'{NCBIFTP.YeastFTP}/{yeast_name}/{NCBIFTP.all_cat}/'

    @staticmethod
    def get_genomes_all_https_url(n):
        yeast_name = Yeast.get_yeast_name(n)
        return {
            'pombe': 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/945/GCF_000002945.1_ASM294v2/',
            'japonicus': 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/149/845/GCF_000149845.2_SJ5/',
            'cryophilus': 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/155/GCF_000004155.1_SCY4/',
            'octosporus': 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/150/505/GCF_000150505.1_SO6/',
        }[yeast_name]


class EntrezQuery:
    url_entrez_efetch = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gene&id={}&rettype=gb&retmode=xml'  # .format(entrez_gene_id)

    @staticmethod
    def fill_efetch(gene_id):
        return EntrezQuery.url_entrez_efetch.format(gene_id)


class BiomartQuery:

    biomary_query_xml = '''<?xml version="1.0" encoding="UTF-8"?>
    <!DOCTYPE Query>
    <Query  virtualSchemaName = "fungi_mart" formatter = "FASTA" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
                    
        <Dataset name = "{query_db}" interface = "default" >
            <Filter name = "upstream_flank" value = "{up_stream}"/>
            <Filter name = "downstream_flank" value = "{down_stream}"/>
            <Filter name = "ensembl_gene_id" value = "{ensembl_geneid}"/>
            <Attribute name = "transcript_exon_intron" />
            <Attribute name = "ensembl_gene_id" />
            <Attribute name = "ensembl_transcript_id" />
            <Attribute name = "external_gene_name" />
            <Attribute name = "description" />
            <Attribute name = "chromosome_name" />
            <Attribute name = "gene_biotype" />
            <Attribute name = "start_position" />
            <Attribute name = "end_position" />
            <Attribute name = "5_utr_start" />
            <Attribute name = "5_utr_end" />
            <Attribute name = "3_utr_start" />
            <Attribute name = "3_utr_end" />
            <Attribute name = "transcript_biotype" />
            <Attribute name = "transcript_start" />
            <Attribute name = "transcript_end" />
            <Attribute name = "transcription_start_site" />
            <Attribute name = "transcript_length" />
            <Attribute name = "cds_length" />
            <Attribute name = "cds_start" />
            <Attribute name = "cds_end" />
            <Attribute name = "exon_chrom_start" />
            <Attribute name = "exon_chrom_end" />
            <Attribute name = "cdna_coding_start" />
            <Attribute name = "cdna_coding_end" />
            <Attribute name = "genomic_coding_start" />
            <Attribute name = "genomic_coding_end" />
            <Attribute name = "rank" />
            <Attribute name = "strand" />
        </Dataset>
    </Query>'''  # .format(query_db, up_stream, down_stream, ensembl_geneid)

    url_biomart = 'http://fungi.ensembl.org/biomart/martservice?query={}'  # .format(biomary_query_xml)

    @staticmethod
    def get_query_db(n):
        fungi_name = Yeast.get_yeast_longname(n)
        split_name = fungi_name.split('_')
        db_ident = f'{split_name[0][0].lower()}{split_name[1]}_eg_gene'
        return db_ident

    @staticmethod
    def get_biomart_query(query_db: str, ensembl_geneid, up_stream: int = 500, down_stream: int = 500):
        if not query_db.endswith('gene'):
            query_db = BiomartQuery.get_query_db(query_db)
        biomart_query_xml = BiomartQuery.biomary_query_xml.format(query_db=query_db, up_stream=up_stream, down_stream=down_stream, ensembl_geneid=ensembl_geneid)
        return BiomartQuery.url_biomart.format(biomart_query_xml)

    @staticmethod
    def biomart_query_attribute():
        attri_list = '''<Attribute name = "ensembl_gene_id" />
            <Attribute name = "ensembl_transcript_id" />
            <Attribute name = "external_gene_name" />
            <Attribute name = "description" />
            <Attribute name = "chromosome_name" />
            <Attribute name = "gene_biotype" />
            <Attribute name = "start_position" />
            <Attribute name = "end_position" />
            <Attribute name = "5_utr_start" />
            <Attribute name = "5_utr_end" />
            <Attribute name = "3_utr_start" />
            <Attribute name = "3_utr_end" />
            <Attribute name = "transcript_biotype" />
            <Attribute name = "transcript_start" />
            <Attribute name = "transcript_end" />
            <Attribute name = "transcription_start_site" />
            <Attribute name = "transcript_length" />
            <Attribute name = "cds_length" />
            <Attribute name = "cds_start" />
            <Attribute name = "cds_end" />
            <Attribute name = "exon_chrom_start" />
            <Attribute name = "exon_chrom_end" />
            <Attribute name = "cdna_coding_start" />
            <Attribute name = "cdna_coding_end" />
            <Attribute name = "genomic_coding_start" />
            <Attribute name = "genomic_coding_end" />
            <Attribute name = "rank" />
            <Attribute name = "strand" />'''.replace('<Attribute name = "', '').replace('" />', '').replace('\t', '').replace(' ', '').split('\n')
        return attri_list


class BaseTrans:
    BaseReverse = {'A': 'T',
                   'C': 'G',
                   'T': 'A',
                   'G': 'C'}
    BaseStart = 'ATG'
    BaseEnd = ['TAA', 'TGA', 'TAG']
    BaseToAA = {'GCT': 'A',
                'GCC': 'A',
                'GCA': 'A',
                'GCG': 'A',
                'CGT': 'R',
                'CGC': 'R',
                'CGA': 'R',
                'CGG': 'R',
                'AGA': 'R',
                'AGG': 'R',
                'AAT': 'N',
                'AAC': 'N',
                'GAT': 'D',
                'GAC': 'D',
                'TGT': 'C',
                'TGC': 'C',
                'CAA': 'Q',
                'CAG': 'Q',
                'GAA': 'E',
                'GAG': 'E',
                'GGT': 'G',
                'GGC': 'G',
                'GGA': 'G',
                'GGG': 'G',
                'CAT': 'H',
                'CAC': 'H',
                'ATT': 'I',
                'ATC': 'I',
                'ATA': 'I',
                'TTA': 'L',
                'TTG': 'L',
                'CTT': 'L',
                'CTC': 'L',
                'CTA': 'L',
                'CTG': 'L',
                'AAA': 'K',
                'AAG': 'K',
                'ATG': 'M',
                'TTT': 'F',
                'TTC': 'F',
                'CCT': 'P',
                'CCC': 'P',
                'CCA': 'P',
                'CCG': 'P',
                'TCT': 'S',
                'TCC': 'S',
                'TCA': 'S',
                'TCG': 'S',
                'AGT': 'S',
                'AGC': 'S',
                'ACT': 'T',
                'ACC': 'T',
                'ACA': 'T',
                'ACG': 'T',
                'TGG': 'W',
                'TAT': 'Y',
                'TAC': 'Y',
                'GTT': 'V',
                'GTC': 'V',
                'GTA': 'V',
                'GTG': 'V'}
