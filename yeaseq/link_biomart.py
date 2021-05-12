from requests import get as requests_get

from .yeaseq_config import BiomartQuery


def query_biomart(stable_gene_id, query_yeast='sjaponicus_eg_gene', up_stream: int = 1000, down_stream: int = 1000):
    biomart_query_url = BiomartQuery.get_biomart_query(
        query_db=query_yeast, ensembl_geneid=stable_gene_id, up_stream=up_stream, down_stream=down_stream)

    biomart_query_handle = requests_get(biomart_query_url)
    biomart_query_content = biomart_query_handle.text
    return biomart_query_content


def biomart_info_parser(c: str) -> list:
    """
    :param c: query result from func query_biomart
    :return: dict contains attributions like the table in /test/biomart_result.md and an additional attribution 'seq'
    """
    info_list = []
    transcript_info_list = c.lstrip('>').split('>')
    for each_transcript in transcript_info_list:
        split_content = each_transcript.split('\n')
        title = split_content[0]
        title_list = title.split('|')
        seq = ''.join(split_content[1:]).replace(' ', '')
        info_dict = dict(zip(BiomartQuery.biomart_query_attribute(), title_list))
        info_dict['seq'] = seq
        info_list.append(info_dict)
    return info_list
