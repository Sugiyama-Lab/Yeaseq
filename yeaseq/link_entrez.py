from requests import get as requests_get
from xml.etree import cElementTree as ElementTree

from .yeaseq_config import EntrezQuery


class XmlListConfig(list):
    def __init__(self, x_list):
        super(XmlListConfig, self).__init__()
        for element in x_list:
            if element:
                if len(element) == 1 or element[0].tag != element[1].tag:
                    self.append(XmlDictConfig(element))
                elif element[0].tag == element[1].tag:
                    self.append(XmlListConfig(element))
            elif element.text:
                text = element.text.strip()
                if text:
                    self.append(text)


class XmlDictConfig(dict):
    def __init__(self, parent_element):
        super(XmlDictConfig, self).__init__()
        if parent_element.items():
            self.update(dict(parent_element.items()))
        for element in parent_element:
            if element:
                if len(element) == 1 or element[0].tag != element[1].tag:
                    x_dict = XmlDictConfig(element)
                else:
                    x_dict = {element[0].tag: XmlListConfig(element)}
                if element.items():
                    x_dict.update(dict(element.items()))
                self.update({element.tag: x_dict})
            elif element.items():
                self.update({element.tag: dict(element.items())})
            else:
                self.update({element.tag: element.text})


def query_entrez(entrez_gene_id):
    """
    With a entrez gene id as input, a gene symbol or gene name as output.
    The gene symbol also indicates ensemble stable gene id.
    Example: input: '2539869', output: 'SPBC11B10.09'
    """
    url_entrez_query = EntrezQuery.fill_efetch(entrez_gene_id)
    handle_entrez_query = requests_get(url_entrez_query)
    entrez_query_content = handle_entrez_query.text
    entrez_query_root = ElementTree.XML(entrez_query_content)
    query_info_dict = XmlDictConfig(entrez_query_root)
    #   str_gene_location = dict_entrez_info['Entrezgene']['Entrezgene_gene-source']['Gene-source']['Gene-source_src-str1']   # located chromosome NC_003423
    gene_symbol = query_info_dict['Entrezgene']['Entrezgene_gene-source']['Gene-source']['Gene-source_src-str2']   # gene symbol SPBC582.03
    return gene_symbol
