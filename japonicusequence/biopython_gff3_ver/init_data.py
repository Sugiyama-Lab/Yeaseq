from urllib import request
import gzip
import os
import re
import pickle


address_work = os.path.dirname(__file__)
address_refseq_data = os.path.join(address_work, 'refseq_data')
address_ftp = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/Schizosaccharomyces_japonicus/all_assembly_versions'


def make_dir_refseq_data():
    list_content_work_address = os.listdir(address_work)
    if 'refseq_data' not in list_content_work_address:
        os.mkdir(address_refseq_data)


def get_dir_refseq_data_list():
    list_content_refseq_data = os.listdir(address_refseq_data)
    return list_content_refseq_data


def search_database_upload():
    url_japonicus_version = address_ftp
    list_url_content = request.urlopen(url_japonicus_version).readlines()
    list_assembly_version = []
    for each_row in list_url_content:
        re_matched_assembly_version = re.match('.*\s*? (GCF_000149845\..*) ->.*', str(each_row), re.I)
        if re_matched_assembly_version:
            inner_str_assembly_version = re_matched_assembly_version.group(1)
            list_assembly_version.append(inner_str_assembly_version)
    return list_assembly_version[-1]


def download_assembly(inner_str_assembly_version):
    url_file_list = address_ftp + '/{}'.format(inner_str_assembly_version)
    list_url_file_content = request.urlopen(url_file_list).readlines()
    for each_file_in_url in list_url_file_content:
        if '.gff.gz' in str(each_file_in_url):
            str_file_gz = str(each_file_in_url).split(' ')[-1][:-5]
    try:
        url_gff3_gz = url_file_list + '/{}'.format(str_file_gz)
    except NameError:
        raise
    else:
        handle_gff3_gz = request.urlopen(url_gff3_gz)
        with open(os.path.join(address_refseq_data, str_file_gz), 'wb') as file_gff3_gz:
            file_gff3_gz.write(handle_gff3_gz.read())

        str_file_name_gff = os.path.splitext(str_file_gz)[0]
        handle_gff = gzip.GzipFile(str_file_gz)
        with open(os.path.join(address_refseq_data, str_file_name_gff), 'wb') as file_gff3:
            file_gff3.write(handle_gff.read())

    return str_file_name_gff


def divide_assembly():
    str_file_name_gff = download_assembly(search_database_upload())
    address_assembly_divided = os.path.join(address_refseq_data, os.path.splitext(str_file_name_gff)[0])

    if os.path.isdir(address_assembly_divided):
        if os.listdir(address_assembly_divided):
            return None

    if not os.path.isdir(address_assembly_divided):
        os.makedirs(address_assembly_divided)

    with open(os.path.join(address_refseq_data, str_file_name_gff), 'r') as handle_gff:
        list_gff_raw_data = handle_gff.readlines()

    list_gff_data_row = []
    for each_data_row in list_gff_raw_data:
        if each_data_row.startswith('#'):
            continue
        if not each_data_row:
            continue

        type_this_row = each_data_row.split('\t')[2]
        if type_this_row == 'region':
            continue

        if type_this_row == 'gene':
            list_gff_data_row.append([each_data_row])
        else:
            list_gff_data_row[-1].append(each_data_row)

    dict_gff_data = dict()
    str_current_chromesome_accession = list_gff_data_row[0][0].split('\t')[0]
    for each_gene_data in list_gff_data_row:
        list_gene_info = each_gene_data[0].split('\t')
        str_this_chromesome_accession = list_gene_info[0]

        if str_this_chromesome_accession != str_current_chromesome_accession:
            with open(os.path.join(address_assembly_divided, str_current_chromesome_accession), 'wb') as handle_each_chromesome:
                pickle.dump(dict_gff_data, handle_each_chromesome)
            str_current_chromesome_accession = str_this_chromesome_accession
            dict_gff_data = dict()

        list_gene_description = list_gene_info[-1].split(';')
        for each_description in list_gene_description:
            if 'name=' in each_description.lower():
                str_genename = each_description.lower().replace('name=', '')
                continue
            re_matched_geneid = re.match('.*GeneID:(\d+).*', each_description, re.I)
            if re_matched_geneid:
                str_geneid = re_matched_geneid.group(1)
                continue

        temp_dict_each_gene = dict()
        temp_dict_each_gene['genename'] = str_genename
        for each_gene_details in each_gene_data:
            list_gene_details = each_gene_details.split('\t')
            if list_gene_details[2] not in temp_dict_each_gene.keys():
                temp_dict_each_gene[list_gene_details[2]] = [(list_gene_details[3], list_gene_details[4])]
            else:
                temp_dict_each_gene[list_gene_details[2]].append((list_gene_details[3], list_gene_details[4]))

        dict_gff_data[str_geneid] = temp_dict_each_gene

