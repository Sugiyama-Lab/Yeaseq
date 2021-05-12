#!/usr/bin/python3
import requests
import wx
import wx.richtext
import webbrowser
import os


class EnsemblQuery(object):
    def __init__(self, search_content=''):
        self.str_gene_name = search_content

        self.xml_query = '''<?xml version="1.0" encoding="UTF-8"?>
        <!DOCTYPE Query>
        <Query  virtualSchemaName = "fungi_mart" formatter = "FASTA" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >

        	<Dataset name = "sjaponicus_eg_gene" interface = "default" >
        		<Filter name = "upstream_flank" value = "100"/>
        		<Filter name = "downstream_flank" value = "100"/>
        		<Filter name = "ensembl_gene_id" value = "{}"/>
        		<Attribute name = "transcript_exon_intron" />
        		<Attribute name = "ensembl_gene_id" />
        		<Attribute name = "start_position" />
        		<Attribute name = "end_position" />
        		<Attribute name = "5_utr_start" />
        		<Attribute name = "5_utr_end" />
        		<Attribute name = "3_utr_start" />
        		<Attribute name = "3_utr_end" />
        		<Attribute name = "exon_chrom_start" />
        		<Attribute name = "exon_chrom_end" />
        		<Attribute name = "strand" />
        	</Dataset>
        </Query>'''.format(self.str_gene_name)

        self.url_biomart = 'http://fungi.ensembl.org/biomart/martservice?query={}'.format(self.xml_query)
        self.dict_info_raw = self.get_ensembl_data()
        self.dict_info_processed = self.raw_dict_process()

    def get_ensembl_data(self):
        handle_biomart_url = requests.request('GET', self.url_biomart)
        str_data_content = handle_biomart_url.text
        str_head = str_data_content.split('\n')[0]
        str_seq = ''.join(str_data_content.split('\n')[1:])

        list_info = []
        dict_info_raw = dict()
        print(str_data_content)
        print(str_head)
        print(str_seq)
        for each_info in str_head.lstrip('>').split('|'):
            list_info.append(each_info)
        dict_info_raw['genename'] = list_info[0]
        dict_info_raw['gene_start'] = list_info[1]
        dict_info_raw['gene_end'] = list_info[2]
        dict_info_raw['utr5_start'] = list_info[3]
        dict_info_raw['utr5_end'] = list_info[4]
        dict_info_raw['utr3_start'] = list_info[5]
        dict_info_raw['utr3_end'] = list_info[6]
        dict_info_raw['exon_start'] = list_info[7]
        dict_info_raw['exon_end'] = list_info[8]
        dict_info_raw['strand'] = list_info[9]
        dict_info_raw['seq'] = str_seq
        print(dict_info_raw)
        return dict_info_raw

    def raw_dict_process(self):
        dict_info_processed = dict()
        dict_info_processed['genename'] = self.dict_info_raw['genename']
        dict_info_processed['gene_site'] = (int(self.dict_info_raw['gene_start']), int(self.dict_info_raw['gene_end']))
        dict_info_processed['gene_length'] = dict_info_processed['gene_site'][1] - dict_info_processed['gene_site'][
            0] + 1
        dict_info_processed['seq'] = self.dict_info_raw['seq']

        if self.dict_info_raw['strand'] == '1':
            int_site_normalization = 100 - int(self.dict_info_raw['gene_start'])

            dict_info_processed['exon_site'] = []
            for each_exon_start_site, each_exon_end_site in zip(
                    sorted([int(s) for s in self.dict_info_raw['exon_start'].split(';')]),
                    sorted([int(e) for e in self.dict_info_raw['exon_end'].split(';')])):
                dict_info_processed['exon_site'].append(
                    (each_exon_start_site + int_site_normalization, each_exon_end_site + int_site_normalization + 1))

            dict_info_processed['utr5_site'] = []
            dict_info_processed['utr3_site'] = []
            if self.dict_info_raw['utr5_start'] and self.dict_info_raw['utr5_end']:
                int_normalized_utr5_start = int(self.dict_info_raw['utr5_start']) + int_site_normalization
                int_normalized_utr5_end = int(self.dict_info_raw['utr5_end']) + int_site_normalization + 1
            if self.dict_info_raw['utr3_start'] and self.dict_info_raw['utr3_end']:
                int_normalized_utr3_start = int(self.dict_info_raw['utr3_start']) + int_site_normalization
                int_normalized_utr3_end = int(self.dict_info_raw['utr3_end']) + 1

            for each_exon_site_tuple in dict_info_processed['exon_site']:
                if self.dict_info_raw['utr5_start'] and self.dict_info_raw['utr5_end']:
                    if each_exon_site_tuple[1] >= int_normalized_utr5_end:
                        dict_info_processed['utr5_site'].append((each_exon_site_tuple[0], int_normalized_utr5_end))
                        break
                    elif each_exon_site_tuple[1] < int_normalized_utr5_end:
                        dict_info_processed['utr5_site'].append(each_exon_site_tuple)
                else:
                    dict_info_processed['utr5_site'] = tuple()

            for each_exon_site_tuple in dict_info_processed['exon_site']:
                if self.dict_info_raw['utr3_start'] and self.dict_info_raw['utr3_end']:
                    if each_exon_site_tuple[1] > int_normalized_utr3_start:
                        if dict_info_processed['utr3_site']:
                            dict_info_processed['utr3_site'].append(each_exon_site_tuple)
                        elif not dict_info_processed['utr3_site']:
                            dict_info_processed['utr3_site'].append(
                                (int_normalized_utr3_start, each_exon_site_tuple[1]))
                    elif each_exon_site_tuple[1] <= int_normalized_utr3_start:
                        continue
                else:
                    dict_info_processed['utr3_site'] = tuple()

            dict_info_processed['cds_site'] = []
            if dict_info_processed['utr5_site'] and dict_info_processed['utr3_site']:
                for each_exon_site_tuple in dict_info_processed['exon_site']:
                    if dict_info_processed['utr5_site'][-1][1] >= each_exon_site_tuple[1]:
                        continue
                    elif dict_info_processed['utr5_site'][-1][1] < each_exon_site_tuple[1]:
                        if dict_info_processed['utr3_site'][0][0] <= each_exon_site_tuple[1]:
                            if dict_info_processed['cds_site']:
                                dict_info_processed['cds_site'].append(
                                    (each_exon_site_tuple[0], dict_info_processed['utr3_site'][0][0]))
                                break
                            elif not dict_info_processed['cds_site']:
                                dict_info_processed['cds_site'].append(
                                    (dict_info_processed['utr5_site'][-1][1], dict_info_processed['utr3_site'][0][0]))
                                break
                        elif dict_info_processed['utr3_site'][0][0] > each_exon_site_tuple[1]:
                            if dict_info_processed['cds_site']:
                                dict_info_processed['cds_site'].append(each_exon_site_tuple)
                            elif not dict_info_processed['cds_site']:
                                dict_info_processed['cds_site'].append(
                                    (dict_info_processed['utr5_site'][-1][1], each_exon_site_tuple[1]))
            elif not dict_info_processed['utr5_site'] and dict_info_processed['utr3_site']:
                for each_exon_site_tuple in dict_info_processed['exon_site']:
                    if dict_info_processed['utr3_site'][0][0] <= each_exon_site_tuple[1]:
                        if dict_info_processed['cds_site']:
                            dict_info_processed['cds_site'].append(
                                (each_exon_site_tuple[0], dict_info_processed['utr3_site'][0][0]))
                            break
                        elif not dict_info_processed['cds_site']:
                            dict_info_processed['cds_site'].append(
                                (each_exon_site_tuple[0], dict_info_processed['utr3_site'][0][0]))
                            break
                    elif dict_info_processed['utr3_site'][0][0] > each_exon_site_tuple[1]:
                        if dict_info_processed['cds_site']:
                            dict_info_processed['cds_site'].append(each_exon_site_tuple)
                        elif not dict_info_processed['cds_site']:
                            dict_info_processed['cds_site'].append((each_exon_site_tuple[0], each_exon_site_tuple[1]))
            elif dict_info_processed['utr5_site'] and not dict_info_processed['utr3_site']:
                for each_exon_site_tuple in dict_info_processed['exon_site']:
                    if dict_info_processed['utr5_site'][-1][1] >= each_exon_site_tuple[1]:
                        continue
                    elif dict_info_processed['utr5_site'][-1][1] < each_exon_site_tuple[1]:
                        if dict_info_processed['cds_site']:
                            dict_info_processed['cds_site'].append(each_exon_site_tuple)
                        elif not dict_info_processed['cds_site']:
                            dict_info_processed['cds_site'].append(
                                (dict_info_processed['utr5_site'][-1][1], each_exon_site_tuple[1]))
            elif not dict_info_processed['utr5_site'] and not dict_info_processed['utr3_site']:
                dict_info_processed['cds_site'] = dict_info_processed['exon_site']

            if len(dict_info_processed['exon_site']) == 1:
                dict_info_processed['intron_site'] = []
            else:
                list_flatten_exon_site = sum(dict_info_processed['exon_site'], ())[1:-1]
                dict_info_processed['intron_site'] = list(
                    zip([list_flatten_exon_site[s] for s in range(0, len(list_flatten_exon_site), 2)],
                        [list_flatten_exon_site[e] for e in range(1, len(list_flatten_exon_site), 2)]))

            dict_info_processed['strand'] = 'plus'

        elif self.dict_info_raw['strand'] == '-1':
            int_site_normalization = int(self.dict_info_raw['gene_end']) + 100

            dict_info_processed['exon_site'] = []
            for each_exon_start_site, each_exon_end_site in zip(
                    sorted([int(e) for e in self.dict_info_raw['exon_end'].split(';')], reverse=True),
                    sorted([int(s) for s in self.dict_info_raw['exon_start'].split(';')], reverse=True)):
                dict_info_processed['exon_site'].append(
                    (int_site_normalization - each_exon_start_site, int_site_normalization - each_exon_end_site + 1))

            dict_info_processed['utr5_site'] = []
            dict_info_processed['utr3_site'] = []
            if self.dict_info_raw['utr5_start'] and self.dict_info_raw['utr5_end']:
                int_normalized_utr5_start = int_site_normalization - int(self.dict_info_raw['utr5_end'])
                int_normalized_utr5_end = int_site_normalization - int(self.dict_info_raw['utr5_start']) + 1
            if self.dict_info_raw['utr3_start'] and self.dict_info_raw['utr3_end']:
                int_normalized_utr3_start = int_site_normalization - int(self.dict_info_raw['utr3_end'])
                int_normalized_utr3_end = int_site_normalization - int(self.dict_info_raw['utr3_start'])

            for each_exon_site_tuple in dict_info_processed['exon_site']:
                if self.dict_info_raw['utr5_start'] and self.dict_info_raw['utr5_end']:
                    if each_exon_site_tuple[1] >= int_normalized_utr5_end:
                        dict_info_processed['utr5_site'].append((each_exon_site_tuple[0], int_normalized_utr5_end))
                        break
                    elif each_exon_site_tuple[1] < int_normalized_utr5_end:
                        dict_info_processed['utr5_site'].append(each_exon_site_tuple)
                else:
                    dict_info_processed['utr5_site'] = tuple()

            for each_exon_site_tuple in dict_info_processed['exon_site']:
                if self.dict_info_raw['utr3_start'] and self.dict_info_raw['utr3_end']:
                    if each_exon_site_tuple[1] > int_normalized_utr3_start:
                        if dict_info_processed['utr3_site']:
                            dict_info_processed['utr3_site'].append(each_exon_site_tuple)
                        elif not dict_info_processed['utr3_site']:
                            dict_info_processed['utr3_site'].append((int_normalized_utr3_start,
                                                                     each_exon_site_tuple[1]))
                    elif each_exon_site_tuple[1] <= int_normalized_utr3_start:
                        continue
                else:
                    dict_info_processed['utr3_site'] = tuple()

            dict_info_processed['cds_site'] = []
            if dict_info_processed['utr5_site'] and dict_info_processed['utr3_site']:
                for each_exon_site_tuple in dict_info_processed['exon_site']:
                    if dict_info_processed['utr5_site'][-1][1] >= each_exon_site_tuple[1]:
                        continue
                    elif dict_info_processed['utr5_site'][-1][1] < each_exon_site_tuple[1]:
                        if dict_info_processed['utr3_site'][0][0] <= each_exon_site_tuple[1]:
                            if dict_info_processed['cds_site']:
                                dict_info_processed['cds_site'].append(
                                    (each_exon_site_tuple[0], dict_info_processed['utr3_site'][0][0]))
                                break
                            elif not dict_info_processed['cds_site']:
                                dict_info_processed['cds_site'].append(
                                    (dict_info_processed['utr5_site'][-1][1], dict_info_processed['utr3_site'][0][0]))
                                break
                        elif dict_info_processed['utr3_site'][0][0] > each_exon_site_tuple[1]:
                            if dict_info_processed['cds_site']:
                                dict_info_processed['cds_site'].append(each_exon_site_tuple)
                            elif not dict_info_processed['cds_site']:
                                dict_info_processed['cds_site'].append(
                                    (dict_info_processed['utr5_site'][-1][1], each_exon_site_tuple[1]))
            elif not dict_info_processed['utr5_site'] and dict_info_processed['utr3_site']:
                for each_exon_site_tuple in dict_info_processed['exon_site']:
                    if dict_info_processed['utr3_site'][0][0] <= each_exon_site_tuple[1]:
                        if dict_info_processed['cds_site']:
                            dict_info_processed['cds_site'].append(
                                (each_exon_site_tuple[0], dict_info_processed['utr3_site'][0][0]))
                            break
                        elif not dict_info_processed['cds_site']:
                            dict_info_processed['cds_site'].append(
                                (each_exon_site_tuple[0], dict_info_processed['utr3_site'][0][0]))
                            break
                    elif dict_info_processed['utr3_site'][0][0] > each_exon_site_tuple[1]:
                        if dict_info_processed['cds_site']:
                            dict_info_processed['cds_site'].append(each_exon_site_tuple)
                        elif not dict_info_processed['cds_site']:
                            dict_info_processed['cds_site'].append((each_exon_site_tuple[0], each_exon_site_tuple[1]))
            elif dict_info_processed['utr5_site'] and not dict_info_processed['utr3_site']:
                for each_exon_site_tuple in dict_info_processed['exon_site']:
                    if dict_info_processed['utr5_site'][-1][1] >= each_exon_site_tuple[1]:
                        continue
                    elif dict_info_processed['utr5_site'][-1][1] < each_exon_site_tuple[1]:
                        if dict_info_processed['cds_site']:
                            dict_info_processed['cds_site'].append(each_exon_site_tuple)
                        elif not dict_info_processed['cds_site']:
                            dict_info_processed['cds_site'].append(
                                (dict_info_processed['utr5_site'][-1][1], each_exon_site_tuple[1]))
            elif not dict_info_processed['utr5_site'] and not dict_info_processed['utr3_site']:
                dict_info_processed['cds_site'] = dict_info_processed['exon_site']

            if len(dict_info_processed['exon_site']) == 1:
                dict_info_processed['intron_site'] = []
            else:
                list_flatten_exon_site = sum(dict_info_processed['exon_site'], ())[1:-1]
                dict_info_processed['intron_site'] = list(
                    zip([list_flatten_exon_site[s] for s in range(0, len(list_flatten_exon_site), 2)],
                        [list_flatten_exon_site[e] for e in range(1, len(list_flatten_exon_site), 2)]))

            dict_info_processed['strand'] = 'minus'

        return dict_info_processed

    def return_data(self):
        return self.dict_info_processed


class FrameSearch(wx.Frame):
    def __init__(self, *args, **kw):
        super(FrameSearch, self).__init__(*args, **kw)
        panel_search = wx.Panel(self)

        self.str_search_gene_name = ''
        self.str_search_organism = 'Schizosaccharomyces japonicus'

        statictext_title = wx.StaticText(panel_search, label='Japonicusequence', pos=(175, 10))
        font_title = statictext_title.GetFont()
        font_title.PointSize += 8
        statictext_title.SetFont(font_title)

        statictext_in_organism = wx.StaticText(panel_search, label='in {}'.format(self.str_search_organism),
                                               pos=(298, 62))
        font_title.PointSize -= 4
        statictext_in_organism.SetFont(font_title)

        self.content_gene_name = wx.TextCtrl(panel_search, value=self.str_search_gene_name,
                                             size=(180, 31), pos=(112, 59), style=wx.TE_NOHIDESEL | wx.TE_PROCESS_ENTER)

        font_input_context = self.content_gene_name.GetFont()
        font_input_context.PointSize += 4
        self.content_gene_name.SetFont(font_input_context)

        self.content_gene_name.Bind(wx.EVT_TEXT_ENTER, self.button_search_gene_name)

        button_search_gene_name = wx.Button(panel_search, -1, '&Search', size=(100, 30), pos=(5, 60))
        button_search_gene_name.Bind(wx.EVT_BUTTON, self.button_search_gene_name)

    def button_search_gene_name(self, event):
        self.str_search_gene_name = self.content_gene_name.GetValue()
        try:
            getseqinfo = EnsemblQuery(self.str_search_gene_name)
            global dict_info
            dict_info = getseqinfo.return_data()
            frame_show_seq = FrameShowSeq(None, title='Japonicusequence', pos=(240, 160), size=(950, 700))
            self.Destroy()
            frame_show_seq.Show()
        except IndexError:
            wx.MessageBox('Not found the gene', 'Error', wx.OK)
        except requests.exceptions.ConnectionError:
            wx.MessageBox('Please check your network', 'Error', wx.OK)


class FrameSearchingNow(wx.Frame):
    def __init__(self, *args, **kw):
        super(FrameSearchingNow, self).__init__(*args, **kw)
        panel_searching_now = wx.Panel(self)
        self.count = 0
        self.gauge = wx.Gauge(panel_searching_now, -1, 50, (20, 50), (250, 25))
        self.gauge.SetBezelFace(3)
        self.gauge.SetShadowWidth(3)
        self.Bind(wx.EVT_IDLE, self.OnIdle)

    def OnIdle(self, event):
        self.count = self.count + 1
        if self.count == 100:
            self.count = 0
        self.gauge.SetValue(self.count)


class FrameShowSeq(wx.Frame):
    def __init__(self, *args, **kw):
        super(FrameShowSeq, self).__init__(*args, **kw)
        panel_show = wx.Panel(self)

        self.dict_show_status = dict()

        statictext_part_description = wx.StaticText(panel_show, label='Display sequence of these gene parts: ',
                                                    pos=(10, 15))
        statictext_upstream = wx.StaticText(panel_show, label='bases upstream of UTR', pos=(115, 50))
        statictext_downstreamn = wx.StaticText(panel_show, label='bases downstream of UTR', pos=(115, 255))
        font_part_description = statictext_part_description.GetFont()
        font_part_description.PointSize += 2
        statictext_part_description.SetFont(font_part_description)
        statictext_upstream.SetFont(font_part_description)
        statictext_downstreamn.SetFont(font_part_description)

        self.spinctrl_upstream = wx.SpinCtrl(panel_show, -1, size=(80, 25), pos=(20, 50))
        self.spinctrl_upstream.SetRange(0, 100)
        self.spinctrl_upstream.SetValue(0)
        self.dict_show_status['upstream'] = 0
        self.spinctrl_downstream = wx.SpinCtrl(panel_show, -1, size=(80, 25), pos=(20, 255))
        self.spinctrl_downstream.SetRange(0, 100)
        self.spinctrl_downstream.SetValue(0)
        self.dict_show_status['downstream'] = 0

        self.checkbox_translation = wx.CheckBox(panel_show, -1, 'Show translation', size=(35, 40), pos=(1000, 1000))
        if dict_info['utr5_site']:
            self.checkbox_utr5 = wx.CheckBox(panel_show, -1, '5\'UTR', size=(100, 50), pos=(20, 85))
            self.dict_show_status['utr5_site'] = False
        else:
            statictext_no_utr5 = wx.StaticText(panel_show, label='(No 5\'UTR)', pos=(20, 85))
        self.checkbox_exon = wx.CheckBox(panel_show, -1, 'Exons', size=(100, 50), pos=(20, 120))
        self.dict_show_status['cds_site'] = True
        self.checkbox_exon.SetValue(True)
        if dict_info['intron_site']:
            self.checkbox_intron = wx.CheckBox(panel_show, -1, 'Introns', size=(100, 50), pos=(20, 155))
            self.dict_show_status['intron_site'] = False
        else:
            statictext_no_intron = wx.StaticText(panel_show, label='(No intron)', pos=(20, 175))
        if dict_info['utr3_site']:
            self.checkbox_utr3 = wx.CheckBox(panel_show, -1, '3\'UTR', size=(100, 50), pos=(20, 195))
            self.dict_show_status['utr3_site'] = False
        else:
            statictext_no_utr3 = wx.StaticText(panel_show, label='(No 3\'UTR)', pos=(20, 210))

        font_upstream = wx.Font(15, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL, False)
        self.style_upstream = wx.TextAttr('#FF34B3', colBack='#FFFFFF', font=font_upstream)
        font_utr5 = wx.Font(15, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL, False)
        self.style_utr5 = wx.TextAttr('#87CEFF', colBack='#FFFFFF', font=font_utr5)
        font_exon = wx.Font(15, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL, False)
        self.style_exon = wx.TextAttr('#000000', colBack='#FFE1FF', font=font_exon)
        font_intron = wx.Font(15, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_ITALIC, wx.FONTWEIGHT_NORMAL, False)
        self.style_intron = wx.TextAttr('#000000', colBack='#EEEE00', font=font_intron)
        font_utr3 = wx.Font(15, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL, False)
        self.style_utr3 = wx.TextAttr('#DA70D6', colBack='#FFFFFF', font=font_utr3)
        font_downstream = wx.Font(15, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL, False)
        self.style_downstream = wx.TextAttr('#76EE00', colBack='#FFFFFF', font=font_downstream)

        init_seq = ''.join(
            [dict_info['seq'][each_cds_site[0]:each_cds_site[1]] for each_cds_site in dict_info['cds_site']])
        str_header = '>' + dict_info['genename'] + ' length:{}'.format(len(init_seq)) + ' strand:{}'.format(
            dict_info['strand']) + '\n'
        self.content_seq = wx.richtext.RichTextCtrl(panel_show,
                                                    value=str_header + init_seq,
                                                    size=(500, 600), pos=(400, 10),
                                                    style=wx.TE_MULTILINE | wx.TE_READONLY)
        self.content_seq.SetStyle(len(str_header), len(str_header) + len(init_seq), self.style_exon)
        try:
            self.checkbox_utr5.Bind(wx.EVT_CHECKBOX, self.action_change_seq)
        except:
            pass
        self.checkbox_exon.Bind(wx.EVT_CHECKBOX, self.action_change_seq)
        try:
            self.checkbox_intron.Bind(wx.EVT_CHECKBOX, self.action_change_seq)
        except:
            pass
        try:
            self.checkbox_utr3.Bind(wx.EVT_CHECKBOX, self.action_change_seq)
        except:
            pass
        self.spinctrl_upstream.Bind(wx.EVT_SPINCTRL, self.action_change_seq)
        self.spinctrl_downstream.Bind(wx.EVT_SPINCTRL, self.action_change_seq)
        self.spinctrl_upstream.Bind(wx.EVT_TEXT, self.action_change_seq)
        self.spinctrl_downstream.Bind(wx.EVT_TEXT, self.action_change_seq)

        button_back = wx.Button(panel_show, -1, 'Back', size=(100, 30), pos=(30, 320))
        button_back.Bind(wx.EVT_BUTTON, self.action_button_back)
        button_save = wx.Button(panel_show, -1, 'Save sequence', size=(135, 30), pos=(180, 320))
        button_save.Bind(wx.EVT_BUTTON, self.action_button_save)
        button_ncbi_blastn = wx.Button(panel_show, -1, 'NCBI BLASTN', size=(120, 30), pos=(30, 420))
        button_ncbi_blastn.Bind(wx.EVT_BUTTON, self.action_button_ncbi_blastn)
        button_ensembl_dna_blast = wx.Button(panel_show, -1, 'Ensembl DNA BLAST', size=(180, 30), pos=(180, 420))
        button_ensembl_dna_blast.Bind(wx.EVT_BUTTON, self.action_button_ensembl_dna_blast)

    def action_button_back(self, event):
        self.Destroy()
        frame_search = FrameSearch(None, title='Japonicusequence', pos=(300, 240), size=(600, 150))
        frame_search.Show()

    def action_button_save(self, event):
        dlg = wx.FileDialog(self, 'Choose a file', os.getcwd(), '', 'All files (*.*)|*.*', wx.FD_SAVE)
        if dlg.ShowModal() == wx.ID_OK:
            try:
                with open(dlg.GetFilename(), 'x') as f:
                    f.write(self.content_seq.GetValue())
                wx.MessageBox('Save successfully', 'Confirm', wx.OK)
            except FileExistsError:
                wx.MessageBox('A file of the same name has already existed.', 'Error', wx.OK)
                self.action_button_save(event)
        dlg.Destroy()

    def action_button_ncbi_blastn(self, event):
        url_ncbi_blastn = 'https://blast.ncbi.nlm.nih.gov/Blast.cgi?QUERY={}&PROGRAM=blastn&PAGE_TYPE=BlastSearch'.format(
            self.content_seq.GetValue())
        webbrowser.open(url_ncbi_blastn)

    def action_button_ensembl_dna_blast(self, event):
        url_ensembl_dna_blast = 'https://fungi.ensembl.org/Schizosaccharomyces_pombe/Tools/Blast?{}'.format(
            self.content_seq.GetValue())
        webbrowser.open(url_ensembl_dna_blast)

    def action_get_parts_value(self, event):
        try:
            self.dict_show_status['utr5_site'] = self.checkbox_utr5.GetValue()
        except:
            pass
        try:
            self.dict_show_status['cds_site'] = self.checkbox_exon.GetValue()
        except:
            pass
        try:
            self.dict_show_status['intron_site'] = self.checkbox_intron.GetValue()
        except:
            pass
        try:
            self.dict_show_status['utr3_site'] = self.checkbox_utr3.GetValue()
        except:
            pass
        self.dict_show_status['upstream'] = int(self.spinctrl_upstream.GetValue())
        self.dict_show_status['downstream'] = int(self.spinctrl_downstream.GetValue())

    def action_change_seq(self, event):
        self.action_get_parts_value(event)
        list_available_part_site = []
        for each_part, status in self.dict_show_status.items():
            if (each_part != 'upstream') and (each_part != 'downstream') and status:
                if len(dict_info[each_part]) >= 1:
                    for each_each_part in dict_info[each_part]:
                        list_available_part_site.append((each_each_part, each_part))

        if self.dict_show_status['upstream'] != 0:
            list_available_part_site.append(((100 - self.dict_show_status['upstream'], 100), 'upstream'))
        if self.dict_show_status['downstream'] != 0:
            list_available_part_site.append(((dict_info['gene_length'] + 100,
                                              dict_info['gene_length'] + self.dict_show_status['downstream'] + 100),
                                             'downstream'))
        print(list_available_part_site)
        list_sorted_part_site = sorted(list_available_part_site, key=lambda l: l[0][0])
        str_seq_show = ''.join([dict_info['seq'][each_site_tuple[0][0]:each_site_tuple[0][1]] for each_site_tuple in
                                list_sorted_part_site])
        str_header = '>' + dict_info['genename'] + ' length:{}'.format(len(str_seq_show)) + ' strand:{}'.format(
            dict_info['strand']) + '\n'
        self.content_seq.SetValue(str_header + str_seq_show)

        list_part_len = []
        for each_part in list_sorted_part_site:
            list_part_len.append((each_part[0][1] - each_part[0][0], each_part[1]))

        int_init_site = len(str_header)
        int_end_site = int_init_site
        for each_len_part in list_part_len:
            int_end_site += each_len_part[0]
            if each_len_part[1] == 'upstream':
                self.content_seq.SetStyle(int_init_site, int_end_site, self.style_upstream)
            elif each_len_part[1] == 'utr5_site':
                self.content_seq.SetStyle(int_init_site, int_end_site, self.style_utr5)
            elif each_len_part[1] == 'cds_site':
                self.content_seq.SetStyle(int_init_site, int_end_site, self.style_exon)
            elif each_len_part[1] == 'intron_site':
                self.content_seq.SetStyle(int_init_site, int_end_site, self.style_intron)
            elif each_len_part[1] == 'utr3_site':
                self.content_seq.SetStyle(int_init_site, int_end_site, self.style_utr3)
            elif each_len_part[1] == 'downstream':
                self.content_seq.SetStyle(int_init_site, int_end_site, self.style_downstream)
            int_init_site = int_end_site


dict_info = dict()
main_app = wx.App(False)
frame_search = FrameSearch(None, title='Japonicusequence', pos=(300, 240), size=(700, 150))
frame_search.Show()
main_app.MainLoop()

