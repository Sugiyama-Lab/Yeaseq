from .yeaseq_config import Yeast
from .yeaseq_config import QueryParam
from .yeaseq_config import GUIConfig
from .yeaseq_config import WorkDir
from .ncbi_gff import OfflineDataManager
from .ncbi_gff import OfflineDataLoader
from .seq_process import OnlineDataProcessor
from .seq_process import OfflineDataProcessor
from .seq_process import SeqProcessor
from .ortho_table import OrthoGenes
from . import link_biomart
from . import link_entrez
from . import __info__

import wx
import wx.richtext
# import webbrowser


class SearchFrame(wx.Frame):
    def __init__(self, *args, **kwargs):
        super(SearchFrame, self).__init__(*args, **kwargs)

        self._search_panel = wx.Panel(self, style=wx.TAB_TRAVERSAL | wx.CLIP_CHILDREN | wx.FULL_REPAINT_ON_RESIZE)

        self._create_menubar()

        # Stable gene ID search
        self._search_content = None
        self._search_yeast = Yeast.YeastList[0]
        self._offline = False
        self._offline_status_button_dict = {True: u'Search online↑',
                                            False: u'Search offline↓'}
        self._offline_status_bar_dict = {True: 'Offline search',
                                         False: 'Online search'}
        self._upstream = QueryParam.MaxUpstream
        self._downstream = QueryParam.MaxDownstream

        # Ortho gene tool
        self._ortho_find_on = False
        self._ortho_find_button_dict = {True: u'Close ortholog finder←',
                                        False: u'Find the orthologs→'}

        # To store show seq frame objects
        self._showseq_frame_list = []

        # Font
        self._font_boxname = wx.Font(15, wx.FONTFAMILY_SWISS, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_BOLD, False)
        self._font_static_text = wx.Font(15, wx.FONTFAMILY_SWISS, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL, False)
        self._font_search_content = wx.Font(15, wx.FONTFAMILY_SWISS, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_BOLD, False)
        self._font_listbox = wx.Font(15, wx.FONTFAMILY_SWISS, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_LIGHT, False)

        # Status bar
        self._status_bar = self._initial_status_bar()

        # Total sizer
        self.boxsizer_full = wx.BoxSizer(wx.HORIZONTAL)  # Contain common and external tool
        self.boxsizer_all = wx.BoxSizer(wx.VERTICAL)  # For common use

        """ Box sizer for search widget sizer and message buttons sizer"""
        search_message_boxsizer = wx.BoxSizer(wx.HORIZONTAL)

        # Box sizer for search widgets
        search_boxsizer = wx.BoxSizer(wx.VERTICAL)

        # Search widgets for ensembl and ncbi
        staticboxsizer_ensembl, self._yeast_ensembl, self._stable_id = self._initial_search_widget(
            'Ensembl ID', self._search_symbol_event)
        staticboxsizer_ncbi, self._yeast_ncbi, self._entrez_id = self._initial_search_widget(
            'Entrez gene ID', self._search_entrez_event)
        search_boxsizer.Add(staticboxsizer_ensembl, 0, wx.ALL, 10)
        search_boxsizer.Add(staticboxsizer_ncbi, 0, wx.ALL, 10)

        search_message_boxsizer.Add(search_boxsizer, 0, wx.RIGHT, 20)

        # Help and about buttons
        message_grid_sizer = self._initial_message_button()
        search_message_boxsizer.Add(message_grid_sizer, 0, wx.TOP, 20)

        self.boxsizer_all.Add(search_message_boxsizer, 0, wx.ALL, 10)

        """ Box sizer for yeast listbox and uptream/downstream """
        listbox_exbase_boxsizer = wx.BoxSizer(wx.HORIZONTAL)

        # Yeast name listbox
        staticboxsizer_listbox, self._yeast_listbox = self._initial_yeast_listbox()
        listbox_exbase_boxsizer.Add(staticboxsizer_listbox, 0, wx.ALL, 10)

        # Entra base num (upstream and downstream)
        staticboxsizer_stream, self._upstream_text, self._downstream_text = self._initial_stream_set()
        listbox_exbase_boxsizer.Add(staticboxsizer_stream, 0, wx.ALL, 10)

        self.boxsizer_all.Add(listbox_exbase_boxsizer, 0, wx.ALL, 10)

        """ Box sizer for offline switch button and ortho gene find button """
        offlinebutton_orthotool_gbsizer = wx.GridBagSizer(hgap=240, vgap=4)

        # Offline switch button
        self._offline_switch_button = self._initial_offline_switch()
        offlinebutton_orthotool_gbsizer.Add(self._offline_switch_button, pos=(0, 0), span=(1, 1), flag=wx.ALIGN_RIGHT | wx.ALIGN_CENTRE_VERTICAL)

        # Orthologous gene tool button
        self._ortho_find_button = self._init_ortho_find_button()
        offlinebutton_orthotool_gbsizer.Add(self._ortho_find_button, pos=(0, 1), span=(1, 1), flag=wx.ALIGN_RIGHT | wx.ALIGN_CENTRE_VERTICAL)

        self.boxsizer_all.Add(offlinebutton_orthotool_gbsizer, 0, wx.ALL, 10)

        # Offline part
        self._staticboxsizer_offline, self._download_one_button = self._initial_offline_part()
        self.boxsizer_all.Add(self._staticboxsizer_offline, 0, wx.ALL, 10)
        self.boxsizer_all.Hide(self._staticboxsizer_offline)

        self.boxsizer_full.Add(self.boxsizer_all, 0, wx.ALL, 10)

        # Ortho tool part
        self._staticboxsizer_ortho_tool, self.ortho_search_content, self.ortho_info_text = self._init_ortho_tool_part()
        self.boxsizer_full.Add(self._staticboxsizer_ortho_tool, 0, wx.ALL, 10)
        self.boxsizer_full.Hide(self._staticboxsizer_ortho_tool)

        # Register boxsizer
        self._search_panel.SetSizerAndFit(self.boxsizer_full)
        self.SetClientSize(self._search_panel.GetBestSize())
        self.boxsizer_full.SetSizeHints(self._search_panel)
        self.boxsizer_full.Layout()

        # Offline data manager
        self._offline_data_manager = OfflineDataManager()

        # Ortho tool
        self._ortho_find_tool = OrthoGenes()

    def _create_menubar(self):
        file_menu = wx.Menu()
        about_item = file_menu.Append(-1, '&About...', 'About this app')
        file_menu.AppendSeparator()
        help_item = file_menu.Append(-1, '&Help...\tCtrl-H', 'Help')
        file_menu.AppendSeparator()
        exit_item = file_menu.Append(wx.ID_EXIT)

        menubar = wx.MenuBar()
        menubar.Append(file_menu, '&File')

        self.SetMenuBar(menubar)
        self.Bind(wx.EVT_MENU, self._help_menu, help_item)
        self.Bind(wx.EVT_MENU, self._exit_menu, exit_item)
        self.Bind(wx.EVT_MENU, self._about_menu, about_item)

    def _help_menu(self, event):
        _ = event
        wx.MessageBox(__info__.HelpMessage, 'Help', wx.OK | wx.ICON_INFORMATION)

    def _exit_menu(self, event):
        _ = event
        self.Close(True)

    def _about_menu(self, event):
        _ = event
        wx.MessageBox(__info__.AboutMessage, 'About', wx.OK | wx.ICON_INFORMATION)

    def _initial_status_bar(self):
        status_bar = self.CreateStatusBar()
        status_bar.SetStatusText(self._offline_status_bar_dict[self._offline])
        return status_bar

    def _initial_search_widget(self, box_name, search_func):
        static_box = wx.StaticBox(self._search_panel, -1, label=box_name)
        static_box.SetFont(self._font_boxname)
        static_box_sizer = wx.StaticBoxSizer(static_box, wx.VERTICAL)
        grid_sizer = wx.GridBagSizer(hgap=10, vgap=10)

        search_button = wx.Button(self._search_panel, -1, '&Search')
        search_button.SetFont(self._font_static_text)
        search_button.Bind(wx.EVT_BUTTON, search_func)
        grid_sizer.Add(search_button, pos=(0, 0), span=(1, 1), flag=wx.ALIGN_RIGHT | wx.ALIGN_CENTRE_VERTICAL)

        _search_content = wx.TextCtrl(self._search_panel, -1, size=(200, 30.5), style=wx.TE_PROCESS_ENTER | wx.TE_HT_ON_TEXT)
        _search_content.SetFont(self._font_search_content)
        _search_content.Bind(wx.EVT_TEXT_ENTER, search_func)
        grid_sizer.Add(_search_content, pos=(0, 1), span=(1, 1), flag=wx.ALIGN_RIGHT | wx.ALIGN_CENTRE_VERTICAL)

        link_preposition = wx.StaticText(self._search_panel, -1, label='in')
        link_preposition.SetFont(self._font_static_text)
        grid_sizer.Add(link_preposition, pos=(0, 2), span=(1, 1), flag=wx.ALIGN_RIGHT | wx.ALIGN_CENTRE_VERTICAL)
        yeast_name = wx.StaticText(self._search_panel, -1, label=self._change_yeast_label())
        yeast_name.SetFont(self._font_static_text)
        grid_sizer.Add(yeast_name, pos=(0, 3), span=(1, 1), flag=wx.ALIGN_RIGHT | wx.ALIGN_CENTRE_VERTICAL)
        static_box_sizer.Add(grid_sizer, proportion=0, flag=wx.ALIGN_CENTRE_VERTICAL | wx.RIGHT, border=50)
        return static_box_sizer, yeast_name, _search_content

    def _change_yeast_label(self):
        return f'S. {self._search_yeast}'

    def _search_symbol_event(self, event):
        _ = event
        self._search_control(search_type='Ensembl')

    def _search_entrez_event(self, event):
        _ = event
        self._search_control(search_type='NCBI')

    def _search_control(self, search_type='Ensembl', ):
        if search_type == 'Ensembl':
            search_content_widget = self._stable_id
        else:
            search_content_widget = self._entrez_id

        search_content = search_content_widget.GetValue().replace(' ', '')
        yeast = self._search_yeast
        upstream = self._get_stream_num(self._upstream_text)
        downstream = self._get_stream_num(self._downstream_text)
        print('Search', search_content, 'in', yeast, f'by {search_type}')

        if self._offline:
            offline_data_loader = OfflineDataLoader(search_content, n=yeast, content_type=search_type)
            try:
                ch_data, gene_id = offline_data_loader.get_ch_data()
                if gene_id is None:
                    wx.MessageBox(f'Not found input gene {search_content} in {yeast}\n'
                                  f'If input gene is in NCBI database, please check offline data (download again).', 'Error', wx.OK | wx.ICON_ERROR)
                    return None
                if ch_data is None:
                    wx.MessageBox(f'No information of input gene, please check offline data and download again', 'Error', wx.OK | wx.ICON_ERROR)
                    return None
            except FileNotFoundError:
                wx.MessageBox('Offline data is incomplete. Download again.', 'Error', wx.OK | wx.ICON_ERROR)
                return None
            data_processor = OfflineDataProcessor(ch_data, gene_id, max_upstream=upstream, max_downstream=downstream)
        else:
            if search_type == 'Ensembl':
                gene_symbol = search_content
            else:
                gene_symbol = link_entrez.query_entrez(search_content)
            biomart_result = link_biomart.query_biomart(
                gene_symbol, query_yeast=yeast, up_stream=upstream, down_stream=downstream)
            biomart_info = link_biomart.biomart_info_parser(biomart_result)
            print(biomart_info)
            try:
                data_processor = OnlineDataProcessor(biomart_info, max_upstream=upstream, max_downstream=downstream)
            except KeyError:
                wx.MessageBox(f'Not find the input gene "{search_content}" in {yeast}', 'No gene found', wx.OK | wx.ICON_ERROR)
                return None

        self._showseq_frame_list.append(
            SeqFrame(None, pos=[int(_/2) for _ in self.GetPosition()],
                     title=f'Yeaseq-Sequence - {self._offline_status_bar_dict[self._offline]}: {search_content}',
                     size=GUIConfig.SeqSize))
        self._showseq_frame_list[-1].set_data_processor(data_processor)
        self._showseq_frame_list[-1].Show(True)

    def _initial_message_button(self):
        grid_sizer = wx.GridBagSizer(hgap=10, vgap=10)

        help_button = wx.Button(self._search_panel, -1, '&Help')
        help_button.SetFont(self._font_static_text)
        help_button.Bind(wx.EVT_BUTTON, self._help_menu)

        about_button = wx.Button(self._search_panel, -1, '&About')
        about_button.SetFont(self._font_static_text)
        about_button.Bind(wx.EVT_BUTTON, self._about_menu)

        grid_sizer.Add(help_button, pos=(0, 0), span=(1, 1), flag=wx.ALIGN_RIGHT | wx.ALIGN_CENTRE_VERTICAL)
        grid_sizer.Add(about_button, pos=(1, 0), span=(1, 1), flag=wx.ALIGN_RIGHT | wx.ALIGN_CENTRE_VERTICAL)
        return grid_sizer

    def _initial_yeast_listbox(self):
        staticbox = wx.StaticBox(self._search_panel, -1, label='Yeast species')
        staticbox.SetFont(self._font_boxname)
        staticbox_sizer = wx.StaticBoxSizer(staticbox, wx.VERTICAL)
        grid_sizer = wx.GridBagSizer(hgap=10, vgap=10)

        yeast_listbox = wx.ListBox(self._search_panel, -1, choices=Yeast.YeastFullName, style=wx.LB_SINGLE)
        yeast_listbox.SetFont(self._font_static_text)
        yeast_listbox.SetSelection(0)
        yeast_listbox.Bind(wx.EVT_LISTBOX, self._yeast_listbox_event)

        grid_sizer.Add(yeast_listbox, pos=(0, 0), span=(1, 1), flag=wx.ALIGN_RIGHT | wx.ALIGN_CENTRE_VERTICAL)
        staticbox_sizer.Add(grid_sizer, proportion=0, flag=wx.ALL, border=10)
        return staticbox_sizer, yeast_listbox

    def _yeast_listbox_event(self, event):
        _ = event
        self._search_yeast = self._yeast_listbox.GetStringSelection().split(' ')[1]
        self._yeast_ensembl.SetLabel(self._change_yeast_label())
        self._yeast_ncbi.SetLabel(self._change_yeast_label())
        self._download_one_button.SetLabel(f'Download data for {self._search_yeast}')

    def _initial_stream_set(self):
        staticbox = wx.StaticBox(self._search_panel, -1, label='Extra base number')
        staticbox.SetFont(self._font_boxname)
        staticbox_sizer = wx.StaticBoxSizer(staticbox, wx.VERTICAL)

        up_grid_sizer, up_num_text = self._single_stream_set('Max upstream: ', self._upstream)
        down_grid_sizer, down_num_text = self._single_stream_set('Max downstream: ', self._downstream)
        staticbox_sizer.Add(up_grid_sizer, proportion=0, flag=wx.ALL, border=10)
        staticbox_sizer.Add(down_grid_sizer, proportion=0, flag=wx.ALL, border=10)
        return staticbox_sizer, up_num_text, down_num_text

    def _single_stream_set(self, desc_text, max_num):
        grid_sizer = wx.GridBagSizer(hgap=10, vgap=10)
        desc_text = wx.StaticText(self._search_panel, -1, label=desc_text)
        desc_text.SetFont(self._font_static_text)
        grid_sizer.Add(desc_text, pos=(0, 0), span=(1, 1), flag=wx.ALIGN_RIGHT | wx.ALIGN_CENTRE_VERTICAL)

        num_text = wx.TextCtrl(self._search_panel, -1, size=(100, 25))
        num_text.SetFont(self._font_search_content)
        num_text.SetValue(str(max_num))
        # Here two lines below were used to change the format of up-/down-stream number to thousand-separated number
        # num_text.SetValue(format(max_num, ','))
        # num_text.Bind(wx.EVT_TEXT, self._num_thousand_format_event)
        grid_sizer.Add(num_text, pos=(0, 1), span=(1, 1), flag=wx.ALIGN_RIGHT | wx.ALIGN_CENTRE_VERTICAL)
        return grid_sizer, num_text

    @staticmethod
    def _get_stream_num(text_widget):
        stream_content = text_widget.GetValue()
        try:
            # int_stream = int(stream_content.replace(',', '').replace(' ', ''))
            int_stream = int(stream_content.replace(' ', ''))
            return int_stream
        except ValueError:
            print('Error when convert upstream or downstream number')
            wx.MessageBox('Error when convert upstream or downstream number')

    def _num_thousand_format_event(self, event):
        event_id = event.GetId()
        event_widget = self._search_panel.FindWindowById(event_id)
        str_num = event_widget.GetValue()
        int_num = int(str_num.replace(',', '').replace(' ', ''))
        event_widget.SetValue(format(int_num, ','))

    def _initial_offline_switch(self):
        offline_switch_button = wx.Button(self._search_panel, -1, label=self._offline_status_button_dict[self._offline])
        offline_switch_button.SetFont(self._font_static_text)
        offline_switch_button.Bind(wx.EVT_BUTTON, self._switch_offline)
        return offline_switch_button

    def _switch_offline(self, event):
        _ = event
        if self._offline:
            self._offline = False

            self._status_bar.SetFieldsCount(1)
            self._status_bar.SetStatusText(self._offline_status_bar_dict[self._offline])

            self.boxsizer_all.Hide(self._staticboxsizer_offline)
            self._offline_switch_button.SetLabel(self._offline_status_button_dict[self._offline])
        else:
            self._offline = True

            self._status_bar.SetFieldsCount(2)
            self._status_bar.SetStatusWidths([-1, -4])
            self._status_bar.SetStatusText(self._offline_status_bar_dict[self._offline], 0)
            self._check_offline_data()

            self.boxsizer_all.Show(self._staticboxsizer_offline)
            self._offline_switch_button.SetLabel(self._offline_status_button_dict[self._offline])
        print('Panel size: ', self._search_panel.GetBestSize())
        self.SetClientSize(self._search_panel.GetBestSize())
        self.boxsizer_all.Layout()
        # self.boxsizer_full.Layout()

    def _initial_offline_part(self):
        staticbox = wx.StaticBox(self._search_panel, -1, label='Offline data')
        staticbox.SetFont(self._font_boxname)
        staticbox_sizer = wx.StaticBoxSizer(staticbox, wx.VERTICAL)

        grid_sizer = wx.GridBagSizer(hgap=4, vgap=4)

        # 生成一个message box，更新status bar
        data_check_button = wx.Button(self._search_panel, -1, label='Check offline data status')
        data_check_button.SetFont(self._font_static_text)
        data_check_button.Bind(wx.EVT_BUTTON, self._check_offline_data_button_event)
        grid_sizer.Add(data_check_button, pos=(0, 0), span=(1, 1), flag=wx.ALIGN_RIGHT | wx.ALIGN_CENTRE_VERTICAL)

        # 下载所有文件，生成一个message box
        download_all_button = wx.Button(self._search_panel, -1, label='Download all')
        download_all_button.SetFont(self._font_static_text)
        download_all_button.Bind(wx.EVT_BUTTON, self._download_all)
        grid_sizer.Add(download_all_button, pos=(0, 1), span=(1, 1), flag=wx.ALIGN_RIGHT | wx.ALIGN_CENTRE_VERTICAL)

        staticbox_sizer.Add(grid_sizer, proportion=0, flag=wx.EXPAND, border=10)

        # 下载其中一个文件，生成一个message box
        download_one_button = wx.Button(self._search_panel, -1, label=f'Download data for {self._search_yeast}')
        download_one_button.SetFont(self._font_static_text)
        download_one_button.Bind(wx.EVT_BUTTON, self._download_one)
        staticbox_sizer.Add(download_one_button, proportion=0, flag=wx.EXPAND, border=10)
        return staticbox_sizer, download_one_button

    def _check_offline_data(self):
        exist_list = self._offline_data_manager.check_data_existence()
        exist_status = ['{}-{}'.format(each_yeast, 'OK' if exist_list[i] else 'None') for i, each_yeast in enumerate(Yeast.YeastList)]
        status_text = 'Data status: {}'.format(' | '.join(exist_status))
        self._status_bar.SetStatusText(status_text, 1)
        return exist_status

    def _check_offline_data_button_event(self, event):
        _ = event
        exist_status = self._check_offline_data()
        wx.MessageBox('''
Data existion status:
    {}

Data folder:
    {}'''.format('\n    '.join(exist_status), WorkDir.WorkDir))

    def _download_all(self, event):
        _ = event
        try:
            for n in set(Yeast.YeastList):
                print(f'Download {n} data')
                self._offline_data_manager.download_data(n)
                self._offline_data_manager.save_data(n)
            wx.MessageBox('Done')
            self._check_offline_data()
        except Exception as r:
            print(r)
            wx.MessageBox('Error')

    def _download_one(self, event):
        _ = event
        n = self._search_yeast
        try:
            # if n == 'pombe':
            #     wx.MessageBox('Offline data for pombe if not supported now')
            #     return None
            print(f'Download {n} data')
            self._offline_data_manager.download_data(n)
            self._offline_data_manager.save_data(n)
            wx.MessageBox('Done')
            self._check_offline_data()
        except TimeoutError as r:
            print(r)
            wx.MessageBox(f'Timeout error when downloading {n}')
        except Exception as r:
            print(r)
            wx.MessageBox(f'Error when downloading {n}')

    def _init_ortho_find_button(self):
        ortho_find_button = wx.Button(self._search_panel, -1, label=self._ortho_find_button_dict[self._ortho_find_on])
        ortho_find_button.SetFont(self._font_static_text)
        ortho_find_button.Bind(wx.EVT_BUTTON, self._switch_ortho_find_tool)
        return ortho_find_button

    def _switch_ortho_find_tool(self, event):
        _ = event
        if self._ortho_find_on:
            self._ortho_find_on = False

            self.boxsizer_full.Hide(self._staticboxsizer_ortho_tool)
            self._ortho_find_button.SetLabel(self._ortho_find_button_dict[self._ortho_find_on])
        else:
            self._ortho_find_on = True

            self.boxsizer_full.Show(self._staticboxsizer_ortho_tool)
            self._ortho_find_button.SetLabel(self._ortho_find_button_dict[self._ortho_find_on])
        print('Panel size: ', self._search_panel.GetBestSize())
        self.SetClientSize(self._search_panel.GetBestSize())
        self.boxsizer_full.Layout()

    def _init_ortho_tool_part(self):
        staticbox = wx.StaticBox(self._search_panel, -1, label='Orthologous gene finder')
        staticbox.SetFont(self._font_boxname)
        staticbox_sizer = wx.StaticBoxSizer(staticbox, wx.VERTICAL)

        grid_sizer = wx.GridBagSizer(hgap=4, vgap=4)

        """ BoxSizer for search button and text """
        search_boxsizer = wx.BoxSizer(wx.HORIZONTAL)

        # Search button
        search_ortho_button = wx.Button(self._search_panel, -1, label='Search orthologous information for')
        search_ortho_button.SetFont(self._font_static_text)
        search_ortho_button.Bind(wx.EVT_BUTTON, self._search_ortho_gene)
        search_boxsizer.Add(search_ortho_button, proportion=0, flag=wx.ALIGN_CENTRE_VERTICAL | wx.ALIGN_LEFT, border=50)

        # Search text
        search_content = wx.TextCtrl(self._search_panel, -1, size=(150, 30.5), style=wx.TE_PROCESS_ENTER | wx.TE_HT_ON_TEXT)
        search_content.SetFont(self._font_search_content)
        search_content.Bind(wx.EVT_TEXT_ENTER, self._search_ortho_gene)
        search_boxsizer.Add(search_content, proportion=0, flag=wx.ALIGN_CENTRE_VERTICAL | wx.RIGHT, border=50)

        grid_sizer.Add(search_boxsizer, pos=(0, 0), span=(1, 1), flag=wx.ALIGN_LEFT | wx.ALIGN_CENTRE_VERTICAL)

        """ Text to show result """
        ortho_info_text = wx.TextCtrl(self._search_panel, -1, size=(525, 325),
                                      value='', style=wx.TE_MULTILINE | wx.TE_READONLY | wx.TE_RICH)
        font_ortho_data_source = wx.Font(15, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL, False)
        style_ortho_data_source = wx.TextAttr('#FF34B3', colBack='#FFFFFF', font=font_ortho_data_source)
        ortho_info_text.SetDefaultStyle(style_ortho_data_source)
        ortho_info_text.AppendText(__info__.OrthoGeneInfo)

        grid_sizer.Add(ortho_info_text, pos=(1, 0), span=(1, 1), flag=wx.ALIGN_LEFT | wx.ALIGN_CENTRE_VERTICAL)

        staticbox_sizer.Add(grid_sizer, proportion=0, flag=wx.ALIGN_CENTRE_VERTICAL | wx.RIGHT, border=50)

        return staticbox_sizer, search_content, ortho_info_text

    def _search_ortho_gene(self, event):
        _ = event
        gene = self.ortho_search_content.GetValue().replace(' ', '')
        record = self._ortho_find_tool.find_gene(gene)
        if record is None:
            wx.MessageBox(f'Input "{gene}" is not found.', 'No gene found', wx.OK | wx.ICON_INFORMATION)
            return None
        texts = self._ortho_find_tool.fill_shown_text(gene, record)

        self.ortho_info_text.SetValue('')

        font_ortho_data_source = wx.Font(10, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL, False)
        style_ortho_data_source = wx.TextAttr('#FF34B3', colBack='#FFFFFF', font=font_ortho_data_source)
        self.ortho_info_text.SetDefaultStyle(style_ortho_data_source)
        self.ortho_info_text.AppendText(__info__.OrthoGeneInfo)

        font_ortho_info = wx.Font(15, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL, False)
        style_ortho_info = wx.TextAttr('#000000', colBack='#FFFFFF', font=font_ortho_info)
        self.ortho_info_text.SetDefaultStyle(style_ortho_info)
        self.ortho_info_text.AppendText(''.join(texts))


class SeqPanel(wx.Panel):
    def __init__(self, notebook, transcript_data, **kwargs):
        super(SeqPanel, self).__init__(notebook, **kwargs)

        self._trans_data = transcript_data
        self._seq_processor = SeqProcessor(self._trans_data)

        self._max_upstream, self._max_downstream = self._seq_processor.get_stream_num()
        self._future_list = self._seq_processor.get_seq_future()
        self._future_num = len(self._future_list)
        self._future_widget_dict = {}  # For storing seq display control widgets
        self._widget_rename_dict = {'UTR_5': '5\'UTR',
                                    'Exon': 'Exon',
                                    'Intron': 'Intron',
                                    'UTR_3': '3\'UTR'}
        self._seq_attr_dict = self._set_seq_attr()
        self._protein_seq = None

        # Font
        self._font_boxname = wx.Font(15, wx.FONTFAMILY_SWISS, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_BOLD, False)
        self._font_static_text = wx.Font(15, wx.FONTFAMILY_SWISS, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL, False)
        self._font_seq = wx.Font(15, wx.FONTFAMILY_SWISS, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_BOLD, False)

        """Whole box sizer"""
        # Total sizer for the sizer with info and operation and the sizer for sequence display
        self.boxsizer_all = wx.BoxSizer(wx.HORIZONTAL)

        """First col of the whole panel"""
        # Vertical sizer for information box and operation sizer
        left_boxsizer = wx.BoxSizer(wx.VERTICAL)

        # Add info box to first vertical sizer
        info_boxsizer = wx.BoxSizer(wx.VERTICAL)
        staticboxsizer_basic_info = self._initial_basic_information()
        info_boxsizer.Add(staticboxsizer_basic_info, 0, wx.ALL | wx.EXPAND, 10)
        left_boxsizer.Add(info_boxsizer, 0, wx.ALL | wx.EXPAND, 10)

        # Horizontal sizer for two operation widgets
        operation_boxsizer = wx.BoxSizer(wx.HORIZONTAL)

        # Nucleotide sequence operations
        staticboxsizer_operations = self._initial_operate_col()
        operation_boxsizer.Add(staticboxsizer_operations, 0, wx.ALL, 10)

        # Change sequence to nucleotide or amino acid
        staticboxsizer_seq_type = self._initial_seq_type()
        operation_boxsizer.Add(staticboxsizer_seq_type, 0, wx.ALL, 10)

        # Add operation sizer to the left sizer
        left_boxsizer.Add(operation_boxsizer, 0, wx.ALL, 10)

        # Add second col to total sizer
        self.boxsizer_all.Add(left_boxsizer, 0, wx.ALL, 10)

        """Second col of the whole panel for showing sequence"""
        staticboxsizer_seq, self._seq_length_text, self._seq_text = self._initial_seq_col()
        self.boxsizer_all.Add(staticboxsizer_seq, 0, wx.ALL | wx.EXPAND, 10)
        self._change_seq()

        # Register boxsizer
        self.SetSizerAndFit(self.boxsizer_all)
        self.boxsizer_all.SetSizeHints(self)
        self.boxsizer_all.Layout()

    @staticmethod
    def _set_seq_attr():
        font_dict = dict()
        font_upstream = wx.Font(16, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL, False)
        style_upstream = wx.TextAttr('#FF34B3', colBack='#FFFFFF', font=font_upstream)
        font_utr5 = wx.Font(16, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL, False)
        style_utr5 = wx.TextAttr('#87CEFF', colBack='#FFFFFF', font=font_utr5)
        font_exon = wx.Font(16, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL, False)
        style_exon = wx.TextAttr('#000000', colBack='#FFE1FF', font=font_exon)
        font_intron = wx.Font(16, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_ITALIC, wx.FONTWEIGHT_NORMAL, False)
        style_intron = wx.TextAttr('#000000', colBack='#EEEE00', font=font_intron)
        font_utr3 = wx.Font(16, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL, False)
        style_utr3 = wx.TextAttr('#DA70D6', colBack='#FFFFFF', font=font_utr3)
        font_downstream = wx.Font(16, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL, False)
        style_downstream = wx.TextAttr('#76EE00', colBack='#FFFFFF', font=font_downstream)

        font_protein = wx.Font(16, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL, False)
        style_protein = wx.TextAttr('#000000', colBack='#FFFFFF', font=font_protein)

        font_dict['upstream'] = style_upstream
        font_dict['UTR_5'] = style_utr5
        font_dict['Exon'] = style_exon
        font_dict['Intron'] = style_intron
        font_dict['UTR_3'] = style_utr3
        font_dict['downstream'] = style_downstream
        font_dict['protein'] = style_protein
        return font_dict

    def _initial_basic_information(self):
        staticbox = wx.StaticBox(self, -1, label='Search information')
        staticbox.SetFont(self._font_boxname)
        staticbox_sizer = wx.StaticBoxSizer(staticbox, wx.VERTICAL)

        grid_sizer = wx.GridBagSizer(vgap=10, hgap=10)
        _i = 0
        for desc_key, desc_value in self._seq_processor.get_description():
            key_text = wx.StaticText(self, -1, label=desc_key)
            key_text.SetFont(self._font_static_text)
            value_text = wx.StaticText(self, -1, label=desc_value)
            value_text.SetFont(self._font_static_text)

            grid_sizer.Add(key_text, pos=(_i, 0), span=(1, 1), flag=wx.ALIGN_LEFT | wx.ALIGN_CENTRE_VERTICAL)
            grid_sizer.Add(value_text, pos=(_i, 1), span=(1, 1), flag=wx.ALIGN_LEFT | wx.ALIGN_CENTRE_VERTICAL)
            _i += 1

            if desc_key == 'GeneDescription':
                # desc_value = ' '.join(([f'word{k}' for k in range(100)]))
                _desc_content = desc_value

                def _gene_desc_event(e):
                    _ = e
                    wx.MessageBox(_desc_content, 'Gene description', wx.OK | wx.ICON_INFORMATION)

                # TODO split this to single func and move the numbers to config
                desc_words = desc_value.split(' ')
                aligned_desc_lines = ['']
                desc_out_flag = False
                for word in desc_words:
                    len_this_word = len(word) + 1
                    if len(aligned_desc_lines) == 2:
                        if (len(aligned_desc_lines[-1]) + len_this_word + 3) > 45:
                            desc_out_flag = True
                            break
                        else:
                            pass
                    else:
                        if (len(aligned_desc_lines[-1]) + len_this_word) > 45:
                            aligned_desc_lines.append('')
                        else:
                            pass
                    aligned_desc_lines[-1] += word + ' '
                aligned_desc = '\n'.join(aligned_desc_lines)
                if desc_out_flag:
                    value_text.SetLabel(aligned_desc + '...')
                    gene_desc_button = wx.Button(self,  -1, '...')
                    gene_desc_button.SetFont(self._font_static_text)
                    gene_desc_button.Bind(wx.EVT_BUTTON, _gene_desc_event)
                    grid_sizer.Add(gene_desc_button, pos=(_i, 1), span=(1, 1), flag=wx.ALIGN_LEFT | wx.ALIGN_CENTRE_VERTICAL)
                else:
                    value_text.SetLabel(aligned_desc)
                _i += 1

        staticbox_sizer.Add(grid_sizer, proportion=0, flag=wx.EXPAND, border=10)
        return staticbox_sizer

    def _initial_operate_col(self):
        staticbox = wx.StaticBox(self, -1, label='Display options')
        staticbox.SetFont(self._font_boxname)
        staticbox_sizer = wx.StaticBoxSizer(staticbox, wx.VERTICAL)

        grid_sizer = wx.GridBagSizer(vgap=10, hgap=10)

        up_spinctrl, up_statictext = self._initial_single_stream_widget('upstream', self._max_upstream)
        grid_sizer.Add(up_spinctrl, pos=(0, 0), span=(1, 1), flag=wx.ALIGN_LEFT | wx.ALIGN_CENTRE_VERTICAL)
        grid_sizer.Add(up_statictext, pos=(0, 1), span=(1, 1), flag=wx.ALIGN_LEFT | wx.ALIGN_CENTRE_VERTICAL)

        for widget_num, widget_name in enumerate(self._future_list):
            checkbox = wx.CheckBox(self, -1, self._widget_rename_dict[widget_name])
            checkbox.SetValue(False)
            checkbox.SetFont(self._font_static_text)
            checkbox.Bind(wx.EVT_CHECKBOX, self._change_seq_event)
            self._future_widget_dict[widget_name] = checkbox

            grid_sizer.Add(checkbox, pos=(widget_num + 1, 0), span=(1, 1), flag=wx.ALIGN_LEFT | wx.ALIGN_CENTRE_VERTICAL)

        down_spinctrl, down_statictext = self._initial_single_stream_widget('downstream', self._max_downstream)
        down_x_pos = self._future_num + 1
        grid_sizer.Add(down_spinctrl, pos=(down_x_pos, 0), span=(1, 1), flag=wx.ALIGN_LEFT | wx.ALIGN_CENTRE_VERTICAL)
        grid_sizer.Add(down_statictext, pos=(down_x_pos, 1), span=(1, 1), flag=wx.ALIGN_LEFT | wx.ALIGN_CENTRE_VERTICAL)

        if 'Exon' in self._future_widget_dict:
            self._future_widget_dict['Exon'].SetValue(True)
        staticbox_sizer.Add(grid_sizer, proportion=0, flag=wx.EXPAND, border=10)
        return staticbox_sizer

    def _initial_single_stream_widget(self, stream_type, max_num):
        spinctrl = wx.SpinCtrl(self, -1)
        spinctrl.SetRange(0, max_num)
        spinctrl.SetValue(0)
        spinctrl.Bind(wx.EVT_SPINCTRL, self._change_seq_event)
        spinctrl.Bind(wx.EVT_TEXT, self._change_seq_event)

        spinctrl.SetFont(self._font_static_text)
        self._future_widget_dict[stream_type] = spinctrl
        spinctrl.GetRange()
        statictext = wx.StaticText(self, label=f'bases {stream_type} of UTR (0-{max_num})')
        statictext.SetFont(self._font_static_text)
        return spinctrl, statictext

    def _initial_seq_type(self):
        staticbox = wx.StaticBox(self, -1, label='Sequence type')
        staticbox.SetFont(self._font_boxname)
        staticbox_sizer = wx.StaticBoxSizer(staticbox, wx.VERTICAL)
        grid_sizer = wx.GridBagSizer(vgap=10, hgap=10)

        radio_box = []
        radio_box.append(wx.RadioButton(self, -1, label='Nucleotide', style=wx.RB_GROUP))
        radio_box.append(wx.RadioButton(self, -1, label='Amino acid'))
        for _i, _box in enumerate(radio_box):
            _box.SetFont(self._font_static_text)
            _box.Bind(wx.EVT_RADIOBUTTON, self._change_seq_type_event)
            grid_sizer.Add(_box, pos=(_i, 0), span=(1, 1), flag=wx.ALIGN_LEFT | wx.ALIGN_CENTRE_VERTICAL)
        radio_box[0].SetValue(True)

        staticbox_sizer.Add(grid_sizer, proportion=0, flag=wx.EXPAND, border=10)
        return staticbox_sizer

    def _change_seq_type_event(self, event):
        event_id = event.GetId()
        event_widget = self.FindWindowById(event_id)
        seq_type = event_widget.GetLabel()
        if seq_type == 'Nucleotide':
            self._change_seq()
        else:
            if not self._protein_seq:
                self._protein_seq = self._seq_processor.get_protein_seq()
            print(self._protein_seq)
            _protein_seq_len = len(self._protein_seq)
            self._seq_text.SetValue(self._protein_seq)
            self._seq_text.SetStyle(0, _protein_seq_len, self._seq_attr_dict['protein'])
            self._seq_length_text.SetLabel(str(_protein_seq_len))

    def _initial_seq_col(self):
        staticbox = wx.StaticBox(self, -1, label='Sequence')
        staticbox.SetFont(self._font_boxname)
        staticbox_sizer = wx.StaticBoxSizer(staticbox, wx.VERTICAL)

        grid_sizer = wx.GridBagSizer(vgap=10, hgap=10)
        seq_length_desc_text = wx.StaticText(self, label='Sequence length: ')
        seq_length_desc_text.SetFont(self._font_static_text)
        seq_length_value_text = wx.StaticText(self, label=str(0))
        seq_length_value_text.SetFont(self._font_static_text)
        grid_sizer.Add(seq_length_desc_text, pos=(0, 0), span=(1, 1), flag=wx.ALIGN_LEFT | wx.ALIGN_CENTRE_VERTICAL)
        grid_sizer.Add(seq_length_value_text, pos=(0, 1), span=(1, 1), flag=wx.ALIGN_LEFT | wx.ALIGN_CENTRE_VERTICAL)
        staticbox_sizer.Add(grid_sizer, proportion=0, flag=wx.EXPAND, border=10)

        seq_text_ctrl = wx.richtext.RichTextCtrl(self, size=(1500, 1500),
                                                 value='', style=wx.TE_MULTILINE | wx.TE_READONLY)
        staticbox_sizer.Add(seq_text_ctrl, proportion=0, flag=wx.EXPAND, border=10)
        return staticbox_sizer, seq_length_value_text, seq_text_ctrl

    def _change_seq_event(self, e):
        _ = e
        self._change_seq()

    def _change_seq(self):
        future_status = [(future_name, _each_widget.GetValue()
                          ) for future_name, _each_widget in self._future_widget_dict.items()]
        _seq_future_dict = dict(future_status)
        if not _seq_future_dict:
            _seq_future_dict = self._seq_processor.get_initial_future_dict()

        self._seq_processor.set_feature(_seq_future_dict)
        _seq, _site_info = self._seq_processor.get_seq_info()
        self._seq_length_text.SetLabel(str(len(_seq)))
        self._seq_text.SetValue(_seq)
        print(_site_info)
        for _site in _site_info:
            self._seq_text.SetStyle(*_site[1:], self._seq_attr_dict[_site[0]])


class SeqFrame(wx.Frame):
    def __init__(self, *args, **kwargs):
        super(SeqFrame, self).__init__(*args, **kwargs)
        self._data_processor = None

        self._notebook = wx.Notebook(self)
        self._transcript_name_list = None
        self._nb_dict = dict()

    def set_data_processor(self, data_processor):
        """
        接收 online 或 offline data precessor 类，初始化 frame
        """
        self._data_processor = data_processor
        self._transcript_name_list, trans_data = self._data_processor.get_transcript_data()
        for each_trans in self._transcript_name_list:
            _p = SeqPanel(self._notebook, trans_data[each_trans], style=wx.TAB_TRAVERSAL | wx.CLIP_CHILDREN | wx.FULL_REPAINT_ON_RESIZE)
            self._nb_dict[each_trans] = _p
            self._notebook.AddPage(_p, each_trans)

