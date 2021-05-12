import wx
from data_process import GetSeqInfo


class FrameInitiatingNow(wx.Frame):
    def __init__(self, *args, **kw):
        super(FrameInitiatingNow, self).__init__(*args, **kw)


class FrameSearch(wx.Frame):
    def __init__(self, *args, **kw):
        super(FrameSearch, self).__init__(*args, **kw)
        panel_search = wx.Panel(self)

        self.str_search_gene_name = ''
        self.str_search_gene_id = ''
        self.str_search_accession = ''
        self.str_search_organism = 'Schizosaccharomyces japonicus'

        statictext_title = wx.StaticText(panel_search, label='Japonicusequence', pos=(175, 10))
        font_title = statictext_title.GetFont()
        font_title.PointSize += 8
        statictext_title.SetFont(font_title)

        statictext_genename = wx.StaticText(panel_search, label='Gene name', pos=(10, 60))
        statictext_in_organism = wx.StaticText(panel_search, label='in', pos=(220, 100))
        statictext_geneid = wx.StaticText(panel_search, label='Gene id', pos=(10, 155))
        statictext_accession = wx.StaticText(panel_search, label='RNA accession', pos=(300, 155))

        font_title.PointSize -= 4
        statictext_genename.SetFont(font_title)
        statictext_in_organism.SetFont(font_title)
        statictext_geneid.SetFont(font_title)
        statictext_accession.SetFont(font_title)

        self.content_gene_name = wx.TextCtrl(panel_search, value=self.str_search_gene_name, size=(200, 40), pos=(10, 100), style=wx.TE_NOHIDESEL)
        self.content_organism = wx.TextCtrl(panel_search, value=self.str_search_organism, size=(315, 40), pos=(250, 100), style=wx.TE_NOHIDESEL)
        self.content_gene_id = wx.TextCtrl(panel_search, value=self.str_search_gene_id, size=(200, 40), pos=(10, 200), style=wx.TE_NOHIDESEL)
        self.content_accession = wx.TextCtrl(panel_search, value=self.str_search_accession, size=(200, 40), pos=(300, 200), style=wx.TE_NOHIDESEL)

        font_input_context = self.content_gene_name.GetFont()
        font_input_context.PointSize += 3
        self.content_gene_name.SetFont(font_input_context)
        self.content_organism.SetFont(font_input_context)
        self.content_gene_id.SetFont(font_input_context)
        self.content_accession.SetFont(font_input_context)

        button_search_gene_name = wx.Button(panel_search, -1, '&Search', size=(100, 30), pos=(135, 63))
        button_search_gene_id = wx.Button(panel_search, -1, '&Search', size=(100, 30), pos=(100, 155))
        button_search_accession = wx.Button(panel_search, -1, '&Search', size=(100, 30), pos=(465, 155))
        button_search_gene_name.Bind(wx.EVT_BUTTON, self.action_button_search_gene_name)
        button_search_gene_id.Bind(wx.EVT_BUTTON, self.action_button_search_gene_id)
        button_search_accession.Bind(wx.EVT_BUTTON, self.action_button_search_accession)

    # 打开一个searching now frame，调用process类，返回字典，关闭searching now frame，将字典全局化，打开显示frame
    def action_button_search_gene_name(self, event):
        getseqinfo = GetSeqInfo('geneid', '22831531')
        wx.MessageBox(str(getseqinfo.return_data()))

    def action_button_search_gene_id(self, event):
        pass

    def action_button_search_accession(self, event):
        pass


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


main_app = wx.App(False)
frame_search = FrameSearch(None, title='Japonicusequence', pos=(240, 160), size=(600, 300))
frame_search.Show()
main_app.MainLoop()
