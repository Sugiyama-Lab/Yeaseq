import wx
from data_process import *


class MyFrame(wx.Frame):
    def __init__(self, parent=None, id=-1, title='Initiating'):
        wx.Frame.__init__(self, parent, id, title, style=wx.FRAME_SHAPED | wx.SIMPLE_BORDER)
        self.SetClientSize((300, 100))
        self.Center()
        self.txtShow = wx.StaticText(self, label='Initiating...', pos=(50, 10))
        font_title = self.txtShow.GetFont()
        font_title.PointSize += 5
        self.txtShow.SetFont(font_title)


class FrameInitiatingNow(wx.Frame):
    def __init__(self, *args, **kw):
        super(FrameInitiatingNow, self).__init__(*args, **kw)
        panel_initiating = wx.Panel(self)
        statictext_title = wx.StaticText(panel_initiating, label='Initiating...', pos=(50, 10))
        font_title = statictext_title.GetFont()
        font_title.PointSize += 5
        statictext_title.SetFont(font_title)


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
        statictext_in_organism = wx.StaticText(panel_search, label='in {}'.format(self.str_search_organism), pos=(220, 100))
        statictext_geneid = wx.StaticText(panel_search, label='Gene id', pos=(10, 155))
        statictext_accession = wx.StaticText(panel_search, label='RNA accession', pos=(300, 155))

        font_title.PointSize -= 4
        statictext_genename.SetFont(font_title)
        statictext_in_organism.SetFont(font_title)
        statictext_geneid.SetFont(font_title)
        statictext_accession.SetFont(font_title)

        self.content_gene_name = wx.TextCtrl(panel_search, value=self.str_search_gene_name, size=(200, 40), pos=(10, 100), style=wx.TE_NOHIDESEL)
        #   self.content_organism = wx.TextCtrl(panel_search, value=self.str_search_organism, size=(315, 40), pos=(250, 100), style=wx.TE_NOHIDESEL)
        self.content_gene_id = wx.TextCtrl(panel_search, value=self.str_search_gene_id, size=(200, 40), pos=(10, 200), style=wx.TE_NOHIDESEL)
        self.content_accession = wx.TextCtrl(panel_search, value=self.str_search_accession, size=(200, 40), pos=(300, 200), style=wx.TE_NOHIDESEL)

        font_input_context = self.content_gene_name.GetFont()
        font_input_context.PointSize += 3
        self.content_gene_name.SetFont(font_input_context)
        #   self.content_organism.SetFont(font_input_context)
        self.content_gene_id.SetFont(font_input_context)
        self.content_accession.SetFont(font_input_context)

        button_search_gene_name = wx.Button(panel_search, -1, '&Search', size=(100, 30), pos=(135, 63))
        button_search_gene_id = wx.Button(panel_search, -1, '&Search', size=(100, 30), pos=(100, 155))
        button_search_accession = wx.Button(panel_search, -1, '&Search', size=(100, 30), pos=(465, 155))
        button_search_gene_name.Bind(wx.EVT_BUTTON, self.button_search_gene_name)
        button_search_gene_id.Bind(wx.EVT_BUTTON, self.button_search_gene_id)
        button_search_accession.Bind(wx.EVT_BUTTON, self.button_search_accession)

    # 打开一个searching now frame，调用process类，返回字典，关闭searching now frame，将字典全局化，打开显示frame
    def button_search_gene_name(self, event):
        frame_searching_now = FrameSearchingNow(None, title='Searching', pos=(400, 260), size=(300, 80))
        frame_searching_now.Show()
        self.str_search_gene_name = self.content_gene_name.GetValue()
        getseqinfo = GetSeqInfo('genename', self.str_search_gene_name)
        global dict_gene_data
        dict_gene_data = getseqinfo.return_data()
        frame_searching_now.Destroy()

    def button_search_gene_id(self, event):
        frame_searching_now = FrameSearchingNow(None, title='Searching', pos=(400, 260), size=(300, 80))
        frame_searching_now.Show()
        self.str_search_gene_id = self.content_gene_id.GetValue()
        getseqinfo = GetSeqInfo('geneid', self.str_search_gene_id)
        global dict_gene_data
        dict_gene_data = getseqinfo.return_data()
        frame_searching_now.Destroy()

    def button_search_accession(self, event):
        frame_searching_now = FrameSearchingNow(None, title='Searching', pos=(400, 260), size=(300, 80))
        frame_searching_now.Show()
        self.str_search_accession = self.content_accession.GetValue()
        getseqinfo = GetSeqInfo('accession', self.str_search_accession)
        global dict_gene_data
        dict_gene_data = getseqinfo.return_data()
        frame_searching_now.Destroy()


class FrameSearchingNow(wx.Frame):
    def __init__(self, *args, **kw):
        super(FrameSearchingNow, self).__init__(*args, **kw)
        panel_searching_now = wx.Panel(self)
        self.count = 0
        self.gauge = wx.Gauge(panel_searching_now, -1, 1000, pos=(20, 20), size=(250, 15))
        self.gauge.SetBezelFace(3)
        self.gauge.SetShadowWidth(3)
        self.Bind(wx.EVT_IDLE, self.OnIdle)

    def OnIdle(self, event):
        self.count = self.count + 1
        if self.count == 300:
            self.count = 350
        if self.count == 500:
            self.count = 700
        if self.count == 900:
            self.count = 970
        self.gauge.SetValue(self.count)


a = wx.App(None)
frame_search = FrameSearch(None, title='Japonicusequence', pos=(300, 240), size=(600, 150))
frame_search.Show()
a.MainLoop()
