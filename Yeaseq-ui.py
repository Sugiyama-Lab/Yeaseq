import wx
from yeaseq import wx_ui
# from yeaseq import yeaseq_config


if __name__ == '__main__':
    app_main = wx.App()  # redirect=True, filename=yeaseq_config.WorkDir.get_log_path())
    screen_size = wx.DisplaySize()
    pos_middle = (screen_size[0] / 4, screen_size[1] / 4)
    frame_search = wx_ui.SearchFrame(None, pos=pos_middle, title='Yeaseq')  # , size=yeaseq_config.GUIConfig.SearchSizeLow)
    frame_search.Show(True)
    app_main.MainLoop()
