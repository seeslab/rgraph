import os
import wx

ID_BUTTON1=110

SA_PARAMETERS = {
    'Iteration factor' : (1.0, 1001),
    'Cooling factor' : (0.995, 1002),
    'Randomizations' : (0, 1003),
    'Initial T' : (1., 1004),
    }

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# The MainWindow class
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
class MainWindow(wx.Frame):
    """The MainWindow class.

    This class takes care of the main window of the application.
    """
    
    # -------------------------------------------------------------------------
    # 
    # -------------------------------------------------------------------------
    def network_sizer(self):
        box = wx.StaticBox(self, -1, 'Network')
        sizer = wx.StaticBoxSizer(box, wx.VERTICAL)
        rows = {}

##         rows[aParam] = wx.BoxSizer(wx.HORIZONTAL)
##         nodesText = wx.StaticText(self, 0, 'Nodes')
##         text.SetFont(wx.Font(10, wx.SWISS, wx.NORMAL, wx.NORMAL))
##         rows[aParam].Add(text, 1, wx.CENTER)
##         textcontrol = wx.TextCtrl(self, 0, value=`PARAMETERS[aParam]`)
##         rows[aParam].Add(textcontrol, 1, wx.CENTER)
##         button = wx.Button(self, 1, "Default")
##         rows[aParam].Add(button, 1, wx.CENTER)
##         sizer.Add(rows[aParam], wx.EXPAND)

        return sizer

    # -------------------------------------------------------------------------
    # 
    # -------------------------------------------------------------------------
    def SA_parameter_sizer(self):

        # The overall frame
        # ---------------------------------------------------------------------
        sizer = wx.StaticBoxSizer(
            wx.StaticBox(self, -1, 'Simulated annealing parameters'),
            wx.VERTICAL
            )

        # Initialize
        # ---------------------------------------------------------------------
        self.SAParametersControl = {}

        # Create all titles, textcontrols, and buttons
        # ---------------------------------------------------------------------
        box = wx.FlexGridSizer(cols=2)
        for aParam in SA_PARAMETERS:
            title = wx.StaticText(self, 0, aParam)
##             title.SetFont(wx.Font(10, wx.DEFAULT, wx.DEFAULT, wx.DEFAULT))
            textControl = wx.TextCtrl(
                self,
                SA_PARAMETERS[aParam][1],
                value=`SA_PARAMETERS[aParam][0]`
                )
            textbox = wx.BoxSizer(wx.HORIZONTAL)
            textbox.Add(title, 1, wx.CENTER)
            textbox.Add(textControl, 1, wx.CENTER)
            box.Add(textbox, 1, wx.EXPAND)
            self.SAParametersControl[aParam] = textControl
            button = wx.Button(
                self,
                SA_PARAMETERS[aParam][1],
                "Set default")
            box.Add(button, 1, wx.CENTER)
        sizer.Add(box)

        # Add all event handlers
        # ---------------------------------------------------------------------
        wx.EVT_BUTTON(
            self,
            SA_PARAMETERS['Iteration factor'][1],
            self.OnIterationFactorDefault
            )

        return sizer

    # -------------------------------------------------------------------------
    # Constructor
    # -------------------------------------------------------------------------
    def __init__(self,parent,id,title):
        """MainWindow constructor."""

        # Init the parent class
        # ---------------------------------------------------------------------
        wx.Frame.__init__(self,parent,wx.ID_ANY, title)

        # Some stuff
        # ---------------------------------------------------------------------
        self.dirName=''
        self.fileName=''
        self.networkNodes=''
        self.networkLinks=''
##         self.control = wx.TextCtrl(self, 1, style=wx.TE_MULTILINE)

        # Statusbar at the bottom of the frame
        # ---------------------------------------------------------------------
        self.CreateStatusBar()

        # Menu bar (and menus) at the top of the frame
        # ---------------------------------------------------------------------
        # The manu bar
        menuBar = wx.MenuBar()
        self.SetMenuBar(menuBar)

        # The File menu
        filemenu= wx.Menu()
        filemenu.Append(wx.ID_OPEN, "&Open", "Open a file to edit")
        wx.EVT_MENU(self, wx.ID_OPEN, self.OnOpen)

        filemenu.AppendSeparator()
        filemenu.Append(wx.ID_ABOUT, "&About", "Information about this program")
        wx.EVT_MENU(self, wx.ID_ABOUT, self.OnAbout)

        filemenu.AppendSeparator()
        filemenu.Append(wx.ID_EXIT, "E&xit", "Terminate the program")
        wx.EVT_MENU(self, wx.ID_EXIT, self.OnExit)

        menuBar.Append(filemenu, "&File") # Adding the "filemenu" to the MenuBar

        # The sizers
        self.sizer = wx.GridSizer(cols=2)
        self.NetworkSizer = self.network_sizer()
        self.sizer.Add(self.NetworkSizer, wx.ALL)
        self.SAParameterSizer = self.SA_parameter_sizer()
        self.sizer.Add(self.SAParameterSizer, wx.ALL)
        
        # Layout sizers
        self.SetSizer(self.sizer)
        self.SetAutoLayout(1)
        self.sizer.Fit(self)
        self.Show(1)

    def OnAbout(self,e):
        d= wx.MessageDialog( self, " A sample editor \n"
                            " in wxPython","About Sample Editor", wx.OK)
                            # Create a message dialog box
        d.ShowModal() # Shows it
        d.Destroy() # finally destroy it when finished.
    def OnExit(self,e):
        self.Close(True)  # Close the frame.
    def OnOpen(self,e):
        """ Open a file"""
        dlg = wx.FileDialog(self, "Choose a file", self.dirName, "", "*.*", wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            self.filename=dlg.GetFilename()
            self.dirName=dlg.GetDirectory()
            f=open(os.path.join(self.dirName, self.filename),'r')
            self.control.SetValue(f.read())
            f.close()
        dlg.Destroy()

    def OnIterationFactorDefault(self, e):
        print self.SAParametersControl['Iteration factor'].GetValue()
##         self.SAParameters['Iteration factor'] = SA_PARAMETERS['Iteration factor'][0]
##         print self.SAParameters['Iteration factor']


app = wx.PySimpleApp()
frame = MainWindow(None, -1, "NetCarto")
app.MainLoop()
