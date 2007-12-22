import os
import wx

import network

SA_PARAMETERS = {
    'Iteration factor' : (1.0, 1001),
    'Cooling factor' : (0.995, 1002),
    'Randomizations' : (0, 1003),
    'Initial T' : (1., 1004),
    }

SA_PARAMETERS_ID = dict([(SA_PARAMETERS[anSAParam][1], anSAParam)
                         for anSAParam in SA_PARAMETERS])

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# The MainWindow class
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
class MainWindow(wx.Frame):
    """ The MainWindow class.

    This class takes care of the main window of the application.
    """
    
    # -------------------------------------------------------------------------
    # 
    # -------------------------------------------------------------------------
    def network_sizer(self):
        box = wx.StaticBox(self, -1, 'Network')
        sizer = wx.StaticBoxSizer(box, wx.VERTICAL)
        rows = {}

        self.NNodesText = wx.StaticText(self, 0, 'Nodes: 0')
        sizer.Add(self.NNodesText)

        self.NLinksText = wx.StaticText(self, 0, 'Links: 0')
        sizer.Add(self.NLinksText)

        return sizer

    # -------------------------------------------------------------------------
    # The SA parameters sizer
    # -------------------------------------------------------------------------
    def SA_parameter_sizer(self):
        """ The SA parameters sizer.

        This sizer contains text control areas to enter the simulated
        annealing parameters, as well as buttons to set the defaults.
        """

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
            # The title
            title = wx.StaticText(self, 0, aParam)
            # Text control
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
            # Button
            button = wx.Button(
                self,
                SA_PARAMETERS[aParam][1],
                "Set default")
            # Button event handler
            button.Bind(wx.EVT_BUTTON, self.OnSAParameterDefault)
            # Add the button
            box.Add(button, 1, wx.CENTER)
        sizer.Add(box)
        
        # Done
        # ---------------------------------------------------------------------
        return sizer

    # -------------------------------------------------------------------------
    # Constructor
    # -------------------------------------------------------------------------
    def __init__(self,parent,id,title):
        """ MainWindow constructor. """

        # Init the parent class
        # ---------------------------------------------------------------------
        wx.Frame.__init__(self,parent,wx.ID_ANY, title)

        # Some stuff
        # ---------------------------------------------------------------------
        self.dirName = ''
        self.fileName = ''
        self.fullFileName = ''
        self.networkNodes = ''
        self.networkLinks = ''

        # Statusbar at the bottom of the frame
        # ---------------------------------------------------------------------
        self.CreateStatusBar()

        # Menu bar (and menus) at the top of the frame
        # ---------------------------------------------------------------------
        # The manu bar
        menuBar = wx.MenuBar()
        self.SetMenuBar(menuBar)

        # The File menu
        filemenu = wx.Menu()
        filemenu.Append(wx.ID_OPEN, "&Open", "Open a network file")
        wx.EVT_MENU(self, wx.ID_OPEN, self.OnOpen)

        filemenu.AppendSeparator()
        filemenu.Append(wx.ID_EXIT, "&Quit", "Terminate NetCarto")
        wx.EVT_MENU(self, wx.ID_EXIT, self.OnExit)

        menuBar.Append(filemenu, "&File") # Adding "filemenu" to MenuBar

        # The Help menu
        helpmenu = wx.Menu()
        helpmenu.Append(wx.ID_ABOUT,
                        "&About",
                        "Information about NetCarto")
        wx.EVT_MENU(self, wx.ID_ABOUT, self.OnAbout)

        menuBar.Append(helpmenu, "&Help") # Adding "help" to MenuBar

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

    # -------------------------------------------------------------------------
    # Set default SA parameter
    # -------------------------------------------------------------------------
    def OnSAParameterDefault(self, event):
        """ Set default SA parameter.

        Which parameter is actually set to the default is determined
        by the ID of the event.
        """
        ID = event.GetId()
        param = SA_PARAMETERS_ID[ID]
        self.SAParametersControl[param].SetValue(`SA_PARAMETERS[param][0]`)

    # -------------------------------------------------------------------------
    # Open network file
    # -------------------------------------------------------------------------
    def OnOpen(self, event):
        """ Open a network file.

        Besides opening the file, we check that it is a valid network
        file, that is, a list of adjacencies.
        """
        # Open the file dialog
        # ---------------------------------------------------------------------
        dlg = wx.FileDialog(self,
                            "Choose a file",
                            self.dirName,
                            "",
                            "*.*",
                            wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            fileName = dlg.GetFilename()
            dirName = dlg.GetDirectory()
            fullFileName = os.path.join(dirName, fileName)
            success, nNodes, nLinks = \
                  network.read_undirected_network(fullFileName)
            if success: # Good!
                self.fileName = fileName
                self.dirName = dirName
                self.fullFileName = fullFileName
                self.networkNodes = `nNodes`
                self.networkLinks = `nLinks`
                self.NNodesText.SetLabel('Nodes: %s' % self.networkNodes)
                self.NLinksText.SetLabel('Links: %s' % self.networkLinks)
            else:  # Show error message
                d = wx.MessageDialog(
                    self,
                    "%s is not a valid network file" % fileName,
                    "Error opening file",
                    wx.OK
                    )
                d.ShowModal()
                d.Destroy()
                
                
        # Done
        # ---------------------------------------------------------------------
        dlg.Destroy()

    # -------------------------------------------------------------------------
    # Quit
    # -------------------------------------------------------------------------
    def OnExit(self, event):
        """ Quit the program. """
        self.Close(True)  # Close the frame.

    # -------------------------------------------------------------------------
    # About
    #--------------------------------------------------------------------------
    def OnAbout(self, event):
        d= wx.MessageDialog(self,
                            "The Network Cartography software\n"
                            "\n"
                            "by Roger Guimera\n"
                            "and the Amaral Lab",
                            "About NetCarto",
                            wx.OK)
        d.ShowModal()
        d.Destroy()


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# The main loop
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
if __name__ == '__main__':
    app = wx.PySimpleApp()
    frame = MainWindow(None, -1, "NetCarto")
    app.MainLoop()
