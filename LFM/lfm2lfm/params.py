import ConfigParser

class params():
    def __init__(self,ConfigFileName):
        config = ConfigParser.ConfigParser()
        config.read(ConfigFileName)
        
        self.ni = config.getint('Dimensions','NI')
        self.nj = config.getint('Dimensions','NJ')
        self.nk = config.getint('Dimensions','NK')

        self.initLFMfile = config.get('Output','Prefix')+'_%dx%dx%d_mhd_0000000.hdf'%(self.ni,self.nj,self.nk)
        self.dirInitLFMfile = config.get('Output','Dir')
    
        self.xscale = config.getfloat('Normalization','Xscale')
        self.nscale = config.getfloat('Normalization','Nscale')
        self.bscale = config.getfloat('Normalization','Bscale')

        self.rmin = config.getfloat('GridSpecs','RMIN')
        self.rmax = config.getfloat('GridSpecs','RMAX')
        self.thetamin = config.getfloat('GridSpecs','THETAMIN')
        
        self.gamma = config.getfloat('Constants','gamma')
        self.NO2   = config.getint('Constants','NO2')
        self.Tsolar = config.getfloat('Constants','Tsolar')
        
        if config.has_section('MAS'):
            self.masdir = config.get('MAS','masdir')
            self.masTimeLabel = config.get('MAS','masTimeLabel')
            self.masFrame = config.get('MAS','masFrame')
            self.masFakeRotation = config.getboolean('MAS','masFakeRotation')

        if config.has_section('WSA'):
            self.wsaFile = config.get('WSA','wsafile')
            self.gaussSmoothWidth = config.getint('WSA','gauss_smooth_width')
            self.plots = config.getboolean('WSA','plots')
            self.densTempInfile = config.getboolean('WSA','density_temperature_infile')
            self.corotatingFrame = config.getboolean('WSA','corotating_frame')

        if config.has_section('LFM'):
            self.inFile = config.get('LFM','infile')

        self.dumpInit = config.getboolean('DUMPS','init')
        self.dumpBC   = config.getboolean('DUMPS','BC')

