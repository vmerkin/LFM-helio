import ConfigParser

class params():
    def __init__(self,ConfigFileName):
        config = ConfigParser.ConfigParser()
        config.read(ConfigFileName)

        self.masdir = config.get('Inputs','masdir')
        self.lfmfile = config.get('Inputs','lfmfile')
        
        self.variable = config.get('Inputs','variable')
        self.mas_time_label = config.getint('Inputs','mas_time_label')
        self.gamma = config.getfloat('Inputs','gamma')

        self.kslice = config.getint('Inputs','kslice')
        self.rslice = config.getfloat('Inputs','rslice')

        self.slice = config.get('Inputs','slice')
        
