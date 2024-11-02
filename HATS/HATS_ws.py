import os
import numpy as np
import warnings
import pdb
import datetime as dt

############ Global Variables ##########
__version__       = "2024-09-09T1700AST"   
input_path       = '/data/HATS/aux'
#######################################

######################################
#
# HATS_ws
#     A class to read HATS weather station
#
# Usage
#     import HATS_ws
#     ws=HATS_ws('hats-YYYY-MM-DD')
#
# Configuration
#     input_path : set default directory
#
# Author
#      @guiguesp, written in Argentina
#
#####################################

class WS(object):
    
    def __str__(self):
        return 'A Class representing HATS Weather Station Data'
    
    def __init__(self,date='',InputPath=input_path):
        
        warnings.filterwarnings('ignore')
        if (InputPath[-1] != '/') :
            InputPath+='/'        
        if (os.path.exists(InputPath)):
            self.MetaData={"InputPath":InputPath}
        else:
            raise ValueError("\n\n No Input Path {} found \n\n".format(InputPath))


        if (len(date) == 10):
            YMD = date[:10]
        else:
            date = dt.datetime.now()
            YMD="{0:4d}-{1:02d}-{2:02d}".format(date.year,date.month,date.day)

        self.MetaData.update({'Date':YMD})
            
        fname = 'hats-'+YMD+'.ws'    # ASCII file
        self.MetaData.update({'Filename':fname})
        self.data = {}
        self.from_file(fname)
        self.Units = {'Temperature':'°C','Humidity':'%','Pressure':'HPa'}
        
        return

    def from_file(self,fname):
        
        fullpathname = self.MetaData['InputPath']+fname
        if (not os.path.exists(str(fullpathname))):
            raise ValueError('\n\n File {0:s} not found \n\n'.format(fullpathname))

        DateTime = []
        Temperature = []
        Humidity = []
        Pressure = []
        
        f = open(fullpathname,'r',errors='ignore')
        for line in f:
            s = line.split(sep=',')
            if (len(s) == 5):
                try:
                    DateTime.append(dt.datetime.fromisoformat(s[0]))
                    Temperature.append(float(s[2].split(sep='=')[1][:-1]))
                    Humidity.append(float(s[3].split(sep='=')[1][:-1]))
                    temp = s[4].split(sep='=')[1]
                    Pressure.append(float(temp[:temp.find('H')]))
                except:
                    break
                    
        f.close()
        self.data.update({'time':np.asarray(DateTime),
                          'temperature':np.asarray(Temperature),
                          'humidity':np.asarray(Humidity),
                          'pressure':np.asarray(Pressure)})    
        return

    def to_csv(self):
        
        csv_fname=self.MetaData['Filename'][:-3]+'.csv'
        f=open(csv_fname,'w')
        f.write("time,temperature [C],Humidity [%],Pressure [Hpa]\n")
        for i in np.arange(self.data['time'].shape[0]):
            line = self.data['time'][i].isoformat()            + ','
            line+='{0:6.1f}'.format(self.data['temperature'][i]) + ','
            line+='{0:6.1f}'.format(self.data['humidity'][i])    + ','
            line+='{0:6.1f}'.format(self.data['pressure'][i])    + '\n'
            f.write(line)

        f.close()

        return
    
