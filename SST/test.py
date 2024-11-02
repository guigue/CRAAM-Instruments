# A simple test for the scripts

import oSST
d=oSST.SST()
d.data_path='./TestData'
d.initial_time=d.str2datetime('2002-12-21 15')
d.final_time=d.str2datetime('2002-12-21 19')
d.data_type='bi'
d.read()

t=oSST.TimeAxis()
t.getTimeAxis(d,'dt')

y=oSST.yAxis()
y.getValues(d,'adc_sigma')

import matplotlib.pyplot as plt
plt.plot(t.time,y.adc_sigma[:,1])
plt.xlabel('UT')
plt.ylabel(y.AxisName[1]+' ('+y.AxisUnits+')')
plt.show()


#import oSSTMap
#m=oSSTMap.Map(d)
#m.getDataMaps(d)
#tp=oSST.TotalPower()
#tp.getTotalPower(d,oSST.Beams.two)



