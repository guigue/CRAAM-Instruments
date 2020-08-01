import oRBD
import numpy as np
import datetime as dt
from CraamTools.filters import RunningMean as rm
import pdb

def get_par():
    d=oRBD.RBD()
    d.readRBDinDictionary('rs1170914.1700')

    # Liquid Nitrogen 1
    i0=7900
    i1=8500
    r=np.arange(i0,i1+1)
    nitro1  = np.zeros(6)
    snitro1 = np.zeros(6)
    for i in np.arange(6):
       nitro1[i]  = np.mean(d.Data['adcval'][r,i])
       snitro1[i] = np.std(d.Data['adcval'][r,i])

    # Liquid Nitrogen 2
    i0 = 9300
    i1 = 10200
    r  = np.arange(i0,i1+1)
    nitro2  = np.zeros(6)
    snitro2 = np.zeros(6)

    for i in np.arange(6):
        nitro2[i]=np.mean(d.Data['adcval'][r,i])
        snitro2[i]=np.std(d.Data['adcval'][r,i])

    # Water+Ice
    i0 = 12475
    i1 = 13000
    r  = np.arange(i0,i1+1)
    hielo1  = np.zeros(6)
    shielo1 = np.zeros(6)
    for i in np.arange(6):
        hielo1[i]=np.mean(d.Data['adcval'][r,i])
        shielo1[i]=np.std(d.Data['adcval'][r,i])

    i0 = 21950
    i1 = 22650
    r  = np.arange(i0,i1+1)
    hielo2  = np.zeros(6)
    shielo2 = np.zeros(6)
    for i in np.arange(6):
        hielo2[i]  = np.mean(d.Data['adcval'][r,i])
        shielo2[i] = np.std(d.Data['adcval'][r,i])

    b  = oRBD.RBD()
    b.readRBDinDictionary('bi1170914')
    b.CorrectAuxiliary()
    br = b.extract([dt.datetime(2017,9,14,17,0),dt.datetime(2017,9,14,18,0)])
    x  = (br.Data['target']//32 == 1)
    cold  = np.zeros(6)
    scold = np.zeros(6)
    for i in np.arange(6):
        cold[i]  = br.Data['adc'][x,i].mean()
        scold[i] = br.Data['adc'][x,i].std()

    tcold = br.Data['amb_temp'][x].mean()
    stcold = br.Data['amb_temp'][x].std()

    x = (br.Data['target']//32 == 2)
    hot  = np.zeros(6)
    shot = np.zeros(6)
    for i in np.arange(6):
        hot[i]  = br.Data['adc'][x,i].mean()
        shot[i] = br.Data['adc'][x,i].std()

    thot = br.Data['hot_temp'][x].mean()
    sthot = br.Data['hot_temp'][x].std()

    temperature=np.array([77.0, 273.15, tcold+273.15,thot+273.15])

    adc0=np.array([nitro2[0],hielo1[0],cold[0],hot[0]] )
    adc1=np.array([nitro2[1],hielo1[1],cold[1],hot[1]] )
    adc2=np.array([nitro2[2],hielo1[2],cold[2],hot[2]] )
    adc3=np.array([nitro2[3],hielo1[3],cold[3],hot[3]] )
    adc4=np.array([nitro2[4],hielo1[4],cold[4],hot[4]] )
    adc5=np.array([nitro2[5],hielo1[5],cold[5],hot[5]] )

    print('\nReceptor 1')
    print(np.polyfit(adc0,temperature,1))

    print('\nReceptor 2')
    print(np.polyfit(adc1,temperature,1))

    print('\nReceptor 3')
    print(np.polyfit(adc2,temperature,1))

    print('\nReceptor 4')
    print(np.polyfit(adc3,temperature,1))

    print('\nReceptor 5')
    print(np.polyfit(adc4,temperature,1))

    print('\nReceptor 6')
    print(np.polyfit(adc5,temperature,1))

    return temperature, nitro2, hielo1, cold, hot
