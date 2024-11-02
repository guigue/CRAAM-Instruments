import oRBD

d1=oRBD.RBD()
d2=oRBD.RBD()

d1.readRBDinDictionary('TestData/rs1150621.1700')
d2.readRBDinDictionary('TestData/rs1150621.1800')

d=oRBD.RBD()
d.concat([d1,d2])    


