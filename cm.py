#!/usr/bin/env python
import networkx as nx
import math
import csv
import random
import random as rand
import sys
import community
import matplotlib.pyplot as plt
import numpy as np
import copy


f = open("cm.txt","w")
k = 10
pop = 200
maxgen = 1
r = 1000
p = 0.1

proMax = 0.0
proMin = 0.0
totMax = []

popMax = 0.0 #Guarda el maxina difusion del final del programa
poPfin = [] #guarda la semilla final que produce la difusion maxima.
#this method just reads the graph structure from the file
def buildG(G, file_, delimiter_):
    global Nodospajek
    Nodospajek = []
    #construct the weighted version of the contact graph from cgraph.dat file
    #reader = csv.reader(open("/home/kazem/Data/UCI/karate.txt"), delimiter=" ")
    reader = csv.reader(open(file_), delimiter=delimiter_)
    Arcos = 0
    cont = 0
    for line in reader:
        if Arcos == 0 and  line[0] != "*Arcs" and line[0] != "*vertices":
            Nodospajek.append(line)
        if Arcos == 1:
           if len(line) >  2:
              if float(line[2]) != 0.0:
                #line format: u,v,w
                G.add_edge(int(line[0]),int(line[1]),weight=float(line[2]))
           else:
                #line format: u,v
                G.add_edge(int(line[0]),int(line[1]),weight=1.0)
        if line[0] == "*Arcs" : Arcos = 1

def independent_cascade(G, seeds, steps=0):

    h0 = seeds
    h1 = []
    h2 = []
    ht = []
    vecinosItem = []
    for i in h0:
        vecinosItem = G.neighbors(i)
        for j in vecinosItem:
            if random.random() <= p:
               h1.append(j)
 
    vecinosItem = []
    for i in h1:
        vecinosItem = G.neighbors(i)
        for j in vecinosItem:
            if random.random() <= p:
               h2.append(j)


    htemp = h0 + h1 + h2
    for i in htemp:
        while i not in ht:
              ht.append(i)

    return ht
            
#defino mi grafo en grupos usando BGLL y los dibujo -> 
def comunidad(G):
   Grupo = {}
   #first compute the best partition

   partition = community.best_partition(G)
   values =  [partition.get(node) for node in G.nodes()]
   
   #drawing
   size = float(len(set(partition.values())))

   count = 0.
   comp = 0
   for com in set(partition.values()):
       comp += 1
       count = count + 1.
       list_nodes = [nodes for nodes in partition.keys()
                                if partition[nodes] == com]
       #print "Grupo:", comp, "tamano:", len(list_nodes) ,"nodos:", list_nodes
       Grupo[comp] = list_nodes
   #print("numero de componentes", comp)
   mod = community.modularity(partition,G)
   #print("modularidad:", mod)
   #Llamamos funcion candidato
   f.write("cm semilla = ")
   f.write(str(k))
   f.write(" poblacion = ")
   f.write(str(pop))
   f.write(" Generacion = ")
   f.write(str(maxgen))
   f.write(" Repetir = ")
   f.write(str(r))
   f.write(" \n")
   for i  in range(r):
       print i
       candidatos(G,comp, Grupo)

   f.write("promedio inicio = ") 
   f.write(str(proMin/r+0.0))
   f.write(" promedio final = ") 
   f.write(str(proMax/r))
   f.write(" maximo = ") 
   f.write(str(max(totMax)))
   f.write(" diferencia = ") 
   f.write(str((proMax - proMin)/r))
   f.close()
   
def candidatos(G,comp, Grupo):
   #Hallamos la comunidad significativa de menor tamano
   Cmin = len(G.nodes())/comp #definimos como el numero total de nodos por el numero total de grupos
   #Hallamos el valor de la comunidad maxima
   Cmax = 0 #se define como el numero maximo de nodos en un grupo
   Csig = 0 #suma todos los grupos que van a portar a la piscina de nodos candidatos
   PoolNodosCandidatos = 0 #NUMERO de nodos que seran candidatos
   ContGrupo = 0
   Cnodomx = 0 #Extrae el maximo de nodos por comunidad
   contnodosPool = 0
   nPool = []
   nodosPool =[]
   for i in Grupo:
       if len(Grupo[i]) > Cmax: 
          Cmax = len(Grupo[i])
       if len(Grupo[i]) >= Cmin:
          Csig += 1
          PoolNodosCandidatos = PoolNodosCandidatos + len(Grupo[i])
   for i in Grupo:
       ContGrupo += 1
       if len(Grupo[i]) >= Cmin:
          Cnodomx = 10*(len(Grupo[i]) - Cmin)/(Cmax - Cmin)+4
          #print("Comunidad No. ", ContGrupo,"TotalNodosGrp ", len(i), "aporta ", Cnodomx, "Nodos")
          nPool.append(pool(G, Grupo[i], Cnodomx)) #Crea la piscina, nodos piscina lo mejor de cada grupo
          contnodosPool = Cnodomx + contnodosPool #lleva el numero de nodos candidatos de una piscina
   for i in nPool:
         for j in i:
             nodosPool.append(j) #Crea el nodo piscina final.
   '''print("Tamano minimo Comunidad Significativa ", Cmin)
   print ("Tamano maximo comunidad Significativa ", Cmax)
   print ("Numero de Comunidades significativas ", Csig) 
   print ("Nodos prePool: ", PoolNodosCandidatos)
   print ("Nodos Pool: ", contnodosPool)
   print ("piscina Candidatos", nodosPool)'''

   poblacionP(G, nodosPool)




def pool(G,i,aport):

    #inicializo el vector piscina de cada grupo
    #vector que contiene todos los nodos piscina
    Pool = [x*0 for x in range(aport)]
    for x in range(aport):
        for y in i:
            if Pool[x] == 0 and x == 0: #Esta iniciando la piscina temporal
               Pool[x] = y
            if Pool[x] == 0 and x <> 0 and not (y in Pool):
               Pool[x] = y
            if not (y in Pool) and G.degree(y) > G.degree(Pool[x]): #Se escoge los mejores candidatos de cada grupo, segun su grado de conexion.
               Pool[x] = y
    return Pool
        

def poblacionP(G, Pool):
    Ppop = [] #Ppop
    Temp = [] #guarda el nodo temporal del candidato al azar
    verTemp = []    #verifica y guarda un nodo temporal para comprobar que no esta repetido
    PoolCandidate = [] #guarda la piscina
    TempCandidate = Pool #guarda el grupo de candidatos y luego se vacia cuando estos se escogen.
    #Escojo candidato al azar de los grupos 
    PoolCandidate = [i for i in Pool]
    for i in range(pop/2):
        for j in range(k):
            verTemp = random.choice(Pool)
            while verTemp in Temp:
                  verTemp = random.choice(Pool)
            Temp.append(verTemp)
        Ppop.append(list(Temp))
        Temp = []
    highgrade = 0
    nodehigh = []
    vecinode = []
    #buscamos la otra poblacion que usa la funcion SHD
    for i in range(pop/2+1, pop+1): 
        for j in range(k):
            if TempCandidate <> []:
                for x in TempCandidate:
                    nodehigh = []
                    if G.degree(x) > highgrade:
                       highgrade = G.degree(x)
                       nodehigh = x
            if TempCandidate == []:
                  ControlWhile = 0
                  nodoVi = 0
                  while ControlWhile == 0:
                        nodoVi = random.choice(PoolCandidate)
                        if not nodoVi in Temp:
                           Temp.append(nodoVi)
                           ControlWhile = 1
            if TempCandidate <> []:                  
                #codigo clave, aqui se escoge el mas alto y se remueve, pero primero hay que remover
                #los vecinos similares a este.
                Temp.append(nodehigh)
                vecinode = G.neighbors(nodehigh)
                for y in vecinode:
                        #Aqui se mira la simililaridad y se borra los similares
                        preds = nx.jaccard_coefficient(G, [(nodehigh, y)])
                        for u,v,p in preds:
                            if p > 0.6 and (y in TempCandidate):
                                TempCandidate.remove(y)
                if nodehigh in TempCandidate:
                    TempCandidate.remove(nodehigh)
            highgrade = 0
        Ppop.append(list(Temp))
        Temp = []

    
    c = 0
    popMayor = 0
    popGrade = [] #variable que contiene el S de la difusion de 2hop para cada semilla
    for i in Ppop:
            contgrade = 0
            diffusion = independent_cascade(G, i, steps = 2)
            contgrade =  len(diffusion)
            if popMayor < contgrade:
                          pPoolini = [] #variable para imprimir el inicio de la secuencia
                          popMayor = contgrade
                          pPoolini.append(i) 
            popGrade.append(contgrade)

    for i in range(0, maxgen+1):
        Ppop = cruce(G, Ppop, PoolCandidate, popGrade)

    global proMax, proMin, totMax

    #creamos un diccionario de nodos en G vs diffusion 2 step.
    proMax = proMax + popMax
    proMin = proMin + popMayor
    totMax.append(popMax)
    
    f.write("inicio ")
    f.write(str(pPoolini))
    f.write(" dif ")
    f.write(str(popMayor))
    f.write(" final ")
    f.write(str(poPfin))
    f.write(" difini ")
    f.write(str(popMax))
    f.write(" \n")
    #f.write("inicio", " ".join(poolPop(pPoolini)), "dif", str(poolMax(pPoolini)), "final", " ".join(poolPop(Ppop)), "difini", str(poolMax(Ppop)), "\n")

def duplicado(valores):
    repetido = []
    unico = []
    
    for x in valores:
	    if x not in unico:
		   unico.append(x)
	    else:
		   if x not in repetido:
			  repetido.append(x)

    if repetido == []:
       return True
    else:
       return False


#busqueda local mejor hijo. vamos a usar funcion de influencia.
def bestSon(G, Ppop, hijo, popGrade):

    #popGrade = [] #variable que contiene el S de la difusion de 2hop para cada semilla
    
    global popMax, poPfin
    #organizamos Ppop
    tempopGrade = []
    tempopGrade = [i for i in popGrade]
    tempopGrade.sort()

    #hallamos la difusion de los hijos
    sonGrade = []
    minSeed = 0
    ctrPop = 0
    itr = 0
    while ctrPop == 0:
          if tempopGrade[0] == popGrade[itr]:
             minSeed = itr
             ctrPop = 1
          itr = itr + 1
    contgrade = 0
    diffusion = independent_cascade(G, hijo, steps = 2)
    contgrade = len(diffusion)
    if tempopGrade[0] < contgrade:
       if not (hijo in Ppop):
          Ppop[minSeed] = hijo
          popGrade[minSeed] = contgrade
    c = 0
    popMax = 0.0
    for i in popGrade:
       c = c + 1
       if i > popMax:
          popMax = i
    c = 0
    for i in popGrade:
        if i == popMax:
            #poPfin = [] 
            poPfin = Ppop[c] 
        c = c + 1
    return Ppop     
    
    
  
def cruce(G, Ppop, Pool, popGrade):
    #Seleccion padres (2) por iteraccion.
    iposicion = 3
    rpadres = 1 #valor rand de cruce entre padres
    pmhijo = 1
    pc = 0.8 #probabilidad de cruce
    pm = 0.2 #probabilidad de mutacion
    ctrDuplicado = 0
    while ctrDuplicado == 0:
        padre1 = random.choice(Ppop)
        padre2 = random.choice(Ppop)
        while padre1 == padre2 or rpadres > pc:
              hijo1 = []
              hijo2 = []
              padre1 = random.choice(Ppop)
              padre2 = random.choice(Ppop)
              hijo1 = [i for i in padre1]
              hijo2 = [i for i in padre2]
              rpadres = random.random()
        #hijo 1 y 2
        for i in range(iposicion, len(padre1)):
            hijo2[i] = padre1[i]
            hijo1[i] = padre2[i]
        if (duplicado(hijo1) and duplicado(hijo2)):
           ctrDuplicado = 1
 
    #Verificacion de los vectores hijos.
    #mutacion
    if duplicado(hijo1) and duplicado(hijo2):
       #mutamos hijo 1
       for i in range(len(hijo1)):
           pmhijo = random.random()
           extraGen = 0
           if pmhijo < pm:
              contPosy = 0
              for y in hijo1:
                  extraGen = random.choice(Pool)
                  ControlWhile = 0
                  while (extraGen in hijo1) and ControlWhile == 0:
                        ControlWhile = 1
                        #Aqui se mira la simililaridad y se borra los similares
                        preds = nx.jaccard_coefficient(G, [(extraGen, y)])
                        contSim = 0
                        for u,v,p in preds:
                               if p > 0.6: 
                                  contSim = contSim + 1
                               if contSim == 0:
                                  hijo1[contPosy] = extraGen
                        contPosy = contPosy + 1
       #mutamos hijo 2
       for i in range(len(hijo2)):
           pmhijo = random.random()
           extraGen = 0
           if pmhijo < pm:
              contPosy = 0
              for y in hijo2:
                  extraGen = random.choice(Pool)
                  ControlWhile = 0
                  while not (extraGen in hijo2) and ControlWhile == 0:
                        ControlWhile = 1
                        #Aqui se mira la simililaridad y se borra los similares
                        preds = nx.jaccard_coefficient(G, [(extraGen, y)])
                        contSim = 0
                        for u,v,p in preds:
                               if p > 0.6: 
                                  contSim = contSim + 1
                               if contSim == 0:
                                  hijo2[contPosy] = extraGen
                        contPosy = contPosy + 1
    hijos = []
    hijos.append(hijo1)
    hijos.append(hijo2)
    #newPpop = []
    for i in hijos:
        newPpop = bestSon(G, Ppop, i, popGrade)
    return newPpop

                                  

def main(argv):
    if len(argv) < 2:
        sys.stderr.write("Usage: %s <input graph>\n" % (argv[0],))
        return 1
    graph_fn = argv[1]
    G = nx.Graph()  #let's create the graph first
    buildG(G, graph_fn, ' ')
    

    
    comunidad(G)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
