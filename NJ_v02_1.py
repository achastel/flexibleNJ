from ete3 import Tree
import time
import sys

_tamID = 1

def descreve(T, otu):
   st="("
   i=len(T)-1
   r=T[i]
   #print("OTU=",otu)
   
   while i>=0 and r[4]!=otu   :
      i=i-1
      r=T[i]

   #print("descrevendo ", r[0])
  
   if i>=0 :
      if len(r[0])>_tamID:
         st = st + descreve(T, r[0])
      else:
         st = st + r[0]
      st = st + ","

      if len(r[2])>_tamID:
         st = st + descreve(T, r[2])
      else:
         st = st + r[2]
      st = st + ")"
   else:
      st = ""
 
   return st


def createNewickTree(T):
   i = len(T)-1
 
   #print ("====================================================", i)
#   while i>0 :
   r = T[i]
   st="("	+ r[0]
#   print("descreve ",r[0])
#   if len(r[0]) > _tamID :
#      st = st + descreve(T, r[0])
#   else:
#      st = st + r[0]
   st =st +","+r[2]
       
#   print("descreve ",r[2])
#   if len(r[2]) > _tamID :
#      st = st + descreve(T, r[2])
#   else:
#      st = st + r[2]
   st = st + ");"
#       i=i-1

   #print (st)   
   return st



class MatrizDistancia:
    def __init__(self, arq):
        self.L = []                               # lista de nomes das OTUs:: L[i] = u <=> u eh a i-esima OTU em L
        self.H = {}                               # lista dos indices de L:: H[u] = i <=> L[i] = u
        self.n = 0                                # quantidade de OTUs  
        self.d = []                               # distancia entre OTUs:: matriz triangular. Para i > j, d[i][j] eh a distancia entre L[i] e L[j] 
        self.D = []                               # soma das distancias entre L[i] e as demais OTUs, i. e, toda a linha:: D[i] = \sum{j < i} d[i][j] + \sum{j > i} d[j][i]
        
        #### LEITURA DO ARQUIVO DE ENTRADA ####
        # ================================================================
        f = open(arq, "r")
        linha = f.readline().split()
        st=""
        st=linha[0]
        #print(linha[0])
        
        # definindo n
        #self.n = len(linha)                       
        self.n = int(st)                       

        # definindo L e H
        for i in range(self.n):
            self.L.append(str(i))   #linha[i])
            self.H[str(i)] = i           #self.H[linha[i]] = i
  
        # definindo a matriz triangular d
        for i in range(self.n): self.d.append([])

        i = 0
        while i < self.n:
            linha = f.readline().split()
            self.d[i] =  list(map(float, linha[0:i+1]))
            i = i + 1

        f.close()
        # ================================================================

        # definindo D
        
        # ================================================================
        self.D = [0] * self.n                     # soma de toda a linha i 

        for i in range (self.n):
            for j in range(0,i):         self.D[i] = self.D[i] + self.d[i][j]
            for j in range(i+1, self.n): self.D[i] = self.D[i] + self.d[j][i]
        # ================================================================

    def vazia(self):
        return (self.n == 0) 
        
    def __str__(self):
        # ================================================================
        # valor de n
        saida = "n = " + str(self.n) + "\n\n"
        # ================================================================

        # ================================================================
        # valor de L, H e d
        saida = saida + "matriz d (completa) = \n"
        for i in range(self.n): saida = saida + "\t" + self.L[i] + " "
        saida = saida + "\n"

        for i in range(len(self.d)):
            saida = saida + self.L[i] + "\t"
            for j in range (len(self.d[i])): saida = saida + str(self.d[i][j]) + "\t"
            saida = saida + "0.0" + "\t"
            j = i + 1
            while j < len(self.d):
                saida = saida + str(self.d[j][i]) + "\t"
                j = j + 1
            saida = saida + "\n"
        # =================================================================
        # matriz d bruta
        saida = saida + "\nMatriz d Bruta (Estrutura de dados usada)\n"
        for i in range(self.n): saida = saida + "\t" + str(self.L[i])
        for i in range(len(self.d)):
            saida = saida + "\n" + str(self.L[i]) + "\t"
            for j in range (len(self.d[i])): saida = saida + str(self.d[i][j]) + "\t"
        # =================================================================
                
        # =================================================================
        # =================================================================
        # hash H
        saida = saida + "\n\nhash H: "
        saida = saida + str(self.H)
        # =================================================================

        # =================================================================
        # =================================================================
        # matriz D
        saida = saida + "\n\nmatriz D: "
        saida = saida + str(self.D)
        # =================================================================


        return saida               

    def neighborJoining (self, u, v):
        nw = "("
        # se n == 2, remove u e v da matriz e devolve (u, v, d[L[u]][L[v]])
        # se n > 2, insere uma nova OTU uv, e remove u e v da matriz e devolve (u, v, uv, d[L[u]][L[uv]], d[L[uv]][L[v]])
        # =================================================================
        # =================================================================
        # caso especial n = 2, nao precisa fazer nada
        if self.n == 2:
            self.n = 0
            if len(u)<2 :
              nw = nw + u 
            nw = nw + "," 
            if len(v)<2:
              nw = nw + v 
            nw = nw +")"
        
            return u, 0.5 * (int(100 * self.d[1][0])/100.), v, 0.5 * (int(100 * self.d[1][0])/100.) , "raiz" , nw
            #return u, v, "raiz", 0.5 * (int(100 * self.d[1][0])/100.), 0.5 * (int(100 * self.d[1][0])/100.) 
        # =================================================================
        # =================================================================
        #print ("--------NEIGHBOR-JOINING------------")
        L = self.insereOTU(u, v)

        self.removeOTU(u)
        self.removeOTU(v)
        return L

    def insereOTU (self, u, v):
        novo = "("+ u +","+ v+")"
        self.d.append([0] * self.n)
        self.L.append(novo)
        h = self.H[novo] = self.n
        #print (u, v, novo)
        
        i = self.H[u]
        j = self.H[v]
        
        if i < j: 
           i, j = j, i
           u, v = v, u

        self.d[h][i] = abs(0.5 * self.d[i][j] + 1.0 * (self.D[i] - self.D[j])/(2 * (self.n - 2)))
        self.d[h][j] = abs(self.d[i][j] - self.d[h][i])

        for k in range(j):      self.d[h][k] = abs(0.5 * (self.d[j][k] + self.d[i][k] - self.d[i][j]))
        for k in range(j+1, i): self.d[h][k] = abs(0.5 * (self.d[k][j] + self.d[i][k] - self.d[i][j]))
        for k in range(i+1, self.n): self.d[h][k] = abs(0.5 * (self.d[k][j] + self.d[k][i] - self.d[i][j]))

        for k in range(self.n): self.D[k] += self.d[h][k]
        self.D.append(0)
        for k in range(self.n): self.D[self.n] += self.d[h][k]
            
        self.n = self.n + 1
        nw = "("
        if len(u)<2 :
          nw = nw + u 
        nw = nw + "," 
        if len(v)<2:
          nw = nw + v 
        nw = nw +")"
      
        return u, int(100*self.d[h][i])/1000 ,v,  int(100*self.d[h][j])/1000 , novo , nw


    def removeOTU(self, u):
        self.n = self.n - 1
        i = self.H[u]
        h = self.n

        for k in range(i):
            self.D[k] = self.D[k] - self.d[i][k]
            self.d[i][k] = self.d[h][k]

        self.D[i] = self.D[h] - self.d[h][i]

        for k in range(i+1, self.n):
            self.D[k] = self.D[k] - self.d[k][i] 
            self.d[k][i] = self.d[h][k]

        del self.H[self.L[i]]    
        self.L[i] = self.L[h]
        self.L.pop()
        self.H[self.L[i]] = i
        self.d.pop()
        self.D.pop()    

    
    
class Floresta:
    pass
    #     # =================================================================
    #     # nos de T ja calculados.
    #     saida = saida + "\nArvore T:\n" + str(self.T)
    #     # =================================================================


class MatrizQ:
    # matriz Q referente a uma MatrizDistancia dist
    def __init__(self, dist):
        self.L = [[]] * dist.n  # lista dos nomes das OTUs da MatrizDistancia dist
        self.Q = []             # matriz Q
        self.n = dist.n         # quantidade de linhas da matriz

        for i in range(dist.n):
            self.L[i] = dist.L[i]
            
            self.Q.append([])
            # calcular Q[i][j]
            for j in range(i): self.Q[i].append((dist.n - 2) * dist.d[i][j] - (dist.D[i] + dist.D[j]))
      

    def __str__(self):
        saida = "Matriz Q\n"
        for i in range(self.n): saida = saida + "\t" + str(self.L[i])
        for i in range(len(self.Q)):
            saida = saida + "\n" + str(self.L[i]) + "\t"
            for j in range (len(self.Q[i])): saida = saida + str(self.Q[i][j]) + "\t"
        return saida

    def minimo(self):
        u = 1
        v = 0
        for i in range (len(self.Q)):
            for j in range (len(self.Q[i])):
                if self.Q[i][j] < self.Q[u][v]: u, v = i, j

        return self.L[u], self.L[v]



T = [] # arvore
arquivo = sys.argv[1]
arq2 = sys.argv[2]
d = MatrizDistancia(arquivo)
t1 = time.time()
c=0
while not d.vazia():
   Q = MatrizQ (d)
   c=c+1
   x,y = Q.minimo()
   #print ("minimo = ", x, y)
   T.append(d.neighborJoining(x, y))

tempo = time.time() - t1
#print("Neighboor-Joining - Original")
#print("Entrada: ",arquivo,"\n")
stNW = createNewickTree(T)
#print(stNW)
a=Tree(stNW)
#print (a.write(format=9)) # (A:1.000000,(B:1.000000,(E:1.000000,D:1.000000)Internal_1:0.500000)Internal_2:0.500000);
#print(a)
#print (T)
#print("\nV0.2 iterações=",c," execução: ", tempo,"\n\n")
print(tempo)
#a.write(format=9)
#print(a)

saida = arq2+".py"
f = open(saida, "a+")
f.write("from ete3 import Tree\n\n")
f.write("#TEMPO = "+str(tempo)+" s\n")
f.write(arq2+"= Tree(\""+a.write(format=9)+"\")\n\n")
f.close()
