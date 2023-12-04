import random
import math
import time

def distancia (p, q):
    return math.sqrt((p[0]-q[0])**2 + (p[1]-q[1])**2 + (p[2]-q[2])**2)

def imprimaMatriz (arq, M):
    fc = arq+ ".cmp"
    f = open (fc, "w")
    n=len(M)
    #print(n)
    #f.write(str(len(M))+"\n")
    for i in range (len(M)):
        for j in range (len(M[i])):
            f.write(str(M[i][j])+" ")
        f.write("\n")
    f.close()
    fc = arq+ ".inc"
    f = open (fc, "w")
    f.write(str(len(M))+"\n")
    i=0
    while i< len(M):
       j=0
       seq=""
       while j<=i :
          seq = seq + str(M[i][j])+" "
          j=j+1
    
       f.write(seq+"\n")   #print(seq)
       i=i+1
    f.close()

    
# definindo quantas OTUs vai gerar e o arquivo onde imprimir
print ("Digite o número de pontos: ")
n = int(input())
print ("Digite o nome do arquivo: ")
arq = input() 


start = time.time()

# fazendo uma lista de pontos no plano
L = []
for i in range(n):
   x = random.uniform(0,100)#n**4)
   y = random.uniform(0,100)#n**4)
   z = random.uniform(0,100)#n**4)
   L.append([x, y, z])

# criando a matriz
M = []
for i in range (n):
    M.append([])
    for j in range(0,n):   #i+1):
        M[i].append([])

for i in range (n): M[i][i] = 0

i = 0
while i < n:
    j = 0
    while j < i:
        M[i][j] = int(distancia (L[i], L[j])+0.5)
        M[j][i] = M[i][j]
        j += 1
    i += 1

# pontos impróprios

i = 0
while i < n:
    j = 0
    while j < i:
        k = i+1
        while k < j:
            if M[i][j] > M[i][k] + M[k][j]:
                print (str([i,j,k]), str([M[i][j],M[i][k],M[k][j]]) + " ")
            k += 1
        j += 1
    i += 1


imprimaMatriz (arq, M)

end = time.time()

print (end - start)
