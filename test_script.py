import numpy as np

nso = 4

e0 = np.random.rand(nso)
eri = np.random.rand(nso,nso,nso,nso)

t1 = np.random.rand(nso,nso)
t2 = np.random.rand(nso,nso,nso,nso)

s1 = np.zeros((nso,nso))
s2 = np.zeros((nso,nso,nso,nso))
s3 = np.zeros((nso,nso,nso,nso,nso,nso))
s4 = np.zeros((nso,nso,nso,nso,nso,nso,nso,nso))

eri2 = eri*0

for a in range(nso):
    for b in range(nso):
        for c in range(nso):
            for d in range(nso):
                eri2[a,b,c,d] = (1/2)*( eri[a,b,c,d] - eri[b,a,c,d] )

for a in range(nso):
    for b in range(nso):
        for c in range(nso):
            for d in range(nso):
                eri[a,b,c,d] = (1/2)*( eri2[a,b,c,d] - eri2[a,b,d,c] )


x = np.random.rand(nso)
y = np.sqrt(1 - x**2)
