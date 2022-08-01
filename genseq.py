# uso: python3 genseq.py (quantidade de caracteres a serem gerados)

import sys, random, string

def random_char(y):
    return ''.join(random.choice(string.ascii_letters) for x in range(y))

fA = open(sys.argv[1] + "_A.in", "x") 
fB = open(sys.argv[1] + "_B.in", "x") 

sA = random_char(int(sys.argv[1]))
sB = random_char(int(sys.argv[1]))

fA.write(sA)
fB.write(sB)

fA.close()
fB.close()