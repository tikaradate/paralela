# uso: python3 genseq.py (quantidade de caracteres a serem gerados)

import sys, random, string

def random_char(y):
    return ''.join(random.choice(string.ascii_letters) for x in range(y))

f = open(sys.argv[1] + ".in", "x") 

s = random_char(int(sys.argv[1]))

f.write(s)

f.close()