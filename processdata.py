import glob
import statistics
import pandas as pd
import matplotlib.pyplot as plt

dirs = []
for i in range(1, 51):
    dirs.append(glob.glob("experimentos/exp"+str(i)+"/*.time"))

files = {}
for _dir in dirs:
    for exp in _dir:
        name = exp.split('/')[-1]
        with open(exp, "r") as f:
            time = float(f.read().strip())
        if name in files.keys(): 
            files[name].append(time)
        else:
            files[name] = [time] 

dict(sorted(files.items()))
table = {}
for f in dict(sorted(files.items())):
    # print("%s avg: %f s stdev: %f s" % (f, statistics.mean(files[f]), statistics.stdev(files[f])))
    mean = statistics.mean(files[f])
    stdev = statistics.stdev(files[f])
    if(f in table.keys()):
        table[f].append(mean)
        table[f].append(stdev)
    else:
        table[f] = [mean]
        table[f].append(stdev)

df = pd.DataFrame(table, index = ["Média", "Desvio Padrão"])
df = df.transpose()
print(df)
df.to_csv("results.csv")
