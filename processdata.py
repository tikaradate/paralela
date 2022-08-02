import glob
import statistics

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
for f in dict(sorted(files.items())):
    print("%s avg: %f s stdev: %f s" % (f, statistics.mean(files[f]), statistics.stdev(files[f])))
