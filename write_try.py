import numpy as np

fn="test.txt"

f=open(fn,"w")

a=1

f.write(str(a)+'\n')
f.write(str(a))
f.write(str(a))
f.write(str(a))
f.write(str(a))
f.write(str(a))
f.close()
