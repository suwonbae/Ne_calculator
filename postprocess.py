import numpy
import subprocess

#change num_chains and num_atoms_per_mol depending on model
num_chains=27
num_atoms_per_mol=2000

lee=numpy.zeros(num_chains)
lpp=numpy.zeros(num_chains)
for i in range(num_chains):
   subprocess.Popen("tail -2 stdout."+str(i+1)+" | head -1 > temp",shell=True).wait()
   result=numpy.loadtxt('temp')
   lee[i]=result[0]**2
   lpp[i]=result[1]

subprocess.Popen("rm temp",shell=True).wait()

bpp=numpy.mean(lpp)/(num_atoms_per_mol-1)
app=numpy.mean(lee)/numpy.mean(lpp)
Ne=app/bpp

print lpp
print Ne
