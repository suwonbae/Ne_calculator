import numpy
import subprocess

lee=numpy.zeros(27)
lpp=numpy.zeros(27)
for i in range(27):
   subprocess.Popen("tail -2 stdout."+str(i+1)+" | head -1 > temp",shell=True).wait()
   result=numpy.loadtxt('temp')
   lee[i]=result[0]**2
   lpp[i]=result[1]

bpp=numpy.mean(lpp)/(2000-1)
app=numpy.mean(lee)/numpy.mean(lpp)
Ne=app/bpp

print lpp
print Ne
