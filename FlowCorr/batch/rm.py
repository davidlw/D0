import subprocess
for i in range(0, 74):
   id = "21.%d" % i
   print id
   subprocess.call(["condor_rm", id])
