import subprocess
name = "submit_c24_d0.jdl"
f = open(name, "w")

command_lines = '''Universe   = vanilla
GetEnv     = True
Executable = submit_c24_d0.sh
Arguments  = 00
Log        = log/submit_c24_d0.$(Process).log
Output     = out/submit_c24_d0.$(Process).out
Error      = err/submit_c24_d0.$(Process).err
Queue
'''

for i in range(1, 11):
   temp = '''
Arguments  = %02d
Queue
   ''' % i
   command_lines += temp

f.write(command_lines)
f.close()
subprocess.call(["condor_submit", name]);
