import subprocess
name = "submit_corr2D_trg_ref.jdl"
f = open(name, "w")

command_lines = '''Universe   = vanilla
GetEnv     = True
Executable = submit_corr2D_trg_ref.sh
Arguments  = small.list.00
Log        = log/submit_corr2D_trg_ref.$(Process).log
Output     = out/submit_corr2D_trg_ref.$(Process).out
Error      = err/submit_corr2D_trg_ref.$(Process).err
Queue
'''

for i in range(1, 80):
   temp = '''
Arguments  = small.list.%02d
Queue
   ''' % i
   command_lines += temp

f.write(command_lines)
f.close()
subprocess.call(["condor_submit", name]);
