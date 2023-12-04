#!/usr/bin/python
import os
import shutil

f=open('exercise4.py','r')                                               #readtemplate
template=f.read()
f.close


#Loop over the atoms number 
for t in [1, 2, 3, 4, 5, 6, 7, 8]:
   t=str(t)                                                     # float to string
   print(t)
   if not os.path.exists(t): os.mkdir(t)                        # Create directory with the xc as name
   shutil.copy('exercise4.py', t+'/exercise4.py')                    # copy the template file into the new directoy
   out=open(t+'/exercise4.py','w')                                      # replace the term xc in temaplte by the current xc
   out.write(template.replace('KPOINT', t))                               #   and copy the new file into the directory
   out.close()

   os.chdir(t)

   shutil.copyfile("../" + 'Si.vbc.UPF',  "Si.vbc.UPF")

   os.system("python3 exercise4.py | tee output.log")

   os.chdir('..')
