import os,glob
cpp_files = glob.glob('*.cpp')
for f in cpp_files:
  comp_str = 'g++ %s -o %s' %(f, f.split('.cpp')[0])
  print('compiling: %s' %comp_str)
  os.system(comp_str)
