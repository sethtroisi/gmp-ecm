
from __future__ import print_function
import os
import sys
import string
import platform
from re import match
from subprocess import Popen, PIPE, STDOUT
from tempfile import *
from time import clock

class Timer() :
  def __enter__(self): self.start = clock()
  def __exit__(self, *args): print(' time {:.3f} milliseconds'.format(1000 * (clock() - self.start)))

test_dir = '..\\bin\\x64\\Release\\'
# test_dir = '..\\bin\\win32\\Release\\'
test_normal = True
test_gpu = True

def get_tests(filename):
  print('running tests in {:s}'.format(filename))
  start, sub, tests = True, dict(), []
  with open(filename) as f:
    lines = f.readlines()
    cnt, lnth = 0, len(lines)
    while cnt < lnth:
      line = lines[cnt].strip()
      cnt += 1
      tkns = line.split()
      if line.startswith('echo') and len(tkns) > 2 and tkns[2] == '|':
        if 'checkcode' not in line:
          while not lines[cnt] and cnt < lnth:
            cnt += 1 
          if cnt < lnth and 'checkcode' in lines[cnt]:
            line += '; ' + lines[cnt]
            cnt += 1
        start = False
      elif start:
        sp = line.split('="')
        if len(sp) == 2:
          if sp[1].startswith('$1'):
            sub[sp[0]] = sp[1][2:-1]
          else:
            sub[sp[0]] = sp[1][:-1]
        continue
      else:
        continue
      tkns = line.split()
      cmd =  ''
      for tok in tkns[3:]:
        if tok.startswith('"') and tok.endswith('"'):
          tok = tok[1:-1]
        if tok[0] == '$' and tok[1:] in sub:
          tok = tok.replace(tok, sub[tok[1:]])
        if tok.endswith(';'):
          cmd += tok[:-1]
          break
        else:
          cmd += tok + ' '
      code = tkns[-1]
      t = (tkns[1].replace('\"', ''), cmd.strip(), int(tkns[-1]))
      tests += [t]
  return tests

def run_exe(exe, args, inp) :
  al = {'stdin' : PIPE, 'stdout' : PIPE, 'stderr' : STDOUT }
  if sys.platform.startswith('win'):
    al['creationflags'] = 0x08000000
  p = Popen([exe] + args.split(' '), **al)
  res = p.communicate(inp.encode())[0].decode()
  ret = p.poll()
  return (ret, res)

def do_tests(tests, gpu=False, out=False):
  exe  = test_dir + ("ecm_gpu.exe" if gpu else "ecm.exe")
  err_cnt = 0
  for ix, tt in enumerate(tests):
    print(tt[1], tt[0], end='')
    rv = run_exe(exe, tt[1], tt[0])
    if type(tt[2]) == int and rv[0] != tt[2]:
      print(" - *** ERROR in test {:d}: {:d} {:d} ***".format(ix, rv[0], tt[2]))
      err_cnt += 1
    elif type(tt[2]) == tuple and rv[0] != tt[2][0] and rv[0] != tt[2][1]:
      print(" - *** ERROR in test {:d}: {:d} {:s} ***".format(ix, rv[0], tt[2]))
      err_cnt += 1
    else:
      print(" - passed")
    if out:
      op = rv[1].rsplit('\r\n')
      for i in op :
        print(i)

  if not err_cnt:
    print('  all tests passed')

with Timer():
  if os.path.exists('test.pm1.save'):
    os.remove('test.pm1.save')
  if test_normal:
    do_tests(get_tests("..\\test.ecm"))
    do_tests(get_tests("..\\test.pm1"))
    do_tests(get_tests("..\\test.pp1"))
    do_tests(get_tests("..\\testlong.pp1"))
  if test_gpu:
    do_tests(get_tests("..\\test.gpuecm"), gpu=True)
