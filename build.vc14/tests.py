
from __future__ import print_function
import os
import sys
import string
import platform
from re import match
from subprocess import Popen, PIPE, STDOUT
from tempfile import *

try:
  from time import perf_counter as timer
except ImportError:
  from time import clock as timer
  

x64 = True
debug = False
test_gpu_version = True
run_non_gpu_tests = True
run_gpu_tests = True

class Timer() :
  def __enter__(self): self.start = timer()
  def __exit__(self, *args): print(' time {:.3f} milliseconds'.format(1000 * (timer() - self.start)))

cpath = os.path.dirname(__file__)
config = 'x64' if x64 else 'Win32'
mode = 'Debug' if debug else 'Release'
test_dir = '..\\bin\\{:s}\\{:s}\\'.format(config, mode)

def get_tests(filename):
  print('running tests in {:s}'.format(filename))
  start, sub, tests, c_tests = True, dict(), [], []
  with open(os.path.join(cpath, filename)) as f:
    lines = f.readlines()
    cnt, lnth = 0, len(lines)
    while cnt < lnth:
      try:
        line = lines[cnt].strip()
        cnt += 1
        tkns = line.split()
        if line.startswith('echo') and len(tkns) > 2 and tkns[2] == '|':
          while cnt < lnth and 'checkcode' not in line:
            while cnt < lnth and not lines[cnt]:
              cnt += 1
            if cnt < lnth:
              line += '|' + lines[cnt]
              cnt += 1
          start = False
        elif start:
          sp = line.split('="')
          if len(sp) == 2:
            if sp[1].startswith('${1:-./ecm}'):
              sub[sp[0]] = sp[1][12:-1]
            else:
              sub[sp[0]] = sp[1][:-1]
          continue
        else:
          continue
        line = line.replace(';', '|')
        sub_tests = []
        for line_part in line.split('|'):
          tkns = line_part.strip().split()
          cmd =  []
          for tok in tkns:
            if tok.startswith('"') and tok.endswith('"'):
              tok = tok[1:-1]
            if tok[0] == '$' and tok[1:] in sub:
              tok = tok.replace(tok, sub[tok[1:]])
            elif tok == './ecm':
              tok = ''
            cmd += [tok]
          cseq = []
          if cmd and cmd[0] == 'echo':
            cseq += [cmd[1]]
            cmd = cmd[2:]
          if len(cmd) >= 3 and cmd[-3] == 'checkcode' and cmd[-2] == '$?':
            cseq += [int(cmd[-1])]
            cmd = cmd[:-3]
          cmd = (' '.join(cmd)).strip()
          if cmd:
            cseq += [cmd]
          sub_tests += [cseq]
        if len(sub_tests) == 3 and all(len(x) == 1 for x in sub_tests):
          tests += [tuple(x[0] for x in sub_tests)]
        else:
          c_tests += [sub_tests]
      except ValueError:
        print('parsing error on line {} in text "{}"'.format(cnt, line))
  return tests, c_tests

def run_exe(exe, args, inp) :
  al = {'stdin' : PIPE, 'stdout' : PIPE, 'stderr' : STDOUT }
  if sys.platform.startswith('win'):
    al['creationflags'] = 0x08000000
  p = Popen([os.path.join(cpath, exe)] + args.split(' '), **al)
  res = p.communicate(inp.encode())[0].decode()
  ret = p.poll()
  return (ret, res)

def output_complex_tests(x):
  print('these tests are too complex:')
  for t in x:
    print(t)

def do_tests(tests, ctests, out=False, gpu=False):
  ecm_exe  = test_dir + ("ecm_gpu.exe" if gpu else "ecm.exe")
  err_cnt = 0
  for ix, tt in enumerate(tests):
    print(tt[1], tt[0], end='')
    rv = run_exe(ecm_exe, tt[1], tt[0])
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

  if ctests:
    output_complex_tests(ctests)
  if not err_cnt:
    if ctests:
      print('all other tests passed')
    else:
      print('all tests passed')

with Timer():
  if os.path.exists('test.pm1.save'):
    os.remove('test.pm1.save')
  if run_non_gpu_tests:
    t, ct = get_tests("..\\test.ecm")
    do_tests(t, ct)
    t, ct = get_tests("..\\test.pm1")
    do_tests(t, ct)
    t, ct = get_tests("..\\test.pp1")
    do_tests(t, ct)
    t, ct = get_tests("..\\testlong.pp1")
    do_tests(t, ct)
    t, ct = get_tests("..\\testlong.pm1")
    do_tests(t, ct)
  if run_gpu_tests:
    t, ct = get_tests("..\\test.gpuecm")
    do_tests(t, ct, gpu=True)
