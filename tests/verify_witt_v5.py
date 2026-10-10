import ast,contextlib,io,itertools
from pathlib import Path
p=str(Path(__file__).resolve().parents[1] / 'Witt_Basis_Sum_Quadratics.py')
s=open(p).read(); s=s[:s.index('\nwhile True:\n\n    construct_witt_basis()')]
ns={}; exec(compile(s,p,'exec'),ns)
classes=['epsilon','pi','epsilon*pi']; results=[]
for minus, a,b in itertools.product([True,False],classes,classes):
 for m1,m2 in itertools.product(['1','pi'] if a=='epsilon' else ['1','epsilon'], ['1','pi'] if b=='epsilon' else ['1','epsilon']):
  ns['minus_one_square']=minus
  h=ns['hilbert']; valid=h(m1,a)*h(m2,b)==h(a,b)
  if not valid:continue
  case='4.1' if a==b else '4.2' if h(a,b)==1 else '4.3'
  # intercept input sequence
  sel=lambda v,cl: str(cl.index(v)+1)
  it=iter(['yes' if minus else 'no',sel(a,classes),sel(b,classes),sel(m1,['1','pi'] if a=='epsilon' else ['1','epsilon']),sel(m2,['1','pi'] if b=='epsilon' else ['1','epsilon'])])
  ns['input']=lambda prompt='':next(it)
  out=io.StringIO()
  with contextlib.redirect_stdout(out): ns['construct_witt_basis']()
  t=out.getvalue(); results.append((minus,a,b,m1,m2,case,t))
from collections import Counter
print('VALID CASE COUNTS',Counter((x[0],x[5]) for x in results))
print('NO BASIS OUTPUT:')
for row in results:
 if 'e1 =' not in row[-1]:print(row[:6])
print('SYMBOLIC WITNESSES:')
for row in results:
 if 'Symbolic witnesses' in row[-1] or 'Explicit norm witness' not in row[-1] and row[5]=='4.1':print(row[:6], [v for v in row[-1].splitlines() if v.startswith(('r =','s =','u ='))])
print('CASE 4.3 ALL:')
for row in results:
 if row[5]=='4.3':print(row[:5], [v for v in row[-1].splitlines() if v.startswith(('r =','s =','e2 =','f2 ='))])
