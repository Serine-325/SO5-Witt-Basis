from pathlib import Path
import contextlib,io,itertools,re,collections
import sympy as S
source=open(Path(__file__).resolve().parents[1] / 'Witt_Basis_Sum_Quadratics.py').read().split('\nwhile True:\n\n    construct_witt_basis()')[0]
ns={};exec(compile(source,'v5','exec'),ns)
U,V=S.symbols('U V'); s0,t0=S.symbols('s0 t0'); classes=['epsilon','pi','epsilon*pi']
counts=collections.Counter(); failures=[]
for square,a,b in itertools.product([True,False],classes,classes):
 for m1,m2 in itertools.product(['1','pi'] if a=='epsilon' else ['1','epsilon'], ['1','pi'] if b=='epsilon' else ['1','epsilon']):
  ns['minus_one_square']=square;h=ns['hilbert']
  if h(m1,a)*h(m2,b)!=h(a,b):continue
  case='4.1' if a==b else '4.2' if h(a,b)==1 else '4.3'
  it=iter(['yes' if square else 'no',str(classes.index(a)+1),str(classes.index(b)+1),str((['1','pi'] if a=='epsilon' else ['1','epsilon']).index(m1)+1),str((['1','pi'] if b=='epsilon' else ['1','epsilon']).index(m2)+1)])
  ns['input']=lambda prompt='':next(it)
  out=io.StringIO()
  with contextlib.redirect_stdout(out):ns['construct_witt_basis']()
  lines=out.getvalue().splitlines(); vals={}
  for line in lines:
   for k in ['e1','e2','f1','f2','g1']:
    if line.startswith(k+' = '):vals[k]=line.split(' = ',1)[1]
  if case!='4.3': continue
  key=(square,a,b,m1,m2)
  if len(vals)!=5:failures.append((key,'missing vectors'));continue
  eps=2 if square else -1; pi=3
  cv=lambda x: {'1':1,'epsilon':eps,'pi':pi,'epsilon*pi':eps*pi}[x]
  alpha,beta,mu1,mu2=cv(a),cv(b),cv(m1),cv(m2);lam=cv(ns['multiply_classes'](a,b))
  def root(z):
   z=S.sympify(z)
   if z==alpha:return U
   if z==beta:return V
   return S.sqrt(z)
  i = S.I if square else (U if a=='epsilon' else V if b=='epsilon' else S.I)
  loc={'epsilon':S.Integer(eps),'pi':S.Integer(pi),'sqrt':root,'i':i,'s0':s0,'t0':t0,'u':S.Symbol('u'),'r':S.Symbol('r'),'s':S.Symbol('s')}
  def parse(expr):
   expr=expr.strip().replace('^','**')
   if expr.startswith('(') and expr.endswith(')'):
    expr=expr[1:-1]
   parts=expr.split(',')
   if len(parts)!=3:raise ValueError(expr)
   return tuple(S.sympify(p,locals=loc) for p in parts)
  try: vectors={k:parse(v) for k,v in vals.items()}
  except Exception as e:failures.append((key,'parse '+str(e)));continue
  def red(z):
   z=S.cancel(z)
   numer,den=S.fraction(z)
   if den.has(U,V,s0,t0):return S.simplify(z)
   numer=S.Poly(S.expand(numer),U,V,s0,t0)
   gb=S.groebner([U**2-alpha,V**2-beta]+([] if square else [s0**4+t0**4+1]),U,V,s0,t0,order='lex',extension=S.I)
   return S.factor(gb.reduce(numer.as_expr())[1]/den)
  def conj(z,w):return red(z.subs(w,-w, simultaneous=True))
  def norm(z,w):return red(z*conj(z,w))
  def bil(z,w,gen):return red((z*conj(w,gen)+w*conj(z,gen))/2)
  def q(v):return red(mu1*norm(v[0],U)+mu2*norm(v[1],V)+lam*v[2]**2)
  def B(x,y):return red(mu1*bil(x[0],y[0],U)+mu2*bil(x[1],y[1],V)+lam*x[2]*y[2])
  tests={}
  for k in ['e1','e2','f1','f2']:tests['q('+k+')']=q(vectors[k])
  tests['q(g1)-1']=red(q(vectors['g1'])-1)
  for x,y in itertools.combinations(['e1','e2','f1','f2','g1'],2):
   tests['B('+x+','+y+')']=red(B(vectors[x],vectors[y])-(1 if (x,y) in [('e1','f1'),('e2','f2')] else 0))
  bad={k:str(v) for k,v in tests.items() if v!=0}
  counts[(square,case)]+=1
  if bad:failures.append((key,bad))
print('COUNTS',dict(counts));print('TOTAL',sum(counts.values()),'FAILURES',len(failures))
for f in failures:print('FAIL',f)
