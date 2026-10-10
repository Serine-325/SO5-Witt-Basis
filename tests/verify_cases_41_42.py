import sympy as S
from pathlib import Path
import contextlib,io,itertools
source=(Path(__file__).resolve().parents[1] / 'Witt_Basis_Sum_Quadratics.py').read_text().split('\nwhile True:\n\n    construct_witt_basis()')[0]
ns={}
with contextlib.redirect_stdout(io.StringIO()):exec(compile(source,'v5','exec'),ns)
classes=['epsilon','pi','epsilon*pi']; rows=[]
def norm(z,d):
 # pair (a,b) represents a+b sqrt(d)
 a,b=z;return S.expand(a*a-d*b*b)
def pair(a,b,d):return S.expand((a[0]*b[0]-d*a[1]*b[1]))
def q(v,m1,m2,a,b,l):return S.expand(m1*norm(v[0],a)+m2*norm(v[1],b)+l*v[2]**2)
def B(v,w,m1,m2,a,b,l):return S.expand(m1*pair(v[0],w[0],a)+m2*pair(v[1],w[1],b)+l*v[2]*w[2])
def check(vec,m1,m2,a,b,l,rels):
 gb=S.groebner(rels,*sorted(set().union(*(z.free_symbols for z in rels)),key=str),domain=S.QQ) if rels else None
 def red(z):
  z=S.cancel(z);num,den=S.fraction(z)
  if gb:
   num=gb.reduce(S.expand(num))[1]
  return S.factor(num/den)
 e1,e2,f1,f2,g=vec
 checks={f'q({i})':red(q(v,m1,m2,a,b,l)-(1 if i=='g' else 0)) for i,v in [('e1',e1),('e2',e2),('f1',f1),('f2',f2),('g',g)]}
 for i,j in itertools.combinations(range(5),2):
  checks[f'B({i},{j})']=red(B(vec[i],vec[j],m1,m2,a,b,l)-(1 if (i,j) in ((0,2),(1,3)) else 0))
 return {k:v for k,v in checks.items() if v!=0}
for sq,a,b in itertools.product((True,False),classes,classes):
 ns['minus_one_square']=sq
 for m1s,m2s in itertools.product(['1','pi'] if a=='epsilon' else ['1','epsilon'], ['1','pi'] if b=='epsilon' else ['1','epsilon']):
  h=ns['hilbert'];
  if h(m1s,a)*h(m2s,b)!=h(a,b):continue
  case='4.1' if a==b else '4.2' if h(a,b)==1 else '4.3'
  if case=='4.3':continue
  eps=S.Integer(2 if sq else -1);pi=S.Integer(3)
  cv=lambda x:{'1':S.Integer(1),'epsilon':eps,'pi':pi,'epsilon*pi':eps*pi}[x]
  aa,bb,m1,m2,l=map(cv,(a,b,m1s,m2s,ns['multiply_classes'](a,b)))
  inputs=iter(['yes' if sq else 'no',str(classes.index(a)+1),str(classes.index(b)+1),str((['1','pi'] if a=='epsilon' else ['1','epsilon']).index(m1s)+1),str((['1','pi'] if b=='epsilon' else ['1','epsilon']).index(m2s)+1)])
  ns['input']=lambda prompt='':next(inputs)
  out=io.StringIO()
  with contextlib.redirect_stdout(out):ns['construct_witt_basis']()
  txt=out.getvalue();key=(sq,a,b,m1s,m2s)
  if case=='4.1':
   u0,u1=S.symbols('u0 u1');u=(u0,u1);nu=norm(u,aa)
   v=[((1,0),u,0),((0,1),(aa*u1,u0),0),((S.Rational(1,2)/m1,0),(-u0/(2*m1),-u1/(2*m1)),0),((0,-1/(2*aa*m1)),(u1/(2*m1),u0/(2*aa*m1)),0),((0,0),(0,0),1)]
   # e2 second sqrt(a)*u=(a*u1,u0); f2 second sqrt(a)*u/(2*a*m1)
   # fix f2 second coefficient u0/(2*a*m1)
   rel=[nu+m1/m2]
   bad=check(v,m1,m2,aa,bb,l,rel)
   witness=ns['case_41_u'](a,m1s,m2s)[0]
   witness_ok=False
   if witness=='i':witness_ok= S.simplify(-m1/m2+1)==0
   elif witness=='1':witness_ok= S.simplify(-m1/m2-1)==0
   elif witness=='s^2 + i*t^2':witness_ok= S.simplify(-m1/m2+1)==0 and aa==-1
   else:witness_ok=False
   if not witness_ok:bad['witness']=f'{witness} does not have norm {-m1/m2}'
  else:
   if m1s==m2s=='1':
    # explicit bases as printed: beta=-alpha since eps=-1
    # first pair (1,0,1),(1/2,0,-1/2)
    # second sqrt(a),sqrt(b) and f2=(-1/(2sqrt(a)),sqrt(b)/(2a))
    # 1/sqrt(a)=sqrt(a)/a
    v=[((1,0),(0,0),1),((0,1),(0,1),0),((S.Rational(1,2),0),(0,0),-S.Rational(1,2)),((0,-1/(2*aa)),(0,1/(2*aa)),0),((0,0),(1,0),0)]
    bad=check(v,m1,m2,aa,bb,l,[])
   else:
    s,t=S.symbols('s t');rel=[s**4+t**4+1]
    # use i sqrt(pi) = sqrt(-pi) = sqrt(epsilon*pi), and other reversed
    if a=='pi':
     v=[((0,1),(0,1),0),((t**2,0),(-s**2,0),1),((0,1/(2*pi)),(0,-1/(2*pi)),0),((t**2/2,0),(-s**2/2,0),-S.Rational(1,2)),((s**2,0),(t**2,0),0)]
    else:
     v=[((0,1),(0,1),0),((-s**2,0),(t**2,0),1),((0,-1/(2*pi)),(0,1/(2*pi)),0),((-s**2/2,0),(t**2/2,0),-S.Rational(1,2)),((t**2,0),(s**2,0),0)]
    bad=check(v,m1,m2,aa,bb,l,rel)
  rows.append((case,key,bad))
from collections import Counter
print('TESTED',dict(Counter(x[0] for x in rows)))
print('PASS',sum(not x[2] for x in rows),'FAIL',sum(bool(x[2]) for x in rows))
for case,key,bad in rows:print('PASS' if not bad else 'FAIL',case,key, bad)
