
load 'csv dates tables/tara'

XL=:readxlworkbook

fdir=:'/users/pietdejong/documents/research/cifr/'
load fdir,'programs/kfs.ijs'
dir=:'~/documents/research/cifr/'

Phi=:pnorm :. qnorm

rddat=: 3 : 0
  dat=:readcsv dir,'/data/cifrdatdaily.csv'
  ({.dat)=:|:x=.>".L:0  '-_'&charsub L:0 }.dat
  date=:1&todayno date
  {.dat
)

NB. rddat''

phi=:(%&(%:&o.2))@^@-@-:@*:     NB.  standard normal density
rhot=:%@%:@>:@%@*:              NB.  

cusph=:3 : 0  NB. cubic spline model with forcing x
   c=.3-%:3['z q s'=.y
   h=:(%mean) 1: P q
   start=.(_;_ _);(0 1 1,c),:0 0 1 1
   D_t=.(1 1 0 0 0),(0 1 1 0,c),:0 0 1 0 1
   D=:start,z;"0 2 (q,.s*h) (0 0;0 3)}"1 2 D_t
)

figx=: 3 : 0
  pd 'reset;sub ',":p,p['p n'=.$u=.P"1>y
  pd 'show'[res=:panel"1/~u 
  header=:<'    b      a0      d0      rho ',:' t-stat  t-stat  t-stat %:1-*:ssig'
  header,res
)

panel=: 4 : 0
  pd 'new;pensize 1;yticpos -3 -2 -1 0 1 2 3;xticpos -3 -2 -1 0 1 2 3'
  if. x-:y do. <'' return. end. 
  rho=.mean*/zsc"1 'q z'=.Phi^:_1 'u v'=.(x,:y)/:"1 x
  'a d'=.|:ad['ad covad int covint'=.|:SMI cusph z;q;s=.^9
  'b a0 d0'=.bad[tstat=.bad%sig*c['bad sig c'=.beta_KFR;sig_KFR;|getd C_KFR
  pd L:_1 (<q;z);~'type dot;color black;'
  pd L:_1 (<q;rho*q);~'type line;color red'
  pd L:_1 (<q;b*q);~'type line;color green'
  pd L:_1 (<q;a+b*q);~'type line;color black'
  rho_t=.rhot (psi=.b+(#d)*(phi q)*d)%s*sig*h
  pd L:_1 (<q;rho_t);~'type line;color purple'
NB.  rho_t=.rhot s*sig*h%tau=.b+((*#)d)*phi q
  <8j4":L:0 (bad,:tstat),.rho,%:1-*:s*sig
)


figx anz;cba;nab;wbc
stop

simnorm=: 3 : 0
  'n rho'=.y
  (chol 2 2 $ 1,rho,rho,1) mp rnorm 2,n
)

fcop=: 3 : 0
  chisq=.+/*:rnorm 'df n'=.y
  chisq,:eps=.chisq+1.5*rnorm n
)

dplot=: 'dot;pensize 3'&plot
pcop=:dplot@(P L:0)  NB. Plot copula
pncop=:dplot@(qnorm@:P L:0)

cuspf=:3 : 0  NB. cubic spline model with forcing x
   h=.3-%:3['y x s'=.y 
   start=.(_;_ _);(0 1 1,h),:0 0 1 1
   D_t=.(1 1 0,s,0),(0 1 1 0,h),:0 0 1 0 1
   D=:start,y;"0 2 x (<0 0)}"0 2 D_t
)

cusp=:3 : 0  NB. cubic spline model
   h=.3-%:3['y s'=.y
   start=.('';_ _);(1 1,h),:0 1 1
   D_t=.(1 0,s,0),(1 1 0,h),:0 1 0 1
   D=:start,y;"0 2 D_t
)

sreg=: 3 : 0
  rho=.{:{.cor|:'q z'=.|: y
  'muq muz sdq sdz'=.(mean,sd) y
  CoV=.sdq%muq
  b=.rho*sdz%sdq
  a=.muz-b*muq
  a,b,rho,CoV
)


treg=: 3 : 0
 'x y'=.P"1 y
 'q z'=.(-@-:@# {. ])"1 Phi^:_1 'u v'=.(x,:y)/:"1 x
  'a b rho CoV'=.|:sreg\.q,.z
   plot   _40}."1 rho,a,:b*q
)

simnorm=: 3 : 0
  'n rho'=.y
  (chol 2 2 $ 1,rho,rho,1) mp rnorm 2,n
)




treg fcop 1 1000
  
NB. treg fcop 1 1000
  
stop


R=:Phi^:_1@P

fig1=: 3 : 0
  'a b'=.y
  q=.Phi^:_1 u=./:~ P runif n=.1000
  ssig=.%:1-(3%~*:a)+*:b
  z=.(Ez=.(a*1-+:u)+(b*q))+ssig*rnorm#q
  wrs b,{:{.cor z,.q
  pd 'reset;sub 3 4;pensize 1'
  pd L:_2 ('type dot';(<q;z));('type line';(<q;Ez))
  pd L:_2 'new';('type line';(<u;Ez));('type dot';(<u;z))
  pd L:_1 'new;type dot';(<u;v=.Phi z)
  pd L:_1 'new;type dot';(<(P z); v)
  pd 'show'
)

NB.fig1 0 0.8
NB. stop


sim=: 3 : 0   NB.  draw randomly from empirical copula
  n=.1000
  u=.pnorm ({., {.+ 20*{:)"1 rnorm n,2
  chi=. (_1++:?n#2)*  %:+/"1 *:rnorm n,3$u
  chi=.chi,.(_1++:?n#2)*  %:+/"1 *:rnorm n,3
  dplot ;/|:pnorm chi*<:+:u
)





finds=: 4 : 0
  'q z'=.Phi^:_1 'u v'=.(P"1 x,:y)/:"1 x
  {:@>{:@KFR@KF@cuspm"1 z;q;s=.^8
)
   
NB. cba finds anz 



fitpair=: 4 : 0
  'q z'=.qnorm 'u v'=.|:(|:/:{.) P"1 x,:y
  sm=.smooth@((z;q)&(,<))
NB.  for_i. 7+-:-:-:i.10 do.
NB.    zhat=.sm wrs s=.^i
NB.    wrs 10j3": ell_KFR,(%>:*:s*sig_KFR%beta),beta=.{.beta_KFR end.
  pd 'reset;new;type dot;pensize 2'
  pd u;z
  pd 'type line'
  pd  u;mu=:sm 1600
  pd 'show'
)

anz fitpair *:cba

stop



  

stop
  



C=: 3 : 0
  'x y p q'=.y
  'u v'=. P L:0 x;y
  mean (u<:p)*.(v<:q)
)

br=: 4 : 0
  h=:x
  x=:/:~rnorm 1000
  Fx=:+/\P x
  fx=: ((h}.Fx) - (-h)}.Fx)%(h}.x)-(-h)}.x
)
   

dydx=:}:@{.,:((%~/)@:(dif"1))  
F=:(,:P)@:(/:~)
f=:dydx@F
If=:({.,:{.*{:)@:f
G=:{.@f ,:(}:@{:@F - {:@If)


cop=: 3 : 0
  'u v'=.P L:0  [ 5&{. L:0 cba;anz
NB.  suv=.u (,"0/)&(/:~) v
NB. c=.(|.+/(u,.v) (*/@:<:)"1"1 3 suv)%#u
  'su sv'=./:~ L:0 u;v
  'm1 m2'=.(u;v) (<:"0 1)L:0 (su;sv) 
  
) 

update=: 4 : 0
   'q x y '=.x [ qij=.y
   (C x;y;(q+qij);q)-*:q
)

s=: 4 : 0
   q=.x
   qij=.((q;y)&update)^:_[0
   1<.qij%q*1-q
)

S=: 4 : 0        NB. Sensitivity matrix:   x=q , y=boxed list of series
   (x&s)@,"0/~y   
)
  

NB. F=:(;P)@/:~                NB.  x;F(x) from data series 
NB. f=:%~/@:>@(dif L:0)@F      NB.  f(x) from data series 
u=:>:@i. % >:              NB.  uniform grid or quantiles
ma=: 100&(mean\)           NB.  moving average

rf=: 3  : 0                NB. estimate geom mean of ratio
   u=.ma"1 u <:#>{. y  
   u; %/ ^ ma"1 ^. >  f L:0 y
)
   
panel1=: 3 : 0
   pd 'new;pensize 1'
   if. ({.-:{:)y do.
      pd 'type line' 
      pd ({.;ln@(1:-{:))>F>{.y
   else. 
      pd 'type dot'
      pd u=.P L:0 y
      u=./:~>{.u
      pd 'type line;color red'
      pd u;u s"0 2 y end.   
)


fig1=: 3 : 0
   wrs p=.":2#$y
   pd 'reset;sub ',p
   (panel1@:,)"0/~y
   pd 'show'
   pd 'pdf 320 300 ',fdir,'fig1.pdf' 
)

NB. fig1 cba;anz;mqg;wbc;banks

panel2=: 3 : 0
   pd 'new;pensize 1'
   if. ({.-:{:)y do.
      pd 'type line' 
      pd ({.;ln@(1:-{:))>rf>{.y
   else. 
      pd 'type dot'
      pd u=.P L:0 y
      u=./:~>{.u
      pd 'type line;color red'
      pd u;u s"0 2 y end.   
)


fig2=: 3 : 0
   p=.":2#$y
   pd 'reset;sub ',p
   (panel2@:,)"0/~y
   pd 'show'
   pd 'pdf 320 300 ',fdir,'fig2.pdf' 
)

NB.  fig2 wbc;anz

sens=: 3 : 0
  q=.0.8+(0.95-0.8)*u 4 
  q S"0 1 y
)

NB. sens anz;wbc;mqg;cba
  
  