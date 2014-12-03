
load 'csv dates tables/tara'

fdir=:'/users/pietdejong/documents/research/cifr/'
load fdir,'programs/kfs.ijs'
dir=:'~/documents/research/cifr/'

Phi=:pnorm :. qnorm

rddat=: 3 : 0
  dat=:readcsv dir,'/data/cifrdat.csv'
  ({.dat)=:|:x=.>".L:0  '-_'&charsub L:0 }.dat
  date=:1&todayno date
  {.dat
)

NB. rddat''

dplot=: 'dot;pensize 3'&plot
pcop=:dplot@(P L:0)  NB. Plot copula
pncop=:dplot@(qnorm@:P L:0)

sim=: 3 : 0   NB.  draw randomly from empirical copula
  n=.1000
  u=.pnorm ({., {.+ 20*{:)"1 rnorm n,2
  chi=. (_1++:?n#2)*  %:+/"1 *:rnorm n,3$u
  chi=.chi,.(_1++:?n#2)*  %:+/"1 *:rnorm n,3
  dplot ;/|:pnorm chi*<:+:u
)




   
cuspm=:3 : 0  NB. cubic spline model
   'y x s'=.y
   y_0=.(_;_ _) [ D_0=.(0 1 1,h=.(-%:)3),:0 0 1 1
   D_t=.(0 1 0,s,0),(0 1 1 0,h),:0 0 1 0 1
   D=:(y_0;D_0),y;"0 2 x (<0 0)}"0 2 D_t
)

muest=:,@:>@:{:@:|:@:(]ESM KF)  NB.  'eps coveps vareps covvareps mu'=.|:(]ESM KF) D


figure=: 3 : 0
  pd 'reset;sub ',":p,p['p n'=.$u=.P"1>y
  pd 'show'[res=:panel"1/~u 
  res,~<'   beta    alpha0   delta0   R',:'     t        t       t    ssig'
)

panel=: 4 : 0
  'q z'=.Phi^:_1 'u v'=.(x,:y)/:"1 x
  pd 'new;pensize 1;type dot;color black'
  pd (u;z)
  pd 'type line;color red' NB.; key R2=',6j2": (1-*:r*sig_KFR),*:sig_KFR  
  pd u;zhat=.muest cuspm z;q;s=.^12
  'beta alpha0 delta0'=.ba['ba sig c'=.beta_KFR;sig_KFR;|getd C_KFR
  tstat=.ba%sig*c[R=.%:1-*:ssig=.s*sig
  results=.(ba,:tstat),.R,ssig
NB.  pd 'type line;color green'
NB.  pd 2&(mean\)"1 L:0 u;1,:*:(z-zhat)%ssig
  'u vhat'=.|:(u,.Phi zhat)/:/:x  NB. return to original order 
  <8j4":L:0 results
)

try=: 3 : 0
  'a b'=.y
  c=.0.283
  q=.Phi^:_1 u=./:~ P runif n=.1000
  ssig=.%:1-(3%~*:a)+(*:b)-+:a*b*c
  z=.(Ez=.(a*(1-+:u))+(b*q))+ssig*rnorm#q
  wrs b,{:{.cor z,.q
  pd 'reset;sub 2 2;pensize 1'
  pd L:_2 ('type dot';(<u;z));('type line';(<u;Ez))
  pd L:_2 'new';('type dot';(<q;z));('type line';(<q;Ez))
  pd L:_1 'new;type dot';(<u;v=.Phi z)
  pd L:_1 'new;type dot';(<(P z);v)
  pd 'show'
)

try 0.5 0.99

stop

finds=: 4 : 0
  'q z'=.Phi^:_1 'u v'=.(P"1 x,:y)/:"1 x
  {:@>{:@KFR@KF@cuspm"1 z;q;s=.^8
)
   
NB. cba finds anz 

 figure anz;cba;wbc

stop

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
  
  