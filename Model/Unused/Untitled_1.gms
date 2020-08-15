


$GDXIN G:\My Drive\Side Projects\Cutting stock (APA)\Code\Model\data5.gdx
sets i,t;
$load i,t
parameters due,M,N,q,h,w,BinH,BinW,Cap;
$load due,M,N,q,h,w,BinH,BinW,Cap
$GDXIN

Alias (i,j,k,l);
sets al(j,i), del(l,j,i) , Aset(l,k,j,i);
al(j,i)$(ord(j) <= ord(i))= yes;
del(l,j,i)$(ord(l) < ord(i) and ord(j) < ord(i) and ord(j)>=ord(l) )=yes;
Aset(l,k,j,i)$(ord(l)<= ord(k)) = yes;
variables BQ(l),objV;

Binary variables
alpha(j,i)
beta(k,j)
gamma(l,k)
a(l,k,j,i)
b(i,l,t)
x(i,t)
;
integer variables
o(j,i)
delta(l,j,i)
y(i)
z(i,t)
tp(i)
tn(i)
;

y.lo(i) = 1 ;
y.up(i) = M(i)  ;
equations
obj
c1(i)
c2(j,i)
c2_2(j,i)
c2_3(j)
c3(i)
c4(j)
c5(k,j)
c6(k,j)
c7(k)
c8(k)
*c9(l,j,i)
*c10(l,j,i)
*c11(l,j,i)
*c12(l)
c9_2(l)
c13(l)
c14(l)
c15(l,k,j,i)
c16(i,l)
c17(i,l,t)
c18(t)
c19(i)
c20(i)
c21(i)
c22(i)
;

obj..objV =e= sum(l,BQ(l));
c1(i)..sum(j$(al(j,i)),o(j,i)) =e= y(i);
c2(j,i)..o(j,i) =l= M(i)*alpha(j,i)   ;

c2_2(j,i)$(w(i) ne w(j))..alpha(j,i)=e= 0;
c2_3(j)..sum(i, o(j,i)*h(i) )=l= BinH;

C3(j)$(ord(j)<=N-1)..sum(i$( ord(i) ne ord(j)),alpha(j,i)) =l= N*alpha(j,j) ;

c4(j)..sum(k, beta(k,j)) =e= alpha(j,j);

c5(k,j)$(ord(j)<ord(k) and ord(k)>=1 ).. sum(i, h(i)*o(j,i) ) =l= sum(i, h(i)*o(k,i)) +(BinH+1)*(1-beta(k,j)) ;
c6(k,j)$(ord(j)>=ord(k) and ord(k)<N-1 ).. sum(i, h(i)*o(j,i) ) =l= sum(i, h(i)*o(k,i)) +(BinH)*(1-beta(k,j)) ;

c7(k).. sum(j,w[j]*beta(k,j)) =l= BinW * beta(k,k)*0.5;
c8(k).. sum(l$al(l,k), gamma(l,k)) =e= beta(k,k);

*c9(del(l,j,i)).. delta(l,j,i) =l= o(j,i);
*c10(del(l,j,i)).. delta(l,j,i) =l= M(i)*gamma(l,j);
*c11(del(l,j,i)).. o(j,i)-M(i)*(1-gamma(l,j)) =l= delta(l,j,i)   ;

*c12(l)$(ord(l) < N-1)..sum( i$(ord(i)>=ord(l)) , h(i)*gamma(l,l)*o(i,i) ) +
*         sum((i,j)$(ord(i)>=ord(l)+1 and del(l,j,i)) , h(i)*delta(l,j,i) ) =l= BinH*gamma(l,l);

c9_2(l)$(ord(l)<=N-1) ..  sum(k$(ord(k)>=ord(l)), sum(i,h(i)*o(k,i))*gamma(l,k)  ) =l= BinH*gamma(l,l)   ;

c13(l)$(ord(l) < N-1)..sum( k$(ord(k)>=ord(l)+1), gamma(l,k)) =l= (N-ord(l))*gamma(l,l);

c14(l).. sum(t, x(l,t)) =e= gamma(l,l);
c15(Aset(l,k,j,i)).. a(l,k,j,i) =e= alpha(j,i)*beta(k,j)*gamma(l,k);

c16(i,l).. sum((k,j)$Aset(l,k,j,i), a(l,k,j,i)*q(i)) =l= y(i)*BQ(l);

c17(i,l,t).. b(i,l,t) =e= sum((k,j)$(Aset(l,k,j,i)) , a(l,k,j,i)) * x(l,t);

c18(t).. sum(l , x(l,t)*BQ(l)) =l= Cap ;

c19(i).. sum( (t,l) , ord(t)*b(i,l,t) )-due(i) =l= tp(i)  ;
c20(i).. due(i)- sum( (t,l) , ord(t)*b(i,l,t) ) =l= tn(i) ;
c21(i).. tp(i) =l= 1 ;
c22(i).. tn(i) =l= 2 ;

model spliting /all/;
solve spliting minimizing objV using MINLP;
display q,N;
Execute_Unload 'Spliting_model_optimal2',alpha,beta,gamma,o,y,BQ;
execute "gdx2sqlite -i Spliting_model_optimal2.gdx -o  Spliting_model_optimal2.db" ;
