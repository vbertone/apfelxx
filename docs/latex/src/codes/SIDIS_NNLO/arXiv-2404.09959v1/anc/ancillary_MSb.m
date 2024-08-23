(* ::Package:: *)

(* 
Paper Name: NNLO QCD corrections to polarized semi-inclusive DIS
Authors: Saurav Goyal, Roman Lee, Sven-Olaf Moch, Vaibhav Pathak, Narayan Rana and V. Ravindran,
Year: April 2024
*)

(* This File Contains Subscript[\[ScriptCapitalG], {1,ab}] Coefficient Function in MSbar's Scheme see eq.(5) *)

(*
Symbols/Functions used:
Q2= Q^2, xp= x', zp= z',  
Li[a,b]= PolyLog[a,b], 
S12[a]= Nielsen Polylogarithm: PolyLog[1,2,a],
ArcTan[a]= Tangent Inverse function,
InvTanInt[a] = Integrate[ ArcTan[u]/u,{u,0,a}],
Theta[a]= HeavisideTheta[a].

zeta2, zeta3, . . . are the usual Zeta function with, zeta2 = Pi^2/6 and so on. 
muF2 is the square of factorisation scale,
as is the Strong Coupling constant at muF2 scale,

Plus Distributions are defined as,
Dxp[a_]= Subscript[[Log^a(1-xp)/(1-xp)], +], Dzp[a_]= Subscript[[Log^a(1-zp)/(1-zp)], +],
Delx= DiracDelta[1-xp], Delz= DiracDelta[1-zp],

CA=Nc and CF=(Nc^2-1)/(2 Nc) in SU(Nc) gauge theory. Nc=3 for QCD. 
nF is the number of active light quark flavours which can be 3 or 5.

Subscript[\[ScriptCapitalG], {1,ab}]= \!\(
\*SubscriptBox[
SuperscriptBox[\(\[Sum]\), \(\[Infinity]\)], \(i = 0\)]\(
\*SuperscriptBox[\((
\*SubscriptBox[\(a\), \(s\)]\((muF2)\))\), \(i\)]\ \[ScriptCapitalG]^
\*SubscriptBox[\((i)\), \({1, ab}\)]\)\);
here ab means `a' initiated and `b' fragmenting parton. 

Subscript[\[ScriptCapitalG]^(2),{1,qq}]= Subscript[\[ScriptCapitalG]^(2), {1,OverBar[q] OverBar[q]}] = (eq^2)*Subscript[G^(2),{1,qq,NS}] + (\!\(
\*SubscriptBox[\(\[Sum]\), \(i\)]\ \(eq_{i}^2\)\))*Subscript[G^(2),{1,qq,PS}]
Subscript[\[ScriptCapitalG]^(2),{1,qg}]= Subscript[\[ScriptCapitalG]^(2),{1,OverBar[q] g}]  = (eq^2)*Subscript[G^(2),{1,qg}]
Subscript[\[ScriptCapitalG]^(2),{1,gq}]= Subscript[\[ScriptCapitalG]^(2),{1,g OverBar[q]}]  = (eq^2)*Subscript[G^(2),{1,gq}]
Subscript[\[ScriptCapitalG]^(2),{1,gg}] = (\!\(
\*SubscriptBox[\(\[Sum]\), \(i\)]\ \(eq_{i}^2\)\))*Subscript[G^(2),{1,gg}]
Subscript[\[ScriptCapitalG]^(2),{1,qq'}]= Subscript[\[ScriptCapitalG]^(2),{1,OverBar[q] OverBar[q]'}]  = (eq^2)*Subscript[G^(2),{1,qq',[1]}]  +(eq'^2)*Subscript[G^(2),{1,qq',[2]}] +(eq*eq')*Subscript[G^(2),{1,qq',[3]}]
Subscript[\[ScriptCapitalG]^(2),{1,q OverBar[q]'}]= Subscript[\[ScriptCapitalG]^(2),{1,OverBar[q] q'}] = (eq^2)*Subscript[G^(2), {1,qq',[1]}] +(eq'^2)*Subscript[G^(2),{1,qq',[2]}] -(eq*eq')*Subscript[G^(2),{1,qq',[3]}]

above eq is the electric charge of quark q, eq' is the electric charge of quark q'.
(\!\(
\*SubscriptBox[\(\[Sum]\), \(i\)]\ \(eq_{i}^2\)\)) means sum over square of electric charge of active quark flavours. 

Notations: Subscript[G, {1,ab,MS}] here MS stands for MSbar's coefficient function.
*)

(*
For Larin to MSbar scheme transformation we used (see arxiv: 1409.5131) :
Zqq is to be multiplied by quark(antiquark) initiated processes,
Zgg = Zqg = Zgq=0

Zqq = 1 + as*( Zq^(1)(x)) + as^2*( Zq^(2)(x))+ ...

where,
Zq^(1)(x) = ZqNS^(1)(x) = -8*Cf*(1-x)
Zq^(2)(x) = ZqNS^(2)(x) + ZqPS^(2)(x)

ZqNS^(2)(x) = ZqNS^(2,+)(x) - ZqNS^(2,-)(x)
ZqNS^(2,+)(x) = ( Cf*Ca*( -592/9 -80/3*Log[(x)] -4*Log[(x)]^2 +592/9*x + 8/3*x*Log[(x)] +4*x*Log[(x)]^2 +8*zeta2 -8*zeta2*x)
              + Cf*Tf*nf*( +80/9 +16/3*Log[(x)] -80/9*x -16/3*x*Log[(x)])
              + Cf^2*( -16 -16*Log[(x)] +16*Log[(x)]*Log[(1-x)] +16*x -8*x*Log[(x)] -16*x*Log[(x)]*Log[(1-x)]) )
ZqNSq^(2,-)(x) = ((CF^2-1/2*CA*CF)*( 8*(1+x)*( 4*Li[2,(-x)] +4*Log[(x)]*Log[(1+x)] +2*zeta2 -Log[(x)]^2 -3*Log[(x)]) -56*(1-x) ))

ZqPS^(2)(x) = 2*Cf*Tf*( +8*(1 -x) +4*(3 -x)*Log[(x)] +2*(2 +x)*Log[(x)]^2)
*)



{"\!\(\*SubscriptBox[SuperscriptBox[\(G\), \((2)\)], \({1,qg,MS}\)]\)" -> 
  CF^2*(-86 + 16/(-1 + xp) + 80/zp - (4*(5 + xp)*zp)/(-1 + xp) - 
     (16*xp^2)/(-1 + xp + zp) + xp*(30 - 40/zp + 4*zp) + 
     (8*zeta2*(-8 + 14*zp - 8*zp^2 + xp^2*(2 - 4*zp + zp^2) + 
        xp*(4 - 8*zp + 6*zp^2)))/((-1 + xp)*zp) + 12*(-2 + 2/zp + zp)*
      Dxp[2] + 4*(-2 + xp*(-2 + zp^(-1)) + zp^(-1) + 2*zp)*Li[2, 1 - xp] - 
     (8*(-6 + 10*zp - 5*zp^2 + xp*(4 - 6*zp + 3*zp^2))*Li[2, 1 - zp])/
      ((-1 + xp)*zp) + (8*(xp^2 + (-1 + zp)^2)*
       Li[2, (2 - 2*xp)^(-1) + zp/(2 - 2*xp) - 
         Sqrt[1 - 2*zp + 4*xp*zp + zp^2]/(2 - 2*xp)])/((-1 + xp)*zp) + 
     (8*(xp^2 + (-1 + zp)^2)*Li[2, (2 - 2*xp)^(-1) - (2*xp)/(2 - 2*xp) - 
         zp/(2 - 2*xp) + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]/(2 - 2*xp)])/
      ((-1 + xp)*zp) + (4*(xp^2 + (-1 + zp)^2)*Log[2]^2)/((-1 + xp)*zp) - 
     (4*(-4 + 6*zp - 5*zp^2 + 4*xp*zp^2 + xp^2*(2 - 6*zp + zp^2))*
       Log[1 - xp]^2)/((-1 + xp)*zp) - 
     (4*(10 - 12*zp - 3*zp^2 + xp^2*(12 - 20*zp + zp^2) + 
        6*xp*(1 - 2*zp + 4*zp^2))*Log[xp]^2)/((-1 + xp)*zp) - 
     (4*(2 - 6*zp - zp^2 + 8*xp*zp^2 + xp^2*(6 - 10*zp + zp^2))*
       Log[1 - zp]^2)/((-1 + xp)*zp) - 
     (8*(-6 + 9*xp^2*(-1 + zp) + 16*zp - zp^2 + xp*(12 - 23*zp + 2*zp^2))*
       Log[zp])/((-1 + xp)*zp) + 
     (4*(-3 + 8*zp + 2*zp^2 - 2*xp*zp*(2 + 5*zp) + xp^2*(-3 + 10*zp + zp^2))*
       Log[zp]^2)/((-1 + xp)*zp) + Dxp[1]*(32 - 24/zp - 6*zp + 
       16*(-2 + 2/zp + zp)*Log[1 - zp] + 4*(-2 + 4/zp + zp)*Log[zp]) + 
     Log[1 - zp]*(2*(16 - 16/(-1 + xp) + xp*(-44 + 35/zp - 8*zp) - 17/zp + 
         14*zp - (8*xp^3)/(-1 + xp + zp)^2 + (8*xp^2)/(-1 + xp + zp)) + 
       (8*(zp*(2 + 3*zp) + xp^2*(-4 + 8*zp) + xp*(-2 + 2*zp - 9*zp^2))*
         Log[zp])/((-1 + xp)*zp)) + Log[Q2/muF2]^2*
      ((-4*(1 + xp)*(2 - 2*zp + zp^2))/zp + 8*(-2 + 2/zp + zp)*Dxp[0] + 
       Delx*(-8 + 12/zp + 5*zp + (-8 + 8/zp + 4*zp)*Log[1 - zp] + 
         (4 - 2*zp)*Log[zp])) + Dxp[0]*(34 + zeta2*(80 - 80/zp - 40*zp) - 
       36/zp - 20*zp + (-40 + 48/zp + 20*zp)*Li[2, 1 - zp] + 
       8*(-2 + 2/zp + zp)*Log[1 - zp]^2 - (12*Log[zp])/zp + 
       (8 - 4*zp)*Log[zp]^2 + Log[1 - zp]*(32 - 24/zp - 6*zp + 
         4*(-6 + 8/zp + 3*zp)*Log[zp])) + Log[Q2/muF2]*
      (2*(-4 + (1 + 5*xp)/zp + (2 - 8*xp)*zp) + 24*(-2 + 2/zp + zp)*Dxp[1] - 
       (8*(3 - 4*zp + 3*zp^2 + xp*(3 - 4*zp + zp^2))*Log[1 - xp])/zp + 
       (4*(3 - 2*zp + 2*xp*zp^2 + xp^2*(5 - 6*zp + 2*zp^2))*Log[xp])/
        ((-1 + xp)*zp) - (8*(1 + xp)*(2 - 2*zp + zp^2)*Log[1 - zp])/zp - 
       (8*(1 + xp - 2*xp*zp + zp^2)*Log[zp])/zp + 
       Dxp[0]*(8 + 6*zp + 16*(-2 + 2/zp + zp)*Log[1 - zp] + 
         4*(-2 + 4/zp + zp)*Log[zp]) + 
       Delx*(34 + zeta2*(80 - 80/zp - 40*zp) - 36/zp - 8*zp + 
         (-40 + 48/zp + 20*zp)*Li[2, 1 - zp] + 8*(-2 + 2/zp + zp)*
          Log[1 - zp]^2 + (12*(-1 + zp)^2*Log[zp])/zp + 
         (8 - 4*zp)*Log[zp]^2 + Log[1 - zp]*(8 + 6*zp + 4*(-6 + 8/zp + 3*zp)*
            Log[zp]))) - (4*(xp^2 + (-1 + zp)^2)*
       Log[1 + zp - Sqrt[1 - 2*zp + 4*xp*zp + zp^2]]^2)/((-1 + xp)*zp) - 
     (4*(xp^2 + (-1 + zp)^2)*Log[1 - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]]^2)/
      ((-1 + xp)*zp) + Log[2]*
      ((-8*(-2*(-1 + zp)^3 + 8*xp^3*zp + xp^2*(2 - 8*zp + 6*zp^2) + 
          xp*(-1 + 11*zp - 11*zp^2 + zp^3)))/((-1 + xp)*zp*
         Sqrt[1 + (-2 + 4*xp)*zp + zp^2]) + 
       (16*(xp^2 + (-1 + zp)^2)*Log[1 - xp])/((-1 + xp)*zp) - 
       (8*(xp^2 + (-1 + zp)^2)*Log[xp])/((-1 + xp)*zp) - 
       (16*(xp^2 + (-1 + zp)^2)*Log[1 + zp - Sqrt[1 - 2*zp + 4*xp*zp + 
             zp^2]])/((-1 + xp)*zp) + (8*(xp^2 + (-1 + zp)^2)*
         Log[1 - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/((-1 + xp)*zp)) + 
     (8*(-2*(-1 + zp)^3 + 8*xp^3*zp + xp^2*(2 - 8*zp + 6*zp^2) + 
        xp*(-1 + 11*zp - 11*zp^2 + zp^3))*
       Log[1 - 2*xp - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/
      ((-1 + xp)*zp*Sqrt[1 + (-2 + 4*xp)*zp + zp^2]) - 
     (4*(xp^2 + (-1 + zp)^2)*Log[1 - 2*xp - zp + Sqrt[1 - 2*zp + 4*xp*zp + 
            zp^2]]^2)/((-1 + xp)*zp) + 
     Log[1 - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]]*
      ((8*(-2*(-1 + zp)^3 + 8*xp^3*zp + xp^2*(2 - 8*zp + 6*zp^2) + 
          xp*(-1 + 11*zp - 11*zp^2 + zp^3)))/((-1 + xp)*zp*
         Sqrt[1 + (-2 + 4*xp)*zp + zp^2]) - 
       (8*(xp^2 + (-1 + zp)^2)*Log[1 - 2*xp - zp + 
           Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/((-1 + xp)*zp)) + 
     Log[1 - xp]*(2*(4 - 8/(-1 + xp) + xp*(-32 + 23/zp - 2*zp) - 5/zp + 
         8*zp) + (4*(11 - 10*zp - 4*zp^2 + xp^2*(13 - 22*zp + 2*zp^2) + 
          xp*(4 - 8*zp + 22*zp^2))*Log[xp])/((-1 + xp)*zp) - 
       (8*(1 - 4*zp - zp^2 + xp^2*(5 - 8*zp + zp^2) + 
          xp*(-2 + 4*zp + 4*zp^2))*Log[1 - zp])/((-1 + xp)*zp) + 
       (8*(-1 + 4*zp + zp^2 - xp*zp*(2 + 5*zp) + xp^2*(-3 + 6*zp))*Log[zp])/
        ((-1 + xp)*zp) - (8*(xp^2 + (-1 + zp)^2)*
         Log[1 + zp - Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/((-1 + xp)*zp) - 
       (8*(xp^2 + (-1 + zp)^2)*Log[1 - 2*xp - zp + 
           Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/((-1 + xp)*zp)) + 
     Log[xp]*(2*((8*xp^3)/(-1 + xp + zp)^2 + xp^2*(-8/(-1 + xp + zp) - 
           32/Sqrt[1 + (-2 + 4*xp)*zp + zp^2]) + 
         xp*(56 - (43 + 8/Sqrt[1 - 2*zp + 4*xp*zp + zp^2])/zp + 
           zp*(2 - 24/Sqrt[1 + (-2 + 4*xp)*zp + zp^2])) + 
         (-8 + 8*zp^3 - 37*Sqrt[1 - 2*zp + 4*xp*zp + zp^2] + 
           zp^2*(-24 + 13*Sqrt[1 - 2*zp + 4*xp*zp + zp^2]) + 
           4*zp*(6 + 13*Sqrt[1 - 2*zp + 4*xp*zp + zp^2]) - 
           xp*(4 + 4*zp^3 - 37*Sqrt[1 - 2*zp + 4*xp*zp + zp^2] + 
             zp^2*(-20 + 6*Sqrt[1 - 2*zp + 4*xp*zp + zp^2]) + 
             4*zp*(11 + 7*Sqrt[1 - 2*zp + 4*xp*zp + zp^2])))/
          ((-1 + xp)*zp*Sqrt[1 + (-2 + 4*xp)*zp + zp^2])) + 
       (4*(13 - 18*zp - 2*zp^2 + xp^2*(15 - 26*zp + 2*zp^2) + 
          xp*(4 - 8*zp + 26*zp^2))*Log[1 - zp])/((-1 + xp)*zp) - 
       (4*(-14 + 26*zp + 3*zp^2 + 2*xp^2*(-7 + 12*zp) + 
          xp*(-4 + 4*zp - 30*zp^2))*Log[zp])/((-1 + xp)*zp) - 
       (8*(xp^2 + (-1 + zp)^2)*Log[1 + zp - Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/
        ((-1 + xp)*zp) + (8*(xp^2 + (-1 + zp)^2)*
         Log[1 - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/((-1 + xp)*zp) + 
       (8*(xp^2 + (-1 + zp)^2)*Log[1 - 2*xp - zp + 
           Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/((-1 + xp)*zp)) + 
     Log[1 + zp - Sqrt[1 - 2*zp + 4*xp*zp + zp^2]]*
      ((-8*(-2*(-1 + zp)^3 + 8*xp^3*zp + xp^2*(2 - 8*zp + 6*zp^2) + 
          xp*(-1 + 11*zp - 11*zp^2 + zp^3)))/((-1 + xp)*zp*
         Sqrt[1 + (-2 + 4*xp)*zp + zp^2]) + 
       (8*(xp^2 + (-1 + zp)^2)*Log[1 - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/
        ((-1 + xp)*zp) + (16*(xp^2 + (-1 + zp)^2)*
         Log[1 - 2*xp - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/
        ((-1 + xp)*zp)) + Delx*(-24 + zeta2*(-32 + 24/zp - 2*zp) - 13*zp + 
       44*zeta3*(-2 + 2/zp + zp) + (-12/zp + 8*zp)*Li[2, 1 - zp] + 
       (16 - 24/zp - 8*zp)*Li[3, 1 - zp] + 
       (10*(2 - 2*zp + zp^2)*Log[1 - zp]^3)/(3*zp) + 
       (67 + zeta2*(8 - 16/zp - 4*zp) - 36/zp - 45*zp - 
         8*(-2 + zp)*Li[2, 1 - zp])*Log[zp] + (-8 + (13*zp)/2)*Log[zp]^2 - 
       (5*(-2 + zp)*Log[zp]^3)/3 + Log[1 - zp]^2*(16 - 12/zp - zp + 
         (-20 + 24/zp + 10*zp)*Log[zp]) + Log[1 - zp]*
        (zeta2*(56 - 56/zp - 28*zp) + 4*(9 - 9/zp - 5*zp) + 
         8*(-4 + 5/zp + 2*zp)*Li[2, 1 - zp] + (-12/zp + 8*zp)*Log[zp] + 
         (8*Log[zp]^2)/zp) + (64 - 48/zp - 32*zp)*S12[1 - zp]) + 
     Theta[1 - xp - zp]*((-16*zeta2*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2))/
        ((-1 + xp)*zp) - (8*(-4 + 4*xp^2*(-1 + zp) + 5*zp + 3*zp^2 + 
          xp*(8 - 7*zp - 6*zp^2)))/((-1 + xp)*zp) + 
       ((16*xp*(-1 + zp)^2 - 8*(-2 + zp)*zp)*Li[2, -(zp/(1 - xp - zp))])/
        ((-1 + xp)*zp) + ((-16*xp*(-1 + zp)^2 + 8*(-2 + zp)*zp)*
         Li[2, -((xp*zp)/(1 - xp - zp))])/((-1 + xp)*zp) - 
       (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
         Li[2, (1 - xp)^(-1) - xp/(1 - xp) - zp/(1 - xp)])/((-1 + xp)*zp) - 
       (4*(-4 + 2*zp + 3*zp^2 + xp^2*(-4 + 8*zp) - 2*xp*(3 - 6*zp + 7*zp^2))*
         Log[xp]^2)/((-1 + xp)*zp) + (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 
          2*xp*zp^2)*Log[1 - zp]^2)/((-1 + xp)*zp) + 
       8*(4 - 2/(-1 + xp) + xp*(-4 + 3/zp) - 3/zp)*Log[1 - xp - zp] + 
       (16*(-4 + 2/(-1 + xp) + xp*(4 - 3/zp) + 3/zp) - 
         (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[1 - xp - zp])/
          ((-1 + xp)*zp))*Log[zp] + (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 
          2*xp*zp^2)*Log[zp]^2)/((-1 + xp)*zp) + 
       Log[1 - xp]*(8*(-4 + 2/(-1 + xp) + xp*(4 - 3/zp) + 3/zp) + 
         (16*(-1 + zp^2 + xp^2*(-1 + 2*zp) + xp*(-2 + 4*zp - 4*zp^2))*
           Log[xp])/((-1 + xp)*zp) + (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 
            2*xp*zp^2)*Log[1 - zp])/((-1 + xp)*zp) + 
         (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[zp])/
          ((-1 + xp)*zp)) + Log[1 - zp]*
        (16*(-4 + 2/(-1 + xp) + xp*(4 - 3/zp) + 3/zp) - 
         (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[1 - xp - zp])/
          ((-1 + xp)*zp) + (32*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
           Log[zp])/((-1 + xp)*zp)) + Log[xp]*
        ((-8*(-6 + 6*xp^2*(-1 + zp) + 8*zp + zp^2 - 2*xp*(-6 + 5*zp)))/
          ((-1 + xp)*zp) + (16*(-2 + 2*zp + zp^2 + xp^2*(-2 + 4*zp) + 
            xp*(-2 + 4*zp - 6*zp^2))*Log[1 - zp])/((-1 + xp)*zp) - 
         (8*(-2 + 2*zp + zp^2 + xp^2*(-2 + 4*zp) + xp*(-2 + 4*zp - 6*zp^2))*
           Log[1 - xp - zp])/((-1 + xp)*zp) + 
         (8*(-4 + 6*zp + zp^2 + xp^2*(-4 + 8*zp) - 2*xp*(1 - 2*zp + 5*zp^2))*
           Log[zp])/((-1 + xp)*zp)) + 
       ((16*zeta2*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2))/((-1 + xp)*zp) + 
         (8*(-4 + 4*xp^2*(-1 + zp) + 5*zp + 3*zp^2 + xp*(8 - 7*zp - 6*zp^2)))/
          ((-1 + xp)*zp) + (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
           Li[2, (1 - xp)^(-1) - xp/(1 - xp) - zp/((1 - xp)*xp) + 
             zp^2/((1 - xp)*xp)])/((-1 + xp)*zp) - 
         (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[1 - zp]^2)/
          ((-1 + xp)*zp) + 8*(-4 + 2/(-1 + xp) + xp*(4 - 3/zp) + 3/zp)*
          Log[1 - xp - zp] + 8*(-4 + 2/(-1 + xp) + xp*(4 - 3/zp) + 3/zp)*
          Log[xp - zp] + (16*(4 - 2/(-1 + xp) + xp*(-4 + 3/zp) - 3/zp) + 
           (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[1 - xp - zp])/
            ((-1 + xp)*zp) + (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
             Log[xp - zp])/((-1 + xp)*zp))*Log[zp] - 
         (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[zp]^2)/
          ((-1 + xp)*zp) + Log[1 - zp]*(16*(4 - 2/(-1 + xp) + 
             xp*(-4 + 3/zp) - 3/zp) + (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 
              2*xp*zp^2)*Log[1 - xp - zp])/((-1 + xp)*zp) + 
           (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[xp - zp])/
            ((-1 + xp)*zp) - (32*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
             Log[zp])/((-1 + xp)*zp)) + Log[1 - xp]*
          (8*(4 - 2/(-1 + xp) + xp*(-4 + 3/zp) - 3/zp) + 
           (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[xp])/
            ((-1 + xp)*zp) - (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
             Log[1 - zp])/((-1 + xp)*zp) - (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 
              2*xp*zp^2)*Log[zp])/((-1 + xp)*zp)) + 
         Log[xp]*(8*(-4 + 2/(-1 + xp) + xp*(4 - 3/zp) + 3/zp) + 
           (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[1 - zp])/
            ((-1 + xp)*zp) - (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
             Log[1 - xp - zp])/((-1 + xp)*zp) - 
           (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[xp - zp])/
            ((-1 + xp)*zp) + (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
             Log[zp])/((-1 + xp)*zp)))*Theta[xp - zp] + 
       ((8*zeta2*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2))/((-1 + xp)*zp) + 
         (8*(-4 + 4*xp^2*(-1 + zp) + 5*zp + 3*zp^2 + xp*(8 - 7*zp - 6*zp^2)))/
          ((-1 + xp)*zp) - (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
           Li[2, -(xp/((1 - xp - zp)*(-xp + zp))) + xp^2/((1 - xp - zp)*(
                -xp + zp))])/((-1 + xp)*zp) + 
         ((4 + xp^2*(4 - 8*zp) - 8*zp + 8*xp*zp^2)*Log[1 - xp]^2)/
          (zp - xp*zp) + ((4 + xp^2*(4 - 8*zp) - 8*zp + 8*xp*zp^2)*Log[xp]^2)/
          (zp - xp*zp) - (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
           Log[1 - zp]^2)/((-1 + xp)*zp) + 
         ((4 + xp^2*(4 - 8*zp) - 8*zp + 8*xp*zp^2)*Log[1 - xp - zp]^2)/
          (zp - xp*zp) - (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
           Log[zp]^2)/((-1 + xp)*zp) + 8*(-4 + 2/(-1 + xp) + xp*(4 - 3/zp) + 
           3/zp)*Log[-xp + zp] + ((4 + xp^2*(4 - 8*zp) - 8*zp + 8*xp*zp^2)*
           Log[-xp + zp]^2)/(zp - xp*zp) + Log[1 - xp - zp]*
          (8*(-4 + 2/(-1 + xp) + xp*(4 - 3/zp) + 3/zp) - 
           (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[-xp + zp])/
            ((-1 + xp)*zp)) + Log[xp]*(8*(-4 + 2/(-1 + xp) + xp*(4 - 3/zp) + 
             3/zp) + (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
             Log[1 - zp])/((-1 + xp)*zp) - (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 
              2*xp*zp^2)*Log[1 - xp - zp])/((-1 + xp)*zp) + 
           (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[zp])/
            ((-1 + xp)*zp) - (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
             Log[-xp + zp])/((-1 + xp)*zp)) + Log[1 - xp]*
          (8*(4 - 2/(-1 + xp) + xp*(-4 + 3/zp) - 3/zp) + 
           (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[xp])/
            ((-1 + xp)*zp) - (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
             Log[1 - zp])/((-1 + xp)*zp) + (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 
              2*xp*zp^2)*Log[1 - xp - zp])/((-1 + xp)*zp) - 
           (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[zp])/
            ((-1 + xp)*zp) + (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
             Log[-xp + zp])/((-1 + xp)*zp)) + 
         Log[zp]*(16*(4 - 2/(-1 + xp) + xp*(-4 + 3/zp) - 3/zp) + 
           (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[1 - xp - zp])/
            ((-1 + xp)*zp) + (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
             Log[-xp + zp])/((-1 + xp)*zp)) + Log[1 - zp]*
          (16*(4 - 2/(-1 + xp) + xp*(-4 + 3/zp) - 3/zp) + 
           (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[1 - xp - zp])/
            ((-1 + xp)*zp) - (32*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
             Log[zp])/((-1 + xp)*zp) + (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 
              2*xp*zp^2)*Log[-xp + zp])/((-1 + xp)*zp)))*Theta[-xp + zp]) + 
     ((-8*zeta2*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2))/((-1 + xp)*zp) - 
       (8*(-4 + 4*xp^2*(-1 + zp) + 5*zp + 3*zp^2 + xp*(8 - 7*zp - 6*zp^2)))/
        ((-1 + xp)*zp) + (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
         Li[2, (1 - xp - zp)^(-1) - xp/(1 - xp - zp)])/((-1 + xp)*zp) + 
       ((16*xp*(-1 + zp)^2 - 8*(-2 + zp)*zp)*Li[2, xp^(-1) + zp^(-1) - 
           1/(xp*zp)])/((-1 + xp)*zp) + 
       ((-16*xp*(-1 + zp)^2 + 8*(-2 + zp)*zp)*Li[2, 1 - zp^(-1) + xp/zp])/
        ((-1 + xp)*zp) + ((4 + xp^2*(4 - 8*zp) - 8*zp + 8*xp*zp^2)*
         Log[1 - xp]^2)/((-1 + xp)*zp) - 
       (16*(-1 + zp^2 + xp^2*(-1 + 2*zp) + xp*(-2 + 4*zp - 4*zp^2))*
         Log[xp]^2)/((-1 + xp)*zp) + (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 
          2*xp*zp^2)*Log[1 - zp]^2)/((-1 + xp)*zp) + 
       (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[zp]^2)/
        ((-1 + xp)*zp) + 8*(4 - 2/(-1 + xp) + xp*(-4 + 3/zp) - 3/zp)*
        Log[-1 + xp + zp] + ((4 + xp^2*(4 - 8*zp) - 8*zp + 8*xp*zp^2)*
         Log[-1 + xp + zp]^2)/((-1 + xp)*zp) + 
       Log[zp]*(16*(-4 + 2/(-1 + xp) + xp*(4 - 3/zp) + 3/zp) - 
         (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[-1 + xp + zp])/
          ((-1 + xp)*zp)) + Log[1 - zp]*
        (16*(-4 + 2/(-1 + xp) + xp*(4 - 3/zp) + 3/zp) + 
         (32*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[zp])/
          ((-1 + xp)*zp) - (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
           Log[-1 + xp + zp])/((-1 + xp)*zp)) + Log[1 - xp]*
        (8*(-4 + 2/(-1 + xp) + xp*(4 - 3/zp) + 3/zp) + 
         (16*(-1 + zp^2 + xp^2*(-1 + 2*zp) + xp*(-2 + 4*zp - 4*zp^2))*
           Log[xp])/((-1 + xp)*zp) + (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 
            2*xp*zp^2)*Log[1 - zp])/((-1 + xp)*zp) + 
         (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[zp])/
          ((-1 + xp)*zp) - (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
           Log[-1 + xp + zp])/((-1 + xp)*zp)) + 
       Log[xp]*((-8*(-6 + 6*xp^2*(-1 + zp) + 8*zp + zp^2 - 2*xp*(-6 + 5*zp)))/
          ((-1 + xp)*zp) + (16*(-2 + 2*zp + zp^2 + xp^2*(-2 + 4*zp) + 
            xp*(-2 + 4*zp - 6*zp^2))*Log[1 - zp])/((-1 + xp)*zp) - 
         (32*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[zp])/
          ((-1 + xp)*zp) + (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
           Log[-1 + xp + zp])/((-1 + xp)*zp)))*Theta[-1 + xp + zp] + 
     Theta[xp - zp]*((-8*zeta2*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2))/
        ((-1 + xp)*zp) - (8*(-4 + 4*xp^2*(-1 + zp) + 5*zp + 3*zp^2 + 
          xp*(8 - 7*zp - 6*zp^2)))/((-1 + xp)*zp) + 
       (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Li[2, zp/xp])/
        ((-1 + xp)*zp) + (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
         Li[2, (-xp + zp)^(-1) - xp/(-xp + zp)])/((-1 + xp)*zp) + 
       ((4 + xp^2*(4 - 8*zp) - 8*zp + 8*xp*zp^2)*Log[1 - xp]^2)/
        ((-1 + xp)*zp) + (24*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
         Log[xp]^2)/((-1 + xp)*zp) + (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 
          2*xp*zp^2)*Log[1 - zp]^2)/((-1 + xp)*zp) + 
       8*(4 - 2/(-1 + xp) + xp*(-4 + 3/zp) - 3/zp)*Log[xp - zp] + 
       ((4 + xp^2*(4 - 8*zp) - 8*zp + 8*xp*zp^2)*Log[xp - zp]^2)/
        ((-1 + xp)*zp) + (16*(-4 + 2/(-1 + xp) + xp*(4 - 3/zp) + 3/zp) - 
         (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[xp - zp])/
          ((-1 + xp)*zp))*Log[zp] + (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 
          2*xp*zp^2)*Log[zp]^2)/((-1 + xp)*zp) + 
       Log[xp]*(16*(4 - 2/(-1 + xp) + xp*(-4 + 3/zp) - 3/zp) - 
         (32*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[1 - zp])/
          ((-1 + xp)*zp) + (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
           Log[xp - zp])/((-1 + xp)*zp) - (40*(1 + xp^2*(1 - 2*zp) - 2*zp + 
            2*xp*zp^2)*Log[zp])/((-1 + xp)*zp)) + 
       Log[1 - xp]*(8*(-4 + 2/(-1 + xp) + xp*(4 - 3/zp) + 3/zp) - 
         (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[xp])/
          ((-1 + xp)*zp) + (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
           Log[1 - zp])/((-1 + xp)*zp) - (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 
            2*xp*zp^2)*Log[xp - zp])/((-1 + xp)*zp) + 
         (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[zp])/
          ((-1 + xp)*zp)) + Log[1 - zp]*
        (16*(-4 + 2/(-1 + xp) + xp*(4 - 3/zp) + 3/zp) - 
         (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[xp - zp])/
          ((-1 + xp)*zp) + (32*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
           Log[zp])/((-1 + xp)*zp)) + 
       ((8*zeta2*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2))/((-1 + xp)*zp) + 
         (8*(-4 + 4*xp^2*(-1 + zp) + 5*zp + 3*zp^2 + xp*(8 - 7*zp - 6*zp^2)))/
          ((-1 + xp)*zp) - (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
           Li[2, -(xp/((1 - xp - zp)*(-xp + zp))) + xp^2/((1 - xp - zp)*(
                -xp + zp))])/((-1 + xp)*zp) + 
         ((4 + xp^2*(4 - 8*zp) - 8*zp + 8*xp*zp^2)*Log[1 - xp]^2)/
          (zp - xp*zp) + ((4 + xp^2*(4 - 8*zp) - 8*zp + 8*xp*zp^2)*Log[xp]^2)/
          (zp - xp*zp) - (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
           Log[1 - zp]^2)/((-1 + xp)*zp) + 
         ((4 + xp^2*(4 - 8*zp) - 8*zp + 8*xp*zp^2)*Log[xp - zp]^2)/
          (zp - xp*zp) - (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
           Log[zp]^2)/((-1 + xp)*zp) + 8*(-4 + 2/(-1 + xp) + xp*(4 - 3/zp) + 
           3/zp)*Log[-1 + xp + zp] + ((4 + xp^2*(4 - 8*zp) - 8*zp + 
            8*xp*zp^2)*Log[-1 + xp + zp]^2)/(zp - xp*zp) + 
         Log[xp - zp]*(8*(-4 + 2/(-1 + xp) + xp*(4 - 3/zp) + 3/zp) - 
           (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[-1 + xp + zp])/
            ((-1 + xp)*zp)) + Log[xp]*(8*(-4 + 2/(-1 + xp) + xp*(4 - 3/zp) + 
             3/zp) + (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
             Log[1 - zp])/((-1 + xp)*zp) - (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 
              2*xp*zp^2)*Log[xp - zp])/((-1 + xp)*zp) + 
           (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[zp])/
            ((-1 + xp)*zp) - (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
             Log[-1 + xp + zp])/((-1 + xp)*zp)) + Log[1 - xp]*
          (8*(4 - 2/(-1 + xp) + xp*(-4 + 3/zp) - 3/zp) + 
           (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[xp])/
            ((-1 + xp)*zp) - (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
             Log[1 - zp])/((-1 + xp)*zp) + (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 
              2*xp*zp^2)*Log[xp - zp])/((-1 + xp)*zp) - 
           (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[zp])/
            ((-1 + xp)*zp) + (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
             Log[-1 + xp + zp])/((-1 + xp)*zp)) + 
         Log[zp]*(16*(4 - 2/(-1 + xp) + xp*(-4 + 3/zp) - 3/zp) + 
           (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[xp - zp])/
            ((-1 + xp)*zp) + (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
             Log[-1 + xp + zp])/((-1 + xp)*zp)) + Log[1 - zp]*
          (16*(4 - 2/(-1 + xp) + xp*(-4 + 3/zp) - 3/zp) + 
           (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[xp - zp])/
            ((-1 + xp)*zp) - (32*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
             Log[zp])/((-1 + xp)*zp) + (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 
              2*xp*zp^2)*Log[-1 + xp + zp])/((-1 + xp)*zp)))*
        Theta[-1 + xp + zp]) + Theta[-xp + zp]*
      ((-8*(-4 + 4*xp^2*(-1 + zp) + 5*zp + 3*zp^2 + xp*(8 - 7*zp - 6*zp^2)))/
        ((-1 + xp)*zp) - (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
         Li[2, xp/zp])/((-1 + xp)*zp) - 
       (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
         Li[2, -(xp/(1 - xp)) + zp/(1 - xp)])/((-1 + xp)*zp) + 
       (20*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[xp]^2)/
        ((-1 + xp)*zp) + (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
         Log[1 - zp]^2)/((-1 + xp)*zp) + 
       (12*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[zp]^2)/
        ((-1 + xp)*zp) + Log[1 - xp]*(8*(-4 + 2/(-1 + xp) + xp*(4 - 3/zp) + 
           3/zp) - (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[xp])/
          ((-1 + xp)*zp) + (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
           Log[1 - zp])/((-1 + xp)*zp) + (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 
            2*xp*zp^2)*Log[zp])/((-1 + xp)*zp)) + 
       8*(4 - 2/(-1 + xp) + xp*(-4 + 3/zp) - 3/zp)*Log[-xp + zp] + 
       Log[1 - zp]*(16*(-4 + 2/(-1 + xp) + xp*(4 - 3/zp) + 3/zp) + 
         (32*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[zp])/
          ((-1 + xp)*zp) - (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
           Log[-xp + zp])/((-1 + xp)*zp)) + 
       Log[zp]*(16*(-4 + 2/(-1 + xp) + xp*(4 - 3/zp) + 3/zp) - 
         (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[-xp + zp])/
          ((-1 + xp)*zp)) + Log[xp]*(16*(4 - 2/(-1 + xp) + xp*(-4 + 3/zp) - 
           3/zp) - (32*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[1 - zp])/
          ((-1 + xp)*zp) - (32*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
           Log[zp])/((-1 + xp)*zp) + (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 
            2*xp*zp^2)*Log[-xp + zp])/((-1 + xp)*zp)) + 
       ((16*zeta2*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2))/((-1 + xp)*zp) + 
         (8*(-4 + 4*xp^2*(-1 + zp) + 5*zp + 3*zp^2 + xp*(8 - 7*zp - 6*zp^2)))/
          ((-1 + xp)*zp) + (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
           Li[2, (1 - xp)^(-1) - xp/(1 - xp) - zp/((1 - xp)*xp) + 
             zp^2/((1 - xp)*xp)])/((-1 + xp)*zp) - 
         (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[1 - zp]^2)/
          ((-1 + xp)*zp) - (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
           Log[zp]^2)/((-1 + xp)*zp) + Log[1 - xp]*
          (8*(4 - 2/(-1 + xp) + xp*(-4 + 3/zp) - 3/zp) + 
           (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[xp])/
            ((-1 + xp)*zp) - (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
             Log[1 - zp])/((-1 + xp)*zp) - (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 
              2*xp*zp^2)*Log[zp])/((-1 + xp)*zp)) + 
         8*(-4 + 2/(-1 + xp) + xp*(4 - 3/zp) + 3/zp)*Log[-xp + zp] + 
         8*(-4 + 2/(-1 + xp) + xp*(4 - 3/zp) + 3/zp)*Log[-1 + xp + zp] + 
         Log[xp]*(8*(-4 + 2/(-1 + xp) + xp*(4 - 3/zp) + 3/zp) + 
           (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[1 - zp])/
            ((-1 + xp)*zp) + (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
             Log[zp])/((-1 + xp)*zp) - (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 
              2*xp*zp^2)*Log[-xp + zp])/((-1 + xp)*zp) - 
           (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[-1 + xp + zp])/
            ((-1 + xp)*zp)) + Log[zp]*(16*(4 - 2/(-1 + xp) + xp*(-4 + 3/zp) - 
             3/zp) + (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
             Log[-xp + zp])/((-1 + xp)*zp) + 
           (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[-1 + xp + zp])/
            ((-1 + xp)*zp)) + Log[1 - zp]*
          (16*(4 - 2/(-1 + xp) + xp*(-4 + 3/zp) - 3/zp) - 
           (32*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[zp])/
            ((-1 + xp)*zp) + (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
             Log[-xp + zp])/((-1 + xp)*zp) + 
           (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[-1 + xp + zp])/
            ((-1 + xp)*zp)))*Theta[-1 + xp + zp])) + 
   CA*CF*(16*zeta2*(4 + (-1 + xp)^(-1) + xp*(2 - zp^(-1)) - 2/zp - 3*zp) + 
     (2*(84 - 36/(-1 + xp) - 37/zp + (-57 + 54/(-1 + xp))*zp + 82*zp^2 + 
        (36*xp^2)/(-1 + xp + zp) + xp*(-204 + 119/zp + 87*zp - 38*zp^2)))/9 + 
     (16*xp^2*Li[2, 1 - xp])/(-1 + xp) + 
     (4*(-8 + 4*xp^2 + 10*zp - 17*zp^2 + 4*xp*(1 - 2*zp + 4*zp^2))*
       Li[2, 1 - zp])/((-1 + xp)*zp) + (8*(2 + 2*zp + zp^2)*Li[2, -zp])/
      ((-1 + xp)*zp) + ((4 - 4*zp^2 + 8*xp*zp^2 + xp^2*(4 + 8*zp))*
       Li[2, (2 - 2*xp)^(-1) + zp/(2 - 2*xp) - 
         Sqrt[1 - 2*zp + 4*xp*zp + zp^2]/(2 - 2*xp)])/((-1 + xp)*zp) + 
     ((4 - 4*zp^2 + 8*xp*zp^2 + xp^2*(4 + 8*zp))*
       Li[2, (2 - 2*xp)^(-1) - (2*xp)/(2 - 2*xp) - zp/(2 - 2*xp) + 
         Sqrt[1 - 2*zp + 4*xp*zp + zp^2]/(2 - 2*xp)])/(zp - xp*zp) + 
     (8*(1 - zp^2 + 2*xp*zp^2 + xp^2*(1 + 2*zp))*
       Li[2, 1 - 1/(2*xp) + zp/(2*xp) - Sqrt[1 - 2*zp + 4*xp*zp + zp^2]/
          (2*xp)])/((-1 + xp)*zp) + 
     ((2 - 2*zp^2 + 4*xp*zp^2 + xp^2*(2 + 4*zp))*Log[2]^2)/((-1 + xp)*zp) + 
     ((2 + xp^2*(2 - 4*zp) - 4*zp + 4*xp*zp^2)*Log[1 - xp]^2)/
      ((-1 + xp)*zp) - (2*(-6 + 10*zp + zp^2 + xp^2*(-6 + 4*zp) - 
        6*xp*(1 - 2*zp + 3*zp^2))*Log[xp]^2)/((-1 + xp)*zp) + 
     (4*(3 + xp^2*(1 - 2*zp) - 6*zp + 2*zp^2 + 2*xp*zp^2)*Log[1 - zp]^2)/
      ((-1 + xp)*zp) + ((6 + 32*zp + 38*zp^2 - 52*xp*zp^2 + 
        2*xp^2*(-9 + 2*zp))*Log[zp]^2)/(zp - xp*zp) + 
     Dxp[1]*((4*(24 - 31/zp + 3*zp + 4*zp^2))/3 + 8*(-2 + 2/zp + zp)*
        Log[1 - zp] - (16*(1 + zp + zp^2)*Log[zp])/zp) + 
     Delx*Log[Q2/muF2]^2*(16 - 62/(3*zp) + 2*zp + (8*zp^2)/3 + 
       (-8 + 8/zp + 4*zp)*Log[1 - zp] - (8*(1 + zp + zp^2)*Log[zp])/zp) + 
     Log[1 - zp]*((4*(-30 + 12/(-1 + xp) + 35/zp - 3*zp - 8*zp^2 + 
          (6*xp^3)/(-1 + xp + zp)^2 - (6*xp^2)/(-1 + xp + zp) + 
          xp*(6 + 5/zp - 9*zp + 4*zp^2)))/3 + 
       (8*(2*xp^2 - 4*zp*(1 + zp) + xp*(1 - 2*zp + 7*zp^2))*Log[zp])/
        ((-1 + xp)*zp)) + Log[zp]*
      ((4*(-47 + 24*zp + 3*zp^2 + 8*zp^3 + xp^2*(11 + 27*zp - 9*zp^2 + 
            4*zp^3) - 3*xp*(-9 + 17*zp - zp^2 + 4*zp^3)))/(3*(-1 + xp)*zp) + 
       (8*(2 + 2*zp + zp^2)*Log[1 + zp])/((-1 + xp)*zp)) + 
     Dxp[0]*(8*zeta2*(-4 + 6/zp + 3*zp) + (4*(33 + 13/zp - 24*zp - 13*zp^2))/
        9 + (16 - 48/zp - 32*zp)*Li[2, 1 - zp] + 8*(2 + 2/zp + zp)*
        Li[2, -zp] + (-8 + 8/zp + 4*zp)*Log[1 - zp]^2 + 
       4*(-4 - 6/zp - 5*zp)*Log[zp]^2 + Log[1 - zp]*
        ((4*(24 - 31/zp + 3*zp + 4*zp^2))/3 - (16*(1 + zp + zp^2)*Log[zp])/
          zp) + Log[zp]*((8*(3 - 20/zp + 3*zp + 2*zp^2))/3 + 
         8*(2 + 2/zp + zp)*Log[1 + zp])) + Log[Q2/muF2]*
      ((4*(-1 + zp)*(-17 - 11*zp - 8*zp^2 + xp*(-23 - 5*zp + 4*zp^2)))/
        (3*zp) - (8*(1 + xp - 2*zp - 2*xp*zp + 2*zp^2)*Log[1 - zp])/zp + 
       (8*(1 + xp + 2*zp + 2*xp*zp + 4*zp^2)*Log[zp])/zp + 
       Dxp[0]*((4*(24 - 31/zp + 3*zp + 4*zp^2))/3 + 8*(-2 + 2/zp + zp)*
          Log[1 - zp] - (16*(1 + zp + zp^2)*Log[zp])/zp) + 
       Delx*(8*zeta2*(-4 + 6/zp + 3*zp) + (4*(33 + 13/zp - 24*zp - 13*zp^2))/
          9 + (16 - 48/zp - 32*zp)*Li[2, 1 - zp] + 8*(2 + 2/zp + zp)*
          Li[2, -zp] + (-8 + 8/zp + 4*zp)*Log[1 - zp]^2 + 
         4*(-4 - 6/zp - 5*zp)*Log[zp]^2 + Log[1 - zp]*
          ((4*(24 - 31/zp + 3*zp + 4*zp^2))/3 - (16*(1 + zp + zp^2)*Log[zp])/
            zp) + Log[zp]*((8*(3 - 20/zp + 3*zp + 2*zp^2))/3 + 
           8*(2 + 2/zp + zp)*Log[1 + zp]))) + 
     ((2 - 2*zp^2 + 4*xp*zp^2 + xp^2*(2 + 4*zp))*
       Log[1 + zp - Sqrt[1 - 2*zp + 4*xp*zp + zp^2]]^2)/((-1 + xp)*zp) + 
     ((2 - 2*zp^2 + 4*xp*zp^2 + xp^2*(2 + 4*zp))*
       Log[1 - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]]^2)/((-1 + xp)*zp) + 
     Log[1 + zp - Sqrt[1 - 2*zp + 4*xp*zp + zp^2]]*
      ((4*(-2*(-1 + zp)^3 + 8*xp^3*zp + 2*xp^2*(1 - 4*zp + 7*zp^2) + 
          xp*(-1 + 13*zp - 15*zp^2 + 3*zp^3)))/((-1 + xp)*zp*
         Sqrt[1 + (-2 + 4*xp)*zp + zp^2]) + 
       ((4 - 4*zp^2 + 8*xp*zp^2 + xp^2*(4 + 8*zp))*
         Log[1 - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/(zp - xp*zp)) + 
     Log[2]*((4*(-2*(-1 + zp)^3 + 8*xp^3*zp + 2*xp^2*(1 - 4*zp + 7*zp^2) + 
          xp*(-1 + 13*zp - 15*zp^2 + 3*zp^3)))/((-1 + xp)*zp*
         Sqrt[1 + (-2 + 4*xp)*zp + zp^2]) + 
       ((4 - 4*zp^2 + 8*xp*zp^2 + xp^2*(4 + 8*zp))*Log[xp])/((-1 + xp)*zp) + 
       ((4 - 4*zp^2 + 8*xp*zp^2 + xp^2*(4 + 8*zp))*
         Log[1 - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/(zp - xp*zp)) - 
     (4*(-2*(-1 + zp)^3 + 8*xp^3*zp + 2*xp^2*(1 - 4*zp + 7*zp^2) + 
        xp*(-1 + 13*zp - 15*zp^2 + 3*zp^3))*
       Log[1 - 2*xp - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/
      ((-1 + xp)*zp*Sqrt[1 + (-2 + 4*xp)*zp + zp^2]) + 
     ((2 - 2*zp^2 + 4*xp*zp^2 + xp^2*(2 + 4*zp))*
       Log[1 - 2*xp - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]]^2)/(zp - xp*zp) + 
     Log[1 - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]]*
      ((-4*(-2*(-1 + zp)^3 + 8*xp^3*zp + 2*xp^2*(1 - 4*zp + 7*zp^2) + 
          xp*(-1 + 13*zp - 15*zp^2 + 3*zp^3)))/((-1 + xp)*zp*
         Sqrt[1 + (-2 + 4*xp)*zp + zp^2]) + 
       ((4 - 4*zp^2 + 8*xp*zp^2 + xp^2*(4 + 8*zp))*
         Log[1 - 2*xp - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/
        ((-1 + xp)*zp)) + Log[1 - xp]*
      ((4*(-21 + 6/(-1 + xp) + 26/zp - 3*zp - 8*zp^2 + 
          xp*(-3 + 14/zp - 9*zp + 4*zp^2)))/3 - 
       (8*(1 + xp + xp^2 - 2*zp - 2*xp*zp - 2*xp^2*zp + 3*xp*zp^2)*Log[xp])/
        ((-1 + xp)*zp) + (4*(5 + xp^2*(1 - 2*zp) - 10*zp + 4*zp^2 + 
          xp*(-2 + 4*zp))*Log[1 - zp])/((-1 + xp)*zp) + 
       (16*(xp^2 + 3*xp*zp^2 - 2*zp*(1 + zp))*Log[zp])/((-1 + xp)*zp) + 
       ((4 - 4*zp^2 + 8*xp*zp^2 + xp^2*(4 + 8*zp))*
         Log[1 + zp - Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/(zp - xp*zp) + 
       ((4 - 4*zp^2 + 8*xp*zp^2 + xp^2*(4 + 8*zp))*
         Log[1 - 2*xp - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/
        ((-1 + xp)*zp)) + Log[xp]*
      ((4*((-6*xp^3)/(-1 + xp + zp)^2 + 6*xp^2*((-1 + xp + zp)^(-1) + 
            4/Sqrt[1 + (-2 + 4*xp)*zp + zp^2]) + 
          (4 - 6/Sqrt[1 + (-2 + 4*xp)*zp + zp^2] + 
            zp^2*(-3 - 18/Sqrt[1 + (-2 + 4*xp)*zp + zp^2]) + 
            18*zp*(1 + 1/Sqrt[1 + (-2 + 4*xp)*zp + zp^2]) + 
            zp^3*(8 + 6/Sqrt[1 + (-2 + 4*xp)*zp + zp^2]))/(zp - xp*zp) + 
          (xp*(xp*(22 + 3*zp + 8*zp^3 - 6/Sqrt[1 - 2*zp + 4*xp*zp + zp^2] + 
               zp^2*(-18 - 42/Sqrt[1 - 2*zp + 4*xp*zp + zp^2])) + 
             3*(9 + 1/Sqrt[1 - 2*zp + 4*xp*zp + zp^2] + zp*(-11 - 
                 13/Sqrt[1 - 2*zp + 4*xp*zp + zp^2]) + zp^3*(-8 - 
                 3/Sqrt[1 - 2*zp + 4*xp*zp + zp^2]) + zp^2*
                (1 + 15/Sqrt[1 - 2*zp + 4*xp*zp + zp^2]))))/(zp - xp*zp)))/
        3 - (8*(1 + xp + xp^2 - 4*zp - 2*xp*zp - 2*xp^2*zp + zp^2 + 
          3*xp*zp^2)*Log[1 - zp])/((-1 + xp)*zp) - 
       (4*(7 - 8*zp - 4*zp^2 + xp^2*(7 + 2*zp) + xp*(2 - 4*zp + 24*zp^2))*
         Log[zp])/((-1 + xp)*zp) + ((4 - 4*zp^2 + 8*xp*zp^2 + 
          xp^2*(4 + 8*zp))*Log[1 + zp - Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/
        ((-1 + xp)*zp) + ((4 - 4*zp^2 + 8*xp*zp^2 + xp^2*(4 + 8*zp))*
         Log[1 - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/(zp - xp*zp) + 
       ((4 - 4*zp^2 + 8*xp*zp^2 + xp^2*(4 + 8*zp))*
         Log[1 - 2*xp - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/
        (zp - xp*zp)) + Delx*(zeta3*(36 - 20/zp - 10*zp) + 
       (4*zeta2*(-24 + 31/zp + 9*zp - 4*zp^2))/3 - 
       (2*(-1169 + 723*zp + 204*zp^2 + 269*zp^3))/(27*zp) + 
       (8*(3 - 20/zp - 3*zp + 2*zp^2)*Li[2, 1 - zp])/3 + 4*zp*Li[2, -zp] + 
       (24 + 40/zp + 28*zp)*Li[3, 1 - zp] - 
       (8*(2 + 2*zp + zp^2)*Li[3, 1/2 - zp/2])/zp - 
       (8*(2 + 2*zp + zp^2)*Li[3, 1/2 + zp/2])/zp + (8 + 8/zp + 4*zp)*
        Li[3, -zp] + 8*(2 + 2/zp + zp)*Li[3, (-2*zp)/(1 - zp)] + 
       8*(2 + 2/zp + zp)*Li[3, (2*zp)/(1 + zp)] + 
       (8*(2 + 2*zp + zp^2)*Log[2]^3)/(3*zp) - 
       (2*(2 + 6*zp + zp^2)*Log[1 - zp]^3)/(3*zp) - 
       (10*(4 + 2*zp + 3*zp^2)*Log[zp]^3)/(3*zp) + 
       Log[1 - zp]^2*(16 - 62/(3*zp) + (8*zp^2)/3 + (4 - 4/zp - 6*zp)*
          Log[zp]) + Log[1 - zp]*(4*zeta2*(-6 + 2/zp + zp) + 
         (2*(57 + 26/zp - 48*zp - 26*zp^2))/9 + (8 - 40/zp - 28*zp)*
          Li[2, 1 - zp] + 8*(2 + 2/zp + zp)*Li[2, -zp] + 
         (4*(6 - 40/zp + 3*zp + 4*zp^2)*Log[zp])/3 + 4*(-4 - 6/zp - 5*zp)*
          Log[zp]^2) + (12*zeta2*(2 + 2/zp + zp) - 
         (8*(2 + 2*zp + zp^2)*Li[2, -zp])/zp)*Log[1 + zp] - 
       (4*(2 + 2*zp + zp^2)*Log[1 + zp]^3)/(3*zp) + 
       Log[2]^2*((-8 - 8/zp - 4*zp)*Log[1 - zp] + (-8 - 8/zp - 4*zp)*
          Log[1 + zp]) + Log[zp]*(90 + 304/(9*zp) + 38*zp + (124*zp^2)/9 + 
         zeta2*(48/zp + 32*zp) - (16*(4 + zp + 3*zp^2)*Li[2, 1 - zp])/zp + 
         (8 + 8/zp + 4*zp)*Li[2, -zp] + 4*zp*Log[1 + zp]) + 
       Log[zp]^2*(-16 - 80/(3*zp) - 5*zp + 6*(2 + 2/zp + zp)*Log[1 + zp]) + 
       Log[2]*((-8*zeta2*(2 + 2*zp + zp^2))/zp + (8 + 8/zp + 4*zp)*
          Log[1 - zp]^2 + (8 + 8/zp + 4*zp)*Log[1 + zp]^2) + 
       (-24 - 56/zp - 44*zp)*S12[1 - zp] - (8*(2 + 2*zp + zp^2)*S12[-zp])/
        zp) + Theta[1 - xp - zp]*
      ((8*zeta2*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2))/((-1 + xp)*zp) + 
       (4*(-4 + 4*xp^2*(-1 + zp) + 5*zp + 3*zp^2 + xp*(8 - 7*zp - 6*zp^2)))/
        ((-1 + xp)*zp) + ((-8*xp*(-1 + zp)^2 + 4*(-2 + zp)*zp)*
         Li[2, -(zp/(1 - xp - zp))])/((-1 + xp)*zp) + 
       ((8*xp*(-1 + zp)^2 - 4*(-2 + zp)*zp)*Li[2, -((xp*zp)/(1 - xp - zp))])/
        ((-1 + xp)*zp) + ((4 + xp^2*(4 - 8*zp) - 8*zp + 8*xp*zp^2)*
         Li[2, (1 - xp)^(-1) - xp/(1 - xp) - zp/(1 - xp)])/((-1 + xp)*zp) + 
       ((-8 + 4*zp + 6*zp^2 + 8*xp^2*(-1 + 2*zp) - 4*xp*(3 - 6*zp + 7*zp^2))*
         Log[xp]^2)/((-1 + xp)*zp) - (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 
          2*xp*zp^2)*Log[1 - zp]^2)/((-1 + xp)*zp) + 
       4*(-4 + 2/(-1 + xp) + xp*(4 - 3/zp) + 3/zp)*Log[1 - xp - zp] + 
       (8*(4 - 2/(-1 + xp) + xp*(-4 + 3/zp) - 3/zp) + 
         (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[1 - xp - zp])/
          ((-1 + xp)*zp))*Log[zp] - (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 
          2*xp*zp^2)*Log[zp]^2)/((-1 + xp)*zp) + 
       Log[1 - zp]*(8*(4 - 2/(-1 + xp) + xp*(-4 + 3/zp) - 3/zp) + 
         (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[1 - xp - zp])/
          ((-1 + xp)*zp) - (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
           Log[zp])/((-1 + xp)*zp)) + Log[1 - xp]*
        (4*(4 - 2/(-1 + xp) + xp*(-4 + 3/zp) - 3/zp) - 
         (8*(-1 + zp^2 + xp^2*(-1 + 2*zp) + xp*(-2 + 4*zp - 4*zp^2))*Log[xp])/
          ((-1 + xp)*zp) - (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
           Log[1 - zp])/((-1 + xp)*zp) - (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 
            2*xp*zp^2)*Log[zp])/((-1 + xp)*zp)) + 
       Log[xp]*((4*(-6 + 6*xp^2*(-1 + zp) + 8*zp + zp^2 - 2*xp*(-6 + 5*zp)))/
          ((-1 + xp)*zp) - (8*(-2 + 2*zp + zp^2 + xp^2*(-2 + 4*zp) + 
            xp*(-2 + 4*zp - 6*zp^2))*Log[1 - zp])/((-1 + xp)*zp) + 
         (4*(-2 + 2*zp + zp^2 + xp^2*(-2 + 4*zp) + xp*(-2 + 4*zp - 6*zp^2))*
           Log[1 - xp - zp])/((-1 + xp)*zp) - 
         (4*(-4 + 6*zp + zp^2 + xp^2*(-4 + 8*zp) - 2*xp*(1 - 2*zp + 5*zp^2))*
           Log[zp])/((-1 + xp)*zp)) + ((-8*zeta2*(-2 + zp))/(-1 + xp) - 
         (4*(13 + 14*xp^2 + 3*zp - xp*(25 + 6*zp)))/(-1 + xp) + 
         ((4 + xp^2*(4 - 8*zp) - 4*zp^2 + 8*xp*zp^2)*
           Li[2, (1 - zp)^(-1) - zp/(xp*(1 - zp))])/((-1 + xp)*zp) + 
         ((4 + xp^2*(4 - 8*zp) - 8*zp + 8*xp*zp^2)*
           Li[2, (1 - xp)^(-1) - xp/(1 - xp) - zp/((1 - xp)*xp) + 
             zp^2/((1 - xp)*xp)])/(zp - xp*zp) + 
         (8*(-1 + zp^2 - 2*xp*zp^2 + xp^2*(-1 + 2*zp))*Log[1 - xp]^2)/
          ((-1 + xp)*zp) + (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
           Log[1 - zp]^2)/((-1 + xp)*zp) + 4*(4 - 2/(-1 + xp) + 
           xp*(-4 + 3/zp) - 3/zp)*Log[1 - xp - zp] - 
         (4*(7 - 9*xp + 4*xp^2 - zp)*Log[xp - zp])/(-1 + xp) + 
         ((8*(7 - 9*xp + 4*xp^2 - zp))/(-1 + xp) - 
           (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[1 - xp - zp])/
            ((-1 + xp)*zp) - (8*(-2 + zp)*Log[xp - zp])/(-1 + xp))*Log[zp] + 
         (8*(-2 + zp)*Log[zp]^2)/(-1 + xp) + 
         Log[xp]*((-4*(7 - 9*xp + 4*xp^2 - zp))/(-1 + xp) - 
           (8*(-2 + zp)*Log[1 - zp])/(-1 + xp) + 
           (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[1 - xp - zp])/
            ((-1 + xp)*zp) + (8*(-2 + zp)*Log[xp - zp])/(-1 + xp) - 
           (8*(-2 + zp)*Log[zp])/(-1 + xp)) + Log[1 - xp]*
          ((4*(3 + 8*zp - 2*zp^2 + xp^2*(3 + 4*zp) - 2*xp*(3 + 5*zp)))/
            ((-1 + xp)*zp) - (8*(-2 + zp)*Log[xp])/(-1 + xp) + 
           (8*(-2 + zp)*Log[1 - zp])/(-1 + xp) + 
           (8*(1 + xp^2*(1 - 2*zp) - zp^2 + 2*xp*zp^2)*Log[xp - zp])/
            ((-1 + xp)*zp) - (8*(1 + xp^2*(1 - 2*zp) + 2*zp - 2*zp^2 + 
              2*xp*zp^2)*Log[zp])/((-1 + xp)*zp)) + Log[1 - zp]*
          ((4*(-3 + xp*(6 - 17*zp) + 13*zp - zp^2 + xp^2*(-3 + 8*zp)))/
            ((-1 + xp)*zp) - (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
             Log[1 - xp - zp])/((-1 + xp)*zp) - 
           (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[xp - zp])/
            ((-1 + xp)*zp) + (8*(1 + xp^2*(1 - 2*zp) - 4*zp + zp^2 + 
              2*xp*zp^2)*Log[zp])/((-1 + xp)*zp)))*Theta[xp - zp] + 
       ((zeta2*(4 + xp^2*(4 - 8*zp) - 8*zp + 8*xp*zp^2))/(zp - xp*zp) - 
         (4*(-4 + 4*xp^2*(-1 + zp) + 5*zp + 3*zp^2 + xp*(8 - 7*zp - 6*zp^2)))/
          ((-1 + xp)*zp) + ((4 + xp^2*(4 - 8*zp) - 8*zp + 8*xp*zp^2)*
           Li[2, -(xp/((1 - xp - zp)*(-xp + zp))) + xp^2/((1 - xp - zp)*(
                -xp + zp))])/((-1 + xp)*zp) + 
         ((2 + xp^2*(2 - 4*zp) - 4*zp + 4*xp*zp^2)*Log[1 - xp]^2)/
          ((-1 + xp)*zp) + ((2 + xp^2*(2 - 4*zp) - 4*zp + 4*xp*zp^2)*
           Log[xp]^2)/((-1 + xp)*zp) + (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 
            2*xp*zp^2)*Log[1 - zp]^2)/((-1 + xp)*zp) + 
         ((2 + xp^2*(2 - 4*zp) - 4*zp + 4*xp*zp^2)*Log[1 - xp - zp]^2)/
          ((-1 + xp)*zp) + (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
           Log[zp]^2)/((-1 + xp)*zp) + 4*(4 - 2/(-1 + xp) + xp*(-4 + 3/zp) - 
           3/zp)*Log[-xp + zp] + ((2 + xp^2*(2 - 4*zp) - 4*zp + 4*xp*zp^2)*
           Log[-xp + zp]^2)/((-1 + xp)*zp) + 
         Log[zp]*(8*(-4 + 2/(-1 + xp) + xp*(4 - 3/zp) + 3/zp) - 
           (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[1 - xp - zp])/
            ((-1 + xp)*zp) - (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
             Log[-xp + zp])/((-1 + xp)*zp)) + Log[1 - zp]*
          (8*(-4 + 2/(-1 + xp) + xp*(4 - 3/zp) + 3/zp) - 
           (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[1 - xp - zp])/
            ((-1 + xp)*zp) + (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
             Log[zp])/((-1 + xp)*zp) - (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 
              2*xp*zp^2)*Log[-xp + zp])/((-1 + xp)*zp)) + 
         Log[1 - xp - zp]*(4*(4 - 2/(-1 + xp) + xp*(-4 + 3/zp) - 3/zp) + 
           ((4 + xp^2*(4 - 8*zp) - 8*zp + 8*xp*zp^2)*Log[-xp + zp])/
            ((-1 + xp)*zp)) + Log[xp]*(4*(4 - 2/(-1 + xp) + xp*(-4 + 3/zp) - 
             3/zp) - (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[1 - zp])/
            ((-1 + xp)*zp) + ((4 + xp^2*(4 - 8*zp) - 8*zp + 8*xp*zp^2)*
             Log[1 - xp - zp])/((-1 + xp)*zp) - 
           (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[zp])/
            ((-1 + xp)*zp) + ((4 + xp^2*(4 - 8*zp) - 8*zp + 8*xp*zp^2)*
             Log[-xp + zp])/((-1 + xp)*zp)) + Log[1 - xp]*
          (4*(-4 + 2/(-1 + xp) + xp*(4 - 3/zp) + 3/zp) + 
           ((4 + xp^2*(4 - 8*zp) - 8*zp + 8*xp*zp^2)*Log[xp])/(zp - xp*zp) + 
           (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[1 - zp])/
            ((-1 + xp)*zp) + ((4 + xp^2*(4 - 8*zp) - 8*zp + 8*xp*zp^2)*
             Log[1 - xp - zp])/(zp - xp*zp) + 
           (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[zp])/
            ((-1 + xp)*zp) + ((4 + xp^2*(4 - 8*zp) - 8*zp + 8*xp*zp^2)*
             Log[-xp + zp])/(zp - xp*zp)))*Theta[-xp + zp]) + 
     ((zeta2*(4 + xp^2*(4 - 8*zp) - 8*zp + 8*xp*zp^2))/((-1 + xp)*zp) + 
       (4*(-4 + 4*xp^2*(-1 + zp) + 5*zp + 3*zp^2 + xp*(8 - 7*zp - 6*zp^2)))/
        ((-1 + xp)*zp) + ((4 + xp^2*(4 - 8*zp) - 8*zp + 8*xp*zp^2)*
         Li[2, (1 - xp - zp)^(-1) - xp/(1 - xp - zp)])/(zp - xp*zp) + 
       ((-8*xp*(-1 + zp)^2 + 4*(-2 + zp)*zp)*Li[2, xp^(-1) + zp^(-1) - 
           1/(xp*zp)])/((-1 + xp)*zp) + ((8*xp*(-1 + zp)^2 - 4*(-2 + zp)*zp)*
         Li[2, 1 - zp^(-1) + xp/zp])/((-1 + xp)*zp) + 
       ((2 + xp^2*(2 - 4*zp) - 4*zp + 4*xp*zp^2)*Log[1 - xp]^2)/
        (zp - xp*zp) + (8*(-1 + zp^2 + xp^2*(-1 + 2*zp) + 
          xp*(-2 + 4*zp - 4*zp^2))*Log[xp]^2)/((-1 + xp)*zp) - 
       (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[1 - zp]^2)/
        ((-1 + xp)*zp) - (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
         Log[zp]^2)/((-1 + xp)*zp) + 4*(-4 + 2/(-1 + xp) + xp*(4 - 3/zp) + 
         3/zp)*Log[-1 + xp + zp] + ((2 + xp^2*(2 - 4*zp) - 4*zp + 4*xp*zp^2)*
         Log[-1 + xp + zp]^2)/(zp - xp*zp) + 
       Log[xp]*((4*(-6 + 6*xp^2*(-1 + zp) + 8*zp + zp^2 - 2*xp*(-6 + 5*zp)))/
          ((-1 + xp)*zp) - (8*(-2 + 2*zp + zp^2 + xp^2*(-2 + 4*zp) + 
            xp*(-2 + 4*zp - 6*zp^2))*Log[1 - zp])/((-1 + xp)*zp) + 
         (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[zp])/
          ((-1 + xp)*zp) - (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
           Log[-1 + xp + zp])/((-1 + xp)*zp)) + 
       Log[zp]*(8*(4 - 2/(-1 + xp) + xp*(-4 + 3/zp) - 3/zp) + 
         (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[-1 + xp + zp])/
          ((-1 + xp)*zp)) + Log[1 - zp]*
        (8*(4 - 2/(-1 + xp) + xp*(-4 + 3/zp) - 3/zp) - 
         (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[zp])/
          ((-1 + xp)*zp) + (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
           Log[-1 + xp + zp])/((-1 + xp)*zp)) + Log[1 - xp]*
        (4*(4 - 2/(-1 + xp) + xp*(-4 + 3/zp) - 3/zp) - 
         (8*(-1 + zp^2 + xp^2*(-1 + 2*zp) + xp*(-2 + 4*zp - 4*zp^2))*Log[xp])/
          ((-1 + xp)*zp) - (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
           Log[1 - zp])/((-1 + xp)*zp) - (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 
            2*xp*zp^2)*Log[zp])/((-1 + xp)*zp) + 
         ((4 + xp^2*(4 - 8*zp) - 8*zp + 8*xp*zp^2)*Log[-1 + xp + zp])/
          ((-1 + xp)*zp)))*Theta[-1 + xp + zp] + 
     Theta[-xp + zp]*((8*zeta2*(1 + xp^2*(1 - 2*zp) - zp^2 + 2*xp*zp^2))/
        ((-1 + xp)*zp) + (4*(-4 + 4*xp^2*(-1 + zp) + 5*zp + 3*zp^2 + 
          xp*(8 - 7*zp - 6*zp^2)))/((-1 + xp)*zp) + 
       (4*(-2 + zp)*Li[2, xp/zp])/(-1 + xp) + 
       ((4 + xp^2*(4 - 8*zp) - 8*zp + 8*xp*zp^2)*
         Li[2, -(xp/(1 - xp)) + zp/(1 - xp)])/((-1 + xp)*zp) + 
       ((4 + xp^2*(4 - 8*zp) - 4*zp^2 + 8*xp*zp^2)*
         Li[2, -(-xp + zp)^(-1) + zp/(-xp + zp)])/((-1 + xp)*zp) + 
       (4*(-1 + zp^2 - 2*xp*zp^2 + xp^2*(-1 + 2*zp))*
         Li[2, -(xp/(-xp + zp)) + (xp*zp)/(-xp + zp)])/((-1 + xp)*zp) + 
       ((2 + xp^2*(2 - 4*zp) - 20*zp + 8*zp^2 + 4*xp*zp^2)*Log[xp]^2)/
        (zp - xp*zp) - (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
         Log[1 - zp]^2)/((-1 + xp)*zp) + 
       (2*(-4 + 6*zp + zp^2 - 8*xp*zp^2 + xp^2*(-4 + 8*zp))*Log[zp]^2)/
        ((-1 + xp)*zp) + Log[1 - xp]*(4*(4 - 2/(-1 + xp) + xp*(-4 + 3/zp) - 
           3/zp) + (8*(-2 + zp)*Log[xp])/(-1 + xp) - 
         (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[1 - zp])/
          ((-1 + xp)*zp) - (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
           Log[zp])/((-1 + xp)*zp)) + 4*(-4 + 2/(-1 + xp) + xp*(4 - 3/zp) + 
         3/zp)*Log[-xp + zp] + Log[1 - zp]*
        (8*(4 - 2/(-1 + xp) + xp*(-4 + 3/zp) - 3/zp) - 
         (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[zp])/
          ((-1 + xp)*zp) + (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
           Log[-xp + zp])/((-1 + xp)*zp)) + 
       Log[xp]*((4*(-3 + xp*(6 - 17*zp) + 13*zp - zp^2 + xp^2*(-3 + 8*zp)))/
          ((-1 + xp)*zp) + (4*(3 + xp^2*(3 - 6*zp) - 8*zp + zp^2 + 6*xp*zp^2)*
           Log[1 - zp])/((-1 + xp)*zp) + (8*(1 + xp^2*(1 - 2*zp) - 4*zp + 
            zp^2 + 2*xp*zp^2)*Log[zp])/((-1 + xp)*zp) + 
         ((4 + xp^2*(4 - 8*zp) - 8*zp + 8*xp*zp^2)*Log[-xp + zp])/
          (zp - xp*zp)) + Log[zp]*(8*(4 - 2/(-1 + xp) + xp*(-4 + 3/zp) - 
           3/zp) - (4*(-2 + 2*zp + zp^2 - 4*xp*zp^2 + xp^2*(-2 + 4*zp))*
           Log[-xp + zp])/((-1 + xp)*zp)) + 
       ((-8*zeta2*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2))/((-1 + xp)*zp) - 
         (4*(-4 + 4*xp^2*(-1 + zp) + 5*zp + 3*zp^2 + xp*(8 - 7*zp - 6*zp^2)))/
          ((-1 + xp)*zp) + ((4 + xp^2*(4 - 8*zp) - 8*zp + 8*xp*zp^2)*
           Li[2, (1 - xp)^(-1) - xp/(1 - xp) - zp/((1 - xp)*xp) + 
             zp^2/((1 - xp)*xp)])/(zp - xp*zp) + 
         (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[1 - zp]^2)/
          ((-1 + xp)*zp) + (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
           Log[zp]^2)/((-1 + xp)*zp) + Log[1 - xp]*
          (4*(-4 + 2/(-1 + xp) + xp*(4 - 3/zp) + 3/zp) - 
           (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[xp])/
            ((-1 + xp)*zp) + (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
             Log[1 - zp])/((-1 + xp)*zp) + (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 
              2*xp*zp^2)*Log[zp])/((-1 + xp)*zp)) + 
         4*(4 - 2/(-1 + xp) + xp*(-4 + 3/zp) - 3/zp)*Log[-xp + zp] + 
         4*(4 - 2/(-1 + xp) + xp*(-4 + 3/zp) - 3/zp)*Log[-1 + xp + zp] + 
         Log[zp]*(8*(-4 + 2/(-1 + xp) + xp*(4 - 3/zp) + 3/zp) - 
           (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[-xp + zp])/
            ((-1 + xp)*zp) - (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
             Log[-1 + xp + zp])/((-1 + xp)*zp)) + Log[1 - zp]*
          (8*(-4 + 2/(-1 + xp) + xp*(4 - 3/zp) + 3/zp) + 
           (16*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[zp])/
            ((-1 + xp)*zp) - (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
             Log[-xp + zp])/((-1 + xp)*zp) - (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 
              2*xp*zp^2)*Log[-1 + xp + zp])/((-1 + xp)*zp)) + 
         Log[xp]*(4*(4 - 2/(-1 + xp) + xp*(-4 + 3/zp) - 3/zp) - 
           (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[1 - zp])/
            ((-1 + xp)*zp) - (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
             Log[zp])/((-1 + xp)*zp) + (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 
              2*xp*zp^2)*Log[-xp + zp])/((-1 + xp)*zp) + 
           (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[-1 + xp + zp])/
            ((-1 + xp)*zp)))*Theta[-1 + xp + zp]) + 
     Theta[xp - zp]*((4*(13 + 14*xp^2 + 3*zp - xp*(25 + 6*zp)))/(-1 + xp) + 
       zeta2*(8 + xp*(8 - 4/zp) - 4/zp - 8*zp + 8/(zp - xp*zp)) - 
       (4*(-2 + zp)*Li[2, zp/xp])/(-1 + xp) + 
       (4*(-1 + zp^2 - 2*xp*zp^2 + xp^2*(-1 + 2*zp))*
         Li[2, xp/(1 - zp) - zp/(1 - zp)])/((-1 + xp)*zp) + 
       ((4 + xp^2*(4 - 8*zp) - 8*zp + 8*xp*zp^2)*
         Li[2, (-xp + zp)^(-1) - xp/(-xp + zp)])/(zp - xp*zp) + 
       ((6 + xp^2*(6 - 12*zp) + 4*zp - 8*zp^2 + 12*xp*zp^2)*Log[1 - xp]^2)/
        ((-1 + xp)*zp) - (12*(-2 + zp)*Log[xp]^2)/(-1 + xp) - 
       (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[1 - zp]^2)/
        ((-1 + xp)*zp) + (4*(7 - 9*xp + 4*xp^2 - zp)*Log[xp - zp])/
        (-1 + xp) + ((2 + xp^2*(2 - 4*zp) - 4*zp + 4*xp*zp^2)*Log[xp - zp]^2)/
        (zp - xp*zp) + ((-8*(7 - 9*xp + 4*xp^2 - zp))/(-1 + xp) + 
         (4*(-2 + zp)*Log[xp - zp])/(-1 + xp))*Log[zp] - 
       (8*(-2 + zp)*Log[zp]^2)/(-1 + xp) + 
       Log[xp]*((8*(7 - 9*xp + 4*xp^2 - zp))/(-1 + xp) + 
         (8*(1 + xp^2*(1 - 2*zp) - 4*zp + zp^2 + 2*xp*zp^2)*Log[1 - zp])/
          ((-1 + xp)*zp) - (4*(-2 + zp)*Log[xp - zp])/(-1 + xp) + 
         (20*(-2 + zp)*Log[zp])/(-1 + xp)) + Log[1 - xp]*
        ((-4*(3 + 8*zp - 2*zp^2 + xp^2*(3 + 4*zp) - 2*xp*(3 + 5*zp)))/
          ((-1 + xp)*zp) - (8*(1 + xp^2*(1 - 2*zp) + 2*zp - 2*zp^2 + 
            2*xp*zp^2)*Log[xp])/((-1 + xp)*zp) - (8*(-2 + zp)*Log[1 - zp])/
          (-1 + xp) + (8 + xp*(8 - 4/zp) - 4/zp - 8*zp + 8/(zp - xp*zp))*
          Log[xp - zp] + (8*(1 + xp^2*(1 - 2*zp) + 2*zp - 2*zp^2 + 2*xp*zp^2)*
           Log[zp])/((-1 + xp)*zp)) + Log[1 - zp]*
        ((-4*(-3 + xp*(6 - 17*zp) + 13*zp - zp^2 + xp^2*(-3 + 8*zp)))/
          ((-1 + xp)*zp) + (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
           Log[xp - zp])/((-1 + xp)*zp) - (8*(1 + xp^2*(1 - 2*zp) - 4*zp + 
            zp^2 + 2*xp*zp^2)*Log[zp])/((-1 + xp)*zp)) + 
       ((-4*(13 + 14*xp^2 + 3*zp - xp*(25 + 6*zp)))/(-1 + xp) + 
         4*zeta2*(-2 + xp*(-2 + zp^(-1)) + zp^(-1) + 2*zp - 2/(zp - xp*zp)) + 
         ((4 + xp^2*(4 - 8*zp) - 4*zp^2 + 8*xp*zp^2)*
           Li[2, (1 - zp)^(-1) - zp/(xp*(1 - zp))])/((-1 + xp)*zp) + 
         ((4 + xp^2*(4 - 8*zp) - 8*zp + 8*xp*zp^2)*
           Li[2, -(xp/((1 - xp - zp)*(-xp + zp))) + xp^2/((1 - xp - zp)*(
                -xp + zp))])/((-1 + xp)*zp) + 
         ((6 + xp^2*(6 - 12*zp) + 4*zp - 8*zp^2 + 12*xp*zp^2)*Log[1 - xp]^2)/
          (zp - xp*zp) + ((2 + xp^2*(2 - 4*zp) - 4*zp + 4*xp*zp^2)*Log[xp]^2)/
          ((-1 + xp)*zp) + (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
           Log[1 - zp]^2)/((-1 + xp)*zp) + 
         ((2 + xp^2*(2 - 4*zp) - 4*zp + 4*xp*zp^2)*Log[xp - zp]^2)/
          ((-1 + xp)*zp) + (8*(-2 + zp)*Log[zp]^2)/(-1 + xp) + 
         4*(4 - 2/(-1 + xp) + xp*(-4 + 3/zp) - 3/zp)*Log[-1 + xp + zp] + 
         ((2 + xp^2*(2 - 4*zp) - 4*zp + 4*xp*zp^2)*Log[-1 + xp + zp]^2)/
          ((-1 + xp)*zp) + Log[zp]*((8*(7 - 9*xp + 4*xp^2 - zp))/(-1 + xp) - 
           (8*(-2 + zp)*Log[xp - zp])/(-1 + xp) - 
           (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[-1 + xp + zp])/
            ((-1 + xp)*zp)) + Log[1 - zp]*
          ((4*(-3 + xp*(6 - 17*zp) + 13*zp - zp^2 + xp^2*(-3 + 8*zp)))/
            ((-1 + xp)*zp) - (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*
             Log[xp - zp])/((-1 + xp)*zp) + (8*(1 + xp^2*(1 - 2*zp) - 4*zp + 
              zp^2 + 2*xp*zp^2)*Log[zp])/((-1 + xp)*zp) - 
           (8*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*zp^2)*Log[-1 + xp + zp])/
            ((-1 + xp)*zp)) + Log[xp - zp]*((-4*(7 - 9*xp + 4*xp^2 - zp))/
            (-1 + xp) + ((4 + xp^2*(4 - 8*zp) - 8*zp + 8*xp*zp^2)*
             Log[-1 + xp + zp])/((-1 + xp)*zp)) + 
         Log[xp]*((-4*(7 - 9*xp + 4*xp^2 - zp))/(-1 + xp) - 
           (8*(-2 + zp)*Log[1 - zp])/(-1 + xp) + (8 + xp*(8 - 4/zp) - 4/zp - 
             8*zp + 8/(zp - xp*zp))*Log[xp - zp] - (8*(-2 + zp)*Log[zp])/
            (-1 + xp) + ((4 + xp^2*(4 - 8*zp) - 8*zp + 8*xp*zp^2)*
             Log[-1 + xp + zp])/((-1 + xp)*zp)) + Log[1 - xp]*
          ((4*(3 + 8*zp - 2*zp^2 + xp^2*(3 + 4*zp) - 2*xp*(3 + 5*zp)))/
            ((-1 + xp)*zp) + 4*(-2 + xp*(-2 + zp^(-1)) + zp^(-1) + 2*zp - 
             2/(zp - xp*zp))*Log[xp] + (8*(-2 + zp)*Log[1 - zp])/(-1 + xp) + 
           4*(-2 + xp*(-2 + zp^(-1)) + zp^(-1) + 2*zp - 2/(zp - xp*zp))*
            Log[xp - zp] - (8*(1 + xp^2*(1 - 2*zp) + 2*zp - 2*zp^2 + 
              2*xp*zp^2)*Log[zp])/((-1 + xp)*zp) + 
           ((4 + xp^2*(4 - 8*zp) - 8*zp + 8*xp*zp^2)*Log[-1 + xp + zp])/
            (zp - xp*zp)))*Theta[-1 + xp + zp])), 
 "\!\(\*SubscriptBox[SuperscriptBox[\(G\), \((2)\)], \({1,gq,MS}\)]\)" -> 
  (2*CA*(10 - 38*zp + 30*zp^2 + xp*(-11 + 29*zp - 20*zp^2)))/((-1 + zp)*zp) + 
   (2*CA*zeta2*(2 + xp^2*(6 - 10*zp) - 2*zp - 3*zp^2 + 
      2*xp*(2 - 4*zp + 5*zp^2)))/((-1 + zp)*zp) - 
   (4*CA*(3 + 2*(-3 - 3*xp + xp^2)*zp + (5 + 2*xp)*zp^2)*Li[2, 1 - xp])/
    ((-1 + zp)*zp) + (4*CA*(-2 - 2*xp^2 + 2*xp*(-2 + zp) + zp)*Li[2, -xp])/
    zp + (4*CA*(1 - 2*xp + 2*xp^2)*(-1 + 2*zp)*Li[2, 1 - zp])/
    ((-1 + zp)*zp) + (2*CA*(-1 + 2*zp - 2*zp^2 + xp^5*(-2 + 4*zp) - 
      4*xp^4*(2 - 7*zp + 7*zp^2) + xp*(-5 + 8*zp + 6*zp^2 - 4*zp^3) + 
      xp^3*(-13 + 46*zp - 60*zp^2 + 40*zp^3) + 
      xp^2*(-11 + 28*zp + 12*zp^2 - 80*zp^3 + 40*zp^4))*
     Li[2, 1/2 - xp/2 - Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]/2])/
    ((1 + xp^2 + xp*(2 - 4*zp))^(3/2)*(-1 + zp)*zp) + 
   (2*CA*(-1 + 2*zp - 2*zp^2 + xp^5*(-2 + 4*zp) - 
      4*xp^4*(2 - 7*zp + 7*zp^2) + xp*(-5 + 8*zp + 6*zp^2 - 4*zp^3) + 
      xp^3*(-13 + 46*zp - 60*zp^2 + 40*zp^3) + 
      xp^2*(-11 + 28*zp + 12*zp^2 - 80*zp^3 + 40*zp^4))*
     Li[2, 1/2 + xp/2 - Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]/2])/
    ((1 + xp^2 + xp*(2 - 4*zp))^(3/2)*(-1 + zp)*zp) - 
   (2*CA*(-1 + 2*zp - 2*zp^2 + xp^5*(-2 + 4*zp) - 
      4*xp^4*(2 - 7*zp + 7*zp^2) + xp*(-5 + 8*zp + 6*zp^2 - 4*zp^3) + 
      xp^3*(-13 + 46*zp - 60*zp^2 + 40*zp^3) + 
      xp^2*(-11 + 28*zp + 12*zp^2 - 80*zp^3 + 40*zp^4))*
     Li[2, 1/2 - 1/(2*xp) - Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]/(2*xp)])/
    ((1 + xp^2 + xp*(2 - 4*zp))^(3/2)*(-1 + zp)*zp) - 
   (2*CA*(-1 + 2*zp - 2*zp^2 + xp^5*(-2 + 4*zp) - 
      4*xp^4*(2 - 7*zp + 7*zp^2) + xp*(-5 + 8*zp + 6*zp^2 - 4*zp^3) + 
      xp^3*(-13 + 46*zp - 60*zp^2 + 40*zp^3) + 
      xp^2*(-11 + 28*zp + 12*zp^2 - 80*zp^3 + 40*zp^4))*
     Li[2, 1/2 + 1/(2*xp) - Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]/(2*xp)])/
    ((1 + xp^2 + xp*(2 - 4*zp))^(3/2)*(-1 + zp)*zp) - 
   (CA*(-1 + 2*xp)*(3 - 8*zp + 6*zp^2)*Log[1 - xp]^2)/((-1 + zp)*zp) + 
   (CA*(2 + xp^5*(4 - 8*zp) + 3*Sqrt[1 + 2*xp + xp^2 - 4*xp*zp] + 
      4*zp^2*(1 + 6*Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]) - 
      4*zp*(1 + 8*Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]) - 
      4*xp^4*(-14*zp^2 + 2*(-2 + Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]) + 
        zp*(14 + Sqrt[1 + 2*xp + xp^2 - 4*xp*zp])) + 
      2*xp^3*(13 - 40*zp^3 - 3*Sqrt[1 + 2*xp + xp^2 - 4*xp*zp] + 
        6*zp^2*(10 + Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]) + 
        2*zp*(-23 + 7*Sqrt[1 + 2*xp + xp^2 - 4*xp*zp])) + 
      xp^2*(22 - 24*zp^2 - 80*zp^4 + 15*Sqrt[1 + 2*xp + xp^2 - 4*xp*zp] + 
        16*zp^3*(10 + Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]) - 
        4*zp*(14 + 17*Sqrt[1 + 2*xp + xp^2 - 4*xp*zp])) - 
      2*xp*(-5 - 8*Sqrt[1 + 2*xp + xp^2 - 4*xp*zp] + 
        zp^2*(6 - 86*Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]) + 
        4*zp*(2 + 9*Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]) + 
        zp^3*(-4 + 48*Sqrt[1 + 2*xp + xp^2 - 4*xp*zp])))*Log[xp]^2)/
    ((1 + xp^2 + xp*(2 - 4*zp))^(3/2)*(-1 + zp)*zp) + 
   Delz*Log[Q2/muF2]^2*(-12*CA*(-1 + xp) + 2*CA*(-1 + 2*xp)*Log[1 - xp] + 
     4*CA*(1 + xp)*Log[xp]) + Dzp[1]*(-24*CA*(-1 + xp) + 
     4*CA*(-1 + 2*xp)*Log[1 - xp] + 8*CA*(1 + xp)*Log[xp]) + 
   Dzp[0]*(2*CA*(-8 + 9*xp) - 8*CA*xp*zeta2 + 8*CA*(1 + xp)*Li[2, 1 - xp] - 
     4*(CA + 2*CA*xp)*Li[2, -xp] + 2*CA*(-1 + 2*xp)*Log[1 - xp]^2 - 
     2*CA*(3 + 4*xp)*Log[xp]^2 + Log[1 - xp]*(-24*CA*(-1 + xp) + 
       12*CA*Log[xp]) + Log[xp]*(2*CA*(-11 + 16*xp) - 
       4*(CA + 2*CA*xp)*Log[1 + xp])) + 
   Log[Q2/muF2]*((24*CA*(-1 + xp)*(-1 + 2*zp))/zp - 
     (4*CA*(-1 + 2*xp)*(-1 + 2*zp)*Log[1 - xp])/zp - 
     (8*CA*(1 + xp)*(-1 + 2*zp)*Log[xp])/zp + 
     Dzp[0]*(-24*CA*(-1 + xp) + 4*CA*(-1 + 2*xp)*Log[1 - xp] + 
       8*CA*(1 + xp)*Log[xp]) + Delz*(2*CA*(-8 + 9*xp) - 8*CA*xp*zeta2 + 
       8*CA*(1 + xp)*Li[2, 1 - xp] - 4*(CA + 2*CA*xp)*Li[2, -xp] + 
       2*CA*(-1 + 2*xp)*Log[1 - xp]^2 - 2*CA*(3 + 4*xp)*Log[xp]^2 + 
       Log[1 - xp]*(-24*CA*(-1 + xp) + 12*CA*Log[xp]) + 
       Log[xp]*(2*CA*(-11 + 16*xp) - 4*(CA + 2*CA*xp)*Log[1 + xp]))) - 
   (4*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*Log[1 - zp]^2)/((-1 + zp)*zp) + 
   (2*CA*(2*xp^4*(8 - 12*zp + 9*zp^2) + zp*(16 - 31*zp + 25*zp^2) - 
      8*xp^3*(-2 + 12*zp - 16*zp^2 + 11*zp^3) - 
      2*xp*(8 - 23*zp + 63*zp^2 - 80*zp^3 + 52*zp^4) + 
      xp^2*(-16 + 86*zp - 73*zp^2 - 13*zp^3 + 76*zp^4))*Log[zp])/
    ((1 + xp^2 + xp*(2 - 4*zp))*(xp - zp)*(-1 + zp)*zp) - 
   (4*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*Log[zp]^2)/((-1 + zp)*zp) + 
   Log[1 - zp]*((2*CA*(16 - 59*zp + 74*zp^2 - 31*zp^3 + 
        4*xp^2*(4 - 9*zp + 6*zp^2) + xp*(-32 + 94*zp - 88*zp^2 + 22*zp^3)))/
      ((-1 + zp)*zp*(-1 + xp + zp)) + 
     (2*CA*(2 + 2*xp^2*(-1 + zp) - 6*zp + 7*zp^2 - 2*xp*(2 - 6*zp + 7*zp^2))*
       Log[zp])/((-1 + zp)*zp)) + Log[1 - xp]*
    ((-4*CA*(7 - 21*zp + xp^2*(-1 + zp)*zp + 15*zp^2 + 
        xp*(-7 + 20*zp - 14*zp^2)))/((-1 + zp)*zp) + 
     (4*CA*(-3 + xp^2 + 11*zp + 4*xp*(-1 + zp)*zp - 8*zp^2)*Log[xp])/
      ((-1 + zp)*zp) + (2*CA*(5 - 2*xp^2*(-1 + zp) - 10*zp + 7*zp^2 - 
        2*xp*(5 - 10*zp + 7*zp^2))*Log[1 - zp])/((-1 + zp)*zp) + 
     (2*CA*(4 - 6*zp + 2*xp^2*zp + 7*zp^2 - 2*xp*(4 - 6*zp + 7*zp^2))*
       Log[zp])/((-1 + zp)*zp)) + 
   (2*CA*(-1 + 2*zp - 2*zp^2 + xp^5*(-2 + 4*zp) - 
      4*xp^4*(2 - 7*zp + 7*zp^2) + xp*(-5 + 8*zp + 6*zp^2 - 4*zp^3) + 
      xp^3*(-13 + 46*zp - 60*zp^2 + 40*zp^3) + 
      xp^2*(-11 + 28*zp + 12*zp^2 - 80*zp^3 + 40*zp^4))*
     Log[1 + xp - Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]]*
     Log[1 - xp + Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]])/
    ((1 + xp^2 + xp*(2 - 4*zp))^(3/2)*(-1 + zp)*zp) - 
   (2*CA*(-1 + 2*zp - 2*zp^2 + xp^5*(-2 + 4*zp) - 
      4*xp^4*(2 - 7*zp + 7*zp^2) + xp*(-5 + 8*zp + 6*zp^2 - 4*zp^3) + 
      xp^3*(-13 + 46*zp - 60*zp^2 + 40*zp^3) + 
      xp^2*(-11 + 28*zp + 12*zp^2 - 80*zp^3 + 40*zp^4))*
     Log[1 + xp - Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]]*
     Log[-1 + xp + Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]])/
    ((1 + xp^2 + xp*(2 - 4*zp))^(3/2)*(-1 + zp)*zp) + 
   Log[2]*((-4*CA*(-1 + 2*zp - 2*zp^2 + xp^5*(-2 + 4*zp) - 
        4*xp^4*(2 - 7*zp + 7*zp^2) + xp*(-5 + 8*zp + 6*zp^2 - 4*zp^3) + 
        xp^3*(-13 + 46*zp - 60*zp^2 + 40*zp^3) + 
        xp^2*(-11 + 28*zp + 12*zp^2 - 80*zp^3 + 40*zp^4))*Log[xp])/
      ((1 + xp^2 + xp*(2 - 4*zp))^(3/2)*(-1 + zp)*zp) - 
     (2*CA*(-1 + 2*zp - 2*zp^2 + xp^5*(-2 + 4*zp) - 
        4*xp^4*(2 - 7*zp + 7*zp^2) + xp*(-5 + 8*zp + 6*zp^2 - 4*zp^3) + 
        xp^3*(-13 + 46*zp - 60*zp^2 + 40*zp^3) + 
        xp^2*(-11 + 28*zp + 12*zp^2 - 80*zp^3 + 40*zp^4))*
       Log[1 - xp + Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]])/
      ((1 + xp^2 + xp*(2 - 4*zp))^(3/2)*(-1 + zp)*zp) + 
     (2*CA*(-1 + 2*zp - 2*zp^2 + xp^5*(-2 + 4*zp) - 
        4*xp^4*(2 - 7*zp + 7*zp^2) + xp*(-5 + 8*zp + 6*zp^2 - 4*zp^3) + 
        xp^3*(-13 + 46*zp - 60*zp^2 + 40*zp^3) + 
        xp^2*(-11 + 28*zp + 12*zp^2 - 80*zp^3 + 40*zp^4))*
       Log[-1 + xp + Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]])/
      ((1 + xp^2 + xp*(2 - 4*zp))^(3/2)*(-1 + zp)*zp)) + 
   Log[xp]*(2*CA*(36 - 11/zp + (1 + 2*xp + xp^2 - 4*xp*zp)^(-1) - 
       xp^3/(1 + 2*xp + xp^2 - 4*xp*zp) + xp*(-29 + 4/(-1 + zp) + 12/zp + 
         (-xp + zp)^(-1) + (-1 + xp + zp)^(-1) - (1 + 2*xp + xp^2 - 4*xp*zp)^
          (-1)) + xp^2*(2 + 2/(xp - zp) - 6/(-1 + zp) + 6/zp - 
         2/(-1 + xp + zp) + (1 + 2*xp + xp^2 - 4*xp*zp)^(-1))) + 
     (4*CA*(-2 - 2*xp^2 + 2*xp*(-2 + zp) + zp)*Log[1 + xp])/zp + 
     (4*CA*(-3 + xp^2 + 10*zp - 8*zp^2 + 2*xp*zp*(-1 + 2*zp))*Log[1 - zp])/
      ((-1 + zp)*zp) + (2*CA*(-7 + 18*zp + 4*xp^2*zp - 14*zp^2 + 
        2*xp*(-1 - 6*zp + 8*zp^2))*Log[zp])/((-1 + zp)*zp) + 
     (2*CA*(-1 + 2*zp - 2*zp^2 + xp^5*(-2 + 4*zp) - 
        4*xp^4*(2 - 7*zp + 7*zp^2) + xp*(-5 + 8*zp + 6*zp^2 - 4*zp^3) + 
        xp^3*(-13 + 46*zp - 60*zp^2 + 40*zp^3) + 
        xp^2*(-11 + 28*zp + 12*zp^2 - 80*zp^3 + 40*zp^4))*
       Log[-1 + xp + Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]])/
      ((1 + xp^2 + xp*(2 - 4*zp))^(3/2)*(-1 + zp)*zp) + 
     (2*CA*(-1 + 2*zp - 2*zp^2 + xp^5*(-2 + 4*zp) - 
        4*xp^4*(2 - 7*zp + 7*zp^2) + xp*(-5 + 8*zp + 6*zp^2 - 4*zp^3) + 
        xp^3*(-13 + 46*zp - 60*zp^2 + 40*zp^3) + 
        xp^2*(-11 + 28*zp + 12*zp^2 - 80*zp^3 + 40*zp^4))*
       Log[1 + xp + Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]])/
      ((1 + xp^2 + xp*(2 - 4*zp))^(3/2)*(-1 + zp)*zp)) + 
   Delz*(CA*(57 - 59*xp) + 4*CA*(-5 + 6*xp)*zeta2 + 7*CA*(1 - 2*xp)*zeta3 + 
     2*CA*(1 + 4*xp)*Li[2, 1 - xp] + 4*CA*(1 + xp)*Li[2, -xp] - 
     20*CA*Li[3, 1 - xp] + 4*CA*(1 + 2*xp)*Li[3, 1/2 - xp/2] + 
     4*CA*(1 + 2*xp)*Li[3, 1/2 + xp/2] - 4*(CA + 2*CA*xp)*
      Li[3, (-2*xp)/(1 - xp)] - 4*(CA + 2*CA*xp)*Li[3, (2*xp)/(1 + xp)] - 
     (4*(CA + 2*CA*xp)*Log[2]^3)/3 + (CA*(1 + 6*xp)*Log[1 - xp]^3)/3 + 
     (CA*(7 + 10*xp)*Log[xp]^3)/3 + Log[1 - xp]^2*(-10*CA*(-1 + xp) - 
       4*CA*(-1 + xp)*Log[xp]) + Log[1 - xp]*(CA*(-12 + 13*xp) + 4*CA*zeta2 + 
       2*CA*(5 + 2*xp)*Li[2, 1 - xp] - 4*(CA + 2*CA*xp)*Li[2, -xp] + 
       2*CA*(-11 + 16*xp)*Log[xp] - 4*CA*(2 + xp)*Log[xp]^2) - 
     8*(CA + 2*CA*xp)*zeta2*Log[1 + xp] + (2*(CA + 2*CA*xp)*Log[1 + xp]^3)/
      3 + Log[2]^2*(2*CA*(1 + 2*xp)*Log[1 - xp] + 2*CA*(1 + 2*xp)*
        Log[1 + xp]) + Log[xp]^2*(CA*(29/2 - 22*xp) + 
       4*CA*(1 + 2*xp)*Log[1 + xp]) + Log[2]*(4*CA*(1 + 2*xp)*zeta2 - 
       2*(CA + 2*CA*xp)*Log[1 - xp]^2 - 2*(CA + 2*CA*xp)*Log[1 + xp]^2) + 
     Log[xp]*(40*CA - 8*CA*zeta2 - 4*(CA + 2*CA*xp)*Li[2, 1 - xp] + 
       4*CA*(1 + 2*xp)*Li[2, -xp] + 4*CA*(1 + xp)*Log[1 + xp] - 
       2*(CA + 2*CA*xp)*Log[1 + xp]^2) + 8*CA*xp*S12[1 - xp]) + 
   Theta[1 - xp - zp]*((-4*CA*(-1 + 2*xp)*zeta2*(1 - 2*zp + 2*zp^2))/
      ((-1 + zp)*zp) - (4*CA*(1 - 2*zp + 2*zp^2 + xp*(-1 - 4*zp + 4*zp^2)))/
      ((-1 + zp)*zp) + (2*CA*(2 - 2*xp^2*(-1 + zp) - 2*zp + zp^2 - 
        2*xp*(2 - 2*zp + zp^2))*Li[2, -(zp/(1 - xp - zp))])/((-1 + zp)*zp) + 
     (2*CA*(-2 + 2*xp^2*(-1 + zp) + 2*zp - zp^2 + 2*xp*(2 - 2*zp + zp^2))*
       Li[2, -((xp*zp)/(1 - xp - zp))])/((-1 + zp)*zp) - 
     (2*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*
       Li[2, (1 - xp)^(-1) - xp/(1 - xp) - zp/(1 - xp)])/((-1 + zp)*zp) + 
     (CA*(2 - 6*xp^2*(-1 + zp) + 2*zp - 5*zp^2 + 2*xp*(-2 - 2*zp + 5*zp^2))*
       Log[xp]^2)/((-1 + zp)*zp) + (4*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*
       Log[1 - zp]^2)/((-1 + zp)*zp) + (4*CA*(-1 + xp + 2*zp - 2*zp^2)*
       Log[1 - xp - zp])/((-1 + zp)*zp) + 
     ((-8*CA*(-1 + xp + 2*zp - 2*zp^2))/((-1 + zp)*zp) - 
       (4*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*Log[1 - xp - zp])/
        ((-1 + zp)*zp))*Log[zp] + (4*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*
       Log[zp]^2)/((-1 + zp)*zp) + Log[1 - xp]*
      ((-4*CA*(-1 + xp + 2*zp - 2*zp^2))/((-1 + zp)*zp) + 
       (4*CA*(1 + 2*xp^2 + zp - 2*xp*(1 + zp))*Log[xp])/zp + 
       (4*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*Log[1 - zp])/((-1 + zp)*zp) + 
       (4*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*Log[zp])/((-1 + zp)*zp)) + 
     Log[1 - zp]*((-8*CA*(-1 + xp + 2*zp - 2*zp^2))/((-1 + zp)*zp) - 
       (4*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*Log[1 - xp - zp])/
        ((-1 + zp)*zp) + (8*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*Log[zp])/
        ((-1 + zp)*zp)) + Log[xp]*
      ((-2*CA*(xp*(4 - 6*zp) + 6*xp^2*(-1 + zp) + zp*(-5 + 7*zp)))/
        ((-1 + zp)*zp) + (4*CA*(2*xp^2*(-1 + zp) + 2*xp*(2 - 3*zp)*zp + 
          zp*(-2 + 3*zp))*Log[1 - zp])/((-1 + zp)*zp) + 
       (2*CA*(-2*xp^2*(-1 + zp) + (2 - 3*zp)*zp + 2*xp*zp*(-2 + 3*zp))*
         Log[1 - xp - zp])/((-1 + zp)*zp) + 
       (2*CA*(2 + 2*xp^2*(-1 + zp) - 6*zp + 7*zp^2 - 
          2*xp*(2 - 6*zp + 7*zp^2))*Log[zp])/((-1 + zp)*zp)) + 
     ((4*CA*zeta2*(2 + 2*xp^2 + 2*xp*(-2 + zp) - zp))/(-1 + zp) + 
       (4*CA*(-7*xp^2 + 2*(-2 + zp) + 4*xp*(1 + zp)))/(-1 + zp) + 
       (2*CA*(1 + 2*xp^2*zp + zp^2 - 2*xp*(1 + zp^2))*
         Li[2, (1 - zp)^(-1) - zp/(xp*(1 - zp))])/((-1 + zp)*zp) + 
       (2*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*Li[2, (1 - xp)^(-1) - 
           xp/(1 - xp) - zp/((1 - xp)*xp) + zp^2/((1 - xp)*xp)])/
        ((-1 + zp)*zp) - (4*CA*(1 + 2*xp^2*zp + zp^2 - 2*xp*(1 + zp^2))*
         Log[1 - xp]^2)/((-1 + zp)*zp) - 
       (4*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*Log[1 - zp]^2)/((-1 + zp)*zp) - 
       (4*CA*(-1 + xp + 2*zp - 2*zp^2)*Log[1 - xp - zp])/((-1 + zp)*zp) - 
       (2*CA*(5 - 6*xp + 6*xp^2 - 3*zp)*Log[xp - zp])/(-1 + zp) + 
       ((4*CA*(5 - 6*xp + 6*xp^2 - 3*zp))/(-1 + zp) + 
         (4*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*Log[1 - xp - zp])/
          ((-1 + zp)*zp) + (4*CA*(2 + 2*xp^2 + 2*xp*(-2 + zp) - zp)*
           Log[xp - zp])/(-1 + zp))*Log[zp] + 
       (4*CA*(-2 - 2*xp^2 - 2*xp*(-2 + zp) + zp)*Log[zp]^2)/(-1 + zp) + 
       Log[1 - zp]*((2*CA*(-2 + xp*(2 - 6*zp) + 9*zp + 6*xp^2*zp - 7*zp^2))/
          ((-1 + zp)*zp) + (4*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*
           Log[1 - xp - zp])/((-1 + zp)*zp) + 
         (4*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*Log[xp - zp])/((-1 + zp)*zp) + 
         4*CA*(3 + xp*(-6 + 2/zp) - (2*xp^2)/(-1 + zp) - zp^(-1))*Log[zp]) + 
       Log[xp]*((-2*CA*(5 - 6*xp + 6*xp^2 - 3*zp))/(-1 + zp) + 
         (4*CA*(2 + 2*xp^2 + 2*xp*(-2 + zp) - zp)*Log[1 - zp])/(-1 + zp) - 
         (4*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*Log[1 - xp - zp])/
          ((-1 + zp)*zp) + (4*CA*(-2 - 2*xp^2 - 2*xp*(-2 + zp) + zp)*
           Log[xp - zp])/(-1 + zp) + (4*CA*(2 + 2*xp^2 + 2*xp*(-2 + zp) - zp)*
           Log[zp])/(-1 + zp)) + Log[1 - xp]*
        ((-4*CA*(-1 + xp - 3*zp + 6*xp*zp - 6*xp^2*zp + zp^2))/
          ((-1 + zp)*zp) + (4*CA*(2 + 2*xp^2 + 2*xp*(-2 + zp) - zp)*Log[xp])/
          (-1 + zp) + (4*CA*(-2 - 2*xp^2 - 2*xp*(-2 + zp) + zp)*Log[1 - zp])/
          (-1 + zp) + (4*CA*(1 + 2*xp^2*zp + zp^2 - 2*xp*(1 + zp^2))*
           Log[xp - zp])/((-1 + zp)*zp) - 
         (4*CA*(1 + 2*zp + 4*xp^2*zp - 2*xp*(1 + 2*zp))*Log[zp])/
          ((-1 + zp)*zp)))*Theta[xp - zp] + 
     ((2*CA*(-1 + 2*xp)*zeta2*(1 - 2*zp + 2*zp^2))/((-1 + zp)*zp) + 
       (4*CA*(1 - 2*zp + 2*zp^2 + xp*(-1 - 4*zp + 4*zp^2)))/((-1 + zp)*zp) - 
       (2*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*
         Li[2, -(xp/((1 - xp - zp)*(-xp + zp))) + 
           xp^2/((1 - xp - zp)*(-xp + zp))])/((-1 + zp)*zp) - 
       (CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*Log[1 - xp]^2)/((-1 + zp)*zp) - 
       (CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*Log[xp]^2)/((-1 + zp)*zp) - 
       (4*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*Log[1 - zp]^2)/((-1 + zp)*zp) - 
       (CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*Log[1 - xp - zp]^2)/
        ((-1 + zp)*zp) - (4*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*Log[zp]^2)/
        ((-1 + zp)*zp) - (4*CA*(-1 + xp + 2*zp - 2*zp^2)*Log[-xp + zp])/
        ((-1 + zp)*zp) - (CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*Log[-xp + zp]^2)/
        ((-1 + zp)*zp) + Log[1 - xp - zp]*((-4*CA*(-1 + xp + 2*zp - 2*zp^2))/
          ((-1 + zp)*zp) - (2*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*
           Log[-xp + zp])/((-1 + zp)*zp)) + 
       Log[xp]*((-4*CA*(-1 + xp + 2*zp - 2*zp^2))/((-1 + zp)*zp) + 
         (4*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*Log[1 - zp])/((-1 + zp)*zp) - 
         (2*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*Log[1 - xp - zp])/
          ((-1 + zp)*zp) + (4*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*Log[zp])/
          ((-1 + zp)*zp) - (2*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*
           Log[-xp + zp])/((-1 + zp)*zp)) + Log[1 - xp]*
        ((4*CA*(-1 + xp + 2*zp - 2*zp^2))/((-1 + zp)*zp) + 
         (2*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*Log[xp])/((-1 + zp)*zp) - 
         (4*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*Log[1 - zp])/((-1 + zp)*zp) + 
         (2*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*Log[1 - xp - zp])/
          ((-1 + zp)*zp) - (4*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*Log[zp])/
          ((-1 + zp)*zp) + (2*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*
           Log[-xp + zp])/((-1 + zp)*zp)) + 
       Log[zp]*((8*CA*(-1 + xp + 2*zp - 2*zp^2))/((-1 + zp)*zp) + 
         (4*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*Log[1 - xp - zp])/
          ((-1 + zp)*zp) + (4*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*
           Log[-xp + zp])/((-1 + zp)*zp)) + Log[1 - zp]*
        ((8*CA*(-1 + xp + 2*zp - 2*zp^2))/((-1 + zp)*zp) + 
         (4*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*Log[1 - xp - zp])/
          ((-1 + zp)*zp) - (8*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*Log[zp])/
          ((-1 + zp)*zp) + (4*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*
           Log[-xp + zp])/((-1 + zp)*zp)))*Theta[-xp + zp]) + 
   ((-2*CA*(-1 + 2*xp)*zeta2*(1 - 2*zp + 2*zp^2))/((-1 + zp)*zp) - 
     (4*CA*(1 - 2*zp + 2*zp^2 + xp*(-1 - 4*zp + 4*zp^2)))/((-1 + zp)*zp) + 
     (2*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*Li[2, (1 - xp - zp)^(-1) - 
         xp/(1 - xp - zp)])/((-1 + zp)*zp) + 
     (2*CA*(2 - 2*xp^2*(-1 + zp) - 2*zp + zp^2 - 2*xp*(2 - 2*zp + zp^2))*
       Li[2, xp^(-1) + zp^(-1) - 1/(xp*zp)])/((-1 + zp)*zp) + 
     (2*CA*(-2 + 2*xp^2*(-1 + zp) + 2*zp - zp^2 + 2*xp*(2 - 2*zp + zp^2))*
       Li[2, 1 - zp^(-1) + xp/zp])/((-1 + zp)*zp) + 
     (CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*Log[1 - xp]^2)/((-1 + zp)*zp) - 
     (4*CA*(1 + 2*xp^2 + zp - 2*xp*(1 + zp))*Log[xp]^2)/zp + 
     (4*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*Log[1 - zp]^2)/((-1 + zp)*zp) + 
     (4*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*Log[zp]^2)/((-1 + zp)*zp) + 
     (4*CA*(-1 + xp + 2*zp - 2*zp^2)*Log[-1 + xp + zp])/((-1 + zp)*zp) + 
     (CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*Log[-1 + xp + zp]^2)/
      ((-1 + zp)*zp) + Log[zp]*((-8*CA*(-1 + xp + 2*zp - 2*zp^2))/
        ((-1 + zp)*zp) - (4*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*
         Log[-1 + xp + zp])/((-1 + zp)*zp)) + 
     Log[1 - zp]*((-8*CA*(-1 + xp + 2*zp - 2*zp^2))/((-1 + zp)*zp) + 
       (8*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*Log[zp])/((-1 + zp)*zp) - 
       (4*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*Log[-1 + xp + zp])/
        ((-1 + zp)*zp)) + Log[1 - xp]*((-4*CA*(-1 + xp + 2*zp - 2*zp^2))/
        ((-1 + zp)*zp) + (4*CA*(1 + 2*xp^2 + zp - 2*xp*(1 + zp))*Log[xp])/
        zp + (4*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*Log[1 - zp])/
        ((-1 + zp)*zp) + (4*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*Log[zp])/
        ((-1 + zp)*zp) - (2*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*
         Log[-1 + xp + zp])/((-1 + zp)*zp)) + 
     Log[xp]*((-2*CA*(xp*(4 - 6*zp) + 6*xp^2*(-1 + zp) + zp*(-5 + 7*zp)))/
        ((-1 + zp)*zp) + (4*CA*(2*xp^2*(-1 + zp) + 2*xp*(2 - 3*zp)*zp + 
          zp*(-2 + 3*zp))*Log[1 - zp])/((-1 + zp)*zp) - 
       (8*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*Log[zp])/((-1 + zp)*zp) + 
       (4*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*Log[-1 + xp + zp])/
        ((-1 + zp)*zp)))*Theta[-1 + xp + zp] + 
   Theta[xp - zp]*((4*CA*(4 + 7*xp^2 - 2*zp - 4*xp*(1 + zp)))/(-1 + zp) - 
     (2*CA*zeta2*(1 + 2*zp + 4*xp^2*zp - 2*xp*(1 + 2*zp)))/((-1 + zp)*zp) + 
     (4*CA*(2 + 2*xp^2 + 2*xp*(-2 + zp) - zp)*Li[2, zp/xp])/(-1 + zp) - 
     (2*CA*(1 + 2*xp^2*zp + zp^2 - 2*xp*(1 + zp^2))*
       Li[2, xp/(1 - zp) - zp/(1 - zp)])/((-1 + zp)*zp) + 
     (2*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*Li[2, (-xp + zp)^(-1) - 
         xp/(-xp + zp)])/((-1 + zp)*zp) + 
     (CA*(3 + 2*zp + 8*xp^2*zp + 2*zp^2 - 2*xp*(3 + 2*zp + 2*zp^2))*
       Log[1 - xp]^2)/((-1 + zp)*zp) + 
     (8*CA*(2 + 2*xp^2 + 2*xp*(-2 + zp) - zp)*Log[xp]^2)/(-1 + zp) + 
     (4*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*Log[1 - zp]^2)/((-1 + zp)*zp) + 
     (2*CA*(5 - 6*xp + 6*xp^2 - 3*zp)*Log[xp - zp])/(-1 + zp) + 
     (CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*Log[xp - zp]^2)/((-1 + zp)*zp) - 
     (4*CA*(5 - 6*xp + 6*xp^2 - 3*zp)*Log[zp])/(-1 + zp) + 
     (4*CA*(2 + 2*xp^2 + 2*xp*(-2 + zp) - zp)*Log[zp]^2)/(-1 + zp) + 
     Log[1 - zp]*((2*CA*(2 - 9*zp - 6*xp^2*zp + 7*zp^2 + xp*(-2 + 6*zp)))/
        ((-1 + zp)*zp) - (4*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*Log[xp - zp])/
        ((-1 + zp)*zp) + 4*CA*(-3 + xp*(6 - 2/zp) + (2*xp^2)/(-1 + zp) + 
         zp^(-1))*Log[zp]) + Log[xp]*((4*CA*(5 - 6*xp + 6*xp^2 - 3*zp))/
        (-1 + zp) + 4*CA*(3 + xp*(-6 + 2/zp) - (2*xp^2)/(-1 + zp) - zp^(-1))*
        Log[1 - zp] - (12*CA*(2 + 2*xp^2 + 2*xp*(-2 + zp) - zp)*Log[zp])/
        (-1 + zp)) + Log[1 - xp]*
      ((4*CA*(-1 + xp - 3*zp + 6*xp*zp - 6*xp^2*zp + zp^2))/((-1 + zp)*zp) - 
       (4*CA*(1 + 2*zp + 4*xp^2*zp - 2*xp*(1 + 2*zp))*Log[xp])/
        ((-1 + zp)*zp) + (4*CA*(2 + 2*xp^2 + 2*xp*(-2 + zp) - zp)*
         Log[1 - zp])/(-1 + zp) - (2*CA*(1 + 2*zp + 4*xp^2*zp - 
          2*xp*(1 + 2*zp))*Log[xp - zp])/((-1 + zp)*zp) + 
       (4*CA*(1 + 2*zp + 4*xp^2*zp - 2*xp*(1 + 2*zp))*Log[zp])/
        ((-1 + zp)*zp)) + ((4*CA*(-7*xp^2 + 2*(-2 + zp) + 4*xp*(1 + zp)))/
        (-1 + zp) + (2*CA*zeta2*(1 + 2*zp + 4*xp^2*zp - 2*xp*(1 + 2*zp)))/
        ((-1 + zp)*zp) + (2*CA*(1 + 2*xp^2*zp + zp^2 - 2*xp*(1 + zp^2))*
         Li[2, (1 - zp)^(-1) - zp/(xp*(1 - zp))])/((-1 + zp)*zp) - 
       (2*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*
         Li[2, -(xp/((1 - xp - zp)*(-xp + zp))) + 
           xp^2/((1 - xp - zp)*(-xp + zp))])/((-1 + zp)*zp) + 
       (CA*(-3 - 2*zp - 8*xp^2*zp - 2*zp^2 + xp*(6 + 4*zp + 4*zp^2))*
         Log[1 - xp]^2)/((-1 + zp)*zp) - (CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*
         Log[xp]^2)/((-1 + zp)*zp) - (4*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*
         Log[1 - zp]^2)/((-1 + zp)*zp) - (CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*
         Log[xp - zp]^2)/((-1 + zp)*zp) + 
       (4*CA*(-2 - 2*xp^2 - 2*xp*(-2 + zp) + zp)*Log[zp]^2)/(-1 + zp) - 
       (4*CA*(-1 + xp + 2*zp - 2*zp^2)*Log[-1 + xp + zp])/((-1 + zp)*zp) - 
       (CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*Log[-1 + xp + zp]^2)/
        ((-1 + zp)*zp) + Log[xp - zp]*((-2*CA*(5 - 6*xp + 6*xp^2 - 3*zp))/
          (-1 + zp) - (2*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*
           Log[-1 + xp + zp])/((-1 + zp)*zp)) + 
       Log[xp]*((-2*CA*(5 - 6*xp + 6*xp^2 - 3*zp))/(-1 + zp) + 
         (4*CA*(2 + 2*xp^2 + 2*xp*(-2 + zp) - zp)*Log[1 - zp])/(-1 + zp) - 
         (2*CA*(1 + 2*zp + 4*xp^2*zp - 2*xp*(1 + 2*zp))*Log[xp - zp])/
          ((-1 + zp)*zp) + (4*CA*(2 + 2*xp^2 + 2*xp*(-2 + zp) - zp)*Log[zp])/
          (-1 + zp) - (2*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*
           Log[-1 + xp + zp])/((-1 + zp)*zp)) + Log[1 - xp]*
        ((-4*CA*(-1 + xp - 3*zp + 6*xp*zp - 6*xp^2*zp + zp^2))/
          ((-1 + zp)*zp) + (2*CA*(1 + 2*zp + 4*xp^2*zp - 2*xp*(1 + 2*zp))*
           Log[xp])/((-1 + zp)*zp) + (4*CA*(-2 - 2*xp^2 - 2*xp*(-2 + zp) + 
            zp)*Log[1 - zp])/(-1 + zp) + (2*CA*(1 + 2*zp + 4*xp^2*zp - 
            2*xp*(1 + 2*zp))*Log[xp - zp])/((-1 + zp)*zp) - 
         (4*CA*(1 + 2*zp + 4*xp^2*zp - 2*xp*(1 + 2*zp))*Log[zp])/
          ((-1 + zp)*zp) + (2*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*
           Log[-1 + xp + zp])/((-1 + zp)*zp)) + 
       Log[zp]*((4*CA*(5 - 6*xp + 6*xp^2 - 3*zp))/(-1 + zp) + 
         (4*CA*(2 + 2*xp^2 + 2*xp*(-2 + zp) - zp)*Log[xp - zp])/(-1 + zp) + 
         (4*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*Log[-1 + xp + zp])/
          ((-1 + zp)*zp)) + Log[1 - zp]*
        ((2*CA*(-2 + xp*(2 - 6*zp) + 9*zp + 6*xp^2*zp - 7*zp^2))/
          ((-1 + zp)*zp) + (4*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*
           Log[xp - zp])/((-1 + zp)*zp) + 4*CA*(3 + xp*(-6 + 2/zp) - 
           (2*xp^2)/(-1 + zp) - zp^(-1))*Log[zp] + 
         (4*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*Log[-1 + xp + zp])/
          ((-1 + zp)*zp)))*Theta[-1 + xp + zp]) + 
   Theta[-xp + zp]*((4*CA*zeta2*(1 + 2*zp + 4*xp^2*zp - 2*xp*(1 + 2*zp)))/
      ((-1 + zp)*zp) - (4*CA*(1 - 2*zp + 2*zp^2 + xp*(-1 - 4*zp + 4*zp^2)))/
      ((-1 + zp)*zp) + (4*CA*(-2 - 2*xp^2 - 2*xp*(-2 + zp) + zp)*
       Li[2, xp/zp])/(-1 + zp) - (2*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*
       Li[2, -(xp/(1 - xp)) + zp/(1 - xp)])/((-1 + zp)*zp) + 
     (2*CA*(1 + 2*xp^2*zp + zp^2 - 2*xp*(1 + zp^2))*
       Li[2, -(-xp + zp)^(-1) + zp/(-xp + zp)])/((-1 + zp)*zp) - 
     (2*CA*(1 + 2*xp^2*zp + zp^2 - 2*xp*(1 + zp^2))*
       Li[2, -(xp/(-xp + zp)) + (xp*zp)/(-xp + zp)])/((-1 + zp)*zp) + 
     (CA*(-1 + 12*zp + 10*xp^2*zp - 7*zp^2 + 2*xp*(1 - 12*zp + 7*zp^2))*
       Log[xp]^2)/((-1 + zp)*zp) + (4*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*
       Log[1 - zp]^2)/((-1 + zp)*zp) - 
     (2*CA*(2 - 2*zp + 2*xp^2*zp + 3*zp^2 + xp*(-4 + 4*zp - 6*zp^2))*
       Log[zp]^2)/((-1 + zp)*zp) + Log[1 - xp]*
      ((-4*CA*(-1 + xp + 2*zp - 2*zp^2))/((-1 + zp)*zp) + 
       (4*CA*(-2 - 2*xp^2 - 2*xp*(-2 + zp) + zp)*Log[xp])/(-1 + zp) + 
       (4*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*Log[1 - zp])/((-1 + zp)*zp) + 
       (4*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*Log[zp])/((-1 + zp)*zp)) + 
     (4*CA*(-1 + xp + 2*zp - 2*zp^2)*Log[-xp + zp])/((-1 + zp)*zp) + 
     Log[1 - zp]*((-8*CA*(-1 + xp + 2*zp - 2*zp^2))/((-1 + zp)*zp) + 
       (8*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*Log[zp])/((-1 + zp)*zp) - 
       (4*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*Log[-xp + zp])/((-1 + zp)*zp)) + 
     Log[xp]*((2*CA*(-2 + xp*(2 - 6*zp) + 9*zp + 6*xp^2*zp - 7*zp^2))/
        ((-1 + zp)*zp) - (2*CA*(-3 + 8*zp + 2*xp^2*zp - 7*zp^2 + 
          2*xp*(3 - 8*zp + 7*zp^2))*Log[1 - zp])/((-1 + zp)*zp) + 
       4*CA*(3 + xp*(-6 + 2/zp) - (2*xp^2)/(-1 + zp) - zp^(-1))*Log[zp] - 
       (2*CA*(1 + 2*xp^2*zp + zp^2 - 2*xp*(1 + zp^2))*Log[-xp + zp])/
        ((-1 + zp)*zp)) + Log[zp]*((-8*CA*(-1 + xp + 2*zp - 2*zp^2))/
        ((-1 + zp)*zp) + (4*CA*(1 + 2*xp^2*zp + zp^2 - 2*xp*(1 + zp^2))*
         Log[-xp + zp])/((-1 + zp)*zp)) + 
     ((4*CA*(-1 + 2*xp)*zeta2*(1 - 2*zp + 2*zp^2))/((-1 + zp)*zp) + 
       (4*CA*(1 - 2*zp + 2*zp^2 + xp*(-1 - 4*zp + 4*zp^2)))/((-1 + zp)*zp) + 
       (2*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*Li[2, (1 - xp)^(-1) - 
           xp/(1 - xp) - zp/((1 - xp)*xp) + zp^2/((1 - xp)*xp)])/
        ((-1 + zp)*zp) - (4*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*Log[1 - zp]^2)/
        ((-1 + zp)*zp) - (4*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*Log[zp]^2)/
        ((-1 + zp)*zp) + Log[1 - xp]*((4*CA*(-1 + xp + 2*zp - 2*zp^2))/
          ((-1 + zp)*zp) + (4*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*Log[xp])/
          ((-1 + zp)*zp) - (4*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*Log[1 - zp])/
          ((-1 + zp)*zp) - (4*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*Log[zp])/
          ((-1 + zp)*zp)) - (4*CA*(-1 + xp + 2*zp - 2*zp^2)*Log[-xp + zp])/
        ((-1 + zp)*zp) - (4*CA*(-1 + xp + 2*zp - 2*zp^2)*Log[-1 + xp + zp])/
        ((-1 + zp)*zp) + Log[xp]*((-4*CA*(-1 + xp + 2*zp - 2*zp^2))/
          ((-1 + zp)*zp) + (4*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*Log[1 - zp])/
          ((-1 + zp)*zp) + (4*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*Log[zp])/
          ((-1 + zp)*zp) - (4*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*
           Log[-xp + zp])/((-1 + zp)*zp) - (4*CA*(-1 + 2*xp)*
           (1 - 2*zp + 2*zp^2)*Log[-1 + xp + zp])/((-1 + zp)*zp)) + 
       Log[zp]*((8*CA*(-1 + xp + 2*zp - 2*zp^2))/((-1 + zp)*zp) + 
         (4*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*Log[-xp + zp])/
          ((-1 + zp)*zp) + (4*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*
           Log[-1 + xp + zp])/((-1 + zp)*zp)) + Log[1 - zp]*
        ((8*CA*(-1 + xp + 2*zp - 2*zp^2))/((-1 + zp)*zp) - 
         (8*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*Log[zp])/((-1 + zp)*zp) + 
         (4*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*Log[-xp + zp])/
          ((-1 + zp)*zp) + (4*CA*(-1 + 2*xp)*(1 - 2*zp + 2*zp^2)*
           Log[-1 + xp + zp])/((-1 + zp)*zp)))*Theta[-1 + xp + zp]) + 
   CF*(26 - 17/(2*zp) - (43*zp)/2 + 2/(-1 + xp + zp) + 
     (1 + 2*xp + xp^2 - 4*xp*zp)^(-1) + xp^3/(1 + 2*xp + xp^2 - 4*xp*zp) + 
     xp*(-9 + 2/zp + 22*zp - 6/(-1 + xp + zp) - (1 + 2*xp + xp^2 - 4*xp*zp)^
        (-1)) + xp^2*(4/(-1 + xp + zp) - (1 + 2*xp + xp^2 - 4*xp*zp)^(-1)) - 
     (4*zeta2*(3 - 6*zp + 2*xp^2*(3 - 6*zp + zp^2) + 
        xp*(2 - 4*zp + 12*zp^2 - 4*zp^3)))/((-1 + zp)*zp) + 
     (-6 + 12*xp)*Dzp[2] + (4*(1 + zp + 2*zp^3 + 2*xp^2*zp*(1 + zp) - 
        2*xp*(1 + 2*zp^2 + zp^3))*Li[2, 1 - xp])/((-1 + zp)*zp) + 
     (8*(2 - 2*xp^2*(-1 + zp) - 2*zp + zp^2 + 2*xp*(2 - 2*zp + zp^2))*
       Li[2, -xp])/zp + (4*(-1 + 2*zp)*(-1 - 4*xp^2 - 2*zp + xp*(2 + 4*zp))*
       Li[2, 1 - zp])/((-1 + zp)*zp) - 
     (4*(-2 + 3*zp - 3*zp^2 + xp^7*(-2 + 4*zp) - 
        2*xp^6*(7 - 23*zp + 23*zp^2) + 7*xp^5*(-6 + 25*zp - 39*zp^2 + 
          26*zp^3) + xp*(-14 + 45*zp - 51*zp^2 + 34*zp^3) - 
        5*xp^4*(14 - 65*zp + 121*zp^2 - 112*zp^3 + 56*zp^4) - 
        2*xp^2*(21 - 89*zp + 146*zp^2 - 114*zp^3 + 57*zp^4) + 
        2*xp^3*(-35 + 164*zp - 311*zp^2 + 304*zp^3 - 145*zp^4 + 58*zp^5))*
       Li[2, 1/2 - xp/2 - Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]/2])/
      ((1 + xp^2 + xp*(2 - 4*zp))^(5/2)*(-1 + zp)*zp) - 
     (4*(-2 + 3*zp - 3*zp^2 + xp^7*(-2 + 4*zp) - 
        2*xp^6*(7 - 23*zp + 23*zp^2) + 7*xp^5*(-6 + 25*zp - 39*zp^2 + 
          26*zp^3) + xp*(-14 + 45*zp - 51*zp^2 + 34*zp^3) - 
        5*xp^4*(14 - 65*zp + 121*zp^2 - 112*zp^3 + 56*zp^4) - 
        2*xp^2*(21 - 89*zp + 146*zp^2 - 114*zp^3 + 57*zp^4) + 
        2*xp^3*(-35 + 164*zp - 311*zp^2 + 304*zp^3 - 145*zp^4 + 58*zp^5))*
       Li[2, 1/2 + xp/2 - Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]/2])/
      ((1 + xp^2 + xp*(2 - 4*zp))^(5/2)*(-1 + zp)*zp) + 
     (4*(-2 + 3*zp - 3*zp^2 + xp^7*(-2 + 4*zp) - 
        2*xp^6*(7 - 23*zp + 23*zp^2) + 7*xp^5*(-6 + 25*zp - 39*zp^2 + 
          26*zp^3) + xp*(-14 + 45*zp - 51*zp^2 + 34*zp^3) - 
        5*xp^4*(14 - 65*zp + 121*zp^2 - 112*zp^3 + 56*zp^4) - 
        2*xp^2*(21 - 89*zp + 146*zp^2 - 114*zp^3 + 57*zp^4) + 
        2*xp^3*(-35 + 164*zp - 311*zp^2 + 304*zp^3 - 145*zp^4 + 58*zp^5))*
       Li[2, 1/2 - 1/(2*xp) - Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]/(2*xp)])/
      ((1 + xp^2 + xp*(2 - 4*zp))^(5/2)*(-1 + zp)*zp) + 
     (4*(-2 + 3*zp - 3*zp^2 + xp^7*(-2 + 4*zp) - 
        2*xp^6*(7 - 23*zp + 23*zp^2) + 7*xp^5*(-6 + 25*zp - 39*zp^2 + 
          26*zp^3) + xp*(-14 + 45*zp - 51*zp^2 + 34*zp^3) - 
        5*xp^4*(14 - 65*zp + 121*zp^2 - 112*zp^3 + 56*zp^4) - 
        2*xp^2*(21 - 89*zp + 146*zp^2 - 114*zp^3 + 57*zp^4) + 
        2*xp^3*(-35 + 164*zp - 311*zp^2 + 304*zp^3 - 145*zp^4 + 58*zp^5))*
       Li[2, 1/2 + 1/(2*xp) - Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]/(2*xp)])/
      ((1 + xp^2 + xp*(2 - 4*zp))^(5/2)*(-1 + zp)*zp) - 
     2*(-1 + 2*xp)*(1 + zp)*Log[1 - xp]^2 + 
     ((96 + 3/(1 + xp^2 + xp*(2 - 4*zp))^(5/2) - 
        (3*xp^6)/(1 + xp^2 + xp*(2 - 4*zp))^(5/2) + 
        5/(1 + xp^2 + xp*(2 - 4*zp))^(3/2) - 
        56/Sqrt[1 + xp^2 + xp*(2 - 4*zp)] + 
        (xp^4*(10 + xp^2 + xp*(2 - 4*zp)))/(1 + xp^2 + xp*(2 - 4*zp))^(5/2) + 
        xp^2*(48 - 9/(1 + xp^2 + xp*(2 - 4*zp))^(5/2) - 
          6/(1 + xp^2 + xp*(2 - 4*zp))^(3/2) - 
          222/Sqrt[1 + xp^2 + xp*(2 - 4*zp)] + 96/(-1 + zp) - 
          96/(Sqrt[1 + xp^2 + xp*(2 - 4*zp)]*(-1 + zp)) + 
          (-64 + 96/Sqrt[1 + xp^2 + xp*(2 - 4*zp)])/zp) + 120/(-1 + zp) - 
        32/(Sqrt[1 + xp^2 + xp*(2 - 4*zp)]*(-1 + zp)) + 
        (32*(-3 + 1/Sqrt[1 + xp^2 + xp*(2 - 4*zp)]))/zp + 
        (32*xp^3*(-1 + 2*zp))/(Sqrt[1 + xp^2 + xp*(2 - 4*zp)]*(-1 + zp)*zp) + 
        xp*((-64 + 116/Sqrt[1 + xp^2 + xp*(2 - 4*zp)])*zp + 
          (128 + 96/Sqrt[1 + 2*xp + xp^2 - 4*xp*zp])/zp + 
          (144 + 96*zp - 154/Sqrt[1 + 2*xp + xp^2 - 4*xp*zp] + 
            (58*zp)/Sqrt[1 + 2*xp + xp^2 - 4*xp*zp])/(1 - zp)))*Log[xp]^2)/
      4 + Dzp[1]*(17 - 20*xp + (-8 + 16*xp)*Log[1 - xp] + 
       (6 - 12*xp)*Log[xp]) + Dzp[0]*(-10 - xp + (8 - 16*xp)*zeta2 + 
       (-2 + 4*xp)*Li[2, 1 - xp] + (-4 + 8*xp)*Log[1 - xp]^2 + 
       2*(-7 + 8*xp)*Log[xp] + (-2 + 4*xp)*Log[xp]^2 + 
       Log[1 - xp]*(17 - 20*xp + (6 - 12*xp)*Log[xp])) + 
     Log[Q2/muF2]^2*(-2*(-1 + 2*xp)*(1 + zp) + (-4 + 8*xp)*Dzp[0] + 
       Delz*(-3/2 + 6*xp + (-2 + 4*xp)*Log[1 - xp] + (1 - 2*xp)*Log[xp])) - 
     (2*(-1 + 2*xp)*(-1 + 3*zp + zp^2)*Log[1 - zp]^2)/zp + 
     ((35 - (6*xp^3)/(1 + xp^2 + xp*(2 - 4*zp))^2 + 
        (3*xp^4)/(1 + xp^2 + xp*(2 - 4*zp))^2 + 
        (3*xp^5)/(1 + xp^2 + xp*(2 - 4*zp))^2 - 16/(-1 + zp) - 6/zp - 18*zp + 
        3/(1 + 2*xp + xp^2 - 4*xp*zp)^2 + 6/(1 + 2*xp + xp^2 - 4*xp*zp) + 
        2*xp^2*(-3/(1 + xp^2 + xp*(2 - 4*zp))^2 + 8/(xp - zp) + 
          (1 + 2*xp + xp^2 - 4*xp*zp)^(-1)) + 
        xp*(-43 + 16/(-1 + zp) + 20*zp + 8/(-xp + zp) + 
          3/(1 + 2*xp + xp^2 - 4*xp*zp)^2 + 8/(1 + 2*xp + xp^2 - 4*xp*zp)))*
       Log[zp])/2 + (2*(-1 + 2*xp)*(-1 + 2*zp + zp^2)*Log[zp]^2)/(-1 + zp) + 
     Log[Q2/muF2]*(3*(-4 - 2*xp*(-2 + zp) + zp) + (-12 + 24*xp)*Dzp[1] - 
       4*(-1 + 2*xp)*(1 + zp)*Log[1 - xp] + (-4 + 8*xp*zp)*Log[xp] + 
       Dzp[0]*(11 - 8*xp + (-8 + 16*xp)*Log[1 - xp] + (6 - 12*xp)*Log[xp]) + 
       Delz*(2 - 13*xp + (8 - 16*xp)*zeta2 + (-2 + 4*xp)*Li[2, 1 - xp] + 
         (-4 + 8*xp)*Log[1 - xp]^2 + 4*(-2 + xp)*Log[xp] + 
         (-2 + 4*xp)*Log[xp]^2 + Log[1 - xp]*(11 - 8*xp + 
           (6 - 12*xp)*Log[xp])) - (4*(-1 + 2*xp)*(-1 + 3*zp + zp^2)*
         Log[1 - zp])/zp + (-4 + 8*xp)*Log[zp]) + 
     Log[1 - zp]*((6 - 40*zp + 59*zp^2 - 28*zp^3 + 3*zp^4 - 
         2*xp^3*(6 - 22*zp + 3*zp^2) + xp^2*(30 - 136*zp + 87*zp^2 - 
           12*zp^3) - 2*xp*(12 - 66*zp + 72*zp^2 - 25*zp^3 + 3*zp^4))/
        (zp*(-1 + xp + zp)^2) + ((8 + 8*xp*(-2 + zp) - 8*xp^2*(-1 + zp) - 
          4*zp)*Log[zp])/((-1 + zp)*zp)) + Log[1 - xp]*
      (-16 + 8*xp^2 + xp*(20 - 12/zp - 6*zp) + 6/zp + 3*zp + 
       ((-8*xp^2 - 12*(1 - zp + zp^2) + 8*xp*(3 - 2*zp + zp^2 + zp^3))*
         Log[xp])/((-1 + zp)*zp) + (4*(2*xp^2 + (1 + zp)^2 - 2*xp*(1 + zp)^2)*
         Log[1 - zp])/zp + (4*(-1 + zp - 2*xp^2*zp - 2*zp^2 + 
          xp*(2 - 2*zp + 4*zp^2))*Log[zp])/((-1 + zp)*zp)) - 
     (4*(-2 + 3*zp - 3*zp^2 + xp^7*(-2 + 4*zp) - 
        2*xp^6*(7 - 23*zp + 23*zp^2) + 7*xp^5*(-6 + 25*zp - 39*zp^2 + 
          26*zp^3) + xp*(-14 + 45*zp - 51*zp^2 + 34*zp^3) - 
        5*xp^4*(14 - 65*zp + 121*zp^2 - 112*zp^3 + 56*zp^4) - 
        2*xp^2*(21 - 89*zp + 146*zp^2 - 114*zp^3 + 57*zp^4) + 
        2*xp^3*(-35 + 164*zp - 311*zp^2 + 304*zp^3 - 145*zp^4 + 58*zp^5))*
       Log[1 + xp - Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]]*
       Log[1 - xp + Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]])/
      ((1 + xp^2 + xp*(2 - 4*zp))^(5/2)*(-1 + zp)*zp) + 
     (4*(-2 + 3*zp - 3*zp^2 + xp^7*(-2 + 4*zp) - 
        2*xp^6*(7 - 23*zp + 23*zp^2) + 7*xp^5*(-6 + 25*zp - 39*zp^2 + 
          26*zp^3) + xp*(-14 + 45*zp - 51*zp^2 + 34*zp^3) - 
        5*xp^4*(14 - 65*zp + 121*zp^2 - 112*zp^3 + 56*zp^4) - 
        2*xp^2*(21 - 89*zp + 146*zp^2 - 114*zp^3 + 57*zp^4) + 
        2*xp^3*(-35 + 164*zp - 311*zp^2 + 304*zp^3 - 145*zp^4 + 58*zp^5))*
       Log[1 + xp - Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]]*
       Log[-1 + xp + Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]])/
      ((1 + xp^2 + xp*(2 - 4*zp))^(5/2)*(-1 + zp)*zp) + 
     Log[2]*((8*(-2 + 3*zp - 3*zp^2 + xp^7*(-2 + 4*zp) - 
          2*xp^6*(7 - 23*zp + 23*zp^2) + 7*xp^5*(-6 + 25*zp - 39*zp^2 + 
            26*zp^3) + xp*(-14 + 45*zp - 51*zp^2 + 34*zp^3) - 
          5*xp^4*(14 - 65*zp + 121*zp^2 - 112*zp^3 + 56*zp^4) - 
          2*xp^2*(21 - 89*zp + 146*zp^2 - 114*zp^3 + 57*zp^4) + 
          2*xp^3*(-35 + 164*zp - 311*zp^2 + 304*zp^3 - 145*zp^4 + 58*zp^5))*
         Log[xp])/((1 + xp^2 + xp*(2 - 4*zp))^(5/2)*(-1 + zp)*zp) + 
       (4*(-2 + 3*zp - 3*zp^2 + xp^7*(-2 + 4*zp) - 
          2*xp^6*(7 - 23*zp + 23*zp^2) + 7*xp^5*(-6 + 25*zp - 39*zp^2 + 
            26*zp^3) + xp*(-14 + 45*zp - 51*zp^2 + 34*zp^3) - 
          5*xp^4*(14 - 65*zp + 121*zp^2 - 112*zp^3 + 56*zp^4) - 
          2*xp^2*(21 - 89*zp + 146*zp^2 - 114*zp^3 + 57*zp^4) + 
          2*xp^3*(-35 + 164*zp - 311*zp^2 + 304*zp^3 - 145*zp^4 + 58*zp^5))*
         Log[1 - xp + Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]])/
        ((1 + xp^2 + xp*(2 - 4*zp))^(5/2)*(-1 + zp)*zp) - 
       (4*(-2 + 3*zp - 3*zp^2 + xp^7*(-2 + 4*zp) - 
          2*xp^6*(7 - 23*zp + 23*zp^2) + 7*xp^5*(-6 + 25*zp - 39*zp^2 + 
            26*zp^3) + xp*(-14 + 45*zp - 51*zp^2 + 34*zp^3) - 
          5*xp^4*(14 - 65*zp + 121*zp^2 - 112*zp^3 + 56*zp^4) - 
          2*xp^2*(21 - 89*zp + 146*zp^2 - 114*zp^3 + 57*zp^4) + 
          2*xp^3*(-35 + 164*zp - 311*zp^2 + 304*zp^3 - 145*zp^4 + 58*zp^5))*
         Log[-1 + xp + Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]])/
        ((1 + xp^2 + xp*(2 - 4*zp))^(5/2)*(-1 + zp)*zp)) + 
     Log[xp]*(((3*xp^4)/(1 + xp^2 + xp*(2 - 4*zp))^2 - 
         (3*xp^5)/(1 + xp^2 + xp*(2 - 4*zp))^2 - 58/zp - 28*zp + 
         xp^3*(6/(1 + xp^2 + xp*(2 - 4*zp))^2 - 8/(-1 + xp + zp)^2) + 
         xp*(3 + 8/(xp - zp) - 80/(-1 + zp) + 84/zp - 4/(-1 + xp + zp)^2 + 
           24/(-1 + xp + zp) - 3/(1 + 2*xp + xp^2 - 4*xp*zp)^2 - 
           8/(1 + 2*xp + xp^2 - 4*xp*zp)) + 2*xp^2*(-8 + 24/(-1 + zp) - 
           24/zp + 8/(-xp + zp) + 6/(-1 + xp + zp)^2 - 8/(-1 + xp + zp) - 
           3/(1 + 2*xp + xp^2 - 4*xp*zp)^2 + (1 + 2*xp + xp^2 - 4*xp*zp)^
            (-1)) + 3*(17 + 16/(-1 + zp) - 4/(-1 + xp + zp) + 
           (1 + 2*xp + xp^2 - 4*xp*zp)^(-2) + 2/(1 + 2*xp + xp^2 - 4*xp*zp)))/
        2 + (8*(2 - 2*xp^2*(-1 + zp) - 2*zp + zp^2 + 2*xp*(2 - 2*zp + zp^2))*
         Log[1 + xp])/zp + (4*(-4 - 2*xp^2 + 5*zp - 5*zp^2 + 
          2*xp*(4 - 4*zp + 3*zp^2 + zp^3))*Log[1 - zp])/((-1 + zp)*zp) - 
       (2*(3 + 7*zp + 8*xp^2*zp - 2*zp^2 + 2*zp^3 - 
          2*xp*(11 - 9*zp + 8*zp^2))*Log[zp])/((-1 + zp)*zp) - 
       (4*(-2 + 3*zp - 3*zp^2 + xp^7*(-2 + 4*zp) - 
          2*xp^6*(7 - 23*zp + 23*zp^2) + 7*xp^5*(-6 + 25*zp - 39*zp^2 + 
            26*zp^3) + xp*(-14 + 45*zp - 51*zp^2 + 34*zp^3) - 
          5*xp^4*(14 - 65*zp + 121*zp^2 - 112*zp^3 + 56*zp^4) - 
          2*xp^2*(21 - 89*zp + 146*zp^2 - 114*zp^3 + 57*zp^4) + 
          2*xp^3*(-35 + 164*zp - 311*zp^2 + 304*zp^3 - 145*zp^4 + 58*zp^5))*
         Log[-1 + xp + Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]])/
        ((1 + xp^2 + xp*(2 - 4*zp))^(5/2)*(-1 + zp)*zp) - 
       (4*(-2 + 3*zp - 3*zp^2 + xp^7*(-2 + 4*zp) - 
          2*xp^6*(7 - 23*zp + 23*zp^2) + 7*xp^5*(-6 + 25*zp - 39*zp^2 + 
            26*zp^3) + xp*(-14 + 45*zp - 51*zp^2 + 34*zp^3) - 
          5*xp^4*(14 - 65*zp + 121*zp^2 - 112*zp^3 + 56*zp^4) - 
          2*xp^2*(21 - 89*zp + 146*zp^2 - 114*zp^3 + 57*zp^4) + 
          2*xp^3*(-35 + 164*zp - 311*zp^2 + 304*zp^3 - 145*zp^4 + 58*zp^5))*
         Log[1 + xp + Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]])/
        ((1 + xp^2 + xp*(2 - 4*zp))^(5/2)*(-1 + zp)*zp)) + 
     Delz*((-5*(-7 + 9*xp))/2 + (-13 + 16*xp)*zeta2 + (-22 + 44*xp)*zeta3 + 
       (3 - 4*xp)*Li[2, 1 - xp] + (12 - 24*xp)*Li[3, 1 - xp] + 
       (5*(-1 + 2*xp)*Log[1 - xp]^3)/3 + ((38 + xp)/2 + (-6 + 12*xp)*zeta2 + 
         (-6 + 12*xp)*Li[2, 1 - xp])*Log[xp] + (37/4 - 6*xp)*Log[xp]^2 + 
       ((1 - 2*xp)*Log[xp]^3)/6 + Log[1 - xp]^2*(21/2 - 12*xp + 
         (4 - 8*xp)*Log[xp]) + Log[1 - xp]*(-14 + 4*xp + (8 - 16*xp)*zeta2 + 
         (-4 + 8*xp)*Li[2, 1 - xp] + (-18 + 20*xp)*Log[xp] + 
         (-4 + 8*xp)*Log[xp]^2) + (-10 + 20*xp)*S12[1 - xp]) + 
     Theta[1 - xp - zp]*(((-8 + 8*xp^2*(-1 + zp) + 8*zp - 4*zp^2 + 
          8*xp*(2 - 2*zp + zp^2))*Li[2, -(zp/(1 - xp - zp))])/
        ((-1 + zp)*zp) + (4*(2 - 2*xp^2*(-1 + zp) - 2*zp + zp^2 - 
          2*xp*(2 - 2*zp + zp^2))*Li[2, -((xp*zp)/(1 - xp - zp))])/
        ((-1 + zp)*zp) - (8*(-2 + 2*xp^2*(-1 + zp) + 2*zp - zp^2 + 
          2*xp*(2 - 2*zp + zp^2))*Log[1 - xp]*Log[xp])/((-1 + zp)*zp) + 
       (6*(-2 + 2*xp^2*(-1 + zp) + 2*zp - zp^2 + 2*xp*(2 - 2*zp + zp^2))*
         Log[xp]^2)/((-1 + zp)*zp) + Log[xp]*
        ((-4*(4 - 6*xp^2*(-1 + zp) - 3*zp + zp^2 + xp*(-8 + 6*zp)))/
          ((-1 + zp)*zp) - (8*(-2 + 2*xp^2*(-1 + zp) + 2*zp - zp^2 + 
            2*xp*(2 - 2*zp + zp^2))*Log[1 - zp])/((-1 + zp)*zp) + 
         ((-8 + 8*xp^2*(-1 + zp) + 8*zp - 4*zp^2 + 8*xp*(2 - 2*zp + zp^2))*
           Log[1 - xp - zp])/((-1 + zp)*zp) + 
         (4*(2 - 2*xp^2*(-1 + zp) - 2*zp + zp^2 - 2*xp*(2 - 2*zp + zp^2))*
           Log[zp])/((-1 + zp)*zp)) + 
       ((8*(1 + 2*zp + 7*xp^2*zp - xp*(1 + 8*zp)))/((-1 + zp)*zp) - 
         (8*zeta2*(1 + 2*xp^2*zp + zp^2 - 2*xp*(1 + zp^2)))/((-1 + zp)*zp) - 
         (4*(1 + 2*xp^2*zp + zp^2 - 2*xp*(1 + zp^2))*
           Li[2, (1 - zp)^(-1) - zp/(xp*(1 - zp))])/((-1 + zp)*zp) + 
         (8*(1 + 2*xp^2*zp + zp^2 - 2*xp*(1 + zp^2))*Log[1 - xp]^2)/
          ((-1 + zp)*zp) + (4*(2 + zp + 6*xp^2*zp + zp^2 - 2*xp*(1 + 3*zp))*
           Log[xp - zp])/((-1 + zp)*zp) + 
         ((-8*(2 + zp + 6*xp^2*zp + zp^2 - 2*xp*(1 + 3*zp)))/((-1 + zp)*zp) - 
           (8*(1 + 2*xp^2*zp + zp^2 - 2*xp*(1 + zp^2))*Log[xp - zp])/
            ((-1 + zp)*zp))*Log[zp] + (8*(1 + 2*xp^2*zp + zp^2 - 
            2*xp*(1 + zp^2))*Log[zp]^2)/((-1 + zp)*zp) + 
         Log[xp]*((4*(2 + zp + 6*xp^2*zp + zp^2 - 2*xp*(1 + 3*zp)))/
            ((-1 + zp)*zp) - (8*(1 + 2*xp^2*zp + zp^2 - 2*xp*(1 + zp^2))*
             Log[1 - zp])/((-1 + zp)*zp) + (8*(1 + 2*xp^2*zp + zp^2 - 
              2*xp*(1 + zp^2))*Log[xp - zp])/((-1 + zp)*zp) - 
           (8*(1 + 2*xp^2*zp + zp^2 - 2*xp*(1 + zp^2))*Log[zp])/
            ((-1 + zp)*zp)) + Log[1 - zp]*
          ((-4*(2 + zp + 6*xp^2*zp + zp^2 - 2*xp*(1 + 3*zp)))/
            ((-1 + zp)*zp) + (8*(1 + 2*xp^2*zp + zp^2 - 2*xp*(1 + zp^2))*
             Log[zp])/((-1 + zp)*zp)) + Log[1 - xp]*
          ((-8*(2 + zp + 6*xp^2*zp + zp^2 - 2*xp*(1 + 3*zp)))/
            ((-1 + zp)*zp) - (8*(1 + 2*xp^2*zp + zp^2 - 2*xp*(1 + zp^2))*
             Log[xp])/((-1 + zp)*zp) + (8*(1 + 2*xp^2*zp + zp^2 - 
              2*xp*(1 + zp^2))*Log[1 - zp])/((-1 + zp)*zp) - 
           (8*(1 + 2*xp^2*zp + zp^2 - 2*xp*(1 + zp^2))*Log[xp - zp])/
            ((-1 + zp)*zp) + (16*(1 + 2*xp^2*zp + zp^2 - 2*xp*(1 + zp^2))*
             Log[zp])/((-1 + zp)*zp)))*Theta[xp - zp]) + 
     ((-8*zeta2*(1 + zp + 4*xp^2*zp + zp^2 - 2*xp*(1 + zp + zp^2)))/
        ((-1 + zp)*zp) + (4*(1 + zp + 4*xp^2*zp + zp^2 - 
          2*xp*(1 + zp + zp^2))*Li[2, xp/zp])/((-1 + zp)*zp) - 
       (4*(1 + 2*xp^2*zp + zp^2 - 2*xp*(1 + zp^2))*
         Li[2, -(-xp + zp)^(-1) + zp/(-xp + zp)])/((-1 + zp)*zp) + 
       (4*(1 + 2*xp^2*zp + zp^2 - 2*xp*(1 + zp^2))*
         Li[2, -(xp/(-xp + zp)) + (xp*zp)/(-xp + zp)])/((-1 + zp)*zp) + 
       (8*(1 + 2*xp^2*zp + zp^2 - 2*xp*(1 + zp^2))*Log[1 - xp]*Log[xp])/
        ((-1 + zp)*zp) - (2*(4 + zp + 10*xp^2*zp + 4*zp^2 - 
          2*xp*(4 + zp + 4*zp^2))*Log[xp]^2)/((-1 + zp)*zp) + 
       (2*(1 + zp + 4*xp^2*zp + zp^2 - 2*xp*(1 + zp + zp^2))*Log[zp]^2)/
        ((-1 + zp)*zp) - (4*(1 + zp + 4*xp^2*zp + zp^2 - 
          2*xp*(1 + zp + zp^2))*Log[zp]*Log[-xp + zp])/((-1 + zp)*zp) + 
       Log[xp]*((-4*(2 + zp + 6*xp^2*zp + zp^2 - 2*xp*(1 + 3*zp)))/
          ((-1 + zp)*zp) + (4*(1 + 2*xp^2*zp + zp^2 - 2*xp*(1 + zp^2))*
           Log[1 - zp])/((-1 + zp)*zp) + 
         (8*(1 + 2*xp^2*zp + zp^2 - 2*xp*(1 + zp^2))*Log[zp])/
          ((-1 + zp)*zp) + ((4 - 8*xp + 8*xp^2)*Log[-xp + zp])/(-1 + zp)))*
      Theta[-xp + zp] + 
     (((-8 + 8*xp^2*(-1 + zp) + 8*zp - 4*zp^2 + 8*xp*(2 - 2*zp + zp^2))*
         Li[2, xp^(-1) + zp^(-1) - 1/(xp*zp)])/((-1 + zp)*zp) + 
       (4*(2 - 2*xp^2*(-1 + zp) - 2*zp + zp^2 - 2*xp*(2 - 2*zp + zp^2))*
         Li[2, 1 - zp^(-1) + xp/zp])/((-1 + zp)*zp) - 
       (8*(-2 + 2*xp^2*(-1 + zp) + 2*zp - zp^2 + 2*xp*(2 - 2*zp + zp^2))*
         Log[1 - xp]*Log[xp])/((-1 + zp)*zp) + 
       (8*(-2 + 2*xp^2*(-1 + zp) + 2*zp - zp^2 + 2*xp*(2 - 2*zp + zp^2))*
         Log[xp]^2)/((-1 + zp)*zp) + Log[xp]*
        ((-4*(4 - 6*xp^2*(-1 + zp) - 3*zp + zp^2 + xp*(-8 + 6*zp)))/
          ((-1 + zp)*zp) - (8*(-2 + 2*xp^2*(-1 + zp) + 2*zp - zp^2 + 
            2*xp*(2 - 2*zp + zp^2))*Log[1 - zp])/((-1 + zp)*zp)))*
      Theta[-1 + xp + zp] + Theta[xp - zp]*
      ((-8*(1 + 2*zp + 7*xp^2*zp - xp*(1 + 8*zp)))/((-1 + zp)*zp) + 
       (8*zeta2*(1 + 2*xp^2*zp + zp^2 - 2*xp*(1 + zp^2)))/((-1 + zp)*zp) - 
       (4*(1 + zp + 4*xp^2*zp + zp^2 - 2*xp*(1 + zp + zp^2))*Li[2, zp/xp])/
        ((-1 + zp)*zp) + (4*(1 + 2*xp^2*zp + zp^2 - 2*xp*(1 + zp^2))*
         Li[2, xp/(1 - zp) - zp/(1 - zp)])/((-1 + zp)*zp) - 
       (8*(1 + 2*xp^2*zp + zp^2 - 2*xp*(1 + zp^2))*Log[1 - xp]^2)/
        ((-1 + zp)*zp) - (4*(3 + zp + 8*xp^2*zp + 3*zp^2 - 
          2*xp*(3 + zp + 3*zp^2))*Log[xp]^2)/((-1 + zp)*zp) - 
       (4*(2 + zp + 6*xp^2*zp + zp^2 - 2*xp*(1 + 3*zp))*Log[xp - zp])/
        ((-1 + zp)*zp) + ((8*(2 + zp + 6*xp^2*zp + zp^2 - 2*xp*(1 + 3*zp)))/
          ((-1 + zp)*zp) - (4*(-1 + 2*xp)*(1 - zp + zp^2)*Log[xp - zp])/
          ((-1 + zp)*zp))*Log[zp] - 
       (8*(1 + 2*xp^2*zp + zp^2 - 2*xp*(1 + zp^2))*Log[zp]^2)/
        ((-1 + zp)*zp) + Log[1 - xp]*
        ((8*(2 + zp + 6*xp^2*zp + zp^2 - 2*xp*(1 + 3*zp)))/((-1 + zp)*zp) + 
         (16*(1 + 2*xp^2*zp + zp^2 - 2*xp*(1 + zp^2))*Log[xp])/
          ((-1 + zp)*zp) - (8*(1 + 2*xp^2*zp + zp^2 - 2*xp*(1 + zp^2))*
           Log[1 - zp])/((-1 + zp)*zp) + 
         (8*(1 + 2*xp^2*zp + zp^2 - 2*xp*(1 + zp^2))*Log[xp - zp])/
          ((-1 + zp)*zp) - (16*(1 + 2*xp^2*zp + zp^2 - 2*xp*(1 + zp^2))*
           Log[zp])/((-1 + zp)*zp)) + Log[1 - zp]*
        ((4*(2 + zp + 6*xp^2*zp + zp^2 - 2*xp*(1 + 3*zp)))/((-1 + zp)*zp) - 
         (8*(1 + 2*xp^2*zp + zp^2 - 2*xp*(1 + zp^2))*Log[zp])/
          ((-1 + zp)*zp)) + Log[xp]*
        ((-8*(2 + zp + 6*xp^2*zp + zp^2 - 2*xp*(1 + 3*zp)))/((-1 + zp)*zp) + 
         (8*(1 + 2*xp^2*zp + zp^2 - 2*xp*(1 + zp^2))*Log[1 - zp])/
          ((-1 + zp)*zp) + (4*(-1 + 2*xp)*(1 - zp + zp^2)*Log[xp - zp])/
          ((-1 + zp)*zp) + (4*(5 + zp + 12*xp^2*zp + 5*zp^2 - 
            2*xp*(5 + zp + 5*zp^2))*Log[zp])/((-1 + zp)*zp)) + 
       ((8*(1 + 2*zp + 7*xp^2*zp - xp*(1 + 8*zp)))/((-1 + zp)*zp) - 
         (8*zeta2*(1 + 2*xp^2*zp + zp^2 - 2*xp*(1 + zp^2)))/((-1 + zp)*zp) - 
         (4*(1 + 2*xp^2*zp + zp^2 - 2*xp*(1 + zp^2))*
           Li[2, (1 - zp)^(-1) - zp/(xp*(1 - zp))])/((-1 + zp)*zp) + 
         (8*(1 + 2*xp^2*zp + zp^2 - 2*xp*(1 + zp^2))*Log[1 - xp]^2)/
          ((-1 + zp)*zp) + (4*(2 + zp + 6*xp^2*zp + zp^2 - 2*xp*(1 + 3*zp))*
           Log[xp - zp])/((-1 + zp)*zp) + 
         ((-8*(2 + zp + 6*xp^2*zp + zp^2 - 2*xp*(1 + 3*zp)))/((-1 + zp)*zp) - 
           (8*(1 + 2*xp^2*zp + zp^2 - 2*xp*(1 + zp^2))*Log[xp - zp])/
            ((-1 + zp)*zp))*Log[zp] + (8*(1 + 2*xp^2*zp + zp^2 - 
            2*xp*(1 + zp^2))*Log[zp]^2)/((-1 + zp)*zp) + 
         Log[xp]*((4*(2 + zp + 6*xp^2*zp + zp^2 - 2*xp*(1 + 3*zp)))/
            ((-1 + zp)*zp) - (8*(1 + 2*xp^2*zp + zp^2 - 2*xp*(1 + zp^2))*
             Log[1 - zp])/((-1 + zp)*zp) + (8*(1 + 2*xp^2*zp + zp^2 - 
              2*xp*(1 + zp^2))*Log[xp - zp])/((-1 + zp)*zp) - 
           (8*(1 + 2*xp^2*zp + zp^2 - 2*xp*(1 + zp^2))*Log[zp])/
            ((-1 + zp)*zp)) + Log[1 - zp]*
          ((-4*(2 + zp + 6*xp^2*zp + zp^2 - 2*xp*(1 + 3*zp)))/
            ((-1 + zp)*zp) + (8*(1 + 2*xp^2*zp + zp^2 - 2*xp*(1 + zp^2))*
             Log[zp])/((-1 + zp)*zp)) + Log[1 - xp]*
          ((-8*(2 + zp + 6*xp^2*zp + zp^2 - 2*xp*(1 + 3*zp)))/
            ((-1 + zp)*zp) - (8*(1 + 2*xp^2*zp + zp^2 - 2*xp*(1 + zp^2))*
             Log[xp])/((-1 + zp)*zp) + (8*(1 + 2*xp^2*zp + zp^2 - 
              2*xp*(1 + zp^2))*Log[1 - zp])/((-1 + zp)*zp) - 
           (8*(1 + 2*xp^2*zp + zp^2 - 2*xp*(1 + zp^2))*Log[xp - zp])/
            ((-1 + zp)*zp) + (16*(1 + 2*xp^2*zp + zp^2 - 2*xp*(1 + zp^2))*
             Log[zp])/((-1 + zp)*zp)))*Theta[-1 + xp + zp])), 
 "\!\(\*SubscriptBox[SuperscriptBox[\(G\), \((2)\)], \({1,gg,MS}\)]\)" -> 
  (-8*CA*(-9 + 11*xp)*(-1 + zp))/zp + 
   (8*CA*zeta2*(1 + zp^2 + 2*xp^2*zp^2 + 2*xp*(zp + zp^3)))/
    (zp*(-1 + zp^2)) - (4*CA*(zp + 5*xp^2*zp + 5*xp*(1 + zp^2))*
     InvTanInt[-(Sqrt[zp]/Sqrt[xp])])/(Sqrt[xp]*zp^(3/2)) + 
   (4*CA*(zp + 5*xp^2*zp + 5*xp*(1 + zp^2))*InvTanInt[Sqrt[zp]/Sqrt[xp]])/
    (Sqrt[xp]*zp^(3/2)) + (4*CA*(zp + 5*xp^2*zp + 5*xp*(1 + zp^2))*
     InvTanInt[-(Sqrt[xp]*Sqrt[zp])])/(Sqrt[xp]*zp^(3/2)) - 
   (4*CA*(zp + 5*xp^2*zp + 5*xp*(1 + zp^2))*InvTanInt[Sqrt[xp]*Sqrt[zp]])/
    (Sqrt[xp]*zp^(3/2)) + (8*CA*(-1 + 2*xp^2*zp - zp^2 + 2*xp*(1 + zp^2))*
     Li[2, 1 - xp])/(zp*(1 + zp)) + 
   (8*CA*(1 + 2*xp^2*zp + zp^2 + 2*xp*(1 + zp^2))*Li[2, -xp])/
    ((-1 + zp)*zp) + (8*CA*(1 + zp^2 + 2*xp^2*zp^2 + 2*xp*(zp + zp^3))*
     Li[2, -(zp/xp)])/(zp*(-1 + zp^2)) - 
   (8*CA*(zp + 2*xp^2*zp + zp^3 + 2*xp*(1 + zp^2))*Li[2, -(xp*zp)])/
    (zp*(-1 + zp^2)) + (8*CA*(-1 + 2*xp^2*zp - zp^2 + 2*xp*(1 + zp^2))*
     Li[2, (2 - 2*xp)^(-1) + zp/(2 - 2*xp) - Sqrt[1 - 2*zp + 4*xp*zp + zp^2]/
        (2 - 2*xp)])/(zp*(1 + zp)) - 
   (8*CA*(-1 + 2*xp^2*zp - zp^2 + 2*xp*(1 + zp^2))*
     Li[2, (2 - 2*xp)^(-1) - (2*xp)/(2 - 2*xp) - zp/(2 - 2*xp) + 
       Sqrt[1 - 2*zp + 4*xp*zp + zp^2]/(2 - 2*xp)])/(zp*(1 + zp)) + 
   (16*CA*(-1 + 2*xp^2*zp - zp^2 + 2*xp*(1 + zp^2))*
     Li[2, 1 - 1/(2*xp) + zp/(2*xp) - Sqrt[1 - 2*zp + 4*xp*zp + zp^2]/
        (2*xp)])/(zp*(1 + zp)) + 
   (4*CA*(-1 + 2*xp^2*zp - zp^2 + 2*xp*(1 + zp^2))*Log[2]^2)/(zp*(1 + zp)) + 
   (4*CA*(1 - 2*zp + 2*xp^2*(-2 + zp)*zp + zp^2 - 2*zp^3 + 
      2*xp*(-2 + zp - 2*zp^2 + zp^3))*Log[xp]^2)/(zp*(-1 + zp^2)) + 
   (4*CA*(1 + 2*xp^2 + zp^2 - 2*xp*(1 + zp^2))*Log[zp]^2)/(-1 + zp^2) + 
   Log[zp]*((28*CA*(-1 + xp)*(1 + zp))/zp - 
     (4*CA*(zp + 5*xp^2*zp + 5*xp*(1 + zp^2))*ArcTan[Sqrt[zp]/Sqrt[xp]])/
      (Sqrt[xp]*zp^(3/2)) + (4*CA*(zp + 5*xp^2*zp + 5*xp*(1 + zp^2))*
       ArcTan[Sqrt[xp]*Sqrt[zp]])/(Sqrt[xp]*zp^(3/2)) + 
     (8*CA*(1 + zp^2 + 2*xp^2*zp^2 + 2*xp*(zp + zp^3))*Log[1 + zp/xp])/
      (zp*(-1 + zp^2)) - (8*CA*(zp + 2*xp^2*zp + zp^3 + 2*xp*(1 + zp^2))*
       Log[1 + xp*zp])/(zp*(-1 + zp^2))) + 
   (4*CA*(-1 + 2*xp^2*zp - zp^2 + 2*xp*(1 + zp^2))*
     Log[1 + zp - Sqrt[1 - 2*zp + 4*xp*zp + zp^2]]^2)/(zp*(1 + zp)) + 
   (4*CA*(-1 + 2*xp^2*zp - zp^2 + 2*xp*(1 + zp^2))*
     Log[1 - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]]^2)/(zp*(1 + zp)) + 
   Log[1 + zp - Sqrt[1 - 2*zp + 4*xp*zp + zp^2]]*
    ((8*CA*(-1 + 2*xp)*Sqrt[1 + (-2 + 4*xp)*zp + zp^2])/zp - 
     (8*CA*(-1 + 2*xp^2*zp - zp^2 + 2*xp*(1 + zp^2))*
       Log[1 - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/(zp*(1 + zp))) + 
   Log[2]*((8*CA*(-1 + 2*xp)*Sqrt[1 + (-2 + 4*xp)*zp + zp^2])/zp + 
     (8*CA*(-1 + 2*xp^2*zp - zp^2 + 2*xp*(1 + zp^2))*Log[xp])/(zp*(1 + zp)) - 
     (8*CA*(-1 + 2*xp^2*zp - zp^2 + 2*xp*(1 + zp^2))*
       Log[1 - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/(zp*(1 + zp))) - 
   (8*CA*(-1 + 2*xp)*Sqrt[1 + (-2 + 4*xp)*zp + zp^2]*
     Log[1 - 2*xp - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/zp + 
   (4*CA*(1 - 2*xp^2*zp + zp^2 - 2*xp*(1 + zp^2))*
     Log[1 - 2*xp - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]]^2)/(zp*(1 + zp)) + 
   Log[xp]*((4*CA*(16*xp^2*zp + 2*(-1 + zp)*(1 - zp + 
          Sqrt[1 - 2*zp + 4*xp*zp + zp^2]) + 
        xp*(4 + 4*zp^2 - 9*Sqrt[1 - 2*zp + 4*xp*zp + zp^2] + 
          zp*(-16 + 9*Sqrt[1 - 2*zp + 4*xp*zp + zp^2]))))/
      (zp*Sqrt[1 + (-2 + 4*xp)*zp + zp^2]) + 
     (4*CA*(zp + 5*xp^2*zp + 5*xp*(1 + zp^2))*ArcTan[Sqrt[zp]/Sqrt[xp]])/
      (Sqrt[xp]*zp^(3/2)) + (4*CA*(zp + 5*xp^2*zp + 5*xp*(1 + zp^2))*
       ArcTan[Sqrt[xp]*Sqrt[zp]])/(Sqrt[xp]*zp^(3/2)) + 
     (8*CA*(1 + 2*xp^2*zp + zp^2 + 2*xp*(1 + zp^2))*Log[1 + xp])/
      ((-1 + zp)*zp) + (8*CA*(1 + 2*zp + 2*xp^2*zp - zp^2 + 
        xp*(4 + 2*zp - 2*zp^3))*Log[zp])/(zp*(-1 + zp^2)) - 
     (8*CA*(1 + zp^2 + 2*xp^2*zp^2 + 2*xp*(zp + zp^3))*Log[1 + zp/xp])/
      (zp*(-1 + zp^2)) - (8*CA*(zp + 2*xp^2*zp + zp^3 + 2*xp*(1 + zp^2))*
       Log[1 + xp*zp])/(zp*(-1 + zp^2)) + 
     (8*CA*(-1 + 2*xp^2*zp - zp^2 + 2*xp*(1 + zp^2))*
       Log[1 + zp - Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/(zp*(1 + zp)) - 
     (8*CA*(-1 + 2*xp^2*zp - zp^2 + 2*xp*(1 + zp^2))*
       Log[1 - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/(zp*(1 + zp)) - 
     (8*CA*(-1 + 2*xp^2*zp - zp^2 + 2*xp*(1 + zp^2))*
       Log[1 - 2*xp - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/(zp*(1 + zp))) + 
   Log[1 - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]]*
    ((-8*CA*(-1 + 2*xp)*Sqrt[1 + (-2 + 4*xp)*zp + zp^2])/zp + 
     (8*CA*(-1 + 2*xp^2*zp - zp^2 + 2*xp*(1 + zp^2))*
       Log[1 - 2*xp - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/(zp*(1 + zp))) + 
   Log[1 - xp]*((-8*CA*(-1 + 2*xp^2*zp - zp^2 + 2*xp*(1 + zp^2))*
       Log[1 + zp - Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/(zp*(1 + zp)) + 
     (8*CA*(-1 + 2*xp^2*zp - zp^2 + 2*xp*(1 + zp^2))*
       Log[1 - 2*xp - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/(zp*(1 + zp))) + 
   CF*(-14 + xp*(40 - 12/zp - 44*zp) - 13/zp + 43*zp - 
     (8*zeta2*(1 + zp^2 + 2*zp^4 - 2*xp^2*zp*(-1 - 2*zp + zp^2) + 
        2*xp*(-1 + 2*zp + zp^2 + 2*zp^3)))/(zp*(-1 + zp^2)) - 
     32*Sqrt[xp]*Sqrt[zp]*InvTanInt[-(Sqrt[zp]/Sqrt[xp])] + 
     32*Sqrt[xp]*Sqrt[zp]*InvTanInt[Sqrt[zp]/Sqrt[xp]] + 
     32*Sqrt[xp]*Sqrt[zp]*InvTanInt[-(Sqrt[xp]*Sqrt[zp])] - 
     32*Sqrt[xp]*Sqrt[zp]*InvTanInt[Sqrt[xp]*Sqrt[zp]] + 
     (4*(1 - zp + 4*xp^2*(-1 + zp)*zp + 2*zp^2 - 4*zp^3 + 
        2*xp*(-1 + zp - 4*zp^2 + 2*zp^3))*Li[2, 1 - xp])/(zp*(1 + zp)) - 
     (16*(1 + 2*xp^2*zp + zp^2 + 2*xp*(1 + zp^2))*Li[2, -xp])/(-1 + zp) + 
     (16 - 32*xp)*Li[2, 1 - zp] - 
     (16*(zp + 2*xp^2*zp + zp^3 + 2*xp*(1 + zp^2))*Li[2, -(zp/xp)])/
      (-1 + zp^2) + (16*(1 + zp^2 + 2*xp^2*zp^2 + 2*xp*(zp + zp^3))*
       Li[2, -(xp*zp)])/(-1 + zp^2) - 
     (16*(-1 - 2*xp*(-1 + zp) + zp + 2*xp^2*zp)*
       Li[2, (2 - 2*xp)^(-1) + zp/(2 - 2*xp) - 
         Sqrt[1 - 2*zp + 4*xp*zp + zp^2]/(2 - 2*xp)])/(zp*(1 + zp)) - 
     (16*(1 + 2*xp^2 + 2*xp*(-1 + zp) - zp)*zp*
       Li[2, (2 - 2*xp)^(-1) - (2*xp)/(2 - 2*xp) - zp/(2 - 2*xp) + 
         Sqrt[1 - 2*zp + 4*xp*zp + zp^2]/(2 - 2*xp)])/(1 + zp) + 
     (16*(-1 + zp)*(-1 + 2*xp^2*zp - zp^2 + 2*xp*(1 + zp^2))*
       Li[2, 1 - 1/(2*xp) + zp/(2*xp) - Sqrt[1 - 2*zp + 4*xp*zp + zp^2]/
          (2*xp)])/(zp*(1 + zp)) - (8*(-1 - 2*xp*(-1 + zp) + zp + 2*xp^2*zp)*
       Log[2]^2)/(zp*(1 + zp)) + (4*(-1 + 2*xp)*(2 - 2*zp + zp^2)*
       Log[Q2/muF2]^2)/zp + 4*(-2 - 4*xp^2 - 2*xp*(-2 + zp) + zp)*
      Log[1 - xp]^2 + (8*(1 - 2*zp + 3*zp^2 + 4*xp*zp*(1 - zp + zp^2) + 
        xp^2*(1 - 2*zp + 3*zp^2))*Log[xp]^2)/(-1 + zp^2) + 
     (4*(-1 + 2*xp)*(2 - 2*zp + zp^2)*Log[1 - zp]^2)/zp - 
     (4*(4 + 2*zp - zp^2 + 2*zp^3 + 4*xp^2*zp^3 - 3*zp^4 + 
        2*xp*(-4 - 2*zp + zp^2 - 2*zp^3 + 3*zp^4))*Log[zp]^2)/
      (zp*(-1 + zp^2)) + Log[1 - zp]*(-44 + 50/zp - 6*zp + 
       xp*(64 - 72/zp + 12*zp) + (16 - 32*xp)*Log[zp]) + 
     Log[Q2/muF2]*(-44 + 50/zp - 6*zp + xp*(64 - 72/zp + 12*zp) + 
       (8*(-1 + 2*xp)*(2 - 2*zp + zp^2)*Log[1 - xp])/zp - 
       (4*(-3 + 2*zp + xp*(6 - 4*zp + 4*zp^2))*Log[xp])/zp + 
       (8*(-1 + 2*xp)*(2 - 2*zp + zp^2)*Log[1 - zp])/zp + 
       (16 - 32*xp)*Log[zp]) + Log[zp]*
      (2*(-6 + (33 - 52*xp)/zp + (9 - 10*xp)*zp) - 32*Sqrt[xp]*Sqrt[zp]*
        ArcTan[Sqrt[zp]/Sqrt[xp]] + 32*Sqrt[xp]*Sqrt[zp]*
        ArcTan[Sqrt[xp]*Sqrt[zp]] - (16*(zp + 2*xp^2*zp + zp^3 + 
          2*xp*(1 + zp^2))*Log[1 + zp/xp])/(-1 + zp^2) + 
       (16*(1 + zp^2 + 2*xp^2*zp^2 + 2*xp*(zp + zp^3))*Log[1 + xp*zp])/
        (-1 + zp^2)) + (8*(1 + 2*xp^2 + 2*xp*(-1 + zp) - zp)*zp*
       Log[1 + zp - Sqrt[1 - 2*zp + 4*xp*zp + zp^2]]^2)/(1 + zp) + 
     (8*(1 + 2*xp^2 + 2*xp*(-1 + zp) - zp)*zp*
       Log[1 - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]]^2)/(1 + zp) + 
     Log[2]*((-16*(-(-1 + zp)^2 + 2*xp*(-1 + zp)^2 + 2*xp^2*zp)*Log[1 - xp])/
        zp + (16*(1 + 2*xp^2 + 2*xp*(-1 + zp) - zp)*zp*Log[xp])/(1 + zp) + 
       (16*(-(-1 + zp)^2 + 2*xp*(-1 + zp)^2 + 2*xp^2*zp)*
         Log[1 + zp - Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/zp - 
       (16*(1 + 2*xp^2 + 2*xp*(-1 + zp) - zp)*zp*
         Log[1 - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/(1 + zp)) + 
     (16*(1 + 2*xp^2 + 2*xp*(-1 + zp) - zp)*zp*
       Log[1 - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]]*
       Log[1 - 2*xp - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/(1 + zp) + 
     (8*(-1 - 2*xp*(-1 + zp) + zp + 2*xp^2*zp)*
       Log[1 - 2*xp - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]]^2)/
      (zp*(1 + zp)) + Log[xp]*((4*(-11 - 16*xp*(-1 + zp) + 6*zp + 7*zp^2))/
        zp + 32*Sqrt[xp]*Sqrt[zp]*ArcTan[Sqrt[zp]/Sqrt[xp]] + 
       32*Sqrt[xp]*Sqrt[zp]*ArcTan[Sqrt[xp]*Sqrt[zp]] - 
       (16*(1 + 2*xp^2*zp + zp^2 + 2*xp*(1 + zp^2))*Log[1 + xp])/(-1 + zp) - 
       (4*(-3 + 2*zp + xp*(6 - 4*zp + 4*zp^2))*Log[1 - zp])/zp - 
       (4*(-1 - 2*zp - zp^2 + 10*zp^3 + 8*xp^2*zp^3 + 2*zp^4 + 
          xp*(2 + 4*zp + 6*zp^2 - 4*zp^3 + 8*zp^4))*Log[zp])/
        (zp*(-1 + zp^2)) + (16*(zp + 2*xp^2*zp + zp^3 + 2*xp*(1 + zp^2))*
         Log[1 + zp/xp])/(-1 + zp^2) + 
       (16*(1 + zp^2 + 2*xp^2*zp^2 + 2*xp*(zp + zp^3))*Log[1 + xp*zp])/
        (-1 + zp^2) + (16*(1 + 2*xp^2 + 2*xp*(-1 + zp) - zp)*zp*
         Log[1 + zp - Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/(1 + zp) - 
       (16*(1 + 2*xp^2 + 2*xp*(-1 + zp) - zp)*zp*
         Log[1 - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/(1 + zp) - 
       (16*(1 + 2*xp^2 + 2*xp*(-1 + zp) - zp)*zp*
         Log[1 - 2*xp - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/(1 + zp)) + 
     Log[1 - xp]*(-44 + 50/zp - 6*zp + xp*(64 - 72/zp + 12*zp) - 
       (4*(-3 + 2*zp + xp*(6 - 4*zp + 4*zp^2))*Log[xp])/zp + 
       (8*(-1 + 2*xp)*(2 - 2*zp + zp^2)*Log[1 - zp])/zp + 
       (16 - 32*xp)*Log[zp] + (16*(-1 - 2*xp*(-1 + zp) + zp + 2*xp^2*zp)*
         Log[1 + zp - Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/(zp*(1 + zp)) + 
       (16*(1 + 2*xp^2 + 2*xp*(-1 + zp) - zp)*zp*
         Log[1 - 2*xp - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/(1 + zp)) + 
     Log[1 + zp - Sqrt[1 - 2*zp + 4*xp*zp + zp^2]]*
      ((-16*(1 + 2*xp^2 + 2*xp*(-1 + zp) - zp)*zp*
         Log[1 - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/(1 + zp) - 
       (16*(-(-1 + zp)^2 + 2*xp*(-1 + zp)^2 + 2*xp^2*zp)*
         Log[1 - 2*xp - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/zp)), 
 "\!\(\*SubscriptBox[SuperscriptBox[\(G\), \((2)\)], \
\({1,q\!\(\*OverscriptBox[\(q\), \(_\)]\),MS}\)]\)" -> 
  CA*CF*(xp*(44 - 14/zp - 6*zp) + 8*(-6 + 2/zp + zp) + 
     (4*zeta2*(1 + zp - 3*zp^2 + zp^3 + 2*xp^3*(-1 + zp^2) + 
        2*xp^2*(zp + zp^3) + xp*(-1 - zp - 5*zp^2 + 3*zp^3)))/
      ((-1 + xp)*(1 + xp)*(-1 + zp)*(1 + zp)) + 
     ((-4 + 8*zp + 8*zp^2 + 8*xp*(-1 + zp)*zp^2 - 4*zp^3 + 
        4*xp^2*(-1 + zp + 2*zp^2))*Li[2, 1 - xp])/((-1 + xp)*zp*(1 + zp)) + 
     (16*(-1 + xp^2*(-1 + zp) + zp + xp*(-1 + 2*zp))*Li[2, -xp])/
      ((1 + xp)*(-1 + zp)) + (8*(1 + zp^2)*Li[2, -zp])/((-1 + xp)*(1 + zp)) + 
     (16*xp*zp^2*(-1 + xp*zp)*Li[2, -(zp/xp)])/((-1 + xp)*(1 + xp)*(-1 + zp)*
       (1 + zp)) + (16*xp*zp^2*(-xp + zp)*Li[2, -(xp*zp)])/
      ((-1 + xp)*(1 + xp)*(-1 + zp)*(1 + zp)) + 
     ((4 + 4*(-1 + 4*xp)*zp^2)*Li[2, (2 - 2*xp)^(-1) + zp/(2 - 2*xp) - 
         Sqrt[1 - 2*zp + 4*xp*zp + zp^2]/(2 - 2*xp)])/((-1 + xp)*(1 + zp)) + 
     ((-4 + (4 - 16*xp)*zp^2)*Li[2, (2 - 2*xp)^(-1) - (2*xp)/(2 - 2*xp) - 
         zp/(2 - 2*xp) + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]/(2 - 2*xp)])/
      ((-1 + xp)*(1 + zp)) + ((8 + 8*(-1 + 4*xp)*zp^2)*
       Li[2, 1 - 1/(2*xp) + zp/(2*xp) - Sqrt[1 - 2*zp + 4*xp*zp + zp^2]/
          (2*xp)])/((-1 + xp)*(1 + zp)) + ((2 + (-2 + 8*xp)*zp^2)*Log[2]^2)/
      ((-1 + xp)*(1 + zp)) + ((-2*(1 + xp))/zp + (8 + 4/(-1 + xp))*zp + 
       4*(-3 + (1 - zp)^(-1) + 1/((1 + xp)*(-1 + zp)) + 3/(1 + zp) + 
         3/((-1 + xp)*(1 + zp)) + (zp - xp*zp)^(-1)))*Log[xp]^2 + 
     Dzp[0]*(8 - 8*xp - (4*(1 + xp^2)*zeta2)/(1 + xp) - 
       (8*(1 + xp^2)*Li[2, -xp])/(1 + xp) + (2*(1 + xp^2)*Log[xp]^2)/
        (1 + xp) + Log[xp]*(4*(1 + xp) - (8*(1 + xp^2)*Log[1 + xp])/
          (1 + xp))) - (2*(1 + (1 + 2*xp)*zp^2)*Log[zp]^2)/
      ((-1 + xp)*(1 + zp)) + Dxp[0]*(8*(-1 + zp) + (4*zeta2*(1 + zp^2))/
        (1 + zp) + (8*(1 + zp^2)*Li[2, -zp])/(1 + zp) - 
       (2*(1 + zp^2)*Log[zp]^2)/(1 + zp) + 
       Log[zp]*(-4*(1 + zp) + (8*(1 + zp^2)*Log[1 + zp])/(1 + zp))) + 
     Log[Q2/muF2]*(Delz*(8 - 8*xp - (4*(1 + xp^2)*zeta2)/(1 + xp) - 
         (8*(1 + xp^2)*Li[2, -xp])/(1 + xp) + (2*(1 + xp^2)*Log[xp]^2)/
          (1 + xp) + Log[xp]*(4*(1 + xp) - (8*(1 + xp^2)*Log[1 + xp])/
            (1 + xp))) + Delx*(8*(-1 + zp) + (4*zeta2*(1 + zp^2))/(1 + zp) + 
         (8*(1 + zp^2)*Li[2, -zp])/(1 + zp) - (2*(1 + zp^2)*Log[zp]^2)/
          (1 + zp) + Log[zp]*(-4*(1 + zp) + (8*(1 + zp^2)*Log[1 + zp])/
            (1 + zp)))) + Log[zp]*
      ((4*(-1 + xp*(4 - 6*zp) + 3*zp - zp^2 + xp^2*(-2 + 3*zp)))/
        ((-1 + xp)*(-1 + zp)) + (8*(1 + zp^2)*Log[1 + zp])/
        ((-1 + xp)*(1 + zp)) + (16*xp*zp^2*(-1 + xp*zp)*Log[1 + zp/xp])/
        ((-1 + xp)*(1 + xp)*(-1 + zp)*(1 + zp)) + 
       (16*xp*zp^2*(-xp + zp)*Log[1 + xp*zp])/((-1 + xp)*(1 + xp)*(-1 + zp)*
         (1 + zp))) + ((2 + (-2 + 8*xp)*zp^2)*
       Log[1 + zp - Sqrt[1 - 2*zp + 4*xp*zp + zp^2]]^2)/
      ((-1 + xp)*(1 + zp)) + ((2 + (-2 + 8*xp)*zp^2)*
       Log[1 - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]]^2)/
      ((-1 + xp)*(1 + zp)) + Log[1 + zp - Sqrt[1 - 2*zp + 4*xp*zp + zp^2]]*
      ((4*Sqrt[1 + (-2 + 4*xp)*zp + zp^2])/(-1 + xp) + 
       ((-4 + (4 - 16*xp)*zp^2)*Log[1 - zp + Sqrt[1 - 2*zp + 4*xp*zp + 
             zp^2]])/((-1 + xp)*(1 + zp))) + 
     Log[2]*((4*Sqrt[1 + (-2 + 4*xp)*zp + zp^2])/(-1 + xp) + 
       ((4 + 4*(-1 + 4*xp)*zp^2)*Log[xp])/((-1 + xp)*(1 + zp)) + 
       ((-4 + (4 - 16*xp)*zp^2)*Log[1 - zp + Sqrt[1 - 2*zp + 4*xp*zp + 
             zp^2]])/((-1 + xp)*(1 + zp))) - 
     (4*Sqrt[1 + (-2 + 4*xp)*zp + zp^2]*Log[1 - 2*xp - zp + 
         Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/(-1 + xp) + 
     ((-2 + (2 - 8*xp)*zp^2)*Log[1 - 2*xp - zp + Sqrt[1 - 2*zp + 4*xp*zp + 
            zp^2]]^2)/((-1 + xp)*(1 + zp)) + 
     Log[xp]*((2*(5 + xp^2*(-2 + 8*zp) - (2*zp^3)/Sqrt[1 - 2*zp + 4*xp*zp + 
             zp^2] - (8*xp*zp^2)/Sqrt[1 + (-2 + 4*xp)*zp + zp^2] + 
          zp*(-12 - 2/Sqrt[1 - 2*zp + 4*xp*zp + zp^2]) + 
          zp^2*(1 + 4/Sqrt[1 - 2*zp + 4*xp*zp + zp^2])))/(zp - xp*zp) + 
       (16*(-1 + xp^2*(-1 + zp) + zp + xp*(-1 + 2*zp))*Log[1 + xp])/
        ((1 + xp)*(-1 + zp)) - (4*(-(zp^2*(1 + zp)) + xp*zp^2*(1 + 5*zp) + 
          xp^3*(-1 + zp + 2*zp^2) + xp^2*(-1 + zp - 4*zp^2 - 2*zp^3))*
         Log[zp])/((-1 + xp)*(1 + xp)*(-1 + zp)*(1 + zp)) - 
       (16*xp*zp^2*(-1 + xp*zp)*Log[1 + zp/xp])/((-1 + xp)*(1 + xp)*(-1 + zp)*
         (1 + zp)) + (16*xp*zp^2*(-xp + zp)*Log[1 + xp*zp])/
        ((-1 + xp)*(1 + xp)*(-1 + zp)*(1 + zp)) + 
       ((4 + 4*(-1 + 4*xp)*zp^2)*Log[1 + zp - Sqrt[1 - 2*zp + 4*xp*zp + 
             zp^2]])/((-1 + xp)*(1 + zp)) + 
       ((-4 + (4 - 16*xp)*zp^2)*Log[1 - zp + Sqrt[1 - 2*zp + 4*xp*zp + 
             zp^2]])/((-1 + xp)*(1 + zp)) + 
       ((-4 + (4 - 16*xp)*zp^2)*Log[1 - 2*xp - zp + 
           Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/((-1 + xp)*(1 + zp))) + 
     Log[1 - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]]*
      ((-4*Sqrt[1 + (-2 + 4*xp)*zp + zp^2])/(-1 + xp) + 
       ((4 + 4*(-1 + 4*xp)*zp^2)*Log[1 - 2*xp - zp + 
           Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/((-1 + xp)*(1 + zp))) + 
     Log[1 - xp]*(((-4 + (4 - 16*xp)*zp^2)*
         Log[1 + zp - Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/
        ((-1 + xp)*(1 + zp)) + ((4 + 4*(-1 + 4*xp)*zp^2)*
         Log[1 - 2*xp - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/
        ((-1 + xp)*(1 + zp))) + Delz*(13*(-1 + xp) + 6*(1 + xp)*zeta2 + 
       4*(1 + xp)*Li[2, 1 - xp] + 12*(1 + xp)*Li[2, -xp] - 
       (8*(1 + xp^2)*Li[3, 1 - xp])/(1 + xp) + 
       (8*(1 + xp^2)*Li[3, 1/2 - xp/2])/(1 + xp) + 
       (8*(1 + xp^2)*Li[3, 1/2 + xp/2])/(1 + xp) - 
       (8*(1 + xp^2)*Li[3, (-2*xp)/(1 - xp)])/(1 + xp) - 
       (8*(1 + xp^2)*Li[3, (2*xp)/(1 + xp)])/(1 + xp) - 
       (8*(1 + xp^2)*Log[2]^3)/(3*(1 + xp)) + (4*(1 + xp^2)*Log[1 - xp]^3)/
        (3*(1 + xp)) - (4*(1 + xp^2)*Log[1 - xp]^2*Log[xp])/(1 + xp) - 
       ((1 + xp^2)*Log[xp]^3)/(1 + xp) + Log[1 - xp]*
        (8 - 8*xp + (4*(1 + xp^2)*zeta2)/(1 + xp) - (8*(1 + xp^2)*Li[2, -xp])/
          (1 + xp) + 4*(1 + xp)*Log[xp] + (2*(1 + xp^2)*Log[xp]^2)/
          (1 + xp)) - (16*(1 + xp^2)*zeta2*Log[1 + xp])/(1 + xp) + 
       (4*(1 + xp^2)*Log[1 + xp]^3)/(3*(1 + xp)) + 
       Log[2]^2*((4*(1 + xp^2)*Log[1 - xp])/(1 + xp) + 
         (4*(1 + xp^2)*Log[1 + xp])/(1 + xp)) + 
       Log[xp]^2*(-4 - 8*xp + (8*(1 + xp^2)*Log[1 + xp])/(1 + xp)) + 
       Log[2]*((8*(1 + xp^2)*zeta2)/(1 + xp) - (4*(1 + xp^2)*Log[1 - xp]^2)/
          (1 + xp) - (4*(1 + xp^2)*Log[1 + xp]^2)/(1 + xp)) + 
       Log[xp]*(-9 + 7*xp + (4*(1 + xp^2)*zeta2)/(1 + xp) + 
         (4*(1 + xp^2)*Li[2, 1 - xp])/(1 + xp) + (8*(1 + xp^2)*Li[2, -xp])/
          (1 + xp) + 12*(1 + xp)*Log[1 + xp] - (4*(1 + xp^2)*Log[1 + xp]^2)/
          (1 + xp)) + (4*(1 + xp^2)*S12[1 - xp])/(1 + xp)) + 
     Delx*(1 - zp + 2*zeta2*(1 + zp) + (4*zeta3*(1 + zp^2))/(1 + zp) - 
       4*(1 + zp)*Li[2, 1 - zp] + 4*(1 + zp)*Li[2, -zp] + 
       (8*(1 + zp^2)*Li[3, 1 - zp])/(1 + zp) - 
       (8*(1 + zp^2)*Li[3, 1/2 - zp/2])/(1 + zp) - 
       (8*(1 + zp^2)*Li[3, 1/2 + zp/2])/(1 + zp) + (4*(1 + zp^2)*Li[3, -zp])/
        (1 + zp) + (8*(1 + zp^2)*Li[3, (-2*zp)/(1 - zp)])/(1 + zp) + 
       (8*(1 + zp^2)*Li[3, (2*zp)/(1 + zp)])/(1 + zp) + 
       (8*(1 + zp^2)*Log[2]^3)/(3*(1 + zp)) - (4*(1 + zp^2)*Log[1 - zp]^3)/
        (3*(1 + zp)) + (4*(1 + zp^2)*Log[1 - zp]^2*Log[zp])/(1 + zp) - 
       (5*(1 + zp^2)*Log[zp]^3)/(3*(1 + zp)) + Log[1 - zp]*
        (8*(-1 + zp) - (4*zeta2*(1 + zp^2))/(1 + zp) + 
         (8*(1 + zp^2)*Li[2, -zp])/(1 + zp) - 4*(1 + zp)*Log[zp] - 
         (2*(1 + zp^2)*Log[zp]^2)/(1 + zp)) + 
       ((12*zeta2*(1 + zp^2))/(1 + zp) - (8*(1 + zp^2)*Li[2, -zp])/(1 + zp))*
        Log[1 + zp] - (4*(1 + zp^2)*Log[1 + zp]^3)/(3*(1 + zp)) + 
       Log[zp]*(-7 + zp + (4*zeta2*(1 + zp^2))/(1 + zp) - 
         (4*(1 + zp^2)*Li[2, 1 - zp])/(1 + zp) + (4*(1 + zp^2)*Li[2, -zp])/
          (1 + zp) + 4*(1 + zp)*Log[1 + zp]) + 
       Log[2]^2*((-4*(1 + zp^2)*Log[1 - zp])/(1 + zp) - 
         (4*(1 + zp^2)*Log[1 + zp])/(1 + zp)) + 
       Log[zp]^2*(-4*(1 + zp) + (6*(1 + zp^2)*Log[1 + zp])/(1 + zp)) + 
       Log[2]*((-8*zeta2*(1 + zp^2))/(1 + zp) + (4*(1 + zp^2)*Log[1 - zp]^2)/
          (1 + zp) + (4*(1 + zp^2)*Log[1 + zp]^2)/(1 + zp)) - 
       (4*(1 + zp^2)*S12[1 - zp])/(1 + zp) - (8*(1 + zp^2)*S12[-zp])/
        (1 + zp))) + 
   CF^2*((4*(-4*(2 - 6*zp + zp^2) + xp*(7 - 22*zp + 3*zp^2)))/zp - 
     (8*zeta2*(1 + zp - 3*zp^2 + zp^3 + 2*xp^3*(-1 + zp^2) + 
        2*xp^2*(zp + zp^3) + xp*(-1 - zp - 5*zp^2 + 3*zp^3)))/
      ((-1 + xp)*(1 + xp)*(-1 + zp)*(1 + zp)) - 
     (8*(-1 + 2*zp + 2*zp^2 + 2*xp*(-1 + zp)*zp^2 - zp^3 + 
        xp^2*(-1 + zp + 2*zp^2))*Li[2, 1 - xp])/((-1 + xp)*zp*(1 + zp)) - 
     (32*(-1 + xp^2*(-1 + zp) + zp + xp*(-1 + 2*zp))*Li[2, -xp])/
      ((1 + xp)*(-1 + zp)) - (16*(1 + zp^2)*Li[2, -zp])/
      ((-1 + xp)*(1 + zp)) - (32*xp*zp^2*(-1 + xp*zp)*Li[2, -(zp/xp)])/
      ((-1 + xp)*(1 + xp)*(-1 + zp)*(1 + zp)) + 
     (32*xp*(xp - zp)*zp^2*Li[2, -(xp*zp)])/((-1 + xp)*(1 + xp)*(-1 + zp)*
       (1 + zp)) + ((-8 + (8 - 32*xp)*zp^2)*
       Li[2, (2 - 2*xp)^(-1) + zp/(2 - 2*xp) - 
         Sqrt[1 - 2*zp + 4*xp*zp + zp^2]/(2 - 2*xp)])/((-1 + xp)*(1 + zp)) + 
     ((8 + 8*(-1 + 4*xp)*zp^2)*Li[2, (2 - 2*xp)^(-1) - (2*xp)/(2 - 2*xp) - 
         zp/(2 - 2*xp) + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]/(2 - 2*xp)])/
      ((-1 + xp)*(1 + zp)) - (16*(1 + (-1 + 4*xp)*zp^2)*
       Li[2, 1 - 1/(2*xp) + zp/(2*xp) - Sqrt[1 - 2*zp + 4*xp*zp + zp^2]/
          (2*xp)])/((-1 + xp)*(1 + zp)) + ((-4 + (4 - 16*xp)*zp^2)*Log[2]^2)/
      ((-1 + xp)*(1 + zp)) + 4*(6 + 2/(-1 + zp) - 2/((1 + xp)*(-1 + zp)) + 
       (1 + xp)/zp + (-4 - 2/(-1 + xp))*zp - 6/(1 + zp) - 
       6/((-1 + xp)*(1 + zp)) - 2/(zp - xp*zp))*Log[xp]^2 + 
     Dzp[0]*(16*(-1 + xp) + (8*(1 + xp^2)*zeta2)/(1 + xp) + 
       (16*(1 + xp^2)*Li[2, -xp])/(1 + xp) - (4*(1 + xp^2)*Log[xp]^2)/
        (1 + xp) + Log[xp]*(-8*(1 + xp) + (16*(1 + xp^2)*Log[1 + xp])/
          (1 + xp))) + ((4 + (4 + 8*xp)*zp^2)*Log[zp]^2)/
      ((-1 + xp)*(1 + zp)) + Dxp[0]*(-16*(-1 + zp) - (8*zeta2*(1 + zp^2))/
        (1 + zp) - (16*(1 + zp^2)*Li[2, -zp])/(1 + zp) + 
       (4*(1 + zp^2)*Log[zp]^2)/(1 + zp) + 
       Log[zp]*(8*(1 + zp) - (16*(1 + zp^2)*Log[1 + zp])/(1 + zp))) + 
     Log[Q2/muF2]*(Delz*(16*(-1 + xp) + (8*(1 + xp^2)*zeta2)/(1 + xp) + 
         (16*(1 + xp^2)*Li[2, -xp])/(1 + xp) - (4*(1 + xp^2)*Log[xp]^2)/
          (1 + xp) + Log[xp]*(-8*(1 + xp) + (16*(1 + xp^2)*Log[1 + xp])/
            (1 + xp))) + Delx*(-16*(-1 + zp) - (8*zeta2*(1 + zp^2))/
          (1 + zp) - (16*(1 + zp^2)*Li[2, -zp])/(1 + zp) + 
         (4*(1 + zp^2)*Log[zp]^2)/(1 + zp) + Log[zp]*(8*(1 + zp) - 
           (16*(1 + zp^2)*Log[1 + zp])/(1 + zp)))) + 
     Log[zp]*((-8*(-1 + xp*(4 - 6*zp) + 3*zp - zp^2 + xp^2*(-2 + 3*zp)))/
        ((-1 + xp)*(-1 + zp)) - (16*(1 + zp^2)*Log[1 + zp])/
        ((-1 + xp)*(1 + zp)) - (32*xp*zp^2*(-1 + xp*zp)*Log[1 + zp/xp])/
        ((-1 + xp)*(1 + xp)*(-1 + zp)*(1 + zp)) + 
       (32*xp*(xp - zp)*zp^2*Log[1 + xp*zp])/((-1 + xp)*(1 + xp)*(-1 + zp)*
         (1 + zp))) + ((-4 + (4 - 16*xp)*zp^2)*
       Log[1 + zp - Sqrt[1 - 2*zp + 4*xp*zp + zp^2]]^2)/
      ((-1 + xp)*(1 + zp)) + ((-4 + (4 - 16*xp)*zp^2)*
       Log[1 - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]]^2)/
      ((-1 + xp)*(1 + zp)) + Log[1 + zp - Sqrt[1 - 2*zp + 4*xp*zp + zp^2]]*
      ((-8*Sqrt[1 + (-2 + 4*xp)*zp + zp^2])/(-1 + xp) + 
       ((8 + 8*(-1 + 4*xp)*zp^2)*Log[1 - zp + Sqrt[1 - 2*zp + 4*xp*zp + 
             zp^2]])/((-1 + xp)*(1 + zp))) + 
     Log[2]*((-8*Sqrt[1 + (-2 + 4*xp)*zp + zp^2])/(-1 + xp) + 
       ((-8 + (8 - 32*xp)*zp^2)*Log[xp])/((-1 + xp)*(1 + zp)) + 
       ((8 + 8*(-1 + 4*xp)*zp^2)*Log[1 - zp + Sqrt[1 - 2*zp + 4*xp*zp + 
             zp^2]])/((-1 + xp)*(1 + zp))) + 
     (8*Sqrt[1 + (-2 + 4*xp)*zp + zp^2]*Log[1 - 2*xp - zp + 
         Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/(-1 + xp) + 
     ((4 + 4*(-1 + 4*xp)*zp^2)*
       Log[1 - 2*xp - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]]^2)/
      ((-1 + xp)*(1 + zp)) + Log[1 - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]]*
      ((8*Sqrt[1 + (-2 + 4*xp)*zp + zp^2])/(-1 + xp) + 
       ((-8 + (8 - 32*xp)*zp^2)*Log[1 - 2*xp - zp + 
           Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/((-1 + xp)*(1 + zp))) + 
     Log[1 - xp]*(((8 + 8*(-1 + 4*xp)*zp^2)*
         Log[1 + zp - Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/
        ((-1 + xp)*(1 + zp)) + ((-8 + (8 - 32*xp)*zp^2)*
         Log[1 - 2*xp - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/
        ((-1 + xp)*(1 + zp))) + 
     Log[xp]*((-4*(5 + xp^2*(-2 + 8*zp) - (2*zp^3)/
           Sqrt[1 - 2*zp + 4*xp*zp + zp^2] - (8*xp*zp^2)/
           Sqrt[1 + (-2 + 4*xp)*zp + zp^2] + 
          zp*(-12 - 2/Sqrt[1 - 2*zp + 4*xp*zp + zp^2]) + 
          zp^2*(1 + 4/Sqrt[1 - 2*zp + 4*xp*zp + zp^2])))/(zp - xp*zp) - 
       (32*(-1 + xp^2*(-1 + zp) + zp + xp*(-1 + 2*zp))*Log[1 + xp])/
        ((1 + xp)*(-1 + zp)) + (8*(-(zp^2*(1 + zp)) + xp*zp^2*(1 + 5*zp) + 
          xp^3*(-1 + zp + 2*zp^2) + xp^2*(-1 + zp - 4*zp^2 - 2*zp^3))*
         Log[zp])/((-1 + xp)*(1 + xp)*(-1 + zp)*(1 + zp)) + 
       (32*xp*zp^2*(-1 + xp*zp)*Log[1 + zp/xp])/((-1 + xp)*(1 + xp)*(-1 + zp)*
         (1 + zp)) + (32*xp*(xp - zp)*zp^2*Log[1 + xp*zp])/
        ((-1 + xp)*(1 + xp)*(-1 + zp)*(1 + zp)) + 
       ((-8 + (8 - 32*xp)*zp^2)*Log[1 + zp - Sqrt[1 - 2*zp + 4*xp*zp + 
             zp^2]])/((-1 + xp)*(1 + zp)) + ((8 + 8*(-1 + 4*xp)*zp^2)*
         Log[1 - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/
        ((-1 + xp)*(1 + zp)) + ((8 + 8*(-1 + 4*xp)*zp^2)*
         Log[1 - 2*xp - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/
        ((-1 + xp)*(1 + zp))) + Delz*(-26*(-1 + xp) - 12*(1 + xp)*zeta2 - 
       8*(1 + xp)*Li[2, 1 - xp] - 24*(1 + xp)*Li[2, -xp] + 
       (16*(1 + xp^2)*Li[3, 1 - xp])/(1 + xp) - 
       (16*(1 + xp^2)*Li[3, 1/2 - xp/2])/(1 + xp) - 
       (16*(1 + xp^2)*Li[3, 1/2 + xp/2])/(1 + xp) + 
       (16*(1 + xp^2)*Li[3, (-2*xp)/(1 - xp)])/(1 + xp) + 
       (16*(1 + xp^2)*Li[3, (2*xp)/(1 + xp)])/(1 + xp) + 
       (16*(1 + xp^2)*Log[2]^3)/(3*(1 + xp)) - (8*(1 + xp^2)*Log[1 - xp]^3)/
        (3*(1 + xp)) + (8*(1 + xp^2)*Log[1 - xp]^2*Log[xp])/(1 + xp) + 
       (2*(1 + xp^2)*Log[xp]^3)/(1 + xp) + Log[1 - xp]*
        (16*(-1 + xp) - (8*(1 + xp^2)*zeta2)/(1 + xp) + 
         (16*(1 + xp^2)*Li[2, -xp])/(1 + xp) - 8*(1 + xp)*Log[xp] - 
         (4*(1 + xp^2)*Log[xp]^2)/(1 + xp)) + 
       (32*(1 + xp^2)*zeta2*Log[1 + xp])/(1 + xp) - 
       (8*(1 + xp^2)*Log[1 + xp]^3)/(3*(1 + xp)) + 
       Log[xp]^2*(8 + 16*xp - (16*(1 + xp^2)*Log[1 + xp])/(1 + xp)) + 
       Log[2]^2*((-8*(1 + xp^2)*Log[1 - xp])/(1 + xp) - 
         (8*(1 + xp^2)*Log[1 + xp])/(1 + xp)) + 
       Log[2]*((-16*(1 + xp^2)*zeta2)/(1 + xp) + (8*(1 + xp^2)*Log[1 - xp]^2)/
          (1 + xp) + (8*(1 + xp^2)*Log[1 + xp]^2)/(1 + xp)) + 
       Log[xp]*(18 - 14*xp - (8*(1 + xp^2)*zeta2)/(1 + xp) - 
         (8*(1 + xp^2)*Li[2, 1 - xp])/(1 + xp) - (16*(1 + xp^2)*Li[2, -xp])/
          (1 + xp) - 24*(1 + xp)*Log[1 + xp] + (8*(1 + xp^2)*Log[1 + xp]^2)/
          (1 + xp)) - (8*(1 + xp^2)*S12[1 - xp])/(1 + xp)) + 
     Delx*(2*(-1 + zp) - 4*zeta2*(1 + zp) - (8*zeta3*(1 + zp^2))/(1 + zp) + 
       8*(1 + zp)*Li[2, 1 - zp] - 8*(1 + zp)*Li[2, -zp] - 
       (16*(1 + zp^2)*Li[3, 1 - zp])/(1 + zp) + 
       (16*(1 + zp^2)*Li[3, 1/2 - zp/2])/(1 + zp) + 
       (16*(1 + zp^2)*Li[3, 1/2 + zp/2])/(1 + zp) - (8*(1 + zp^2)*Li[3, -zp])/
        (1 + zp) - (16*(1 + zp^2)*Li[3, (-2*zp)/(1 - zp)])/(1 + zp) - 
       (16*(1 + zp^2)*Li[3, (2*zp)/(1 + zp)])/(1 + zp) - 
       (16*(1 + zp^2)*Log[2]^3)/(3*(1 + zp)) + (8*(1 + zp^2)*Log[1 - zp]^3)/
        (3*(1 + zp)) - (8*(1 + zp^2)*Log[1 - zp]^2*Log[zp])/(1 + zp) + 
       (10*(1 + zp^2)*Log[zp]^3)/(3*(1 + zp)) + Log[1 - zp]*
        (-16*(-1 + zp) + (8*zeta2*(1 + zp^2))/(1 + zp) - 
         (16*(1 + zp^2)*Li[2, -zp])/(1 + zp) + 8*(1 + zp)*Log[zp] + 
         (4*(1 + zp^2)*Log[zp]^2)/(1 + zp)) + 
       ((-24*zeta2*(1 + zp^2))/(1 + zp) + (16*(1 + zp^2)*Li[2, -zp])/
          (1 + zp))*Log[1 + zp] + (8*(1 + zp^2)*Log[1 + zp]^3)/(3*(1 + zp)) + 
       Log[zp]*(-2*(-7 + zp) - (8*zeta2*(1 + zp^2))/(1 + zp) + 
         (8*(1 + zp^2)*Li[2, 1 - zp])/(1 + zp) - (8*(1 + zp^2)*Li[2, -zp])/
          (1 + zp) - 8*(1 + zp)*Log[1 + zp]) + 
       Log[zp]^2*(8*(1 + zp) - (12*(1 + zp^2)*Log[1 + zp])/(1 + zp)) + 
       Log[2]^2*((8*(1 + zp^2)*Log[1 - zp])/(1 + zp) + 
         (8*(1 + zp^2)*Log[1 + zp])/(1 + zp)) + 
       Log[2]*((16*zeta2*(1 + zp^2))/(1 + zp) - (8*(1 + zp^2)*Log[1 - zp]^2)/
          (1 + zp) - (8*(1 + zp^2)*Log[1 + zp]^2)/(1 + zp)) + 
       (8*(1 + zp^2)*S12[1 - zp])/(1 + zp) + (16*(1 + zp^2)*S12[-zp])/
        (1 + zp))) + 
   CF*((8*zeta2*(-1 + xp^2 + 4*zp + xp*(2 - 4*zp + 4*zp^2)))/(-1 + xp) + 
     (-113 + 204*zp + 255*zp^2 - 130*zp^3 + xp*(127 - 300*zp - 141*zp^2 + 
         62*zp^3))/(9*zp) - 16*Sqrt[xp]*Sqrt[zp]*
      InvTanInt[-(Sqrt[zp]/Sqrt[xp])] + 16*Sqrt[xp]*Sqrt[zp]*
      InvTanInt[Sqrt[zp]/Sqrt[xp]] + 16*Sqrt[xp]*Sqrt[zp]*
      InvTanInt[-(Sqrt[xp]*Sqrt[zp])] - 16*Sqrt[xp]*Sqrt[zp]*
      InvTanInt[Sqrt[xp]*Sqrt[zp]] + 
     (4*(-1 + 2*zp + 4*zp^2 + xp^2*(1 - 2*zp + 4*zp^2))*Li[2, 1 - xp])/
      ((-1 + xp)*zp) - 16*(1 + xp)*(-1 + 2*zp)*Li[2, -xp] + 
     ((8 - 4*xp^2 + 8*xp*(-1 + zp))*Li[2, 1 - zp])/(-1 + xp) + 
     (8*(xp + xp^3 + 2*zp + 2*xp^2*zp + 4*xp*zp^2)*Li[2, -(zp/xp)])/
      (-1 + xp^2) + (8*(1 + 2*xp*zp + 2*xp^3*zp + xp^2*(1 + 4*zp^2))*
       Li[2, -(xp*zp)])/(-1 + xp^2) - 
     (8*(1 + xp^4 + xp*(3 - 6*zp) + xp^3*(3 - 6*zp) + 
        4*xp^2*(1 - 3*zp + 3*zp^2))*Li[2, 1/2 - xp/2 - 
         Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]/2])/(1 + xp^2 + xp*(2 - 4*zp))^
       (3/2) - (8*(1 + xp^4 + xp*(3 - 6*zp) + xp^3*(3 - 6*zp) + 
        4*xp^2*(1 - 3*zp + 3*zp^2))*Li[2, 1/2 + xp/2 - 
         Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]/2])/(1 + xp^2 + xp*(2 - 4*zp))^
       (3/2) + (8*(1 + xp^4 + xp*(3 - 6*zp) + xp^3*(3 - 6*zp) + 
        4*xp^2*(1 - 3*zp + 3*zp^2))*Li[2, 1/2 - 1/(2*xp) - 
         Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]/(2*xp)])/(1 + xp^2 + xp*(2 - 4*zp))^
       (3/2) + (8*(1 + xp^4 + xp*(3 - 6*zp) + xp^3*(3 - 6*zp) + 
        4*xp^2*(1 - 3*zp + 3*zp^2))*Li[2, 1/2 + 1/(2*xp) - 
         Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]/(2*xp)])/(1 + xp^2 + xp*(2 - 4*zp))^
       (3/2) + (4*(1 + 2*xp^2 + 2*zp + xp*(-2 + 4*zp))*
       Li[2, (2 - 2*xp)^(-1) + zp/(2 - 2*xp) - 
         Sqrt[1 - 2*zp + 4*xp*zp + zp^2]/(2 - 2*xp)])/(-1 + xp) - 
     (4*(1 + 2*zp + 4*xp^2*zp + xp*(2 - 4*zp + 8*zp^2))*
       Li[2, (2 - 2*xp)^(-1) - (2*xp)/(2 - 2*xp) - zp/(2 - 2*xp) + 
         Sqrt[1 - 2*zp + 4*xp*zp + zp^2]/(2 - 2*xp)])/(-1 + xp) + 
     (8*(1 + 2*zp + 4*xp*zp^2 + xp^2*(1 + 2*zp))*
       Li[2, 1 - 1/(2*xp) + zp/(2*xp) - Sqrt[1 - 2*zp + 4*xp*zp + zp^2]/
          (2*xp)])/(-1 + xp) + ((2 + 4*xp^2 + 4*zp + xp*(-4 + 8*zp))*
       Log[2]^2)/(-1 + xp) - (4*xp*(2 - 4*zp + 4*zp^2 + xp*(-1 + 2*zp))*
       Log[1 - xp]^2)/(-1 + xp) + 2*(2 + (1 - xp)^(-1) + 
       (1 + xp^2 + xp*(2 - 4*zp))^(-3/2) + xp^4/(1 + xp^2 + xp*(2 - 4*zp))^
         (3/2) + (xp^2*(1 + 3*xp^2 + xp*(6 - 12*zp)))/
        (1 + xp^2 + xp*(2 - 4*zp))^(3/2) + 3/Sqrt[1 + xp^2 + xp*(2 - 4*zp)] + 
       (6*xp)/Sqrt[1 + xp^2 + xp*(2 - 4*zp)] - (2*(1 + xp))/zp + 
       2*(4 + 5/(-1 + xp) + xp*(6 - 6/Sqrt[1 + xp^2 + xp*(2 - 4*zp)]))*zp)*
      Log[xp]^2 + Dzp[1]*(-10*(-1 + xp) + 4*(1 + xp)*Log[xp]) + 
     Dzp[0]*(12*(-1 + xp) + 4*(1 + xp)*Li[2, 1 - xp] + 12*(-1 + xp)*Log[xp] - 
       4*(1 + xp)*Log[xp]^2 + Log[1 - xp]*(-10*(-1 + xp) + 
         4*(1 + xp)*Log[xp])) + Delz*(-36*(-1 + xp) + 10*(-1 + xp)*zeta2 + 
       2*(-1 + xp)*Li[2, 1 - xp] - 4*(1 + xp)*Li[3, 1 - xp] + 
       (25 - 3*xp - 4*(1 + xp)*zeta2 - 4*(1 + xp)*Li[2, 1 - xp])*Log[xp] + 
       ((17 - 15*xp)*Log[xp]^2)/2 + (5*(1 + xp)*Log[xp]^3)/3 + 
       Log[1 - xp]^2*(5 - 5*xp + 2*(1 + xp)*Log[xp]) + 
       Log[1 - xp]*(12*(-1 + xp) + 4*(1 + xp)*Li[2, 1 - xp] + 
         12*(-1 + xp)*Log[xp] - 4*(1 + xp)*Log[xp]^2)) + 
     ((2 + 4*zp - 4*xp^2*(1 + 2*zp) - 4*xp*(1 + 4*zp^2))*Log[zp]^2)/
      (-1 + xp) + Dxp[1]*(2 + 8/(3*zp) - 2*zp - (8*zp^2)/3 + 
       4*(1 + zp)*Log[zp]) + Log[1 - zp]*
      ((26 - 78*zp + 6*zp^2 + 16*zp^3 + xp*(-34 + 54*zp + 18*zp^2 - 8*zp^3))/
        (3*zp) - 4*(1 + xp + 2*zp)*Log[zp]) + Log[Q2/muF2]^2*
      (Delz*(5 - 5*xp + 2*(1 + xp)*Log[xp]) + 
       Delx*(1 + 4/(3*zp) - zp - (4*zp^2)/3 + 2*(1 + zp)*Log[zp])) + 
     Dxp[0]*((2*(-7 - 60*zp + 42*zp^2 + 25*zp^3))/(9*zp) + 
       4*(1 + zp)*Li[2, 1 - zp] + (-2 + 8/(3*zp) - 8*zp - (8*zp^2)/3)*
        Log[zp] + 4*(1 + zp)*Log[zp]^2 + Log[1 - zp]*(2 + 8/(3*zp) - 2*zp - 
         (8*zp^2)/3 + 4*(1 + zp)*Log[zp])) + Log[Q2/muF2]*
      ((26 - 78*zp + 6*zp^2 + 16*zp^3 + xp*(-34 + 54*zp + 18*zp^2 - 8*zp^3))/
        (3*zp) - (4*(1 + xp)*(-1 + 2*zp)*Log[xp])/zp + 
       Dzp[0]*(-10*(-1 + xp) + 4*(1 + xp)*Log[xp]) + 
       Delz*(12*(-1 + xp) + 4*(1 + xp)*Li[2, 1 - xp] + 12*(-1 + xp)*Log[xp] - 
         4*(1 + xp)*Log[xp]^2 + Log[1 - xp]*(-10*(-1 + xp) + 
           4*(1 + xp)*Log[xp])) - 4*(1 + xp + 2*zp)*Log[zp] + 
       Dxp[0]*(2 + 8/(3*zp) - 2*zp - (8*zp^2)/3 + 4*(1 + zp)*Log[zp]) + 
       Delx*((2*(-7 - 60*zp + 42*zp^2 + 25*zp^3))/(9*zp) + 
         4*(1 + zp)*Li[2, 1 - zp] + (-2 + 8/(3*zp) - 8*zp - (8*zp^2)/3)*
          Log[zp] + 4*(1 + zp)*Log[zp]^2 + Log[1 - zp]*(2 + 8/(3*zp) - 2*zp - 
           (8*zp^2)/3 + 4*(1 + zp)*Log[zp]))) + 
     Log[zp]*((2*(6*(-9 - 3/(-1 + xp) + (1 + xp^2 + xp*(2 - 4*zp))^(-1) + 
            (1 - zp)^(-1)) + 13/zp - (3*(5 + xp)*zp)/(-1 + xp) + 8*zp^2 - 
          (6*xp^3)/(1 + 2*xp + xp^2 - 4*xp*zp) + 
          xp^2*(12/(xp - zp) - 6/(1 + 2*xp + xp^2 - 4*xp*zp)) + 
          xp*(-17/zp + 9*zp - 4*zp^2 + 6*(1 + (-1 + zp)^(-1) + 
              (-xp + zp)^(-1) + (1 + 2*xp + xp^2 - 4*xp*zp)^(-1)))))/3 - 
       16*Sqrt[xp]*Sqrt[zp]*ArcTan[Sqrt[zp]/Sqrt[xp]] + 
       16*Sqrt[xp]*Sqrt[zp]*ArcTan[Sqrt[xp]*Sqrt[zp]] + 
       (8*(xp + xp^3 + 2*zp + 2*xp^2*zp + 4*xp*zp^2)*Log[1 + zp/xp])/
        (-1 + xp^2) + (8*(1 + 2*xp*zp + 2*xp^3*zp + xp^2*(1 + 4*zp^2))*
         Log[1 + xp*zp])/(-1 + xp^2)) - 
     (8*(1 + xp^4 + xp*(3 - 6*zp) + xp^3*(3 - 6*zp) + 
        4*xp^2*(1 - 3*zp + 3*zp^2))*
       Log[1 + xp - Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]]*
       Log[1 - xp + Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]])/
      (1 + xp^2 + xp*(2 - 4*zp))^(3/2) + 
     (8*(1 + xp^4 + xp*(3 - 6*zp) + xp^3*(3 - 6*zp) + 
        4*xp^2*(1 - 3*zp + 3*zp^2))*
       Log[1 + xp - Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]]*
       Log[-1 + xp + Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]])/
      (1 + xp^2 + xp*(2 - 4*zp))^(3/2) + 
     ((2 + 4*zp + 8*xp^2*zp + 4*xp*(1 - 2*zp + 4*zp^2))*
       Log[1 + zp - Sqrt[1 - 2*zp + 4*xp*zp + zp^2]]^2)/(-1 + xp) + 
     ((2 + 4*zp + 8*xp^2*zp + 4*xp*(1 - 2*zp + 4*zp^2))*
       Log[1 - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]]^2)/(-1 + xp) + 
     Log[2]*((4*(1 + 2*xp)*Sqrt[1 + (-2 + 4*xp)*zp + zp^2])/(-1 + xp) - 
       (8*xp*(2 - 4*zp + 4*zp^2 + xp*(-1 + 2*zp))*Log[1 - xp])/(-1 + xp) + 
       4*(2 + 3/(-1 + xp) + (1 + xp^2 + xp*(2 - 4*zp))^(-3/2) + 
         xp^4/(1 + xp^2 + xp*(2 - 4*zp))^(3/2) + 
         (xp^2*(1 + 3*xp^2 + xp*(6 - 12*zp)))/(1 + xp^2 + xp*(2 - 4*zp))^
           (3/2) + 3/Sqrt[1 + xp^2 + xp*(2 - 4*zp)] + 
         (6*xp)/Sqrt[1 + xp^2 + xp*(2 - 4*zp)] + 
         2*((-1 + xp)^(-1) + xp*(2 - 6/Sqrt[1 + xp^2 + xp*(2 - 4*zp)]))*zp + 
         (8*xp*zp^2)/(-1 + xp))*Log[xp] + 
       (8*(1 + xp^4 + xp*(3 - 6*zp) + xp^3*(3 - 6*zp) + 
          4*xp^2*(1 - 3*zp + 3*zp^2))*
         Log[1 - xp + Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]])/
        (1 + xp^2 + xp*(2 - 4*zp))^(3/2) - 
       (8*(1 + xp^4 + xp*(3 - 6*zp) + xp^3*(3 - 6*zp) + 
          4*xp^2*(1 - 3*zp + 3*zp^2))*Log[-1 + xp + 
           Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]])/(1 + xp^2 + xp*(2 - 4*zp))^
         (3/2) + (8*xp*(2 - 4*zp + 4*zp^2 + xp*(-1 + 2*zp))*
         Log[1 + zp - Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/(-1 + xp) - 
       (4*(1 + 2*zp + 4*xp^2*zp + xp*(2 - 4*zp + 8*zp^2))*
         Log[1 - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/(-1 + xp)) - 
     (4*(1 + 2*xp)*Sqrt[1 + (-2 + 4*xp)*zp + zp^2]*
       Log[1 - 2*xp - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/(-1 + xp) + 
     ((2 - 4*xp + 4*xp^2 + 4*zp + 8*xp*zp)*
       Log[1 - 2*xp - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]]^2)/(1 - xp) + 
     Log[1 + zp - Sqrt[1 - 2*zp + 4*xp*zp + zp^2]]*
      ((4*(1 + 2*xp)*Sqrt[1 + (-2 + 4*xp)*zp + zp^2])/(-1 + xp) - 
       (4*(1 + 2*zp + 4*xp^2*zp + xp*(2 - 4*zp + 8*zp^2))*
         Log[1 - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/(-1 + xp) - 
       (8*xp*(2 - 4*zp + 4*zp^2 + xp*(-1 + 2*zp))*
         Log[1 - 2*xp - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/(-1 + xp)) + 
     Log[xp]*((2*(57 + 3/(-1 + xp) - 14/zp + 6/(1 + 2*xp + xp^2 - 4*xp*zp) + 
          (6*xp^3)/(1 + 2*xp + xp^2 - 4*xp*zp) - 8/(zp - xp*zp) + 
          12/Sqrt[1 + (-2 + 4*xp)*zp + zp^2] + 
          18/((-1 + xp)*Sqrt[1 + (-2 + 4*xp)*zp + zp^2]) + 
          xp^2*(12/(-xp + zp) - 6/(1 + 2*xp + xp^2 - 4*xp*zp)) + 
          (2*zp^2*(4 + 3/Sqrt[1 + (-2 + 4*xp)*zp + zp^2]))/(-1 + xp) + 
          (zp*(3 + 12/Sqrt[1 + (-2 + 4*xp)*zp + zp^2]))/(1 - xp) + 
          xp*(-9 + 6/(xp - zp) + 22/zp - 6/(1 + 2*xp + xp^2 - 4*xp*zp) + 
            (4*zp^2*(6 - 2*xp - 3/Sqrt[1 - 2*zp + 4*xp*zp + zp^2]))/
             (1 - xp) + (6*zp*(-4 + xp*(3 - 8/Sqrt[1 - 2*zp + 4*xp*zp + 
                    zp^2])))/(1 - xp))))/3 + 16*Sqrt[xp]*Sqrt[zp]*
        ArcTan[Sqrt[zp]/Sqrt[xp]] + 16*Sqrt[xp]*Sqrt[zp]*
        ArcTan[Sqrt[xp]*Sqrt[zp]] - 16*(1 + xp)*(-1 + 2*zp)*Log[1 + xp] - 
       (4*(1 + xp)*(-1 + 2*zp)*Log[1 - zp])/zp - 
       (4*((-1 + zp)^2 + xp^3*(-1 - zp + 4*zp^2) + 
          xp*(1 - 6*zp + 5*zp^2 - 8*zp^3) + xp^2*(-1 + 3*zp - 8*zp^2 + 
            8*zp^3))*Log[zp])/((-1 + xp)*(1 + xp)*zp) - 
       (8*(xp + xp^3 + 2*zp + 2*xp^2*zp + 4*xp*zp^2)*Log[1 + zp/xp])/
        (-1 + xp^2) + (8*(1 + 2*xp*zp + 2*xp^3*zp + xp^2*(1 + 4*zp^2))*
         Log[1 + xp*zp])/(-1 + xp^2) - 
       (8*(1 + xp^4 + xp*(3 - 6*zp) + xp^3*(3 - 6*zp) + 
          4*xp^2*(1 - 3*zp + 3*zp^2))*Log[-1 + xp + 
           Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]])/(1 + xp^2 + xp*(2 - 4*zp))^
         (3/2) - (8*(1 + xp^4 + xp*(3 - 6*zp) + xp^3*(3 - 6*zp) + 
          4*xp^2*(1 - 3*zp + 3*zp^2))*
         Log[1 + xp + Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]])/
        (1 + xp^2 + xp*(2 - 4*zp))^(3/2) + 
       (4*(1 + 2*zp + 4*xp^2*zp + xp*(2 - 4*zp + 8*zp^2))*
         Log[1 + zp - Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/(-1 + xp) - 
       (4*(1 + 2*zp + 4*xp^2*zp + xp*(2 - 4*zp + 8*zp^2))*
         Log[1 - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/(-1 + xp) - 
       (4*(1 + 2*zp + 4*xp^2*zp + xp*(2 - 4*zp + 8*zp^2))*
         Log[1 - 2*xp - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/(-1 + xp)) + 
     Log[1 - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]]*
      ((-4*(1 + 2*xp)*Sqrt[1 + (-2 + 4*xp)*zp + zp^2])/(-1 + xp) + 
       (4*(1 + 2*zp + 4*xp^2*zp + xp*(2 - 4*zp + 8*zp^2))*
         Log[1 - 2*xp - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/(-1 + xp)) + 
     Log[1 - xp]*((26 - 78*zp + 6*zp^2 + 16*zp^3 + 
         xp*(-34 + 54*zp + 18*zp^2 - 8*zp^3))/(3*zp) - 
       (4*(1 + xp)*(-1 + 2*zp)*Log[xp])/zp - 4*(1 + xp + 2*zp)*Log[zp] - 
       (4*(1 + 2*xp^2 + 2*zp + xp*(-2 + 4*zp))*
         Log[1 + zp - Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/(-1 + xp) + 
       (4*(1 + 2*zp + 4*xp^2*zp + xp*(2 - 4*zp + 8*zp^2))*
         Log[1 - 2*xp - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/(-1 + xp)) + 
     Delx*(zeta2*(-2 - 8/(3*zp) + 2*zp + (8*zp^2)/3) - 
       (107 + 132*zp + 48*zp^2 - 287*zp^3)/(27*zp) + 
       (-2 + 8/(3*zp) - 8*zp - (8*zp^2)/3)*Li[2, 1 - zp] - 
       4*(1 + zp)*Li[3, 1 - zp] + (-4*zeta2*(1 + zp) - 
         (2*(7 + 90*zp + 81*zp^2 + 31*zp^3))/(9*zp) + 
         8*(1 + zp)*Li[2, 1 - zp])*Log[zp] + ((-3 + 8/zp - 15*zp)*Log[zp]^2)/
        6 + (5*(1 + zp)*Log[zp]^3)/3 + Log[1 - zp]^2*(1 + 4/(3*zp) - zp - 
         (4*zp^2)/3 + 2*(1 + zp)*Log[zp]) + Log[1 - zp]*
        ((2*(-7 - 60*zp + 42*zp^2 + 25*zp^3))/(9*zp) + 
         4*(1 + zp)*Li[2, 1 - zp] + (-2 + 8/(3*zp) - 8*zp - (8*zp^2)/3)*
          Log[zp] + 4*(1 + zp)*Log[zp]^2) + 8*(1 + zp)*S12[1 - zp]) + 
     ((4*(-1 + 2*xp)*(-1 + 2*zp)*Li[2, zp/xp])/(-1 + xp) + 
       (4*(-1 + 2*xp)*(-1 + 2*zp)*Log[xp]^2)/(-1 + xp) + 
       (4*(-1 + 2*xp)*(-1 + 2*zp)*Log[xp - zp]*Log[zp])/(-1 + xp) + 
       Log[xp]*((-4*(-1 + 2*xp)*(-1 + 2*zp)*Log[xp - zp])/(-1 + xp) - 
         (4*(-1 + 2*xp)*(-1 + 2*zp)*Log[zp])/(-1 + xp)))*Theta[xp - zp] + 
     ((8*(-1 + 2*xp)*zeta2*(-1 + 2*zp))/(-1 + xp) - 
       (4*(-1 + 2*xp)*(-1 + 2*zp)*Li[2, xp/zp])/(-1 + xp) + 
       (2*(-1 + 2*xp)*(-1 + 2*zp)*Log[xp]^2)/(-1 + xp) - 
       (2*(-1 + 2*xp)*(-1 + 2*zp)*Log[zp]^2)/(-1 + xp) - 
       (4*(-1 + 2*xp)*(-1 + 2*zp)*Log[xp]*Log[-xp + zp])/(-1 + xp) + 
       (4*(-1 + 2*xp)*(-1 + 2*zp)*Log[zp]*Log[-xp + zp])/(-1 + xp))*
      Theta[-xp + zp]), "\!\(\*SubscriptBox[SuperscriptBox[\(G\), \((2)\)], \
\({1,qq,NS,MS}\)]\)" -> -(Delz* Zq2xpNS) + 
   CF*nF*((-8*(-3 + 8*xp + 5*zp))/9 + Dxp[1]*((-4*(1 + zp))/3 + 
       (8*Dzp[0])/3) - (4*(1 + xp)*Dzp[1])/3 + 
     (Delz*((-2*(1 + xp))/3 + (4*Dxp[0])/3) + 
       Delx*(2*Delz - (2*(1 + zp))/3 + (4*Dzp[0])/3))*Log[Q2/muF2]^2 + 
     (8*(xp + zp)*Log[1 - xp])/3 - (4*(3 + 4*xp^2 + 4*xp*(-1 + zp) - zp)*
       Log[xp])/(3*(-1 + xp)) + Dzp[0]*((8*(1 + 4*xp))/9 - 
       (4*(1 + xp)*Log[1 - xp])/3 + (8*(1 + xp^2)*Log[xp])/(3*(-1 + xp))) + 
     Delz*((-2*(-23 + 79*xp))/27 + (4*(1 + xp)*zeta2)/3 + 
       (112/27 - (8*zeta2)/3)*Dxp[0] - (40*Dxp[1])/9 + (4*Dxp[2])/3 + 
       (4*(1 + xp^2)*Li[2, 1 - xp])/(3*(-1 + xp)) - 
       (2*(1 + xp)*Log[1 - xp]^2)/3 + ((18 - 20*xp + 22*xp^2)*Log[xp])/
        (3 - 3*xp) + ((5 + 5*xp^2)*Log[xp]^2)/(3 - 3*xp) + 
       Log[1 - xp]*((8*(1 + 4*xp))/9 + (8*(1 + xp^2)*Log[xp])/
          (3*(-1 + xp)))) + (8*(xp + zp)*Log[1 - zp])/3 + 
     Dxp[0]*((8*(1 + 4*zp))/9 - (40*Dzp[0])/9 + (8*Dzp[1])/3 - 
       (4*(1 + zp)*Log[1 - zp])/3) + Log[Q2/muF2]*((8*(xp + zp))/3 - 
       (4*(1 + xp)*Dzp[0])/3 + Dxp[0]*((-4*(1 + zp))/3 + (8*Dzp[0])/3) + 
       Delz*((8*(1 + 4*xp))/9 - (40*Dxp[0])/9 + (8*Dxp[1])/3 - 
         (4*(1 + xp)*Log[1 - xp])/3 + (8*(1 + xp^2)*Log[xp])/(3*(-1 + xp))) + 
       Delx*(Delz*(-34/3 - (16*zeta2)/3) + (8*(1 + 4*zp))/9 - (40*Dzp[0])/9 + 
         (8*Dzp[1])/3 - (4*(1 + zp)*Log[1 - zp])/3)) + 
     Delx*(Delz*(127/6 + (76*zeta2)/9 + (8*zeta3)/3) + (4*zeta2*(1 + zp))/3 - 
       (2*(19 + 37*zp))/27 + (112/27 - (8*zeta2)/3)*Dzp[0] - (40*Dzp[1])/9 + 
       (4*Dzp[2])/3 + (8*(1 + 4*zp)*Log[1 - zp])/9 - 
       (2*(1 + zp)*Log[1 - zp]^2)/3 + (2*(11 - 12*zp + 11*zp^2)*Log[zp])/
        (9*(-1 + zp)) + ((1 + zp^2)*Log[zp]^2)/(3*(-1 + zp)))) + 
   CF*((-8*zeta2*(-1 + xp^2 + 4*zp + xp*(2 - 4*zp + 4*zp^2)))/(-1 + xp) + 
     (-113 + 564*zp - 33*zp^2 - 130*zp^3 + xp*(127 - 300*zp - 141*zp^2 + 
         62*zp^3))/(9*zp) + 16*Sqrt[xp]*Sqrt[zp]*
      InvTanInt[-(Sqrt[zp]/Sqrt[xp])] - 16*Sqrt[xp]*Sqrt[zp]*
      InvTanInt[Sqrt[zp]/Sqrt[xp]] - 16*Sqrt[xp]*Sqrt[zp]*
      InvTanInt[-(Sqrt[xp]*Sqrt[zp])] + 16*Sqrt[xp]*Sqrt[zp]*
      InvTanInt[Sqrt[xp]*Sqrt[zp]] - 
     (4*(1 - 2*zp + 4*zp^2 + xp^2*(-1 + 2*zp + 4*zp^2))*Li[2, 1 - xp])/
      ((-1 + xp)*zp) + 16*(1 + xp)*(-1 + 2*zp)*Li[2, -xp] - 
     (4*(xp^2 - 4*zp + xp*(-2 + 6*zp))*Li[2, 1 - zp])/(-1 + xp) - 
     (8*(xp + xp^3 + 2*zp + 2*xp^2*zp + 4*xp*zp^2)*Li[2, -(zp/xp)])/
      (-1 + xp^2) - (8*(1 + 2*xp*zp + 2*xp^3*zp + xp^2*(1 + 4*zp^2))*
       Li[2, -(xp*zp)])/(-1 + xp^2) - (8*(1 + xp^2 + xp*(-1 + 2*zp))*
       Li[2, 1/2 - xp/2 - Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]/2])/
      Sqrt[1 + xp^2 + xp*(2 - 4*zp)] - (8*(1 + xp^2 + xp*(-1 + 2*zp))*
       Li[2, 1/2 + xp/2 - Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]/2])/
      Sqrt[1 + xp^2 + xp*(2 - 4*zp)] + (8*(1 + xp^2 + xp*(-1 + 2*zp))*
       Li[2, 1/2 - 1/(2*xp) - Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]/(2*xp)])/
      Sqrt[1 + xp^2 + xp*(2 - 4*zp)] + (8*(1 + xp^2 + xp*(-1 + 2*zp))*
       Li[2, 1/2 + 1/(2*xp) - Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]/(2*xp)])/
      Sqrt[1 + xp^2 + xp*(2 - 4*zp)] - 
     (4*(1 + 2*xp^2 + 2*zp + xp*(-2 + 4*zp))*
       Li[2, (2 - 2*xp)^(-1) + zp/(2 - 2*xp) - 
         Sqrt[1 - 2*zp + 4*xp*zp + zp^2]/(2 - 2*xp)])/(-1 + xp) + 
     (4*(1 + 2*zp + 4*xp^2*zp + xp*(2 - 4*zp + 8*zp^2))*
       Li[2, (2 - 2*xp)^(-1) - (2*xp)/(2 - 2*xp) - zp/(2 - 2*xp) + 
         Sqrt[1 - 2*zp + 4*xp*zp + zp^2]/(2 - 2*xp)])/(-1 + xp) - 
     (8*(1 + 2*zp + 4*xp*zp^2 + xp^2*(1 + 2*zp))*
       Li[2, 1 - 1/(2*xp) + zp/(2*xp) - Sqrt[1 - 2*zp + 4*xp*zp + zp^2]/
          (2*xp)])/(-1 + xp) + ((2 - 4*xp + 4*xp^2 + 4*zp + 8*xp*zp)*
       Log[2]^2)/(1 - xp) + (4*xp*(2 - 4*zp + 4*zp^2 + xp*(-1 + 2*zp))*
       Log[1 - xp]^2)/(-1 + xp) + 2*(6 + (-1 + xp)^(-1) + 
       4/Sqrt[1 + xp^2 + xp*(2 - 4*zp)] + 
       (4*xp^2)/Sqrt[1 + xp^2 + xp*(2 - 4*zp)] - 2/zp - 
       (2*(1 + 4*xp)*zp)/(-1 + xp) + 
       xp*(8 - 4/Sqrt[1 + xp^2 + xp*(2 - 4*zp)] - 2/zp + 
         (-12 + 8/Sqrt[1 + xp^2 + xp*(2 - 4*zp)])*zp))*Log[xp]^2 + 
     Dzp[1]*(-10*(-1 + xp) + 4*(1 + xp)*Log[xp]) + 
     Dzp[0]*(12*(-1 + xp) + 4*(1 + xp)*Li[2, 1 - xp] + 12*(-1 + xp)*Log[xp] - 
       4*(1 + xp)*Log[xp]^2 + Log[1 - xp]*(-10*(-1 + xp) + 
         4*(1 + xp)*Log[xp])) + Delz*(-36*(-1 + xp) + 10*(-1 + xp)*zeta2 + 
       2*(-1 + xp)*Li[2, 1 - xp] - 4*(1 + xp)*Li[3, 1 - xp] + 
       (25 - 3*xp - 4*(1 + xp)*zeta2 - 4*(1 + xp)*Li[2, 1 - xp])*Log[xp] + 
       ((17 - 15*xp)*Log[xp]^2)/2 + (5*(1 + xp)*Log[xp]^3)/3 + 
       Log[1 - xp]^2*(5 - 5*xp + 2*(1 + xp)*Log[xp]) + 
       Log[1 - xp]*(12*(-1 + xp) + 4*(1 + xp)*Li[2, 1 - xp] + 
         12*(-1 + xp)*Log[xp] - 4*(1 + xp)*Log[xp]^2)) + 
     (2*(3 + 2*xp*(1 - 2*zp)^2 + 6*zp + xp^2*(-2 + 4*zp))*Log[zp]^2)/
      (-1 + xp) + Dxp[1]*(2 + 8/(3*zp) - 2*zp - (8*zp^2)/3 + 
       4*(1 + zp)*Log[zp]) + Log[1 - zp]*
      ((26 - 78*zp + 6*zp^2 + 16*zp^3 + xp*(-34 + 54*zp + 18*zp^2 - 8*zp^3))/
        (3*zp) - 4*(1 + xp + 2*zp)*Log[zp]) + Log[Q2/muF2]^2*
      (Delz*(5 - 5*xp + 2*(1 + xp)*Log[xp]) + 
       Delx*(1 + 4/(3*zp) - zp - (4*zp^2)/3 + 2*(1 + zp)*Log[zp])) + 
     Dxp[0]*((2*(-7 - 60*zp + 42*zp^2 + 25*zp^3))/(9*zp) + 
       4*(1 + zp)*Li[2, 1 - zp] + (-2 + 8/(3*zp) - 8*zp - (8*zp^2)/3)*
        Log[zp] + 4*(1 + zp)*Log[zp]^2 + Log[1 - zp]*(2 + 8/(3*zp) - 2*zp - 
         (8*zp^2)/3 + 4*(1 + zp)*Log[zp])) + Log[Q2/muF2]*
      ((26 - 78*zp + 6*zp^2 + 16*zp^3 + xp*(-34 + 54*zp + 18*zp^2 - 8*zp^3))/
        (3*zp) - (4*(1 + xp)*(-1 + 2*zp)*Log[xp])/zp + 
       Dzp[0]*(-10*(-1 + xp) + 4*(1 + xp)*Log[xp]) + 
       Delz*(12*(-1 + xp) + 4*(1 + xp)*Li[2, 1 - xp] + 12*(-1 + xp)*Log[xp] - 
         4*(1 + xp)*Log[xp]^2 + Log[1 - xp]*(-10*(-1 + xp) + 
           4*(1 + xp)*Log[xp])) - 4*(1 + xp + 2*zp)*Log[zp] + 
       Dxp[0]*(2 + 8/(3*zp) - 2*zp - (8*zp^2)/3 + 4*(1 + zp)*Log[zp]) + 
       Delx*((2*(-7 - 60*zp + 42*zp^2 + 25*zp^3))/(9*zp) + 
         4*(1 + zp)*Li[2, 1 - zp] + (-2 + 8/(3*zp) - 8*zp - (8*zp^2)/3)*
          Log[zp] + 4*(1 + zp)*Log[zp]^2 + Log[1 - zp]*(2 + 8/(3*zp) - 2*zp - 
           (8*zp^2)/3 + 4*(1 + zp)*Log[zp]))) + 
     Log[zp]*((2*(18/(-1 + xp) - 6/(-1 + zp) + 13/zp + (21 + 18/(-1 + xp))*
           zp + 8*zp^2 + (12*xp^2)/(-xp + zp) + 
          xp*(6*(4 + (xp - zp)^(-1) + (-1 + zp)^(-1)) - 17/zp + 9*zp - 
            4*zp^2)))/3 + 16*Sqrt[xp]*Sqrt[zp]*ArcTan[Sqrt[zp]/Sqrt[xp]] - 
       16*Sqrt[xp]*Sqrt[zp]*ArcTan[Sqrt[xp]*Sqrt[zp]] - 
       (8*(xp + xp^3 + 2*zp + 2*xp^2*zp + 4*xp*zp^2)*Log[1 + zp/xp])/
        (-1 + xp^2) - (8*(1 + 2*xp*zp + 2*xp^3*zp + xp^2*(1 + 4*zp^2))*
         Log[1 + xp*zp])/(-1 + xp^2)) - (8*(1 + xp^2 + xp*(-1 + 2*zp))*
       Log[1 + xp - Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]]*
       Log[1 - xp + Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]])/
      Sqrt[1 + xp^2 + xp*(2 - 4*zp)] + (8*(1 + xp^2 + xp*(-1 + 2*zp))*
       Log[1 + xp - Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]]*
       Log[-1 + xp + Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]])/
      Sqrt[1 + xp^2 + xp*(2 - 4*zp)] - 
     (2*(1 + 2*zp + 4*xp^2*zp + xp*(2 - 4*zp + 8*zp^2))*
       Log[1 + zp - Sqrt[1 - 2*zp + 4*xp*zp + zp^2]]^2)/(-1 + xp) - 
     (2*(1 + 2*zp + 4*xp^2*zp + xp*(2 - 4*zp + 8*zp^2))*
       Log[1 - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]]^2)/(-1 + xp) + 
     Log[2]*((-4*(1 + 2*xp)*Sqrt[1 + (-2 + 4*xp)*zp + zp^2])/(-1 + xp) + 
       (8*xp*(2 - 4*zp + 4*zp^2 + xp*(-1 + 2*zp))*Log[1 - xp])/(-1 + xp) + 
       4*(-2 - 3/(-1 + xp) + 4/Sqrt[1 + xp^2 + xp*(2 - 4*zp)] - 
         (4*xp)/Sqrt[1 + xp^2 + xp*(2 - 4*zp)] + 
         (4*xp^2)/Sqrt[1 + xp^2 + xp*(2 - 4*zp)] - (2*zp)/(-1 + xp) + 
         xp*(-4 + 8/Sqrt[1 + xp^2 + xp*(2 - 4*zp)])*zp - 
         (8*xp*zp^2)/(-1 + xp))*Log[xp] + (8*(1 + xp^2 + xp*(-1 + 2*zp))*
         Log[1 - xp + Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]])/
        Sqrt[1 + xp^2 + xp*(2 - 4*zp)] - (8*(1 + xp^2 + xp*(-1 + 2*zp))*
         Log[-1 + xp + Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]])/
        Sqrt[1 + xp^2 + xp*(2 - 4*zp)] - 
       (8*xp*(2 - 4*zp + 4*zp^2 + xp*(-1 + 2*zp))*
         Log[1 + zp - Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/(-1 + xp) + 
       (4*(1 + 2*zp + 4*xp^2*zp + xp*(2 - 4*zp + 8*zp^2))*
         Log[1 - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/(-1 + xp)) + 
     (4*(1 + 2*xp)*Sqrt[1 + (-2 + 4*xp)*zp + zp^2]*
       Log[1 - 2*xp - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/(-1 + xp) + 
     ((2 + 4*xp^2 + 4*zp + xp*(-4 + 8*zp))*
       Log[1 - 2*xp - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]]^2)/(-1 + xp) + 
     Log[1 + zp - Sqrt[1 - 2*zp + 4*xp*zp + zp^2]]*
      ((-4*(1 + 2*xp)*Sqrt[1 + (-2 + 4*xp)*zp + zp^2])/(-1 + xp) + 
       (4*(1 + 2*zp + 4*xp^2*zp + xp*(2 - 4*zp + 8*zp^2))*
         Log[1 - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/(-1 + xp) + 
       (8*xp*(2 - 4*zp + 4*zp^2 + xp*(-1 + 2*zp))*
         Log[1 - 2*xp - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/(-1 + xp)) + 
     Log[1 - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]]*
      ((4*(1 + 2*xp)*Sqrt[1 + (-2 + 4*xp)*zp + zp^2])/(-1 + xp) - 
       (4*(1 + 2*zp + 4*xp^2*zp + xp*(2 - 4*zp + 8*zp^2))*
         Log[1 - 2*xp - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/(-1 + xp)) + 
     Log[1 - xp]*((26 - 78*zp + 6*zp^2 + 16*zp^3 + 
         xp*(-34 + 54*zp + 18*zp^2 - 8*zp^3))/(3*zp) - 
       (4*(1 + xp)*(-1 + 2*zp)*Log[xp])/zp - 4*(1 + xp + 2*zp)*Log[zp] + 
       (4*(1 + 2*xp^2 + 2*zp + xp*(-2 + 4*zp))*
         Log[1 + zp - Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/(-1 + xp) - 
       (4*(1 + 2*zp + 4*xp^2*zp + xp*(2 - 4*zp + 8*zp^2))*
         Log[1 - 2*xp - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/(-1 + xp)) + 
     Log[xp]*((2*((12*xp^2)/(xp - zp) + xp*(-27 + 22/zp + 8*zp^2 + 
            6/(-xp + zp) + zp*(-18 - 48/Sqrt[1 + (-2 + 4*xp)*zp + zp^2])) - 
          (22 + zp*(-48 - 6/Sqrt[1 - 2*zp + 4*xp*zp + zp^2]) + 
            zp^3*(8 - 6/Sqrt[1 - 2*zp + 4*xp*zp + zp^2]) + 
            zp^2*(-15 + 12/Sqrt[1 - 2*zp + 4*xp*zp + zp^2]) + 
            xp*(-14 + zp^2*(-18 - 48/Sqrt[1 - 2*zp + 4*xp*zp + zp^2]) + 
              zp^3*(-16 - 12/Sqrt[1 - 2*zp + 4*xp*zp + zp^2]) + 
              zp*(75 - 12/Sqrt[1 - 2*zp + 4*xp*zp + zp^2])))/(zp - xp*zp)))/
        3 - 16*Sqrt[xp]*Sqrt[zp]*ArcTan[Sqrt[zp]/Sqrt[xp]] - 
       16*Sqrt[xp]*Sqrt[zp]*ArcTan[Sqrt[xp]*Sqrt[zp]] + 
       16*(1 + xp)*(-1 + 2*zp)*Log[1 + xp] - 
       (4*(1 + xp)*(-1 + 2*zp)*Log[1 - zp])/zp + 
       (4*(-(-1 + zp)^2 + xp^3*(1 + zp + 4*zp^2) + xp^2*(1 + 5*zp + 8*zp^3) - 
          xp*(1 + 2*zp - 11*zp^2 + 8*zp^3))*Log[zp])/((-1 + xp)*(1 + xp)*
         zp) + (8*(xp + xp^3 + 2*zp + 2*xp^2*zp + 4*xp*zp^2)*Log[1 + zp/xp])/
        (-1 + xp^2) - (8*(1 + 2*xp*zp + 2*xp^3*zp + xp^2*(1 + 4*zp^2))*
         Log[1 + xp*zp])/(-1 + xp^2) - (8*(1 + xp^2 + xp*(-1 + 2*zp))*
         Log[-1 + xp + Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]])/
        Sqrt[1 + xp^2 + xp*(2 - 4*zp)] - (8*(1 + xp^2 + xp*(-1 + 2*zp))*
         Log[1 + xp + Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]])/
        Sqrt[1 + xp^2 + xp*(2 - 4*zp)] - 
       (4*(1 + 2*zp + 4*xp^2*zp + xp*(2 - 4*zp + 8*zp^2))*
         Log[1 + zp - Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/(-1 + xp) + 
       (4*(1 + 2*zp + 4*xp^2*zp + xp*(2 - 4*zp + 8*zp^2))*
         Log[1 - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/(-1 + xp) + 
       (4*(1 + 2*zp + 4*xp^2*zp + xp*(2 - 4*zp + 8*zp^2))*
         Log[1 - 2*xp - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/(-1 + xp)) + 
     Delx*(zeta2*(-2 - 8/(3*zp) + 2*zp + (8*zp^2)/3) - 
       (107 + 132*zp + 48*zp^2 - 287*zp^3)/(27*zp) + 
       (-2 + 8/(3*zp) - 8*zp - (8*zp^2)/3)*Li[2, 1 - zp] - 
       4*(1 + zp)*Li[3, 1 - zp] + (-4*zeta2*(1 + zp) - 
         (2*(7 + 90*zp + 81*zp^2 + 31*zp^3))/(9*zp) + 
         8*(1 + zp)*Li[2, 1 - zp])*Log[zp] + ((-3 + 8/zp - 15*zp)*Log[zp]^2)/
        6 + (5*(1 + zp)*Log[zp]^3)/3 + Log[1 - zp]^2*(1 + 4/(3*zp) - zp - 
         (4*zp^2)/3 + 2*(1 + zp)*Log[zp]) + Log[1 - zp]*
        ((2*(-7 - 60*zp + 42*zp^2 + 25*zp^3))/(9*zp) + 
         4*(1 + zp)*Li[2, 1 - zp] + (-2 + 8/(3*zp) - 8*zp - (8*zp^2)/3)*
          Log[zp] + 4*(1 + zp)*Log[zp]^2) + 8*(1 + zp)*S12[1 - zp]) + 
     ((-4*(-1 + 2*xp)*(-1 + 2*zp)*Li[2, zp/xp])/(-1 + xp) - 
       (4*(-1 + 2*xp)*(-1 + 2*zp)*Log[xp]^2)/(-1 + xp) - 
       (4*(-1 + 2*xp)*(-1 + 2*zp)*Log[xp - zp]*Log[zp])/(-1 + xp) + 
       Log[xp]*((4*(-1 + 2*xp)*(-1 + 2*zp)*Log[xp - zp])/(-1 + xp) + 
         (4*(-1 + 2*xp)*(-1 + 2*zp)*Log[zp])/(-1 + xp)))*Theta[xp - zp] + 
     ((-8*(-1 + 2*xp)*zeta2*(-1 + 2*zp))/(-1 + xp) + 
       (4*(-1 + 2*xp)*(-1 + 2*zp)*Li[2, xp/zp])/(-1 + xp) - 
       (2*(-1 + 2*xp)*(-1 + 2*zp)*Log[xp]^2)/(-1 + xp) + 
       (2*(-1 + 2*xp)*(-1 + 2*zp)*Log[zp]^2)/(-1 + xp) + 
       (4*(-1 + 2*xp)*(-1 + 2*zp)*Log[xp]*Log[-xp + zp])/(-1 + xp) - 
       (4*(-1 + 2*xp)*(-1 + 2*zp)*Log[zp]*Log[-xp + zp])/(-1 + xp))*
      Theta[-xp + zp]) + 
   CF^2*((-8*zeta2*(-2 + (3 - 2*xp - 3*xp^2)*zp + 2*(-1 + xp + xp^2)*zp^2 + 
        (1 - 4*xp + xp^2)*zp^3))/((-1 + xp)*(-1 + zp)*zp) + 
     4*(5 - 2/(-1 + xp) - 8/(-1 + zp) - 12/zp + (-2 + 6/(-1 + xp))*zp - 
       (-1 + xp + zp)^(-1) + 2/(1 + 2*xp + xp^2 - 4*xp*zp) - 
       (2*xp^3)/(1 + 2*xp + xp^2 - 4*xp*zp) + 
       xp*(-32 + 8/(-1 + zp) + 11/zp - 4*zp + (-xp + zp)^(-1) + 
         (-1 + xp + zp)^(-1) - 2/(1 + 2*xp + xp^2 - 4*xp*zp)) + 
       2*xp^2*((xp - zp)^(-1) + (1 + 2*xp + xp^2 - 4*xp*zp)^(-1))) + 
     Dxp[2]*(-12*(1 + zp) + 24*Dzp[0]) - 12*(1 + xp)*Dzp[2] - 
     (8*(3*zp + 3*xp^2*zp + xp*(-1 + 3*zp^2))*
       InvTanInt[-(Sqrt[zp]/Sqrt[xp])])/(Sqrt[xp]*zp^(3/2)) + 
     (8*(3*zp + 3*xp^2*zp + xp*(-1 + 3*zp^2))*InvTanInt[Sqrt[zp]/Sqrt[xp]])/
      (Sqrt[xp]*zp^(3/2)) + (8*(3*zp + 3*xp^2*zp + xp*(-1 + 3*zp^2))*
       InvTanInt[-(Sqrt[xp]*Sqrt[zp])])/(Sqrt[xp]*zp^(3/2)) - 
     (8*(3*zp + 3*xp^2*zp + xp*(-1 + 3*zp^2))*InvTanInt[Sqrt[xp]*Sqrt[zp]])/
      (Sqrt[xp]*zp^(3/2)) + 
     (8*(1 + 2*zp - 3*zp^2 + xp^2*(1 + 2*zp - 3*zp^2) + 
        xp*zp*(-1 + 4*zp + zp^2))*Li[2, 1 - xp])/((-1 + xp)*(-1 + zp)*zp) - 
     (16*(1 + xp)*(-1 + 2*zp)*Li[2, -xp])/zp - 
     (8*(-1 + 5*zp - 3*zp^2 + xp^2*(-4 + 5*zp) + xp*(5 - 10*zp + 9*zp^2))*
       Li[2, 1 - zp])/((-1 + xp)*(-1 + zp)) - 
     (16*(1 + 2*xp*zp + 2*xp^3*zp + xp^2*(1 + 2*zp^2))*Li[2, -(zp/xp)])/
      ((-1 + xp^2)*zp) - (16*(xp + xp^3 + 2*zp + 2*xp^2*zp + 2*xp*zp^2)*
       Li[2, -(xp*zp)])/((-1 + xp)*(1 + xp)*zp) - 
     (8*((1 - 2*zp)^2 + xp^6*(1 - 2*zp)^2 + xp*(2 - 20*zp + 48*zp^2 - 
          32*zp^3) + xp^5*(2 - 20*zp + 48*zp^2 - 32*zp^3) + 
        xp^2*(-1 - 16*zp + 104*zp^2 - 176*zp^3 + 88*zp^4) + 
        xp^4*(-1 - 16*zp + 104*zp^2 - 176*zp^3 + 88*zp^4) - 
        4*xp^3*(1 - 12*zp^2 + 28*zp^3 - 30*zp^4 + 12*zp^5))*
       Li[2, 1/2 - xp/2 - Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]/2])/
      ((1 + xp^2 + xp*(2 - 4*zp))^(5/2)*(-1 + zp)*zp) - 
     (8*((1 - 2*zp)^2 + xp^6*(1 - 2*zp)^2 + xp*(2 - 20*zp + 48*zp^2 - 
          32*zp^3) + xp^5*(2 - 20*zp + 48*zp^2 - 32*zp^3) + 
        xp^2*(-1 - 16*zp + 104*zp^2 - 176*zp^3 + 88*zp^4) + 
        xp^4*(-1 - 16*zp + 104*zp^2 - 176*zp^3 + 88*zp^4) - 
        4*xp^3*(1 - 12*zp^2 + 28*zp^3 - 30*zp^4 + 12*zp^5))*
       Li[2, 1/2 + xp/2 - Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]/2])/
      ((1 + xp^2 + xp*(2 - 4*zp))^(5/2)*(-1 + zp)*zp) + 
     (8*((1 - 2*zp)^2 + xp^6*(1 - 2*zp)^2 + xp*(2 - 20*zp + 48*zp^2 - 
          32*zp^3) + xp^5*(2 - 20*zp + 48*zp^2 - 32*zp^3) + 
        xp^2*(-1 - 16*zp + 104*zp^2 - 176*zp^3 + 88*zp^4) + 
        xp^4*(-1 - 16*zp + 104*zp^2 - 176*zp^3 + 88*zp^4) - 
        4*xp^3*(1 - 12*zp^2 + 28*zp^3 - 30*zp^4 + 12*zp^5))*
       Li[2, 1/2 - 1/(2*xp) - Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]/(2*xp)])/
      ((1 + xp^2 + xp*(2 - 4*zp))^(5/2)*(-1 + zp)*zp) + 
     (8*((1 - 2*zp)^2 + xp^6*(1 - 2*zp)^2 + xp*(2 - 20*zp + 48*zp^2 - 
          32*zp^3) + xp^5*(2 - 20*zp + 48*zp^2 - 32*zp^3) + 
        xp^2*(-1 - 16*zp + 104*zp^2 - 176*zp^3 + 88*zp^4) + 
        xp^4*(-1 - 16*zp + 104*zp^2 - 176*zp^3 + 88*zp^4) - 
        4*xp^3*(1 - 12*zp^2 + 28*zp^3 - 30*zp^4 + 12*zp^5))*
       Li[2, 1/2 + 1/(2*xp) - Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]/(2*xp)])/
      ((1 + xp^2 + xp*(2 - 4*zp))^(5/2)*(-1 + zp)*zp) - 
     (8*(1 + 2*zp + 2*xp*zp^2 + xp^2*(1 + 2*zp))*
       Li[2, (2 - 2*xp)^(-1) + zp/(2 - 2*xp) - 
         Sqrt[1 - 2*zp + 4*xp*zp + zp^2]/(2 - 2*xp)])/((-1 + xp)*zp) + 
     (8*(1 + 2*zp + 2*xp*zp^2 + xp^2*(1 + 2*zp))*
       Li[2, (2 - 2*xp)^(-1) - (2*xp)/(2 - 2*xp) - zp/(2 - 2*xp) + 
         Sqrt[1 - 2*zp + 4*xp*zp + zp^2]/(2 - 2*xp)])/((-1 + xp)*zp) - 
     (16*(1 + 2*zp + 2*xp*zp^2 + xp^2*(1 + 2*zp))*
       Li[2, 1 - 1/(2*xp) + zp/(2*xp) - Sqrt[1 - 2*zp + 4*xp*zp + zp^2]/
          (2*xp)])/((-1 + xp)*zp) + ((4 + 8*zp + 8*xp*zp^2 + xp^2*(4 + 8*zp))*
       Log[2]^2)/(zp - xp*zp) + (4*(4*xp*(-1 + zp)^2 + (4 - 3*zp)*zp + 
        xp^2*(-4 + 4*zp + zp^2))*Log[1 - xp]^2)/((-1 + xp)*(-1 + zp)) + 
     2*(26 + 48/(-1 + xp) + 3/(1 + xp^2 + xp*(2 - 4*zp))^(5/2) + 
       (3*xp^6)/(1 + xp^2 + xp*(2 - 4*zp))^(5/2) - 
       3/(1 + xp^2 + xp*(2 - 4*zp))^(3/2) + 
       16/Sqrt[1 + xp^2 + xp*(2 - 4*zp)] - 
       (3*xp^4*(2 + xp^2 + xp*(2 - 4*zp)))/(1 + xp^2 + xp*(2 - 4*zp))^(5/2) + 
       (xp^2*(-3 - 2*(1 + xp^2 + xp*(2 - 4*zp)) + 
          16*(1 + xp^2 + xp*(2 - 4*zp))^2 + (4*(1 + xp^2 + xp*(2 - 4*zp))^2)/
           (-1 + zp) - (4*(1 + xp^2 + xp*(2 - 4*zp))^2)/zp))/
        (1 + xp^2 + xp*(2 - 4*zp))^(5/2) + 44/(-1 + zp) + 
       64/((-1 + xp)*(-1 + zp)) + 4/(Sqrt[1 + xp^2 + xp*(2 - 4*zp)]*
         (-1 + zp)) + (62 + 52/(-1 + xp))*zp + 4/(zp - xp*zp) + 
       2*xp*(21 + 3/Sqrt[1 + xp^2 + xp*(2 - 4*zp)] + 10/(-1 + zp) + 
         4/(Sqrt[1 + xp^2 + xp*(2 - 4*zp)]*(-1 + zp)) + 
         (-2 + 4/Sqrt[1 + xp^2 + xp*(2 - 4*zp)])/zp + zp - 
         (6*zp)/Sqrt[1 + xp^2 + xp*(2 - 4*zp)]) - 
       (4 + 4/Sqrt[1 + 2*xp + xp^2 - 4*xp*zp])/zp)*Log[xp]^2 + 
     Dzp[1]*(-24*(1 + xp)*Log[1 - xp] + (4*(3 + 5*xp^2)*Log[xp])/(-1 + xp)) + 
     Dzp[0]*(8 + 56*xp + 16*(1 + xp)*zeta2 + 4*(1 + xp)*Li[2, 1 - xp] - 
       12*(1 + xp)*Log[1 - xp]^2 + 12*((-1 + xp)^(-1) - xp)*Log[xp] + 
       (4*(5 + 7*xp^2)*Log[1 - xp]*Log[xp])/(-1 + xp) - 
       (8*xp^2*Log[xp]^2)/(-1 + xp)) + 
     (4*(-3 + 10*xp*(-1 + zp)^2 + 10*zp - 3*zp^2 + xp^2*(-7 + 10*zp + zp^2))*
       Log[1 - zp]^2)/((-1 + xp)*(-1 + zp)) - 
     (4*(1 + 4*zp - 8*zp^2 - 5*zp^3 - 2*xp*zp*(3 - 7*zp + 3*zp^2) + 
        xp^2*(1 + 4*zp - 8*zp^2 + zp^3))*Log[zp]^2)/((-1 + xp)*(-1 + zp)*
       zp) + Dxp[1]*(48*Dzp[1] - 24*(1 + zp)*Log[1 - zp] + 
       4*(1 + zp)*Log[zp]) + Log[1 - zp]*
      (4*(18 - 8/(-1 + xp) + 12/(-1 + zp) - 4*zp + xp^2/(-1 + xp + zp)^2 - 
         3/(-1 + xp + zp) + xp*(-20 - 12/(-1 + zp) + 4*zp - 
           (-1 + xp + zp)^(-2) + (-1 + xp + zp)^(-1))) + 
       (8*(-4 + 7*xp*(-1 + zp)^2 + 7*zp - zp^2 + xp^2*(-3 + 7*zp))*Log[zp])/
        ((-1 + xp)*(-1 + zp))) + Dxp[0]*(8 + 56*zp + 16*zeta2*(1 + zp) + 
       (-64 - 32*zeta2)*Dzp[0] + 24*Dzp[2] + 4*(1 + zp)*Li[2, 1 - zp] - 
       12*(1 + zp)*Log[1 - zp]^2 + (-28 - 12/(-1 + zp) - 8*zp)*Log[zp] - 
       (4*(3 + zp^2)*Log[1 - zp]*Log[zp])/(-1 + zp) + 
       (8*(1 + 2*zp^2)*Log[zp]^2)/(-1 + zp)) + Log[Q2/muF2]^2*
      (4*(1 + xp)*(1 + zp) - 8*(1 + xp)*Dzp[0] + 
       Dxp[0]*(-8*(1 + zp) + 16*Dzp[0]) + Delz*(-8*(2 + xp) + 24*Dxp[0] + 
         16*Dxp[1] - 8*(1 + xp)*Log[1 - xp] + ((2 + 6*xp^2)*Log[xp])/
          (-1 + xp)) + Delx*(Delz*(18 - 16*zeta2) - 8*(2 + zp) + 24*Dzp[0] + 
         16*Dzp[1] - 8*(1 + zp)*Log[1 - zp] + ((2 + 6*zp^2)*Log[zp])/
          (-1 + zp))) + Log[Q2/muF2]*(8*(2 + xp + zp + 2*xp*zp) + 
       Dxp[1]*(-24*(1 + zp) + 48*Dzp[0]) - 24*(1 + xp)*Dzp[1] + 
       8*(1 + 3*zp + xp*(3 + zp))*Log[1 - xp] - 
       (8*(1 + xp*(-1 + zp) + xp^2*(2 + zp))*Log[xp])/(-1 + xp) + 
       Dzp[0]*(-12*(1 + xp) - 24*(1 + xp)*Log[1 - xp] + 
         (4*(3 + 5*xp^2)*Log[xp])/(-1 + xp)) + 
       Delz*(20 + 44*xp + 16*(1 + xp)*zeta2 + (-64 - 32*zeta2)*Dxp[0] + 
         24*Dxp[1] + 24*Dxp[2] + 4*(1 + xp)*Li[2, 1 - xp] - 
         12*(1 + xp)*Log[1 - xp]^2 + (12*(2 + xp)*Log[xp])/(-1 + xp) - 
         (8*xp^2*Log[xp]^2)/(-1 + xp) + Log[1 - xp]*(-12*(1 + xp) + 
           (4*(5 + 7*xp^2)*Log[xp])/(-1 + xp))) + 8*(1 + 3*zp + xp*(3 + zp))*
        Log[1 - zp] - 8*(xp + zp)*Log[zp] + 
       Dxp[0]*(-12*(1 + zp) + 24*Dzp[0] + 48*Dzp[1] - 
         24*(1 + zp)*Log[1 - zp] + 4*(1 + zp)*Log[zp]) + 
       Delx*(20 + Delz*(-93 - 24*zeta2 + 80*zeta3) + 44*zp + 
         16*zeta2*(1 + zp) + (-64 - 32*zeta2)*Dzp[0] + 24*Dzp[1] + 
         24*Dzp[2] + 4*(1 + zp)*Li[2, 1 - zp] - 12*(1 + zp)*Log[1 - zp]^2 - 
         (4*(-1 + 5*zp + 5*zp^2)*Log[zp])/(-1 + zp) + 
         (8*(1 + 2*zp^2)*Log[zp]^2)/(-1 + zp) + Log[1 - zp]*
          (-12*(1 + zp) - (4*(3 + zp^2)*Log[zp])/(-1 + zp)))) + 
     Log[zp]*(4*(33 + 4/(-1 + xp) - (3*xp^4)/(1 + xp^2 + xp*(2 - 4*zp))^2 - 
         (3*xp^5)/(1 + xp^2 + xp*(2 - 4*zp))^2 - 3/((-1 + xp)*(xp - zp)) + 
         12/(-1 + zp) - 3/((-1 + xp)*(-1 + zp)) + 6*(2 + (-1 + xp)^(-1))*zp + 
         3/(-xp + zp) + (xp^2*(-1 + 6*zp))/(xp - zp)^2 + 
         3/(1 + 2*xp + xp^2 - 4*xp*zp)^2 - 2/(1 + 2*xp + xp^2 - 4*xp*zp) - 
         6/(zp - xp*zp) + xp*(-17 - 12/(-1 + zp) + 6/zp + 3/(-xp + zp) + 
           3/(1 + 2*xp + xp^2 - 4*xp*zp)^2) + xp^3*(-4/(xp - zp)^2 + 
           2/(1 + 2*xp + xp^2 - 4*xp*zp))) - 
       (8*(3*zp + 3*xp^2*zp + xp*(-1 + 3*zp^2))*ArcTan[Sqrt[zp]/Sqrt[xp]])/
        (Sqrt[xp]*zp^(3/2)) + (8*(3*zp + 3*xp^2*zp + xp*(-1 + 3*zp^2))*
         ArcTan[Sqrt[xp]*Sqrt[zp]])/(Sqrt[xp]*zp^(3/2)) - 
       (16*(1 + 2*xp*zp + 2*xp^3*zp + xp^2*(1 + 2*zp^2))*Log[1 + zp/xp])/
        ((-1 + xp^2)*zp) - (16*(xp + xp^3 + 2*zp + 2*xp^2*zp + 2*xp*zp^2)*
         Log[1 + xp*zp])/((-1 + xp)*(1 + xp)*zp)) - 
     (8*((1 - 2*zp)^2 + xp^6*(1 - 2*zp)^2 + xp*(2 - 20*zp + 48*zp^2 - 
          32*zp^3) + xp^5*(2 - 20*zp + 48*zp^2 - 32*zp^3) + 
        xp^2*(-1 - 16*zp + 104*zp^2 - 176*zp^3 + 88*zp^4) + 
        xp^4*(-1 - 16*zp + 104*zp^2 - 176*zp^3 + 88*zp^4) - 
        4*xp^3*(1 - 12*zp^2 + 28*zp^3 - 30*zp^4 + 12*zp^5))*
       Log[1 + xp - Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]]*
       Log[1 - xp + Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]])/
      ((1 + xp^2 + xp*(2 - 4*zp))^(5/2)*(-1 + zp)*zp) + 
     (8*((1 - 2*zp)^2 + xp^6*(1 - 2*zp)^2 + xp*(2 - 20*zp + 48*zp^2 - 
          32*zp^3) + xp^5*(2 - 20*zp + 48*zp^2 - 32*zp^3) + 
        xp^2*(-1 - 16*zp + 104*zp^2 - 176*zp^3 + 88*zp^4) + 
        xp^4*(-1 - 16*zp + 104*zp^2 - 176*zp^3 + 88*zp^4) - 
        4*xp^3*(1 - 12*zp^2 + 28*zp^3 - 30*zp^4 + 12*zp^5))*
       Log[1 + xp - Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]]*
       Log[-1 + xp + Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]])/
      ((1 + xp^2 + xp*(2 - 4*zp))^(5/2)*(-1 + zp)*zp) + 
     ((4 + 8*zp + 8*xp*zp^2 + xp^2*(4 + 8*zp))*
       Log[1 + zp - Sqrt[1 - 2*zp + 4*xp*zp + zp^2]]^2)/(zp - xp*zp) + 
     ((4 + 8*zp + 8*xp*zp^2 + xp^2*(4 + 8*zp))*
       Log[1 - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]]^2)/(zp - xp*zp) + 
     Log[1 + zp - Sqrt[1 - 2*zp + 4*xp*zp + zp^2]]*
      ((-8*(2*(-1 + zp)^2 + 8*xp^3*zp + 2*xp^2*(1 - 4*zp + 7*zp^2) + 
          xp*(-1 + 13*zp - 7*zp^2 + 3*zp^3)))/((-1 + xp)*zp*
         Sqrt[1 + (-2 + 4*xp)*zp + zp^2]) + 
       (8*(1 + 2*zp + 2*xp*zp^2 + xp^2*(1 + 2*zp))*
         Log[1 - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/((-1 + xp)*zp)) + 
     Log[2]*((-8*(2*(-1 + zp)^2 + 8*xp^3*zp + 2*xp^2*(1 - 4*zp + 7*zp^2) + 
          xp*(-1 + 13*zp - 7*zp^2 + 3*zp^3)))/((-1 + xp)*zp*
         Sqrt[1 + (-2 + 4*xp)*zp + zp^2]) + 
       4*(-4 - 8/(-1 + xp) + 3/(1 + xp^2 + xp*(2 - 4*zp))^(5/2) + 
         (3*xp^6)/(1 + xp^2 + xp*(2 - 4*zp))^(5/2) - 
         3/(1 + xp^2 + xp*(2 - 4*zp))^(3/2) + 
         16/Sqrt[1 + xp^2 + xp*(2 - 4*zp)] - 
         (3*xp^4*(2 + xp^2 + xp*(2 - 4*zp)))/(1 + xp^2 + xp*(2 - 4*zp))^
           (5/2) + (xp^2*(-3 - 2*(1 + xp^2 + xp*(2 - 4*zp)) + 
            16*(1 + xp^2 + xp*(2 - 4*zp))^2 + (4*(1 + xp^2 + xp*(2 - 4*zp))^
               2)/(-1 + zp) - (4*(1 + xp^2 + xp*(2 - 4*zp))^2)/zp))/
          (1 + xp^2 + xp*(2 - 4*zp))^(5/2) + 
         4/(Sqrt[1 + xp^2 + xp*(2 - 4*zp)]*(-1 + zp)) - (4*xp*zp)/(-1 + xp) + 
         4/(zp - xp*zp) + xp*(-4 + 6/Sqrt[1 + xp^2 + xp*(2 - 4*zp)] + 
           8/(Sqrt[1 + xp^2 + xp*(2 - 4*zp)]*(-1 + zp)) + 
           (-2 + 8/Sqrt[1 + xp^2 + xp*(2 - 4*zp)])/zp - 
           (12*zp)/Sqrt[1 + xp^2 + xp*(2 - 4*zp)]) - 
         (2 + 4/Sqrt[1 + 2*xp + xp^2 - 4*xp*zp])/zp)*Log[xp] + 
       (8*((1 - 2*zp)^2 + xp^6*(1 - 2*zp)^2 + xp*(2 - 20*zp + 48*zp^2 - 
            32*zp^3) + xp^5*(2 - 20*zp + 48*zp^2 - 32*zp^3) + 
          xp^2*(-1 - 16*zp + 104*zp^2 - 176*zp^3 + 88*zp^4) + 
          xp^4*(-1 - 16*zp + 104*zp^2 - 176*zp^3 + 88*zp^4) - 
          4*xp^3*(1 - 12*zp^2 + 28*zp^3 - 30*zp^4 + 12*zp^5))*
         Log[1 - xp + Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]])/
        ((1 + xp^2 + xp*(2 - 4*zp))^(5/2)*(-1 + zp)*zp) - 
       (8*((1 - 2*zp)^2 + xp^6*(1 - 2*zp)^2 + xp*(2 - 20*zp + 48*zp^2 - 
            32*zp^3) + xp^5*(2 - 20*zp + 48*zp^2 - 32*zp^3) + 
          xp^2*(-1 - 16*zp + 104*zp^2 - 176*zp^3 + 88*zp^4) + 
          xp^4*(-1 - 16*zp + 104*zp^2 - 176*zp^3 + 88*zp^4) - 
          4*xp^3*(1 - 12*zp^2 + 28*zp^3 - 30*zp^4 + 12*zp^5))*
         Log[-1 + xp + Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]])/
        ((1 + xp^2 + xp*(2 - 4*zp))^(5/2)*(-1 + zp)*zp) + 
       (8*(1 + 2*zp + 2*xp*zp^2 + xp^2*(1 + 2*zp))*
         Log[1 - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/((-1 + xp)*zp)) + 
     (8*(2*(-1 + zp)^2 + 8*xp^3*zp + 2*xp^2*(1 - 4*zp + 7*zp^2) + 
        xp*(-1 + 13*zp - 7*zp^2 + 3*zp^3))*
       Log[1 - 2*xp - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/
      ((-1 + xp)*zp*Sqrt[1 + (-2 + 4*xp)*zp + zp^2]) + 
     ((4 + 8*zp + 8*xp*zp^2 + xp^2*(4 + 8*zp))*
       Log[1 - 2*xp - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]]^2)/
      ((-1 + xp)*zp) + Log[1 - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]]*
      ((8*(2*(-1 + zp)^2 + 8*xp^3*zp + 2*xp^2*(1 - 4*zp + 7*zp^2) + 
          xp*(-1 + 13*zp - 7*zp^2 + 3*zp^3)))/((-1 + xp)*zp*
         Sqrt[1 + (-2 + 4*xp)*zp + zp^2]) - 
       (8*(1 + 2*zp + 2*xp*zp^2 + xp^2*(1 + 2*zp))*
         Log[1 - 2*xp - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/
        ((-1 + xp)*zp)) + Log[1 - xp]*
      (4*(11 - 4/(-1 + xp) + 6/(-1 + zp) - zp + 
         xp*(-11 - 6/(-1 + zp) + zp)) - 
       (8*(-4 + 9*zp - 2*zp^2 + xp^2*(-7 + 9*zp + zp^2) + 
          xp*(9 - 18*zp + 11*zp^2))*Log[xp])/((-1 + xp)*(-1 + zp)) + 
       (8*(-1 + 6*xp*(-1 + zp)^2 + 6*zp - 3*zp^2 + xp^2*(-5 + 6*zp + zp^2))*
         Log[1 - zp])/((-1 + xp)*(-1 + zp)) + 
       (8*(-3 + 5*zp + zp^2 + xp^2*(-2 + 5*zp) + xp*(5 - 10*zp + 3*zp^2))*
         Log[zp])/((-1 + xp)*(-1 + zp)) + 
       (8*(1 + 2*zp + 2*xp*zp^2 + xp^2*(1 + 2*zp))*
         Log[1 + zp - Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/((-1 + xp)*zp) - 
       (8*(1 + 2*zp + 2*xp*zp^2 + xp^2*(1 + 2*zp))*
         Log[1 - 2*xp - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/
        ((-1 + xp)*zp)) + Log[xp]*
      (4*(-4 + 10/(-1 + xp) + 3/(1 + xp^2 + xp*(2 - 4*zp))^2 - 
         (3*xp^4)/(1 + xp^2 + xp*(2 - 4*zp))^2 + 
         (3*xp^5)/(1 + xp^2 + xp*(2 - 4*zp))^2 + 3/(xp - zp) + 
         3/((-1 + xp)*(xp - zp)) - 18/(-1 + zp) + 3/(-1 + xp + zp) - 
         2/(1 + 2*xp + xp^2 - 4*xp*zp) - 9/(zp - xp*zp) - 
         26/Sqrt[1 + (-2 + 4*xp)*zp + zp^2] - 
         18/((-1 + xp)*Sqrt[1 + (-2 + 4*xp)*zp + zp^2]) - 
         (6*xp*zp^2)/((-1 + xp)*Sqrt[1 + (-2 + 4*xp)*zp + zp^2]) + 
         6/((zp - xp*zp)*Sqrt[1 + (-2 + 4*xp)*zp + zp^2]) + 
         xp^3*(-2/(xp - zp)^2 - 2/(1 + 2*xp + xp^2 - 4*xp*zp)) - 
         (2 + 2/Sqrt[1 - 2*zp + 4*xp*zp + zp^2])/zp + 
         xp^2*((xp - zp)^(-2) + 6/(xp - zp) - (-1 + xp + zp)^(-2) - 
           16/Sqrt[1 + (-2 + 4*xp)*zp + zp^2]) + 
         xp*(26 - 3/(1 + xp^2 + xp*(2 - 4*zp))^2 + 3/(xp - zp) + 
           18/(-1 + zp) + (-1 + xp + zp)^(-2) - (-1 + xp + zp)^(-1) - 
           4/(zp*Sqrt[1 + (-2 + 4*xp)*zp + zp^2]) + 
           zp*(-1 - 28/Sqrt[1 + (-2 + 4*xp)*zp + zp^2])) + 
         (zp*(4 + 4/Sqrt[1 + (-2 + 4*xp)*zp + zp^2] + 
            xp*(11 + 14/Sqrt[1 + (-2 + 4*xp)*zp + zp^2])))/(1 - xp)) + 
       (8*(3*zp + 3*xp^2*zp + xp*(-1 + 3*zp^2))*ArcTan[Sqrt[zp]/Sqrt[xp]])/
        (Sqrt[xp]*zp^(3/2)) + (8*(3*zp + 3*xp^2*zp + xp*(-1 + 3*zp^2))*
         ArcTan[Sqrt[xp]*Sqrt[zp]])/(Sqrt[xp]*zp^(3/2)) - 
       (16*(1 + xp)*(-1 + 2*zp)*Log[1 + xp])/zp - 
       (8*(-5 + 13*zp - 2*zp^2 + xp^2*(-8 + 13*zp + zp^2) + 
          xp*(13 - 26*zp + 15*zp^2))*Log[1 - zp])/((-1 + xp)*(-1 + zp)) - 
       (4*(zp*(3 + 14*zp - 3*zp^2) + 2*xp^3*(2 - 8*zp + 13*zp^2) + 
          xp*(8 + 5*zp - 10*zp^2 + 19*zp^3) + xp^2*(-4 + 34*zp - 46*zp^2 + 
            38*zp^3))*Log[zp])/((-1 + xp)*(1 + xp)*(-1 + zp)*zp) + 
       (16*(1 + 2*xp*zp + 2*xp^3*zp + xp^2*(1 + 2*zp^2))*Log[1 + zp/xp])/
        ((-1 + xp)*(1 + xp)*zp) - (16*(xp + xp^3 + 2*zp + 2*xp^2*zp + 
          2*xp*zp^2)*Log[1 + xp*zp])/((-1 + xp)*(1 + xp)*zp) - 
       (8*((1 - 2*zp)^2 + xp^6*(1 - 2*zp)^2 + xp*(2 - 20*zp + 48*zp^2 - 
            32*zp^3) + xp^5*(2 - 20*zp + 48*zp^2 - 32*zp^3) + 
          xp^2*(-1 - 16*zp + 104*zp^2 - 176*zp^3 + 88*zp^4) + 
          xp^4*(-1 - 16*zp + 104*zp^2 - 176*zp^3 + 88*zp^4) - 
          4*xp^3*(1 - 12*zp^2 + 28*zp^3 - 30*zp^4 + 12*zp^5))*
         Log[-1 + xp + Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]])/
        ((1 + xp^2 + xp*(2 - 4*zp))^(5/2)*(-1 + zp)*zp) - 
       (8*((1 - 2*zp)^2 + xp^6*(1 - 2*zp)^2 + xp*(2 - 20*zp + 48*zp^2 - 
            32*zp^3) + xp^5*(2 - 20*zp + 48*zp^2 - 32*zp^3) + 
          xp^2*(-1 - 16*zp + 104*zp^2 - 176*zp^3 + 88*zp^4) + 
          xp^4*(-1 - 16*zp + 104*zp^2 - 176*zp^3 + 88*zp^4) - 
          4*xp^3*(1 - 12*zp^2 + 28*zp^3 - 30*zp^4 + 12*zp^5))*
         Log[1 + xp + Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]])/
        ((1 + xp^2 + xp*(2 - 4*zp))^(5/2)*(-1 + zp)*zp) - 
       (8*(1 + 2*zp + 2*xp*zp^2 + xp^2*(1 + 2*zp))*
         Log[1 + zp - Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/((-1 + xp)*zp) + 
       (8*(1 + 2*zp + 2*xp*zp^2 + xp^2*(1 + 2*zp))*
         Log[1 - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/((-1 + xp)*zp) + 
       (8*(1 + 2*zp + 2*xp*zp^2 + xp^2*(1 + 2*zp))*
         Log[1 - 2*xp - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/
        ((-1 + xp)*zp)) + Delz*(-74 + 78*xp + (8 - 8*xp)*zeta2 - 
       16*(1 + xp)*zeta3 + 32*zeta3*Dxp[0] + (-64 - 32*zeta2)*Dxp[1] + 
       8*Dxp[3] + 12*((-1 + xp)^(-1) - xp)*Li[2, 1 - xp] - 
       (8*(1 + 2*xp^2)*Li[3, 1 - xp])/(-1 + xp) - 4*(1 + xp)*Log[1 - xp]^3 + 
       (-58 - 80/(-1 + xp) - 96*xp - (4*(3 + 5*xp^2)*zeta2)/(-1 + xp) - 
         (4*(1 + 3*xp^2)*Li[2, 1 - xp])/(-1 + xp))*Log[xp] + 
       (4*(3 + 4*xp^2)*Log[1 - xp]^2*Log[xp])/(-1 + xp) + 
       (-2 - 9/(-1 + xp) + 12*xp)*Log[xp]^2 + (5*(1 + xp)*Log[xp]^3)/3 + 
       Log[1 - xp]*(8 + 54*xp + 16*(1 + xp)*zeta2 + (8*xp^2*Li[2, 1 - xp])/
          (-1 + xp) + (12*(3 - 2*xp)*xp*Log[xp])/(-1 + xp) - 
         (8*(1 + 2*xp^2)*Log[xp]^2)/(-1 + xp)) - (28*(1 + xp^2)*S12[1 - xp])/
        (-1 + xp)) + Delx*(2 + Delz*(511/4 + 58*zeta2 - (56*zeta2^2)/5 - 
         60*zeta3) + zeta2*(8 - 8*zp) - 2*zp - 16*zeta3*(1 + zp) + 
       32*zeta3*Dzp[0] + (-64 - 32*zeta2)*Dzp[1] + 8*Dzp[3] + 
       12*(-3 + (1 - zp)^(-1))*Li[2, 1 - zp] + (8*(2 + zp^2)*Li[3, 1 - zp])/
        (-1 + zp) - 4*(1 + zp)*Log[1 - zp]^3 + 
       (30 + 76/(-1 + zp) + 58*zp - (8*zeta2*(1 + 2*zp^2))/(-1 + zp) + 
         (16*(1 + 2*zp^2)*Li[2, 1 - zp])/(-1 + zp))*Log[zp] - 
       (4*(2 + zp^2)*Log[1 - zp]^2*Log[zp])/(-1 + zp) + 
       (-21 - 9/(-1 + zp) - 7*zp)*Log[zp]^2 + (5*(3 + 5*zp^2)*Log[zp]^3)/
        (3*(-1 + zp)) + Log[1 - zp]*(10 + 56*zp + 16*zeta2*(1 + zp) - 
         (8*Li[2, 1 - zp])/(-1 + zp) - (12*(-1 + zp + zp^2)*Log[zp])/
          (-1 + zp) + (2*(1 + 5*zp^2)*Log[zp]^2)/(-1 + zp)) + 
       (16*(2 + 3*zp^2)*S12[1 - zp])/(-1 + zp)) + Theta[1 - xp - zp]*
      ((16*zeta2*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp)))/
        ((-1 + xp)*(-1 + zp)) - (8*(-4 + 11*zp + 4*xp^2*zp - 3*zp^2 + 
          xp*(5 - 19*zp + 6*zp^2)))/((-1 + xp)*(-1 + zp)) + 
       (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*
         Li[2, (1 - xp)^(-1) - xp/(1 - xp) - zp/(1 - xp)])/
        ((-1 + xp)*(-1 + zp)) - (16*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + 
          xp^2*(-1 + 2*zp))*Log[xp]^2)/((-1 + xp)*(-1 + zp)) - 
       (16*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*Log[1 - zp]^2)/
        ((-1 + xp)*(-1 + zp)) - (8*(-3 + xp*(2 - 8*zp) + 6*zp + 
          xp^2*(-1 + 4*zp))*Log[1 - xp - zp])/((-1 + xp)*(-1 + zp)) + 
       ((16*(-3 + xp*(2 - 8*zp) + 6*zp + xp^2*(-1 + 4*zp)))/
          ((-1 + xp)*(-1 + zp)) + (16*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + 
            xp^2*(-1 + 2*zp))*Log[1 - xp - zp])/((-1 + xp)*(-1 + zp)))*
        Log[zp] - (16*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*
         Log[zp]^2)/((-1 + xp)*(-1 + zp)) + Log[1 - zp]*
        ((16*(-3 + xp*(2 - 8*zp) + 6*zp + xp^2*(-1 + 4*zp)))/
          ((-1 + xp)*(-1 + zp)) + (16*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + 
            xp^2*(-1 + 2*zp))*Log[1 - xp - zp])/((-1 + xp)*(-1 + zp)) - 
         (32*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*Log[zp])/
          ((-1 + xp)*(-1 + zp))) + Log[1 - xp]*
        ((8*(-3 + xp*(2 - 8*zp) + 6*zp + xp^2*(-1 + 4*zp)))/
          ((-1 + xp)*(-1 + zp)) + (16*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + 
            xp^2*(-1 + 2*zp))*Log[xp])/((-1 + xp)*(-1 + zp)) - 
         (16*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*Log[1 - zp])/
          ((-1 + xp)*(-1 + zp)) - (16*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + 
            xp^2*(-1 + 2*zp))*Log[zp])/((-1 + xp)*(-1 + zp))) + 
       Log[xp]*((-16*(-3 + xp*(2 - 8*zp) + 6*zp + xp^2*(-1 + 4*zp)))/
          ((-1 + xp)*(-1 + zp)) + (32*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + 
            xp^2*(-1 + 2*zp))*Log[1 - zp])/((-1 + xp)*(-1 + zp)) - 
         (16*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*
           Log[1 - xp - zp])/((-1 + xp)*(-1 + zp)) + 
         (32*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*Log[zp])/
          ((-1 + xp)*(-1 + zp))) + 
       ((-8*(3*zp^2 + xp^2*(-5 + 8*zp) + xp*(4 - 4*zp - 6*zp^2)))/
          ((-1 + xp)*(-1 + zp)) - (16*zeta2*(-((-2 + zp)*zp) + 
            xp^2*(-1 + 2*zp) + xp*(2 - 4*zp + 4*zp^2)))/
          ((-1 + xp)*(-1 + zp)) + ((-8 + (8 - 16*xp)*zp^2)*
           Li[2, (1 - zp)^(-1) - zp/(xp*(1 - zp))])/((-1 + xp)*(-1 + zp)) - 
         (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*
           Li[2, (1 - xp)^(-1) - xp/(1 - xp) - zp/((1 - xp)*xp) + 
             zp^2/((1 - xp)*xp)])/((-1 + xp)*(-1 + zp)) + 
         (16*(1 + (-1 + 2*xp)*zp^2)*Log[1 - xp]^2)/((-1 + xp)*(-1 + zp)) + 
         (16*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*Log[1 - zp]^2)/
          ((-1 + xp)*(-1 + zp)) + (8*(-3 + xp*(2 - 8*zp) + 6*zp + 
            xp^2*(-1 + 4*zp))*Log[1 - xp - zp])/((-1 + xp)*(-1 + zp)) + 
         (8*(-((-4 + zp)*zp) - 2*xp*(2 + zp) + xp^2*(1 + 2*zp))*Log[xp - zp])/
          ((-1 + xp)*(-1 + zp)) + 
         ((16*((-4 + zp)*zp + 2*xp*(2 + zp) - xp^2*(1 + 2*zp)))/
            ((-1 + xp)*(-1 + zp)) - (16*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + 
              xp^2*(-1 + 2*zp))*Log[1 - xp - zp])/((-1 + xp)*(-1 + zp)) - 
           (16*(-((-2 + zp)*zp) + xp^2*(-1 + 2*zp) + xp*(2 - 4*zp + 4*zp^2))*
             Log[xp - zp])/((-1 + xp)*(-1 + zp)))*Log[zp] + 
         (16*(-((-2 + zp)*zp) + xp^2*(-1 + 2*zp) + xp*(2 - 4*zp + 4*zp^2))*
           Log[zp]^2)/((-1 + xp)*(-1 + zp)) + 
         Log[xp]*((8*(-((-4 + zp)*zp) - 2*xp*(2 + zp) + xp^2*(1 + 2*zp)))/
            ((-1 + xp)*(-1 + zp)) - (16*(-((-2 + zp)*zp) + xp^2*(-1 + 2*zp) + 
              xp*(2 - 4*zp + 4*zp^2))*Log[1 - zp])/((-1 + xp)*(-1 + zp)) + 
           (16*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*
             Log[1 - xp - zp])/((-1 + xp)*(-1 + zp)) + 
           (16*(-((-2 + zp)*zp) + xp^2*(-1 + 2*zp) + xp*(2 - 4*zp + 4*zp^2))*
             Log[xp - zp])/((-1 + xp)*(-1 + zp)) - 
           (16*(-((-2 + zp)*zp) + xp^2*(-1 + 2*zp) + xp*(2 - 4*zp + 4*zp^2))*
             Log[zp])/((-1 + xp)*(-1 + zp))) + Log[1 - zp]*
          ((8*(3 - 10*zp - 6*xp^2*zp + zp^2 + 2*xp*(1 + 5*zp)))/
            ((-1 + xp)*(-1 + zp)) - (16*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + 
              xp^2*(-1 + 2*zp))*Log[1 - xp - zp])/((-1 + xp)*(-1 + zp)) - 
           (16*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*
             Log[xp - zp])/((-1 + xp)*(-1 + zp)) + 
           (16*(-1 + 4*zp - zp^2 + xp^2*(-2 + 4*zp) + xp*(4 - 8*zp + 6*zp^2))*
             Log[zp])/((-1 + xp)*(-1 + zp))) + Log[1 - xp]*
          ((-8*(3 + 3*xp^2 + 2*zp - 2*zp^2 + 2*xp*(-5 + 2*zp)))/
            ((-1 + xp)*(-1 + zp)) - (16*(-((-2 + zp)*zp) + xp^2*(-1 + 2*zp) + 
              xp*(2 - 4*zp + 4*zp^2))*Log[xp])/((-1 + xp)*(-1 + zp)) + 
           (16*(-((-2 + zp)*zp) + xp^2*(-1 + 2*zp) + xp*(2 - 4*zp + 4*zp^2))*
             Log[1 - zp])/((-1 + xp)*(-1 + zp)) - (16*(1 + (-1 + 2*xp)*zp^2)*
             Log[xp - zp])/((-1 + xp)*(-1 + zp)) + 
           (16*(1 + 2*zp - 2*zp^2 + xp^2*(-1 + 2*zp) + xp*(2 - 4*zp + 
                6*zp^2))*Log[zp])/((-1 + xp)*(-1 + zp))))*Theta[xp - zp] + 
       ((-8*zeta2*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp)))/
          ((-1 + xp)*(-1 + zp)) + (8*(-4 + 11*zp + 4*xp^2*zp - 3*zp^2 + 
            xp*(5 - 19*zp + 6*zp^2)))/((-1 + xp)*(-1 + zp)) + 
         (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*
           Li[2, -(xp/((1 - xp - zp)*(-xp + zp))) + xp^2/((1 - xp - zp)*(
                -xp + zp))])/((-1 + xp)*(-1 + zp)) + 
         ((-4 + 8*xp*(-1 + zp)^2 + 8*zp + xp^2*(-4 + 8*zp))*Log[1 - xp]^2)/
          ((-1 + xp)*(-1 + zp)) + ((-4 + 8*xp*(-1 + zp)^2 + 8*zp + 
            xp^2*(-4 + 8*zp))*Log[xp]^2)/((-1 + xp)*(-1 + zp)) + 
         (16*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*Log[1 - zp]^2)/
          ((-1 + xp)*(-1 + zp)) + ((-4 + 8*xp*(-1 + zp)^2 + 8*zp + 
            xp^2*(-4 + 8*zp))*Log[1 - xp - zp]^2)/((-1 + xp)*(-1 + zp)) + 
         (16*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*Log[zp]^2)/
          ((-1 + xp)*(-1 + zp)) + (8*(-3 + xp*(2 - 8*zp) + 6*zp + 
            xp^2*(-1 + 4*zp))*Log[-xp + zp])/((-1 + xp)*(-1 + zp)) + 
         ((-4 + 8*xp*(-1 + zp)^2 + 8*zp + xp^2*(-4 + 8*zp))*Log[-xp + zp]^2)/
          ((-1 + xp)*(-1 + zp)) + Log[zp]*
          ((-16*(-3 + xp*(2 - 8*zp) + 6*zp + xp^2*(-1 + 4*zp)))/
            ((-1 + xp)*(-1 + zp)) - (16*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + 
              xp^2*(-1 + 2*zp))*Log[1 - xp - zp])/((-1 + xp)*(-1 + zp)) - 
           (16*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*
             Log[-xp + zp])/((-1 + xp)*(-1 + zp))) + 
         Log[1 - zp]*((-16*(-3 + xp*(2 - 8*zp) + 6*zp + xp^2*(-1 + 4*zp)))/
            ((-1 + xp)*(-1 + zp)) - (16*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + 
              xp^2*(-1 + 2*zp))*Log[1 - xp - zp])/((-1 + xp)*(-1 + zp)) + 
           (32*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*Log[zp])/
            ((-1 + xp)*(-1 + zp)) - (16*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + 
              xp^2*(-1 + 2*zp))*Log[-xp + zp])/((-1 + xp)*(-1 + zp))) + 
         Log[1 - xp]*((-8*(-3 + xp*(2 - 8*zp) + 6*zp + xp^2*(-1 + 4*zp)))/
            ((-1 + xp)*(-1 + zp)) - (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + 
              xp^2*(-1 + 2*zp))*Log[xp])/((-1 + xp)*(-1 + zp)) + 
           (16*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*Log[1 - zp])/
            ((-1 + xp)*(-1 + zp)) - (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + 
              xp^2*(-1 + 2*zp))*Log[1 - xp - zp])/((-1 + xp)*(-1 + zp)) + 
           (16*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*Log[zp])/
            ((-1 + xp)*(-1 + zp)) - (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + 
              xp^2*(-1 + 2*zp))*Log[-xp + zp])/((-1 + xp)*(-1 + zp))) + 
         Log[1 - xp - zp]*((8*(-3 + xp*(2 - 8*zp) + 6*zp + xp^2*(-1 + 4*zp)))/
            ((-1 + xp)*(-1 + zp)) + (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + 
              xp^2*(-1 + 2*zp))*Log[-xp + zp])/((-1 + xp)*(-1 + zp))) + 
         Log[xp]*((8*(-3 + xp*(2 - 8*zp) + 6*zp + xp^2*(-1 + 4*zp)))/
            ((-1 + xp)*(-1 + zp)) - (16*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + 
              xp^2*(-1 + 2*zp))*Log[1 - zp])/((-1 + xp)*(-1 + zp)) + 
           (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*
             Log[1 - xp - zp])/((-1 + xp)*(-1 + zp)) - 
           (16*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*Log[zp])/
            ((-1 + xp)*(-1 + zp)) + (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + 
              xp^2*(-1 + 2*zp))*Log[-xp + zp])/((-1 + xp)*(-1 + zp))))*
        Theta[-xp + zp]) + 
     ((8*zeta2*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp)))/
        ((-1 + xp)*(-1 + zp)) - (8*(-4 + 11*zp + 4*xp^2*zp - 3*zp^2 + 
          xp*(5 - 19*zp + 6*zp^2)))/((-1 + xp)*(-1 + zp)) - 
       (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*
         Li[2, (1 - xp - zp)^(-1) - xp/(1 - xp - zp)])/
        ((-1 + xp)*(-1 + zp)) + ((4 + xp^2*(4 - 8*zp) - 8*xp*(-1 + zp)^2 - 
          8*zp)*Log[1 - xp]^2)/((-1 + xp)*(-1 + zp)) - 
       (16*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*Log[xp]^2)/
        ((-1 + xp)*(-1 + zp)) - (16*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + 
          xp^2*(-1 + 2*zp))*Log[1 - zp]^2)/((-1 + xp)*(-1 + zp)) - 
       (16*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*Log[zp]^2)/
        ((-1 + xp)*(-1 + zp)) - (8*(-3 + xp*(2 - 8*zp) + 6*zp + 
          xp^2*(-1 + 4*zp))*Log[-1 + xp + zp])/((-1 + xp)*(-1 + zp)) + 
       ((4 + xp^2*(4 - 8*zp) - 8*xp*(-1 + zp)^2 - 8*zp)*Log[-1 + xp + zp]^2)/
        ((-1 + xp)*(-1 + zp)) + Log[xp]*
        ((-16*(-3 + xp*(2 - 8*zp) + 6*zp + xp^2*(-1 + 4*zp)))/
          ((-1 + xp)*(-1 + zp)) + (32*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + 
            xp^2*(-1 + 2*zp))*Log[1 - zp])/((-1 + xp)*(-1 + zp)) + 
         (32*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*Log[zp])/
          ((-1 + xp)*(-1 + zp)) - (16*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + 
            xp^2*(-1 + 2*zp))*Log[-1 + xp + zp])/((-1 + xp)*(-1 + zp))) + 
       Log[1 - xp]*((8*(-3 + xp*(2 - 8*zp) + 6*zp + xp^2*(-1 + 4*zp)))/
          ((-1 + xp)*(-1 + zp)) + (16*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + 
            xp^2*(-1 + 2*zp))*Log[xp])/((-1 + xp)*(-1 + zp)) - 
         (16*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*Log[1 - zp])/
          ((-1 + xp)*(-1 + zp)) - (16*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + 
            xp^2*(-1 + 2*zp))*Log[zp])/((-1 + xp)*(-1 + zp)) + 
         (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*
           Log[-1 + xp + zp])/((-1 + xp)*(-1 + zp))) + 
       Log[zp]*((16*(-3 + xp*(2 - 8*zp) + 6*zp + xp^2*(-1 + 4*zp)))/
          ((-1 + xp)*(-1 + zp)) + (16*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + 
            xp^2*(-1 + 2*zp))*Log[-1 + xp + zp])/((-1 + xp)*(-1 + zp))) + 
       Log[1 - zp]*((16*(-3 + xp*(2 - 8*zp) + 6*zp + xp^2*(-1 + 4*zp)))/
          ((-1 + xp)*(-1 + zp)) - (32*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + 
            xp^2*(-1 + 2*zp))*Log[zp])/((-1 + xp)*(-1 + zp)) + 
         (16*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*
           Log[-1 + xp + zp])/((-1 + xp)*(-1 + zp))))*Theta[-1 + xp + zp] + 
     Theta[xp - zp]*((8*(3*zp^2 + xp^2*(-5 + 8*zp) + xp*(4 - 4*zp - 6*zp^2)))/
        ((-1 + xp)*(-1 + zp)) + (8*zeta2*(1 + 2*zp - 2*zp^2 + 
          xp^2*(-1 + 2*zp) + xp*(2 - 4*zp + 6*zp^2)))/((-1 + xp)*(-1 + zp)) - 
       (8*(-1 + 4*zp - 2*zp^2 + xp^2*(-3 + 4*zp) + xp*(4 - 8*zp + 8*zp^2))*
         Li[2, zp/xp])/((-1 + xp)*(-1 + zp)) + 
       ((8 + 8*(-1 + 2*xp)*zp^2)*Li[2, xp/(1 - zp) - zp/(1 - zp)])/
        ((-1 + xp)*(-1 + zp)) - (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + 
          xp^2*(-1 + 2*zp))*Li[2, (-xp + zp)^(-1) - xp/(-xp + zp)])/
        ((-1 + xp)*(-1 + zp)) - (4*(3 + 2*zp - 4*zp^2 + xp^2*(-1 + 2*zp) + 
          2*xp*(1 - 2*zp + 5*zp^2))*Log[1 - xp]^2)/((-1 + xp)*(-1 + zp)) - 
       (8*(-1 + 8*zp - 4*zp^2 + xp^2*(-5 + 8*zp) + 8*xp*(1 - 2*zp + 2*zp^2))*
         Log[xp]^2)/((-1 + xp)*(-1 + zp)) - 
       (16*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*Log[1 - zp]^2)/
        ((-1 + xp)*(-1 + zp)) + ((8*(-4 + zp)*zp + 16*xp*(2 + zp) - 
          8*xp^2*(1 + 2*zp))*Log[xp - zp])/((-1 + xp)*(-1 + zp)) + 
       ((4 + xp^2*(4 - 8*zp) - 8*xp*(-1 + zp)^2 - 8*zp)*Log[xp - zp]^2)/
        ((-1 + xp)*(-1 + zp)) + 
       ((16*(-((-4 + zp)*zp) - 2*xp*(2 + zp) + xp^2*(1 + 2*zp)))/
          ((-1 + xp)*(-1 + zp)) + (8*(1 + xp^2)*Log[xp - zp])/
          ((-1 + xp)*(-1 + zp)))*Log[zp] - 
       (16*(-((-2 + zp)*zp) + xp^2*(-1 + 2*zp) + xp*(2 - 4*zp + 4*zp^2))*
         Log[zp]^2)/((-1 + xp)*(-1 + zp)) + 
       Log[xp]*((16*((-4 + zp)*zp + 2*xp*(2 + zp) - xp^2*(1 + 2*zp)))/
          ((-1 + xp)*(-1 + zp)) + (16*(-1 + 4*zp - zp^2 + xp^2*(-2 + 4*zp) + 
            xp*(4 - 8*zp + 6*zp^2))*Log[1 - zp])/((-1 + xp)*(-1 + zp)) - 
         (8*(1 + xp^2)*Log[xp - zp])/((-1 + xp)*(-1 + zp)) + 
         (8*(-1 + 12*zp - 6*zp^2 + xp^2*(-7 + 12*zp) + 
            12*xp*(1 - 2*zp + 2*zp^2))*Log[zp])/((-1 + xp)*(-1 + zp))) + 
       Log[1 - zp]*((-8*(3 - 10*zp - 6*xp^2*zp + zp^2 + 2*xp*(1 + 5*zp)))/
          ((-1 + xp)*(-1 + zp)) + (16*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + 
            xp^2*(-1 + 2*zp))*Log[xp - zp])/((-1 + xp)*(-1 + zp)) - 
         (16*(-1 + 4*zp - zp^2 + xp^2*(-2 + 4*zp) + xp*(4 - 8*zp + 6*zp^2))*
           Log[zp])/((-1 + xp)*(-1 + zp))) + Log[1 - xp]*
        ((8*(3 + 3*xp^2 + 2*zp - 2*zp^2 + 2*xp*(-5 + 2*zp)))/
          ((-1 + xp)*(-1 + zp)) + (16*(1 + 2*zp - 2*zp^2 + xp^2*(-1 + 2*zp) + 
            xp*(2 - 4*zp + 6*zp^2))*Log[xp])/((-1 + xp)*(-1 + zp)) - 
         (16*(-((-2 + zp)*zp) + xp^2*(-1 + 2*zp) + xp*(2 - 4*zp + 4*zp^2))*
           Log[1 - zp])/((-1 + xp)*(-1 + zp)) + 
         (8*(1 + 2*zp - 2*zp^2 + xp^2*(-1 + 2*zp) + xp*(2 - 4*zp + 6*zp^2))*
           Log[xp - zp])/((-1 + xp)*(-1 + zp)) - 
         (16*(1 + 2*zp - 2*zp^2 + xp^2*(-1 + 2*zp) + xp*(2 - 4*zp + 6*zp^2))*
           Log[zp])/((-1 + xp)*(-1 + zp))) + 
       ((-8*(3*zp^2 + xp^2*(-5 + 8*zp) + xp*(4 - 4*zp - 6*zp^2)))/
          ((-1 + xp)*(-1 + zp)) - (8*zeta2*(1 + 2*zp - 2*zp^2 + 
            xp^2*(-1 + 2*zp) + xp*(2 - 4*zp + 6*zp^2)))/
          ((-1 + xp)*(-1 + zp)) + ((-8 + (8 - 16*xp)*zp^2)*
           Li[2, (1 - zp)^(-1) - zp/(xp*(1 - zp))])/((-1 + xp)*(-1 + zp)) + 
         (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*
           Li[2, -(xp/((1 - xp - zp)*(-xp + zp))) + xp^2/((1 - xp - zp)*(
                -xp + zp))])/((-1 + xp)*(-1 + zp)) + 
         (4*(3 + 2*zp - 4*zp^2 + xp^2*(-1 + 2*zp) + 2*xp*(1 - 2*zp + 5*zp^2))*
           Log[1 - xp]^2)/((-1 + xp)*(-1 + zp)) + 
         ((-4 + 8*xp*(-1 + zp)^2 + 8*zp + xp^2*(-4 + 8*zp))*Log[xp]^2)/
          ((-1 + xp)*(-1 + zp)) + (16*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + 
            xp^2*(-1 + 2*zp))*Log[1 - zp]^2)/((-1 + xp)*(-1 + zp)) + 
         ((-4 + 8*xp*(-1 + zp)^2 + 8*zp + xp^2*(-4 + 8*zp))*Log[xp - zp]^2)/
          ((-1 + xp)*(-1 + zp)) + (16*(-((-2 + zp)*zp) + xp^2*(-1 + 2*zp) + 
            xp*(2 - 4*zp + 4*zp^2))*Log[zp]^2)/((-1 + xp)*(-1 + zp)) + 
         (8*(-3 + xp*(2 - 8*zp) + 6*zp + xp^2*(-1 + 4*zp))*Log[-1 + xp + zp])/
          ((-1 + xp)*(-1 + zp)) + ((-4 + 8*xp*(-1 + zp)^2 + 8*zp + 
            xp^2*(-4 + 8*zp))*Log[-1 + xp + zp]^2)/((-1 + xp)*(-1 + zp)) + 
         Log[zp]*((16*((-4 + zp)*zp + 2*xp*(2 + zp) - xp^2*(1 + 2*zp)))/
            ((-1 + xp)*(-1 + zp)) - (16*(-((-2 + zp)*zp) + xp^2*(-1 + 2*zp) + 
              xp*(2 - 4*zp + 4*zp^2))*Log[xp - zp])/((-1 + xp)*(-1 + zp)) - 
           (16*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*
             Log[-1 + xp + zp])/((-1 + xp)*(-1 + zp))) + 
         Log[1 - zp]*((8*(3 - 10*zp - 6*xp^2*zp + zp^2 + 2*xp*(1 + 5*zp)))/
            ((-1 + xp)*(-1 + zp)) - (16*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + 
              xp^2*(-1 + 2*zp))*Log[xp - zp])/((-1 + xp)*(-1 + zp)) + 
           (16*(-1 + 4*zp - zp^2 + xp^2*(-2 + 4*zp) + xp*(4 - 8*zp + 6*zp^2))*
             Log[zp])/((-1 + xp)*(-1 + zp)) - 
           (16*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*
             Log[-1 + xp + zp])/((-1 + xp)*(-1 + zp))) + 
         Log[1 - xp]*((-8*(3 + 3*xp^2 + 2*zp - 2*zp^2 + 2*xp*(-5 + 2*zp)))/
            ((-1 + xp)*(-1 + zp)) - (8*(1 + 2*zp - 2*zp^2 + xp^2*(-1 + 
                2*zp) + xp*(2 - 4*zp + 6*zp^2))*Log[xp])/
            ((-1 + xp)*(-1 + zp)) + (16*(-((-2 + zp)*zp) + xp^2*(-1 + 2*zp) + 
              xp*(2 - 4*zp + 4*zp^2))*Log[1 - zp])/((-1 + xp)*(-1 + zp)) - 
           (8*(1 + 2*zp - 2*zp^2 + xp^2*(-1 + 2*zp) + xp*(2 - 4*zp + 6*zp^2))*
             Log[xp - zp])/((-1 + xp)*(-1 + zp)) + 
           (16*(1 + 2*zp - 2*zp^2 + xp^2*(-1 + 2*zp) + xp*(2 - 4*zp + 
                6*zp^2))*Log[zp])/((-1 + xp)*(-1 + zp)) - 
           (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*
             Log[-1 + xp + zp])/((-1 + xp)*(-1 + zp))) + 
         Log[xp - zp]*((8*(-((-4 + zp)*zp) - 2*xp*(2 + zp) + 
              xp^2*(1 + 2*zp)))/((-1 + xp)*(-1 + zp)) + 
           (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*
             Log[-1 + xp + zp])/((-1 + xp)*(-1 + zp))) + 
         Log[xp]*((8*(-((-4 + zp)*zp) - 2*xp*(2 + zp) + xp^2*(1 + 2*zp)))/
            ((-1 + xp)*(-1 + zp)) - (16*(-((-2 + zp)*zp) + xp^2*(-1 + 2*zp) + 
              xp*(2 - 4*zp + 4*zp^2))*Log[1 - zp])/((-1 + xp)*(-1 + zp)) + 
           (8*(1 + 2*zp - 2*zp^2 + xp^2*(-1 + 2*zp) + xp*(2 - 4*zp + 6*zp^2))*
             Log[xp - zp])/((-1 + xp)*(-1 + zp)) - 
           (16*(-((-2 + zp)*zp) + xp^2*(-1 + 2*zp) + xp*(2 - 4*zp + 4*zp^2))*
             Log[zp])/((-1 + xp)*(-1 + zp)) + 
           (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*
             Log[-1 + xp + zp])/((-1 + xp)*(-1 + zp))))*
        Theta[-1 + xp + zp]) + Theta[-xp + zp]*
      ((-32*zeta2*(xp^2*(-1 + zp) - (-1 + zp)*zp + xp*(1 - 2*zp + 3*zp^2)))/
        ((-1 + xp)*(-1 + zp)) - (8*(-4 + 11*zp + 4*xp^2*zp - 3*zp^2 + 
          xp*(5 - 19*zp + 6*zp^2)))/((-1 + xp)*(-1 + zp)) + 
       (8*(-1 + 4*zp - 2*zp^2 + xp^2*(-3 + 4*zp) + xp*(4 - 8*zp + 8*zp^2))*
         Li[2, xp/zp])/((-1 + xp)*(-1 + zp)) + 
       (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*
         Li[2, -(xp/(1 - xp)) + zp/(1 - xp)])/((-1 + xp)*(-1 + zp)) + 
       ((-8 + (8 - 16*xp)*zp^2)*Li[2, -(-xp + zp)^(-1) + zp/(-xp + zp)])/
        ((-1 + xp)*(-1 + zp)) + ((8 + 8*(-1 + 2*xp)*zp^2)*
         Li[2, -(xp/(-xp + zp)) + (xp*zp)/(-xp + zp)])/
        ((-1 + xp)*(-1 + zp)) + ((8 + xp^2*(28 - 48*zp) - 48*zp + 20*zp^2 - 
          8*xp*(6 - 12*zp + 11*zp^2))*Log[xp]^2)/((-1 + xp)*(-1 + zp)) - 
       (16*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*Log[1 - zp]^2)/
        ((-1 + xp)*(-1 + zp)) - (4*(-3 + xp*(4 - 8*zp) + 4*zp + 2*zp^2 + 
          xp^2*(-1 + 4*zp))*Log[zp]^2)/((-1 + xp)*(-1 + zp)) + 
       Log[1 - xp]*((8*(-3 + xp*(2 - 8*zp) + 6*zp + xp^2*(-1 + 4*zp)))/
          ((-1 + xp)*(-1 + zp)) + (16*(-((-2 + zp)*zp) + xp^2*(-1 + 2*zp) + 
            xp*(2 - 4*zp + 4*zp^2))*Log[xp])/((-1 + xp)*(-1 + zp)) - 
         (16*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*Log[1 - zp])/
          ((-1 + xp)*(-1 + zp)) - (16*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + 
            xp^2*(-1 + 2*zp))*Log[zp])/((-1 + xp)*(-1 + zp))) - 
       (8*(-3 + xp*(2 - 8*zp) + 6*zp + xp^2*(-1 + 4*zp))*Log[-xp + zp])/
        ((-1 + xp)*(-1 + zp)) + Log[zp]*
        ((16*(-3 + xp*(2 - 8*zp) + 6*zp + xp^2*(-1 + 4*zp)))/
          ((-1 + xp)*(-1 + zp)) + (8*(-1 + xp^2 + 2*zp^2 - 4*xp*zp^2)*
           Log[-xp + zp])/((-1 + xp)*(-1 + zp))) + 
       Log[xp]*((8*(3 - 10*zp - 6*xp^2*zp + zp^2 + 2*xp*(1 + 5*zp)))/
          ((-1 + xp)*(-1 + zp)) + (8*(-3 + 8*zp - zp^2 + xp^2*(-4 + 8*zp) + 
            2*xp*(4 - 8*zp + 5*zp^2))*Log[1 - zp])/((-1 + xp)*(-1 + zp)) + 
         (16*(-1 + 4*zp - zp^2 + xp^2*(-2 + 4*zp) + xp*(4 - 8*zp + 6*zp^2))*
           Log[zp])/((-1 + xp)*(-1 + zp)) - (8*(xp^2 + zp^2 - 2*xp*zp^2)*
           Log[-xp + zp])/((-1 + xp)*(-1 + zp))) + 
       Log[1 - zp]*((16*(-3 + xp*(2 - 8*zp) + 6*zp + xp^2*(-1 + 4*zp)))/
          ((-1 + xp)*(-1 + zp)) - (32*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + 
            xp^2*(-1 + 2*zp))*Log[zp])/((-1 + xp)*(-1 + zp)) + 
         (16*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*Log[-xp + zp])/
          ((-1 + xp)*(-1 + zp))) + 
       ((-16*zeta2*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp)))/
          ((-1 + xp)*(-1 + zp)) + (8*(-4 + 11*zp + 4*xp^2*zp - 3*zp^2 + 
            xp*(5 - 19*zp + 6*zp^2)))/((-1 + xp)*(-1 + zp)) - 
         (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*
           Li[2, (1 - xp)^(-1) - xp/(1 - xp) - zp/((1 - xp)*xp) + 
             zp^2/((1 - xp)*xp)])/((-1 + xp)*(-1 + zp)) + 
         (16*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*Log[1 - zp]^2)/
          ((-1 + xp)*(-1 + zp)) + (16*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + 
            xp^2*(-1 + 2*zp))*Log[zp]^2)/((-1 + xp)*(-1 + zp)) + 
         Log[1 - xp]*((-8*(-3 + xp*(2 - 8*zp) + 6*zp + xp^2*(-1 + 4*zp)))/
            ((-1 + xp)*(-1 + zp)) - (16*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + 
              xp^2*(-1 + 2*zp))*Log[xp])/((-1 + xp)*(-1 + zp)) + 
           (16*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*Log[1 - zp])/
            ((-1 + xp)*(-1 + zp)) + (16*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + 
              xp^2*(-1 + 2*zp))*Log[zp])/((-1 + xp)*(-1 + zp))) + 
         (8*(-3 + xp*(2 - 8*zp) + 6*zp + xp^2*(-1 + 4*zp))*Log[-xp + zp])/
          ((-1 + xp)*(-1 + zp)) + (8*(-3 + xp*(2 - 8*zp) + 6*zp + 
            xp^2*(-1 + 4*zp))*Log[-1 + xp + zp])/((-1 + xp)*(-1 + zp)) + 
         Log[zp]*((-16*(-3 + xp*(2 - 8*zp) + 6*zp + xp^2*(-1 + 4*zp)))/
            ((-1 + xp)*(-1 + zp)) - (16*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + 
              xp^2*(-1 + 2*zp))*Log[-xp + zp])/((-1 + xp)*(-1 + zp)) - 
           (16*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*
             Log[-1 + xp + zp])/((-1 + xp)*(-1 + zp))) + 
         Log[1 - zp]*((-16*(-3 + xp*(2 - 8*zp) + 6*zp + xp^2*(-1 + 4*zp)))/
            ((-1 + xp)*(-1 + zp)) + (32*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + 
              xp^2*(-1 + 2*zp))*Log[zp])/((-1 + xp)*(-1 + zp)) - 
           (16*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*
             Log[-xp + zp])/((-1 + xp)*(-1 + zp)) - 
           (16*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*
             Log[-1 + xp + zp])/((-1 + xp)*(-1 + zp))) + 
         Log[xp]*((8*(-3 + xp*(2 - 8*zp) + 6*zp + xp^2*(-1 + 4*zp)))/
            ((-1 + xp)*(-1 + zp)) - (16*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + 
              xp^2*(-1 + 2*zp))*Log[1 - zp])/((-1 + xp)*(-1 + zp)) - 
           (16*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*Log[zp])/
            ((-1 + xp)*(-1 + zp)) + (16*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + 
              xp^2*(-1 + 2*zp))*Log[-xp + zp])/((-1 + xp)*(-1 + zp)) + 
           (16*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*
             Log[-1 + xp + zp])/((-1 + xp)*(-1 + zp))))*
        Theta[-1 + xp + zp])) + CA*CF*(-56/3 + 4/(-1 + xp) + 16/(-1 + zp) + 
     24/zp + (52/9 - 12/(-1 + xp))*zp - 4/(1 + 2*xp + xp^2 - 4*xp*zp) + 
     (4*xp^3)/(1 + 2*xp + xp^2 - 4*xp*zp) - 
     (4*zeta2*(2 + (-2 + 4*xp + xp^2)*zp + (4 - 6*xp)*zp^2 + 
        (-3 + 6*xp)*zp^3))/((-1 + xp)*(-1 + zp)*zp) + 
     xp^2*(8/(-xp + zp) - 4/(1 + 2*xp + xp^2 - 4*xp*zp)) + 
     xp*(688/9 + 4/(xp - zp) - 16/(-1 + zp) - 22/zp + 6*zp + 
       4/(1 + 2*xp + xp^2 - 4*xp*zp)) + 
     Dxp[1]*((22*(1 + zp))/3 - (44*Dzp[0])/3) + (22*(1 + xp)*Dzp[1])/3 + 
     (4*(3*zp + 3*xp^2*zp + xp*(-1 + 3*zp^2))*
       InvTanInt[-(Sqrt[zp]/Sqrt[xp])])/(Sqrt[xp]*zp^(3/2)) - 
     (4*(3*zp + 3*xp^2*zp + xp*(-1 + 3*zp^2))*InvTanInt[Sqrt[zp]/Sqrt[xp]])/
      (Sqrt[xp]*zp^(3/2)) - (4*(3*zp + 3*xp^2*zp + xp*(-1 + 3*zp^2))*
       InvTanInt[-(Sqrt[xp]*Sqrt[zp])])/(Sqrt[xp]*zp^(3/2)) + 
     (4*(3*zp + 3*xp^2*zp + xp*(-1 + 3*zp^2))*InvTanInt[Sqrt[xp]*Sqrt[zp]])/
      (Sqrt[xp]*zp^(3/2)) + (4*(-1 - 2*zp + 2*zp^2 + zp^3 - 
        2*xp*zp^2*(1 + zp) + xp^2*(-1 - zp + 2*zp^2))*Li[2, 1 - xp])/
      ((-1 + xp)*(-1 + zp)*zp) + (8*(1 + xp)*(-1 + 2*zp)*Li[2, -xp])/zp + 
     (4*(-1 + 4*zp - zp^2 + xp^2*(-2 + 4*zp) + xp*(4 - 8*zp + 8*zp^2))*
       Li[2, 1 - zp])/((-1 + xp)*(-1 + zp)) + 
     (8*(1 + 2*xp*zp + 2*xp^3*zp + xp^2*(1 + 2*zp^2))*Li[2, -(zp/xp)])/
      ((-1 + xp)*(1 + xp)*zp) + (8*(xp + xp^3 + 2*zp + 2*xp^2*zp + 2*xp*zp^2)*
       Li[2, -(xp*zp)])/((-1 + xp)*(1 + xp)*zp) + 
     (4*((1 - 2*zp)^2 + xp^6*(1 - 2*zp)^2 + xp*(2 - 20*zp + 48*zp^2 - 
          32*zp^3) + xp^5*(2 - 20*zp + 48*zp^2 - 32*zp^3) + 
        xp^2*(-1 - 16*zp + 104*zp^2 - 176*zp^3 + 88*zp^4) + 
        xp^4*(-1 - 16*zp + 104*zp^2 - 176*zp^3 + 88*zp^4) - 
        4*xp^3*(1 - 12*zp^2 + 28*zp^3 - 30*zp^4 + 12*zp^5))*
       Li[2, 1/2 - xp/2 - Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]/2])/
      ((1 + xp^2 + xp*(2 - 4*zp))^(5/2)*(-1 + zp)*zp) + 
     (4*((1 - 2*zp)^2 + xp^6*(1 - 2*zp)^2 + xp*(2 - 20*zp + 48*zp^2 - 
          32*zp^3) + xp^5*(2 - 20*zp + 48*zp^2 - 32*zp^3) + 
        xp^2*(-1 - 16*zp + 104*zp^2 - 176*zp^3 + 88*zp^4) + 
        xp^4*(-1 - 16*zp + 104*zp^2 - 176*zp^3 + 88*zp^4) - 
        4*xp^3*(1 - 12*zp^2 + 28*zp^3 - 30*zp^4 + 12*zp^5))*
       Li[2, 1/2 + xp/2 - Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]/2])/
      ((1 + xp^2 + xp*(2 - 4*zp))^(5/2)*(-1 + zp)*zp) - 
     (4*((1 - 2*zp)^2 + xp^6*(1 - 2*zp)^2 + xp*(2 - 20*zp + 48*zp^2 - 
          32*zp^3) + xp^5*(2 - 20*zp + 48*zp^2 - 32*zp^3) + 
        xp^2*(-1 - 16*zp + 104*zp^2 - 176*zp^3 + 88*zp^4) + 
        xp^4*(-1 - 16*zp + 104*zp^2 - 176*zp^3 + 88*zp^4) - 
        4*xp^3*(1 - 12*zp^2 + 28*zp^3 - 30*zp^4 + 12*zp^5))*
       Li[2, 1/2 - 1/(2*xp) - Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]/(2*xp)])/
      ((1 + xp^2 + xp*(2 - 4*zp))^(5/2)*(-1 + zp)*zp) - 
     (4*((1 - 2*zp)^2 + xp^6*(1 - 2*zp)^2 + xp*(2 - 20*zp + 48*zp^2 - 
          32*zp^3) + xp^5*(2 - 20*zp + 48*zp^2 - 32*zp^3) + 
        xp^2*(-1 - 16*zp + 104*zp^2 - 176*zp^3 + 88*zp^4) + 
        xp^4*(-1 - 16*zp + 104*zp^2 - 176*zp^3 + 88*zp^4) - 
        4*xp^3*(1 - 12*zp^2 + 28*zp^3 - 30*zp^4 + 12*zp^5))*
       Li[2, 1/2 + 1/(2*xp) - Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]/(2*xp)])/
      ((1 + xp^2 + xp*(2 - 4*zp))^(5/2)*(-1 + zp)*zp) + 
     ((4 + 8*zp + 8*xp*zp^2 + xp^2*(4 + 8*zp))*
       Li[2, (2 - 2*xp)^(-1) + zp/(2 - 2*xp) - 
         Sqrt[1 - 2*zp + 4*xp*zp + zp^2]/(2 - 2*xp)])/((-1 + xp)*zp) + 
     ((4 + 8*zp + 8*xp*zp^2 + xp^2*(4 + 8*zp))*
       Li[2, (2 - 2*xp)^(-1) - (2*xp)/(2 - 2*xp) - zp/(2 - 2*xp) + 
         Sqrt[1 - 2*zp + 4*xp*zp + zp^2]/(2 - 2*xp)])/(zp - xp*zp) + 
     (8*(1 + 2*zp + 2*xp*zp^2 + xp^2*(1 + 2*zp))*
       Li[2, 1 - 1/(2*xp) + zp/(2*xp) - Sqrt[1 - 2*zp + 4*xp*zp + zp^2]/
          (2*xp)])/((-1 + xp)*zp) + ((2 + 4*zp + 4*xp*zp^2 + xp^2*(2 + 4*zp))*
       Log[2]^2)/((-1 + xp)*zp) + (Delz*((11*(1 + xp))/3 - (22*Dxp[0])/3) + 
       Delx*(-11*Delz + (11*(1 + zp))/3 - (22*Dzp[0])/3))*Log[Q2/muF2]^2 + 
     ((2 + xp^2*(2 - 4*zp) - 4*xp*(-1 + zp)^2 - 4*zp)*Log[1 - xp]^2)/
      ((-1 + xp)*(-1 + zp)) + (-24 - 38/(-1 + xp) - 
       3/(1 + xp^2 + xp*(2 - 4*zp))^(5/2) - 
       (3*xp^6)/(1 + xp^2 + xp*(2 - 4*zp))^(5/2) + 
       3/(1 + xp^2 + xp*(2 - 4*zp))^(3/2) - 
       16/Sqrt[1 + xp^2 + xp*(2 - 4*zp)] + 
       (3*xp^4*(2 + xp^2 + xp*(2 - 4*zp)))/(1 + xp^2 + xp*(2 - 4*zp))^(5/2) + 
       (xp^2*(5 + 2*xp^2 - 16*(1 + xp^2 + xp*(2 - 4*zp))^2 + 
          2*xp*(2 - 4*zp) - (4*(1 + xp^2 + xp*(2 - 4*zp))^2)/(-1 + zp) + 
          (4*(1 + xp^2 + xp*(2 - 4*zp))^2)/zp))/(1 + xp^2 + xp*(2 - 4*zp))^
         (5/2) - 38/(-1 + zp) - 52/((-1 + xp)*(-1 + zp)) - 
       4/(Sqrt[1 + xp^2 + xp*(2 - 4*zp)]*(-1 + zp)) + 
       ((2 - 44*xp)*zp)/(-1 + xp) - 4/(zp - xp*zp) - 
       (2*xp*(2*(-2 + Sqrt[1 + xp^2 + xp*(2 - 4*zp)]) + 
          (5 - 7*Sqrt[1 + xp^2 + xp*(2 - 4*zp)])*zp + 
          3*(3 + 4*Sqrt[1 + xp^2 + xp*(2 - 4*zp)])*zp^2 - 6*zp^3))/
        (Sqrt[1 + xp^2 + xp*(2 - 4*zp)]*(-1 + zp)*zp) + 
       (4 + 4/Sqrt[1 + 2*xp + xp^2 - 4*xp*zp])/zp)*Log[xp]^2 + 
     Dzp[0]*((4*(10 - 77*xp))/9 + 4*(1 + xp)*zeta2 + 
       (22*(1 + xp)*Log[1 - xp])/3 + ((56 + 32*xp^2)*Log[xp])/(3 - 3*xp) - 
       (2*(1 + xp^2)*Log[xp]^2)/(-1 + xp)) - 
     (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*Log[1 - zp]^2)/
      ((-1 + xp)*(-1 + zp)) - (2*(-1 - 5*zp + 10*zp^2 + 2*zp^3 + 
        2*xp*zp*(4 - 9*zp + 4*zp^2) + xp^2*(-1 - 5*zp + 10*zp^2))*Log[zp]^2)/
      ((-1 + xp)*(-1 + zp)*zp) + Log[1 - zp]*
      (4*(-7 + 4/(-1 + xp) - 6/(-1 + zp)) - (44*zp)/3 + 
       xp*(52/3 + 24/(-1 + zp) - 4/(-1 + xp + zp)) - 
       (12*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*Log[zp])/
        ((-1 + xp)*(-1 + zp))) + Dxp[0]*((4*(10 - 77*zp))/9 + 
       4*zeta2*(1 + zp) + (268/9 - 8*zeta2)*Dzp[0] - (44*Dzp[1])/3 + 
       (22*(1 + zp)*Log[1 - zp])/3 + 4*(1 + zp)*Log[zp] - 
       (2*(1 + zp^2)*Log[zp]^2)/(-1 + zp)) + Log[Q2/muF2]*
      ((-44*(xp + zp))/3 + Dxp[0]*((22*(1 + zp))/3 - (44*Dzp[0])/3) + 
       (22*(1 + xp)*Dzp[0])/3 + Delz*((4*(10 - 77*xp))/9 + 4*(1 + xp)*zeta2 + 
         (268/9 - 8*zeta2)*Dxp[0] - (44*Dxp[1])/3 + (22*(1 + xp)*Log[1 - xp])/
          3 + ((56 + 32*xp^2)*Log[xp])/(3 - 3*xp) - (2*(1 + xp^2)*Log[xp]^2)/
          (-1 + xp)) + Delx*(Delz*(193/3 + (88*zeta2)/3 - 24*zeta3) + 
         (4*(10 - 77*zp))/9 + 4*zeta2*(1 + zp) + (268/9 - 8*zeta2)*Dzp[0] - 
         (44*Dzp[1])/3 + (22*(1 + zp)*Log[1 - zp])/3 + 4*(1 + zp)*Log[zp] - 
         (2*(1 + zp^2)*Log[zp]^2)/(-1 + zp))) + 
     Log[zp]*(2*(-29 - 4/(-1 + xp) - 3/(1 + xp^2 + xp*(2 - 4*zp))^2 + 
         (3*xp^4)/(1 + xp^2 + xp*(2 - 4*zp))^2 + 
         (3*xp^5)/(1 + xp^2 + xp*(2 - 4*zp))^2 - 14/(-1 + zp) - 
         (2*xp^2*(-1 + 2*zp))/(xp - zp)^2 + 2/(1 + 2*xp + xp^2 - 4*xp*zp) - 
         (2*xp^3)/(1 + 2*xp + xp^2 - 4*xp*zp) + 6/(zp - xp*zp) + 
         xp*(21 + 14/(-1 + zp) - 6/zp - (6*zp)/(-1 + xp) - 
           3/(1 + 2*xp + xp^2 - 4*xp*zp)^2)) + 
       (4*(3*zp + 3*xp^2*zp + xp*(-1 + 3*zp^2))*ArcTan[Sqrt[zp]/Sqrt[xp]])/
        (Sqrt[xp]*zp^(3/2)) - (4*(3*zp + 3*xp^2*zp + xp*(-1 + 3*zp^2))*
         ArcTan[Sqrt[xp]*Sqrt[zp]])/(Sqrt[xp]*zp^(3/2)) + 
       (8*(1 + 2*xp*zp + 2*xp^3*zp + xp^2*(1 + 2*zp^2))*Log[1 + zp/xp])/
        ((-1 + xp)*(1 + xp)*zp) + (8*(xp + xp^3 + 2*zp + 2*xp^2*zp + 
          2*xp*zp^2)*Log[1 + xp*zp])/((-1 + xp)*(1 + xp)*zp)) + 
     (4*((1 - 2*zp)^2 + xp^6*(1 - 2*zp)^2 + xp*(2 - 20*zp + 48*zp^2 - 
          32*zp^3) + xp^5*(2 - 20*zp + 48*zp^2 - 32*zp^3) + 
        xp^2*(-1 - 16*zp + 104*zp^2 - 176*zp^3 + 88*zp^4) + 
        xp^4*(-1 - 16*zp + 104*zp^2 - 176*zp^3 + 88*zp^4) - 
        4*xp^3*(1 - 12*zp^2 + 28*zp^3 - 30*zp^4 + 12*zp^5))*
       Log[1 + xp - Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]]*
       Log[1 - xp + Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]])/
      ((1 + xp^2 + xp*(2 - 4*zp))^(5/2)*(-1 + zp)*zp) - 
     (4*((1 - 2*zp)^2 + xp^6*(1 - 2*zp)^2 + xp*(2 - 20*zp + 48*zp^2 - 
          32*zp^3) + xp^5*(2 - 20*zp + 48*zp^2 - 32*zp^3) + 
        xp^2*(-1 - 16*zp + 104*zp^2 - 176*zp^3 + 88*zp^4) + 
        xp^4*(-1 - 16*zp + 104*zp^2 - 176*zp^3 + 88*zp^4) - 
        4*xp^3*(1 - 12*zp^2 + 28*zp^3 - 30*zp^4 + 12*zp^5))*
       Log[1 + xp - Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]]*
       Log[-1 + xp + Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]])/
      ((1 + xp^2 + xp*(2 - 4*zp))^(5/2)*(-1 + zp)*zp) + 
     ((2 + 4*zp + 4*xp*zp^2 + xp^2*(2 + 4*zp))*
       Log[1 + zp - Sqrt[1 - 2*zp + 4*xp*zp + zp^2]]^2)/((-1 + xp)*zp) + 
     ((2 + 4*zp + 4*xp*zp^2 + xp^2*(2 + 4*zp))*
       Log[1 - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]]^2)/((-1 + xp)*zp) + 
     Log[1 + zp - Sqrt[1 - 2*zp + 4*xp*zp + zp^2]]*
      ((4*(2*(-1 + zp)^2 + 8*xp^3*zp + 2*xp^2*(1 - 4*zp + 7*zp^2) + 
          xp*(-1 + 13*zp - 7*zp^2 + 3*zp^3)))/((-1 + xp)*zp*
         Sqrt[1 + (-2 + 4*xp)*zp + zp^2]) + 
       ((4 + 8*zp + 8*xp*zp^2 + xp^2*(4 + 8*zp))*
         Log[1 - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/(zp - xp*zp)) + 
     Log[2]*((4*(2*(-1 + zp)^2 + 8*xp^3*zp + 2*xp^2*(1 - 4*zp + 7*zp^2) + 
          xp*(-1 + 13*zp - 7*zp^2 + 3*zp^3)))/((-1 + xp)*zp*
         Sqrt[1 + (-2 + 4*xp)*zp + zp^2]) + 
       2*(4 + 8/(-1 + xp) - 3/(1 + xp^2 + xp*(2 - 4*zp))^(5/2) - 
         (3*xp^6)/(1 + xp^2 + xp*(2 - 4*zp))^(5/2) + 
         3/(1 + xp^2 + xp*(2 - 4*zp))^(3/2) - 
         16/Sqrt[1 + xp^2 + xp*(2 - 4*zp)] + 
         (3*xp^4*(2 + xp^2 + xp*(2 - 4*zp)))/(1 + xp^2 + xp*(2 - 4*zp))^
           (5/2) + (xp^2*(5 + 2*xp^2 - 16*(1 + xp^2 + xp*(2 - 4*zp))^2 + 
            2*xp*(2 - 4*zp) - (4*(1 + xp^2 + xp*(2 - 4*zp))^2)/(-1 + zp) + 
            (4*(1 + xp^2 + xp*(2 - 4*zp))^2)/zp))/(1 + xp^2 + xp*(2 - 4*zp))^
           (5/2) - 4/(Sqrt[1 + xp^2 + xp*(2 - 4*zp)]*(-1 + zp)) + 
         (4*xp*zp)/(-1 + xp) - 4/(zp - xp*zp) + 
         (2 + 4/Sqrt[1 + 2*xp + xp^2 - 4*xp*zp])/zp + 
         2*xp*(2 - 3/Sqrt[1 + xp^2 + xp*(2 - 4*zp)] - 
           4/(Sqrt[1 + xp^2 + xp*(2 - 4*zp)]*(-1 + zp)) + 
           (6*zp)/Sqrt[1 + xp^2 + xp*(2 - 4*zp)] + 
           (1 - 4/Sqrt[1 + 2*xp + xp^2 - 4*xp*zp])/zp))*Log[xp] - 
       (4*((1 - 2*zp)^2 + xp^6*(1 - 2*zp)^2 + xp*(2 - 20*zp + 48*zp^2 - 
            32*zp^3) + xp^5*(2 - 20*zp + 48*zp^2 - 32*zp^3) + 
          xp^2*(-1 - 16*zp + 104*zp^2 - 176*zp^3 + 88*zp^4) + 
          xp^4*(-1 - 16*zp + 104*zp^2 - 176*zp^3 + 88*zp^4) - 
          4*xp^3*(1 - 12*zp^2 + 28*zp^3 - 30*zp^4 + 12*zp^5))*
         Log[1 - xp + Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]])/
        ((1 + xp^2 + xp*(2 - 4*zp))^(5/2)*(-1 + zp)*zp) + 
       (4*((1 - 2*zp)^2 + xp^6*(1 - 2*zp)^2 + xp*(2 - 20*zp + 48*zp^2 - 
            32*zp^3) + xp^5*(2 - 20*zp + 48*zp^2 - 32*zp^3) + 
          xp^2*(-1 - 16*zp + 104*zp^2 - 176*zp^3 + 88*zp^4) + 
          xp^4*(-1 - 16*zp + 104*zp^2 - 176*zp^3 + 88*zp^4) - 
          4*xp^3*(1 - 12*zp^2 + 28*zp^3 - 30*zp^4 + 12*zp^5))*
         Log[-1 + xp + Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]])/
        ((1 + xp^2 + xp*(2 - 4*zp))^(5/2)*(-1 + zp)*zp) + 
       ((4 + 8*zp + 8*xp*zp^2 + xp^2*(4 + 8*zp))*
         Log[1 - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/(zp - xp*zp)) - 
     (4*(2*(-1 + zp)^2 + 8*xp^3*zp + 2*xp^2*(1 - 4*zp + 7*zp^2) + 
        xp*(-1 + 13*zp - 7*zp^2 + 3*zp^3))*
       Log[1 - 2*xp - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/
      ((-1 + xp)*zp*Sqrt[1 + (-2 + 4*xp)*zp + zp^2]) + 
     ((2 + 4*zp + 4*xp*zp^2 + xp^2*(2 + 4*zp))*
       Log[1 - 2*xp - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]]^2)/(zp - xp*zp) + 
     Log[1 - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]]*
      ((-4*(2*(-1 + zp)^2 + 8*xp^3*zp + 2*xp^2*(1 - 4*zp + 7*zp^2) + 
          xp*(-1 + 13*zp - 7*zp^2 + 3*zp^3)))/((-1 + xp)*zp*
         Sqrt[1 + (-2 + 4*xp)*zp + zp^2]) + 
       ((4 + 8*zp + 8*xp*zp^2 + xp^2*(4 + 8*zp))*
         Log[1 - 2*xp - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/
        ((-1 + xp)*zp)) + Log[1 - xp]*
      ((4*(-15 + 6/(-1 + xp) + xp*(4 + 9/(-1 + zp)) - 9/(-1 + zp) - 11*zp))/
        3 + (8*(-1 + 2*zp + xp^2*(-1 + 2*zp) + xp*(2 - 4*zp + 3*zp^2))*
         Log[xp])/((-1 + xp)*(-1 + zp)) - 
       (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*Log[1 - zp])/
        ((-1 + xp)*(-1 + zp)) - (4*(-3 + 6*zp + xp^2*(-3 + 6*zp) + 
          2*xp*(3 - 6*zp + 2*zp^2))*Log[zp])/((-1 + xp)*(-1 + zp)) + 
       ((4 + 8*zp + 8*xp*zp^2 + xp^2*(4 + 8*zp))*
         Log[1 + zp - Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/(zp - xp*zp) + 
       ((4 + 8*zp + 8*xp*zp^2 + xp^2*(4 + 8*zp))*
         Log[1 - 2*xp - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/
        ((-1 + xp)*zp)) + Log[xp]*(6 - 6/(-1 + xp) - 
       6/(1 + xp^2 + xp*(2 - 4*zp))^2 + (6*xp^4)/(1 + xp^2 + xp*(2 - 4*zp))^
         2 - (6*xp^5)/(1 + xp^2 + xp*(2 - 4*zp))^2 + 24/(-1 + zp) + 
       4/(1 + 2*xp + xp^2 - 4*xp*zp) + 18/(zp - xp*zp) - 
       12/((zp - xp*zp)*Sqrt[1 - 2*zp + 4*xp*zp + zp^2]) + 
       52/Sqrt[1 + (-2 + 4*xp)*zp + zp^2] + 
       36/((-1 + xp)*Sqrt[1 + (-2 + 4*xp)*zp + zp^2]) + 
       (12*xp*zp^2)/((-1 + xp)*Sqrt[1 + (-2 + 4*xp)*zp + zp^2]) + 
       xp^3*(8/(xp - zp)^2 + 4/(1 + 2*xp + xp^2 - 4*xp*zp)) + 
       (4 + 4/Sqrt[1 - 2*zp + 4*xp*zp + zp^2])/zp + 
       xp^2*(-4/(xp - zp)^2 + 8/(-xp + zp) + 
         32/Sqrt[1 + (-2 + 4*xp)*zp + zp^2]) + 
       xp*(-74/3 + 6/(1 + xp^2 + xp*(2 - 4*zp))^2 - 24/(-1 + zp) + 
         4/(-1 + xp + zp) + 8/(zp*Sqrt[1 + (-2 + 4*xp)*zp + zp^2]) + 
         (56*zp)/Sqrt[1 + (-2 + 4*xp)*zp + zp^2]) + 
       (4*zp*(1 - 6/Sqrt[1 + (-2 + 4*xp)*zp + zp^2] + 
          xp*(-31 - 21/Sqrt[1 + (-2 + 4*xp)*zp + zp^2])))/(3 - 3*xp) - 
       (4*(3*zp + 3*xp^2*zp + xp*(-1 + 3*zp^2))*ArcTan[Sqrt[zp]/Sqrt[xp]])/
        (Sqrt[xp]*zp^(3/2)) - (4*(3*zp + 3*xp^2*zp + xp*(-1 + 3*zp^2))*
         ArcTan[Sqrt[xp]*Sqrt[zp]])/(Sqrt[xp]*zp^(3/2)) + 
       (8*(1 + xp)*(-1 + 2*zp)*Log[1 + xp])/zp + 
       (8*(-2 + 4*zp + xp^2*(-2 + 4*zp) + xp*(4 - 8*zp + 5*zp^2))*
         Log[1 - zp])/((-1 + xp)*(-1 + zp)) + 
       (4*(zp + 6*zp^2 - zp^3 + 2*xp^3*(1 - 4*zp + 6*zp^2) + 
          2*xp^2*(-1 + 8*zp - 11*zp^2 + 9*zp^3) + 
          xp*(4 + zp - 4*zp^2 + 9*zp^3))*Log[zp])/((-1 + xp)*(1 + xp)*
         (-1 + zp)*zp) - (8*(1 + 2*xp*zp + 2*xp^3*zp + xp^2*(1 + 2*zp^2))*
         Log[1 + zp/xp])/((-1 + xp^2)*zp) + 
       (8*(xp + xp^3 + 2*zp + 2*xp^2*zp + 2*xp*zp^2)*Log[1 + xp*zp])/
        ((-1 + xp)*(1 + xp)*zp) + (4*((1 - 2*zp)^2 + xp^6*(1 - 2*zp)^2 + 
          xp*(2 - 20*zp + 48*zp^2 - 32*zp^3) + xp^5*(2 - 20*zp + 48*zp^2 - 
            32*zp^3) + xp^2*(-1 - 16*zp + 104*zp^2 - 176*zp^3 + 88*zp^4) + 
          xp^4*(-1 - 16*zp + 104*zp^2 - 176*zp^3 + 88*zp^4) - 
          4*xp^3*(1 - 12*zp^2 + 28*zp^3 - 30*zp^4 + 12*zp^5))*
         Log[-1 + xp + Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]])/
        ((1 + xp^2 + xp*(2 - 4*zp))^(5/2)*(-1 + zp)*zp) + 
       (4*((1 - 2*zp)^2 + xp^6*(1 - 2*zp)^2 + xp*(2 - 20*zp + 48*zp^2 - 
            32*zp^3) + xp^5*(2 - 20*zp + 48*zp^2 - 32*zp^3) + 
          xp^2*(-1 - 16*zp + 104*zp^2 - 176*zp^3 + 88*zp^4) + 
          xp^4*(-1 - 16*zp + 104*zp^2 - 176*zp^3 + 88*zp^4) - 
          4*xp^3*(1 - 12*zp^2 + 28*zp^3 - 30*zp^4 + 12*zp^5))*
         Log[1 + xp + Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]])/
        ((1 + xp^2 + xp*(2 - 4*zp))^(5/2)*(-1 + zp)*zp) + 
       ((4 + 8*zp + 8*xp*zp^2 + xp^2*(4 + 8*zp))*
         Log[1 + zp - Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/((-1 + xp)*zp) + 
       ((4 + 8*zp + 8*xp*zp^2 + xp^2*(4 + 8*zp))*
         Log[1 - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/(zp - xp*zp) + 
       ((4 + 8*zp + 8*xp*zp^2 + xp^2*(4 + 8*zp))*
         Log[1 - 2*xp - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/
        (zp - xp*zp)) + Delz*((4*(-271 + 446*xp))/27 - 
       (4*(1 + 10*xp)*zeta2)/3 - 14*(1 + xp)*zeta3 + 
       (-808/27 + (44*zeta2)/3 + 28*zeta3)*Dxp[0] + (268/9 - 8*zeta2)*
        Dxp[1] - (22*Dxp[2])/3 + ((34 + 10*xp^2)*Li[2, 1 - xp])/(3 - 3*xp) + 
       (16*(1 + xp^2)*Li[3, 1 - xp])/(-1 + xp) + (11*(1 + xp)*Log[1 - xp]^2)/
        3 + ((2*(78 - 91*xp + 83*xp^2))/(3*(-1 + xp)) - (4*(1 + xp^2)*zeta2)/
          (-1 + xp) - (4*(1 + xp^2)*Li[2, 1 - xp])/(-1 + xp))*Log[xp] + 
       ((79 - 12*xp + 43*xp^2)*Log[xp]^2)/(6*(-1 + xp)) + 
       ((1 + xp^2)*Log[xp]^3)/(-1 + xp) + Log[1 - xp]*((10*(4 - 29*xp))/9 + 
         4*(1 + xp)*zeta2 - (4*(1 + xp^2)*Li[2, 1 - xp])/(-1 + xp) + 
         ((56 + 32*xp^2)*Log[xp])/(3 - 3*xp) - (2*(1 + xp^2)*Log[xp]^2)/
          (-1 + xp)) + (4*(1 + xp^2)*S12[1 - xp])/(-1 + xp)) + 
     Delx*(Delz*(-1535/12 - (538*zeta2)/9 + (68*zeta2^2)/5 + (172*zeta3)/3) - 
       14*zeta3*(1 + zp) - (4*zeta2*(7 + 4*zp))/3 + (2*(85 + 319*zp))/27 + 
       (-808/27 + (44*zeta2)/3 + 28*zeta3)*Dzp[0] + (268/9 - 8*zeta2)*
        Dzp[1] - (22*Dzp[2])/3 + 4*(1 + zp)*Li[2, 1 - zp] - 
       (12*(1 + zp^2)*Li[3, 1 - zp])/(-1 + zp) + (11*(1 + zp)*Log[1 - zp]^2)/
        3 + ((2*(49 - 58/(-1 + zp) - 80*zp))/9 + (4*(1 + zp^2)*Li[2, 1 - zp])/
          (-1 + zp))*Log[zp] + ((-35 + 12*zp + zp^2)*Log[zp]^2)/
        (6*(-1 + zp)) + ((5 + 5*zp^2)*Log[zp]^3)/(3 - 3*zp) + 
       Log[1 - zp]*((22*(1 - 14*zp))/9 + 4*zeta2*(1 + zp) + 
         (4*(1 + zp^2)*Li[2, 1 - zp])/(-1 + zp) + 4*(1 + zp)*Log[zp]) + 
       (8*(1 + zp^2)*S12[1 - zp])/(-1 + zp)) + Theta[1 - xp - zp]*
      ((-8*zeta2*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp)))/
        ((-1 + xp)*(-1 + zp)) + (4*(-4 + 11*zp + 4*xp^2*zp - 3*zp^2 + 
          xp*(5 - 19*zp + 6*zp^2)))/((-1 + xp)*(-1 + zp)) + 
       (4*(xp^2*(1 - 2*zp) - 2*xp*(-1 + zp)^2 + (-2 + zp)*zp)*
         Li[2, -(zp/(1 - xp - zp))])/((-1 + xp)*(-1 + zp)) + 
       ((8*xp*(-1 + zp)^2 - 4*(-2 + zp)*zp + xp^2*(-4 + 8*zp))*
         Li[2, -((xp*zp)/(1 - xp - zp))])/((-1 + xp)*(-1 + zp)) + 
       ((4 + xp^2*(4 - 8*zp) - 8*xp*(-1 + zp)^2 - 8*zp)*
         Li[2, (1 - xp)^(-1) - xp/(1 - xp) - zp/(1 - xp)])/
        ((-1 + xp)*(-1 + zp)) + ((-8 + 4*xp*(-1 + zp)^2 + 4*zp + 6*zp^2 + 
          xp^2*(-2 + 4*zp))*Log[xp]^2)/((-1 + xp)*(-1 + zp)) + 
       (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*Log[1 - zp]^2)/
        ((-1 + xp)*(-1 + zp)) + (4*(-3 + xp*(2 - 8*zp) + 6*zp + 
          xp^2*(-1 + 4*zp))*Log[1 - xp - zp])/((-1 + xp)*(-1 + zp)) + 
       ((-8*(-3 + xp*(2 - 8*zp) + 6*zp + xp^2*(-1 + 4*zp)))/
          ((-1 + xp)*(-1 + zp)) - (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + 
            xp^2*(-1 + 2*zp))*Log[1 - xp - zp])/((-1 + xp)*(-1 + zp)))*
        Log[zp] + (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*
         Log[zp]^2)/((-1 + xp)*(-1 + zp)) + Log[1 - xp]*
        ((-4*(-3 + xp*(2 - 8*zp) + 6*zp + xp^2*(-1 + 4*zp)))/
          ((-1 + xp)*(-1 + zp)) - (8*(1 + zp)*Log[xp])/(-1 + xp) + 
         (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*Log[1 - zp])/
          ((-1 + xp)*(-1 + zp)) + (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + 
            xp^2*(-1 + 2*zp))*Log[zp])/((-1 + xp)*(-1 + zp))) + 
       Log[1 - zp]*((-8*(-3 + xp*(2 - 8*zp) + 6*zp + xp^2*(-1 + 4*zp)))/
          ((-1 + xp)*(-1 + zp)) - (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + 
            xp^2*(-1 + 2*zp))*Log[1 - xp - zp])/((-1 + xp)*(-1 + zp)) + 
         (16*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*Log[zp])/
          ((-1 + xp)*(-1 + zp))) + Log[xp]*
        ((4*(-9 + xp*(11 - 17*zp) + 11*zp + zp^2 + xp^2*(-5 + 8*zp)))/
          ((-1 + xp)*(-1 + zp)) - (8*(-2 + 2*xp*(-1 + zp)^2 + 2*zp + zp^2 + 
            xp^2*(-1 + 2*zp))*Log[1 - zp])/((-1 + xp)*(-1 + zp)) + 
         (4*(-2 + 2*xp*(-1 + zp)^2 + 2*zp + zp^2 + xp^2*(-1 + 2*zp))*
           Log[1 - xp - zp])/((-1 + xp)*(-1 + zp)) - 
         (4*(-4 + 6*xp*(-1 + zp)^2 + 6*zp + zp^2 + xp^2*(-3 + 6*zp))*Log[zp])/
          ((-1 + xp)*(-1 + zp))) + 
       ((4*(3*zp^2 + xp^2*(-5 + 8*zp) + xp*(4 - 4*zp - 6*zp^2)))/
          ((-1 + xp)*(-1 + zp)) + (8*zeta2*(-((-2 + zp)*zp) + 
            xp^2*(-1 + 2*zp) + xp*(2 - 4*zp + 4*zp^2)))/
          ((-1 + xp)*(-1 + zp)) + ((4 + (-4 + 8*xp)*zp^2)*
           Li[2, (1 - zp)^(-1) - zp/(xp*(1 - zp))])/((-1 + xp)*(-1 + zp)) + 
         ((-4 + 8*xp*(-1 + zp)^2 + 8*zp + xp^2*(-4 + 8*zp))*
           Li[2, (1 - xp)^(-1) - xp/(1 - xp) - zp/((1 - xp)*xp) + 
             zp^2/((1 - xp)*xp)])/((-1 + xp)*(-1 + zp)) + 
         ((-8 + (8 - 16*xp)*zp^2)*Log[1 - xp]^2)/((-1 + xp)*(-1 + zp)) - 
         (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*Log[1 - zp]^2)/
          ((-1 + xp)*(-1 + zp)) - (4*(-3 + xp*(2 - 8*zp) + 6*zp + 
            xp^2*(-1 + 4*zp))*Log[1 - xp - zp])/((-1 + xp)*(-1 + zp)) + 
         ((4*(-4 + zp)*zp + 8*xp*(2 + zp) - 4*xp^2*(1 + 2*zp))*Log[xp - zp])/
          ((-1 + xp)*(-1 + zp)) + 
         ((8*(-((-4 + zp)*zp) - 2*xp*(2 + zp) + xp^2*(1 + 2*zp)))/
            ((-1 + xp)*(-1 + zp)) + (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + 
              xp^2*(-1 + 2*zp))*Log[1 - xp - zp])/((-1 + xp)*(-1 + zp)) + 
           (8*(-((-2 + zp)*zp) + xp^2*(-1 + 2*zp) + xp*(2 - 4*zp + 4*zp^2))*
             Log[xp - zp])/((-1 + xp)*(-1 + zp)))*Log[zp] - 
         (8*(-((-2 + zp)*zp) + xp^2*(-1 + 2*zp) + xp*(2 - 4*zp + 4*zp^2))*
           Log[zp]^2)/((-1 + xp)*(-1 + zp)) + 
         Log[xp]*((4*(-4 + zp)*zp + 8*xp*(2 + zp) - 4*xp^2*(1 + 2*zp))/
            ((-1 + xp)*(-1 + zp)) + (8*(-((-2 + zp)*zp) + xp^2*(-1 + 2*zp) + 
              xp*(2 - 4*zp + 4*zp^2))*Log[1 - zp])/((-1 + xp)*(-1 + zp)) - 
           (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*
             Log[1 - xp - zp])/((-1 + xp)*(-1 + zp)) - 
           (8*(-((-2 + zp)*zp) + xp^2*(-1 + 2*zp) + xp*(2 - 4*zp + 4*zp^2))*
             Log[xp - zp])/((-1 + xp)*(-1 + zp)) + 
           (8*(-((-2 + zp)*zp) + xp^2*(-1 + 2*zp) + xp*(2 - 4*zp + 4*zp^2))*
             Log[zp])/((-1 + xp)*(-1 + zp))) + Log[1 - zp]*
          ((-4*(3 - 10*zp - 6*xp^2*zp + zp^2 + 2*xp*(1 + 5*zp)))/
            ((-1 + xp)*(-1 + zp)) + (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + 
              xp^2*(-1 + 2*zp))*Log[1 - xp - zp])/((-1 + xp)*(-1 + zp)) + 
           (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*Log[xp - zp])/
            ((-1 + xp)*(-1 + zp)) - (8*(-1 + 4*zp - zp^2 + xp^2*(-2 + 4*zp) + 
              xp*(4 - 8*zp + 6*zp^2))*Log[zp])/((-1 + xp)*(-1 + zp))) + 
         Log[1 - xp]*((4*(3 + 3*xp^2 + 2*zp - 2*zp^2 + 2*xp*(-5 + 2*zp)))/
            ((-1 + xp)*(-1 + zp)) + (8*(-((-2 + zp)*zp) + xp^2*(-1 + 2*zp) + 
              xp*(2 - 4*zp + 4*zp^2))*Log[xp])/((-1 + xp)*(-1 + zp)) - 
           (8*(-((-2 + zp)*zp) + xp^2*(-1 + 2*zp) + xp*(2 - 4*zp + 4*zp^2))*
             Log[1 - zp])/((-1 + xp)*(-1 + zp)) + ((8 + 8*(-1 + 2*xp)*zp^2)*
             Log[xp - zp])/((-1 + xp)*(-1 + zp)) - 
           (8*(1 + 2*zp - 2*zp^2 + xp^2*(-1 + 2*zp) + xp*(2 - 4*zp + 6*zp^2))*
             Log[zp])/((-1 + xp)*(-1 + zp))))*Theta[xp - zp] + 
       ((zeta2*(-4 + 8*xp*(-1 + zp)^2 + 8*zp + xp^2*(-4 + 8*zp)))/
          ((-1 + xp)*(-1 + zp)) - (4*(-4 + 11*zp + 4*xp^2*zp - 3*zp^2 + 
            xp*(5 - 19*zp + 6*zp^2)))/((-1 + xp)*(-1 + zp)) + 
         ((4 + xp^2*(4 - 8*zp) - 8*xp*(-1 + zp)^2 - 8*zp)*
           Li[2, -(xp/((1 - xp - zp)*(-xp + zp))) + xp^2/((1 - xp - zp)*(
                -xp + zp))])/((-1 + xp)*(-1 + zp)) + 
         ((2 + xp^2*(2 - 4*zp) - 4*xp*(-1 + zp)^2 - 4*zp)*Log[1 - xp]^2)/
          ((-1 + xp)*(-1 + zp)) + ((2 + xp^2*(2 - 4*zp) - 4*xp*(-1 + zp)^2 - 
            4*zp)*Log[xp]^2)/((-1 + xp)*(-1 + zp)) - 
         (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*Log[1 - zp]^2)/
          ((-1 + xp)*(-1 + zp)) + ((2 + xp^2*(2 - 4*zp) - 4*xp*(-1 + zp)^2 - 
            4*zp)*Log[1 - xp - zp]^2)/((-1 + xp)*(-1 + zp)) - 
         (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*Log[zp]^2)/
          ((-1 + xp)*(-1 + zp)) - (4*(-3 + xp*(2 - 8*zp) + 6*zp + 
            xp^2*(-1 + 4*zp))*Log[-xp + zp])/((-1 + xp)*(-1 + zp)) + 
         ((2 + xp^2*(2 - 4*zp) - 4*xp*(-1 + zp)^2 - 4*zp)*Log[-xp + zp]^2)/
          ((-1 + xp)*(-1 + zp)) + Log[1 - xp - zp]*
          ((-4*(-3 + xp*(2 - 8*zp) + 6*zp + xp^2*(-1 + 4*zp)))/
            ((-1 + xp)*(-1 + zp)) + ((4 + xp^2*(4 - 8*zp) - 8*xp*(-1 + zp)^
                2 - 8*zp)*Log[-xp + zp])/((-1 + xp)*(-1 + zp))) + 
         Log[xp]*((-4*(-3 + xp*(2 - 8*zp) + 6*zp + xp^2*(-1 + 4*zp)))/
            ((-1 + xp)*(-1 + zp)) + (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + 
              xp^2*(-1 + 2*zp))*Log[1 - zp])/((-1 + xp)*(-1 + zp)) + 
           ((4 + xp^2*(4 - 8*zp) - 8*xp*(-1 + zp)^2 - 8*zp)*Log[1 - xp - zp])/
            ((-1 + xp)*(-1 + zp)) + (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + 
              xp^2*(-1 + 2*zp))*Log[zp])/((-1 + xp)*(-1 + zp)) + 
           ((4 + xp^2*(4 - 8*zp) - 8*xp*(-1 + zp)^2 - 8*zp)*Log[-xp + zp])/
            ((-1 + xp)*(-1 + zp))) + Log[zp]*
          ((8*(-3 + xp*(2 - 8*zp) + 6*zp + xp^2*(-1 + 4*zp)))/
            ((-1 + xp)*(-1 + zp)) + (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + 
              xp^2*(-1 + 2*zp))*Log[1 - xp - zp])/((-1 + xp)*(-1 + zp)) + 
           (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*
             Log[-xp + zp])/((-1 + xp)*(-1 + zp))) + 
         Log[1 - zp]*((8*(-3 + xp*(2 - 8*zp) + 6*zp + xp^2*(-1 + 4*zp)))/
            ((-1 + xp)*(-1 + zp)) + (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + 
              xp^2*(-1 + 2*zp))*Log[1 - xp - zp])/((-1 + xp)*(-1 + zp)) - 
           (16*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*Log[zp])/
            ((-1 + xp)*(-1 + zp)) + (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + 
              xp^2*(-1 + 2*zp))*Log[-xp + zp])/((-1 + xp)*(-1 + zp))) + 
         Log[1 - xp]*((4*(-3 + xp*(2 - 8*zp) + 6*zp + xp^2*(-1 + 4*zp)))/
            ((-1 + xp)*(-1 + zp)) + ((-4 + 8*xp*(-1 + zp)^2 + 8*zp + 
              xp^2*(-4 + 8*zp))*Log[xp])/((-1 + xp)*(-1 + zp)) - 
           (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*Log[1 - zp])/
            ((-1 + xp)*(-1 + zp)) + ((-4 + 8*xp*(-1 + zp)^2 + 8*zp + 
              xp^2*(-4 + 8*zp))*Log[1 - xp - zp])/((-1 + xp)*(-1 + zp)) - 
           (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*Log[zp])/
            ((-1 + xp)*(-1 + zp)) + ((-4 + 8*xp*(-1 + zp)^2 + 8*zp + 
              xp^2*(-4 + 8*zp))*Log[-xp + zp])/((-1 + xp)*(-1 + zp))))*
        Theta[-xp + zp]) + 
     ((zeta2*(4 + xp^2*(4 - 8*zp) - 8*xp*(-1 + zp)^2 - 8*zp))/
        ((-1 + xp)*(-1 + zp)) + (4*(-4 + 11*zp + 4*xp^2*zp - 3*zp^2 + 
          xp*(5 - 19*zp + 6*zp^2)))/((-1 + xp)*(-1 + zp)) + 
       ((-4 + 8*xp*(-1 + zp)^2 + 8*zp + xp^2*(-4 + 8*zp))*
         Li[2, (1 - xp - zp)^(-1) - xp/(1 - xp - zp)])/
        ((-1 + xp)*(-1 + zp)) + (4*(xp^2*(1 - 2*zp) - 2*xp*(-1 + zp)^2 + 
          (-2 + zp)*zp)*Li[2, xp^(-1) + zp^(-1) - 1/(xp*zp)])/
        ((-1 + xp)*(-1 + zp)) + ((8*xp*(-1 + zp)^2 - 4*(-2 + zp)*zp + 
          xp^2*(-4 + 8*zp))*Li[2, 1 - zp^(-1) + xp/zp])/
        ((-1 + xp)*(-1 + zp)) + ((-2 + 4*xp*(-1 + zp)^2 + 4*zp + 
          xp^2*(-2 + 4*zp))*Log[1 - xp]^2)/((-1 + xp)*(-1 + zp)) + 
       (8*(1 + zp)*Log[xp]^2)/(-1 + xp) + 
       (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*Log[1 - zp]^2)/
        ((-1 + xp)*(-1 + zp)) + (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + 
          xp^2*(-1 + 2*zp))*Log[zp]^2)/((-1 + xp)*(-1 + zp)) + 
       (4*(-3 + xp*(2 - 8*zp) + 6*zp + xp^2*(-1 + 4*zp))*Log[-1 + xp + zp])/
        ((-1 + xp)*(-1 + zp)) + ((-2 + 4*xp*(-1 + zp)^2 + 4*zp + 
          xp^2*(-2 + 4*zp))*Log[-1 + xp + zp]^2)/((-1 + xp)*(-1 + zp)) + 
       Log[1 - xp]*((-4*(-3 + xp*(2 - 8*zp) + 6*zp + xp^2*(-1 + 4*zp)))/
          ((-1 + xp)*(-1 + zp)) - (8*(1 + zp)*Log[xp])/(-1 + xp) + 
         (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*Log[1 - zp])/
          ((-1 + xp)*(-1 + zp)) + (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + 
            xp^2*(-1 + 2*zp))*Log[zp])/((-1 + xp)*(-1 + zp)) + 
         ((4 + xp^2*(4 - 8*zp) - 8*xp*(-1 + zp)^2 - 8*zp)*Log[-1 + xp + zp])/
          ((-1 + xp)*(-1 + zp))) + Log[zp]*
        ((-8*(-3 + xp*(2 - 8*zp) + 6*zp + xp^2*(-1 + 4*zp)))/
          ((-1 + xp)*(-1 + zp)) - (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + 
            xp^2*(-1 + 2*zp))*Log[-1 + xp + zp])/((-1 + xp)*(-1 + zp))) + 
       Log[1 - zp]*((-8*(-3 + xp*(2 - 8*zp) + 6*zp + xp^2*(-1 + 4*zp)))/
          ((-1 + xp)*(-1 + zp)) + (16*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + 
            xp^2*(-1 + 2*zp))*Log[zp])/((-1 + xp)*(-1 + zp)) - 
         (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*
           Log[-1 + xp + zp])/((-1 + xp)*(-1 + zp))) + 
       Log[xp]*((4*(-9 + xp*(11 - 17*zp) + 11*zp + zp^2 + xp^2*(-5 + 8*zp)))/
          ((-1 + xp)*(-1 + zp)) - (8*(-2 + 2*xp*(-1 + zp)^2 + 2*zp + zp^2 + 
            xp^2*(-1 + 2*zp))*Log[1 - zp])/((-1 + xp)*(-1 + zp)) - 
         (16*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*Log[zp])/
          ((-1 + xp)*(-1 + zp)) + (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + 
            xp^2*(-1 + 2*zp))*Log[-1 + xp + zp])/((-1 + xp)*(-1 + zp))))*
      Theta[-1 + xp + zp] + Theta[-xp + zp]*
      ((4*(-4 + 11*zp + 4*xp^2*zp - 3*zp^2 + xp*(5 - 19*zp + 6*zp^2)))/
        ((-1 + xp)*(-1 + zp)) + (8*zeta2*(-((-2 + zp)*zp) + 
          xp^2*(-1 + 2*zp) + xp*(2 - 4*zp + 6*zp^2)))/((-1 + xp)*(-1 + zp)) - 
       (4*(-1 + 4*zp - zp^2 + xp^2*(-2 + 4*zp) + xp*(4 - 8*zp + 8*zp^2))*
         Li[2, xp/zp])/((-1 + xp)*(-1 + zp)) + 
       ((4 + xp^2*(4 - 8*zp) - 8*xp*(-1 + zp)^2 - 8*zp)*
         Li[2, -(xp/(1 - xp)) + zp/(1 - xp)])/((-1 + xp)*(-1 + zp)) + 
       ((4 + (-4 + 8*xp)*zp^2)*Li[2, -(-xp + zp)^(-1) + zp/(-xp + zp)])/
        ((-1 + xp)*(-1 + zp)) + ((-4 + (4 - 8*xp)*zp^2)*
         Li[2, -(xp/(-xp + zp)) + (xp*zp)/(-xp + zp)])/
        ((-1 + xp)*(-1 + zp)) + (4*(-1 + 6*zp - 2*zp^2 + xp^2*(-3 + 6*zp) + 
          xp*(6 - 12*zp + 11*zp^2))*Log[xp]^2)/((-1 + xp)*(-1 + zp)) + 
       (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*Log[1 - zp]^2)/
        ((-1 + xp)*(-1 + zp)) + (2*(-3 + xp*(4 - 8*zp) + 4*zp + zp^2 + 
          xp^2*(-2 + 4*zp))*Log[zp]^2)/((-1 + xp)*(-1 + zp)) + 
       Log[1 - xp]*((-4*(-3 + xp*(2 - 8*zp) + 6*zp + xp^2*(-1 + 4*zp)))/
          ((-1 + xp)*(-1 + zp)) - (8*(-((-2 + zp)*zp) + xp^2*(-1 + 2*zp) + 
            xp*(2 - 4*zp + 4*zp^2))*Log[xp])/((-1 + xp)*(-1 + zp)) + 
         (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*Log[1 - zp])/
          ((-1 + xp)*(-1 + zp)) + (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + 
            xp^2*(-1 + 2*zp))*Log[zp])/((-1 + xp)*(-1 + zp))) + 
       (4*(-3 + xp*(2 - 8*zp) + 6*zp + xp^2*(-1 + 4*zp))*Log[-xp + zp])/
        ((-1 + xp)*(-1 + zp)) + Log[xp]*
        ((-4*(3 - 10*zp - 6*xp^2*zp + zp^2 + 2*xp*(1 + 5*zp)))/
          ((-1 + xp)*(-1 + zp)) - (4*(-3 + 8*zp - zp^2 + xp^2*(-4 + 8*zp) + 
            2*xp*(4 - 8*zp + 5*zp^2))*Log[1 - zp])/((-1 + xp)*(-1 + zp)) - 
         (8*(-1 + 4*zp - zp^2 + xp^2*(-2 + 4*zp) + xp*(4 - 8*zp + 6*zp^2))*
           Log[zp])/((-1 + xp)*(-1 + zp)) - (8*xp*zp^2*Log[-xp + zp])/
          ((-1 + xp)*(-1 + zp))) + Log[zp]*
        ((-8*(-3 + xp*(2 - 8*zp) + 6*zp + xp^2*(-1 + 4*zp)))/
          ((-1 + xp)*(-1 + zp)) + ((4 + 4*(-1 + 4*xp)*zp^2)*Log[-xp + zp])/
          ((-1 + xp)*(-1 + zp))) + Log[1 - zp]*
        ((-8*(-3 + xp*(2 - 8*zp) + 6*zp + xp^2*(-1 + 4*zp)))/
          ((-1 + xp)*(-1 + zp)) + (16*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + 
            xp^2*(-1 + 2*zp))*Log[zp])/((-1 + xp)*(-1 + zp)) - 
         (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*Log[-xp + zp])/
          ((-1 + xp)*(-1 + zp))) + 
       ((8*zeta2*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp)))/
          ((-1 + xp)*(-1 + zp)) - (4*(-4 + 11*zp + 4*xp^2*zp - 3*zp^2 + 
            xp*(5 - 19*zp + 6*zp^2)))/((-1 + xp)*(-1 + zp)) + 
         ((-4 + 8*xp*(-1 + zp)^2 + 8*zp + xp^2*(-4 + 8*zp))*
           Li[2, (1 - xp)^(-1) - xp/(1 - xp) - zp/((1 - xp)*xp) + 
             zp^2/((1 - xp)*xp)])/((-1 + xp)*(-1 + zp)) - 
         (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*Log[1 - zp]^2)/
          ((-1 + xp)*(-1 + zp)) - (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + 
            xp^2*(-1 + 2*zp))*Log[zp]^2)/((-1 + xp)*(-1 + zp)) + 
         Log[1 - xp]*((4*(-3 + xp*(2 - 8*zp) + 6*zp + xp^2*(-1 + 4*zp)))/
            ((-1 + xp)*(-1 + zp)) + (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + 
              xp^2*(-1 + 2*zp))*Log[xp])/((-1 + xp)*(-1 + zp)) - 
           (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*Log[1 - zp])/
            ((-1 + xp)*(-1 + zp)) - (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + 
              xp^2*(-1 + 2*zp))*Log[zp])/((-1 + xp)*(-1 + zp))) - 
         (4*(-3 + xp*(2 - 8*zp) + 6*zp + xp^2*(-1 + 4*zp))*Log[-xp + zp])/
          ((-1 + xp)*(-1 + zp)) - (4*(-3 + xp*(2 - 8*zp) + 6*zp + 
            xp^2*(-1 + 4*zp))*Log[-1 + xp + zp])/((-1 + xp)*(-1 + zp)) + 
         Log[xp]*((-4*(-3 + xp*(2 - 8*zp) + 6*zp + xp^2*(-1 + 4*zp)))/
            ((-1 + xp)*(-1 + zp)) + (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + 
              xp^2*(-1 + 2*zp))*Log[1 - zp])/((-1 + xp)*(-1 + zp)) + 
           (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*Log[zp])/
            ((-1 + xp)*(-1 + zp)) - (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + 
              xp^2*(-1 + 2*zp))*Log[-xp + zp])/((-1 + xp)*(-1 + zp)) - 
           (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*
             Log[-1 + xp + zp])/((-1 + xp)*(-1 + zp))) + 
         Log[zp]*((8*(-3 + xp*(2 - 8*zp) + 6*zp + xp^2*(-1 + 4*zp)))/
            ((-1 + xp)*(-1 + zp)) + (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + 
              xp^2*(-1 + 2*zp))*Log[-xp + zp])/((-1 + xp)*(-1 + zp)) + 
           (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*
             Log[-1 + xp + zp])/((-1 + xp)*(-1 + zp))) + 
         Log[1 - zp]*((8*(-3 + xp*(2 - 8*zp) + 6*zp + xp^2*(-1 + 4*zp)))/
            ((-1 + xp)*(-1 + zp)) - (16*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + 
              xp^2*(-1 + 2*zp))*Log[zp])/((-1 + xp)*(-1 + zp)) + 
           (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*
             Log[-xp + zp])/((-1 + xp)*(-1 + zp)) + 
           (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*
             Log[-1 + xp + zp])/((-1 + xp)*(-1 + zp))))*
        Theta[-1 + xp + zp]) + Theta[xp - zp]*
      ((-4*(3*zp^2 + xp^2*(-5 + 8*zp) + xp*(4 - 4*zp - 6*zp^2)))/
        ((-1 + xp)*(-1 + zp)) - (4*zeta2*(1 + 2*zp - 2*zp^2 + 
          xp^2*(-1 + 2*zp) + xp*(2 - 4*zp + 6*zp^2)))/((-1 + xp)*(-1 + zp)) + 
       (4*(-1 + 4*zp - zp^2 + xp^2*(-2 + 4*zp) + xp*(4 - 8*zp + 8*zp^2))*
         Li[2, zp/xp])/((-1 + xp)*(-1 + zp)) + 
       ((-4 + (4 - 8*xp)*zp^2)*Li[2, xp/(1 - zp) - zp/(1 - zp)])/
        ((-1 + xp)*(-1 + zp)) + ((-4 + 8*xp*(-1 + zp)^2 + 8*zp + 
          xp^2*(-4 + 8*zp))*Li[2, (-xp + zp)^(-1) - xp/(-xp + zp)])/
        ((-1 + xp)*(-1 + zp)) + ((6 + 4*zp - 8*zp^2 + xp^2*(-2 + 4*zp) + 
          4*xp*(1 - 2*zp + 5*zp^2))*Log[1 - xp]^2)/((-1 + xp)*(-1 + zp)) + 
       (4*(-1 + 8*zp - 3*zp^2 + xp^2*(-4 + 8*zp) + 8*xp*(1 - 2*zp + 2*zp^2))*
         Log[xp]^2)/((-1 + xp)*(-1 + zp)) + 
       (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*Log[1 - zp]^2)/
        ((-1 + xp)*(-1 + zp)) + ((-4*(-4 + zp)*zp - 8*xp*(2 + zp) + 
          xp^2*(4 + 8*zp))*Log[xp - zp])/((-1 + xp)*(-1 + zp)) + 
       ((-2 + 4*xp*(-1 + zp)^2 + 4*zp + xp^2*(-2 + 4*zp))*Log[xp - zp]^2)/
        ((-1 + xp)*(-1 + zp)) + ((8*(-4 + zp)*zp + 16*xp*(2 + zp) - 
           8*xp^2*(1 + 2*zp))/((-1 + xp)*(-1 + zp)) + 
         (4*(1 + zp)*Log[xp - zp])/(-1 + xp))*Log[zp] + 
       (8*(-((-2 + zp)*zp) + xp^2*(-1 + 2*zp) + xp*(2 - 4*zp + 4*zp^2))*
         Log[zp]^2)/((-1 + xp)*(-1 + zp)) + 
       Log[xp]*((8*(-((-4 + zp)*zp) - 2*xp*(2 + zp) + xp^2*(1 + 2*zp)))/
          ((-1 + xp)*(-1 + zp)) - (8*(-1 + 4*zp - zp^2 + xp^2*(-2 + 4*zp) + 
            xp*(4 - 8*zp + 6*zp^2))*Log[1 - zp])/((-1 + xp)*(-1 + zp)) - 
         (4*(1 + zp)*Log[xp - zp])/(-1 + xp) + 
         ((4 + xp^2*(24 - 48*zp) - 48*zp + 20*zp^2 - 
            48*xp*(1 - 2*zp + 2*zp^2))*Log[zp])/((-1 + xp)*(-1 + zp))) + 
       Log[1 - zp]*((4*(3 - 10*zp - 6*xp^2*zp + zp^2 + 2*xp*(1 + 5*zp)))/
          ((-1 + xp)*(-1 + zp)) - (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + 
            xp^2*(-1 + 2*zp))*Log[xp - zp])/((-1 + xp)*(-1 + zp)) + 
         (8*(-1 + 4*zp - zp^2 + xp^2*(-2 + 4*zp) + xp*(4 - 8*zp + 6*zp^2))*
           Log[zp])/((-1 + xp)*(-1 + zp))) + Log[1 - xp]*
        ((-4*(3 + 3*xp^2 + 2*zp - 2*zp^2 + 2*xp*(-5 + 2*zp)))/
          ((-1 + xp)*(-1 + zp)) - (8*(1 + 2*zp - 2*zp^2 + xp^2*(-1 + 2*zp) + 
            xp*(2 - 4*zp + 6*zp^2))*Log[xp])/((-1 + xp)*(-1 + zp)) + 
         (8*(-((-2 + zp)*zp) + xp^2*(-1 + 2*zp) + xp*(2 - 4*zp + 4*zp^2))*
           Log[1 - zp])/((-1 + xp)*(-1 + zp)) - 
         (4*(1 + 2*zp - 2*zp^2 + xp^2*(-1 + 2*zp) + xp*(2 - 4*zp + 6*zp^2))*
           Log[xp - zp])/((-1 + xp)*(-1 + zp)) + 
         (8*(1 + 2*zp - 2*zp^2 + xp^2*(-1 + 2*zp) + xp*(2 - 4*zp + 6*zp^2))*
           Log[zp])/((-1 + xp)*(-1 + zp))) + 
       ((4*(3*zp^2 + xp^2*(-5 + 8*zp) + xp*(4 - 4*zp - 6*zp^2)))/
          ((-1 + xp)*(-1 + zp)) + (4*zeta2*(1 + 2*zp - 2*zp^2 + 
            xp^2*(-1 + 2*zp) + xp*(2 - 4*zp + 6*zp^2)))/
          ((-1 + xp)*(-1 + zp)) + ((4 + (-4 + 8*xp)*zp^2)*
           Li[2, (1 - zp)^(-1) - zp/(xp*(1 - zp))])/((-1 + xp)*(-1 + zp)) + 
         ((4 + xp^2*(4 - 8*zp) - 8*xp*(-1 + zp)^2 - 8*zp)*
           Li[2, -(xp/((1 - xp - zp)*(-xp + zp))) + xp^2/((1 - xp - zp)*(
                -xp + zp))])/((-1 + xp)*(-1 + zp)) + 
         ((-6 + xp^2*(2 - 4*zp) - 4*zp + 8*zp^2 - 4*xp*(1 - 2*zp + 5*zp^2))*
           Log[1 - xp]^2)/((-1 + xp)*(-1 + zp)) + 
         ((2 + xp^2*(2 - 4*zp) - 4*xp*(-1 + zp)^2 - 4*zp)*Log[xp]^2)/
          ((-1 + xp)*(-1 + zp)) - (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + 
            xp^2*(-1 + 2*zp))*Log[1 - zp]^2)/((-1 + xp)*(-1 + zp)) + 
         ((2 + xp^2*(2 - 4*zp) - 4*xp*(-1 + zp)^2 - 4*zp)*Log[xp - zp]^2)/
          ((-1 + xp)*(-1 + zp)) - (8*(-((-2 + zp)*zp) + xp^2*(-1 + 2*zp) + 
            xp*(2 - 4*zp + 4*zp^2))*Log[zp]^2)/((-1 + xp)*(-1 + zp)) - 
         (4*(-3 + xp*(2 - 8*zp) + 6*zp + xp^2*(-1 + 4*zp))*Log[-1 + xp + zp])/
          ((-1 + xp)*(-1 + zp)) + ((2 + xp^2*(2 - 4*zp) - 4*xp*(-1 + zp)^2 - 
            4*zp)*Log[-1 + xp + zp]^2)/((-1 + xp)*(-1 + zp)) + 
         Log[xp - zp]*((4*(-4 + zp)*zp + 8*xp*(2 + zp) - 4*xp^2*(1 + 2*zp))/
            ((-1 + xp)*(-1 + zp)) + ((4 + xp^2*(4 - 8*zp) - 8*xp*(-1 + zp)^
                2 - 8*zp)*Log[-1 + xp + zp])/((-1 + xp)*(-1 + zp))) + 
         Log[xp]*((4*(-4 + zp)*zp + 8*xp*(2 + zp) - 4*xp^2*(1 + 2*zp))/
            ((-1 + xp)*(-1 + zp)) + (8*(-((-2 + zp)*zp) + xp^2*(-1 + 2*zp) + 
              xp*(2 - 4*zp + 4*zp^2))*Log[1 - zp])/((-1 + xp)*(-1 + zp)) - 
           (4*(1 + 2*zp - 2*zp^2 + xp^2*(-1 + 2*zp) + xp*(2 - 4*zp + 6*zp^2))*
             Log[xp - zp])/((-1 + xp)*(-1 + zp)) + 
           (8*(-((-2 + zp)*zp) + xp^2*(-1 + 2*zp) + xp*(2 - 4*zp + 4*zp^2))*
             Log[zp])/((-1 + xp)*(-1 + zp)) + 
           ((4 + xp^2*(4 - 8*zp) - 8*xp*(-1 + zp)^2 - 8*zp)*
             Log[-1 + xp + zp])/((-1 + xp)*(-1 + zp))) + 
         Log[zp]*((8*(-((-4 + zp)*zp) - 2*xp*(2 + zp) + xp^2*(1 + 2*zp)))/
            ((-1 + xp)*(-1 + zp)) + (8*(-((-2 + zp)*zp) + xp^2*(-1 + 2*zp) + 
              xp*(2 - 4*zp + 4*zp^2))*Log[xp - zp])/((-1 + xp)*(-1 + zp)) + 
           (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*
             Log[-1 + xp + zp])/((-1 + xp)*(-1 + zp))) + 
         Log[1 - zp]*((-4*(3 - 10*zp - 6*xp^2*zp + zp^2 + 2*xp*(1 + 5*zp)))/
            ((-1 + xp)*(-1 + zp)) + (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + 
              xp^2*(-1 + 2*zp))*Log[xp - zp])/((-1 + xp)*(-1 + zp)) - 
           (8*(-1 + 4*zp - zp^2 + xp^2*(-2 + 4*zp) + xp*(4 - 8*zp + 6*zp^2))*
             Log[zp])/((-1 + xp)*(-1 + zp)) + 
           (8*(-1 + 2*xp*(-1 + zp)^2 + 2*zp + xp^2*(-1 + 2*zp))*
             Log[-1 + xp + zp])/((-1 + xp)*(-1 + zp))) + 
         Log[1 - xp]*((4*(3 + 3*xp^2 + 2*zp - 2*zp^2 + 2*xp*(-5 + 2*zp)))/
            ((-1 + xp)*(-1 + zp)) + (4*(1 + 2*zp - 2*zp^2 + xp^2*(-1 + 
                2*zp) + xp*(2 - 4*zp + 6*zp^2))*Log[xp])/
            ((-1 + xp)*(-1 + zp)) - (8*(-((-2 + zp)*zp) + xp^2*(-1 + 2*zp) + 
              xp*(2 - 4*zp + 4*zp^2))*Log[1 - zp])/((-1 + xp)*(-1 + zp)) + 
           (4*(1 + 2*zp - 2*zp^2 + xp^2*(-1 + 2*zp) + xp*(2 - 4*zp + 6*zp^2))*
             Log[xp - zp])/((-1 + xp)*(-1 + zp)) - 
           (8*(1 + 2*zp - 2*zp^2 + xp^2*(-1 + 2*zp) + xp*(2 - 4*zp + 6*zp^2))*
             Log[zp])/((-1 + xp)*(-1 + zp)) + 
           ((-4 + 8*xp*(-1 + zp)^2 + 8*zp + xp^2*(-4 + 8*zp))*
             Log[-1 + xp + zp])/((-1 + xp)*(-1 + zp))))*
        Theta[-1 + xp + zp])), "\!\(\*SubscriptBox[SuperscriptBox[\(G\), \
\((2)\)], \({1,qq,PS,MS}\)]\)" -> -(Delz* Zq2xpPS) + 
   CF*((-32*(-1 + xp)*(-1 + zp))/zp - (8*(xp + zp + xp^2*zp + xp*zp^2)*
       InvTanInt[-(Sqrt[zp]/Sqrt[xp])])/(Sqrt[xp]*zp^(3/2)) + 
     (8*(xp + zp + xp^2*zp + xp*zp^2)*InvTanInt[Sqrt[zp]/Sqrt[xp]])/
      (Sqrt[xp]*zp^(3/2)) + (8*(xp + zp + xp^2*zp + xp*zp^2)*
       InvTanInt[-(Sqrt[xp]*Sqrt[zp])])/(Sqrt[xp]*zp^(3/2)) - 
     (8*(xp + zp + xp^2*zp + xp*zp^2)*InvTanInt[Sqrt[xp]*Sqrt[zp]])/
      (Sqrt[xp]*zp^(3/2)) + ((12*(-1 + xp)*(1 + zp))/zp - 
       (8*(xp + zp + xp^2*zp + xp*zp^2)*ArcTan[Sqrt[zp]/Sqrt[xp]])/
        (Sqrt[xp]*zp^(3/2)) + (8*(xp + zp + xp^2*zp + xp*zp^2)*
         ArcTan[Sqrt[xp]*Sqrt[zp]])/(Sqrt[xp]*zp^(3/2)))*Log[zp] + 
     Log[xp]*((12*(1 + xp)*(-1 + zp))/zp + (8*(xp + zp + xp^2*zp + xp*zp^2)*
         ArcTan[Sqrt[zp]/Sqrt[xp]])/(Sqrt[xp]*zp^(3/2)) + 
       (8*(xp + zp + xp^2*zp + xp*zp^2)*ArcTan[Sqrt[xp]*Sqrt[zp]])/
        (Sqrt[xp]*zp^(3/2)) - (8*(1 + xp)*(1 + zp)*Log[zp])/zp)), 
 "\!\(\*SubscriptBox[SuperscriptBox[\(G\), \((2)\)], \
\({1,qq',\([1]\),MS}\)]\)" -> 
  CF*(((-1 + zp)*(5 - 19*zp - 130*zp^2 + xp*(-19 - 79*zp + 62*zp^2)))/
     (9*zp) - 4*(1 + xp + 2*zp)*Li[2, 1 - zp] - 
    (2*(2 + 15*zp - 9*zp^2 - 8*zp^3 + xp*(2 - 3*zp - 9*zp^2 + 4*zp^3))*
      Log[zp])/(3*zp) - 4*(1 + xp + 2*zp)*Log[zp]^2 + 
    Delx*Log[Q2/muF2]^2*(1 + 4/(3*zp) - zp - (4*zp^2)/3 + 
      2*(1 + zp)*Log[zp]) + Dxp[1]*(2 + 8/(3*zp) - 2*zp - (8*zp^2)/3 + 
      4*(1 + zp)*Log[zp]) + Log[1 - xp]*
     ((-2*(-1 + zp)*(-2 - 11*zp - 8*zp^2 + xp*(-2 - 5*zp + 4*zp^2)))/(3*zp) - 
      4*(1 + xp + 2*zp)*Log[zp]) + Log[1 - zp]*
     ((-2*(-1 + zp)*(-2 - 11*zp - 8*zp^2 + xp*(-2 - 5*zp + 4*zp^2)))/(3*zp) - 
      4*(1 + xp + 2*zp)*Log[zp]) + 
    Log[xp]*((2*(-1 + zp)*(-4 - zp + 8*zp^2 - 12*xp*zp*(1 + 2*zp) + 
         2*xp^2*(-2 - 5*zp + 4*zp^2)))/(3*(-1 + xp)*zp) + 
      (4*(1 + 2*xp^2 - zp + 4*xp*zp)*Log[zp])/(-1 + xp)) + 
    Dxp[0]*((2*(-7 - 60*zp + 42*zp^2 + 25*zp^3))/(9*zp) + 
      4*(1 + zp)*Li[2, 1 - zp] + (-2 + 8/(3*zp) - 8*zp - (8*zp^2)/3)*
       Log[zp] + 4*(1 + zp)*Log[zp]^2 + Log[1 - zp]*(2 + 8/(3*zp) - 2*zp - 
        (8*zp^2)/3 + 4*(1 + zp)*Log[zp])) + 
    Log[Q2/muF2]*((-2*(-1 + zp)*(-2 - 11*zp - 8*zp^2 + 
         xp*(-2 - 5*zp + 4*zp^2)))/(3*zp) - 4*(1 + xp + 2*zp)*Log[zp] + 
      Dxp[0]*(2 + 8/(3*zp) - 2*zp - (8*zp^2)/3 + 4*(1 + zp)*Log[zp]) + 
      Delx*((2*(-7 - 60*zp + 42*zp^2 + 25*zp^3))/(9*zp) + 
        4*(1 + zp)*Li[2, 1 - zp] + (-2 + 8/(3*zp) - 8*zp - (8*zp^2)/3)*
         Log[zp] + 4*(1 + zp)*Log[zp]^2 + Log[1 - zp]*(2 + 8/(3*zp) - 2*zp - 
          (8*zp^2)/3 + 4*(1 + zp)*Log[zp]))) + 
    Delx*(zeta2*(-2 - 8/(3*zp) + 2*zp + (8*zp^2)/3) - 
      (107 + 132*zp + 48*zp^2 - 287*zp^3)/(27*zp) + 
      (-2 + 8/(3*zp) - 8*zp - (8*zp^2)/3)*Li[2, 1 - zp] - 
      4*(1 + zp)*Li[3, 1 - zp] + (-4*zeta2*(1 + zp) - 
        (2*(7 + 90*zp + 81*zp^2 + 31*zp^3))/(9*zp) + 
        8*(1 + zp)*Li[2, 1 - zp])*Log[zp] + ((-3 + 8/zp - 15*zp)*Log[zp]^2)/
       6 + (5*(1 + zp)*Log[zp]^3)/3 + Log[1 - zp]^2*(1 + 4/(3*zp) - zp - 
        (4*zp^2)/3 + 2*(1 + zp)*Log[zp]) + Log[1 - zp]*
       ((2*(-7 - 60*zp + 42*zp^2 + 25*zp^3))/(9*zp) + 
        4*(1 + zp)*Li[2, 1 - zp] + (-2 + 8/(3*zp) - 8*zp - (8*zp^2)/3)*
         Log[zp] + 4*(1 + zp)*Log[zp]^2) + 8*(1 + zp)*S12[1 - zp])), 
 "\!\(\*SubscriptBox[SuperscriptBox[\(G\), \((2)\)], \
\({1,qq',\([2]\),MS}\)]\)" -> CF*((-4*(-1 + xp)*(-3 + 10*zp))/zp - 
    (4*(1 + xp)*(-1 + 2*zp)*Li[2, 1 - xp])/zp - 
    (8*(1 + xp^4 + xp*(2 - 4*zp) + xp^3*(2 - 4*zp) + 2*xp^2*(1 - zp + zp^2))*
      Li[2, 1/2 - xp/2 - Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]/2])/
     (1 + xp^2 + xp*(2 - 4*zp))^(3/2) - 
    (8*(1 + xp^4 + xp*(2 - 4*zp) + xp^3*(2 - 4*zp) + 2*xp^2*(1 - zp + zp^2))*
      Li[2, 1/2 + xp/2 - Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]/2])/
     (1 + xp^2 + xp*(2 - 4*zp))^(3/2) + 
    (8*(1 + xp^4 + xp*(2 - 4*zp) + xp^3*(2 - 4*zp) + 2*xp^2*(1 - zp + zp^2))*
      Li[2, 1/2 - 1/(2*xp) - Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]/(2*xp)])/
     (1 + xp^2 + xp*(2 - 4*zp))^(3/2) + 
    (8*(1 + xp^4 + xp*(2 - 4*zp) + xp^3*(2 - 4*zp) + 2*xp^2*(1 - zp + zp^2))*
      Li[2, 1/2 + 1/(2*xp) - Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]/(2*xp)])/
     (1 + xp^2 + xp*(2 - 4*zp))^(3/2) + 
    ((-4*(1 + xp^2 + xp*(2 - 4*zp))^(3/2) + zp + 5*xp^2*zp + 8*xp^4*zp + 
       7*(1 + xp^2 + xp*(2 - 4*zp))*zp + 8*(1 + xp^2 + xp*(2 - 4*zp))^(3/2)*
        zp + xp^3*(14 - 28*zp)*zp + 2*xp*(1 + xp^2 + xp*(2 - 4*zp))*
        (1 - 2*zp)*(-2*Sqrt[1 + xp^2 + xp*(2 - 4*zp)] + zp))*Log[xp]^2)/
     ((1 + xp^2 + xp*(2 - 4*zp))^(3/2)*zp) + Delz*Log[Q2/muF2]^2*
     (5 - 5*xp + 2*(1 + xp)*Log[xp]) + 
    Dzp[1]*(-10*(-1 + xp) + 4*(1 + xp)*Log[xp]) + 
    Log[1 - xp]*((10*(-1 + xp)*(-1 + 2*zp))/zp - 
      (4*(1 + xp)*(-1 + 2*zp)*Log[xp])/zp) + 
    Dzp[0]*(12*(-1 + xp) + 4*(1 + xp)*Li[2, 1 - xp] + 12*(-1 + xp)*Log[xp] - 
      4*(1 + xp)*Log[xp]^2 + Log[1 - xp]*(-10*(-1 + xp) + 
        4*(1 + xp)*Log[xp])) + Delz*(-36*(-1 + xp) + 10*(-1 + xp)*zeta2 + 
      2*(-1 + xp)*Li[2, 1 - xp] - 4*(1 + xp)*Li[3, 1 - xp] + 
      (25 - 3*xp - 4*(1 + xp)*zeta2 - 4*(1 + xp)*Li[2, 1 - xp])*Log[xp] + 
      ((17 - 15*xp)*Log[xp]^2)/2 + (5*(1 + xp)*Log[xp]^3)/3 + 
      Log[1 - xp]^2*(5 - 5*xp + 2*(1 + xp)*Log[xp]) + 
      Log[1 - xp]*(12*(-1 + xp) + 4*(1 + xp)*Li[2, 1 - xp] + 
        12*(-1 + xp)*Log[xp] - 4*(1 + xp)*Log[xp]^2)) + 
    Log[Q2/muF2]*((10*(-1 + xp)*(-1 + 2*zp))/zp - 
      (4*(1 + xp)*(-1 + 2*zp)*Log[xp])/zp + 
      Dzp[0]*(-10*(-1 + xp) + 4*(1 + xp)*Log[xp]) + 
      Delz*(12*(-1 + xp) + 4*(1 + xp)*Li[2, 1 - xp] + 12*(-1 + xp)*Log[xp] - 
        4*(1 + xp)*Log[xp]^2 + Log[1 - xp]*(-10*(-1 + xp) + 
          4*(1 + xp)*Log[xp]))) + (10*(-1 + xp)*(-1 + 2*zp)*Log[1 - zp])/zp + 
    (2*(-1 + xp)*(5 - 6*zp + 3*zp^2 + xp^2*(5 - 6*zp + 3*zp^2) - 
       2*xp*(-5 + 16*zp - 17*zp^2 + 8*zp^3))*Log[zp])/
     ((1 + xp^2 + xp*(2 - 4*zp))*(-1 + zp)*zp) - 
    (8*(1 + xp^4 + xp*(2 - 4*zp) + xp^3*(2 - 4*zp) + 2*xp^2*(1 - zp + zp^2))*
      Log[1 + xp - Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]]*
      Log[1 - xp + Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]])/
     (1 + xp^2 + xp*(2 - 4*zp))^(3/2) + 
    (8*(1 + xp^4 + xp*(2 - 4*zp) + xp^3*(2 - 4*zp) + 2*xp^2*(1 - zp + zp^2))*
      Log[1 + xp - Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]]*
      Log[-1 + xp + Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]])/
     (1 + xp^2 + xp*(2 - 4*zp))^(3/2) + 
    Log[2]*((8*(1 + xp^4 + xp*(2 - 4*zp) + xp^3*(2 - 4*zp) + 
         2*xp^2*(1 - zp + zp^2))*Log[1 - xp + Sqrt[1 + 2*xp + xp^2 - 
            4*xp*zp]])/(1 + xp^2 + xp*(2 - 4*zp))^(3/2) - 
      (8*(1 + xp^4 + xp*(2 - 4*zp) + xp^3*(2 - 4*zp) + 
         2*xp^2*(1 - zp + zp^2))*Log[-1 + xp + Sqrt[1 + 2*xp + xp^2 - 
            4*xp*zp]])/(1 + xp^2 + xp*(2 - 4*zp))^(3/2)) + 
    Log[xp]*((-2*(6 - 17*zp + xp^3*(-6 + 7*zp) + 
         xp^2*(-6 + 25*zp - 32*zp^2) + xp*(6 - 47*zp + 64*zp^2)))/
       ((1 + xp^2 + xp*(2 - 4*zp))*zp) + 
      (16*(1 + xp^4 + xp*(2 - 4*zp) + xp^3*(2 - 4*zp) + 
         2*xp^2*(1 - zp + zp^2))*Log[2])/(1 + xp^2 + xp*(2 - 4*zp))^(3/2) - 
      (4*(1 + xp)*(-1 + 2*zp)*Log[1 - zp])/zp - 
      (4*(1 + xp)*(-1 + zp)*Log[zp])/zp - 
      (8*(1 + xp^4 + xp*(2 - 4*zp) + xp^3*(2 - 4*zp) + 
         2*xp^2*(1 - zp + zp^2))*Log[-1 + xp + Sqrt[1 + 2*xp + xp^2 - 
            4*xp*zp]])/(1 + xp^2 + xp*(2 - 4*zp))^(3/2) - 
      (8*(1 + xp^4 + xp*(2 - 4*zp) + xp^3*(2 - 4*zp) + 
         2*xp^2*(1 - zp + zp^2))*Log[1 + xp + Sqrt[1 + 2*xp + xp^2 - 
            4*xp*zp]])/(1 + xp^2 + xp*(2 - 4*zp))^(3/2))), 
 "\!\(\*SubscriptBox[SuperscriptBox[\(G\), \((2)\)], \
\({1,qq',\([3]\),MS}\)]\)" -> 
  CF*(20 - 16*zp - (8*zeta2*(-1 + xp^2 + 4*zp + xp*(2 - 4*zp + 4*zp^2)))/
     (-1 + xp) + 16*Sqrt[xp]*Sqrt[zp]*InvTanInt[-(Sqrt[zp]/Sqrt[xp])] - 
    16*Sqrt[xp]*Sqrt[zp]*InvTanInt[Sqrt[zp]/Sqrt[xp]] - 
    16*Sqrt[xp]*Sqrt[zp]*InvTanInt[-(Sqrt[xp]*Sqrt[zp])] + 
    16*Sqrt[xp]*Sqrt[zp]*InvTanInt[Sqrt[xp]*Sqrt[zp]] - 
    (16*(1 + xp^2)*zp*Li[2, 1 - xp])/(-1 + xp) + 16*(1 + xp)*(-1 + 2*zp)*
     Li[2, -xp] - (4*(-1 + 2*xp)*(-1 + 2*zp)*Li[2, 1 - zp])/(-1 + xp) - 
    (8*(xp + xp^3 + 2*zp + 2*xp^2*zp + 4*xp*zp^2)*Li[2, -(zp/xp)])/
     (-1 + xp^2) - (8*(1 + 2*xp*zp + 2*xp^3*zp + xp^2*(1 + 4*zp^2))*
      Li[2, -(xp*zp)])/(-1 + xp^2) + 
    (8*xp*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*(1 - 5*zp + 5*zp^2))*
      Li[2, 1/2 - xp/2 - Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]/2])/
     (1 + xp^2 + xp*(2 - 4*zp))^(3/2) + 
    (8*xp*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*(1 - 5*zp + 5*zp^2))*
      Li[2, 1/2 + xp/2 - Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]/2])/
     (1 + xp^2 + xp*(2 - 4*zp))^(3/2) + 
    (8*xp*(-1 + 2*zp + xp^2*(-1 + 2*zp) - 2*xp*(1 - 5*zp + 5*zp^2))*
      Li[2, 1/2 - 1/(2*xp) - Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]/(2*xp)])/
     (1 + xp^2 + xp*(2 - 4*zp))^(3/2) + 
    (8*xp*(-1 + 2*zp + xp^2*(-1 + 2*zp) - 2*xp*(1 - 5*zp + 5*zp^2))*
      Li[2, 1/2 + 1/(2*xp) - Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]/(2*xp)])/
     (1 + xp^2 + xp*(2 - 4*zp))^(3/2) - 
    (4*(1 + 2*xp^2 + 2*zp + xp*(-2 + 4*zp))*
      Li[2, (2 - 2*xp)^(-1) + zp/(2 - 2*xp) - Sqrt[1 - 2*zp + 4*xp*zp + zp^2]/
         (2 - 2*xp)])/(-1 + xp) + 
    (4*(1 + 2*zp + 4*xp^2*zp + xp*(2 - 4*zp + 8*zp^2))*
      Li[2, (2 - 2*xp)^(-1) - (2*xp)/(2 - 2*xp) - zp/(2 - 2*xp) + 
        Sqrt[1 - 2*zp + 4*xp*zp + zp^2]/(2 - 2*xp)])/(-1 + xp) - 
    (8*(1 + 2*zp + 4*xp*zp^2 + xp^2*(1 + 2*zp))*
      Li[2, 1 - 1/(2*xp) + zp/(2*xp) - Sqrt[1 - 2*zp + 4*xp*zp + zp^2]/
         (2*xp)])/(-1 + xp) + ((2 - 4*xp + 4*xp^2 + 4*zp + 8*xp*zp)*Log[2]^2)/
     (1 - xp) + (4*xp*(2 - 4*zp + 4*zp^2 + xp*(-1 + 2*zp))*Log[1 - xp]^2)/
     (-1 + xp) + (4 + 2/(-1 + xp) + 8*xp - (1 + xp^2 + xp*(2 - 4*zp))^
       (-3/2) + (3*xp^2)/(1 + xp^2 + xp*(2 - 4*zp))^(3/2) + 
      1/Sqrt[1 + xp^2 + xp*(2 - 4*zp)] - 
      (10*xp)/Sqrt[1 + xp^2 + xp*(2 - 4*zp)] + (xp^3*(2 - 4*zp))/
       (1 + xp^2 + xp*(2 - 4*zp))^(3/2) - (4*(1 + 4*xp)*zp)/(-1 + xp) + 
      xp*(-24 + 20/Sqrt[1 + xp^2 + xp*(2 - 4*zp)])*zp)*Log[xp]^2 + 
    ((2 + 4*zp + 8*xp^2*zp + 4*xp*(1 - 2*zp + 4*zp^2))*Log[zp]^2)/(-1 + xp) + 
    Log[zp]*((4*(-((-1 + zp)*zp) + xp^4*(5 + 2*zp) + 
         xp^3*(8 - 25*zp - 4*zp^2) + xp*(-2 + zp - 10*zp^2 + 4*zp^3) + 
         xp^2*(1 + 9*zp + 3*zp^2 + 8*zp^3)))/((-1 + xp)*
        (1 + xp^2 + xp*(2 - 4*zp))*(xp - zp)) + 16*Sqrt[xp]*Sqrt[zp]*
       ArcTan[Sqrt[zp]/Sqrt[xp]] - 16*Sqrt[xp]*Sqrt[zp]*
       ArcTan[Sqrt[xp]*Sqrt[zp]] - (8*(xp + xp^3 + 2*zp + 2*xp^2*zp + 
         4*xp*zp^2)*Log[1 + zp/xp])/(-1 + xp^2) - 
      (8*(1 + 2*xp*zp + 2*xp^3*zp + xp^2*(1 + 4*zp^2))*Log[1 + xp*zp])/
       (-1 + xp^2)) + (8*xp*(1 + xp^2*(1 - 2*zp) - 2*zp + 
       2*xp*(1 - 5*zp + 5*zp^2))*
      Log[1 + xp - Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]]*
      Log[1 - xp + Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]])/
     (1 + xp^2 + xp*(2 - 4*zp))^(3/2) + 
    (8*xp*(-1 + 2*zp + xp^2*(-1 + 2*zp) - 2*xp*(1 - 5*zp + 5*zp^2))*
      Log[1 + xp - Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]]*
      Log[-1 + xp + Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]])/
     (1 + xp^2 + xp*(2 - 4*zp))^(3/2) - 
    (2*(1 + 2*zp + 4*xp^2*zp + xp*(2 - 4*zp + 8*zp^2))*
      Log[1 + zp - Sqrt[1 - 2*zp + 4*xp*zp + zp^2]]^2)/(-1 + xp) - 
    (2*(1 + 2*zp + 4*xp^2*zp + xp*(2 - 4*zp + 8*zp^2))*
      Log[1 - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]]^2)/(-1 + xp) + 
    Log[2]*((-4*(1 + 2*xp)*Sqrt[1 + (-2 + 4*xp)*zp + zp^2])/(-1 + xp) + 
      (8*xp*(2 - 4*zp + 4*zp^2 + xp*(-1 + 2*zp))*Log[1 - xp])/(-1 + xp) + 
      2*(-4 - 6/(-1 + xp) - (1 + xp^2 + xp*(2 - 4*zp))^(-3/2) + 
        (3*xp^2)/(1 + xp^2 + xp*(2 - 4*zp))^(3/2) + 
        1/Sqrt[1 + xp^2 + xp*(2 - 4*zp)] - 
        (10*xp)/Sqrt[1 + xp^2 + xp*(2 - 4*zp)] - (4*zp)/(-1 + xp) + 
        xp*(-8 + 20/Sqrt[1 + xp^2 + xp*(2 - 4*zp)])*zp - 
        (16*xp*zp^2)/(-1 + xp) - (2*xp^3*(-1 + 2*zp))/
         (1 + xp^2 + xp*(2 - 4*zp))^(3/2))*Log[xp] + 
      (8*xp*(-1 + 2*zp + xp^2*(-1 + 2*zp) - 2*xp*(1 - 5*zp + 5*zp^2))*
        Log[1 - xp + Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]])/
       (1 + xp^2 + xp*(2 - 4*zp))^(3/2) + 
      (8*xp*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*(1 - 5*zp + 5*zp^2))*
        Log[-1 + xp + Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]])/
       (1 + xp^2 + xp*(2 - 4*zp))^(3/2) - 
      (8*xp*(2 - 4*zp + 4*zp^2 + xp*(-1 + 2*zp))*
        Log[1 + zp - Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/(-1 + xp) + 
      (4*(1 + 2*zp + 4*xp^2*zp + xp*(2 - 4*zp + 8*zp^2))*
        Log[1 - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/(-1 + xp)) + 
    (4*(1 + 2*xp)*Sqrt[1 + (-2 + 4*xp)*zp + zp^2]*
      Log[1 - 2*xp - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/(-1 + xp) + 
    ((2 + 4*xp^2 + 4*zp + xp*(-4 + 8*zp))*
      Log[1 - 2*xp - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]]^2)/(-1 + xp) + 
    Log[1 + zp - Sqrt[1 - 2*zp + 4*xp*zp + zp^2]]*
     ((-4*(1 + 2*xp)*Sqrt[1 + (-2 + 4*xp)*zp + zp^2])/(-1 + xp) + 
      (4*(1 + 2*zp + 4*xp^2*zp + xp*(2 - 4*zp + 8*zp^2))*
        Log[1 - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/(-1 + xp) + 
      (8*xp*(2 - 4*zp + 4*zp^2 + xp*(-1 + 2*zp))*
        Log[1 - 2*xp - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/(-1 + xp)) + 
    Log[1 - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]]*
     ((4*(1 + 2*xp)*Sqrt[1 + (-2 + 4*xp)*zp + zp^2])/(-1 + xp) - 
      (4*(1 + 2*zp + 4*xp^2*zp + xp*(2 - 4*zp + 8*zp^2))*
        Log[1 - 2*xp - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/(-1 + xp)) + 
    Log[1 - xp]*((4*(1 + 2*xp^2 + 2*zp + xp*(-2 + 4*zp))*
        Log[1 + zp - Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/(-1 + xp) - 
      (4*(1 + 2*zp + 4*xp^2*zp + xp*(2 - 4*zp + 8*zp^2))*
        Log[1 - 2*xp - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/(-1 + xp)) + 
    Log[xp]*(2*(3 + 4/(-1 + xp) - (1 + xp^2 + xp*(2 - 4*zp))^(-1) + 
        xp^2*((1 + xp^2 + xp*(2 - 4*zp))^(-1) + 4/(xp - zp)) - 
        xp^3/(1 + 2*xp + xp^2 - 4*xp*zp) - 
        4/Sqrt[1 + (-2 + 4*xp)*zp + zp^2] - 
        6/((-1 + xp)*Sqrt[1 + (-2 + 4*xp)*zp + zp^2]) - 
        (2*(1 + 2*xp)*zp^2)/((-1 + xp)*Sqrt[1 + (-2 + 4*xp)*zp + zp^2]) + 
        xp*(-3 + 2/(-xp + zp) + (1 + 2*xp + xp^2 - 4*xp*zp)^(-1) - 
          (16*zp)/Sqrt[1 + (-2 + 4*xp)*zp + zp^2]) + 
        (2*zp*(1 - 2/Sqrt[1 + (-2 + 4*xp)*zp + zp^2] + 
           xp*(2 + 8/Sqrt[1 + (-2 + 4*xp)*zp + zp^2])))/(1 - xp)) - 
      16*Sqrt[xp]*Sqrt[zp]*ArcTan[Sqrt[zp]/Sqrt[xp]] - 
      16*Sqrt[xp]*Sqrt[zp]*ArcTan[Sqrt[xp]*Sqrt[zp]] + 
      16*(1 + xp)*(-1 + 2*zp)*Log[1 + xp] + 
      (16*xp*(-1 + xp + 2*zp - xp*zp + xp^2*zp - 2*zp^2 + 2*xp*zp^2)*Log[zp])/
       (-1 + xp^2) + (8*(xp + xp^3 + 2*zp + 2*xp^2*zp + 4*xp*zp^2)*
        Log[1 + zp/xp])/(-1 + xp^2) - 
      (8*(1 + 2*xp*zp + 2*xp^3*zp + xp^2*(1 + 4*zp^2))*Log[1 + xp*zp])/
       (-1 + xp^2) + (8*xp*(1 + xp^2*(1 - 2*zp) - 2*zp + 
         2*xp*(1 - 5*zp + 5*zp^2))*Log[-1 + xp + 
          Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]])/(1 + xp^2 + xp*(2 - 4*zp))^
        (3/2) + (8*xp*(1 + xp^2*(1 - 2*zp) - 2*zp + 2*xp*(1 - 5*zp + 5*zp^2))*
        Log[1 + xp + Sqrt[1 + 2*xp + xp^2 - 4*xp*zp]])/
       (1 + xp^2 + xp*(2 - 4*zp))^(3/2) - 
      (4*(1 + 2*zp + 4*xp^2*zp + xp*(2 - 4*zp + 8*zp^2))*
        Log[1 + zp - Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/(-1 + xp) + 
      (4*(1 + 2*zp + 4*xp^2*zp + xp*(2 - 4*zp + 8*zp^2))*
        Log[1 - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/(-1 + xp) + 
      (4*(1 + 2*zp + 4*xp^2*zp + xp*(2 - 4*zp + 8*zp^2))*
        Log[1 - 2*xp - zp + Sqrt[1 - 2*zp + 4*xp*zp + zp^2]])/(-1 + xp)) + 
    ((-4*(-1 + 2*xp)*(-1 + 2*zp)*Li[2, zp/xp])/(-1 + xp) - 
      (4*(-1 + 2*xp)*(-1 + 2*zp)*Log[xp]^2)/(-1 + xp) - 
      (4*(-1 + 2*xp)*(-1 + 2*zp)*Log[xp - zp]*Log[zp])/(-1 + xp) + 
      Log[xp]*((4*(-1 + 2*xp)*(-1 + 2*zp)*Log[xp - zp])/(-1 + xp) + 
        (4*(-1 + 2*xp)*(-1 + 2*zp)*Log[zp])/(-1 + xp)))*Theta[xp - zp] + 
    ((-8*(-1 + 2*xp)*zeta2*(-1 + 2*zp))/(-1 + xp) + 
      (4*(-1 + 2*xp)*(-1 + 2*zp)*Li[2, xp/zp])/(-1 + xp) - 
      (2*(-1 + 2*xp)*(-1 + 2*zp)*Log[xp]^2)/(-1 + xp) + 
      (2*(-1 + 2*xp)*(-1 + 2*zp)*Log[zp]^2)/(-1 + xp) + 
      (4*(-1 + 2*xp)*(-1 + 2*zp)*Log[xp]*Log[-xp + zp])/(-1 + xp) - 
      (4*(-1 + 2*xp)*(-1 + 2*zp)*Log[zp]*Log[-xp + zp])/(-1 + xp))*
     Theta[-xp + zp]), 
 "\!\(\*SubscriptBox[SuperscriptBox[\(G\), \((1)\)], \({1,qg,MS}\)]\)" -> 
  CF*(4 + xp*(4 - 2/zp) - 2/zp - 4*zp + (-4 + 4/zp + 2*zp)*Dxp[0] + 
    Delx*(-4 + 4/zp + 2*zp)*Log[Q2/muF2] + 
    Delx*(2*zp + (-4 + 4/zp + 2*zp)*Log[1 - zp] + (-4 + 4/zp + 2*zp)*
       Log[zp])), 
 "\!\(\*SubscriptBox[SuperscriptBox[\(G\), \((1)\)], \({1,gq,MS}\)]\)" -> 
  -(((-1 + 2*xp)*(-1 + 2*zp))/zp) + (-1 + 2*xp)*Dzp[0] + 
   Delz*(-1 + 2*xp)*Log[Q2/muF2] + Delz*(2 - 2*xp + (-1 + 2*xp)*Log[1 - xp] + 
     (1 - 2*xp)*Log[xp]), 
 "\!\(\*SubscriptBox[SuperscriptBox[\(G\), \((1)\)], \({1,qq,MS}\)]\)" -> 
  CF*(4*(xp + zp) - 2*(1 + xp)*Dzp[0] + Dxp[0]*(-2*(1 + zp) + 4*Dzp[0]) + 
    (Delz*(-2*(1 + xp) + 4*Dxp[0]) + Delx*(6*Delz - 2*(1 + zp) + 4*Dzp[0]))*
     Log[Q2/muF2] + Delz*(2 - 2*xp + 4*Dxp[1] - 2*(1 + xp)*Log[1 - xp] + 
      (2*(1 + xp^2)*Log[xp])/(-1 + xp)) + 
    Delx*(2 - 16*Delz - 2*zp + 4*Dzp[1] - 2*(1 + zp)*Log[1 - zp] - 
      (2*(1 + zp^2)*Log[zp])/(-1 + zp))), 
 "\!\(\*SubscriptBox[SuperscriptBox[\(G\), \((0)\)], \({1,qq,MS}\)]\)" -> 
  Delx*Delz}


MSb=%;

Zq2xpPS  = Cf*1/2*( + 16 + 24*Log[(xp)] + 8*Log[(xp)]^2 - 16*xp - 8*xp*Log[(xp)] + 4*xp*Log[(xp)]^2 ); 
Zq2xpNSp = ( Cf*Ca*( -592/9 -80/3*Log[(xp)] -4*Log[(xp)]^2 +592/9*xp + 8/3*xp*Log[(xp)] +4*xp*Log[(xp)]^2 +8*zeta2 -8*zeta2*xp)
              + Cf*1/2*nf*( +80/9 +16/3*Log[(xp)] -80/9*xp -16/3*xp*Log[(xp)])
              + Cf^2*( -16 -16*Log[(xp)] +16*Log[(xp)]*Log[(1-xp)] +16*xp -8*xp*Log[(xp)] -16*xp*Log[(xp)]*Log[(1-xp)]) );
Zq2xpNSm = ((CF^2-1/2*CA*CF)*( 8*(1+xp)*( 4*Li[2,(-xp)] +4*Log[(xp)]*Log[(1+xp)] +2*zeta2 -Log[(xp)]^2 -3*Log[(xp)]) -56*(1-xp) ));

Zq2xpNS = Zq2xpNSp - Zq2xpNSm ;


Print[" ---------------------------------------------------------------------- "];
Print[" ---------------------------------------------------------------------- "];
Print[" For Reading the G1 Coefficient Functions in MSbar scheme: "]

(* NNLO *)
G12qgMS  =  "\!\(\*SubscriptBox[SuperscriptBox[\(G\), \((2)\)], \({1,qg,MS}\)]\)";
G12gqMS  =  "\!\(\*SubscriptBox[SuperscriptBox[\(G\), \((2)\)], \({1,gq,MS}\)]\)";
G12qqNSMS="\!\(\*SubscriptBox[SuperscriptBox[\(G\), \((2)\)], \({1,qq,NS,MS}\)]\)";
G12qqPSMS="\!\(\*SubscriptBox[SuperscriptBox[\(G\), \((2)\)], \({1,qq,NS,MS}\)]\)";
G12qQMS  = "\!\(\*SubscriptBox[SuperscriptBox[\(G\), \((2)\)], \({1,q\!\(\*OverscriptBox[\(q\), \(_\)]\),MS}\)]\)";
G12ggMS  = "\!\(\*SubscriptBox[SuperscriptBox[\(G\), \((2)\)], \({1,gg,MS}\)]\)";
G12qqp1MS="\!\(\*SubscriptBox[SuperscriptBox[\(G\), \((2)\)], \({1,qq',\([1]\),MS}\)]\)";
G12qqp2MS="\!\(\*SubscriptBox[SuperscriptBox[\(G\), \((2)\)], \({1,qq',\([2]\),MS}\)]\)";
G12qqp3MS="\!\(\*SubscriptBox[SuperscriptBox[\(G\), \((2)\)], \({1,qq',\([3]\),MS}\)]\)";
(* NLO *)
G11qgMS = "\!\(\*SubscriptBox[SuperscriptBox[\(G\), \((1)\)], \({1,qg,MS}\)]\)";
G11gqMS = "\!\(\*SubscriptBox[SuperscriptBox[\(G\), \((1)\)], \({1,gq,MS}\)]\)";
G11qqMS = "\!\(\*SubscriptBox[SuperscriptBox[\(G\), \((1)\)], \({1,qq,MS}\)]\)";
(* LO *)
G10qqMS = "\!\(\*SubscriptBox[SuperscriptBox[\(G\), \((0)\)], \({1,qq,MS}\)]\)";


(* NNLO qg channel in MSbar, "Subscript[G^(2), {1,qg,MS}]" *)
G12qgMS/.MSb
(* NLO qq channel in MSbar, "Subscript[G^(1), {1,qq,MS}]" *)
G11qqMS/.MSb



Print[" ---------------------------------------------------------------------- "];
