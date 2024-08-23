*************************************************************************************
** Polarized SIDIS coefficient functions up to NNLO from:                            
**                                                                                   
** Polarized semi-inclusive deep-inelastic scattering at NNLO in QCD                 
** L. Bonino, T. Gehrmann, M. Loechner, K. Schoenwald and G. Stagnitto               
**                                                                                   
** FORM readable format                                                              
**                                                                                   
** Notation, according to eq.(8) of the paper:                                       
**           DC[order][a2b][label] with                                              
**           - order: 1 = NLO, 2 = NNLO                                              
**           - a2b: means a -> b, for a and b partons                                
**           - label (it can be none, NS, PS, 1, 2, 3)                               
**                                                                                   
** Symbols:                                                                          
**          NC = 3: number of colours                                                
**          NF = 5: number of active flavours                                        
**                                                                                   
** Scales:                                                                           
**         LMUR = ln(muR^2/Q2)                                                       
**         LMUF = ln(muF^2/Q2)                                                       
**         LMUA = ln(muA^2/Q2)                                                       
**         with Q2 = -q2, invariant mass of the photon                               
**              muR: renormalisation scale                                           
**              muF: initial-state factorisation scale                               
**              muA: final-state factorisation scale                                 
**                                                                                   
** Functions:                                                                        
**            Li2(a) = PolyLog(2,a)                                                  
**            Li3(a) = PolyLog(3,a)                                                  
**            sqrtxz1 = sqrt(1 - 2*z + z*z + 4*x*z)                                  
**            poly2   = 1 + 2*x + x*x - 4*x*z                                        
**            sqrtxz2 = sqrt(poly2)                                                  
**            sqrtxz3 = sqrt(x/z)                                                    
**            InvTanInt(x) = int_0^x dt arctan(t)/t : Arctangent integral            
**            T(region): Heaviside Theta function                                    
**                                                                                   
** Distributions:                                                                    
**           Dd([1-x]) is the Dirac delta function of argument [1-x]                 
**           Dn(a,[1-x]) = (ln^a(1-x)/(1-x))_+ (plus-prescription) for a = 0,1,2,3   
**           same for z                                                              
**                                                                                   
** Kinematic regions in (x,z)-plane: as defined in 2201.06982                        
**           ui = Ui for i = 1,2,3,4                                                 
**           ri = Ri, ti = Ti for i = 1,2                                            
**           Ri, Ti and Ui defined in eq. (5.9), (5.12) and (5.16)                   
**                                                                                   
** Constants: pi, zeta3 = Zeta(3) with Zeta Riemann Zeta function                    
**                                                                                   
*************************************************************************************

id DC1Q2Q      = (  + NC^-1 * (  - z
       - x
       )

       + NC * ( z
       + x
       )

       + Dd([1-x])*NC^-1 * (  - 1/2
       + 1/2*z
       + 1/2*ln(z)
       - ln(z)*[1-z]^-1
       + 1/2*ln(z)*z
       + 1/2*ln([1-z])
       + 1/2*ln([1-z])*z
       )

       + Dd([1-x])*NC * ( 1/2
       - 1/2*z
       - 1/2*ln(z)
       + ln(z)*[1-z]^-1
       - 1/2*ln(z)*z
       - 1/2*ln([1-z])
       - 1/2*ln([1-z])*z
       )

       + Dd([1-x])*LMUA*NC^-1 * (  - 1/2
       - 1/2*z
       )

       + Dd([1-x])*LMUA*NC * ( 1/2
       + 1/2*z
       )

       + Dd([1-x])*Dd([1-z])*NC^-1 * ( 4
       )

       + Dd([1-x])*Dd([1-z])*NC * (  - 4
       )

       + Dd([1-x])*Dd([1-z])*LMUF*NC^-1 * ( 3/4
       )

       + Dd([1-x])*Dd([1-z])*LMUF*NC * (  - 3/4
       )

       + Dd([1-x])*Dd([1-z])*LMUA*NC^-1 * ( 3/4
       )

       + Dd([1-x])*Dd([1-z])*LMUA*NC * (  - 3/4
       )

       + Dd([1-x])*Dn(0,[1-z])*LMUA*NC^-1 * ( 1
       )

       + Dd([1-x])*Dn(0,[1-z])*LMUA*NC * (  - 1
       )

       + Dd([1-x])*Dn(1,[1-z])*NC^-1 * (  - 1
       )

       + Dd([1-x])*Dn(1,[1-z])*NC * ( 1
       )

       + Dd([1-z])*NC^-1 * (  - 1/2
       + 1/2*x
       - 1/2*ln(x)
       + ln(x)*[1-x]^-1
       - 1/2*ln(x)*x
       + 1/2*ln([1-x])
       + 1/2*ln([1-x])*x
       )

       + Dd([1-z])*NC * ( 1/2
       - 1/2*x
       + 1/2*ln(x)
       - ln(x)*[1-x]^-1
       + 1/2*ln(x)*x
       - 1/2*ln([1-x])
       - 1/2*ln([1-x])*x
       )

       + Dd([1-z])*LMUF*NC^-1 * (  - 1/2
       - 1/2*x
       )

       + Dd([1-z])*LMUF*NC * ( 1/2
       + 1/2*x
       )

       + Dd([1-z])*Dn(0,[1-x])*LMUF*NC^-1 * ( 1
       )

       + Dd([1-z])*Dn(0,[1-x])*LMUF*NC * (  - 1
       )

       + Dd([1-z])*Dn(1,[1-x])*NC^-1 * (  - 1
       )

       + Dd([1-z])*Dn(1,[1-x])*NC * ( 1
       )

       + Dn(0,[1-x])*NC^-1 * ( 1/2
       + 1/2*z
       )

       + Dn(0,[1-x])*NC * (  - 1/2
       - 1/2*z
       )

       + Dn(0,[1-x])*Dn(0,[1-z])*NC^-1 * (  - 1
       )

       + Dn(0,[1-x])*Dn(0,[1-z])*NC * ( 1
       )

       + Dn(0,[1-z])*NC^-1 * ( 1/2
       + 1/2*x
       )

       + Dn(0,[1-z])*NC * (  - 1/2
       - 1/2*x ) );

id DC1Q2G      = (  + NC^-1 * (  - 1
       + 1/2*z^-1
       + z
       + 1/2*x*z^-1
       - x
       )

       + NC * ( 1
       - 1/2*z^-1
       - z
       - 1/2*x*z^-1
       + x
       )

       + Dd([1-x])*NC^-1 * (  - 1/2*z
       + ln(z)
       - ln(z)*z^-1
       - 1/2*ln(z)*z
       + ln([1-z])
       - ln([1-z])*z^-1
       - 1/2*ln([1-z])*z
       )

       + Dd([1-x])*NC * ( 1/2*z
       - ln(z)
       + ln(z)*z^-1
       + 1/2*ln(z)*z
       - ln([1-z])
       + ln([1-z])*z^-1
       + 1/2*ln([1-z])*z
       )

       + Dd([1-x])*LMUA*NC^-1 * (  - 1
       + z^-1
       + 1/2*z
       )

       + Dd([1-x])*LMUA*NC * ( 1
       - z^-1
       - 1/2*z
       )

       + Dn(0,[1-x])*NC^-1 * ( 1
       - z^-1
       - 1/2*z
       )

       + Dn(0,[1-x])*NC * (  - 1
       + z^-1
       + 1/2*z ) );

id DC1G2Q      = (  + Dd([1-z])*LMUF * ( 1/2
       - x
       )

       + Dd([1-z]) * ( 1
       - x
       + 1/2*ln(x)
       - ln(x)*x
       - 1/2*ln([1-x])
       + ln([1-x])*x
       )

       + Dn(0,[1-z]) * (  - 1/2
       + x
       )

       + 1
       - 1/2*z^-1
       + x*z^-1
       - 2*x );

id DC2Q2G      = (  + NC^-2 * ( 41/8
       + 4*[1-x]^-1*z^-1
       + 2*[1-x]^-1*z^-1*rln2^2
       + 3*[1-x]^-1*z^-1*sqrtxz1*rln2
       - 4*[1-x]^-1
       - 2*[1-x]^-1*rln2^2
       - [1-x]^-1*sqrtxz1*rln2
       + 6*[1-x]^-1*z
       + [1-x]^-1*z*rln2^2
       - 6*z^-1
       - z^-1*rln2^2
       - z^-1*sqrtxz1*rln2
       + [1-x-z]^-1
       - sqrtxz1*rln2
       - 2*z*[1-x-z]^-1
       - 13/4*z
       + z^2*[1-x-z]^-1
       + 9/2*x*z^-1
       - x*z^-1*rln2^2
       - 2*x*z^-1*sqrtxz1*rln2
       - 49/8*x
       + 1/4*x*z
       - 4/3*pi^2*[1-x]^-1*z^-1
       + 5/3*pi^2*[1-x]^-1
       - 5/6*pi^2*[1-x]^-1*z
       + 29/24*pi^2*z^-1
       - 25/12*pi^2
       + 5/4*pi^2*z
       + 13/24*pi^2*x*z^-1
       - 7/12*pi^2*x
       + 1/12*pi^2*x*z
       - 4*ln(1 + sqrtxz1 - z)*[1-x]^-1*z^-1*rln2
       - 3*ln(1 + sqrtxz1 - z)*[1-x]^-1*z^-1*sqrtxz1
       + 4*ln(1 + sqrtxz1 - z)*[1-x]^-1*rln2
       + ln(1 + sqrtxz1 - z)*[1-x]^-1*sqrtxz1
       - 2*ln(1 + sqrtxz1 - z)*[1-x]^-1*z*rln2
       + 2*ln(1 + sqrtxz1 - z)*z^-1*rln2
       + ln(1 + sqrtxz1 - z)*z^-1*sqrtxz1
       + ln(1 + sqrtxz1 - z)*sqrtxz1
       + 2*ln(1 + sqrtxz1 - z)*x*z^-1*rln2
       + 2*ln(1 + sqrtxz1 - z)*x*z^-1*sqrtxz1
       + 2*ln(1 + sqrtxz1 - z)^2*[1-x]^-1*z^-1
       - 2*ln(1 + sqrtxz1 - z)^2*[1-x]^-1
       + ln(1 + sqrtxz1 - z)^2*[1-x]^-1*z
       - ln(1 + sqrtxz1 - z)^2*z^-1
       - ln(1 + sqrtxz1 - z)^2*x*z^-1
       - 55/2*ln(x)
       - 18*ln(x)*[1-x]^-1*z^-1
       + 3*ln(x)*[1-x]^-1*z^-1*rln2
       + 3/2*ln(x)*[1-x]^-1*z^-1*sqrtxz1
       + 24*ln(x)*[1-x]^-1
       - 3*ln(x)*[1-x]^-1*rln2
       - 1/2*ln(x)*[1-x]^-1*sqrtxz1
       - 123/8*ln(x)*[1-x]^-1*z
       + 3/2*ln(x)*[1-x]^-1*z*rln2
       + 145/8*ln(x)*z^-1
       - 3/2*ln(x)*z^-1*rln2
       - 1/2*ln(x)*z^-1*sqrtxz1
       + ln(x)*[1-x-z]^-1
       - 1/2*ln(x)*sqrtxz1
       - 3*ln(x)*z*[1-x-z]^-2
       + 4*ln(x)*z*[1-x-z]^-1
       + 61/4*ln(x)*z
       + 3*ln(x)*z^2*[1-x-z]^-2
       - 2*ln(x)*z^2*[1-x-z]^-1
       - 7/8*ln(x)*x*z^-1
       - 3/2*ln(x)*x*z^-1*rln2
       - ln(x)*x*z^-1*sqrtxz1
       + 3*ln(x)*x*[1-x-z]^-2
       - 6*ln(x)*x*[1-x-z]^-1
       - 2*ln(x)*x
       + 1/4*ln(x)*x*z
       - 3*ln(x)*x^2*[1-x-z]^-2
       + 3*ln(x)*x^2*[1-x-z]^-1
       + ln(x)*x^3*[1-x-z]^-2
       - 2*ln(x)*ln(1 + sqrtxz1 - z)*[1-x]^-1*z^-1
       + 2*ln(x)*ln(1 + sqrtxz1 - z)*[1-x]^-1
       - ln(x)*ln(1 + sqrtxz1 - z)*[1-x]^-1*z
       + ln(x)*ln(1 + sqrtxz1 - z)*z^-1
       + ln(x)*ln(1 + sqrtxz1 - z)*x*z^-1
       - ln(x)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^-1
       + ln(x)*ln(1 + sqrtxz1 + z)*[1-x]^-1
       - 1/2*ln(x)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z
       + 1/2*ln(x)*ln(1 + sqrtxz1 + z)*z^-1
       + 1/2*ln(x)*ln(1 + sqrtxz1 + z)*x*z^-1
       + 23/2*ln(x)^2
       + 17/2*ln(x)^2*[1-x]^-1*z^-1
       - 12*ln(x)^2*[1-x]^-1
       + 6*ln(x)^2*[1-x]^-1*z
       - 25/4*ln(x)^2*z^-1
       - 31/4*ln(x)^2*z
       - 11/4*ln(x)^2*x*z^-1
       + 9/2*ln(x)^2*x
       - 1/4*ln(x)^2*x*z
       - 15*ln(x)*ln([1-x])
       - 10*ln(x)*ln([1-x])*[1-x]^-1*z^-1
       + 13*ln(x)*ln([1-x])*[1-x]^-1
       - 13/2*ln(x)*ln([1-x])*[1-x]^-1*z
       + 15/2*ln(x)*ln([1-x])*z^-1
       + 19/2*ln(x)*ln([1-x])*z
       + 5/2*ln(x)*ln([1-x])*x*z^-1
       - 5*ln(x)*ln([1-x])*x
       + 1/2*ln(x)*ln([1-x])*x*z
       - 10*ln(x)*ln(z)
       - 8*ln(x)*ln(z)*[1-x]^-1*z^-1
       + 25/2*ln(x)*ln(z)*[1-x]^-1
       - 25/4*ln(x)*ln(z)*[1-x]^-1*z
       + 11/2*ln(x)*ln(z)*z^-1
       + 17/2*ln(x)*ln(z)*z
       + 5/2*ln(x)*ln(z)*x*z^-1
       - 5*ln(x)*ln(z)*x
       - 33/2*ln(x)*ln([1-z])
       - 12*ln(x)*ln([1-z])*[1-x]^-1*z^-1
       + 17*ln(x)*ln([1-z])*[1-x]^-1
       - 17/2*ln(x)*ln([1-z])*[1-x]^-1*z
       + 35/4*ln(x)*ln([1-z])*z^-1
       + 11*ln(x)*ln([1-z])*z
       + 15/4*ln(x)*ln([1-z])*x*z^-1
       - 13/2*ln(x)*ln([1-z])*x
       + 1/2*ln(x)*ln([1-z])*x*z
       + 23/2*ln([1-x])
       + 9*ln([1-x])*[1-x]^-1*z^-1
       - 11*ln([1-x])*[1-x]^-1
       + 7*ln([1-x])*[1-x]^-1*z
       - 65/8*ln([1-x])*z^-1
       - 8*ln([1-x])*z
       + 11/8*ln([1-x])*x*z^-1
       + ln([1-x])*x
       - 1/4*ln([1-x])*x*z
       + 11/2*ln([1-x])^2
       + 5/2*ln([1-x])^2*[1-x]^-1*z^-1
       - 3*ln([1-x])^2*[1-x]^-1
       + 3/2*ln([1-x])^2*[1-x]^-1*z
       - 3*ln([1-x])^2*z^-1
       - 13/4*ln([1-x])^2*z
       - ln([1-x])^2*x*z^-1
       + 3/2*ln([1-x])^2*x
       - 1/4*ln([1-x])^2*x*z
       + 6*ln([1-x])*ln(z)
       + 4*ln([1-x])*ln(z)*[1-x]^-1*z^-1
       - 6*ln([1-x])*ln(z)*[1-x]^-1
       + 3*ln([1-x])*ln(z)*[1-x]^-1*z
       - 7/2*ln([1-x])*ln(z)*z^-1
       - 9/2*ln([1-x])*ln(z)*z
       - 3/2*ln([1-x])*ln(z)*x*z^-1
       + 3*ln([1-x])*ln(z)*x
       + 10*ln([1-x])*ln([1-z])
       + 6*ln([1-x])*ln([1-z])*[1-x]^-1*z^-1
       - 8*ln([1-x])*ln([1-z])*[1-x]^-1
       + 4*ln([1-x])*ln([1-z])*[1-x]^-1*z
       - 11/2*ln([1-x])*ln([1-z])*z^-1
       - 13/2*ln([1-x])*ln([1-z])*z
       - 5/2*ln([1-x])*ln([1-z])*x*z^-1
       + 4*ln([1-x])*ln([1-z])*x
       - 1/2*ln([1-x])*ln([1-z])*x*z
       + 17*ln(z)
       + 15/2*ln(z)*[1-x]^-1*z^-1
       + ln(z)*[1-x]^-1*z^-1*rln2
       + 3/2*ln(z)*[1-x]^-1*z^-1*sqrtxz1
       - 14*ln(z)*[1-x]^-1
       - ln(z)*[1-x]^-1*rln2
       - 1/2*ln(z)*[1-x]^-1*sqrtxz1
       + 17/2*ln(z)*[1-x]^-1*z
       + 1/2*ln(z)*[1-x]^-1*z*rln2
       - 15/2*ln(z)*z^-1
       - 1/2*ln(z)*z^-1*rln2
       - 1/2*ln(z)*z^-1*sqrtxz1
       - 1/2*ln(z)*sqrtxz1
       - 10*ln(z)*z
       + 3/2*ln(z)*x*z^-1
       - 1/2*ln(z)*x*z^-1*rln2
       - ln(z)*x*z^-1*sqrtxz1
       + 5/2*ln(z)*x
       - ln(z)*ln(1 + sqrtxz1 - z)*[1-x]^-1*z^-1
       + ln(z)*ln(1 + sqrtxz1 - z)*[1-x]^-1
       - 1/2*ln(z)*ln(1 + sqrtxz1 - z)*[1-x]^-1*z
       + 1/2*ln(z)*ln(1 + sqrtxz1 - z)*z^-1
       + 1/2*ln(z)*ln(1 + sqrtxz1 - z)*x*z^-1
       + 2*ln(z)^2
       + 3/2*ln(z)^2*[1-x]^-1*z^-1
       - 3*ln(z)^2*[1-x]^-1
       + 3/2*ln(z)^2*[1-x]^-1*z
       - ln(z)^2*z^-1
       - 9/4*ln(z)^2*z
       - 1/2*ln(z)^2*x*z^-1
       + 2*ln(z)^2*x
       + 1/4*ln(z)^2*x*z
       + 6*ln(z)*ln([1-z])
       + 6*ln(z)*ln([1-z])*[1-x]^-1*z^-1
       - 10*ln(z)*ln([1-z])*[1-x]^-1
       + 5*ln(z)*ln([1-z])*[1-x]^-1*z
       - 3*ln(z)*ln([1-z])*z^-1
       - 5*ln(z)*ln([1-z])*z
       - 2*ln(z)*ln([1-z])*x*z^-1
       + 4*ln(z)*ln([1-z])*x
       + 17*ln([1-z])
       + 12*ln([1-z])*[1-x]^-1*z^-1
       - 16*ln([1-z])*[1-x]^-1
       + 10*ln([1-z])*[1-x]^-1*z
       - 89/8*ln([1-z])*z^-1
       - ln([1-z])*[1-x-z]^-1
       + 3*ln([1-z])*z*[1-x-z]^-2
       - 4*ln([1-z])*z*[1-x-z]^-1
       - 33/4*ln([1-z])*z
       - 3*ln([1-z])*z^2*[1-x-z]^-2
       + 2*ln([1-z])*z^2*[1-x-z]^-1
       + 11/8*ln([1-z])*x*z^-1
       - 3*ln([1-z])*x*[1-x-z]^-2
       + 6*ln([1-z])*x*[1-x-z]^-1
       + 1/2*ln([1-z])*x
       - ln([1-z])*x*z
       + 3*ln([1-z])*x^2*[1-x-z]^-2
       - 3*ln([1-z])*x^2*[1-x-z]^-1
       - ln([1-z])*x^3*[1-x-z]^-2
       + 13/2*ln([1-z])^2
       + 4*ln([1-z])^2*[1-x]^-1*z^-1
       - 6*ln([1-z])^2*[1-x]^-1
       + 3*ln([1-z])^2*[1-x]^-1*z
       - 7/2*ln([1-z])^2*z^-1
       - 17/4*ln([1-z])^2*z
       - 3/2*ln([1-z])^2*x*z^-1
       + 5/2*ln([1-z])^2*x
       - 1/4*ln([1-z])^2*x*z
       - Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*z^-1
       + Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1
       - 1/2*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*z
       + 1/2*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z^-1
       + 1/2*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x*z^-1
       + Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*z^-1
       - Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1
       + 1/2*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*z
       - 1/2*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z^-1
       - 1/2*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x*z^-1
       - Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*[1-x]^-1*z^-1
       + Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*[1-x]^-1
       - 1/2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*[1-x]^-1*z
       + 1/2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*z^-1
       + 1/2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*x*z^-1
       + Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*[1-x]^-1*z^-1
       - Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*[1-x]^-1
       + 1/2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*[1-x]^-1*z
       - 1/2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*z^-1
       - 1/2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*x*z^-1
       - Li2(1 - x*z^-1)
       - Li2(1 - x*z^-1)*[1-x]^-1*z^-1
       + 2*Li2(1 - x*z^-1)*[1-x]^-1
       - Li2(1 - x*z^-1)*[1-x]^-1*z
       + 1/2*Li2(1 - x*z^-1)*z^-1
       + Li2(1 - x*z^-1)*z
       + 1/2*Li2(1 - x*z^-1)*x*z^-1
       - Li2(1 - x*z^-1)*x
       + 1/2*Li2(x)
       + Li2(x)*[1-x]^-1*z^-1
       - Li2(x)*[1-x]^-1
       + 1/2*Li2(x)*[1-x]^-1*z
       - 3/4*Li2(x)*z^-1
       - 1/2*Li2(x)*z
       - 3/4*Li2(x)*x*z^-1
       + 1/2*Li2(x)*x
       - Li2(z)
       + 2*Li2(z)*[1-x]^-1*z^-1
       - 3*Li2(z)*[1-x]^-1
       + 3/2*Li2(z)*[1-x]^-1*z
       + Li2(z)*z^-1
       + 1/2*Li2(z)*z
       )

       + NC^-1*NF * ( 1/6*ln(x*[1-x]*z*[1-z])*[1-x]^-1*z^-1
       - 1/6*ln(x)*[1-x]^-1*z^-1
       - 1/6*ln([1-x])*[1-x]^-1*z^-1
       - 1/6*ln(z)*[1-x]^-1*z^-1
       - 1/6*ln([1-z])*[1-x]^-1*z^-1
       )

       + NC*NF * (  - 1/6*ln(x*[1-x]*z*[1-z])*[1-x]^-1*z^-1
       + 1/6*ln(x)*[1-x]^-1*z^-1
       + 1/6*ln([1-x])*[1-x]^-1*z^-1
       + 1/6*ln(z)*[1-x]^-1*z^-1
       + 1/6*ln([1-z])*[1-x]^-1*z^-1
       )

       + NC^2 * ( 35/24
       + 2*[1-x]^-1*z^-1
       - 2*[1-x]^-1
       - 4*[1-x]^-1*rln2^2
       - 2*[1-x]^-1*sqrtxz1*rln2
       + 5/2*[1-x]^-1*z
       + 17/36*z^-1
       + 2*rln2^2
       + 2*sqrtxz1*rln2
       - 23/6*z
       + 2*z*rln2^2
       + 41/18*z^2
       + 83/36*x*z^-1
       - 163/24*x
       + 2*x*rln2^2
       + 8/3*x*z
       - 19/18*x*z^2
       - pi^2*[1-x]^-1*z^-1
       + 2/3*pi^2*[1-x]^-1
       - 1/2*pi^2*[1-x]^-1*z
       + 17/24*pi^2*z^-1
       - 5/12*pi^2
       + 5/4*pi^2*z
       + 17/24*pi^2*x*z^-1
       - 1/4*pi^2*x
       + 1/12*pi^2*x*z
       - 2*ln(1 + sqrtxz1 - z)*[1-x]^-1*z^-1*rln2
       + 6*ln(1 + sqrtxz1 - z)*[1-x]^-1*rln2
       + 2*ln(1 + sqrtxz1 - z)*[1-x]^-1*sqrtxz1
       - ln(1 + sqrtxz1 - z)*[1-x]^-1*z*rln2
       + ln(1 + sqrtxz1 - z)*z^-1*rln2
       - 2*ln(1 + sqrtxz1 - z)*rln2
       - 2*ln(1 + sqrtxz1 - z)*sqrtxz1
       - 2*ln(1 + sqrtxz1 - z)*z*rln2
       + ln(1 + sqrtxz1 - z)*x*z^-1*rln2
       - 2*ln(1 + sqrtxz1 - z)*x*rln2
       + 2*ln(1 + sqrtxz1 - z)^2*[1-x]^-1*z^-1
       - 2*ln(1 + sqrtxz1 - z)^2*[1-x]^-1
       + ln(1 + sqrtxz1 - z)^2*[1-x]^-1*z
       - ln(1 + sqrtxz1 - z)^2*z^-1
       - ln(1 + sqrtxz1 - z)^2*x*z^-1
       + 2*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)
       - 2*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^-1
       - 2*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*[1-x]^-1
       - ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z
       + ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*z^-1
       + 2*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*z
       + ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*x*z^-1
       + 2*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*x
       + 2*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^-1*rln2
       + 2*ln(1 + sqrtxz1 + z)*[1-x]^-1*rln2
       + ln(1 + sqrtxz1 + z)*[1-x]^-1*z*rln2
       - ln(1 + sqrtxz1 + z)*z^-1*rln2
       - 2*ln(1 + sqrtxz1 + z)*rln2
       - 2*ln(1 + sqrtxz1 + z)*z*rln2
       - ln(1 + sqrtxz1 + z)*x*z^-1*rln2
       - 2*ln(1 + sqrtxz1 + z)*x*rln2
       - 13/2*ln(x)
       - 1/6*ln(x)*[1-x]^-1*z^-1
       + 2*ln(x)*[1-x]^-1*z^-1*rln2
       + 4*ln(x)*[1-x]^-1
       - 4*ln(x)*[1-x]^-1*rln2
       - ln(x)*[1-x]^-1*sqrtxz1
       - 75/8*ln(x)*[1-x]^-1*z
       + ln(x)*[1-x]^-1*z*rln2
       - 4/3*ln(x)*[1-x]^-1*z^2
       + 95/24*ln(x)*z^-1
       - ln(x)*z^-1*rln2
       + ln(x)*rln2
       + ln(x)*sqrtxz1
       + 43/4*ln(x)*z
       + ln(x)*z*rln2
       + 8/3*ln(x)*z^2
       - 181/24*ln(x)*x*z^-1
       - ln(x)*x*z^-1*rln2
       - 5/2*ln(x)*x
       + ln(x)*x*rln2
       + 13/4*ln(x)*x*z
       - 4/3*ln(x)*x*z^2
       - 2*ln(x)*ln(1 + sqrtxz1 - z)*[1-x]^-1*z^-1
       + 2*ln(x)*ln(1 + sqrtxz1 - z)*[1-x]^-1
       - ln(x)*ln(1 + sqrtxz1 - z)*[1-x]^-1*z
       + ln(x)*ln(1 + sqrtxz1 - z)*z^-1
       + ln(x)*ln(1 + sqrtxz1 - z)*x*z^-1
       - ln(x)*ln(1 + sqrtxz1 + z)
       + 2*ln(x)*ln(1 + sqrtxz1 + z)*[1-x]^-1
       - ln(x)*ln(1 + sqrtxz1 + z)*z
       - ln(x)*ln(1 + sqrtxz1 + z)*x
       + 9/2*ln(x)^2
       + 5*ln(x)^2*[1-x]^-1*z^-1
       - 5*ln(x)^2*[1-x]^-1
       + 5/2*ln(x)^2*[1-x]^-1*z
       - 11/4*ln(x)^2*z^-1
       - 17/4*ln(x)^2*z
       - 11/4*ln(x)^2*x*z^-1
       + 9/2*ln(x)^2*x
       - 1/4*ln(x)^2*x*z
       - 8*ln(x)*ln([1-x])
       - 8*ln(x)*ln([1-x])*[1-x]^-1*z^-1
       + 8*ln(x)*ln([1-x])*[1-x]^-1
       - 4*ln(x)*ln([1-x])*[1-x]^-1*z
       + 4*ln(x)*ln([1-x])*z^-1
       + 15/2*ln(x)*ln([1-x])*z
       + 4*ln(x)*ln([1-x])*x*z^-1
       - 8*ln(x)*ln([1-x])*x
       + 1/2*ln(x)*ln([1-x])*x*z
       - 8*ln(x)*ln(z)
       - 3*ln(x)*ln(z)*[1-x]^-1*z^-1
       + 19/2*ln(x)*ln(z)*[1-x]^-1
       + 7/4*ln(x)*ln(z)*[1-x]^-1*z
       + 3/2*ln(x)*ln(z)*z^-1
       - 1/2*ln(x)*ln(z)*z
       + 3/2*ln(x)*ln(z)*x*z^-1
       - 9*ln(x)*ln(z)*x
       - 13/2*ln(x)*ln([1-z])
       - 7*ln(x)*ln([1-z])*[1-x]^-1*z^-1
       + 7*ln(x)*ln([1-z])*[1-x]^-1
       - 7/2*ln(x)*ln([1-z])*[1-x]^-1*z
       + 15/4*ln(x)*ln([1-z])*z^-1
       + 6*ln(x)*ln([1-z])*z
       + 15/4*ln(x)*ln([1-z])*x*z^-1
       - 13/2*ln(x)*ln([1-z])*x
       + 1/2*ln(x)*ln([1-z])*x*z
       + 2*ln([1-x])
       + 6*ln([1-x])*[1-x]^-1*z^-1
       - 6*ln([1-x])*[1-x]^-1
       + 4*ln([1-x])*[1-x]^-1*z
       - 55/24*ln([1-x])*z^-1
       - 11/2*ln([1-x])*z
       - 4/3*ln([1-x])*z^2
       + 125/24*ln([1-x])*x*z^-1
       + 3/2*ln([1-x])*x
       - 7/4*ln([1-x])*x*z
       + 2/3*ln([1-x])*x*z^2
       + 3*ln([1-x])^2
       + 2*ln([1-x])^2*[1-x]^-1*z^-1
       - 2*ln([1-x])^2*[1-x]^-1
       + ln([1-x])^2*[1-x]^-1*z
       - 7/4*ln([1-x])^2*z^-1
       - 11/4*ln([1-x])^2*z
       - 7/4*ln([1-x])^2*x*z^-1
       + 3*ln([1-x])^2*x
       - 1/4*ln([1-x])^2*x*z
       + 6*ln([1-x])*ln(z)
       + 4*ln([1-x])*ln(z)*[1-x]^-1*z^-1
       - 4*ln([1-x])*ln(z)*[1-x]^-1
       + 2*ln([1-x])*ln(z)*[1-x]^-1*z
       - 3/2*ln([1-x])*ln(z)*z^-1
       - 1/2*ln([1-x])*ln(z)*z
       - 3/2*ln([1-x])*ln(z)*x*z^-1
       + 7*ln([1-x])*ln(z)*x
       + 5*ln([1-x])*ln([1-z])
       + 2*ln([1-x])*ln([1-z])*[1-x]^-1*z^-1
       - 2*ln([1-x])*ln([1-z])*[1-x]^-1
       + ln([1-x])*ln([1-z])*[1-x]^-1*z
       - 3*ln([1-x])*ln([1-z])*z^-1
       - 9/2*ln([1-x])*ln([1-z])*z
       - 3*ln([1-x])*ln([1-z])*x*z^-1
       + 5*ln([1-x])*ln([1-z])*x
       - 1/2*ln([1-x])*ln([1-z])*x*z
       + 8*ln(z)
       + 6*ln(z)*[1-x]^-1*z^-1
       - 2*ln(z)*[1-x]^-1*z^-1*rln2
       - 5*ln(z)*[1-x]^-1
       - 4*ln(z)*[1-x]^-1*rln2
       - ln(z)*[1-x]^-1*sqrtxz1
       + 5*ln(z)*[1-x]^-1*z
       - ln(z)*[1-x]^-1*z*rln2
       - 7/6*ln(z)*z^-1
       + ln(z)*z^-1*rln2
       + 3*ln(z)*rln2
       + ln(z)*sqrtxz1
       - 8*ln(z)*z
       + 3*ln(z)*z*rln2
       - 4/3*ln(z)*z^2
       + 19/3*ln(z)*x*z^-1
       + ln(z)*x*z^-1*rln2
       + 6*ln(z)*x
       + 3*ln(z)*x*rln2
       - 3/2*ln(z)*x*z
       + 2/3*ln(z)*x*z^2
       - ln(z)*ln(1 + sqrtxz1 - z)
       + 2*ln(z)*ln(1 + sqrtxz1 - z)*[1-x]^-1
       - ln(z)*ln(1 + sqrtxz1 - z)*z
       - ln(z)*ln(1 + sqrtxz1 - z)*x
       - 2*ln(z)*ln(1 + sqrtxz1 + z)
       + 2*ln(z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^-1
       + 2*ln(z)*ln(1 + sqrtxz1 + z)*[1-x]^-1
       + ln(z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z
       - ln(z)*ln(1 + sqrtxz1 + z)*z^-1
       - 2*ln(z)*ln(1 + sqrtxz1 + z)*z
       - ln(z)*ln(1 + sqrtxz1 + z)*x*z^-1
       - 2*ln(z)*ln(1 + sqrtxz1 + z)*x
       - 2*ln(z)*ln(1 + z)*[1-x]^-1*z^-1
       - 2*ln(z)*ln(1 + z)*[1-x]^-1
       - ln(z)*ln(1 + z)*[1-x]^-1*z
       + 7/2*ln(z)^2
       + 1/2*ln(z)^2*[1-x]^-1*z^-1
       - 3/2*ln(z)^2*[1-x]^-1
       + 1/4*ln(z)^2*[1-x]^-1*z
       + 5/4*ln(z)^2*z^-1
       + 15/4*ln(z)^2*z
       + 5/4*ln(z)^2*x*z^-1
       + 9/2*ln(z)^2*x
       + 1/4*ln(z)^2*x*z
       + 5*ln(z)*ln([1-z])
       + 2*ln(z)*ln([1-z])*[1-x]^-1*z^-1
       - 2*ln(z)*ln([1-z])*[1-x]^-1
       + ln(z)*ln([1-z])*[1-x]^-1*z
       - 5/2*ln(z)*ln([1-z])*z^-1
       - 5*ln(z)*ln([1-z])*z
       - 5/2*ln(z)*ln([1-z])*x*z^-1
       + 5*ln(z)*ln([1-z])*x
       - 1/2*ln([1-z])
       + 3*ln([1-z])*[1-x]^-1*z^-1
       - 3*ln([1-z])*[1-x]^-1
       + 2*ln([1-z])*[1-x]^-1*z
       + 17/24*ln([1-z])*z^-1
       - 7/4*ln([1-z])*z
       - 4/3*ln([1-z])*z^2
       + 125/24*ln([1-z])*x*z^-1
       - 3/2*ln([1-z])*x
       - 5/2*ln([1-z])*x*z
       + 2/3*ln([1-z])*x*z^2
       + 2*ln([1-z])^2
       + 1/2*ln([1-z])^2*[1-x]^-1*z^-1
       - 1/2*ln([1-z])^2*[1-x]^-1
       + 1/4*ln([1-z])^2*[1-x]^-1*z
       - 5/4*ln([1-z])^2*z^-1
       - 7/4*ln([1-z])^2*z
       - 5/4*ln([1-z])^2*x*z^-1
       + 2*ln([1-z])^2*x
       - 1/4*ln([1-z])^2*x*z
       + Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)
       - 2*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*z^-1
       - Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*z
       + Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z^-1
       + Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z
       + Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x*z^-1
       + Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x
       - Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)
       + 2*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*z^-1
       + Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*z
       - Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z^-1
       - Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z
       - Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x*z^-1
       - Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x
       - Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)
       + 2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*[1-x]^-1
       - Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*z
       - Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*x
       + Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)
       - 2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*[1-x]^-1
       + Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*z
       + Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*x
       - Li2(1 - x*z^-1)
       - Li2(1 - x*z^-1)*[1-x]^-1*z^-1
       + Li2(1 - x*z^-1)*[1-x]^-1
       - 1/2*Li2(1 - x*z^-1)*[1-x]^-1*z
       + 1/2*Li2(1 - x*z^-1)*z^-1
       + Li2(1 - x*z^-1)*z
       + 1/2*Li2(1 - x*z^-1)*x*z^-1
       - Li2(1 - x*z^-1)*x
       - 2*Li2( - z)*[1-x]^-1*z^-1
       - 2*Li2( - z)*[1-x]^-1
       - Li2( - z)*[1-x]^-1*z
       - 1/2*Li2(x)
       - 1/4*Li2(x)*z^-1
       + 1/2*Li2(x)*z
       - 1/4*Li2(x)*x*z^-1
       - 1/2*Li2(x)*x
       - 3/2*Li2(z)*z^-1
       - 11/2*Li2(z)*z
       - 3/2*Li2(z)*x*z^-1
       - Li2(z)*x
       )

       + LMUR*NC^-1*NF * ( 1/3
       - 1/6*z^-1
       - 1/3*z
       - 1/6*x*z^-1
       + 1/3*x
       )

       + LMUR * (  - 11/6
       + 11/12*z^-1
       + 11/6*z
       + 11/12*x*z^-1
       - 11/6*x
       )

       + LMUR*NC*NF * (  - 1/3
       + 1/6*z^-1
       + 1/3*z
       + 1/6*x*z^-1
       - 1/3*x
       )

       + LMUR*NC^2 * ( 11/6
       - 11/12*z^-1
       - 11/6*z
       - 11/12*x*z^-1
       + 11/6*x
       )

       + LMUF*NC^-2 * (  - 5/4
       + 7/8*z^-1
       + 3/4*z
       - 1/8*x*z^-1
       - 1/4*x
       + 1/2*x*z
       + ln(x)
       + ln(x)*[1-x]^-1*z^-1
       - ln(x)*[1-x]^-1
       + 1/2*ln(x)*[1-x]^-1*z
       - 3/4*ln(x)*z^-1
       - 3/4*ln(x)*z
       - 3/4*ln(x)*x*z^-1
       + ln(x)*x
       - 1/4*ln(x)*x*z
       - 3/2*ln([1-x])
       + ln([1-x])*z^-1
       + 5/4*ln([1-x])*z
       + ln([1-x])*x*z^-1
       - 3/2*ln([1-x])*x
       + 1/4*ln([1-x])*x*z
       - 1/2*ln(z)
       + 1/2*ln(z)*z^-1
       + 1/4*ln(z)*z
       + 1/2*ln(z)*x*z^-1
       - 1/2*ln(z)*x
       + 1/4*ln(z)*x*z
       - 1/2*ln([1-z])
       + 1/2*ln([1-z])*z^-1
       + 1/4*ln([1-z])*z
       + 1/2*ln([1-z])*x*z^-1
       - 1/2*ln([1-z])*x
       + 1/4*ln([1-z])*x*z
       )

       + LMUF * ( 5/2
       - 7/4*z^-1
       - 3/2*z
       + 1/4*x*z^-1
       + 1/2*x
       - x*z
       - 2*ln(x)
       - 2*ln(x)*[1-x]^-1*z^-1
       + 2*ln(x)*[1-x]^-1
       - ln(x)*[1-x]^-1*z
       + 3/2*ln(x)*z^-1
       + 3/2*ln(x)*z
       + 3/2*ln(x)*x*z^-1
       - 2*ln(x)*x
       + 1/2*ln(x)*x*z
       + 3*ln([1-x])
       - 2*ln([1-x])*z^-1
       - 5/2*ln([1-x])*z
       - 2*ln([1-x])*x*z^-1
       + 3*ln([1-x])*x
       - 1/2*ln([1-x])*x*z
       + ln(z)
       - ln(z)*z^-1
       - 1/2*ln(z)*z
       - ln(z)*x*z^-1
       + ln(z)*x
       - 1/2*ln(z)*x*z
       + ln([1-z])
       - ln([1-z])*z^-1
       - 1/2*ln([1-z])*z
       - ln([1-z])*x*z^-1
       + ln([1-z])*x
       - 1/2*ln([1-z])*x*z
       )

       + LMUF*NC^2 * (  - 5/4
       + 7/8*z^-1
       + 3/4*z
       - 1/8*x*z^-1
       - 1/4*x
       + 1/2*x*z
       + ln(x)
       + ln(x)*[1-x]^-1*z^-1
       - ln(x)*[1-x]^-1
       + 1/2*ln(x)*[1-x]^-1*z
       - 3/4*ln(x)*z^-1
       - 3/4*ln(x)*z
       - 3/4*ln(x)*x*z^-1
       + ln(x)*x
       - 1/4*ln(x)*x*z
       - 3/2*ln([1-x])
       + ln([1-x])*z^-1
       + 5/4*ln([1-x])*z
       + ln([1-x])*x*z^-1
       - 3/2*ln([1-x])*x
       + 1/4*ln([1-x])*x*z
       - 1/2*ln(z)
       + 1/2*ln(z)*z^-1
       + 1/4*ln(z)*z
       + 1/2*ln(z)*x*z^-1
       - 1/2*ln(z)*x
       + 1/4*ln(z)*x*z
       - 1/2*ln([1-z])
       + 1/2*ln([1-z])*z^-1
       + 1/4*ln([1-z])*z
       + 1/2*ln([1-z])*x*z^-1
       - 1/2*ln([1-z])*x
       + 1/4*ln([1-z])*x*z
       )

       + LMUA*NC^-2 * ( 7/4
       - z^-1
       - z
       - 1/2*x*z^-1
       + 1/4*x
       + 1/2*x*z
       + 1/2*ln(x)
       + ln(x)*[1-x]^-1*z^-1
       - ln(x)*[1-x]^-1
       + 1/2*ln(x)*[1-x]^-1*z
       - 1/2*ln(x)*z^-1
       - 1/4*ln(x)*z
       - 1/2*ln(x)*x*z^-1
       + 1/2*ln(x)*x
       - 1/4*ln(x)*x*z
       - 1/2*ln([1-x])
       + 1/2*ln([1-x])*z^-1
       + 1/4*ln([1-x])*z
       + 1/2*ln([1-x])*x*z^-1
       - 1/2*ln([1-x])*x
       + 1/4*ln([1-x])*x*z
       + 1/2*ln(z)
       + 1/4*ln(z)*z
       - 1/2*ln(z)*x
       - 1/4*ln(z)*x*z
       - 1/2*ln([1-z])
       + 1/2*ln([1-z])*z^-1
       + 1/4*ln([1-z])*z
       + 1/2*ln([1-z])*x*z^-1
       - 1/2*ln([1-z])*x
       + 1/4*ln([1-z])*x*z
       )

       + LMUA*NC^-1*NF * (  - 1/3
       + 1/6*z^-1
       + 1/3*z
       + 1/6*x*z^-1
       - 1/3*x
       )

       + LMUA * (  - 8/3
       + 47/12*z^-1
       - 1/3*z
       - 4/3*z^2
       + 47/12*x*z^-1
       - 5/3*x
       - 5/2*x*z
       + 2/3*x*z^2
       - ln(x)
       - 2*ln(x)*[1-x]^-1*z^-1
       + 2*ln(x)*[1-x]^-1
       - ln(x)*[1-x]^-1*z
       + ln(x)*z^-1
       + 1/2*ln(x)*z
       + ln(x)*x*z^-1
       - ln(x)*x
       + 1/2*ln(x)*x*z
       + ln([1-x])
       - ln([1-x])*z^-1
       - 1/2*ln([1-x])*z
       - ln([1-x])*x*z^-1
       + ln([1-x])*x
       - 1/2*ln([1-x])*x*z
       + ln(z)
       + ln(z)*z^-1
       + 7/2*ln(z)*z
       + ln(z)*x*z^-1
       + 3*ln(z)*x
       + 1/2*ln(z)*x*z
       + 3*ln([1-z])
       - 2*ln([1-z])*z^-1
       - 5/2*ln([1-z])*z
       - 2*ln([1-z])*x*z^-1
       + 3*ln([1-z])*x
       - 1/2*ln([1-z])*x*z
       )

       + LMUA*NC*NF * ( 1/3
       - 1/6*z^-1
       - 1/3*z
       - 1/6*x*z^-1
       + 1/3*x
       )

       + LMUA*NC^2 * ( 11/12
       - 35/12*z^-1
       + 4/3*z
       + 4/3*z^2
       - 41/12*x*z^-1
       + 17/12*x
       + 2*x*z
       - 2/3*x*z^2
       + 1/2*ln(x)
       + ln(x)*[1-x]^-1*z^-1
       - ln(x)*[1-x]^-1
       + 1/2*ln(x)*[1-x]^-1*z
       - 1/2*ln(x)*z^-1
       - 1/4*ln(x)*z
       - 1/2*ln(x)*x*z^-1
       + 1/2*ln(x)*x
       - 1/4*ln(x)*x*z
       - 1/2*ln([1-x])
       + 1/2*ln([1-x])*z^-1
       + 1/4*ln([1-x])*z
       + 1/2*ln([1-x])*x*z^-1
       - 1/2*ln([1-x])*x
       + 1/4*ln([1-x])*x*z
       - 3/2*ln(z)
       - ln(z)*z^-1
       - 15/4*ln(z)*z
       - ln(z)*x*z^-1
       - 5/2*ln(z)*x
       - 1/4*ln(z)*x*z
       - 5/2*ln([1-z])
       + 3/2*ln([1-z])*z^-1
       + 9/4*ln([1-z])*z
       + 3/2*ln([1-z])*x*z^-1
       - 5/2*ln([1-z])*x
       + 1/4*ln([1-z])*x*z
       )

       + LMUA*LMUF*NC^-2 * ( 1/2
       - 1/2*z^-1
       - 1/4*z
       - 1/2*x*z^-1
       + 1/2*x
       - 1/4*x*z
       )

       + LMUA*LMUF * (  - 1
       + z^-1
       + 1/2*z
       + x*z^-1
       - x
       + 1/2*x*z
       )

       + LMUA*LMUF*NC^2 * ( 1/2
       - 1/2*z^-1
       - 1/4*z
       - 1/2*x*z^-1
       + 1/2*x
       - 1/4*x*z
       )

       + Dd([1-x])*NC^-2 * (  - 3/2
       - 13/16*z
       + 5/2*zeta3*z^-1
       - 3/2*zeta3
       + 3/4*zeta3*z
       + 1/8*pi^2*z^-1
       - 1/3*pi^2
       + 1/16*pi^2*z
       + 67/16*ln(z)
       - 9/4*ln(z)*z^-1
       - 45/16*ln(z)*z
       - 1/6*ln(z)*pi^2*z^-1
       + 1/4*ln(z)*pi^2
       - 1/8*ln(z)*pi^2*z
       - 1/2*ln(z)^2
       + 13/32*ln(z)^2*z
       + 5/24*ln(z)^3
       - 5/48*ln(z)^3*z
       + ln(z)^2*ln([1-z])
       - ln(z)^2*ln([1-z])*z^-1
       - 1/2*ln(z)^2*ln([1-z])*z
       - 3/4*ln(z)*ln([1-z])^2
       + 3/8*ln(z)*ln([1-z])^2*z
       + 3*ln(z)*Li2(z)
       - 3*ln(z)*Li2(z)*z^-1
       - 3/2*ln(z)*Li2(z)*z
       + 9/4*ln([1-z])
       - 9/4*ln([1-z])*z^-1
       - 5/4*ln([1-z])*z
       - 1/3*ln([1-z])*pi^2*z^-1
       + 1/2*ln([1-z])*pi^2
       - 1/4*ln([1-z])*pi^2*z
       + ln([1-z])^2
       - 3/4*ln([1-z])^2*z^-1
       - 1/16*ln([1-z])^2*z
       - 5/12*ln([1-z])^3
       + 5/12*ln([1-z])^3*z^-1
       + 5/24*ln([1-z])^3*z
       - 3/2*ln([1-z])*Li2(1 - z)
       + ln([1-z])*Li2(1 - z)*z^-1
       + 3/4*ln([1-z])*Li2(1 - z)*z
       + 1/2*ln([1-z])*Li2(z)
       - 3/2*ln([1-z])*Li2(z)*z^-1
       - 1/4*ln([1-z])*Li2(z)*z
       + Li3(1 - z)
       - 3/2*Li3(1 - z)*z^-1
       - 1/2*Li3(1 - z)*z
       - 4*Li3(z)
       + 3*Li3(z)*z^-1
       + 2*Li3(z)*z
       + 3/4*Li2(z)*z^-1
       - 1/2*Li2(z)*z
       )

       + Dd([1-x])*NC^2 * (  - 295/36
       + 1169/108*z^-1
       + 2/3*z^-1*rln2^3
       + 2/3*rln2^3
       - 389/144*z
       + 1/3*z*rln2^3
       - 269/108*z^2
       - 7*zeta3*z^-1
       - 6*zeta3*z
       - 1/8*pi^2*z^-1
       - 1/3*pi^2*z^-1*rln2
       - 5/6*pi^2
       - 1/3*pi^2*rln2
       + 7/48*pi^2*z
       - 1/6*pi^2*z*rln2
       - ln(1 + z)*z^-1*rln2^2
       - ln(1 + z)*rln2^2
       - 1/2*ln(1 + z)*z*rln2^2
       + 1/6*ln(1 + z)*pi^2*z^-1
       + 1/6*ln(1 + z)*pi^2
       + 1/12*ln(1 + z)*pi^2*z
       + 247/16*ln(z)
       + 71/36*ln(z)*z^-1
       + 31/16*ln(z)*z
       + 31/18*ln(z)*z^2
       - 1/2*ln(z)*pi^2*z^-1
       - 1/12*ln(z)*pi^2
       - 11/24*ln(z)*pi^2*z
       + 1/2*ln(z)*ln(1 + z)*z
       - 5/2*ln(z)^2
       - 10/3*ln(z)^2*z^-1
       - 7/32*ln(z)^2*z
       + 3/2*ln(z)^2*ln(1 + z)
       + 3/2*ln(z)^2*ln(1 + z)*z^-1
       + 3/4*ln(z)^2*ln(1 + z)*z
       - 5/8*ln(z)^3
       - 5/3*ln(z)^3*z^-1
       - 65/48*ln(z)^3*z
       - 1/2*ln(z)^2*ln([1-z])
       + 1/2*ln(z)^2*ln([1-z])*z^-1
       + 1/4*ln(z)^2*ln([1-z])*z
       + 3/2*ln(z)*ln([1-z])*z
       + 2*ln(z)*ln([1-z])*ln(1 + z)
       + 2*ln(z)*ln([1-z])*ln(1 + z)*z^-1
       + ln(z)*ln([1-z])*ln(1 + z)*z
       + 15/4*ln(z)*ln([1-z])^2
       + 3/2*ln(z)*ln([1-z])^2*z^-1
       + 21/8*ln(z)*ln([1-z])^2*z
       + ln(z)*Li2( - z)
       + ln(z)*Li2( - z)*z^-1
       + 1/2*ln(z)*Li2( - z)*z
       + 2*ln(z)*Li2(z)
       - 2*ln(z)*Li2(z)*z^-1
       - ln(z)*Li2(z)*z
       + 23/6*ln([1-z])
       - 55/36*ln([1-z])*z^-1
       - ln([1-z])*z^-1*rln2^2
       - ln([1-z])*rln2^2
       - 31/12*ln([1-z])*z
       - 1/2*ln([1-z])*z*rln2^2
       - 13/18*ln([1-z])*z^2
       - 1/3*ln([1-z])*pi^2*z^-1
       - 1/2*ln([1-z])*pi^2
       - 7/12*ln([1-z])*pi^2*z
       + 2*ln([1-z])*ln(1 + z)*z^-1*rln2
       + 2*ln([1-z])*ln(1 + z)*rln2
       + ln([1-z])*ln(1 + z)*z*rln2
       - ln([1-z])*ln(1 + z)^2
       - ln([1-z])*ln(1 + z)^2*z^-1
       - 1/2*ln([1-z])*ln(1 + z)^2*z
       + 3*ln([1-z])^2
       - 10/3*ln([1-z])^2*z^-1
       - 1/16*ln([1-z])^2*z
       + 1/3*ln([1-z])^2*z^2
       - 7/12*ln([1-z])^3
       + 7/12*ln([1-z])^3*z^-1
       + 7/24*ln([1-z])^3*z
       + 9/2*ln([1-z])*Li2(1 - z)
       - ln([1-z])*Li2(1 - z)*z^-1
       + 3/4*ln([1-z])*Li2(1 - z)*z
       + 2*ln([1-z])*Li2( - z)
       + 2*ln([1-z])*Li2( - z)*z^-1
       + ln([1-z])*Li2( - z)*z
       + 11/2*ln([1-z])*Li2(z)
       + 3/2*ln([1-z])*Li2(z)*z^-1
       + 13/4*ln([1-z])*Li2(z)*z
       + 1/3*ln([1+z])*pi^2*z^-1
       + 1/3*ln([1+z])*pi^2
       + 1/6*ln([1+z])*pi^2*z
       - 2*Li3(1/2 - 1/2*z)
       - 2*Li3(1/2 - 1/2*z)*z^-1
       - Li3(1/2 - 1/2*z)*z
       - 2*Li3(1/2 + 1/2*z)
       - 2*Li3(1/2 + 1/2*z)*z^-1
       - Li3(1/2 + 1/2*z)*z
       + 4*Li3(1 - z)
       + 7/2*Li3(1 - z)*z^-1
       + 3*Li3(1 - z)*z
       + Li3( - z)
       + Li3( - z)*z^-1
       + 1/2*Li3( - z)*z
       + 2*Li3(1/(1 + z))
       + 2*Li3(1/(1 + z))*z^-1
       + Li3(1/(1 + z))*z
       - 2*Li3(1/(1 + z) - 1/(1 + z)*z)
       - 2*Li3(1/(1 + z) - 1/(1 + z)*z)*z^-1
       - Li3(1/(1 + z) - 1/(1 + z)*z)*z
       - Li3(z)
       + 10*Li3(z)*z^-1
       + 15/2*Li3(z)*z
       + 1/2*Li2( - z)*z
       - Li2(z)
       + 89/12*Li2(z)*z^-1
       + 1/2*Li2(z)*z
       - 2/3*Li2(z)*z^2
       )

       + Dd([1-x])*LMUR*NC^-1*NF * ( 1/6*z
       - 1/3*ln(z)
       + 1/3*ln(z)*z^-1
       + 1/6*ln(z)*z
       - 1/3*ln([1-z])
       + 1/3*ln([1-z])*z^-1
       + 1/6*ln([1-z])*z
       )

       + Dd([1-x])*LMUR * (  - 11/12*z
       + 11/6*ln(z)
       - 11/6*ln(z)*z^-1
       - 11/12*ln(z)*z
       + 11/6*ln([1-z])
       - 11/6*ln([1-z])*z^-1
       - 11/12*ln([1-z])*z
       )

       + Dd([1-x])*LMUR*NC*NF * (  - 1/6*z
       + 1/3*ln(z)
       - 1/3*ln(z)*z^-1
       - 1/6*ln(z)*z
       + 1/3*ln([1-z])
       - 1/3*ln([1-z])*z^-1
       - 1/6*ln([1-z])*z
       )

       + Dd([1-x])*LMUR*NC^2 * ( 11/12*z
       - 11/6*ln(z)
       + 11/6*ln(z)*z^-1
       + 11/12*ln(z)*z
       - 11/6*ln([1-z])
       + 11/6*ln([1-z])*z^-1
       + 11/12*ln([1-z])*z
       )

       + Dd([1-x])*LMUF*NC^-2 * (  - 3/8*z
       + 1/6*pi^2*z^-1
       - 1/6*pi^2
       + 1/12*pi^2*z
       + 3/4*ln(z)
       - 3/4*ln(z)*z^-1
       - 3/8*ln(z)*z
       + 3/4*ln([1-z])
       - 3/4*ln([1-z])*z^-1
       - 3/8*ln([1-z])*z
       )

       + Dd([1-x])*LMUF * ( 3/4*z
       - 1/3*pi^2*z^-1
       + 1/3*pi^2
       - 1/6*pi^2*z
       - 3/2*ln(z)
       + 3/2*ln(z)*z^-1
       + 3/4*ln(z)*z
       - 3/2*ln([1-z])
       + 3/2*ln([1-z])*z^-1
       + 3/4*ln([1-z])*z
       )

       + Dd([1-x])*LMUF*NC^2 * (  - 3/8*z
       + 1/6*pi^2*z^-1
       - 1/6*pi^2
       + 1/12*pi^2*z
       + 3/4*ln(z)
       - 3/4*ln(z)*z^-1
       - 3/8*ln(z)*z
       + 3/4*ln([1-z])
       - 3/4*ln([1-z])*z^-1
       - 3/8*ln([1-z])*z
       )

       + Dd([1-x])*LMUA*NC^-2 * (  - 17/8
       + 9/4*z^-1
       + 7/8*z
       + 1/6*pi^2*z^-1
       - 1/4*pi^2
       + 1/8*pi^2*z
       + 3/4*ln(z)
       - 3/8*ln(z)*z
       - 1/2*ln(z)^2
       + 1/4*ln(z)^2*z
       - ln(z)*ln([1-z])
       + ln(z)*ln([1-z])*z^-1
       + 1/2*ln(z)*ln([1-z])*z
       - 5/4*ln([1-z])
       + 3/4*ln([1-z])*z^-1
       + ln([1-z])^2
       - ln([1-z])^2*z^-1
       - 1/2*ln([1-z])^2*z
       - 5/2*Li2(z)
       + 3*Li2(z)*z^-1
       + 5/4*Li2(z)*z
       )

       + Dd([1-x])*LMUA*NC^-1*NF * (  - 1/6*z
       + 1/3*ln(z)
       - 1/3*ln(z)*z^-1
       - 1/6*ln(z)*z
       + 1/3*ln([1-z])
       - 1/3*ln([1-z])*z^-1
       - 1/6*ln([1-z])*z
       )

       + Dd([1-x])*LMUA * ( 73/12
       - 34/9*z^-1
       - 13/6*z
       - 13/18*z^2
       - 1/3*pi^2*z^-1
       + 1/6*pi^2
       - 5/12*pi^2*z
       - 7/3*ln(z)
       - 29/6*ln(z)*z^-1
       + 8/3*ln(z)*z
       + 2/3*ln(z)*z^2
       + 2*ln(z)*ln(1 + z)
       + 2*ln(z)*ln(1 + z)*z^-1
       + ln(z)*ln(1 + z)*z
       - ln(z)^2
       - 3*ln(z)^2*z^-1
       - 3*ln(z)^2*z
       - 2*ln(z)*ln([1-z])
       + 2*ln(z)*ln([1-z])*z^-1
       + ln(z)*ln([1-z])*z
       + 14/3*ln([1-z])
       - 29/6*ln([1-z])*z^-1
       + 17/12*ln([1-z])*z
       + 2/3*ln([1-z])*z^2
       - 3*ln([1-z])^2
       + 3*ln([1-z])^2*z^-1
       + 3/2*ln([1-z])^2*z
       + 2*Li2( - z)
       + 2*Li2( - z)*z^-1
       + Li2( - z)*z
       + 3*Li2(z)
       + 3/2*Li2(z)*z
       )

       + Dd([1-x])*LMUA*NC*NF * ( 1/6*z
       - 1/3*ln(z)
       + 1/3*ln(z)*z^-1
       + 1/6*ln(z)*z
       - 1/3*ln([1-z])
       + 1/3*ln([1-z])*z^-1
       + 1/6*ln([1-z])*z
       )

       + Dd([1-x])*LMUA*NC^2 * (  - 95/24
       + 55/36*z^-1
       + 31/24*z
       + 13/18*z^2
       + 1/6*pi^2*z^-1
       + 1/12*pi^2
       + 7/24*pi^2*z
       + 19/12*ln(z)
       + 29/6*ln(z)*z^-1
       - 55/24*ln(z)*z
       - 2/3*ln(z)*z^2
       - 2*ln(z)*ln(1 + z)
       - 2*ln(z)*ln(1 + z)*z^-1
       - ln(z)*ln(1 + z)*z
       + 3/2*ln(z)^2
       + 3*ln(z)^2*z^-1
       + 11/4*ln(z)^2*z
       + 3*ln(z)*ln([1-z])
       - 3*ln(z)*ln([1-z])*z^-1
       - 3/2*ln(z)*ln([1-z])*z
       - 41/12*ln([1-z])
       + 49/12*ln([1-z])*z^-1
       - 17/12*ln([1-z])*z
       - 2/3*ln([1-z])*z^2
       + 2*ln([1-z])^2
       - 2*ln([1-z])^2*z^-1
       - ln([1-z])^2*z
       - 2*Li2( - z)
       - 2*Li2( - z)*z^-1
       - Li2( - z)*z
       - 1/2*Li2(z)
       - 3*Li2(z)*z^-1
       - 11/4*Li2(z)*z
       )

       + Dd([1-x])*LMUA*LMUR*NC^-1*NF * ( 1/3
       - 1/3*z^-1
       - 1/6*z
       )

       + Dd([1-x])*LMUA*LMUR * (  - 11/6
       + 11/6*z^-1
       + 11/12*z
       )

       + Dd([1-x])*LMUA*LMUR*NC*NF * (  - 1/3
       + 1/3*z^-1
       + 1/6*z
       )

       + Dd([1-x])*LMUA*LMUR*NC^2 * ( 11/6
       - 11/6*z^-1
       - 11/12*z
       )

       + Dd([1-x])*LMUA*LMUF*NC^-2 * (  - 3/4
       + 3/4*z^-1
       + 3/8*z
       )

       + Dd([1-x])*LMUA*LMUF * ( 3/2
       - 3/2*z^-1
       - 3/4*z
       )

       + Dd([1-x])*LMUA*LMUF*NC^2 * (  - 3/4
       + 3/4*z^-1
       + 3/8*z
       )

       + Dd([1-x])*LMUA^2*NC^-2 * ( 1/4
       - 1/16*z
       + 1/4*ln(z)
       - 1/8*ln(z)*z
       - 1/2*ln([1-z])
       + 1/2*ln([1-z])*z^-1
       + 1/4*ln([1-z])*z
       )

       + Dd([1-x])*LMUA^2*NC^-1*NF * (  - 1/3
       + 1/3*z^-1
       + 1/6*z
       )

       + Dd([1-x])*LMUA^2 * (  - 2/3
       + 3/4*z^-1
       - 25/24*z
       - 1/3*z^2
       + 1/2*ln(z)
       + ln(z)*z^-1
       + 5/4*ln(z)*z
       + 2*ln([1-z])
       - 2*ln([1-z])*z^-1
       - ln([1-z])*z
       )

       + Dd([1-x])*LMUA^2*NC*NF * ( 1/3
       - 1/3*z^-1
       - 1/6*z
       )

       + Dd([1-x])*LMUA^2*NC^2 * ( 5/12
       - 3/4*z^-1
       + 53/48*z
       + 1/3*z^2
       - 3/4*ln(z)
       - ln(z)*z^-1
       - 9/8*ln(z)*z
       - 3/2*ln([1-z])
       + 3/2*ln([1-z])*z^-1
       + 3/4*ln([1-z])*z
       )

       + Dd([1-x]) * ( 349/36
       - 1169/108*z^-1
       - 2/3*z^-1*rln2^3
       - 2/3*rln2^3
       + 253/72*z
       - 1/3*z*rln2^3
       + 269/108*z^2
       + 9/2*zeta3*z^-1
       + 3/2*zeta3
       + 21/4*zeta3*z
       + 1/3*pi^2*z^-1*rln2
       + 7/6*pi^2
       + 1/3*pi^2*rln2
       - 5/24*pi^2*z
       + 1/6*pi^2*z*rln2
       + ln(1 + z)*z^-1*rln2^2
       + ln(1 + z)*rln2^2
       + 1/2*ln(1 + z)*z*rln2^2
       - 1/6*ln(1 + z)*pi^2*z^-1
       - 1/6*ln(1 + z)*pi^2
       - 1/12*ln(1 + z)*pi^2*z
       - 157/8*ln(z)
       + 5/18*ln(z)*z^-1
       + 7/8*ln(z)*z
       - 31/18*ln(z)*z^2
       + 2/3*ln(z)*pi^2*z^-1
       - 1/6*ln(z)*pi^2
       + 7/12*ln(z)*pi^2*z
       - 1/2*ln(z)*ln(1 + z)*z
       + 3*ln(z)^2
       + 10/3*ln(z)^2*z^-1
       - 3/16*ln(z)^2*z
       - 3/2*ln(z)^2*ln(1 + z)
       - 3/2*ln(z)^2*ln(1 + z)*z^-1
       - 3/4*ln(z)^2*ln(1 + z)*z
       + 5/12*ln(z)^3
       + 5/3*ln(z)^3*z^-1
       + 35/24*ln(z)^3*z
       - 1/2*ln(z)^2*ln([1-z])
       + 1/2*ln(z)^2*ln([1-z])*z^-1
       + 1/4*ln(z)^2*ln([1-z])*z
       - 3/2*ln(z)*ln([1-z])*z
       - 2*ln(z)*ln([1-z])*ln(1 + z)
       - 2*ln(z)*ln([1-z])*ln(1 + z)*z^-1
       - ln(z)*ln([1-z])*ln(1 + z)*z
       - 3*ln(z)*ln([1-z])^2
       - 3/2*ln(z)*ln([1-z])^2*z^-1
       - 3*ln(z)*ln([1-z])^2*z
       - ln(z)*Li2( - z)
       - ln(z)*Li2( - z)*z^-1
       - 1/2*ln(z)*Li2( - z)*z
       - 5*ln(z)*Li2(z)
       + 5*ln(z)*Li2(z)*z^-1
       + 5/2*ln(z)*Li2(z)*z
       - 73/12*ln([1-z])
       + 34/9*ln([1-z])*z^-1
       + ln([1-z])*z^-1*rln2^2
       + ln([1-z])*rln2^2
       + 23/6*ln([1-z])*z
       + 1/2*ln([1-z])*z*rln2^2
       + 13/18*ln([1-z])*z^2
       + 2/3*ln([1-z])*pi^2*z^-1
       + 5/6*ln([1-z])*pi^2*z
       - 2*ln([1-z])*ln(1 + z)*z^-1*rln2
       - 2*ln([1-z])*ln(1 + z)*rln2
       - ln([1-z])*ln(1 + z)*z*rln2
       + ln([1-z])*ln(1 + z)^2
       + ln([1-z])*ln(1 + z)^2*z^-1
       + 1/2*ln([1-z])*ln(1 + z)^2*z
       - 4*ln([1-z])^2
       + 49/12*ln([1-z])^2*z^-1
       + 1/8*ln([1-z])^2*z
       - 1/3*ln([1-z])^2*z^2
       + ln([1-z])^3
       - ln([1-z])^3*z^-1
       - 1/2*ln([1-z])^3*z
       - 3*ln([1-z])*Li2(1 - z)
       - 3/2*ln([1-z])*Li2(1 - z)*z
       - 2*ln([1-z])*Li2( - z)
       - 2*ln([1-z])*Li2( - z)*z^-1
       - ln([1-z])*Li2( - z)*z
       - 6*ln([1-z])*Li2(z)
       - 3*ln([1-z])*Li2(z)*z
       - 1/3*ln([1+z])*pi^2*z^-1
       - 1/3*ln([1+z])*pi^2
       - 1/6*ln([1+z])*pi^2*z
       + 2*Li3(1/2 - 1/2*z)
       + 2*Li3(1/2 - 1/2*z)*z^-1
       + Li3(1/2 - 1/2*z)*z
       + 2*Li3(1/2 + 1/2*z)
       + 2*Li3(1/2 + 1/2*z)*z^-1
       + Li3(1/2 + 1/2*z)*z
       - 5*Li3(1 - z)
       - 2*Li3(1 - z)*z^-1
       - 5/2*Li3(1 - z)*z
       - Li3( - z)
       - Li3( - z)*z^-1
       - 1/2*Li3( - z)*z
       - 2*Li3(1/(1 + z))
       - 2*Li3(1/(1 + z))*z^-1
       - Li3(1/(1 + z))*z
       + 2*Li3(1/(1 + z) - 1/(1 + z)*z)
       + 2*Li3(1/(1 + z) - 1/(1 + z)*z)*z^-1
       + Li3(1/(1 + z) - 1/(1 + z)*z)*z
       + 5*Li3(z)
       - 13*Li3(z)*z^-1
       - 19/2*Li3(z)*z
       - 1/2*Li2( - z)*z
       + Li2(z)
       - 49/6*Li2(z)*z^-1
       + 2/3*Li2(z)*z^2
       )

       + Dn(0,[1-x])*NC^-2 * ( 17/8
       - 9/4*z^-1
       - 5/4*z
       - 1/3*pi^2*z^-1
       + 5/12*pi^2
       - 5/24*pi^2*z
       - 3/4*ln(z)*z^-1
       + 1/2*ln(z)^2
       - 1/4*ln(z)^2*z
       + ln(z)*ln([1-z])
       - ln(z)*ln([1-z])*z^-1
       - 1/2*ln(z)*ln([1-z])*z
       + 2*ln([1-z])
       - 3/2*ln([1-z])*z^-1
       - 3/8*ln([1-z])*z
       - ln([1-z])^2
       + ln([1-z])^2*z^-1
       + 1/2*ln([1-z])^2*z
       + 5/2*Li2(z)
       - 3*Li2(z)*z^-1
       - 5/4*Li2(z)*z
       )

       + Dn(0,[1-x])*NC^2 * ( 95/24
       - 55/36*z^-1
       - 31/12*z
       - 13/18*z^2
       - 1/3*pi^2*z^-1
       + 1/12*pi^2
       - 3/8*pi^2*z
       + ln(z)
       - 89/12*ln(z)*z^-1
       + ln(z)*z
       + 2/3*ln(z)*z^2
       + 2*ln(z)*ln(1 + z)
       + 2*ln(z)*ln(1 + z)*z^-1
       + ln(z)*ln(1 + z)*z
       - 3/2*ln(z)^2
       - 3*ln(z)^2*z^-1
       - 11/4*ln(z)^2*z
       - 3*ln(z)*ln([1-z])
       + 3*ln(z)*ln([1-z])*z^-1
       + 3/2*ln(z)*ln([1-z])*z
       + 6*ln([1-z])
       - 20/3*ln([1-z])*z^-1
       + 1/8*ln([1-z])*z
       + 2/3*ln([1-z])*z^2
       - 2*ln([1-z])^2
       + 2*ln([1-z])^2*z^-1
       + ln([1-z])^2*z
       + 2*Li2( - z)
       + 2*Li2( - z)*z^-1
       + Li2( - z)*z
       + 1/2*Li2(z)
       + 3*Li2(z)*z^-1
       + 11/4*Li2(z)*z
       )

       + Dn(0,[1-x])*LMUR*NC^-1*NF * (  - 1/3
       + 1/3*z^-1
       + 1/6*z
       )

       + Dn(0,[1-x])*LMUR * ( 11/6
       - 11/6*z^-1
       - 11/12*z
       )

       + Dn(0,[1-x])*LMUR*NC*NF * ( 1/3
       - 1/3*z^-1
       - 1/6*z
       )

       + Dn(0,[1-x])*LMUR*NC^2 * (  - 11/6
       + 11/6*z^-1
       + 11/12*z
       )

       + Dn(0,[1-x])*LMUF*NC^-2 * ( 3/4
       - 3/4*z^-1
       - 7/8*z
       + ln(z)
       - ln(z)*z^-1
       - 1/2*ln(z)*z
       + ln([1-z])
       - ln([1-z])*z^-1
       - 1/2*ln([1-z])*z
       )

       + Dn(0,[1-x])*LMUF * (  - 3/2
       + 3/2*z^-1
       + 7/4*z
       - 2*ln(z)
       + 2*ln(z)*z^-1
       + ln(z)*z
       - 2*ln([1-z])
       + 2*ln([1-z])*z^-1
       + ln([1-z])*z
       )

       + Dn(0,[1-x])*LMUF*NC^2 * ( 3/4
       - 3/4*z^-1
       - 7/8*z
       + ln(z)
       - ln(z)*z^-1
       - 1/2*ln(z)*z
       + ln([1-z])
       - ln([1-z])*z^-1
       - 1/2*ln([1-z])*z
       )

       + Dn(0,[1-x])*LMUA*NC^-2 * (  - 5/4
       + 3/4*z^-1
       + 1/2*z
       - 1/2*ln(z)
       + 1/4*ln(z)*z
       + ln([1-z])
       - ln([1-z])*z^-1
       - 1/2*ln([1-z])*z
       )

       + Dn(0,[1-x])*LMUA*NC^-1*NF * ( 1/3
       - 1/3*z^-1
       - 1/6*z
       )

       + Dn(0,[1-x])*LMUA * ( 14/3
       - 29/6*z^-1
       + 5/12*z
       + 2/3*z^2
       - ln(z)
       - 2*ln(z)*z^-1
       - 5/2*ln(z)*z
       - 4*ln([1-z])
       + 4*ln([1-z])*z^-1
       + 2*ln([1-z])*z
       )

       + Dn(0,[1-x])*LMUA*NC*NF * (  - 1/3
       + 1/3*z^-1
       + 1/6*z
       )

       + Dn(0,[1-x])*LMUA*NC^2 * (  - 41/12
       + 49/12*z^-1
       - 11/12*z
       - 2/3*z^2
       + 3/2*ln(z)
       + 2*ln(z)*z^-1
       + 9/4*ln(z)*z
       + 3*ln([1-z])
       - 3*ln([1-z])*z^-1
       - 3/2*ln([1-z])*z
       )

       + Dn(0,[1-x])*LMUA*LMUF*NC^-2 * (  - 1
       + z^-1
       + 1/2*z
       )

       + Dn(0,[1-x])*LMUA*LMUF * ( 2
       - 2*z^-1
       - z
       )

       + Dn(0,[1-x])*LMUA*LMUF*NC^2 * (  - 1
       + z^-1
       + 1/2*z
       )

       + Dn(0,[1-x]) * (  - 73/12
       + 34/9*z^-1
       + 23/6*z
       + 13/18*z^2
       + 2/3*pi^2*z^-1
       - 1/2*pi^2
       + 7/12*pi^2*z
       - ln(z)
       + 49/6*ln(z)*z^-1
       - ln(z)*z
       - 2/3*ln(z)*z^2
       - 2*ln(z)*ln(1 + z)
       - 2*ln(z)*ln(1 + z)*z^-1
       - ln(z)*ln(1 + z)*z
       + ln(z)^2
       + 3*ln(z)^2*z^-1
       + 3*ln(z)^2*z
       + 2*ln(z)*ln([1-z])
       - 2*ln(z)*ln([1-z])*z^-1
       - ln(z)*ln([1-z])*z
       - 8*ln([1-z])
       + 49/6*ln([1-z])*z^-1
       + 1/4*ln([1-z])*z
       - 2/3*ln([1-z])*z^2
       + 3*ln([1-z])^2
       - 3*ln([1-z])^2*z^-1
       - 3/2*ln([1-z])^2*z
       - 2*Li2( - z)
       - 2*Li2( - z)*z^-1
       - Li2( - z)*z
       - 3*Li2(z)
       - 3/2*Li2(z)*z
       )

       + Dn(1,[1-x])*NC^-2 * ( 2
       - 3/2*z^-1
       - 3/8*z
       - 1/2*ln(z)
       + ln(z)*z^-1
       + 1/4*ln(z)*z
       - 2*ln([1-z])
       + 2*ln([1-z])*z^-1
       + ln([1-z])*z
       )

       + Dn(1,[1-x])*NC^2 * ( 6
       - 20/3*z^-1
       + 1/8*z
       + 2/3*z^2
       - 5/2*ln(z)
       - ln(z)*z^-1
       - 7/4*ln(z)*z
       - 4*ln([1-z])
       + 4*ln([1-z])*z^-1
       + 2*ln([1-z])*z
       )

       + Dn(1,[1-x])*LMUF*NC^-2 * ( 2
       - 2*z^-1
       - z
       )

       + Dn(1,[1-x])*LMUF * (  - 4
       + 4*z^-1
       + 2*z
       )

       + Dn(1,[1-x])*LMUF*NC^2 * ( 2
       - 2*z^-1
       - z
       )

       + Dn(1,[1-x])*LMUA*NC^-2 * ( 1
       - z^-1
       - 1/2*z
       )

       + Dn(1,[1-x])*LMUA * (  - 2
       + 2*z^-1
       + z
       )

       + Dn(1,[1-x])*LMUA*NC^2 * ( 1
       - z^-1
       - 1/2*z
       )

       + Dn(1,[1-x]) * (  - 8
       + 49/6*z^-1
       + 1/4*z
       - 2/3*z^2
       + 3*ln(z)
       + 3/2*ln(z)*z
       + 6*ln([1-z])
       - 6*ln([1-z])*z^-1
       - 3*ln([1-z])*z
       )

       + Dn(2,[1-x])*NC^-2 * (  - 3/2
       + 3/2*z^-1
       + 3/4*z
       )

       + Dn(2,[1-x])*NC^2 * (  - 3/2
       + 3/2*z^-1
       + 3/4*z
       )

       + Dn(2,[1-x]) * ( 3
       - 3*z^-1
       - 3/2*z
       )

       + T(u1)*NC^-2 * (  - 10
       - 4*[1-x]^-1*z^-1
       + 4*[1-x]^-1
       - 6*[1-x]^-1*z
       + 9*z^-1
       + 7*z
       - 5*x*z^-1
       + 5*x
       + 5/6*pi^2*[1-x]^-1*z^-1
       - pi^2*[1-x]^-1
       + 1/2*pi^2*[1-x]^-1*z
       - 3/4*pi^2*z^-1
       + 3/2*pi^2
       - 5/6*pi^2*z
       - 1/12*pi^2*x*z^-1
       + 1/6*pi^2*x
       + 27*ln(x)
       + 18*ln(x)*[1-x]^-1*z^-1
       - 24*ln(x)*[1-x]^-1
       + 15*ln(x)*[1-x]^-1*z
       - 18*ln(x)*z^-1
       - 18*ln(x)*z
       + 6*ln(x)*x
       - 10*ln(x)^2
       - 13/2*ln(x)^2*[1-x]^-1*z^-1
       + 19/2*ln(x)^2*[1-x]^-1
       - 19/4*ln(x)^2*[1-x]^-1*z
       + 5*ln(x)^2*z^-1
       + 13/2*ln(x)^2*z
       + 3/2*ln(x)^2*x*z^-1
       - 3*ln(x)^2*x
       + 16*ln(x)*ln([1-x])
       + 10*ln(x)*ln([1-x])*[1-x]^-1*z^-1
       - 14*ln(x)*ln([1-x])*[1-x]^-1
       + 7*ln(x)*ln([1-x])*[1-x]^-1*z
       - 8*ln(x)*ln([1-x])*z^-1
       - 10*ln(x)*ln([1-x])*z
       - 2*ln(x)*ln([1-x])*x*z^-1
       + 4*ln(x)*ln([1-x])*x
       + 11*ln(x)*ln(z)
       + 8*ln(x)*ln(z)*[1-x]^-1*z^-1
       - 13*ln(x)*ln(z)*[1-x]^-1
       + 13/2*ln(x)*ln(z)*[1-x]^-1*z
       - 11/2*ln(x)*ln(z)*z^-1
       - 8*ln(x)*ln(z)*z
       - 5/2*ln(x)*ln(z)*x*z^-1
       + 5*ln(x)*ln(z)*x
       + 17*ln(x)*ln([1-z])
       + 11*ln(x)*ln([1-z])*[1-x]^-1*z^-1
       - 16*ln(x)*ln([1-z])*[1-x]^-1
       + 8*ln(x)*ln([1-z])*[1-x]^-1*z
       - 17/2*ln(x)*ln([1-z])*z^-1
       - 11*ln(x)*ln([1-z])*z
       - 5/2*ln(x)*ln([1-z])*x*z^-1
       + 5*ln(x)*ln([1-z])*x
       - 3*ln(x)*ln([1-x-z])
       - 2*ln(x)*ln([1-x-z])*[1-x]^-1*z^-1
       + 3*ln(x)*ln([1-x-z])*[1-x]^-1
       - 3/2*ln(x)*ln([1-x-z])*[1-x]^-1*z
       + 3/2*ln(x)*ln([1-x-z])*z^-1
       + 2*ln(x)*ln([1-x-z])*z
       + 1/2*ln(x)*ln([1-x-z])*x*z^-1
       - ln(x)*ln([1-x-z])*x
       - ln(x)*ln([x-z])
       - ln(x)*ln([x-z])*[1-x]^-1*z^-1
       + 2*ln(x)*ln([x-z])*[1-x]^-1
       - ln(x)*ln([x-z])*[1-x]^-1*z
       + 1/2*ln(x)*ln([x-z])*z^-1
       + ln(x)*ln([x-z])*z
       + 1/2*ln(x)*ln([x-z])*x*z^-1
       - ln(x)*ln([x-z])*x
       - 13*ln([1-x])
       - 9*ln([1-x])*[1-x]^-1*z^-1
       + 11*ln([1-x])*[1-x]^-1
       - 7*ln([1-x])*[1-x]^-1*z
       + 9*ln([1-x])*z^-1
       + 9*ln([1-x])*z
       - 3*ln([1-x])*x
       - 9/2*ln([1-x])^2
       - 5/2*ln([1-x])^2*[1-x]^-1*z^-1
       + 3*ln([1-x])^2*[1-x]^-1
       - 3/2*ln([1-x])^2*[1-x]^-1*z
       + 9/4*ln([1-x])^2*z^-1
       + 5/2*ln([1-x])^2*z
       + 1/4*ln([1-x])^2*x*z^-1
       - 1/2*ln([1-x])^2*x
       - 6*ln([1-x])*ln(z)
       - 4*ln([1-x])*ln(z)*[1-x]^-1*z^-1
       + 6*ln([1-x])*ln(z)*[1-x]^-1
       - 3*ln([1-x])*ln(z)*[1-x]^-1*z
       + 3*ln([1-x])*ln(z)*z^-1
       + 4*ln([1-x])*ln(z)*z
       + ln([1-x])*ln(z)*x*z^-1
       - 2*ln([1-x])*ln(z)*x
       - 13*ln([1-x])*ln([1-z])
       - 8*ln([1-x])*ln([1-z])*[1-x]^-1*z^-1
       + 11*ln([1-x])*ln([1-z])*[1-x]^-1
       - 11/2*ln([1-x])*ln([1-z])*[1-x]^-1*z
       + 13/2*ln([1-x])*ln([1-z])*z^-1
       + 8*ln([1-x])*ln([1-z])*z
       + 3/2*ln([1-x])*ln([1-z])*x*z^-1
       - 3*ln([1-x])*ln([1-z])*x
       - 14*ln(z)
       - 9*ln(z)*[1-x]^-1*z^-1
       + 13*ln(z)*[1-x]^-1
       - 8*ln(z)*[1-x]^-1*z
       + 9*ln(z)*z^-1
       + 9*ln(z)*z
       - 3*ln(z)*x
       - 3*ln(z)^2
       - 5/2*ln(z)^2*[1-x]^-1*z^-1
       + 9/2*ln(z)^2*[1-x]^-1
       - 9/4*ln(z)^2*[1-x]^-1*z
       + 3/2*ln(z)^2*z^-1
       + 5/2*ln(z)^2*z
       + ln(z)^2*x*z^-1
       - 2*ln(z)^2*x
       - 8*ln(z)*ln([1-z])
       - 6*ln(z)*ln([1-z])*[1-x]^-1*z^-1
       + 10*ln(z)*ln([1-z])*[1-x]^-1
       - 5*ln(z)*ln([1-z])*[1-x]^-1*z
       + 4*ln(z)*ln([1-z])*z^-1
       + 6*ln(z)*ln([1-z])*z
       + 2*ln(z)*ln([1-z])*x*z^-1
       - 4*ln(z)*ln([1-z])*x
       + ln(z)*ln([x-z])
       + ln(z)*ln([x-z])*[1-x]^-1*z^-1
       - 2*ln(z)*ln([x-z])*[1-x]^-1
       + ln(z)*ln([x-z])*[1-x]^-1*z
       - 1/2*ln(z)*ln([x-z])*z^-1
       - ln(z)*ln([x-z])*z
       - 1/2*ln(z)*ln([x-z])*x*z^-1
       + ln(z)*ln([x-z])*x
       - 18*ln([1-z])
       - 12*ln([1-z])*[1-x]^-1*z^-1
       + 16*ln([1-z])*[1-x]^-1
       - 10*ln([1-z])*[1-x]^-1*z
       + 12*ln([1-z])*z^-1
       + 12*ln([1-z])*z
       - 4*ln([1-z])*x
       - 13/2*ln([1-z])^2
       - 4*ln([1-z])^2*[1-x]^-1*z^-1
       + 11/2*ln([1-z])^2*[1-x]^-1
       - 11/4*ln([1-z])^2*[1-x]^-1*z
       + 13/4*ln([1-z])^2*z^-1
       + 4*ln([1-z])^2*z
       + 3/4*ln([1-z])^2*x*z^-1
       - 3/2*ln([1-z])^2*x
       + 3*ln([1-z])*ln([1-x-z])
       + 2*ln([1-z])*ln([1-x-z])*[1-x]^-1*z^-1
       - 3*ln([1-z])*ln([1-x-z])*[1-x]^-1
       + 3/2*ln([1-z])*ln([1-x-z])*[1-x]^-1*z
       - 3/2*ln([1-z])*ln([1-x-z])*z^-1
       - 2*ln([1-z])*ln([1-x-z])*z
       - 1/2*ln([1-z])*ln([1-x-z])*x*z^-1
       + ln([1-z])*ln([1-x-z])*x
       + Li2(x^-1*[1-x]^-1*z*[1-z])
       + Li2(x^-1*[1-x]^-1*z*[1-z])*[1-x]^-1*z^-1
       - 2*Li2(x^-1*[1-x]^-1*z*[1-z])*[1-x]^-1
       + Li2(x^-1*[1-x]^-1*z*[1-z])*[1-x]^-1*z
       - 1/2*Li2(x^-1*[1-x]^-1*z*[1-z])*z^-1
       - Li2(x^-1*[1-x]^-1*z*[1-z])*z
       - 1/2*Li2(x^-1*[1-x]^-1*z*[1-z])*x*z^-1
       + Li2(x^-1*[1-x]^-1*z*[1-z])*x
       + Li2([1-x]^-1*z)
       + Li2([1-x]^-1*z)*[1-x]^-1
       - 1/2*Li2([1-x]^-1*z)*[1-x]^-1*z
       - 1/2*Li2([1-x]^-1*z)*z^-1
       + 1/2*Li2([1-x]^-1*z)*x*z^-1
       - Li2([1-x]^-1*z)*x
       + Li2([1-x]*[1-z]^-1)
       + Li2([1-x]*[1-z]^-1)*[1-x]^-1*z^-1
       - 2*Li2([1-x]*[1-z]^-1)*[1-x]^-1
       + Li2([1-x]*[1-z]^-1)*[1-x]^-1*z
       - 1/2*Li2([1-x]*[1-z]^-1)*z^-1
       - Li2([1-x]*[1-z]^-1)*z
       - 1/2*Li2([1-x]*[1-z]^-1)*x*z^-1
       + Li2([1-x]*[1-z]^-1)*x
       - 2*Li2(x*[1-x]^-1*z*[1-z]^-1)
       - Li2(x*[1-x]^-1*z*[1-z]^-1)*[1-x]^-1*z^-1
       + Li2(x*[1-x]^-1*z*[1-z]^-1)*[1-x]^-1
       - 1/2*Li2(x*[1-x]^-1*z*[1-z]^-1)*[1-x]^-1*z
       + Li2(x*[1-x]^-1*z*[1-z]^-1)*z^-1
       + Li2(x*[1-x]^-1*z*[1-z]^-1)*z
       - 2*Li2(z)
       - Li2(z)*[1-x]^-1*z^-1
       + Li2(z)*[1-x]^-1
       - 1/2*Li2(z)*[1-x]^-1*z
       + Li2(z)*z^-1
       + Li2(z)*z
       )

       + T(u1) * ( 10
       + 4*[1-x]^-1*z^-1
       - 4*[1-x]^-1
       + 6*[1-x]^-1*z
       - 9*z^-1
       - 7*z
       + 5*x*z^-1
       - 5*x
       - 5/6*pi^2*[1-x]^-1*z^-1
       + pi^2*[1-x]^-1
       - 1/2*pi^2*[1-x]^-1*z
       + 3/4*pi^2*z^-1
       - 3/2*pi^2
       + 5/6*pi^2*z
       + 1/12*pi^2*x*z^-1
       - 1/6*pi^2*x
       - 27*ln(x)
       - 18*ln(x)*[1-x]^-1*z^-1
       + 24*ln(x)*[1-x]^-1
       - 15*ln(x)*[1-x]^-1*z
       + 18*ln(x)*z^-1
       + 18*ln(x)*z
       - 6*ln(x)*x
       + 10*ln(x)^2
       + 13/2*ln(x)^2*[1-x]^-1*z^-1
       - 19/2*ln(x)^2*[1-x]^-1
       + 19/4*ln(x)^2*[1-x]^-1*z
       - 5*ln(x)^2*z^-1
       - 13/2*ln(x)^2*z
       - 3/2*ln(x)^2*x*z^-1
       + 3*ln(x)^2*x
       - 16*ln(x)*ln([1-x])
       - 10*ln(x)*ln([1-x])*[1-x]^-1*z^-1
       + 14*ln(x)*ln([1-x])*[1-x]^-1
       - 7*ln(x)*ln([1-x])*[1-x]^-1*z
       + 8*ln(x)*ln([1-x])*z^-1
       + 10*ln(x)*ln([1-x])*z
       + 2*ln(x)*ln([1-x])*x*z^-1
       - 4*ln(x)*ln([1-x])*x
       - 11*ln(x)*ln(z)
       - 8*ln(x)*ln(z)*[1-x]^-1*z^-1
       + 13*ln(x)*ln(z)*[1-x]^-1
       - 13/2*ln(x)*ln(z)*[1-x]^-1*z
       + 11/2*ln(x)*ln(z)*z^-1
       + 8*ln(x)*ln(z)*z
       + 5/2*ln(x)*ln(z)*x*z^-1
       - 5*ln(x)*ln(z)*x
       - 17*ln(x)*ln([1-z])
       - 11*ln(x)*ln([1-z])*[1-x]^-1*z^-1
       + 16*ln(x)*ln([1-z])*[1-x]^-1
       - 8*ln(x)*ln([1-z])*[1-x]^-1*z
       + 17/2*ln(x)*ln([1-z])*z^-1
       + 11*ln(x)*ln([1-z])*z
       + 5/2*ln(x)*ln([1-z])*x*z^-1
       - 5*ln(x)*ln([1-z])*x
       + 3*ln(x)*ln([1-x-z])
       + 2*ln(x)*ln([1-x-z])*[1-x]^-1*z^-1
       - 3*ln(x)*ln([1-x-z])*[1-x]^-1
       + 3/2*ln(x)*ln([1-x-z])*[1-x]^-1*z
       - 3/2*ln(x)*ln([1-x-z])*z^-1
       - 2*ln(x)*ln([1-x-z])*z
       - 1/2*ln(x)*ln([1-x-z])*x*z^-1
       + ln(x)*ln([1-x-z])*x
       + ln(x)*ln([x-z])
       + ln(x)*ln([x-z])*[1-x]^-1*z^-1
       - 2*ln(x)*ln([x-z])*[1-x]^-1
       + ln(x)*ln([x-z])*[1-x]^-1*z
       - 1/2*ln(x)*ln([x-z])*z^-1
       - ln(x)*ln([x-z])*z
       - 1/2*ln(x)*ln([x-z])*x*z^-1
       + ln(x)*ln([x-z])*x
       + 13*ln([1-x])
       + 9*ln([1-x])*[1-x]^-1*z^-1
       - 11*ln([1-x])*[1-x]^-1
       + 7*ln([1-x])*[1-x]^-1*z
       - 9*ln([1-x])*z^-1
       - 9*ln([1-x])*z
       + 3*ln([1-x])*x
       + 9/2*ln([1-x])^2
       + 5/2*ln([1-x])^2*[1-x]^-1*z^-1
       - 3*ln([1-x])^2*[1-x]^-1
       + 3/2*ln([1-x])^2*[1-x]^-1*z
       - 9/4*ln([1-x])^2*z^-1
       - 5/2*ln([1-x])^2*z
       - 1/4*ln([1-x])^2*x*z^-1
       + 1/2*ln([1-x])^2*x
       + 6*ln([1-x])*ln(z)
       + 4*ln([1-x])*ln(z)*[1-x]^-1*z^-1
       - 6*ln([1-x])*ln(z)*[1-x]^-1
       + 3*ln([1-x])*ln(z)*[1-x]^-1*z
       - 3*ln([1-x])*ln(z)*z^-1
       - 4*ln([1-x])*ln(z)*z
       - ln([1-x])*ln(z)*x*z^-1
       + 2*ln([1-x])*ln(z)*x
       + 13*ln([1-x])*ln([1-z])
       + 8*ln([1-x])*ln([1-z])*[1-x]^-1*z^-1
       - 11*ln([1-x])*ln([1-z])*[1-x]^-1
       + 11/2*ln([1-x])*ln([1-z])*[1-x]^-1*z
       - 13/2*ln([1-x])*ln([1-z])*z^-1
       - 8*ln([1-x])*ln([1-z])*z
       - 3/2*ln([1-x])*ln([1-z])*x*z^-1
       + 3*ln([1-x])*ln([1-z])*x
       + 14*ln(z)
       + 9*ln(z)*[1-x]^-1*z^-1
       - 13*ln(z)*[1-x]^-1
       + 8*ln(z)*[1-x]^-1*z
       - 9*ln(z)*z^-1
       - 9*ln(z)*z
       + 3*ln(z)*x
       + 3*ln(z)^2
       + 5/2*ln(z)^2*[1-x]^-1*z^-1
       - 9/2*ln(z)^2*[1-x]^-1
       + 9/4*ln(z)^2*[1-x]^-1*z
       - 3/2*ln(z)^2*z^-1
       - 5/2*ln(z)^2*z
       - ln(z)^2*x*z^-1
       + 2*ln(z)^2*x
       + 8*ln(z)*ln([1-z])
       + 6*ln(z)*ln([1-z])*[1-x]^-1*z^-1
       - 10*ln(z)*ln([1-z])*[1-x]^-1
       + 5*ln(z)*ln([1-z])*[1-x]^-1*z
       - 4*ln(z)*ln([1-z])*z^-1
       - 6*ln(z)*ln([1-z])*z
       - 2*ln(z)*ln([1-z])*x*z^-1
       + 4*ln(z)*ln([1-z])*x
       - ln(z)*ln([x-z])
       - ln(z)*ln([x-z])*[1-x]^-1*z^-1
       + 2*ln(z)*ln([x-z])*[1-x]^-1
       - ln(z)*ln([x-z])*[1-x]^-1*z
       + 1/2*ln(z)*ln([x-z])*z^-1
       + ln(z)*ln([x-z])*z
       + 1/2*ln(z)*ln([x-z])*x*z^-1
       - ln(z)*ln([x-z])*x
       + 18*ln([1-z])
       + 12*ln([1-z])*[1-x]^-1*z^-1
       - 16*ln([1-z])*[1-x]^-1
       + 10*ln([1-z])*[1-x]^-1*z
       - 12*ln([1-z])*z^-1
       - 12*ln([1-z])*z
       + 4*ln([1-z])*x
       + 13/2*ln([1-z])^2
       + 4*ln([1-z])^2*[1-x]^-1*z^-1
       - 11/2*ln([1-z])^2*[1-x]^-1
       + 11/4*ln([1-z])^2*[1-x]^-1*z
       - 13/4*ln([1-z])^2*z^-1
       - 4*ln([1-z])^2*z
       - 3/4*ln([1-z])^2*x*z^-1
       + 3/2*ln([1-z])^2*x
       - 3*ln([1-z])*ln([1-x-z])
       - 2*ln([1-z])*ln([1-x-z])*[1-x]^-1*z^-1
       + 3*ln([1-z])*ln([1-x-z])*[1-x]^-1
       - 3/2*ln([1-z])*ln([1-x-z])*[1-x]^-1*z
       + 3/2*ln([1-z])*ln([1-x-z])*z^-1
       + 2*ln([1-z])*ln([1-x-z])*z
       + 1/2*ln([1-z])*ln([1-x-z])*x*z^-1
       - ln([1-z])*ln([1-x-z])*x
       - Li2(x^-1*[1-x]^-1*z*[1-z])
       - Li2(x^-1*[1-x]^-1*z*[1-z])*[1-x]^-1*z^-1
       + 2*Li2(x^-1*[1-x]^-1*z*[1-z])*[1-x]^-1
       - Li2(x^-1*[1-x]^-1*z*[1-z])*[1-x]^-1*z
       + 1/2*Li2(x^-1*[1-x]^-1*z*[1-z])*z^-1
       + Li2(x^-1*[1-x]^-1*z*[1-z])*z
       + 1/2*Li2(x^-1*[1-x]^-1*z*[1-z])*x*z^-1
       - Li2(x^-1*[1-x]^-1*z*[1-z])*x
       - Li2([1-x]^-1*z)
       - Li2([1-x]^-1*z)*[1-x]^-1
       + 1/2*Li2([1-x]^-1*z)*[1-x]^-1*z
       + 1/2*Li2([1-x]^-1*z)*z^-1
       - 1/2*Li2([1-x]^-1*z)*x*z^-1
       + Li2([1-x]^-1*z)*x
       - Li2([1-x]*[1-z]^-1)
       - Li2([1-x]*[1-z]^-1)*[1-x]^-1*z^-1
       + 2*Li2([1-x]*[1-z]^-1)*[1-x]^-1
       - Li2([1-x]*[1-z]^-1)*[1-x]^-1*z
       + 1/2*Li2([1-x]*[1-z]^-1)*z^-1
       + Li2([1-x]*[1-z]^-1)*z
       + 1/2*Li2([1-x]*[1-z]^-1)*x*z^-1
       - Li2([1-x]*[1-z]^-1)*x
       + 2*Li2(x*[1-x]^-1*z*[1-z]^-1)
       + Li2(x*[1-x]^-1*z*[1-z]^-1)*[1-x]^-1*z^-1
       - Li2(x*[1-x]^-1*z*[1-z]^-1)*[1-x]^-1
       + 1/2*Li2(x*[1-x]^-1*z*[1-z]^-1)*[1-x]^-1*z
       - Li2(x*[1-x]^-1*z*[1-z]^-1)*z^-1
       - Li2(x*[1-x]^-1*z*[1-z]^-1)*z
       + 2*Li2(z)
       + Li2(z)*[1-x]^-1*z^-1
       - Li2(z)*[1-x]^-1
       + 1/2*Li2(z)*[1-x]^-1*z
       - Li2(z)*z^-1
       - Li2(z)*z
       )

       + T(u2)*NC^-2 * (  - 10
       - 4*[1-x]^-1*z^-1
       + 4*[1-x]^-1
       - 6*[1-x]^-1*z
       + 9*z^-1
       + 7*z
       - 5*x*z^-1
       + 5*x
       + 5/6*pi^2*[1-x]^-1*z^-1
       - pi^2*[1-x]^-1
       + 1/2*pi^2*[1-x]^-1*z
       - 3/4*pi^2*z^-1
       + 3/2*pi^2
       - 5/6*pi^2*z
       - 1/12*pi^2*x*z^-1
       + 1/6*pi^2*x
       + 27*ln(x)
       + 18*ln(x)*[1-x]^-1*z^-1
       - 24*ln(x)*[1-x]^-1
       + 15*ln(x)*[1-x]^-1*z
       - 18*ln(x)*z^-1
       - 18*ln(x)*z
       + 6*ln(x)*x
       - 3*ln(x)*ln( - [1-x-z])
       - 2*ln(x)*ln( - [1-x-z])*[1-x]^-1*z^-1
       + 3*ln(x)*ln( - [1-x-z])*[1-x]^-1
       - 3/2*ln(x)*ln( - [1-x-z])*[1-x]^-1*z
       + 3/2*ln(x)*ln( - [1-x-z])*z^-1
       + 2*ln(x)*ln( - [1-x-z])*z
       + 1/2*ln(x)*ln( - [1-x-z])*x*z^-1
       - ln(x)*ln( - [1-x-z])*x
       - 19/2*ln(x)^2
       - 13/2*ln(x)^2*[1-x]^-1*z^-1
       + 10*ln(x)^2*[1-x]^-1
       - 5*ln(x)^2*[1-x]^-1*z
       + 19/4*ln(x)^2*z^-1
       + 13/2*ln(x)^2*z
       + 7/4*ln(x)^2*x*z^-1
       - 7/2*ln(x)^2*x
       + 13*ln(x)*ln([1-x])
       + 8*ln(x)*ln([1-x])*[1-x]^-1*z^-1
       - 11*ln(x)*ln([1-x])*[1-x]^-1
       + 11/2*ln(x)*ln([1-x])*[1-x]^-1*z
       - 13/2*ln(x)*ln([1-x])*z^-1
       - 8*ln(x)*ln([1-x])*z
       - 3/2*ln(x)*ln([1-x])*x*z^-1
       + 3*ln(x)*ln([1-x])*x
       + 14*ln(x)*ln(z)
       + 10*ln(x)*ln(z)*[1-x]^-1*z^-1
       - 16*ln(x)*ln(z)*[1-x]^-1
       + 8*ln(x)*ln(z)*[1-x]^-1*z
       - 7*ln(x)*ln(z)*z^-1
       - 10*ln(x)*ln(z)*z
       - 3*ln(x)*ln(z)*x*z^-1
       + 6*ln(x)*ln(z)*x
       + 16*ln(x)*ln([1-z])
       + 11*ln(x)*ln([1-z])*[1-x]^-1*z^-1
       - 17*ln(x)*ln([1-z])*[1-x]^-1
       + 17/2*ln(x)*ln([1-z])*[1-x]^-1*z
       - 8*ln(x)*ln([1-z])*z^-1
       - 11*ln(x)*ln([1-z])*z
       - 3*ln(x)*ln([1-z])*x*z^-1
       + 6*ln(x)*ln([1-z])*x
       - ln(x)*ln([x-z])
       - ln(x)*ln([x-z])*[1-x]^-1*z^-1
       + 2*ln(x)*ln([x-z])*[1-x]^-1
       - ln(x)*ln([x-z])*[1-x]^-1*z
       + 1/2*ln(x)*ln([x-z])*z^-1
       + ln(x)*ln([x-z])*z
       + 1/2*ln(x)*ln([x-z])*x*z^-1
       - ln(x)*ln([x-z])*x
       - 13*ln([1-x])
       - 9*ln([1-x])*[1-x]^-1*z^-1
       + 11*ln([1-x])*[1-x]^-1
       - 7*ln([1-x])*[1-x]^-1*z
       + 9*ln([1-x])*z^-1
       + 9*ln([1-x])*z
       - 3*ln([1-x])*x
       - 9/2*ln([1-x])^2
       - 5/2*ln([1-x])^2*[1-x]^-1*z^-1
       + 3*ln([1-x])^2*[1-x]^-1
       - 3/2*ln([1-x])^2*[1-x]^-1*z
       + 9/4*ln([1-x])^2*z^-1
       + 5/2*ln([1-x])^2*z
       + 1/4*ln([1-x])^2*x*z^-1
       - 1/2*ln([1-x])^2*x
       - 6*ln([1-x])*ln(z)
       - 4*ln([1-x])*ln(z)*[1-x]^-1*z^-1
       + 6*ln([1-x])*ln(z)*[1-x]^-1
       - 3*ln([1-x])*ln(z)*[1-x]^-1*z
       + 3*ln([1-x])*ln(z)*z^-1
       + 4*ln([1-x])*ln(z)*z
       + ln([1-x])*ln(z)*x*z^-1
       - 2*ln([1-x])*ln(z)*x
       - 10*ln([1-x])*ln([1-z])
       - 6*ln([1-x])*ln([1-z])*[1-x]^-1*z^-1
       + 8*ln([1-x])*ln([1-z])*[1-x]^-1
       - 4*ln([1-x])*ln([1-z])*[1-x]^-1*z
       + 5*ln([1-x])*ln([1-z])*z^-1
       + 6*ln([1-x])*ln([1-z])*z
       + ln([1-x])*ln([1-z])*x*z^-1
       - 2*ln([1-x])*ln([1-z])*x
       - 14*ln(z)
       - 9*ln(z)*[1-x]^-1*z^-1
       + 13*ln(z)*[1-x]^-1
       - 8*ln(z)*[1-x]^-1*z
       + 9*ln(z)*z^-1
       + 9*ln(z)*z
       - 3*ln(z)*x
       - 3*ln(z)^2
       - 5/2*ln(z)^2*[1-x]^-1*z^-1
       + 9/2*ln(z)^2*[1-x]^-1
       - 9/4*ln(z)^2*[1-x]^-1*z
       + 3/2*ln(z)^2*z^-1
       + 5/2*ln(z)^2*z
       + ln(z)^2*x*z^-1
       - 2*ln(z)^2*x
       - 11*ln(z)*ln([1-z])
       - 8*ln(z)*ln([1-z])*[1-x]^-1*z^-1
       + 13*ln(z)*ln([1-z])*[1-x]^-1
       - 13/2*ln(z)*ln([1-z])*[1-x]^-1*z
       + 11/2*ln(z)*ln([1-z])*z^-1
       + 8*ln(z)*ln([1-z])*z
       + 5/2*ln(z)*ln([1-z])*x*z^-1
       - 5*ln(z)*ln([1-z])*x
       + ln(z)*ln([x-z])
       + ln(z)*ln([x-z])*[1-x]^-1*z^-1
       - 2*ln(z)*ln([x-z])*[1-x]^-1
       + ln(z)*ln([x-z])*[1-x]^-1*z
       - 1/2*ln(z)*ln([x-z])*z^-1
       - ln(z)*ln([x-z])*z
       - 1/2*ln(z)*ln([x-z])*x*z^-1
       + ln(z)*ln([x-z])*x
       - 18*ln([1-z])
       - 12*ln([1-z])*[1-x]^-1*z^-1
       + 16*ln([1-z])*[1-x]^-1
       - 10*ln([1-z])*[1-x]^-1*z
       + 12*ln([1-z])*z^-1
       + 12*ln([1-z])*z
       - 4*ln([1-z])*x
       + 3*ln([1-z])*ln( - [1-x-z])
       + 2*ln([1-z])*ln( - [1-x-z])*[1-x]^-1*z^-1
       - 3*ln([1-z])*ln( - [1-x-z])*[1-x]^-1
       + 3/2*ln([1-z])*ln( - [1-x-z])*[1-x]^-1*z
       - 3/2*ln([1-z])*ln( - [1-x-z])*z^-1
       - 2*ln([1-z])*ln( - [1-x-z])*z
       - 1/2*ln([1-z])*ln( - [1-x-z])*x*z^-1
       + ln([1-z])*ln( - [1-x-z])*x
       - 6*ln([1-z])^2
       - 4*ln([1-z])^2*[1-x]^-1*z^-1
       + 6*ln([1-z])^2*[1-x]^-1
       - 3*ln([1-z])^2*[1-x]^-1*z
       + 3*ln([1-z])^2*z^-1
       + 4*ln([1-z])^2*z
       + ln([1-z])^2*x*z^-1
       - 2*ln([1-z])^2*x
       + 2*Li2(x^-1*[1-x]*z^-1*[1-z])
       + Li2(x^-1*[1-x]*z^-1*[1-z])*[1-x]^-1*z^-1
       - Li2(x^-1*[1-x]*z^-1*[1-z])*[1-x]^-1
       + 1/2*Li2(x^-1*[1-x]*z^-1*[1-z])*[1-x]^-1*z
       - Li2(x^-1*[1-x]*z^-1*[1-z])*z^-1
       - Li2(x^-1*[1-x]*z^-1*[1-z])*z
       - Li2([1-x]*z^-1)
       - Li2([1-x]*z^-1)*[1-x]^-1
       + 1/2*Li2([1-x]*z^-1)*[1-x]^-1*z
       + 1/2*Li2([1-x]*z^-1)*z^-1
       - 1/2*Li2([1-x]*z^-1)*x*z^-1
       + Li2([1-x]*z^-1)*x
       + Li2([1-x]*[1-z]^-1)
       + Li2([1-x]*[1-z]^-1)*[1-x]^-1*z^-1
       - 2*Li2([1-x]*[1-z]^-1)*[1-x]^-1
       + Li2([1-x]*[1-z]^-1)*[1-x]^-1*z
       - 1/2*Li2([1-x]*[1-z]^-1)*z^-1
       - Li2([1-x]*[1-z]^-1)*z
       - 1/2*Li2([1-x]*[1-z]^-1)*x*z^-1
       + Li2([1-x]*[1-z]^-1)*x
       - Li2(x*[1-x]*z^-1*[1-z]^-1)
       - Li2(x*[1-x]*z^-1*[1-z]^-1)*[1-x]^-1*z^-1
       + 2*Li2(x*[1-x]*z^-1*[1-z]^-1)*[1-x]^-1
       - Li2(x*[1-x]*z^-1*[1-z]^-1)*[1-x]^-1*z
       + 1/2*Li2(x*[1-x]*z^-1*[1-z]^-1)*z^-1
       + Li2(x*[1-x]*z^-1*[1-z]^-1)*z
       + 1/2*Li2(x*[1-x]*z^-1*[1-z]^-1)*x*z^-1
       - Li2(x*[1-x]*z^-1*[1-z]^-1)*x
       - 2*Li2(z)
       - Li2(z)*[1-x]^-1*z^-1
       + Li2(z)*[1-x]^-1
       - 1/2*Li2(z)*[1-x]^-1*z
       + Li2(z)*z^-1
       + Li2(z)*z
       )

       + T(u2) * ( 10
       + 4*[1-x]^-1*z^-1
       - 4*[1-x]^-1
       + 6*[1-x]^-1*z
       - 9*z^-1
       - 7*z
       + 5*x*z^-1
       - 5*x
       - 5/6*pi^2*[1-x]^-1*z^-1
       + pi^2*[1-x]^-1
       - 1/2*pi^2*[1-x]^-1*z
       + 3/4*pi^2*z^-1
       - 3/2*pi^2
       + 5/6*pi^2*z
       + 1/12*pi^2*x*z^-1
       - 1/6*pi^2*x
       - 27*ln(x)
       - 18*ln(x)*[1-x]^-1*z^-1
       + 24*ln(x)*[1-x]^-1
       - 15*ln(x)*[1-x]^-1*z
       + 18*ln(x)*z^-1
       + 18*ln(x)*z
       - 6*ln(x)*x
       + 3*ln(x)*ln( - [1-x-z])
       + 2*ln(x)*ln( - [1-x-z])*[1-x]^-1*z^-1
       - 3*ln(x)*ln( - [1-x-z])*[1-x]^-1
       + 3/2*ln(x)*ln( - [1-x-z])*[1-x]^-1*z
       - 3/2*ln(x)*ln( - [1-x-z])*z^-1
       - 2*ln(x)*ln( - [1-x-z])*z
       - 1/2*ln(x)*ln( - [1-x-z])*x*z^-1
       + ln(x)*ln( - [1-x-z])*x
       + 19/2*ln(x)^2
       + 13/2*ln(x)^2*[1-x]^-1*z^-1
       - 10*ln(x)^2*[1-x]^-1
       + 5*ln(x)^2*[1-x]^-1*z
       - 19/4*ln(x)^2*z^-1
       - 13/2*ln(x)^2*z
       - 7/4*ln(x)^2*x*z^-1
       + 7/2*ln(x)^2*x
       - 13*ln(x)*ln([1-x])
       - 8*ln(x)*ln([1-x])*[1-x]^-1*z^-1
       + 11*ln(x)*ln([1-x])*[1-x]^-1
       - 11/2*ln(x)*ln([1-x])*[1-x]^-1*z
       + 13/2*ln(x)*ln([1-x])*z^-1
       + 8*ln(x)*ln([1-x])*z
       + 3/2*ln(x)*ln([1-x])*x*z^-1
       - 3*ln(x)*ln([1-x])*x
       - 14*ln(x)*ln(z)
       - 10*ln(x)*ln(z)*[1-x]^-1*z^-1
       + 16*ln(x)*ln(z)*[1-x]^-1
       - 8*ln(x)*ln(z)*[1-x]^-1*z
       + 7*ln(x)*ln(z)*z^-1
       + 10*ln(x)*ln(z)*z
       + 3*ln(x)*ln(z)*x*z^-1
       - 6*ln(x)*ln(z)*x
       - 16*ln(x)*ln([1-z])
       - 11*ln(x)*ln([1-z])*[1-x]^-1*z^-1
       + 17*ln(x)*ln([1-z])*[1-x]^-1
       - 17/2*ln(x)*ln([1-z])*[1-x]^-1*z
       + 8*ln(x)*ln([1-z])*z^-1
       + 11*ln(x)*ln([1-z])*z
       + 3*ln(x)*ln([1-z])*x*z^-1
       - 6*ln(x)*ln([1-z])*x
       + ln(x)*ln([x-z])
       + ln(x)*ln([x-z])*[1-x]^-1*z^-1
       - 2*ln(x)*ln([x-z])*[1-x]^-1
       + ln(x)*ln([x-z])*[1-x]^-1*z
       - 1/2*ln(x)*ln([x-z])*z^-1
       - ln(x)*ln([x-z])*z
       - 1/2*ln(x)*ln([x-z])*x*z^-1
       + ln(x)*ln([x-z])*x
       + 13*ln([1-x])
       + 9*ln([1-x])*[1-x]^-1*z^-1
       - 11*ln([1-x])*[1-x]^-1
       + 7*ln([1-x])*[1-x]^-1*z
       - 9*ln([1-x])*z^-1
       - 9*ln([1-x])*z
       + 3*ln([1-x])*x
       + 9/2*ln([1-x])^2
       + 5/2*ln([1-x])^2*[1-x]^-1*z^-1
       - 3*ln([1-x])^2*[1-x]^-1
       + 3/2*ln([1-x])^2*[1-x]^-1*z
       - 9/4*ln([1-x])^2*z^-1
       - 5/2*ln([1-x])^2*z
       - 1/4*ln([1-x])^2*x*z^-1
       + 1/2*ln([1-x])^2*x
       + 6*ln([1-x])*ln(z)
       + 4*ln([1-x])*ln(z)*[1-x]^-1*z^-1
       - 6*ln([1-x])*ln(z)*[1-x]^-1
       + 3*ln([1-x])*ln(z)*[1-x]^-1*z
       - 3*ln([1-x])*ln(z)*z^-1
       - 4*ln([1-x])*ln(z)*z
       - ln([1-x])*ln(z)*x*z^-1
       + 2*ln([1-x])*ln(z)*x
       + 10*ln([1-x])*ln([1-z])
       + 6*ln([1-x])*ln([1-z])*[1-x]^-1*z^-1
       - 8*ln([1-x])*ln([1-z])*[1-x]^-1
       + 4*ln([1-x])*ln([1-z])*[1-x]^-1*z
       - 5*ln([1-x])*ln([1-z])*z^-1
       - 6*ln([1-x])*ln([1-z])*z
       - ln([1-x])*ln([1-z])*x*z^-1
       + 2*ln([1-x])*ln([1-z])*x
       + 14*ln(z)
       + 9*ln(z)*[1-x]^-1*z^-1
       - 13*ln(z)*[1-x]^-1
       + 8*ln(z)*[1-x]^-1*z
       - 9*ln(z)*z^-1
       - 9*ln(z)*z
       + 3*ln(z)*x
       + 3*ln(z)^2
       + 5/2*ln(z)^2*[1-x]^-1*z^-1
       - 9/2*ln(z)^2*[1-x]^-1
       + 9/4*ln(z)^2*[1-x]^-1*z
       - 3/2*ln(z)^2*z^-1
       - 5/2*ln(z)^2*z
       - ln(z)^2*x*z^-1
       + 2*ln(z)^2*x
       + 11*ln(z)*ln([1-z])
       + 8*ln(z)*ln([1-z])*[1-x]^-1*z^-1
       - 13*ln(z)*ln([1-z])*[1-x]^-1
       + 13/2*ln(z)*ln([1-z])*[1-x]^-1*z
       - 11/2*ln(z)*ln([1-z])*z^-1
       - 8*ln(z)*ln([1-z])*z
       - 5/2*ln(z)*ln([1-z])*x*z^-1
       + 5*ln(z)*ln([1-z])*x
       - ln(z)*ln([x-z])
       - ln(z)*ln([x-z])*[1-x]^-1*z^-1
       + 2*ln(z)*ln([x-z])*[1-x]^-1
       - ln(z)*ln([x-z])*[1-x]^-1*z
       + 1/2*ln(z)*ln([x-z])*z^-1
       + ln(z)*ln([x-z])*z
       + 1/2*ln(z)*ln([x-z])*x*z^-1
       - ln(z)*ln([x-z])*x
       + 18*ln([1-z])
       + 12*ln([1-z])*[1-x]^-1*z^-1
       - 16*ln([1-z])*[1-x]^-1
       + 10*ln([1-z])*[1-x]^-1*z
       - 12*ln([1-z])*z^-1
       - 12*ln([1-z])*z
       + 4*ln([1-z])*x
       - 3*ln([1-z])*ln( - [1-x-z])
       - 2*ln([1-z])*ln( - [1-x-z])*[1-x]^-1*z^-1
       + 3*ln([1-z])*ln( - [1-x-z])*[1-x]^-1
       - 3/2*ln([1-z])*ln( - [1-x-z])*[1-x]^-1*z
       + 3/2*ln([1-z])*ln( - [1-x-z])*z^-1
       + 2*ln([1-z])*ln( - [1-x-z])*z
       + 1/2*ln([1-z])*ln( - [1-x-z])*x*z^-1
       - ln([1-z])*ln( - [1-x-z])*x
       + 6*ln([1-z])^2
       + 4*ln([1-z])^2*[1-x]^-1*z^-1
       - 6*ln([1-z])^2*[1-x]^-1
       + 3*ln([1-z])^2*[1-x]^-1*z
       - 3*ln([1-z])^2*z^-1
       - 4*ln([1-z])^2*z
       - ln([1-z])^2*x*z^-1
       + 2*ln([1-z])^2*x
       - 2*Li2(x^-1*[1-x]*z^-1*[1-z])
       - Li2(x^-1*[1-x]*z^-1*[1-z])*[1-x]^-1*z^-1
       + Li2(x^-1*[1-x]*z^-1*[1-z])*[1-x]^-1
       - 1/2*Li2(x^-1*[1-x]*z^-1*[1-z])*[1-x]^-1*z
       + Li2(x^-1*[1-x]*z^-1*[1-z])*z^-1
       + Li2(x^-1*[1-x]*z^-1*[1-z])*z
       + Li2([1-x]*z^-1)
       + Li2([1-x]*z^-1)*[1-x]^-1
       - 1/2*Li2([1-x]*z^-1)*[1-x]^-1*z
       - 1/2*Li2([1-x]*z^-1)*z^-1
       + 1/2*Li2([1-x]*z^-1)*x*z^-1
       - Li2([1-x]*z^-1)*x
       - Li2([1-x]*[1-z]^-1)
       - Li2([1-x]*[1-z]^-1)*[1-x]^-1*z^-1
       + 2*Li2([1-x]*[1-z]^-1)*[1-x]^-1
       - Li2([1-x]*[1-z]^-1)*[1-x]^-1*z
       + 1/2*Li2([1-x]*[1-z]^-1)*z^-1
       + Li2([1-x]*[1-z]^-1)*z
       + 1/2*Li2([1-x]*[1-z]^-1)*x*z^-1
       - Li2([1-x]*[1-z]^-1)*x
       + Li2(x*[1-x]*z^-1*[1-z]^-1)
       + Li2(x*[1-x]*z^-1*[1-z]^-1)*[1-x]^-1*z^-1
       - 2*Li2(x*[1-x]*z^-1*[1-z]^-1)*[1-x]^-1
       + Li2(x*[1-x]*z^-1*[1-z]^-1)*[1-x]^-1*z
       - 1/2*Li2(x*[1-x]*z^-1*[1-z]^-1)*z^-1
       - Li2(x*[1-x]*z^-1*[1-z]^-1)*z
       - 1/2*Li2(x*[1-x]*z^-1*[1-z]^-1)*x*z^-1
       + Li2(x*[1-x]*z^-1*[1-z]^-1)*x
       + 2*Li2(z)
       + Li2(z)*[1-x]^-1*z^-1
       - Li2(z)*[1-x]^-1
       + 1/2*Li2(z)*[1-x]^-1*z
       - Li2(z)*z^-1
       - Li2(z)*z
       )

       + T(u3)*NC^-2 * (  - 10
       - 4*[1-x]^-1*z^-1
       + 4*[1-x]^-1
       - 6*[1-x]^-1*z
       + 9*z^-1
       + 7*z
       - 5*x*z^-1
       + 5*x
       + 3/2*pi^2*[1-x]^-1*z^-1
       - 7/3*pi^2*[1-x]^-1
       + 7/6*pi^2*[1-x]^-1*z
       - 13/12*pi^2*z^-1
       + 13/6*pi^2
       - 3/2*pi^2*z
       - 5/12*pi^2*x*z^-1
       + 5/6*pi^2*x
       + 27*ln(x)
       + 18*ln(x)*[1-x]^-1*z^-1
       - 24*ln(x)*[1-x]^-1
       + 15*ln(x)*[1-x]^-1*z
       - 18*ln(x)*z^-1
       - 18*ln(x)*z
       + 6*ln(x)*x
       - ln(x)*ln( - [x-z])
       - ln(x)*ln( - [x-z])*[1-x]^-1*z^-1
       + 2*ln(x)*ln( - [x-z])*[1-x]^-1
       - ln(x)*ln( - [x-z])*[1-x]^-1*z
       + 1/2*ln(x)*ln( - [x-z])*z^-1
       + ln(x)*ln( - [x-z])*z
       + 1/2*ln(x)*ln( - [x-z])*x*z^-1
       - ln(x)*ln( - [x-z])*x
       - 21/2*ln(x)^2
       - 7*ln(x)^2*[1-x]^-1*z^-1
       + 21/2*ln(x)^2*[1-x]^-1
       - 21/4*ln(x)^2*[1-x]^-1*z
       + 21/4*ln(x)^2*z^-1
       + 7*ln(x)^2*z
       + 7/4*ln(x)^2*x*z^-1
       - 7/2*ln(x)^2*x
       + 15*ln(x)*ln([1-x])
       + 9*ln(x)*ln([1-x])*[1-x]^-1*z^-1
       - 12*ln(x)*ln([1-x])*[1-x]^-1
       + 6*ln(x)*ln([1-x])*[1-x]^-1*z
       - 15/2*ln(x)*ln([1-x])*z^-1
       - 9*ln(x)*ln([1-x])*z
       - 3/2*ln(x)*ln([1-x])*x*z^-1
       + 3*ln(x)*ln([1-x])*x
       + 12*ln(x)*ln(z)
       + 9*ln(x)*ln(z)*[1-x]^-1*z^-1
       - 15*ln(x)*ln(z)*[1-x]^-1
       + 15/2*ln(x)*ln(z)*[1-x]^-1*z
       - 6*ln(x)*ln(z)*z^-1
       - 9*ln(x)*ln(z)*z
       - 3*ln(x)*ln(z)*x*z^-1
       + 6*ln(x)*ln(z)*x
       + 18*ln(x)*ln([1-z])
       + 12*ln(x)*ln([1-z])*[1-x]^-1*z^-1
       - 18*ln(x)*ln([1-z])*[1-x]^-1
       + 9*ln(x)*ln([1-z])*[1-x]^-1*z
       - 9*ln(x)*ln([1-z])*z^-1
       - 12*ln(x)*ln([1-z])*z
       - 3*ln(x)*ln([1-z])*x*z^-1
       + 6*ln(x)*ln([1-z])*x
       - 3*ln(x)*ln([1-x-z])
       - 2*ln(x)*ln([1-x-z])*[1-x]^-1*z^-1
       + 3*ln(x)*ln([1-x-z])*[1-x]^-1
       - 3/2*ln(x)*ln([1-x-z])*[1-x]^-1*z
       + 3/2*ln(x)*ln([1-x-z])*z^-1
       + 2*ln(x)*ln([1-x-z])*z
       + 1/2*ln(x)*ln([1-x-z])*x*z^-1
       - ln(x)*ln([1-x-z])*x
       - 13*ln([1-x])
       - 9*ln([1-x])*[1-x]^-1*z^-1
       + 11*ln([1-x])*[1-x]^-1
       - 7*ln([1-x])*[1-x]^-1*z
       + 9*ln([1-x])*z^-1
       + 9*ln([1-x])*z
       - 3*ln([1-x])*x
       - 11/2*ln([1-x])^2
       - 7/2*ln([1-x])^2*[1-x]^-1*z^-1
       + 5*ln([1-x])^2*[1-x]^-1
       - 5/2*ln([1-x])^2*[1-x]^-1*z
       + 11/4*ln([1-x])^2*z^-1
       + 7/2*ln([1-x])^2*z
       + 3/4*ln([1-x])^2*x*z^-1
       - 3/2*ln([1-x])^2*x
       - 5*ln([1-x])*ln(z)
       - 3*ln([1-x])*ln(z)*[1-x]^-1*z^-1
       + 4*ln([1-x])*ln(z)*[1-x]^-1
       - 2*ln([1-x])*ln(z)*[1-x]^-1*z
       + 5/2*ln([1-x])*ln(z)*z^-1
       + 3*ln([1-x])*ln(z)*z
       + 1/2*ln([1-x])*ln(z)*x*z^-1
       - ln([1-x])*ln(z)*x
       - 11*ln([1-x])*ln([1-z])
       - 6*ln([1-x])*ln([1-z])*[1-x]^-1*z^-1
       + 7*ln([1-x])*ln([1-z])*[1-x]^-1
       - 7/2*ln([1-x])*ln([1-z])*[1-x]^-1*z
       + 11/2*ln([1-x])*ln([1-z])*z^-1
       + 6*ln([1-x])*ln([1-z])*z
       + 1/2*ln([1-x])*ln([1-z])*x*z^-1
       - ln([1-x])*ln([1-z])*x
       - 14*ln(z)
       - 9*ln(z)*[1-x]^-1*z^-1
       + 13*ln(z)*[1-x]^-1
       - 8*ln(z)*[1-x]^-1*z
       + 9*ln(z)*z^-1
       + 9*ln(z)*z
       - 3*ln(z)*x
       + ln(z)*ln( - [x-z])
       + ln(z)*ln( - [x-z])*[1-x]^-1*z^-1
       - 2*ln(z)*ln( - [x-z])*[1-x]^-1
       + ln(z)*ln( - [x-z])*[1-x]^-1*z
       - 1/2*ln(z)*ln( - [x-z])*z^-1
       - ln(z)*ln( - [x-z])*z
       - 1/2*ln(z)*ln( - [x-z])*x*z^-1
       + ln(z)*ln( - [x-z])*x
       - 7/2*ln(z)^2
       - 3*ln(z)^2*[1-x]^-1*z^-1
       + 11/2*ln(z)^2*[1-x]^-1
       - 11/4*ln(z)^2*[1-x]^-1*z
       + 7/4*ln(z)^2*z^-1
       + 3*ln(z)^2*z
       + 5/4*ln(z)^2*x*z^-1
       - 5/2*ln(z)^2*x
       - 9*ln(z)*ln([1-z])
       - 7*ln(z)*ln([1-z])*[1-x]^-1*z^-1
       + 12*ln(z)*ln([1-z])*[1-x]^-1
       - 6*ln(z)*ln([1-z])*[1-x]^-1*z
       + 9/2*ln(z)*ln([1-z])*z^-1
       + 7*ln(z)*ln([1-z])*z
       + 5/2*ln(z)*ln([1-z])*x*z^-1
       - 5*ln(z)*ln([1-z])*x
       - 18*ln([1-z])
       - 12*ln([1-z])*[1-x]^-1*z^-1
       + 16*ln([1-z])*[1-x]^-1
       - 10*ln([1-z])*[1-x]^-1*z
       + 12*ln([1-z])*z^-1
       + 12*ln([1-z])*z
       - 4*ln([1-z])*x
       - 15/2*ln([1-z])^2
       - 5*ln([1-z])^2*[1-x]^-1*z^-1
       + 15/2*ln([1-z])^2*[1-x]^-1
       - 15/4*ln([1-z])^2*[1-x]^-1*z
       + 15/4*ln([1-z])^2*z^-1
       + 5*ln([1-z])^2*z
       + 5/4*ln([1-z])^2*x*z^-1
       - 5/2*ln([1-z])^2*x
       + 3*ln([1-z])*ln([1-x-z])
       + 2*ln([1-z])*ln([1-x-z])*[1-x]^-1*z^-1
       - 3*ln([1-z])*ln([1-x-z])*[1-x]^-1
       + 3/2*ln([1-z])*ln([1-x-z])*[1-x]^-1*z
       - 3/2*ln([1-z])*ln([1-x-z])*z^-1
       - 2*ln([1-z])*ln([1-x-z])*z
       - 1/2*ln([1-z])*ln([1-x-z])*x*z^-1
       + ln([1-z])*ln([1-x-z])*x
       - Li2([1-x]^-1*[1-z])
       - Li2([1-x]^-1*[1-z])*[1-x]^-1*z^-1
       + 2*Li2([1-x]^-1*[1-z])*[1-x]^-1
       - Li2([1-x]^-1*[1-z])*[1-x]^-1*z
       + 1/2*Li2([1-x]^-1*[1-z])*z^-1
       + Li2([1-x]^-1*[1-z])*z
       + 1/2*Li2([1-x]^-1*[1-z])*x*z^-1
       - Li2([1-x]^-1*[1-z])*x
       + Li2([1-x]^-1*z)
       + Li2([1-x]^-1*z)*[1-x]^-1
       - 1/2*Li2([1-x]^-1*z)*[1-x]^-1*z
       - 1/2*Li2([1-x]^-1*z)*z^-1
       + 1/2*Li2([1-x]^-1*z)*x*z^-1
       - Li2([1-x]^-1*z)*x
       - 2*Li2(x*[1-x]^-1*z*[1-z]^-1)
       - Li2(x*[1-x]^-1*z*[1-z]^-1)*[1-x]^-1*z^-1
       + Li2(x*[1-x]^-1*z*[1-z]^-1)*[1-x]^-1
       - 1/2*Li2(x*[1-x]^-1*z*[1-z]^-1)*[1-x]^-1*z
       + Li2(x*[1-x]^-1*z*[1-z]^-1)*z^-1
       + Li2(x*[1-x]^-1*z*[1-z]^-1)*z
       - Li2(x*[1-x]*z^-1*[1-z]^-1)
       - Li2(x*[1-x]*z^-1*[1-z]^-1)*[1-x]^-1*z^-1
       + 2*Li2(x*[1-x]*z^-1*[1-z]^-1)*[1-x]^-1
       - Li2(x*[1-x]*z^-1*[1-z]^-1)*[1-x]^-1*z
       + 1/2*Li2(x*[1-x]*z^-1*[1-z]^-1)*z^-1
       + Li2(x*[1-x]*z^-1*[1-z]^-1)*z
       + 1/2*Li2(x*[1-x]*z^-1*[1-z]^-1)*x*z^-1
       - Li2(x*[1-x]*z^-1*[1-z]^-1)*x
       - 2*Li2(z)
       - Li2(z)*[1-x]^-1*z^-1
       + Li2(z)*[1-x]^-1
       - 1/2*Li2(z)*[1-x]^-1*z
       + Li2(z)*z^-1
       + Li2(z)*z
       )

       + T(u3) * ( 10
       + 4*[1-x]^-1*z^-1
       - 4*[1-x]^-1
       + 6*[1-x]^-1*z
       - 9*z^-1
       - 7*z
       + 5*x*z^-1
       - 5*x
       - 3/2*pi^2*[1-x]^-1*z^-1
       + 7/3*pi^2*[1-x]^-1
       - 7/6*pi^2*[1-x]^-1*z
       + 13/12*pi^2*z^-1
       - 13/6*pi^2
       + 3/2*pi^2*z
       + 5/12*pi^2*x*z^-1
       - 5/6*pi^2*x
       - 27*ln(x)
       - 18*ln(x)*[1-x]^-1*z^-1
       + 24*ln(x)*[1-x]^-1
       - 15*ln(x)*[1-x]^-1*z
       + 18*ln(x)*z^-1
       + 18*ln(x)*z
       - 6*ln(x)*x
       + ln(x)*ln( - [x-z])
       + ln(x)*ln( - [x-z])*[1-x]^-1*z^-1
       - 2*ln(x)*ln( - [x-z])*[1-x]^-1
       + ln(x)*ln( - [x-z])*[1-x]^-1*z
       - 1/2*ln(x)*ln( - [x-z])*z^-1
       - ln(x)*ln( - [x-z])*z
       - 1/2*ln(x)*ln( - [x-z])*x*z^-1
       + ln(x)*ln( - [x-z])*x
       + 21/2*ln(x)^2
       + 7*ln(x)^2*[1-x]^-1*z^-1
       - 21/2*ln(x)^2*[1-x]^-1
       + 21/4*ln(x)^2*[1-x]^-1*z
       - 21/4*ln(x)^2*z^-1
       - 7*ln(x)^2*z
       - 7/4*ln(x)^2*x*z^-1
       + 7/2*ln(x)^2*x
       - 15*ln(x)*ln([1-x])
       - 9*ln(x)*ln([1-x])*[1-x]^-1*z^-1
       + 12*ln(x)*ln([1-x])*[1-x]^-1
       - 6*ln(x)*ln([1-x])*[1-x]^-1*z
       + 15/2*ln(x)*ln([1-x])*z^-1
       + 9*ln(x)*ln([1-x])*z
       + 3/2*ln(x)*ln([1-x])*x*z^-1
       - 3*ln(x)*ln([1-x])*x
       - 12*ln(x)*ln(z)
       - 9*ln(x)*ln(z)*[1-x]^-1*z^-1
       + 15*ln(x)*ln(z)*[1-x]^-1
       - 15/2*ln(x)*ln(z)*[1-x]^-1*z
       + 6*ln(x)*ln(z)*z^-1
       + 9*ln(x)*ln(z)*z
       + 3*ln(x)*ln(z)*x*z^-1
       - 6*ln(x)*ln(z)*x
       - 18*ln(x)*ln([1-z])
       - 12*ln(x)*ln([1-z])*[1-x]^-1*z^-1
       + 18*ln(x)*ln([1-z])*[1-x]^-1
       - 9*ln(x)*ln([1-z])*[1-x]^-1*z
       + 9*ln(x)*ln([1-z])*z^-1
       + 12*ln(x)*ln([1-z])*z
       + 3*ln(x)*ln([1-z])*x*z^-1
       - 6*ln(x)*ln([1-z])*x
       + 3*ln(x)*ln([1-x-z])
       + 2*ln(x)*ln([1-x-z])*[1-x]^-1*z^-1
       - 3*ln(x)*ln([1-x-z])*[1-x]^-1
       + 3/2*ln(x)*ln([1-x-z])*[1-x]^-1*z
       - 3/2*ln(x)*ln([1-x-z])*z^-1
       - 2*ln(x)*ln([1-x-z])*z
       - 1/2*ln(x)*ln([1-x-z])*x*z^-1
       + ln(x)*ln([1-x-z])*x
       + 13*ln([1-x])
       + 9*ln([1-x])*[1-x]^-1*z^-1
       - 11*ln([1-x])*[1-x]^-1
       + 7*ln([1-x])*[1-x]^-1*z
       - 9*ln([1-x])*z^-1
       - 9*ln([1-x])*z
       + 3*ln([1-x])*x
       + 11/2*ln([1-x])^2
       + 7/2*ln([1-x])^2*[1-x]^-1*z^-1
       - 5*ln([1-x])^2*[1-x]^-1
       + 5/2*ln([1-x])^2*[1-x]^-1*z
       - 11/4*ln([1-x])^2*z^-1
       - 7/2*ln([1-x])^2*z
       - 3/4*ln([1-x])^2*x*z^-1
       + 3/2*ln([1-x])^2*x
       + 5*ln([1-x])*ln(z)
       + 3*ln([1-x])*ln(z)*[1-x]^-1*z^-1
       - 4*ln([1-x])*ln(z)*[1-x]^-1
       + 2*ln([1-x])*ln(z)*[1-x]^-1*z
       - 5/2*ln([1-x])*ln(z)*z^-1
       - 3*ln([1-x])*ln(z)*z
       - 1/2*ln([1-x])*ln(z)*x*z^-1
       + ln([1-x])*ln(z)*x
       + 11*ln([1-x])*ln([1-z])
       + 6*ln([1-x])*ln([1-z])*[1-x]^-1*z^-1
       - 7*ln([1-x])*ln([1-z])*[1-x]^-1
       + 7/2*ln([1-x])*ln([1-z])*[1-x]^-1*z
       - 11/2*ln([1-x])*ln([1-z])*z^-1
       - 6*ln([1-x])*ln([1-z])*z
       - 1/2*ln([1-x])*ln([1-z])*x*z^-1
       + ln([1-x])*ln([1-z])*x
       + 14*ln(z)
       + 9*ln(z)*[1-x]^-1*z^-1
       - 13*ln(z)*[1-x]^-1
       + 8*ln(z)*[1-x]^-1*z
       - 9*ln(z)*z^-1
       - 9*ln(z)*z
       + 3*ln(z)*x
       - ln(z)*ln( - [x-z])
       - ln(z)*ln( - [x-z])*[1-x]^-1*z^-1
       + 2*ln(z)*ln( - [x-z])*[1-x]^-1
       - ln(z)*ln( - [x-z])*[1-x]^-1*z
       + 1/2*ln(z)*ln( - [x-z])*z^-1
       + ln(z)*ln( - [x-z])*z
       + 1/2*ln(z)*ln( - [x-z])*x*z^-1
       - ln(z)*ln( - [x-z])*x
       + 7/2*ln(z)^2
       + 3*ln(z)^2*[1-x]^-1*z^-1
       - 11/2*ln(z)^2*[1-x]^-1
       + 11/4*ln(z)^2*[1-x]^-1*z
       - 7/4*ln(z)^2*z^-1
       - 3*ln(z)^2*z
       - 5/4*ln(z)^2*x*z^-1
       + 5/2*ln(z)^2*x
       + 9*ln(z)*ln([1-z])
       + 7*ln(z)*ln([1-z])*[1-x]^-1*z^-1
       - 12*ln(z)*ln([1-z])*[1-x]^-1
       + 6*ln(z)*ln([1-z])*[1-x]^-1*z
       - 9/2*ln(z)*ln([1-z])*z^-1
       - 7*ln(z)*ln([1-z])*z
       - 5/2*ln(z)*ln([1-z])*x*z^-1
       + 5*ln(z)*ln([1-z])*x
       + 18*ln([1-z])
       + 12*ln([1-z])*[1-x]^-1*z^-1
       - 16*ln([1-z])*[1-x]^-1
       + 10*ln([1-z])*[1-x]^-1*z
       - 12*ln([1-z])*z^-1
       - 12*ln([1-z])*z
       + 4*ln([1-z])*x
       + 15/2*ln([1-z])^2
       + 5*ln([1-z])^2*[1-x]^-1*z^-1
       - 15/2*ln([1-z])^2*[1-x]^-1
       + 15/4*ln([1-z])^2*[1-x]^-1*z
       - 15/4*ln([1-z])^2*z^-1
       - 5*ln([1-z])^2*z
       - 5/4*ln([1-z])^2*x*z^-1
       + 5/2*ln([1-z])^2*x
       - 3*ln([1-z])*ln([1-x-z])
       - 2*ln([1-z])*ln([1-x-z])*[1-x]^-1*z^-1
       + 3*ln([1-z])*ln([1-x-z])*[1-x]^-1
       - 3/2*ln([1-z])*ln([1-x-z])*[1-x]^-1*z
       + 3/2*ln([1-z])*ln([1-x-z])*z^-1
       + 2*ln([1-z])*ln([1-x-z])*z
       + 1/2*ln([1-z])*ln([1-x-z])*x*z^-1
       - ln([1-z])*ln([1-x-z])*x
       + Li2([1-x]^-1*[1-z])
       + Li2([1-x]^-1*[1-z])*[1-x]^-1*z^-1
       - 2*Li2([1-x]^-1*[1-z])*[1-x]^-1
       + Li2([1-x]^-1*[1-z])*[1-x]^-1*z
       - 1/2*Li2([1-x]^-1*[1-z])*z^-1
       - Li2([1-x]^-1*[1-z])*z
       - 1/2*Li2([1-x]^-1*[1-z])*x*z^-1
       + Li2([1-x]^-1*[1-z])*x
       - Li2([1-x]^-1*z)
       - Li2([1-x]^-1*z)*[1-x]^-1
       + 1/2*Li2([1-x]^-1*z)*[1-x]^-1*z
       + 1/2*Li2([1-x]^-1*z)*z^-1
       - 1/2*Li2([1-x]^-1*z)*x*z^-1
       + Li2([1-x]^-1*z)*x
       + 2*Li2(x*[1-x]^-1*z*[1-z]^-1)
       + Li2(x*[1-x]^-1*z*[1-z]^-1)*[1-x]^-1*z^-1
       - Li2(x*[1-x]^-1*z*[1-z]^-1)*[1-x]^-1
       + 1/2*Li2(x*[1-x]^-1*z*[1-z]^-1)*[1-x]^-1*z
       - Li2(x*[1-x]^-1*z*[1-z]^-1)*z^-1
       - Li2(x*[1-x]^-1*z*[1-z]^-1)*z
       + Li2(x*[1-x]*z^-1*[1-z]^-1)
       + Li2(x*[1-x]*z^-1*[1-z]^-1)*[1-x]^-1*z^-1
       - 2*Li2(x*[1-x]*z^-1*[1-z]^-1)*[1-x]^-1
       + Li2(x*[1-x]*z^-1*[1-z]^-1)*[1-x]^-1*z
       - 1/2*Li2(x*[1-x]*z^-1*[1-z]^-1)*z^-1
       - Li2(x*[1-x]*z^-1*[1-z]^-1)*z
       - 1/2*Li2(x*[1-x]*z^-1*[1-z]^-1)*x*z^-1
       + Li2(x*[1-x]*z^-1*[1-z]^-1)*x
       + 2*Li2(z)
       + Li2(z)*[1-x]^-1*z^-1
       - Li2(z)*[1-x]^-1
       + 1/2*Li2(z)*[1-x]^-1*z
       - Li2(z)*z^-1
       - Li2(z)*z
       )

       + T(u4)*NC^-2 * (  - 10
       - 4*[1-x]^-1*z^-1
       + 4*[1-x]^-1
       - 6*[1-x]^-1*z
       + 9*z^-1
       + 7*z
       - 5*x*z^-1
       + 5*x
       + 5/6*pi^2*[1-x]^-1*z^-1
       - pi^2*[1-x]^-1
       + 1/2*pi^2*[1-x]^-1*z
       - 3/4*pi^2*z^-1
       + 3/2*pi^2
       - 5/6*pi^2*z
       - 1/12*pi^2*x*z^-1
       + 1/6*pi^2*x
       + 27*ln(x)
       + 18*ln(x)*[1-x]^-1*z^-1
       - 24*ln(x)*[1-x]^-1
       + 15*ln(x)*[1-x]^-1*z
       - 18*ln(x)*z^-1
       - 18*ln(x)*z
       + 6*ln(x)*x
       - ln(x)*ln( - [x-z])
       - ln(x)*ln( - [x-z])*[1-x]^-1*z^-1
       + 2*ln(x)*ln( - [x-z])*[1-x]^-1
       - ln(x)*ln( - [x-z])*[1-x]^-1*z
       + 1/2*ln(x)*ln( - [x-z])*z^-1
       + ln(x)*ln( - [x-z])*z
       + 1/2*ln(x)*ln( - [x-z])*x*z^-1
       - ln(x)*ln( - [x-z])*x
       - 3*ln(x)*ln( - [1-x-z])
       - 2*ln(x)*ln( - [1-x-z])*[1-x]^-1*z^-1
       + 3*ln(x)*ln( - [1-x-z])*[1-x]^-1
       - 3/2*ln(x)*ln( - [1-x-z])*[1-x]^-1*z
       + 3/2*ln(x)*ln( - [1-x-z])*z^-1
       + 2*ln(x)*ln( - [1-x-z])*z
       + 1/2*ln(x)*ln( - [1-x-z])*x*z^-1
       - ln(x)*ln( - [1-x-z])*x
       - 9*ln(x)^2
       - 6*ln(x)^2*[1-x]^-1*z^-1
       + 9*ln(x)^2*[1-x]^-1
       - 9/2*ln(x)^2*[1-x]^-1*z
       + 9/2*ln(x)^2*z^-1
       + 6*ln(x)^2*z
       + 3/2*ln(x)^2*x*z^-1
       - 3*ln(x)^2*x
       + 14*ln(x)*ln([1-x])
       + 9*ln(x)*ln([1-x])*[1-x]^-1*z^-1
       - 13*ln(x)*ln([1-x])*[1-x]^-1
       + 13/2*ln(x)*ln([1-x])*[1-x]^-1*z
       - 7*ln(x)*ln([1-x])*z^-1
       - 9*ln(x)*ln([1-x])*z
       - 2*ln(x)*ln([1-x])*x*z^-1
       + 4*ln(x)*ln([1-x])*x
       + 13*ln(x)*ln(z)
       + 9*ln(x)*ln(z)*[1-x]^-1*z^-1
       - 14*ln(x)*ln(z)*[1-x]^-1
       + 7*ln(x)*ln(z)*[1-x]^-1*z
       - 13/2*ln(x)*ln(z)*z^-1
       - 9*ln(x)*ln(z)*z
       - 5/2*ln(x)*ln(z)*x*z^-1
       + 5*ln(x)*ln(z)*x
       + 15*ln(x)*ln([1-z])
       + 10*ln(x)*ln([1-z])*[1-x]^-1*z^-1
       - 15*ln(x)*ln([1-z])*[1-x]^-1
       + 15/2*ln(x)*ln([1-z])*[1-x]^-1*z
       - 15/2*ln(x)*ln([1-z])*z^-1
       - 10*ln(x)*ln([1-z])*z
       - 5/2*ln(x)*ln([1-z])*x*z^-1
       + 5*ln(x)*ln([1-z])*x
       - 13*ln([1-x])
       - 9*ln([1-x])*[1-x]^-1*z^-1
       + 11*ln([1-x])*[1-x]^-1
       - 7*ln([1-x])*[1-x]^-1*z
       + 9*ln([1-x])*z^-1
       + 9*ln([1-x])*z
       - 3*ln([1-x])*x
       - 9/2*ln([1-x])^2
       - 5/2*ln([1-x])^2*[1-x]^-1*z^-1
       + 3*ln([1-x])^2*[1-x]^-1
       - 3/2*ln([1-x])^2*[1-x]^-1*z
       + 9/4*ln([1-x])^2*z^-1
       + 5/2*ln([1-x])^2*z
       + 1/4*ln([1-x])^2*x*z^-1
       - 1/2*ln([1-x])^2*x
       - 7*ln([1-x])*ln(z)
       - 5*ln([1-x])*ln(z)*[1-x]^-1*z^-1
       + 8*ln([1-x])*ln(z)*[1-x]^-1
       - 4*ln([1-x])*ln(z)*[1-x]^-1*z
       + 7/2*ln([1-x])*ln(z)*z^-1
       + 5*ln([1-x])*ln(z)*z
       + 3/2*ln([1-x])*ln(z)*x*z^-1
       - 3*ln([1-x])*ln(z)*x
       - 10*ln([1-x])*ln([1-z])
       - 6*ln([1-x])*ln([1-z])*[1-x]^-1*z^-1
       + 8*ln([1-x])*ln([1-z])*[1-x]^-1
       - 4*ln([1-x])*ln([1-z])*[1-x]^-1*z
       + 5*ln([1-x])*ln([1-z])*z^-1
       + 6*ln([1-x])*ln([1-z])*z
       + ln([1-x])*ln([1-z])*x*z^-1
       - 2*ln([1-x])*ln([1-z])*x
       - 14*ln(z)
       - 9*ln(z)*[1-x]^-1*z^-1
       + 13*ln(z)*[1-x]^-1
       - 8*ln(z)*[1-x]^-1*z
       + 9*ln(z)*z^-1
       + 9*ln(z)*z
       - 3*ln(z)*x
       + ln(z)*ln( - [x-z])
       + ln(z)*ln( - [x-z])*[1-x]^-1*z^-1
       - 2*ln(z)*ln( - [x-z])*[1-x]^-1
       + ln(z)*ln( - [x-z])*[1-x]^-1*z
       - 1/2*ln(z)*ln( - [x-z])*z^-1
       - ln(z)*ln( - [x-z])*z
       - 1/2*ln(z)*ln( - [x-z])*x*z^-1
       + ln(z)*ln( - [x-z])*x
       - 5/2*ln(z)^2
       - 2*ln(z)^2*[1-x]^-1*z^-1
       + 7/2*ln(z)^2*[1-x]^-1
       - 7/4*ln(z)^2*[1-x]^-1*z
       + 5/4*ln(z)^2*z^-1
       + 2*ln(z)^2*z
       + 3/4*ln(z)^2*x*z^-1
       - 3/2*ln(z)^2*x
       - 10*ln(z)*ln([1-z])
       - 7*ln(z)*ln([1-z])*[1-x]^-1*z^-1
       + 11*ln(z)*ln([1-z])*[1-x]^-1
       - 11/2*ln(z)*ln([1-z])*[1-x]^-1*z
       + 5*ln(z)*ln([1-z])*z^-1
       + 7*ln(z)*ln([1-z])*z
       + 2*ln(z)*ln([1-z])*x*z^-1
       - 4*ln(z)*ln([1-z])*x
       - 18*ln([1-z])
       - 12*ln([1-z])*[1-x]^-1*z^-1
       + 16*ln([1-z])*[1-x]^-1
       - 10*ln([1-z])*[1-x]^-1*z
       + 12*ln([1-z])*z^-1
       + 12*ln([1-z])*z
       - 4*ln([1-z])*x
       + 3*ln([1-z])*ln( - [1-x-z])
       + 2*ln([1-z])*ln( - [1-x-z])*[1-x]^-1*z^-1
       - 3*ln([1-z])*ln( - [1-x-z])*[1-x]^-1
       + 3/2*ln([1-z])*ln( - [1-x-z])*[1-x]^-1*z
       - 3/2*ln([1-z])*ln( - [1-x-z])*z^-1
       - 2*ln([1-z])*ln( - [1-x-z])*z
       - 1/2*ln([1-z])*ln( - [1-x-z])*x*z^-1
       + ln([1-z])*ln( - [1-x-z])*x
       - 6*ln([1-z])^2
       - 4*ln([1-z])^2*[1-x]^-1*z^-1
       + 6*ln([1-z])^2*[1-x]^-1
       - 3*ln([1-z])^2*[1-x]^-1*z
       + 3*ln([1-z])^2*z^-1
       + 4*ln([1-z])^2*z
       + ln([1-z])^2*x*z^-1
       - 2*ln([1-z])^2*x
       + Li2(x^-1*[1-x]^-1*z*[1-z])
       + Li2(x^-1*[1-x]^-1*z*[1-z])*[1-x]^-1*z^-1
       - 2*Li2(x^-1*[1-x]^-1*z*[1-z])*[1-x]^-1
       + Li2(x^-1*[1-x]^-1*z*[1-z])*[1-x]^-1*z
       - 1/2*Li2(x^-1*[1-x]^-1*z*[1-z])*z^-1
       - Li2(x^-1*[1-x]^-1*z*[1-z])*z
       - 1/2*Li2(x^-1*[1-x]^-1*z*[1-z])*x*z^-1
       + Li2(x^-1*[1-x]^-1*z*[1-z])*x
       + 2*Li2(x^-1*[1-x]*z^-1*[1-z])
       + Li2(x^-1*[1-x]*z^-1*[1-z])*[1-x]^-1*z^-1
       - Li2(x^-1*[1-x]*z^-1*[1-z])*[1-x]^-1
       + 1/2*Li2(x^-1*[1-x]*z^-1*[1-z])*[1-x]^-1*z
       - Li2(x^-1*[1-x]*z^-1*[1-z])*z^-1
       - Li2(x^-1*[1-x]*z^-1*[1-z])*z
       - Li2([1-x]^-1*[1-z])
       - Li2([1-x]^-1*[1-z])*[1-x]^-1*z^-1
       + 2*Li2([1-x]^-1*[1-z])*[1-x]^-1
       - Li2([1-x]^-1*[1-z])*[1-x]^-1*z
       + 1/2*Li2([1-x]^-1*[1-z])*z^-1
       + Li2([1-x]^-1*[1-z])*z
       + 1/2*Li2([1-x]^-1*[1-z])*x*z^-1
       - Li2([1-x]^-1*[1-z])*x
       - Li2([1-x]*z^-1)
       - Li2([1-x]*z^-1)*[1-x]^-1
       + 1/2*Li2([1-x]*z^-1)*[1-x]^-1*z
       + 1/2*Li2([1-x]*z^-1)*z^-1
       - 1/2*Li2([1-x]*z^-1)*x*z^-1
       + Li2([1-x]*z^-1)*x
       - 2*Li2(z)
       - Li2(z)*[1-x]^-1*z^-1
       + Li2(z)*[1-x]^-1
       - 1/2*Li2(z)*[1-x]^-1*z
       + Li2(z)*z^-1
       + Li2(z)*z
       )

       + T(u4) * ( 10
       + 4*[1-x]^-1*z^-1
       - 4*[1-x]^-1
       + 6*[1-x]^-1*z
       - 9*z^-1
       - 7*z
       + 5*x*z^-1
       - 5*x
       - 5/6*pi^2*[1-x]^-1*z^-1
       + pi^2*[1-x]^-1
       - 1/2*pi^2*[1-x]^-1*z
       + 3/4*pi^2*z^-1
       - 3/2*pi^2
       + 5/6*pi^2*z
       + 1/12*pi^2*x*z^-1
       - 1/6*pi^2*x
       - 27*ln(x)
       - 18*ln(x)*[1-x]^-1*z^-1
       + 24*ln(x)*[1-x]^-1
       - 15*ln(x)*[1-x]^-1*z
       + 18*ln(x)*z^-1
       + 18*ln(x)*z
       - 6*ln(x)*x
       + ln(x)*ln( - [x-z])
       + ln(x)*ln( - [x-z])*[1-x]^-1*z^-1
       - 2*ln(x)*ln( - [x-z])*[1-x]^-1
       + ln(x)*ln( - [x-z])*[1-x]^-1*z
       - 1/2*ln(x)*ln( - [x-z])*z^-1
       - ln(x)*ln( - [x-z])*z
       - 1/2*ln(x)*ln( - [x-z])*x*z^-1
       + ln(x)*ln( - [x-z])*x
       + 3*ln(x)*ln( - [1-x-z])
       + 2*ln(x)*ln( - [1-x-z])*[1-x]^-1*z^-1
       - 3*ln(x)*ln( - [1-x-z])*[1-x]^-1
       + 3/2*ln(x)*ln( - [1-x-z])*[1-x]^-1*z
       - 3/2*ln(x)*ln( - [1-x-z])*z^-1
       - 2*ln(x)*ln( - [1-x-z])*z
       - 1/2*ln(x)*ln( - [1-x-z])*x*z^-1
       + ln(x)*ln( - [1-x-z])*x
       + 9*ln(x)^2
       + 6*ln(x)^2*[1-x]^-1*z^-1
       - 9*ln(x)^2*[1-x]^-1
       + 9/2*ln(x)^2*[1-x]^-1*z
       - 9/2*ln(x)^2*z^-1
       - 6*ln(x)^2*z
       - 3/2*ln(x)^2*x*z^-1
       + 3*ln(x)^2*x
       - 14*ln(x)*ln([1-x])
       - 9*ln(x)*ln([1-x])*[1-x]^-1*z^-1
       + 13*ln(x)*ln([1-x])*[1-x]^-1
       - 13/2*ln(x)*ln([1-x])*[1-x]^-1*z
       + 7*ln(x)*ln([1-x])*z^-1
       + 9*ln(x)*ln([1-x])*z
       + 2*ln(x)*ln([1-x])*x*z^-1
       - 4*ln(x)*ln([1-x])*x
       - 13*ln(x)*ln(z)
       - 9*ln(x)*ln(z)*[1-x]^-1*z^-1
       + 14*ln(x)*ln(z)*[1-x]^-1
       - 7*ln(x)*ln(z)*[1-x]^-1*z
       + 13/2*ln(x)*ln(z)*z^-1
       + 9*ln(x)*ln(z)*z
       + 5/2*ln(x)*ln(z)*x*z^-1
       - 5*ln(x)*ln(z)*x
       - 15*ln(x)*ln([1-z])
       - 10*ln(x)*ln([1-z])*[1-x]^-1*z^-1
       + 15*ln(x)*ln([1-z])*[1-x]^-1
       - 15/2*ln(x)*ln([1-z])*[1-x]^-1*z
       + 15/2*ln(x)*ln([1-z])*z^-1
       + 10*ln(x)*ln([1-z])*z
       + 5/2*ln(x)*ln([1-z])*x*z^-1
       - 5*ln(x)*ln([1-z])*x
       + 13*ln([1-x])
       + 9*ln([1-x])*[1-x]^-1*z^-1
       - 11*ln([1-x])*[1-x]^-1
       + 7*ln([1-x])*[1-x]^-1*z
       - 9*ln([1-x])*z^-1
       - 9*ln([1-x])*z
       + 3*ln([1-x])*x
       + 9/2*ln([1-x])^2
       + 5/2*ln([1-x])^2*[1-x]^-1*z^-1
       - 3*ln([1-x])^2*[1-x]^-1
       + 3/2*ln([1-x])^2*[1-x]^-1*z
       - 9/4*ln([1-x])^2*z^-1
       - 5/2*ln([1-x])^2*z
       - 1/4*ln([1-x])^2*x*z^-1
       + 1/2*ln([1-x])^2*x
       + 7*ln([1-x])*ln(z)
       + 5*ln([1-x])*ln(z)*[1-x]^-1*z^-1
       - 8*ln([1-x])*ln(z)*[1-x]^-1
       + 4*ln([1-x])*ln(z)*[1-x]^-1*z
       - 7/2*ln([1-x])*ln(z)*z^-1
       - 5*ln([1-x])*ln(z)*z
       - 3/2*ln([1-x])*ln(z)*x*z^-1
       + 3*ln([1-x])*ln(z)*x
       + 10*ln([1-x])*ln([1-z])
       + 6*ln([1-x])*ln([1-z])*[1-x]^-1*z^-1
       - 8*ln([1-x])*ln([1-z])*[1-x]^-1
       + 4*ln([1-x])*ln([1-z])*[1-x]^-1*z
       - 5*ln([1-x])*ln([1-z])*z^-1
       - 6*ln([1-x])*ln([1-z])*z
       - ln([1-x])*ln([1-z])*x*z^-1
       + 2*ln([1-x])*ln([1-z])*x
       + 14*ln(z)
       + 9*ln(z)*[1-x]^-1*z^-1
       - 13*ln(z)*[1-x]^-1
       + 8*ln(z)*[1-x]^-1*z
       - 9*ln(z)*z^-1
       - 9*ln(z)*z
       + 3*ln(z)*x
       - ln(z)*ln( - [x-z])
       - ln(z)*ln( - [x-z])*[1-x]^-1*z^-1
       + 2*ln(z)*ln( - [x-z])*[1-x]^-1
       - ln(z)*ln( - [x-z])*[1-x]^-1*z
       + 1/2*ln(z)*ln( - [x-z])*z^-1
       + ln(z)*ln( - [x-z])*z
       + 1/2*ln(z)*ln( - [x-z])*x*z^-1
       - ln(z)*ln( - [x-z])*x
       + 5/2*ln(z)^2
       + 2*ln(z)^2*[1-x]^-1*z^-1
       - 7/2*ln(z)^2*[1-x]^-1
       + 7/4*ln(z)^2*[1-x]^-1*z
       - 5/4*ln(z)^2*z^-1
       - 2*ln(z)^2*z
       - 3/4*ln(z)^2*x*z^-1
       + 3/2*ln(z)^2*x
       + 10*ln(z)*ln([1-z])
       + 7*ln(z)*ln([1-z])*[1-x]^-1*z^-1
       - 11*ln(z)*ln([1-z])*[1-x]^-1
       + 11/2*ln(z)*ln([1-z])*[1-x]^-1*z
       - 5*ln(z)*ln([1-z])*z^-1
       - 7*ln(z)*ln([1-z])*z
       - 2*ln(z)*ln([1-z])*x*z^-1
       + 4*ln(z)*ln([1-z])*x
       + 18*ln([1-z])
       + 12*ln([1-z])*[1-x]^-1*z^-1
       - 16*ln([1-z])*[1-x]^-1
       + 10*ln([1-z])*[1-x]^-1*z
       - 12*ln([1-z])*z^-1
       - 12*ln([1-z])*z
       + 4*ln([1-z])*x
       - 3*ln([1-z])*ln( - [1-x-z])
       - 2*ln([1-z])*ln( - [1-x-z])*[1-x]^-1*z^-1
       + 3*ln([1-z])*ln( - [1-x-z])*[1-x]^-1
       - 3/2*ln([1-z])*ln( - [1-x-z])*[1-x]^-1*z
       + 3/2*ln([1-z])*ln( - [1-x-z])*z^-1
       + 2*ln([1-z])*ln( - [1-x-z])*z
       + 1/2*ln([1-z])*ln( - [1-x-z])*x*z^-1
       - ln([1-z])*ln( - [1-x-z])*x
       + 6*ln([1-z])^2
       + 4*ln([1-z])^2*[1-x]^-1*z^-1
       - 6*ln([1-z])^2*[1-x]^-1
       + 3*ln([1-z])^2*[1-x]^-1*z
       - 3*ln([1-z])^2*z^-1
       - 4*ln([1-z])^2*z
       - ln([1-z])^2*x*z^-1
       + 2*ln([1-z])^2*x
       - Li2(x^-1*[1-x]^-1*z*[1-z])
       - Li2(x^-1*[1-x]^-1*z*[1-z])*[1-x]^-1*z^-1
       + 2*Li2(x^-1*[1-x]^-1*z*[1-z])*[1-x]^-1
       - Li2(x^-1*[1-x]^-1*z*[1-z])*[1-x]^-1*z
       + 1/2*Li2(x^-1*[1-x]^-1*z*[1-z])*z^-1
       + Li2(x^-1*[1-x]^-1*z*[1-z])*z
       + 1/2*Li2(x^-1*[1-x]^-1*z*[1-z])*x*z^-1
       - Li2(x^-1*[1-x]^-1*z*[1-z])*x
       - 2*Li2(x^-1*[1-x]*z^-1*[1-z])
       - Li2(x^-1*[1-x]*z^-1*[1-z])*[1-x]^-1*z^-1
       + Li2(x^-1*[1-x]*z^-1*[1-z])*[1-x]^-1
       - 1/2*Li2(x^-1*[1-x]*z^-1*[1-z])*[1-x]^-1*z
       + Li2(x^-1*[1-x]*z^-1*[1-z])*z^-1
       + Li2(x^-1*[1-x]*z^-1*[1-z])*z
       + Li2([1-x]^-1*[1-z])
       + Li2([1-x]^-1*[1-z])*[1-x]^-1*z^-1
       - 2*Li2([1-x]^-1*[1-z])*[1-x]^-1
       + Li2([1-x]^-1*[1-z])*[1-x]^-1*z
       - 1/2*Li2([1-x]^-1*[1-z])*z^-1
       - Li2([1-x]^-1*[1-z])*z
       - 1/2*Li2([1-x]^-1*[1-z])*x*z^-1
       + Li2([1-x]^-1*[1-z])*x
       + Li2([1-x]*z^-1)
       + Li2([1-x]*z^-1)*[1-x]^-1
       - 1/2*Li2([1-x]*z^-1)*[1-x]^-1*z
       - 1/2*Li2([1-x]*z^-1)*z^-1
       + 1/2*Li2([1-x]*z^-1)*x*z^-1
       - Li2([1-x]*z^-1)*x
       + 2*Li2(z)
       + Li2(z)*[1-x]^-1*z^-1
       - Li2(z)*[1-x]^-1
       + 1/2*Li2(z)*[1-x]^-1*z
       - Li2(z)*z^-1
       - Li2(z)*z
       )

       + T(t1)*NC^2 * (  - 9/2
       - 2*[1-x]^-1*z^-1
       + 2*[1-x]^-1
       - 5/2*[1-x]^-1*z
       + 7/2*z^-1
       + 2*z
       - 3/2*x*z^-1
       + 3*x
       + 1/2*pi^2*[1-x]^-1*z^-1
       - 1/2*pi^2*[1-x]^-1
       + 1/4*pi^2*[1-x]^-1*z
       - 1/4*pi^2*z^-1
       + 1/2*pi^2
       - 1/2*pi^2*z
       - 1/4*pi^2*x*z^-1
       + 1/2*pi^2*x
       + 15/2*ln(x)
       + 9*ln(x)*[1-x]^-1*z^-1
       - 9*ln(x)*[1-x]^-1
       + 6*ln(x)*[1-x]^-1*z
       - 9*ln(x)*z^-1
       - 9*ln(x)*z
       + 9*ln(x)*x
       - ln(x)*ln( - [x-z])
       - ln(x)*ln( - [x-z])*[1-x]^-1*z^-1
       + ln(x)*ln( - [x-z])*[1-x]^-1
       - 1/2*ln(x)*ln( - [x-z])*[1-x]^-1*z
       + 1/2*ln(x)*ln( - [x-z])*z^-1
       + ln(x)*ln( - [x-z])*z
       + 1/2*ln(x)*ln( - [x-z])*x*z^-1
       - ln(x)*ln( - [x-z])*x
       - 7/2*ln(x)^2
       - 7/2*ln(x)^2*[1-x]^-1*z^-1
       + 7/2*ln(x)^2*[1-x]^-1
       - 7/4*ln(x)^2*[1-x]^-1*z
       + 7/4*ln(x)^2*z^-1
       + 7/2*ln(x)^2*z
       + 7/4*ln(x)^2*x*z^-1
       - 7/2*ln(x)^2*x
       + 6*ln(x)*ln([1-x])
       + 6*ln(x)*ln([1-x])*[1-x]^-1*z^-1
       - 6*ln(x)*ln([1-x])*[1-x]^-1
       + 3*ln(x)*ln([1-x])*[1-x]^-1*z
       - 3*ln(x)*ln([1-x])*z^-1
       - 6*ln(x)*ln([1-x])*z
       - 3*ln(x)*ln([1-x])*x*z^-1
       + 6*ln(x)*ln([1-x])*x
       + 6*ln(x)*ln(z)
       + 6*ln(x)*ln(z)*[1-x]^-1*z^-1
       - 6*ln(x)*ln(z)*[1-x]^-1
       + 3*ln(x)*ln(z)*[1-x]^-1*z
       - 3*ln(x)*ln(z)*z^-1
       - 6*ln(x)*ln(z)*z
       - 3*ln(x)*ln(z)*x*z^-1
       + 6*ln(x)*ln(z)*x
       + 3*ln(x)*ln([1-z])
       + 3*ln(x)*ln([1-z])*[1-x]^-1*z^-1
       - 3*ln(x)*ln([1-z])*[1-x]^-1
       + 3/2*ln(x)*ln([1-z])*[1-x]^-1*z
       - 3/2*ln(x)*ln([1-z])*z^-1
       - 3*ln(x)*ln([1-z])*z
       - 3/2*ln(x)*ln([1-z])*x*z^-1
       + 3*ln(x)*ln([1-z])*x
       - 5*ln([1-x])
       - 6*ln([1-x])*[1-x]^-1*z^-1
       + 6*ln([1-x])*[1-x]^-1
       - 4*ln([1-x])*[1-x]^-1*z
       + 6*ln([1-x])*z^-1
       + 6*ln([1-x])*z
       - 6*ln([1-x])*x
       - 2*ln([1-x])^2
       - 2*ln([1-x])^2*[1-x]^-1*z^-1
       + 2*ln([1-x])^2*[1-x]^-1
       - ln([1-x])^2*[1-x]^-1*z
       + ln([1-x])^2*z^-1
       + 2*ln([1-x])^2*z
       + ln([1-x])^2*x*z^-1
       - 2*ln([1-x])^2*x
       - 5*ln([1-x])*ln(z)
       - 5*ln([1-x])*ln(z)*[1-x]^-1*z^-1
       + 5*ln([1-x])*ln(z)*[1-x]^-1
       - 5/2*ln([1-x])*ln(z)*[1-x]^-1*z
       + 5/2*ln([1-x])*ln(z)*z^-1
       + 5*ln([1-x])*ln(z)*z
       + 5/2*ln([1-x])*ln(z)*x*z^-1
       - 5*ln([1-x])*ln(z)*x
       - 2*ln([1-x])*ln([1-z])
       - 2*ln([1-x])*ln([1-z])*[1-x]^-1*z^-1
       + 2*ln([1-x])*ln([1-z])*[1-x]^-1
       - ln([1-x])*ln([1-z])*[1-x]^-1*z
       + ln([1-x])*ln([1-z])*z^-1
       + 2*ln([1-x])*ln([1-z])*z
       + ln([1-x])*ln([1-z])*x*z^-1
       - 2*ln([1-x])*ln([1-z])*x
       - 5*ln(z)
       - 6*ln(z)*[1-x]^-1*z^-1
       + 6*ln(z)*[1-x]^-1
       - 4*ln(z)*[1-x]^-1*z
       + 6*ln(z)*z^-1
       + 6*ln(z)*z
       - 6*ln(z)*x
       + ln(z)*ln( - [x-z])
       + ln(z)*ln( - [x-z])*[1-x]^-1*z^-1
       - ln(z)*ln( - [x-z])*[1-x]^-1
       + 1/2*ln(z)*ln( - [x-z])*[1-x]^-1*z
       - 1/2*ln(z)*ln( - [x-z])*z^-1
       - ln(z)*ln( - [x-z])*z
       - 1/2*ln(z)*ln( - [x-z])*x*z^-1
       + ln(z)*ln( - [x-z])*x
       - 5/2*ln(z)^2
       - 5/2*ln(z)^2*[1-x]^-1*z^-1
       + 5/2*ln(z)^2*[1-x]^-1
       - 5/4*ln(z)^2*[1-x]^-1*z
       + 5/4*ln(z)^2*z^-1
       + 5/2*ln(z)^2*z
       + 5/4*ln(z)^2*x*z^-1
       - 5/2*ln(z)^2*x
       - ln(z)*ln([1-z])
       - ln(z)*ln([1-z])*[1-x]^-1*z^-1
       + ln(z)*ln([1-z])*[1-x]^-1
       - 1/2*ln(z)*ln([1-z])*[1-x]^-1*z
       + 1/2*ln(z)*ln([1-z])*z^-1
       + ln(z)*ln([1-z])*z
       + 1/2*ln(z)*ln([1-z])*x*z^-1
       - ln(z)*ln([1-z])*x
       - 5/2*ln([1-z])
       - 3*ln([1-z])*[1-x]^-1*z^-1
       + 3*ln([1-z])*[1-x]^-1
       - 2*ln([1-z])*[1-x]^-1*z
       + 3*ln([1-z])*z^-1
       + 3*ln([1-z])*z
       - 3*ln([1-z])*x
       - 1/2*ln([1-z])^2
       - 1/2*ln([1-z])^2*[1-x]^-1*z^-1
       + 1/2*ln([1-z])^2*[1-x]^-1
       - 1/4*ln([1-z])^2*[1-x]^-1*z
       + 1/4*ln([1-z])^2*z^-1
       + 1/2*ln([1-z])^2*z
       + 1/4*ln([1-z])^2*x*z^-1
       - 1/2*ln([1-z])^2*x
       + Li2([1-x]^-1*[1-z])
       + Li2([1-x]^-1*[1-z])*[1-x]^-1*z^-1
       - Li2([1-x]^-1*[1-z])*[1-x]^-1
       + 1/2*Li2([1-x]^-1*[1-z])*[1-x]^-1*z
       - 1/2*Li2([1-x]^-1*[1-z])*z^-1
       - Li2([1-x]^-1*[1-z])*z
       - 1/2*Li2([1-x]^-1*[1-z])*x*z^-1
       + Li2([1-x]^-1*[1-z])*x
       - Li2(x*[1-x]^-1*z^-1*[1-z])
       - Li2(x*[1-x]^-1*z^-1*[1-z])*[1-x]^-1*z^-1
       + Li2(x*[1-x]^-1*z^-1*[1-z])*[1-x]^-1
       - 1/2*Li2(x*[1-x]^-1*z^-1*[1-z])*[1-x]^-1*z
       + 1/2*Li2(x*[1-x]^-1*z^-1*[1-z])*z^-1
       + Li2(x*[1-x]^-1*z^-1*[1-z])*z
       + 1/2*Li2(x*[1-x]^-1*z^-1*[1-z])*x*z^-1
       - Li2(x*[1-x]^-1*z^-1*[1-z])*x
       + Li2(z)
       + Li2(z)*[1-x]^-1*z^-1
       - Li2(z)*[1-x]^-1
       + 1/2*Li2(z)*[1-x]^-1*z
       - 1/2*Li2(z)*z^-1
       - Li2(z)*z
       - 1/2*Li2(z)*x*z^-1
       + Li2(z)*x
       )

       + T(t1) * ( 9/2
       + 2*[1-x]^-1*z^-1
       - 2*[1-x]^-1
       + 5/2*[1-x]^-1*z
       - 7/2*z^-1
       - 2*z
       + 3/2*x*z^-1
       - 3*x
       - 1/2*pi^2*[1-x]^-1*z^-1
       + 1/2*pi^2*[1-x]^-1
       - 1/4*pi^2*[1-x]^-1*z
       + 1/4*pi^2*z^-1
       - 1/2*pi^2
       + 1/2*pi^2*z
       + 1/4*pi^2*x*z^-1
       - 1/2*pi^2*x
       - 15/2*ln(x)
       - 9*ln(x)*[1-x]^-1*z^-1
       + 9*ln(x)*[1-x]^-1
       - 6*ln(x)*[1-x]^-1*z
       + 9*ln(x)*z^-1
       + 9*ln(x)*z
       - 9*ln(x)*x
       + ln(x)*ln( - [x-z])
       + ln(x)*ln( - [x-z])*[1-x]^-1*z^-1
       - ln(x)*ln( - [x-z])*[1-x]^-1
       + 1/2*ln(x)*ln( - [x-z])*[1-x]^-1*z
       - 1/2*ln(x)*ln( - [x-z])*z^-1
       - ln(x)*ln( - [x-z])*z
       - 1/2*ln(x)*ln( - [x-z])*x*z^-1
       + ln(x)*ln( - [x-z])*x
       + 7/2*ln(x)^2
       + 7/2*ln(x)^2*[1-x]^-1*z^-1
       - 7/2*ln(x)^2*[1-x]^-1
       + 7/4*ln(x)^2*[1-x]^-1*z
       - 7/4*ln(x)^2*z^-1
       - 7/2*ln(x)^2*z
       - 7/4*ln(x)^2*x*z^-1
       + 7/2*ln(x)^2*x
       - 6*ln(x)*ln([1-x])
       - 6*ln(x)*ln([1-x])*[1-x]^-1*z^-1
       + 6*ln(x)*ln([1-x])*[1-x]^-1
       - 3*ln(x)*ln([1-x])*[1-x]^-1*z
       + 3*ln(x)*ln([1-x])*z^-1
       + 6*ln(x)*ln([1-x])*z
       + 3*ln(x)*ln([1-x])*x*z^-1
       - 6*ln(x)*ln([1-x])*x
       - 6*ln(x)*ln(z)
       - 6*ln(x)*ln(z)*[1-x]^-1*z^-1
       + 6*ln(x)*ln(z)*[1-x]^-1
       - 3*ln(x)*ln(z)*[1-x]^-1*z
       + 3*ln(x)*ln(z)*z^-1
       + 6*ln(x)*ln(z)*z
       + 3*ln(x)*ln(z)*x*z^-1
       - 6*ln(x)*ln(z)*x
       - 3*ln(x)*ln([1-z])
       - 3*ln(x)*ln([1-z])*[1-x]^-1*z^-1
       + 3*ln(x)*ln([1-z])*[1-x]^-1
       - 3/2*ln(x)*ln([1-z])*[1-x]^-1*z
       + 3/2*ln(x)*ln([1-z])*z^-1
       + 3*ln(x)*ln([1-z])*z
       + 3/2*ln(x)*ln([1-z])*x*z^-1
       - 3*ln(x)*ln([1-z])*x
       + 5*ln([1-x])
       + 6*ln([1-x])*[1-x]^-1*z^-1
       - 6*ln([1-x])*[1-x]^-1
       + 4*ln([1-x])*[1-x]^-1*z
       - 6*ln([1-x])*z^-1
       - 6*ln([1-x])*z
       + 6*ln([1-x])*x
       + 2*ln([1-x])^2
       + 2*ln([1-x])^2*[1-x]^-1*z^-1
       - 2*ln([1-x])^2*[1-x]^-1
       + ln([1-x])^2*[1-x]^-1*z
       - ln([1-x])^2*z^-1
       - 2*ln([1-x])^2*z
       - ln([1-x])^2*x*z^-1
       + 2*ln([1-x])^2*x
       + 5*ln([1-x])*ln(z)
       + 5*ln([1-x])*ln(z)*[1-x]^-1*z^-1
       - 5*ln([1-x])*ln(z)*[1-x]^-1
       + 5/2*ln([1-x])*ln(z)*[1-x]^-1*z
       - 5/2*ln([1-x])*ln(z)*z^-1
       - 5*ln([1-x])*ln(z)*z
       - 5/2*ln([1-x])*ln(z)*x*z^-1
       + 5*ln([1-x])*ln(z)*x
       + 2*ln([1-x])*ln([1-z])
       + 2*ln([1-x])*ln([1-z])*[1-x]^-1*z^-1
       - 2*ln([1-x])*ln([1-z])*[1-x]^-1
       + ln([1-x])*ln([1-z])*[1-x]^-1*z
       - ln([1-x])*ln([1-z])*z^-1
       - 2*ln([1-x])*ln([1-z])*z
       - ln([1-x])*ln([1-z])*x*z^-1
       + 2*ln([1-x])*ln([1-z])*x
       + 5*ln(z)
       + 6*ln(z)*[1-x]^-1*z^-1
       - 6*ln(z)*[1-x]^-1
       + 4*ln(z)*[1-x]^-1*z
       - 6*ln(z)*z^-1
       - 6*ln(z)*z
       + 6*ln(z)*x
       - ln(z)*ln( - [x-z])
       - ln(z)*ln( - [x-z])*[1-x]^-1*z^-1
       + ln(z)*ln( - [x-z])*[1-x]^-1
       - 1/2*ln(z)*ln( - [x-z])*[1-x]^-1*z
       + 1/2*ln(z)*ln( - [x-z])*z^-1
       + ln(z)*ln( - [x-z])*z
       + 1/2*ln(z)*ln( - [x-z])*x*z^-1
       - ln(z)*ln( - [x-z])*x
       + 5/2*ln(z)^2
       + 5/2*ln(z)^2*[1-x]^-1*z^-1
       - 5/2*ln(z)^2*[1-x]^-1
       + 5/4*ln(z)^2*[1-x]^-1*z
       - 5/4*ln(z)^2*z^-1
       - 5/2*ln(z)^2*z
       - 5/4*ln(z)^2*x*z^-1
       + 5/2*ln(z)^2*x
       + ln(z)*ln([1-z])
       + ln(z)*ln([1-z])*[1-x]^-1*z^-1
       - ln(z)*ln([1-z])*[1-x]^-1
       + 1/2*ln(z)*ln([1-z])*[1-x]^-1*z
       - 1/2*ln(z)*ln([1-z])*z^-1
       - ln(z)*ln([1-z])*z
       - 1/2*ln(z)*ln([1-z])*x*z^-1
       + ln(z)*ln([1-z])*x
       + 5/2*ln([1-z])
       + 3*ln([1-z])*[1-x]^-1*z^-1
       - 3*ln([1-z])*[1-x]^-1
       + 2*ln([1-z])*[1-x]^-1*z
       - 3*ln([1-z])*z^-1
       - 3*ln([1-z])*z
       + 3*ln([1-z])*x
       + 1/2*ln([1-z])^2
       + 1/2*ln([1-z])^2*[1-x]^-1*z^-1
       - 1/2*ln([1-z])^2*[1-x]^-1
       + 1/4*ln([1-z])^2*[1-x]^-1*z
       - 1/4*ln([1-z])^2*z^-1
       - 1/2*ln([1-z])^2*z
       - 1/4*ln([1-z])^2*x*z^-1
       + 1/2*ln([1-z])^2*x
       - Li2([1-x]^-1*[1-z])
       - Li2([1-x]^-1*[1-z])*[1-x]^-1*z^-1
       + Li2([1-x]^-1*[1-z])*[1-x]^-1
       - 1/2*Li2([1-x]^-1*[1-z])*[1-x]^-1*z
       + 1/2*Li2([1-x]^-1*[1-z])*z^-1
       + Li2([1-x]^-1*[1-z])*z
       + 1/2*Li2([1-x]^-1*[1-z])*x*z^-1
       - Li2([1-x]^-1*[1-z])*x
       + Li2(x*[1-x]^-1*z^-1*[1-z])
       + Li2(x*[1-x]^-1*z^-1*[1-z])*[1-x]^-1*z^-1
       - Li2(x*[1-x]^-1*z^-1*[1-z])*[1-x]^-1
       + 1/2*Li2(x*[1-x]^-1*z^-1*[1-z])*[1-x]^-1*z
       - 1/2*Li2(x*[1-x]^-1*z^-1*[1-z])*z^-1
       - Li2(x*[1-x]^-1*z^-1*[1-z])*z
       - 1/2*Li2(x*[1-x]^-1*z^-1*[1-z])*x*z^-1
       + Li2(x*[1-x]^-1*z^-1*[1-z])*x
       - Li2(z)
       - Li2(z)*[1-x]^-1*z^-1
       + Li2(z)*[1-x]^-1
       - 1/2*Li2(z)*[1-x]^-1*z
       + 1/2*Li2(z)*z^-1
       + Li2(z)*z
       + 1/2*Li2(z)*x*z^-1
       - Li2(z)*x
       )

       + T(t2)*NC^2 * (  - 9/2
       - 2*[1-x]^-1*z^-1
       + 2*[1-x]^-1
       - 5/2*[1-x]^-1*z
       + 7/2*z^-1
       + 2*z
       - 3/2*x*z^-1
       + 3*x
       + 1/2*pi^2*[1-x]^-1*z^-1
       - 1/2*pi^2*[1-x]^-1
       + 1/4*pi^2*[1-x]^-1*z
       - 1/4*pi^2*z^-1
       + 1/2*pi^2
       - 1/2*pi^2*z
       - 1/4*pi^2*x*z^-1
       + 1/2*pi^2*x
       + 15/2*ln(x)
       + 9*ln(x)*[1-x]^-1*z^-1
       - 9*ln(x)*[1-x]^-1
       + 6*ln(x)*[1-x]^-1*z
       - 9*ln(x)*z^-1
       - 9*ln(x)*z
       + 9*ln(x)*x
       - 3*ln(x)^2
       - 3*ln(x)^2*[1-x]^-1*z^-1
       + 3*ln(x)^2*[1-x]^-1
       - 3/2*ln(x)^2*[1-x]^-1*z
       + 3/2*ln(x)^2*z^-1
       + 3*ln(x)^2*z
       + 3/2*ln(x)^2*x*z^-1
       - 3*ln(x)^2*x
       + 5*ln(x)*ln([1-x])
       + 5*ln(x)*ln([1-x])*[1-x]^-1*z^-1
       - 5*ln(x)*ln([1-x])*[1-x]^-1
       + 5/2*ln(x)*ln([1-x])*[1-x]^-1*z
       - 5/2*ln(x)*ln([1-x])*z^-1
       - 5*ln(x)*ln([1-x])*z
       - 5/2*ln(x)*ln([1-x])*x*z^-1
       + 5*ln(x)*ln([1-x])*x
       + 5*ln(x)*ln(z)
       + 5*ln(x)*ln(z)*[1-x]^-1*z^-1
       - 5*ln(x)*ln(z)*[1-x]^-1
       + 5/2*ln(x)*ln(z)*[1-x]^-1*z
       - 5/2*ln(x)*ln(z)*z^-1
       - 5*ln(x)*ln(z)*z
       - 5/2*ln(x)*ln(z)*x*z^-1
       + 5*ln(x)*ln(z)*x
       + 4*ln(x)*ln([1-z])
       + 4*ln(x)*ln([1-z])*[1-x]^-1*z^-1
       - 4*ln(x)*ln([1-z])*[1-x]^-1
       + 2*ln(x)*ln([1-z])*[1-x]^-1*z
       - 2*ln(x)*ln([1-z])*z^-1
       - 4*ln(x)*ln([1-z])*z
       - 2*ln(x)*ln([1-z])*x*z^-1
       + 4*ln(x)*ln([1-z])*x
       - ln(x)*ln([x-z])
       - ln(x)*ln([x-z])*[1-x]^-1*z^-1
       + ln(x)*ln([x-z])*[1-x]^-1
       - 1/2*ln(x)*ln([x-z])*[1-x]^-1*z
       + 1/2*ln(x)*ln([x-z])*z^-1
       + ln(x)*ln([x-z])*z
       + 1/2*ln(x)*ln([x-z])*x*z^-1
       - ln(x)*ln([x-z])*x
       - 5*ln([1-x])
       - 6*ln([1-x])*[1-x]^-1*z^-1
       + 6*ln([1-x])*[1-x]^-1
       - 4*ln([1-x])*[1-x]^-1*z
       + 6*ln([1-x])*z^-1
       + 6*ln([1-x])*z
       - 6*ln([1-x])*x
       - 2*ln([1-x])^2
       - 2*ln([1-x])^2*[1-x]^-1*z^-1
       + 2*ln([1-x])^2*[1-x]^-1
       - ln([1-x])^2*[1-x]^-1*z
       + ln([1-x])^2*z^-1
       + 2*ln([1-x])^2*z
       + ln([1-x])^2*x*z^-1
       - 2*ln([1-x])^2*x
       - 4*ln([1-x])*ln(z)
       - 4*ln([1-x])*ln(z)*[1-x]^-1*z^-1
       + 4*ln([1-x])*ln(z)*[1-x]^-1
       - 2*ln([1-x])*ln(z)*[1-x]^-1*z
       + 2*ln([1-x])*ln(z)*z^-1
       + 4*ln([1-x])*ln(z)*z
       + 2*ln([1-x])*ln(z)*x*z^-1
       - 4*ln([1-x])*ln(z)*x
       - 2*ln([1-x])*ln([1-z])
       - 2*ln([1-x])*ln([1-z])*[1-x]^-1*z^-1
       + 2*ln([1-x])*ln([1-z])*[1-x]^-1
       - ln([1-x])*ln([1-z])*[1-x]^-1*z
       + ln([1-x])*ln([1-z])*z^-1
       + 2*ln([1-x])*ln([1-z])*z
       + ln([1-x])*ln([1-z])*x*z^-1
       - 2*ln([1-x])*ln([1-z])*x
       - 5*ln(z)
       - 6*ln(z)*[1-x]^-1*z^-1
       + 6*ln(z)*[1-x]^-1
       - 4*ln(z)*[1-x]^-1*z
       + 6*ln(z)*z^-1
       + 6*ln(z)*z
       - 6*ln(z)*x
       - 2*ln(z)^2
       - 2*ln(z)^2*[1-x]^-1*z^-1
       + 2*ln(z)^2*[1-x]^-1
       - ln(z)^2*[1-x]^-1*z
       + ln(z)^2*z^-1
       + 2*ln(z)^2*z
       + ln(z)^2*x*z^-1
       - 2*ln(z)^2*x
       - 2*ln(z)*ln([1-z])
       - 2*ln(z)*ln([1-z])*[1-x]^-1*z^-1
       + 2*ln(z)*ln([1-z])*[1-x]^-1
       - ln(z)*ln([1-z])*[1-x]^-1*z
       + ln(z)*ln([1-z])*z^-1
       + 2*ln(z)*ln([1-z])*z
       + ln(z)*ln([1-z])*x*z^-1
       - 2*ln(z)*ln([1-z])*x
       + ln(z)*ln([x-z])
       + ln(z)*ln([x-z])*[1-x]^-1*z^-1
       - ln(z)*ln([x-z])*[1-x]^-1
       + 1/2*ln(z)*ln([x-z])*[1-x]^-1*z
       - 1/2*ln(z)*ln([x-z])*z^-1
       - ln(z)*ln([x-z])*z
       - 1/2*ln(z)*ln([x-z])*x*z^-1
       + ln(z)*ln([x-z])*x
       - 5/2*ln([1-z])
       - 3*ln([1-z])*[1-x]^-1*z^-1
       + 3*ln([1-z])*[1-x]^-1
       - 2*ln([1-z])*[1-x]^-1*z
       + 3*ln([1-z])*z^-1
       + 3*ln([1-z])*z
       - 3*ln([1-z])*x
       - 1/2*ln([1-z])^2
       - 1/2*ln([1-z])^2*[1-x]^-1*z^-1
       + 1/2*ln([1-z])^2*[1-x]^-1
       - 1/4*ln([1-z])^2*[1-x]^-1*z
       + 1/4*ln([1-z])^2*z^-1
       + 1/2*ln([1-z])^2*z
       + 1/4*ln([1-z])^2*x*z^-1
       - 1/2*ln([1-z])^2*x
       + Li2(x^-1*[1-x]*z*[1-z]^-1)
       + Li2(x^-1*[1-x]*z*[1-z]^-1)*[1-x]^-1*z^-1
       - Li2(x^-1*[1-x]*z*[1-z]^-1)*[1-x]^-1
       + 1/2*Li2(x^-1*[1-x]*z*[1-z]^-1)*[1-x]^-1*z
       - 1/2*Li2(x^-1*[1-x]*z*[1-z]^-1)*z^-1
       - Li2(x^-1*[1-x]*z*[1-z]^-1)*z
       - 1/2*Li2(x^-1*[1-x]*z*[1-z]^-1)*x*z^-1
       + Li2(x^-1*[1-x]*z*[1-z]^-1)*x
       - Li2([1-x]*[1-z]^-1)
       - Li2([1-x]*[1-z]^-1)*[1-x]^-1*z^-1
       + Li2([1-x]*[1-z]^-1)*[1-x]^-1
       - 1/2*Li2([1-x]*[1-z]^-1)*[1-x]^-1*z
       + 1/2*Li2([1-x]*[1-z]^-1)*z^-1
       + Li2([1-x]*[1-z]^-1)*z
       + 1/2*Li2([1-x]*[1-z]^-1)*x*z^-1
       - Li2([1-x]*[1-z]^-1)*x
       + Li2(z)
       + Li2(z)*[1-x]^-1*z^-1
       - Li2(z)*[1-x]^-1
       + 1/2*Li2(z)*[1-x]^-1*z
       - 1/2*Li2(z)*z^-1
       - Li2(z)*z
       - 1/2*Li2(z)*x*z^-1
       + Li2(z)*x
       )

       + T(t2) * ( 9/2
       + 2*[1-x]^-1*z^-1
       - 2*[1-x]^-1
       + 5/2*[1-x]^-1*z
       - 7/2*z^-1
       - 2*z
       + 3/2*x*z^-1
       - 3*x
       - 1/2*pi^2*[1-x]^-1*z^-1
       + 1/2*pi^2*[1-x]^-1
       - 1/4*pi^2*[1-x]^-1*z
       + 1/4*pi^2*z^-1
       - 1/2*pi^2
       + 1/2*pi^2*z
       + 1/4*pi^2*x*z^-1
       - 1/2*pi^2*x
       - 15/2*ln(x)
       - 9*ln(x)*[1-x]^-1*z^-1
       + 9*ln(x)*[1-x]^-1
       - 6*ln(x)*[1-x]^-1*z
       + 9*ln(x)*z^-1
       + 9*ln(x)*z
       - 9*ln(x)*x
       + 3*ln(x)^2
       + 3*ln(x)^2*[1-x]^-1*z^-1
       - 3*ln(x)^2*[1-x]^-1
       + 3/2*ln(x)^2*[1-x]^-1*z
       - 3/2*ln(x)^2*z^-1
       - 3*ln(x)^2*z
       - 3/2*ln(x)^2*x*z^-1
       + 3*ln(x)^2*x
       - 5*ln(x)*ln([1-x])
       - 5*ln(x)*ln([1-x])*[1-x]^-1*z^-1
       + 5*ln(x)*ln([1-x])*[1-x]^-1
       - 5/2*ln(x)*ln([1-x])*[1-x]^-1*z
       + 5/2*ln(x)*ln([1-x])*z^-1
       + 5*ln(x)*ln([1-x])*z
       + 5/2*ln(x)*ln([1-x])*x*z^-1
       - 5*ln(x)*ln([1-x])*x
       - 5*ln(x)*ln(z)
       - 5*ln(x)*ln(z)*[1-x]^-1*z^-1
       + 5*ln(x)*ln(z)*[1-x]^-1
       - 5/2*ln(x)*ln(z)*[1-x]^-1*z
       + 5/2*ln(x)*ln(z)*z^-1
       + 5*ln(x)*ln(z)*z
       + 5/2*ln(x)*ln(z)*x*z^-1
       - 5*ln(x)*ln(z)*x
       - 4*ln(x)*ln([1-z])
       - 4*ln(x)*ln([1-z])*[1-x]^-1*z^-1
       + 4*ln(x)*ln([1-z])*[1-x]^-1
       - 2*ln(x)*ln([1-z])*[1-x]^-1*z
       + 2*ln(x)*ln([1-z])*z^-1
       + 4*ln(x)*ln([1-z])*z
       + 2*ln(x)*ln([1-z])*x*z^-1
       - 4*ln(x)*ln([1-z])*x
       + ln(x)*ln([x-z])
       + ln(x)*ln([x-z])*[1-x]^-1*z^-1
       - ln(x)*ln([x-z])*[1-x]^-1
       + 1/2*ln(x)*ln([x-z])*[1-x]^-1*z
       - 1/2*ln(x)*ln([x-z])*z^-1
       - ln(x)*ln([x-z])*z
       - 1/2*ln(x)*ln([x-z])*x*z^-1
       + ln(x)*ln([x-z])*x
       + 5*ln([1-x])
       + 6*ln([1-x])*[1-x]^-1*z^-1
       - 6*ln([1-x])*[1-x]^-1
       + 4*ln([1-x])*[1-x]^-1*z
       - 6*ln([1-x])*z^-1
       - 6*ln([1-x])*z
       + 6*ln([1-x])*x
       + 2*ln([1-x])^2
       + 2*ln([1-x])^2*[1-x]^-1*z^-1
       - 2*ln([1-x])^2*[1-x]^-1
       + ln([1-x])^2*[1-x]^-1*z
       - ln([1-x])^2*z^-1
       - 2*ln([1-x])^2*z
       - ln([1-x])^2*x*z^-1
       + 2*ln([1-x])^2*x
       + 4*ln([1-x])*ln(z)
       + 4*ln([1-x])*ln(z)*[1-x]^-1*z^-1
       - 4*ln([1-x])*ln(z)*[1-x]^-1
       + 2*ln([1-x])*ln(z)*[1-x]^-1*z
       - 2*ln([1-x])*ln(z)*z^-1
       - 4*ln([1-x])*ln(z)*z
       - 2*ln([1-x])*ln(z)*x*z^-1
       + 4*ln([1-x])*ln(z)*x
       + 2*ln([1-x])*ln([1-z])
       + 2*ln([1-x])*ln([1-z])*[1-x]^-1*z^-1
       - 2*ln([1-x])*ln([1-z])*[1-x]^-1
       + ln([1-x])*ln([1-z])*[1-x]^-1*z
       - ln([1-x])*ln([1-z])*z^-1
       - 2*ln([1-x])*ln([1-z])*z
       - ln([1-x])*ln([1-z])*x*z^-1
       + 2*ln([1-x])*ln([1-z])*x
       + 5*ln(z)
       + 6*ln(z)*[1-x]^-1*z^-1
       - 6*ln(z)*[1-x]^-1
       + 4*ln(z)*[1-x]^-1*z
       - 6*ln(z)*z^-1
       - 6*ln(z)*z
       + 6*ln(z)*x
       + 2*ln(z)^2
       + 2*ln(z)^2*[1-x]^-1*z^-1
       - 2*ln(z)^2*[1-x]^-1
       + ln(z)^2*[1-x]^-1*z
       - ln(z)^2*z^-1
       - 2*ln(z)^2*z
       - ln(z)^2*x*z^-1
       + 2*ln(z)^2*x
       + 2*ln(z)*ln([1-z])
       + 2*ln(z)*ln([1-z])*[1-x]^-1*z^-1
       - 2*ln(z)*ln([1-z])*[1-x]^-1
       + ln(z)*ln([1-z])*[1-x]^-1*z
       - ln(z)*ln([1-z])*z^-1
       - 2*ln(z)*ln([1-z])*z
       - ln(z)*ln([1-z])*x*z^-1
       + 2*ln(z)*ln([1-z])*x
       - ln(z)*ln([x-z])
       - ln(z)*ln([x-z])*[1-x]^-1*z^-1
       + ln(z)*ln([x-z])*[1-x]^-1
       - 1/2*ln(z)*ln([x-z])*[1-x]^-1*z
       + 1/2*ln(z)*ln([x-z])*z^-1
       + ln(z)*ln([x-z])*z
       + 1/2*ln(z)*ln([x-z])*x*z^-1
       - ln(z)*ln([x-z])*x
       + 5/2*ln([1-z])
       + 3*ln([1-z])*[1-x]^-1*z^-1
       - 3*ln([1-z])*[1-x]^-1
       + 2*ln([1-z])*[1-x]^-1*z
       - 3*ln([1-z])*z^-1
       - 3*ln([1-z])*z
       + 3*ln([1-z])*x
       + 1/2*ln([1-z])^2
       + 1/2*ln([1-z])^2*[1-x]^-1*z^-1
       - 1/2*ln([1-z])^2*[1-x]^-1
       + 1/4*ln([1-z])^2*[1-x]^-1*z
       - 1/4*ln([1-z])^2*z^-1
       - 1/2*ln([1-z])^2*z
       - 1/4*ln([1-z])^2*x*z^-1
       + 1/2*ln([1-z])^2*x
       - Li2(x^-1*[1-x]*z*[1-z]^-1)
       - Li2(x^-1*[1-x]*z*[1-z]^-1)*[1-x]^-1*z^-1
       + Li2(x^-1*[1-x]*z*[1-z]^-1)*[1-x]^-1
       - 1/2*Li2(x^-1*[1-x]*z*[1-z]^-1)*[1-x]^-1*z
       + 1/2*Li2(x^-1*[1-x]*z*[1-z]^-1)*z^-1
       + Li2(x^-1*[1-x]*z*[1-z]^-1)*z
       + 1/2*Li2(x^-1*[1-x]*z*[1-z]^-1)*x*z^-1
       - Li2(x^-1*[1-x]*z*[1-z]^-1)*x
       + Li2([1-x]*[1-z]^-1)
       + Li2([1-x]*[1-z]^-1)*[1-x]^-1*z^-1
       - Li2([1-x]*[1-z]^-1)*[1-x]^-1
       + 1/2*Li2([1-x]*[1-z]^-1)*[1-x]^-1*z
       - 1/2*Li2([1-x]*[1-z]^-1)*z^-1
       - Li2([1-x]*[1-z]^-1)*z
       - 1/2*Li2([1-x]*[1-z]^-1)*x*z^-1
       + Li2([1-x]*[1-z]^-1)*x
       - Li2(z)
       - Li2(z)*[1-x]^-1*z^-1
       + Li2(z)*[1-x]^-1
       - 1/2*Li2(z)*[1-x]^-1*z
       + 1/2*Li2(z)*z^-1
       + Li2(z)*z
       + 1/2*Li2(z)*x*z^-1
       - Li2(z)*x
       )

       - 79/12
       - 6*[1-x]^-1*z^-1
       - 2*[1-x]^-1*z^-1*rln2^2
       - 3*[1-x]^-1*z^-1*sqrtxz1*rln2
       + 6*[1-x]^-1
       + 6*[1-x]^-1*rln2^2
       + 3*[1-x]^-1*sqrtxz1*rln2
       - 17/2*[1-x]^-1*z
       - [1-x]^-1*z*rln2^2
       + 199/36*z^-1
       + z^-1*rln2^2
       + z^-1*sqrtxz1*rln2
       - [1-x-z]^-1
       - 2*rln2^2
       - sqrtxz1*rln2
       + 2*z*[1-x-z]^-1
       + 85/12*z
       - 2*z*rln2^2
       - z^2*[1-x-z]^-1
       - 41/18*z^2
       - 245/36*x*z^-1
       + x*z^-1*rln2^2
       + 2*x*z^-1*sqrtxz1*rln2
       + 155/12*x
       - 2*x*rln2^2
       - 35/12*x*z
       + 19/18*x*z^2
       + 7/3*pi^2*[1-x]^-1*z^-1
       - 7/3*pi^2*[1-x]^-1
       + 4/3*pi^2*[1-x]^-1*z
       - 23/12*pi^2*z^-1
       + 5/2*pi^2
       - 5/2*pi^2*z
       - 5/4*pi^2*x*z^-1
       + 5/6*pi^2*x
       - 1/6*pi^2*x*z
       + 6*ln(1 + sqrtxz1 - z)*[1-x]^-1*z^-1*rln2
       + 3*ln(1 + sqrtxz1 - z)*[1-x]^-1*z^-1*sqrtxz1
       - 10*ln(1 + sqrtxz1 - z)*[1-x]^-1*rln2
       - 3*ln(1 + sqrtxz1 - z)*[1-x]^-1*sqrtxz1
       + 3*ln(1 + sqrtxz1 - z)*[1-x]^-1*z*rln2
       - 3*ln(1 + sqrtxz1 - z)*z^-1*rln2
       - ln(1 + sqrtxz1 - z)*z^-1*sqrtxz1
       + 2*ln(1 + sqrtxz1 - z)*rln2
       + ln(1 + sqrtxz1 - z)*sqrtxz1
       + 2*ln(1 + sqrtxz1 - z)*z*rln2
       - 3*ln(1 + sqrtxz1 - z)*x*z^-1*rln2
       - 2*ln(1 + sqrtxz1 - z)*x*z^-1*sqrtxz1
       + 2*ln(1 + sqrtxz1 - z)*x*rln2
       - 4*ln(1 + sqrtxz1 - z)^2*[1-x]^-1*z^-1
       + 4*ln(1 + sqrtxz1 - z)^2*[1-x]^-1
       - 2*ln(1 + sqrtxz1 - z)^2*[1-x]^-1*z
       + 2*ln(1 + sqrtxz1 - z)^2*z^-1
       + 2*ln(1 + sqrtxz1 - z)^2*x*z^-1
       - 2*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)
       + 2*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^-1
       + 2*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*[1-x]^-1
       + ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z
       - ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*z^-1
       - 2*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*z
       - ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*x*z^-1
       - 2*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*x
       - 2*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^-1*rln2
       - 2*ln(1 + sqrtxz1 + z)*[1-x]^-1*rln2
       - ln(1 + sqrtxz1 + z)*[1-x]^-1*z*rln2
       + ln(1 + sqrtxz1 + z)*z^-1*rln2
       + 2*ln(1 + sqrtxz1 + z)*rln2
       + 2*ln(1 + sqrtxz1 + z)*z*rln2
       + ln(1 + sqrtxz1 + z)*x*z^-1*rln2
       + 2*ln(1 + sqrtxz1 + z)*x*rln2
       + 34*ln(x)
       + 109/6*ln(x)*[1-x]^-1*z^-1
       - 5*ln(x)*[1-x]^-1*z^-1*rln2
       - 3/2*ln(x)*[1-x]^-1*z^-1*sqrtxz1
       - 28*ln(x)*[1-x]^-1
       + 7*ln(x)*[1-x]^-1*rln2
       + 3/2*ln(x)*[1-x]^-1*sqrtxz1
       + 99/4*ln(x)*[1-x]^-1*z
       - 5/2*ln(x)*[1-x]^-1*z*rln2
       + 4/3*ln(x)*[1-x]^-1*z^2
       - 265/12*ln(x)*z^-1
       + 5/2*ln(x)*z^-1*rln2
       + 1/2*ln(x)*z^-1*sqrtxz1
       - ln(x)*[1-x-z]^-1
       - ln(x)*rln2
       - 1/2*ln(x)*sqrtxz1
       + 3*ln(x)*z*[1-x-z]^-2
       - 4*ln(x)*z*[1-x-z]^-1
       - 26*ln(x)*z
       - ln(x)*z*rln2
       - 3*ln(x)*z^2*[1-x-z]^-2
       + 2*ln(x)*z^2*[1-x-z]^-1
       - 8/3*ln(x)*z^2
       + 101/12*ln(x)*x*z^-1
       + 5/2*ln(x)*x*z^-1*rln2
       + ln(x)*x*z^-1*sqrtxz1
       - 3*ln(x)*x*[1-x-z]^-2
       + 6*ln(x)*x*[1-x-z]^-1
       + 9/2*ln(x)*x
       - ln(x)*x*rln2
       - 7/2*ln(x)*x*z
       + 4/3*ln(x)*x*z^2
       + 3*ln(x)*x^2*[1-x-z]^-2
       - 3*ln(x)*x^2*[1-x-z]^-1
       - ln(x)*x^3*[1-x-z]^-2
       + 4*ln(x)*ln(1 + sqrtxz1 - z)*[1-x]^-1*z^-1
       - 4*ln(x)*ln(1 + sqrtxz1 - z)*[1-x]^-1
       + 2*ln(x)*ln(1 + sqrtxz1 - z)*[1-x]^-1*z
       - 2*ln(x)*ln(1 + sqrtxz1 - z)*z^-1
       - 2*ln(x)*ln(1 + sqrtxz1 - z)*x*z^-1
       + ln(x)*ln(1 + sqrtxz1 + z)
       + ln(x)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^-1
       - 3*ln(x)*ln(1 + sqrtxz1 + z)*[1-x]^-1
       + 1/2*ln(x)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z
       - 1/2*ln(x)*ln(1 + sqrtxz1 + z)*z^-1
       + ln(x)*ln(1 + sqrtxz1 + z)*z
       - 1/2*ln(x)*ln(1 + sqrtxz1 + z)*x*z^-1
       + ln(x)*ln(1 + sqrtxz1 + z)*x
       - 16*ln(x)^2
       - 27/2*ln(x)^2*[1-x]^-1*z^-1
       + 17*ln(x)^2*[1-x]^-1
       - 17/2*ln(x)^2*[1-x]^-1*z
       + 9*ln(x)^2*z^-1
       + 12*ln(x)^2*z
       + 11/2*ln(x)^2*x*z^-1
       - 9*ln(x)^2*x
       + 1/2*ln(x)^2*x*z
       + 23*ln(x)*ln([1-x])
       + 18*ln(x)*ln([1-x])*[1-x]^-1*z^-1
       - 21*ln(x)*ln([1-x])*[1-x]^-1
       + 21/2*ln(x)*ln([1-x])*[1-x]^-1*z
       - 23/2*ln(x)*ln([1-x])*z^-1
       - 17*ln(x)*ln([1-x])*z
       - 13/2*ln(x)*ln([1-x])*x*z^-1
       + 13*ln(x)*ln([1-x])*x
       - ln(x)*ln([1-x])*x*z
       + 18*ln(x)*ln(z)
       + 11*ln(x)*ln(z)*[1-x]^-1*z^-1
       - 22*ln(x)*ln(z)*[1-x]^-1
       + 9/2*ln(x)*ln(z)*[1-x]^-1*z
       - 7*ln(x)*ln(z)*z^-1
       - 8*ln(x)*ln(z)*z
       - 4*ln(x)*ln(z)*x*z^-1
       + 14*ln(x)*ln(z)*x
       + 23*ln(x)*ln([1-z])
       + 19*ln(x)*ln([1-z])*[1-x]^-1*z^-1
       - 24*ln(x)*ln([1-z])*[1-x]^-1
       + 12*ln(x)*ln([1-z])*[1-x]^-1*z
       - 25/2*ln(x)*ln([1-z])*z^-1
       - 17*ln(x)*ln([1-z])*z
       - 15/2*ln(x)*ln([1-z])*x*z^-1
       + 13*ln(x)*ln([1-z])*x
       - ln(x)*ln([1-z])*x*z
       - 27/2*ln([1-x])
       - 15*ln([1-x])*[1-x]^-1*z^-1
       + 17*ln([1-x])*[1-x]^-1
       - 11*ln([1-x])*[1-x]^-1*z
       + 125/12*ln([1-x])*z^-1
       + 27/2*ln([1-x])*z
       + 4/3*ln([1-x])*z^2
       - 79/12*ln([1-x])*x*z^-1
       - 5/2*ln([1-x])*x
       + 2*ln([1-x])*x*z
       - 2/3*ln([1-x])*x*z^2
       - 17/2*ln([1-x])^2
       - 9/2*ln([1-x])^2*[1-x]^-1*z^-1
       + 5*ln([1-x])^2*[1-x]^-1
       - 5/2*ln([1-x])^2*[1-x]^-1*z
       + 19/4*ln([1-x])^2*z^-1
       + 6*ln([1-x])^2*z
       + 11/4*ln([1-x])^2*x*z^-1
       - 9/2*ln([1-x])^2*x
       + 1/2*ln([1-x])^2*x*z
       - 12*ln([1-x])*ln(z)
       - 8*ln([1-x])*ln(z)*[1-x]^-1*z^-1
       + 10*ln([1-x])*ln(z)*[1-x]^-1
       - 5*ln([1-x])*ln(z)*[1-x]^-1*z
       + 5*ln([1-x])*ln(z)*z^-1
       + 5*ln([1-x])*ln(z)*z
       + 3*ln([1-x])*ln(z)*x*z^-1
       - 10*ln([1-x])*ln(z)*x
       - 15*ln([1-x])*ln([1-z])
       - 8*ln([1-x])*ln([1-z])*[1-x]^-1*z^-1
       + 10*ln([1-x])*ln([1-z])*[1-x]^-1
       - 5*ln([1-x])*ln([1-z])*[1-x]^-1*z
       + 17/2*ln([1-x])*ln([1-z])*z^-1
       + 11*ln([1-x])*ln([1-z])*z
       + 11/2*ln([1-x])*ln([1-z])*x*z^-1
       - 9*ln([1-x])*ln([1-z])*x
       + ln([1-x])*ln([1-z])*x*z
       - 25*ln(z)
       - 27/2*ln(z)*[1-x]^-1*z^-1
       + ln(z)*[1-x]^-1*z^-1*rln2
       - 3/2*ln(z)*[1-x]^-1*z^-1*sqrtxz1
       + 19*ln(z)*[1-x]^-1
       + 5*ln(z)*[1-x]^-1*rln2
       + 3/2*ln(z)*[1-x]^-1*sqrtxz1
       - 27/2*ln(z)*[1-x]^-1*z
       + 1/2*ln(z)*[1-x]^-1*z*rln2
       + 26/3*ln(z)*z^-1
       - 1/2*ln(z)*z^-1*rln2
       + 1/2*ln(z)*z^-1*sqrtxz1
       - 3*ln(z)*rln2
       - 1/2*ln(z)*sqrtxz1
       + 18*ln(z)*z
       - 3*ln(z)*z*rln2
       + 4/3*ln(z)*z^2
       - 47/6*ln(z)*x*z^-1
       - 1/2*ln(z)*x*z^-1*rln2
       + ln(z)*x*z^-1*sqrtxz1
       - 17/2*ln(z)*x
       - 3*ln(z)*x*rln2
       + 3/2*ln(z)*x*z
       - 2/3*ln(z)*x*z^2
       + ln(z)*ln(1 + sqrtxz1 - z)
       + ln(z)*ln(1 + sqrtxz1 - z)*[1-x]^-1*z^-1
       - 3*ln(z)*ln(1 + sqrtxz1 - z)*[1-x]^-1
       + 1/2*ln(z)*ln(1 + sqrtxz1 - z)*[1-x]^-1*z
       - 1/2*ln(z)*ln(1 + sqrtxz1 - z)*z^-1
       + ln(z)*ln(1 + sqrtxz1 - z)*z
       - 1/2*ln(z)*ln(1 + sqrtxz1 - z)*x*z^-1
       + ln(z)*ln(1 + sqrtxz1 - z)*x
       + 2*ln(z)*ln(1 + sqrtxz1 + z)
       - 2*ln(z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^-1
       - 2*ln(z)*ln(1 + sqrtxz1 + z)*[1-x]^-1
       - ln(z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z
       + ln(z)*ln(1 + sqrtxz1 + z)*z^-1
       + 2*ln(z)*ln(1 + sqrtxz1 + z)*z
       + ln(z)*ln(1 + sqrtxz1 + z)*x*z^-1
       + 2*ln(z)*ln(1 + sqrtxz1 + z)*x
       + 2*ln(z)*ln(1 + z)*[1-x]^-1*z^-1
       + 2*ln(z)*ln(1 + z)*[1-x]^-1
       + ln(z)*ln(1 + z)*[1-x]^-1*z
       - 11/2*ln(z)^2
       - 2*ln(z)^2*[1-x]^-1*z^-1
       + 9/2*ln(z)^2*[1-x]^-1
       - 7/4*ln(z)^2*[1-x]^-1*z
       - 1/4*ln(z)^2*z^-1
       - 3/2*ln(z)^2*z
       - 3/4*ln(z)^2*x*z^-1
       - 13/2*ln(z)^2*x
       - 1/2*ln(z)^2*x*z
       - 11*ln(z)*ln([1-z])
       - 8*ln(z)*ln([1-z])*[1-x]^-1*z^-1
       + 12*ln(z)*ln([1-z])*[1-x]^-1
       - 6*ln(z)*ln([1-z])*[1-x]^-1*z
       + 11/2*ln(z)*ln([1-z])*z^-1
       + 10*ln(z)*ln([1-z])*z
       + 9/2*ln(z)*ln([1-z])*x*z^-1
       - 9*ln(z)*ln([1-z])*x
       - 33/2*ln([1-z])
       - 15*ln([1-z])*[1-x]^-1*z^-1
       + 19*ln([1-z])*[1-x]^-1
       - 12*ln([1-z])*[1-x]^-1*z
       + 125/12*ln([1-z])*z^-1
       + ln([1-z])*[1-x-z]^-1
       - 3*ln([1-z])*z*[1-x-z]^-2
       + 4*ln([1-z])*z*[1-x-z]^-1
       + 10*ln([1-z])*z
       + 3*ln([1-z])*z^2*[1-x-z]^-2
       - 2*ln([1-z])*z^2*[1-x-z]^-1
       + 4/3*ln([1-z])*z^2
       - 79/12*ln([1-z])*x*z^-1
       + 3*ln([1-z])*x*[1-x-z]^-2
       - 6*ln([1-z])*x*[1-x-z]^-1
       + ln([1-z])*x
       + 7/2*ln([1-z])*x*z
       - 2/3*ln([1-z])*x*z^2
       - 3*ln([1-z])*x^2*[1-x-z]^-2
       + 3*ln([1-z])*x^2*[1-x-z]^-1
       + ln([1-z])*x^3*[1-x-z]^-2
       - 17/2*ln([1-z])^2
       - 9/2*ln([1-z])^2*[1-x]^-1*z^-1
       + 13/2*ln([1-z])^2*[1-x]^-1
       - 13/4*ln([1-z])^2*[1-x]^-1*z
       + 19/4*ln([1-z])^2*z^-1
       + 6*ln([1-z])^2*z
       + 11/4*ln([1-z])^2*x*z^-1
       - 9/2*ln([1-z])^2*x
       + 1/2*ln([1-z])^2*x*z
       - Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)
       + 3*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*z^-1
       - Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1
       + 3/2*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*z
       - 3/2*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z^-1
       - Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z
       - 3/2*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x*z^-1
       - Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x
       + Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)
       - 3*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*z^-1
       + Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1
       - 3/2*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*z
       + 3/2*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z^-1
       + Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z
       + 3/2*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x*z^-1
       + Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x
       + Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)
       + Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*[1-x]^-1*z^-1
       - 3*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*[1-x]^-1
       + 1/2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*[1-x]^-1*z
       - 1/2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*z^-1
       + Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*z
       - 1/2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*x*z^-1
       + Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*x
       - Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)
       - Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*[1-x]^-1*z^-1
       + 3*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*[1-x]^-1
       - 1/2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*[1-x]^-1*z
       + 1/2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*z^-1
       - Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*z
       + 1/2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*x*z^-1
       - Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*x
       + 2*Li2(1 - x*z^-1)
       + 2*Li2(1 - x*z^-1)*[1-x]^-1*z^-1
       - 3*Li2(1 - x*z^-1)*[1-x]^-1
       + 3/2*Li2(1 - x*z^-1)*[1-x]^-1*z
       - Li2(1 - x*z^-1)*z^-1
       - 2*Li2(1 - x*z^-1)*z
       - Li2(1 - x*z^-1)*x*z^-1
       + 2*Li2(1 - x*z^-1)*x
       + 2*Li2( - z)*[1-x]^-1*z^-1
       + 2*Li2( - z)*[1-x]^-1
       + Li2( - z)*[1-x]^-1*z
       - Li2(x)*[1-x]^-1*z^-1
       + Li2(x)*[1-x]^-1
       - 1/2*Li2(x)*[1-x]^-1*z
       + Li2(x)*z^-1
       + Li2(x)*x*z^-1
       + Li2(z)
       - 2*Li2(z)*[1-x]^-1*z^-1
       + 3*Li2(z)*[1-x]^-1
       - 3/2*Li2(z)*[1-x]^-1*z
       + 1/2*Li2(z)*z^-1
       + 5*Li2(z)*z
       + 3/2*Li2(z)*x*z^-1
       + Li2(z)*x );

id DC2G2G      = (  + NC^-1 * ( 7/4
       - 8*z^-1*[1+z]^-1*rln2^2
       - 16*z^-1*[1+z]^-1*sqrtxz1*rln2
       - 12*z^-1*[1+z]^-1*sqrtxz1*rln2^2
       + 13/8*z^-1
       + 8*z^-1*rln2^2
       + 16*z^-1*sqrtxz1*rln2
       + 12*z^-1*sqrtxz1*rln2^2
       - 48*[1+z]^-1*sqrtxz1^-1*rln2
       - 36*[1+z]^-1*sqrtxz1^-1*rln2^2
       + 32*sqrtxz1^-1*rln2
       + 24*sqrtxz1^-1*rln2^2
       - 8*rln2^2
       + 16*z*[1+z]^-1*sqrtxz1^-1*rln2
       + 12*z*[1+z]^-1*sqrtxz1^-1*rln2^2
       - 16*z*sqrtxz1^-1*rln2
       - 12*z*sqrtxz1^-1*rln2^2
       - 43/8*z
       + 4*z*rln2^2
       + 16*x*z^-1*[1+z]^-1*rln2^2
       + 32*x*z^-1*[1+z]^-1*sqrtxz1*rln2
       + 24*x*z^-1*[1+z]^-1*sqrtxz1*rln2^2
       + 3/2*x*z^-1
       - 16*x*z^-1*rln2^2
       - 32*x*z^-1*sqrtxz1*rln2
       - 24*x*z^-1*sqrtxz1*rln2^2
       + 160*x*[1+z]^-1*sqrtxz1^-1*rln2
       + 120*x*[1+z]^-1*sqrtxz1^-1*rln2^2
       - 128*x*sqrtxz1^-1*rln2
       - 96*x*sqrtxz1^-1*rln2^2
       - 5*x
       + 16*x*rln2^2
       - 32*x*z*[1+z]^-1*sqrtxz1^-1*rln2
       - 24*x*z*[1+z]^-1*sqrtxz1^-1*rln2^2
       + 32*x*z*sqrtxz1^-1*rln2
       + 24*x*z*sqrtxz1^-1*rln2^2
       + 11/2*x*z
       - 8*x*z*rln2^2
       - 8*x^2*z^-1*[1+z]^-1*rln2^2
       - 16*x^2*z^-1*[1+z]^-1*sqrtxz1*rln2
       - 12*x^2*z^-1*[1+z]^-1*sqrtxz1*rln2^2
       + 8*x^2*z^-1*rln2^2
       + 16*x^2*z^-1*sqrtxz1*rln2
       + 12*x^2*z^-1*sqrtxz1*rln2^2
       - 176*x^2*[1+z]^-1*sqrtxz1^-1*rln2
       - 132*x^2*[1+z]^-1*sqrtxz1^-1*rln2^2
       + 160*x^2*sqrtxz1^-1*rln2
       + 120*x^2*sqrtxz1^-1*rln2^2
       - 8*x^2*rln2^2
       + 16*x^2*z*[1+z]^-1*sqrtxz1^-1*rln2
       + 12*x^2*z*[1+z]^-1*sqrtxz1^-1*rln2^2
       - 16*x^2*z*sqrtxz1^-1*rln2
       - 12*x^2*z*sqrtxz1^-1*rln2^2
       + 64*x^3*[1+z]^-1*sqrtxz1^-1*rln2
       + 48*x^3*[1+z]^-1*sqrtxz1^-1*rln2^2
       - 64*x^3*sqrtxz1^-1*rln2
       - 48*x^3*sqrtxz1^-1*rln2^2
       - 2/3*pi^2*x^-2*[1+x]^-1*z^-1*[1-z]^-1
       + 2/3*pi^2*x^-2*[1+x]^-1*z^-1
       + 1/3*pi^2*x^-2*[1+x]^-1
       + 1/3*pi^2*x^-2*[1+x]^-1*z
       + 2/3*pi^2*x^-2*z^-1*[1-z]^-1
       - 2/3*pi^2*x^-2*z^-1
       - 1/3*pi^2*x^-2
       - 1/3*pi^2*x^-2*z
       - 2*pi^2*x^-1*[1+x]^-1*z^-1*[1-z]^-1
       + 2*pi^2*x^-1*[1+x]^-1*z^-1
       + pi^2*x^-1*[1+x]^-1
       + pi^2*x^-1*[1+x]^-1*z
       + 4/3*pi^2*x^-1*z^-1*[1-z]^-1
       - 4/3*pi^2*x^-1*z^-1
       - 2/3*pi^2*x^-1
       - 2/3*pi^2*x^-1*z
       - 2*pi^2*[1+x]^-1*z^-1*[1-z]^-1
       + 2*pi^2*[1+x]^-1*z^-1
       + 4/3*pi^2*[1+x]^-1
       + 2/3*pi^2*[1+x]^-1*z
       - 2/3*pi^2*z^-1*[1+z]^-1*sqrtxz1
       - 7/12*pi^2*z^-1
       + 2/3*pi^2*z^-1*sqrtxz1
       + 2/3*pi^2*[1-z]^-1
       - 2*pi^2*[1+z]^-1*sqrtxz1^-1
       + 4/3*pi^2*sqrtxz1^-1
       - 1/6*pi^2
       + 2/3*pi^2*z*[1+z]^-1*sqrtxz1^-1
       - 2/3*pi^2*z*sqrtxz1^-1
       - 1/3*pi^2*z
       - 2/3*pi^2*x*[1+x]^-1*z^-1*[1-z]^-1
       + 2/3*pi^2*x*[1+x]^-1*z^-1
       + 2/3*pi^2*x*[1+x]^-1
       - 4/3*pi^2*x*z^-1*[1-z]^-1
       + 4/3*pi^2*x*z^-1*[1+z]^-1*sqrtxz1
       + 5/2*pi^2*x*z^-1
       - 4/3*pi^2*x*z^-1*sqrtxz1
       + 4/3*pi^2*x*[1-z]^-1
       + 20/3*pi^2*x*[1+z]^-1*sqrtxz1^-1
       - 16/3*pi^2*x*sqrtxz1^-1
       - pi^2*x
       - 4/3*pi^2*x*z*[1+z]^-1*sqrtxz1^-1
       + 4/3*pi^2*x*z*sqrtxz1^-1
       + pi^2*x*z
       - 2/3*pi^2*x^2*z^-1*[1-z]^-1
       - 2/3*pi^2*x^2*z^-1*[1+z]^-1*sqrtxz1
       + 2/3*pi^2*x^2*z^-1
       + 2/3*pi^2*x^2*z^-1*sqrtxz1
       + 2/3*pi^2*x^2*[1-z]^-1
       - 22/3*pi^2*x^2*[1+z]^-1*sqrtxz1^-1
       + 20/3*pi^2*x^2*sqrtxz1^-1
       + 2/3*pi^2*x^2
       + 2/3*pi^2*x^2*z*[1+z]^-1*sqrtxz1^-1
       - 2/3*pi^2*x^2*z*sqrtxz1^-1
       + 8/3*pi^2*x^3*[1+z]^-1*sqrtxz1^-1
       - 8/3*pi^2*x^3*sqrtxz1^-1
       + 8*ln(1 + sqrtxz1 - z)*z^-1*[1+z]^-1*rln2
       + 16*ln(1 + sqrtxz1 - z)*z^-1*[1+z]^-1*sqrtxz1
       + 16*ln(1 + sqrtxz1 - z)*z^-1*[1+z]^-1*sqrtxz1*rln2
       - 10*ln(1 + sqrtxz1 - z)*z^-1*rln2
       - 16*ln(1 + sqrtxz1 - z)*z^-1*sqrtxz1
       - 16*ln(1 + sqrtxz1 - z)*z^-1*sqrtxz1*rln2
       + 48*ln(1 + sqrtxz1 - z)*[1+z]^-1*sqrtxz1^-1
       + 48*ln(1 + sqrtxz1 - z)*[1+z]^-1*sqrtxz1^-1*rln2
       - 32*ln(1 + sqrtxz1 - z)*sqrtxz1^-1
       - 32*ln(1 + sqrtxz1 - z)*sqrtxz1^-1*rln2
       + 12*ln(1 + sqrtxz1 - z)*rln2
       - 16*ln(1 + sqrtxz1 - z)*z*[1+z]^-1*sqrtxz1^-1
       - 16*ln(1 + sqrtxz1 - z)*z*[1+z]^-1*sqrtxz1^-1*rln2
       + 16*ln(1 + sqrtxz1 - z)*z*sqrtxz1^-1
       + 16*ln(1 + sqrtxz1 - z)*z*sqrtxz1^-1*rln2
       - 6*ln(1 + sqrtxz1 - z)*z*rln2
       - 16*ln(1 + sqrtxz1 - z)*x*z^-1*[1+z]^-1*rln2
       - 32*ln(1 + sqrtxz1 - z)*x*z^-1*[1+z]^-1*sqrtxz1
       - 32*ln(1 + sqrtxz1 - z)*x*z^-1*[1+z]^-1*sqrtxz1*rln2
       + 20*ln(1 + sqrtxz1 - z)*x*z^-1*rln2
       + 32*ln(1 + sqrtxz1 - z)*x*z^-1*sqrtxz1
       + 32*ln(1 + sqrtxz1 - z)*x*z^-1*sqrtxz1*rln2
       - 160*ln(1 + sqrtxz1 - z)*x*[1+z]^-1*sqrtxz1^-1
       - 160*ln(1 + sqrtxz1 - z)*x*[1+z]^-1*sqrtxz1^-1*rln2
       + 128*ln(1 + sqrtxz1 - z)*x*sqrtxz1^-1
       + 128*ln(1 + sqrtxz1 - z)*x*sqrtxz1^-1*rln2
       - 24*ln(1 + sqrtxz1 - z)*x*rln2
       + 32*ln(1 + sqrtxz1 - z)*x*z*[1+z]^-1*sqrtxz1^-1
       + 32*ln(1 + sqrtxz1 - z)*x*z*[1+z]^-1*sqrtxz1^-1*rln2
       - 32*ln(1 + sqrtxz1 - z)*x*z*sqrtxz1^-1
       - 32*ln(1 + sqrtxz1 - z)*x*z*sqrtxz1^-1*rln2
       + 12*ln(1 + sqrtxz1 - z)*x*z*rln2
       + 8*ln(1 + sqrtxz1 - z)*x^2*z^-1*[1+z]^-1*rln2
       + 16*ln(1 + sqrtxz1 - z)*x^2*z^-1*[1+z]^-1*sqrtxz1
       + 16*ln(1 + sqrtxz1 - z)*x^2*z^-1*[1+z]^-1*sqrtxz1*rln2
       - 8*ln(1 + sqrtxz1 - z)*x^2*z^-1*rln2
       - 16*ln(1 + sqrtxz1 - z)*x^2*z^-1*sqrtxz1
       - 16*ln(1 + sqrtxz1 - z)*x^2*z^-1*sqrtxz1*rln2
       + 176*ln(1 + sqrtxz1 - z)*x^2*[1+z]^-1*sqrtxz1^-1
       + 176*ln(1 + sqrtxz1 - z)*x^2*[1+z]^-1*sqrtxz1^-1*rln2
       - 160*ln(1 + sqrtxz1 - z)*x^2*sqrtxz1^-1
       - 160*ln(1 + sqrtxz1 - z)*x^2*sqrtxz1^-1*rln2
       + 12*ln(1 + sqrtxz1 - z)*x^2*rln2
       - 16*ln(1 + sqrtxz1 - z)*x^2*z*[1+z]^-1*sqrtxz1^-1
       - 16*ln(1 + sqrtxz1 - z)*x^2*z*[1+z]^-1*sqrtxz1^-1*rln2
       + 16*ln(1 + sqrtxz1 - z)*x^2*z*sqrtxz1^-1
       + 16*ln(1 + sqrtxz1 - z)*x^2*z*sqrtxz1^-1*rln2
       - 64*ln(1 + sqrtxz1 - z)*x^3*[1+z]^-1*sqrtxz1^-1
       - 64*ln(1 + sqrtxz1 - z)*x^3*[1+z]^-1*sqrtxz1^-1*rln2
       + 64*ln(1 + sqrtxz1 - z)*x^3*sqrtxz1^-1
       + 64*ln(1 + sqrtxz1 - z)*x^3*sqrtxz1^-1*rln2
       - 4*ln(1 + sqrtxz1 - z)^2
       - 4*ln(1 + sqrtxz1 - z)^2*z^-1*[1+z]^-1*sqrtxz1
       + 2*ln(1 + sqrtxz1 - z)^2*z^-1
       + 4*ln(1 + sqrtxz1 - z)^2*z^-1*sqrtxz1
       - 12*ln(1 + sqrtxz1 - z)^2*[1+z]^-1*sqrtxz1^-1
       + 8*ln(1 + sqrtxz1 - z)^2*sqrtxz1^-1
       + 4*ln(1 + sqrtxz1 - z)^2*z*[1+z]^-1*sqrtxz1^-1
       - 4*ln(1 + sqrtxz1 - z)^2*z*sqrtxz1^-1
       + 2*ln(1 + sqrtxz1 - z)^2*z
       + 8*ln(1 + sqrtxz1 - z)^2*x*z^-1*[1+z]^-1*sqrtxz1
       - 4*ln(1 + sqrtxz1 - z)^2*x*z^-1
       - 8*ln(1 + sqrtxz1 - z)^2*x*z^-1*sqrtxz1
       + 40*ln(1 + sqrtxz1 - z)^2*x*[1+z]^-1*sqrtxz1^-1
       - 32*ln(1 + sqrtxz1 - z)^2*x*sqrtxz1^-1
       + 8*ln(1 + sqrtxz1 - z)^2*x
       - 8*ln(1 + sqrtxz1 - z)^2*x*z*[1+z]^-1*sqrtxz1^-1
       + 8*ln(1 + sqrtxz1 - z)^2*x*z*sqrtxz1^-1
       - 4*ln(1 + sqrtxz1 - z)^2*x*z
       - 4*ln(1 + sqrtxz1 - z)^2*x^2*z^-1*[1+z]^-1*sqrtxz1
       + 4*ln(1 + sqrtxz1 - z)^2*x^2*z^-1*sqrtxz1
       - 44*ln(1 + sqrtxz1 - z)^2*x^2*[1+z]^-1*sqrtxz1^-1
       + 40*ln(1 + sqrtxz1 - z)^2*x^2*sqrtxz1^-1
       - 4*ln(1 + sqrtxz1 - z)^2*x^2
       + 4*ln(1 + sqrtxz1 - z)^2*x^2*z*[1+z]^-1*sqrtxz1^-1
       - 4*ln(1 + sqrtxz1 - z)^2*x^2*z*sqrtxz1^-1
       + 16*ln(1 + sqrtxz1 - z)^2*x^3*[1+z]^-1*sqrtxz1^-1
       - 16*ln(1 + sqrtxz1 - z)^2*x^3*sqrtxz1^-1
       - 4*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)
       - 8*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*z^-1*[1+z]^-1
       - 8*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*z^-1*[1+z]^-1*sqrtxz1
       + 6*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*z^-1
       + 8*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*z^-1*sqrtxz1
       - 24*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*[1+z]^-1*sqrtxz1^-1
       + 16*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*sqrtxz1^-1
       + 8*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*z*[1+z]^-1*sqrtxz1^-1
       - 8*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*z*sqrtxz1^-1
       + 2*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*z
       + 16*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*x*z^-1*[1+z]^-1
       + 16*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*x*z^-1*[1+z]^-1*sqrtxz1
       - 12*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*x*z^-1
       - 16*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*x*z^-1*sqrtxz1
       + 80*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*x*[1+z]^-1*sqrtxz1^-1
       - 64*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*x*sqrtxz1^-1
       + 8*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*x
       - 16*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*x*z*[1+z]^-1*sqrtxz1^-1
       + 16*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*x*z*sqrtxz1^-1
       - 4*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*x*z
       - 8*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*x^2*z^-1*[1+z]^-1
       - 8*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*x^2*z^-1*[1+z]^-1*sqrtxz1
       + 8*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*x^2*z^-1
       + 8*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*x^2*z^-1*sqrtxz1
       - 88*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*x^2*[1+z]^-1*sqrtxz1^-1
       + 80*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*x^2*sqrtxz1^-1
       - 4*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*x^2
       + 8*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*x^2*z*[1+z]^-1*sqrtxz1^-1
       - 8*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*x^2*z*sqrtxz1^-1
       + 32*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*x^3*[1+z]^-1*sqrtxz1^-1
       - 32*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*x^3*sqrtxz1^-1
       + 8*ln(1 + sqrtxz1 + z)*z^-1*[1+z]^-1*rln2
       + 8*ln(1 + sqrtxz1 + z)*z^-1*[1+z]^-1*sqrtxz1*rln2
       - 6*ln(1 + sqrtxz1 + z)*z^-1*rln2
       - 8*ln(1 + sqrtxz1 + z)*z^-1*sqrtxz1*rln2
       + 24*ln(1 + sqrtxz1 + z)*[1+z]^-1*sqrtxz1^-1*rln2
       - 16*ln(1 + sqrtxz1 + z)*sqrtxz1^-1*rln2
       + 4*ln(1 + sqrtxz1 + z)*rln2
       - 8*ln(1 + sqrtxz1 + z)*z*[1+z]^-1*sqrtxz1^-1*rln2
       + 8*ln(1 + sqrtxz1 + z)*z*sqrtxz1^-1*rln2
       - 2*ln(1 + sqrtxz1 + z)*z*rln2
       - 16*ln(1 + sqrtxz1 + z)*x*z^-1*[1+z]^-1*rln2
       - 16*ln(1 + sqrtxz1 + z)*x*z^-1*[1+z]^-1*sqrtxz1*rln2
       + 12*ln(1 + sqrtxz1 + z)*x*z^-1*rln2
       + 16*ln(1 + sqrtxz1 + z)*x*z^-1*sqrtxz1*rln2
       - 80*ln(1 + sqrtxz1 + z)*x*[1+z]^-1*sqrtxz1^-1*rln2
       + 64*ln(1 + sqrtxz1 + z)*x*sqrtxz1^-1*rln2
       - 8*ln(1 + sqrtxz1 + z)*x*rln2
       + 16*ln(1 + sqrtxz1 + z)*x*z*[1+z]^-1*sqrtxz1^-1*rln2
       - 16*ln(1 + sqrtxz1 + z)*x*z*sqrtxz1^-1*rln2
       + 4*ln(1 + sqrtxz1 + z)*x*z*rln2
       + 8*ln(1 + sqrtxz1 + z)*x^2*z^-1*[1+z]^-1*rln2
       + 8*ln(1 + sqrtxz1 + z)*x^2*z^-1*[1+z]^-1*sqrtxz1*rln2
       - 8*ln(1 + sqrtxz1 + z)*x^2*z^-1*rln2
       - 8*ln(1 + sqrtxz1 + z)*x^2*z^-1*sqrtxz1*rln2
       + 88*ln(1 + sqrtxz1 + z)*x^2*[1+z]^-1*sqrtxz1^-1*rln2
       - 80*ln(1 + sqrtxz1 + z)*x^2*sqrtxz1^-1*rln2
       + 4*ln(1 + sqrtxz1 + z)*x^2*rln2
       - 8*ln(1 + sqrtxz1 + z)*x^2*z*[1+z]^-1*sqrtxz1^-1*rln2
       + 8*ln(1 + sqrtxz1 + z)*x^2*z*sqrtxz1^-1*rln2
       - 32*ln(1 + sqrtxz1 + z)*x^3*[1+z]^-1*sqrtxz1^-1*rln2
       + 32*ln(1 + sqrtxz1 + z)*x^3*sqrtxz1^-1*rln2
       - ln(1 - 2*z + z^2 + 4*x*z)^2*z^-1*[1+z]^-1*sqrtxz1
       + ln(1 - 2*z + z^2 + 4*x*z)^2*z^-1*sqrtxz1
       - 3*ln(1 - 2*z + z^2 + 4*x*z)^2*[1+z]^-1*sqrtxz1^-1
       + 2*ln(1 - 2*z + z^2 + 4*x*z)^2*sqrtxz1^-1
       + ln(1 - 2*z + z^2 + 4*x*z)^2*z*[1+z]^-1*sqrtxz1^-1
       - ln(1 - 2*z + z^2 + 4*x*z)^2*z*sqrtxz1^-1
       + 2*ln(1 - 2*z + z^2 + 4*x*z)^2*x*z^-1*[1+z]^-1*sqrtxz1
       - 2*ln(1 - 2*z + z^2 + 4*x*z)^2*x*z^-1*sqrtxz1
       + 10*ln(1 - 2*z + z^2 + 4*x*z)^2*x*[1+z]^-1*sqrtxz1^-1
       - 8*ln(1 - 2*z + z^2 + 4*x*z)^2*x*sqrtxz1^-1
       - 2*ln(1 - 2*z + z^2 + 4*x*z)^2*x*z*[1+z]^-1*sqrtxz1^-1
       + 2*ln(1 - 2*z + z^2 + 4*x*z)^2*x*z*sqrtxz1^-1
       - ln(1 - 2*z + z^2 + 4*x*z)^2*x^2*z^-1*[1+z]^-1*sqrtxz1
       + ln(1 - 2*z + z^2 + 4*x*z)^2*x^2*z^-1*sqrtxz1
       - 11*ln(1 - 2*z + z^2 + 4*x*z)^2*x^2*[1+z]^-1*sqrtxz1^-1
       + 10*ln(1 - 2*z + z^2 + 4*x*z)^2*x^2*sqrtxz1^-1
       + ln(1 - 2*z + z^2 + 4*x*z)^2*x^2*z*[1+z]^-1*sqrtxz1^-1
       - ln(1 - 2*z + z^2 + 4*x*z)^2*x^2*z*sqrtxz1^-1
       + 4*ln(1 - 2*z + z^2 + 4*x*z)^2*x^3*[1+z]^-1*sqrtxz1^-1
       - 4*ln(1 - 2*z + z^2 + 4*x*z)^2*x^3*sqrtxz1^-1
       - 8*ln(z*sqrtxz3)*ArcTan(z*sqrtxz3)*z*sqrtxz3
       - 11*ln(x)
       - 8*ln(x)*x^-2*[1+x]^-1*z^-1*[1-z]^-1
       + 8*ln(x)*x^-2*[1+x]^-1*z^-1
       + 4*ln(x)*x^-2*[1+x]^-1
       + 4*ln(x)*x^-2*[1+x]^-1*z
       + 8*ln(x)*x^-2*z^-1*[1-z]^-1
       - 8*ln(x)*x^-2*z^-1
       - 4*ln(x)*x^-2
       - 4*ln(x)*x^-2*z
       - 24*ln(x)*x^-1*[1+x]^-1*z^-1*[1-z]^-1
       + 24*ln(x)*x^-1*[1+x]^-1*z^-1
       + 12*ln(x)*x^-1*[1+x]^-1
       + 12*ln(x)*x^-1*[1+x]^-1*z
       + 16*ln(x)*x^-1*z^-1*[1-z]^-1
       - 16*ln(x)*x^-1*z^-1
       - 8*ln(x)*x^-1
       - 8*ln(x)*x^-1*z
       - 24*ln(x)*[1+x]^-1*z^-1*[1-z]^-1
       + 24*ln(x)*[1+x]^-1*z^-1
       + 16*ln(x)*[1+x]^-1
       + 8*ln(x)*[1+x]^-1*z
       - 4*ln(x)*z^-1*[1+z]^-1*rln2
       - 8*ln(x)*z^-1*[1+z]^-1*sqrtxz1
       + 11/2*ln(x)*z^-1
       + 6*ln(x)*z^-1*rln2
       + 8*ln(x)*z^-1*sqrtxz1
       + 8*ln(x)*[1-z]^-1
       - 24*ln(x)*[1+z]^-1*sqrtxz1^-1
       + 16*ln(x)*sqrtxz1^-1
       - 8*ln(x)*rln2
       + 8*ln(x)*z*[1+z]^-1*sqrtxz1^-1
       - 8*ln(x)*z*sqrtxz1^-1
       - 7/2*ln(x)*z
       + 4*ln(x)*z*rln2
       - 8*ln(x)*x*[1+x]^-1*z^-1*[1-z]^-1
       + 8*ln(x)*x*[1+x]^-1*z^-1
       + 8*ln(x)*x*[1+x]^-1
       - 16*ln(x)*x*z^-1*[1-z]^-1
       + 8*ln(x)*x*z^-1*[1+z]^-1*rln2
       + 16*ln(x)*x*z^-1*[1+z]^-1*sqrtxz1
       + 8*ln(x)*x*z^-1
       - 12*ln(x)*x*z^-1*rln2
       - 16*ln(x)*x*z^-1*sqrtxz1
       + 16*ln(x)*x*[1-z]^-1
       + 80*ln(x)*x*[1+z]^-1*sqrtxz1^-1
       - 64*ln(x)*x*sqrtxz1^-1
       + 8*ln(x)*x
       + 16*ln(x)*x*rln2
       - 16*ln(x)*x*z*[1+z]^-1*sqrtxz1^-1
       + 16*ln(x)*x*z*sqrtxz1^-1
       - 8*ln(x)*x*z*rln2
       - 8*ln(x)*x^2*z^-1*[1-z]^-1
       - 4*ln(x)*x^2*z^-1*[1+z]^-1*rln2
       - 8*ln(x)*x^2*z^-1*[1+z]^-1*sqrtxz1
       + 8*ln(x)*x^2*z^-1
       + 4*ln(x)*x^2*z^-1*rln2
       + 8*ln(x)*x^2*z^-1*sqrtxz1
       + 8*ln(x)*x^2*[1-z]^-1
       - 88*ln(x)*x^2*[1+z]^-1*sqrtxz1^-1
       + 80*ln(x)*x^2*sqrtxz1^-1
       - 8*ln(x)*x^2*rln2
       + 8*ln(x)*x^2*z*[1+z]^-1*sqrtxz1^-1
       - 8*ln(x)*x^2*z*sqrtxz1^-1
       + 32*ln(x)*x^3*[1+z]^-1*sqrtxz1^-1
       - 32*ln(x)*x^3*sqrtxz1^-1
       + 4*ln(x)*ln(1 + sqrtxz1 - z)
       - 4*ln(x)*ln(1 + sqrtxz1 - z)*z^-1*[1+z]^-1*sqrtxz1
       - 2*ln(x)*ln(1 + sqrtxz1 - z)*z^-1
       + 4*ln(x)*ln(1 + sqrtxz1 - z)*z^-1*sqrtxz1
       - 12*ln(x)*ln(1 + sqrtxz1 - z)*[1+z]^-1*sqrtxz1^-1
       + 8*ln(x)*ln(1 + sqrtxz1 - z)*sqrtxz1^-1
       + 4*ln(x)*ln(1 + sqrtxz1 - z)*z*[1+z]^-1*sqrtxz1^-1
       - 4*ln(x)*ln(1 + sqrtxz1 - z)*z*sqrtxz1^-1
       - 2*ln(x)*ln(1 + sqrtxz1 - z)*z
       + 8*ln(x)*ln(1 + sqrtxz1 - z)*x*z^-1*[1+z]^-1*sqrtxz1
       + 4*ln(x)*ln(1 + sqrtxz1 - z)*x*z^-1
       - 8*ln(x)*ln(1 + sqrtxz1 - z)*x*z^-1*sqrtxz1
       + 40*ln(x)*ln(1 + sqrtxz1 - z)*x*[1+z]^-1*sqrtxz1^-1
       - 32*ln(x)*ln(1 + sqrtxz1 - z)*x*sqrtxz1^-1
       - 8*ln(x)*ln(1 + sqrtxz1 - z)*x
       - 8*ln(x)*ln(1 + sqrtxz1 - z)*x*z*[1+z]^-1*sqrtxz1^-1
       + 8*ln(x)*ln(1 + sqrtxz1 - z)*x*z*sqrtxz1^-1
       + 4*ln(x)*ln(1 + sqrtxz1 - z)*x*z
       - 4*ln(x)*ln(1 + sqrtxz1 - z)*x^2*z^-1*[1+z]^-1*sqrtxz1
       + 4*ln(x)*ln(1 + sqrtxz1 - z)*x^2*z^-1*sqrtxz1
       - 44*ln(x)*ln(1 + sqrtxz1 - z)*x^2*[1+z]^-1*sqrtxz1^-1
       + 40*ln(x)*ln(1 + sqrtxz1 - z)*x^2*sqrtxz1^-1
       + 4*ln(x)*ln(1 + sqrtxz1 - z)*x^2
       + 4*ln(x)*ln(1 + sqrtxz1 - z)*x^2*z*[1+z]^-1*sqrtxz1^-1
       - 4*ln(x)*ln(1 + sqrtxz1 - z)*x^2*z*sqrtxz1^-1
       + 16*ln(x)*ln(1 + sqrtxz1 - z)*x^3*[1+z]^-1*sqrtxz1^-1
       - 16*ln(x)*ln(1 + sqrtxz1 - z)*x^3*sqrtxz1^-1
       + 4*ln(x)*ln(1 + sqrtxz1 + z)
       + 4*ln(x)*ln(1 + sqrtxz1 + z)*z^-1*[1+z]^-1
       + 4*ln(x)*ln(1 + sqrtxz1 + z)*z^-1*[1+z]^-1*sqrtxz1
       - 4*ln(x)*ln(1 + sqrtxz1 + z)*z^-1
       - 4*ln(x)*ln(1 + sqrtxz1 + z)*z^-1*sqrtxz1
       + 12*ln(x)*ln(1 + sqrtxz1 + z)*[1+z]^-1*sqrtxz1^-1
       - 8*ln(x)*ln(1 + sqrtxz1 + z)*sqrtxz1^-1
       - 4*ln(x)*ln(1 + sqrtxz1 + z)*z*[1+z]^-1*sqrtxz1^-1
       + 4*ln(x)*ln(1 + sqrtxz1 + z)*z*sqrtxz1^-1
       - 2*ln(x)*ln(1 + sqrtxz1 + z)*z
       - 8*ln(x)*ln(1 + sqrtxz1 + z)*x*z^-1*[1+z]^-1
       - 8*ln(x)*ln(1 + sqrtxz1 + z)*x*z^-1*[1+z]^-1*sqrtxz1
       + 8*ln(x)*ln(1 + sqrtxz1 + z)*x*z^-1
       + 8*ln(x)*ln(1 + sqrtxz1 + z)*x*z^-1*sqrtxz1
       - 40*ln(x)*ln(1 + sqrtxz1 + z)*x*[1+z]^-1*sqrtxz1^-1
       + 32*ln(x)*ln(1 + sqrtxz1 + z)*x*sqrtxz1^-1
       - 8*ln(x)*ln(1 + sqrtxz1 + z)*x
       + 8*ln(x)*ln(1 + sqrtxz1 + z)*x*z*[1+z]^-1*sqrtxz1^-1
       - 8*ln(x)*ln(1 + sqrtxz1 + z)*x*z*sqrtxz1^-1
       + 4*ln(x)*ln(1 + sqrtxz1 + z)*x*z
       + 4*ln(x)*ln(1 + sqrtxz1 + z)*x^2*z^-1*[1+z]^-1
       + 4*ln(x)*ln(1 + sqrtxz1 + z)*x^2*z^-1*[1+z]^-1*sqrtxz1
       - 4*ln(x)*ln(1 + sqrtxz1 + z)*x^2*z^-1
       - 4*ln(x)*ln(1 + sqrtxz1 + z)*x^2*z^-1*sqrtxz1
       + 44*ln(x)*ln(1 + sqrtxz1 + z)*x^2*[1+z]^-1*sqrtxz1^-1
       - 40*ln(x)*ln(1 + sqrtxz1 + z)*x^2*sqrtxz1^-1
       + 4*ln(x)*ln(1 + sqrtxz1 + z)*x^2
       - 4*ln(x)*ln(1 + sqrtxz1 + z)*x^2*z*[1+z]^-1*sqrtxz1^-1
       + 4*ln(x)*ln(1 + sqrtxz1 + z)*x^2*z*sqrtxz1^-1
       - 16*ln(x)*ln(1 + sqrtxz1 + z)*x^3*[1+z]^-1*sqrtxz1^-1
       + 16*ln(x)*ln(1 + sqrtxz1 + z)*x^3*sqrtxz1^-1
       + ln(x)*ln(1 + x*z^-1)
       + 2*ln(x)*ln(1 + x*z^-1)*x^-2*[1+x]^-1*z^-1*[1-z]^-1
       - 2*ln(x)*ln(1 + x*z^-1)*x^-2*[1+x]^-1*z^-1
       - ln(x)*ln(1 + x*z^-1)*x^-2*[1+x]^-1
       - ln(x)*ln(1 + x*z^-1)*x^-2*[1+x]^-1*z
       - 2*ln(x)*ln(1 + x*z^-1)*x^-2*z^-1*[1-z]^-1
       + 2*ln(x)*ln(1 + x*z^-1)*x^-2*z^-1
       + ln(x)*ln(1 + x*z^-1)*x^-2
       + ln(x)*ln(1 + x*z^-1)*x^-2*z
       + 6*ln(x)*ln(1 + x*z^-1)*x^-1*[1+x]^-1*z^-1*[1-z]^-1
       - 6*ln(x)*ln(1 + x*z^-1)*x^-1*[1+x]^-1*z^-1
       - 3*ln(x)*ln(1 + x*z^-1)*x^-1*[1+x]^-1
       - 3*ln(x)*ln(1 + x*z^-1)*x^-1*[1+x]^-1*z
       - 4*ln(x)*ln(1 + x*z^-1)*x^-1*z^-1*[1-z]^-1
       + 4*ln(x)*ln(1 + x*z^-1)*x^-1*z^-1
       + 2*ln(x)*ln(1 + x*z^-1)*x^-1
       + 2*ln(x)*ln(1 + x*z^-1)*x^-1*z
       + 6*ln(x)*ln(1 + x*z^-1)*[1+x]^-1*z^-1*[1-z]^-1
       - 6*ln(x)*ln(1 + x*z^-1)*[1+x]^-1*z^-1
       - 4*ln(x)*ln(1 + x*z^-1)*[1+x]^-1
       - 2*ln(x)*ln(1 + x*z^-1)*[1+x]^-1*z
       - ln(x)*ln(1 + x*z^-1)*z
       + 2*ln(x)*ln(1 + x*z^-1)*x*[1+x]^-1*z^-1*[1-z]^-1
       - 2*ln(x)*ln(1 + x*z^-1)*x*[1+x]^-1*z^-1
       - 2*ln(x)*ln(1 + x*z^-1)*x*[1+x]^-1
       + 4*ln(x)*ln(1 + x*z^-1)*x*z^-1*[1-z]^-1
       - 4*ln(x)*ln(1 + x*z^-1)*x*z^-1
       - 2*ln(x)*ln(1 + x*z^-1)*x
       - 2*ln(x)*ln(1 + x*z^-1)*x*z
       + 2*ln(x)*ln(1 + x*z^-1)*x^2*z^-1*[1-z]^-1
       - 2*ln(x)*ln(1 + x*z^-1)*x^2*z^-1
       - 2*ln(x)*ln(1 + x*z^-1)*x^2
       + 2*ln(x)*ln(1 + x*z)*x^-2*[1+x]^-1*z^-1*[1-z]^-1
       - 2*ln(x)*ln(1 + x*z)*x^-2*[1+x]^-1*z^-1
       - ln(x)*ln(1 + x*z)*x^-2*[1+x]^-1
       - ln(x)*ln(1 + x*z)*x^-2*[1+x]^-1*z
       - 2*ln(x)*ln(1 + x*z)*x^-2*z^-1*[1-z]^-1
       + 2*ln(x)*ln(1 + x*z)*x^-2*z^-1
       + ln(x)*ln(1 + x*z)*x^-2
       + ln(x)*ln(1 + x*z)*x^-2*z
       + 6*ln(x)*ln(1 + x*z)*x^-1*[1+x]^-1*z^-1*[1-z]^-1
       - 6*ln(x)*ln(1 + x*z)*x^-1*[1+x]^-1*z^-1
       - 3*ln(x)*ln(1 + x*z)*x^-1*[1+x]^-1
       - 3*ln(x)*ln(1 + x*z)*x^-1*[1+x]^-1*z
       - 4*ln(x)*ln(1 + x*z)*x^-1*z^-1*[1-z]^-1
       + 4*ln(x)*ln(1 + x*z)*x^-1*z^-1
       + 2*ln(x)*ln(1 + x*z)*x^-1
       + 2*ln(x)*ln(1 + x*z)*x^-1*z
       + 6*ln(x)*ln(1 + x*z)*[1+x]^-1*z^-1*[1-z]^-1
       - 6*ln(x)*ln(1 + x*z)*[1+x]^-1*z^-1
       - 4*ln(x)*ln(1 + x*z)*[1+x]^-1
       - 2*ln(x)*ln(1 + x*z)*[1+x]^-1*z
       - 2*ln(x)*ln(1 + x*z)*z^-1*[1+z]^-1
       + 2*ln(x)*ln(1 + x*z)*z^-1
       + 2*ln(x)*ln(1 + x*z)*x*[1+x]^-1*z^-1*[1-z]^-1
       - 2*ln(x)*ln(1 + x*z)*x*[1+x]^-1*z^-1
       - 2*ln(x)*ln(1 + x*z)*x*[1+x]^-1
       + 4*ln(x)*ln(1 + x*z)*x*z^-1*[1-z]^-1
       + 4*ln(x)*ln(1 + x*z)*x*z^-1*[1+z]^-1
       - 8*ln(x)*ln(1 + x*z)*x*z^-1
       - 4*ln(x)*ln(1 + x*z)*x*z
       + 2*ln(x)*ln(1 + x*z)*x^2*z^-1*[1-z]^-1
       - 2*ln(x)*ln(1 + x*z)*x^2*z^-1*[1+z]^-1
       - 4*ln(x)*ln(1 + x*z)*x^2
       + ln(x)*ln(z + x)
       + 2*ln(x)*ln(z + x)*z^-1*[1+z]^-1
       - 2*ln(x)*ln(z + x)*z^-1
       - ln(x)*ln(z + x)*z
       - 4*ln(x)*ln(z + x)*x*z^-1*[1+z]^-1
       + 4*ln(x)*ln(z + x)*x*z^-1
       - 2*ln(x)*ln(z + x)*x
       + 2*ln(x)*ln(z + x)*x*z
       + 2*ln(x)*ln(z + x)*x^2*z^-1*[1+z]^-1
       - 2*ln(x)*ln(z + x)*x^2*z^-1
       + 2*ln(x)*ln(z + x)*x^2
       + 2*ln(x)^2
       + 3*ln(x)^2*x^-2*[1+x]^-1*z^-1*[1-z]^-1
       - 3*ln(x)^2*x^-2*[1+x]^-1*z^-1
       - 3/2*ln(x)^2*x^-2*[1+x]^-1
       - 3/2*ln(x)^2*x^-2*[1+x]^-1*z
       - 3*ln(x)^2*x^-2*z^-1*[1-z]^-1
       + 3*ln(x)^2*x^-2*z^-1
       + 3/2*ln(x)^2*x^-2
       + 3/2*ln(x)^2*x^-2*z
       + 9*ln(x)^2*x^-1*[1+x]^-1*z^-1*[1-z]^-1
       - 9*ln(x)^2*x^-1*[1+x]^-1*z^-1
       - 9/2*ln(x)^2*x^-1*[1+x]^-1
       - 9/2*ln(x)^2*x^-1*[1+x]^-1*z
       - 6*ln(x)^2*x^-1*z^-1*[1-z]^-1
       + 6*ln(x)^2*x^-1*z^-1
       + 3*ln(x)^2*x^-1
       + 3*ln(x)^2*x^-1*z
       + 9*ln(x)^2*[1+x]^-1*z^-1*[1-z]^-1
       - 9*ln(x)^2*[1+x]^-1*z^-1
       - 6*ln(x)^2*[1+x]^-1
       - 3*ln(x)^2*[1+x]^-1*z
       + 3*ln(x)^2*z^-1*[1+z]^-1*sqrtxz1
       + ln(x)^2*z^-1
       - 3*ln(x)^2*z^-1*sqrtxz1
       - 3*ln(x)^2*[1-z]^-1
       + 9*ln(x)^2*[1+z]^-1*sqrtxz1^-1
       - 6*ln(x)^2*sqrtxz1^-1
       - 3*ln(x)^2*z*[1+z]^-1*sqrtxz1^-1
       + 3*ln(x)^2*z*sqrtxz1^-1
       + 3*ln(x)^2*x*[1+x]^-1*z^-1*[1-z]^-1
       - 3*ln(x)^2*x*[1+x]^-1*z^-1
       - 3*ln(x)^2*x*[1+x]^-1
       + 6*ln(x)^2*x*z^-1*[1-z]^-1
       - 6*ln(x)^2*x*z^-1*[1+z]^-1*sqrtxz1
       - 8*ln(x)^2*x*z^-1
       + 6*ln(x)^2*x*z^-1*sqrtxz1
       - 6*ln(x)^2*x*[1-z]^-1
       - 30*ln(x)^2*x*[1+z]^-1*sqrtxz1^-1
       + 24*ln(x)^2*x*sqrtxz1^-1
       + 2*ln(x)^2*x
       + 6*ln(x)^2*x*z*[1+z]^-1*sqrtxz1^-1
       - 6*ln(x)^2*x*z*sqrtxz1^-1
       - 2*ln(x)^2*x*z
       + 3*ln(x)^2*x^2*z^-1*[1-z]^-1
       + 3*ln(x)^2*x^2*z^-1*[1+z]^-1*sqrtxz1
       - 3*ln(x)^2*x^2*z^-1
       - 3*ln(x)^2*x^2*z^-1*sqrtxz1
       - 3*ln(x)^2*x^2*[1-z]^-1
       + 33*ln(x)^2*x^2*[1+z]^-1*sqrtxz1^-1
       - 30*ln(x)^2*x^2*sqrtxz1^-1
       - ln(x)^2*x^2
       - 3*ln(x)^2*x^2*z*[1+z]^-1*sqrtxz1^-1
       + 3*ln(x)^2*x^2*z*sqrtxz1^-1
       - 12*ln(x)^2*x^3*[1+z]^-1*sqrtxz1^-1
       + 12*ln(x)^2*x^3*sqrtxz1^-1
       + 4*ln(x)*ln([1-x])
       + 4*ln(x)*ln([1-x])*x^-2*[1+x]^-1*z^-1*[1-z]^-1
       - 4*ln(x)*ln([1-x])*x^-2*[1+x]^-1*z^-1
       - 2*ln(x)*ln([1-x])*x^-2*[1+x]^-1
       - 2*ln(x)*ln([1-x])*x^-2*[1+x]^-1*z
       - 4*ln(x)*ln([1-x])*x^-2*z^-1*[1-z]^-1
       + 4*ln(x)*ln([1-x])*x^-2*z^-1
       + 2*ln(x)*ln([1-x])*x^-2
       + 2*ln(x)*ln([1-x])*x^-2*z
       + 12*ln(x)*ln([1-x])*x^-1*[1+x]^-1*z^-1*[1-z]^-1
       - 12*ln(x)*ln([1-x])*x^-1*[1+x]^-1*z^-1
       - 6*ln(x)*ln([1-x])*x^-1*[1+x]^-1
       - 6*ln(x)*ln([1-x])*x^-1*[1+x]^-1*z
       - 8*ln(x)*ln([1-x])*x^-1*z^-1*[1-z]^-1
       + 8*ln(x)*ln([1-x])*x^-1*z^-1
       + 4*ln(x)*ln([1-x])*x^-1
       + 4*ln(x)*ln([1-x])*x^-1*z
       + 12*ln(x)*ln([1-x])*[1+x]^-1*z^-1*[1-z]^-1
       - 12*ln(x)*ln([1-x])*[1+x]^-1*z^-1
       - 8*ln(x)*ln([1-x])*[1+x]^-1
       - 4*ln(x)*ln([1-x])*[1+x]^-1*z
       - ln(x)*ln([1-x])*z^-1
       - 4*ln(x)*ln([1-x])*[1-z]^-1
       + 4*ln(x)*ln([1-x])*x*[1+x]^-1*z^-1*[1-z]^-1
       - 4*ln(x)*ln([1-x])*x*[1+x]^-1*z^-1
       - 4*ln(x)*ln([1-x])*x*[1+x]^-1
       + 8*ln(x)*ln([1-x])*x*z^-1*[1-z]^-1
       - 6*ln(x)*ln([1-x])*x*z^-1
       - 8*ln(x)*ln([1-x])*x*[1-z]^-1
       + 4*ln(x)*ln([1-x])*x^2*z^-1*[1-z]^-1
       - 4*ln(x)*ln([1-x])*x^2*z^-1
       - 4*ln(x)*ln([1-x])*x^2*[1-z]^-1
       - 2*ln(x)*ln([1-x])*x^2
       - 2*ln(x)*ln([1+x])
       - 4*ln(x)*ln([1+x])*x^-2*[1+x]^-1*z^-1*[1-z]^-1
       + 4*ln(x)*ln([1+x])*x^-2*[1+x]^-1*z^-1
       + 2*ln(x)*ln([1+x])*x^-2*[1+x]^-1
       + 2*ln(x)*ln([1+x])*x^-2*[1+x]^-1*z
       + 4*ln(x)*ln([1+x])*x^-2*z^-1*[1-z]^-1
       - 4*ln(x)*ln([1+x])*x^-2*z^-1
       - 2*ln(x)*ln([1+x])*x^-2
       - 2*ln(x)*ln([1+x])*x^-2*z
       - 12*ln(x)*ln([1+x])*x^-1*[1+x]^-1*z^-1*[1-z]^-1
       + 12*ln(x)*ln([1+x])*x^-1*[1+x]^-1*z^-1
       + 6*ln(x)*ln([1+x])*x^-1*[1+x]^-1
       + 6*ln(x)*ln([1+x])*x^-1*[1+x]^-1*z
       + 8*ln(x)*ln([1+x])*x^-1*z^-1*[1-z]^-1
       - 8*ln(x)*ln([1+x])*x^-1*z^-1
       - 4*ln(x)*ln([1+x])*x^-1
       - 4*ln(x)*ln([1+x])*x^-1*z
       - 12*ln(x)*ln([1+x])*[1+x]^-1*z^-1*[1-z]^-1
       + 12*ln(x)*ln([1+x])*[1+x]^-1*z^-1
       + 8*ln(x)*ln([1+x])*[1+x]^-1
       + 4*ln(x)*ln([1+x])*[1+x]^-1*z
       + 2*ln(x)*ln([1+x])*z
       - 4*ln(x)*ln([1+x])*x*[1+x]^-1*z^-1*[1-z]^-1
       + 4*ln(x)*ln([1+x])*x*[1+x]^-1*z^-1
       + 4*ln(x)*ln([1+x])*x*[1+x]^-1
       - 8*ln(x)*ln([1+x])*x*z^-1*[1-z]^-1
       + 8*ln(x)*ln([1+x])*x*z^-1
       + 4*ln(x)*ln([1+x])*x
       + 4*ln(x)*ln([1+x])*x*z
       - 4*ln(x)*ln([1+x])*x^2*z^-1*[1-z]^-1
       + 4*ln(x)*ln([1+x])*x^2*z^-1
       + 4*ln(x)*ln([1+x])*x^2
       - 2*ln(x)*ln(z)
       - 2*ln(x)*ln(z)*x^-2*[1+x]^-1*z^-1*[1-z]^-1
       + 2*ln(x)*ln(z)*x^-2*[1+x]^-1*z^-1
       + ln(x)*ln(z)*x^-2*[1+x]^-1
       + ln(x)*ln(z)*x^-2*[1+x]^-1*z
       + 2*ln(x)*ln(z)*x^-2*z^-1*[1-z]^-1
       - 2*ln(x)*ln(z)*x^-2*z^-1
       - ln(x)*ln(z)*x^-2
       - ln(x)*ln(z)*x^-2*z
       - 6*ln(x)*ln(z)*x^-1*[1+x]^-1*z^-1*[1-z]^-1
       + 6*ln(x)*ln(z)*x^-1*[1+x]^-1*z^-1
       + 3*ln(x)*ln(z)*x^-1*[1+x]^-1
       + 3*ln(x)*ln(z)*x^-1*[1+x]^-1*z
       + 4*ln(x)*ln(z)*x^-1*z^-1*[1-z]^-1
       - 4*ln(x)*ln(z)*x^-1*z^-1
       - 2*ln(x)*ln(z)*x^-1
       - 2*ln(x)*ln(z)*x^-1*z
       - 6*ln(x)*ln(z)*[1+x]^-1*z^-1*[1-z]^-1
       + 6*ln(x)*ln(z)*[1+x]^-1*z^-1
       + 4*ln(x)*ln(z)*[1+x]^-1
       + 2*ln(x)*ln(z)*[1+x]^-1*z
       - 2*ln(x)*ln(z)*z^-1*[1+z]^-1
       + 2*ln(x)*ln(z)*z^-1*[1+z]^-1*sqrtxz1
       + 5/2*ln(x)*ln(z)*z^-1
       - 2*ln(x)*ln(z)*z^-1*sqrtxz1
       + 2*ln(x)*ln(z)*[1-z]^-1
       + 6*ln(x)*ln(z)*[1+z]^-1*sqrtxz1^-1
       - 4*ln(x)*ln(z)*sqrtxz1^-1
       - 2*ln(x)*ln(z)*z*[1+z]^-1*sqrtxz1^-1
       + 2*ln(x)*ln(z)*z*sqrtxz1^-1
       + 2*ln(x)*ln(z)*z
       - 2*ln(x)*ln(z)*x*[1+x]^-1*z^-1*[1-z]^-1
       + 2*ln(x)*ln(z)*x*[1+x]^-1*z^-1
       + 2*ln(x)*ln(z)*x*[1+x]^-1
       - 4*ln(x)*ln(z)*x*z^-1*[1-z]^-1
       + 4*ln(x)*ln(z)*x*z^-1*[1+z]^-1
       - 4*ln(x)*ln(z)*x*z^-1*[1+z]^-1*sqrtxz1
       - ln(x)*ln(z)*x*z^-1
       + 4*ln(x)*ln(z)*x*z^-1*sqrtxz1
       + 4*ln(x)*ln(z)*x*[1-z]^-1
       - 20*ln(x)*ln(z)*x*[1+z]^-1*sqrtxz1^-1
       + 16*ln(x)*ln(z)*x*sqrtxz1^-1
       + 4*ln(x)*ln(z)*x
       + 4*ln(x)*ln(z)*x*z*[1+z]^-1*sqrtxz1^-1
       - 4*ln(x)*ln(z)*x*z*sqrtxz1^-1
       - 2*ln(x)*ln(z)*x*z
       - 2*ln(x)*ln(z)*x^2*z^-1*[1-z]^-1
       - 2*ln(x)*ln(z)*x^2*z^-1*[1+z]^-1
       + 2*ln(x)*ln(z)*x^2*z^-1*[1+z]^-1*sqrtxz1
       + 4*ln(x)*ln(z)*x^2*z^-1
       - 2*ln(x)*ln(z)*x^2*z^-1*sqrtxz1
       + 2*ln(x)*ln(z)*x^2*[1-z]^-1
       + 22*ln(x)*ln(z)*x^2*[1+z]^-1*sqrtxz1^-1
       - 20*ln(x)*ln(z)*x^2*sqrtxz1^-1
       - 2*ln(x)*ln(z)*x^2
       - 2*ln(x)*ln(z)*x^2*z*[1+z]^-1*sqrtxz1^-1
       + 2*ln(x)*ln(z)*x^2*z*sqrtxz1^-1
       - 8*ln(x)*ln(z)*x^3*[1+z]^-1*sqrtxz1^-1
       + 8*ln(x)*ln(z)*x^3*sqrtxz1^-1
       + ln(x)*ln([1-z])
       - 4*ln(x)*ln([1-z])*z^-1*[1+z]^-1*sqrtxz1
       - 3/2*ln(x)*ln([1-z])*z^-1
       + 4*ln(x)*ln([1-z])*z^-1*sqrtxz1
       - 12*ln(x)*ln([1-z])*[1+z]^-1*sqrtxz1^-1
       + 8*ln(x)*ln([1-z])*sqrtxz1^-1
       + 4*ln(x)*ln([1-z])*z*[1+z]^-1*sqrtxz1^-1
       - 4*ln(x)*ln([1-z])*z*sqrtxz1^-1
       + 8*ln(x)*ln([1-z])*x*z^-1*[1+z]^-1*sqrtxz1
       + 3*ln(x)*ln([1-z])*x*z^-1
       - 8*ln(x)*ln([1-z])*x*z^-1*sqrtxz1
       + 40*ln(x)*ln([1-z])*x*[1+z]^-1*sqrtxz1^-1
       - 32*ln(x)*ln([1-z])*x*sqrtxz1^-1
       - 2*ln(x)*ln([1-z])*x
       - 8*ln(x)*ln([1-z])*x*z*[1+z]^-1*sqrtxz1^-1
       + 8*ln(x)*ln([1-z])*x*z*sqrtxz1^-1
       + 2*ln(x)*ln([1-z])*x*z
       - 4*ln(x)*ln([1-z])*x^2*z^-1*[1+z]^-1*sqrtxz1
       + 4*ln(x)*ln([1-z])*x^2*z^-1*sqrtxz1
       - 44*ln(x)*ln([1-z])*x^2*[1+z]^-1*sqrtxz1^-1
       + 40*ln(x)*ln([1-z])*x^2*sqrtxz1^-1
       + 4*ln(x)*ln([1-z])*x^2*z*[1+z]^-1*sqrtxz1^-1
       - 4*ln(x)*ln([1-z])*x^2*z*sqrtxz1^-1
       + 16*ln(x)*ln([1-z])*x^3*[1+z]^-1*sqrtxz1^-1
       - 16*ln(x)*ln([1-z])*x^3*sqrtxz1^-1
       + 11/2*ln([1-x])
       - 25/4*ln([1-x])*z^-1
       + 3/4*ln([1-x])*z
       + 9*ln([1-x])*x*z^-1
       - 8*ln([1-x])*x
       - 3/2*ln([1-x])*x*z
       - ln([1-x])^2
       + ln([1-x])^2*z^-1
       + 1/2*ln([1-x])^2*z
       - 2*ln([1-x])^2*x*z^-1
       + 2*ln([1-x])^2*x
       - ln([1-x])^2*x*z
       - 2*ln([1-x])*ln(z)
       + 4*ln([1-x])*ln(z)*x
       - 2*ln([1-x])*ln([1-z])
       + 2*ln([1-x])*ln([1-z])*z^-1
       + ln([1-x])*ln([1-z])*z
       - 4*ln([1-x])*ln([1-z])*x*z^-1
       + 4*ln([1-x])*ln([1-z])*x
       - 2*ln([1-x])*ln([1-z])*x*z
       + 3/2*ln(z)
       - 12*ln(z)*z^-1*[1+z]^-1*rln2
       - 8*ln(z)*z^-1*[1+z]^-1*sqrtxz1
       - 16*ln(z)*z^-1*[1+z]^-1*sqrtxz1*rln2
       - 33/4*ln(z)*z^-1
       + 10*ln(z)*z^-1*rln2
       + 8*ln(z)*z^-1*sqrtxz1
       + 16*ln(z)*z^-1*sqrtxz1*rln2
       - 24*ln(z)*[1+z]^-1*sqrtxz1^-1
       - 48*ln(z)*[1+z]^-1*sqrtxz1^-1*rln2
       + 16*ln(z)*sqrtxz1^-1
       + 32*ln(z)*sqrtxz1^-1*rln2
       - 8*ln(z)*rln2
       + 8*ln(z)*z*[1+z]^-1*sqrtxz1^-1
       + 16*ln(z)*z*[1+z]^-1*sqrtxz1^-1*rln2
       - 8*ln(z)*z*sqrtxz1^-1
       - 16*ln(z)*z*sqrtxz1^-1*rln2
       - 9/4*ln(z)*z
       + 4*ln(z)*z*rln2
       + 24*ln(z)*x*z^-1*[1+z]^-1*rln2
       + 16*ln(z)*x*z^-1*[1+z]^-1*sqrtxz1
       + 32*ln(z)*x*z^-1*[1+z]^-1*sqrtxz1*rln2
       + 13*ln(z)*x*z^-1
       - 20*ln(z)*x*z^-1*rln2
       - 16*ln(z)*x*z^-1*sqrtxz1
       - 32*ln(z)*x*z^-1*sqrtxz1*rln2
       + 80*ln(z)*x*[1+z]^-1*sqrtxz1^-1
       + 160*ln(z)*x*[1+z]^-1*sqrtxz1^-1*rln2
       - 64*ln(z)*x*sqrtxz1^-1
       - 128*ln(z)*x*sqrtxz1^-1*rln2
       + 16*ln(z)*x*rln2
       - 16*ln(z)*x*z*[1+z]^-1*sqrtxz1^-1
       - 32*ln(z)*x*z*[1+z]^-1*sqrtxz1^-1*rln2
       + 16*ln(z)*x*z*sqrtxz1^-1
       + 32*ln(z)*x*z*sqrtxz1^-1*rln2
       + 5/2*ln(z)*x*z
       - 8*ln(z)*x*z*rln2
       - 12*ln(z)*x^2*z^-1*[1+z]^-1*rln2
       - 8*ln(z)*x^2*z^-1*[1+z]^-1*sqrtxz1
       - 16*ln(z)*x^2*z^-1*[1+z]^-1*sqrtxz1*rln2
       + 12*ln(z)*x^2*z^-1*rln2
       + 8*ln(z)*x^2*z^-1*sqrtxz1
       + 16*ln(z)*x^2*z^-1*sqrtxz1*rln2
       - 88*ln(z)*x^2*[1+z]^-1*sqrtxz1^-1
       - 176*ln(z)*x^2*[1+z]^-1*sqrtxz1^-1*rln2
       + 80*ln(z)*x^2*sqrtxz1^-1
       + 160*ln(z)*x^2*sqrtxz1^-1*rln2
       - 8*ln(z)*x^2*rln2
       + 8*ln(z)*x^2*z*[1+z]^-1*sqrtxz1^-1
       + 16*ln(z)*x^2*z*[1+z]^-1*sqrtxz1^-1*rln2
       - 8*ln(z)*x^2*z*sqrtxz1^-1
       - 16*ln(z)*x^2*z*sqrtxz1^-1*rln2
       + 32*ln(z)*x^3*[1+z]^-1*sqrtxz1^-1
       + 64*ln(z)*x^3*[1+z]^-1*sqrtxz1^-1*rln2
       - 32*ln(z)*x^3*sqrtxz1^-1
       - 64*ln(z)*x^3*sqrtxz1^-1*rln2
       + 4*ln(z)*ln(1 + sqrtxz1 - z)
       + 4*ln(z)*ln(1 + sqrtxz1 - z)*z^-1*[1+z]^-1
       + 8*ln(z)*ln(1 + sqrtxz1 - z)*z^-1*[1+z]^-1*sqrtxz1
       - 4*ln(z)*ln(1 + sqrtxz1 - z)*z^-1
       - 8*ln(z)*ln(1 + sqrtxz1 - z)*z^-1*sqrtxz1
       + 24*ln(z)*ln(1 + sqrtxz1 - z)*[1+z]^-1*sqrtxz1^-1
       - 16*ln(z)*ln(1 + sqrtxz1 - z)*sqrtxz1^-1
       - 8*ln(z)*ln(1 + sqrtxz1 - z)*z*[1+z]^-1*sqrtxz1^-1
       + 8*ln(z)*ln(1 + sqrtxz1 - z)*z*sqrtxz1^-1
       - 2*ln(z)*ln(1 + sqrtxz1 - z)*z
       - 8*ln(z)*ln(1 + sqrtxz1 - z)*x*z^-1*[1+z]^-1
       - 16*ln(z)*ln(1 + sqrtxz1 - z)*x*z^-1*[1+z]^-1*sqrtxz1
       + 8*ln(z)*ln(1 + sqrtxz1 - z)*x*z^-1
       + 16*ln(z)*ln(1 + sqrtxz1 - z)*x*z^-1*sqrtxz1
       - 80*ln(z)*ln(1 + sqrtxz1 - z)*x*[1+z]^-1*sqrtxz1^-1
       + 64*ln(z)*ln(1 + sqrtxz1 - z)*x*sqrtxz1^-1
       - 8*ln(z)*ln(1 + sqrtxz1 - z)*x
       + 16*ln(z)*ln(1 + sqrtxz1 - z)*x*z*[1+z]^-1*sqrtxz1^-1
       - 16*ln(z)*ln(1 + sqrtxz1 - z)*x*z*sqrtxz1^-1
       + 4*ln(z)*ln(1 + sqrtxz1 - z)*x*z
       + 4*ln(z)*ln(1 + sqrtxz1 - z)*x^2*z^-1*[1+z]^-1
       + 8*ln(z)*ln(1 + sqrtxz1 - z)*x^2*z^-1*[1+z]^-1*sqrtxz1
       - 4*ln(z)*ln(1 + sqrtxz1 - z)*x^2*z^-1
       - 8*ln(z)*ln(1 + sqrtxz1 - z)*x^2*z^-1*sqrtxz1
       + 88*ln(z)*ln(1 + sqrtxz1 - z)*x^2*[1+z]^-1*sqrtxz1^-1
       - 80*ln(z)*ln(1 + sqrtxz1 - z)*x^2*sqrtxz1^-1
       + 4*ln(z)*ln(1 + sqrtxz1 - z)*x^2
       - 8*ln(z)*ln(1 + sqrtxz1 - z)*x^2*z*[1+z]^-1*sqrtxz1^-1
       + 8*ln(z)*ln(1 + sqrtxz1 - z)*x^2*z*sqrtxz1^-1
       - 32*ln(z)*ln(1 + sqrtxz1 - z)*x^3*[1+z]^-1*sqrtxz1^-1
       + 32*ln(z)*ln(1 + sqrtxz1 - z)*x^3*sqrtxz1^-1
       + 4*ln(z)*ln(1 + sqrtxz1 + z)
       + 8*ln(z)*ln(1 + sqrtxz1 + z)*z^-1*[1+z]^-1
       + 8*ln(z)*ln(1 + sqrtxz1 + z)*z^-1*[1+z]^-1*sqrtxz1
       - 6*ln(z)*ln(1 + sqrtxz1 + z)*z^-1
       - 8*ln(z)*ln(1 + sqrtxz1 + z)*z^-1*sqrtxz1
       + 24*ln(z)*ln(1 + sqrtxz1 + z)*[1+z]^-1*sqrtxz1^-1
       - 16*ln(z)*ln(1 + sqrtxz1 + z)*sqrtxz1^-1
       - 8*ln(z)*ln(1 + sqrtxz1 + z)*z*[1+z]^-1*sqrtxz1^-1
       + 8*ln(z)*ln(1 + sqrtxz1 + z)*z*sqrtxz1^-1
       - 2*ln(z)*ln(1 + sqrtxz1 + z)*z
       - 16*ln(z)*ln(1 + sqrtxz1 + z)*x*z^-1*[1+z]^-1
       - 16*ln(z)*ln(1 + sqrtxz1 + z)*x*z^-1*[1+z]^-1*sqrtxz1
       + 12*ln(z)*ln(1 + sqrtxz1 + z)*x*z^-1
       + 16*ln(z)*ln(1 + sqrtxz1 + z)*x*z^-1*sqrtxz1
       - 80*ln(z)*ln(1 + sqrtxz1 + z)*x*[1+z]^-1*sqrtxz1^-1
       + 64*ln(z)*ln(1 + sqrtxz1 + z)*x*sqrtxz1^-1
       - 8*ln(z)*ln(1 + sqrtxz1 + z)*x
       + 16*ln(z)*ln(1 + sqrtxz1 + z)*x*z*[1+z]^-1*sqrtxz1^-1
       - 16*ln(z)*ln(1 + sqrtxz1 + z)*x*z*sqrtxz1^-1
       + 4*ln(z)*ln(1 + sqrtxz1 + z)*x*z
       + 8*ln(z)*ln(1 + sqrtxz1 + z)*x^2*z^-1*[1+z]^-1
       + 8*ln(z)*ln(1 + sqrtxz1 + z)*x^2*z^-1*[1+z]^-1*sqrtxz1
       - 8*ln(z)*ln(1 + sqrtxz1 + z)*x^2*z^-1
       - 8*ln(z)*ln(1 + sqrtxz1 + z)*x^2*z^-1*sqrtxz1
       + 88*ln(z)*ln(1 + sqrtxz1 + z)*x^2*[1+z]^-1*sqrtxz1^-1
       - 80*ln(z)*ln(1 + sqrtxz1 + z)*x^2*sqrtxz1^-1
       + 4*ln(z)*ln(1 + sqrtxz1 + z)*x^2
       - 8*ln(z)*ln(1 + sqrtxz1 + z)*x^2*z*[1+z]^-1*sqrtxz1^-1
       + 8*ln(z)*ln(1 + sqrtxz1 + z)*x^2*z*sqrtxz1^-1
       - 32*ln(z)*ln(1 + sqrtxz1 + z)*x^3*[1+z]^-1*sqrtxz1^-1
       + 32*ln(z)*ln(1 + sqrtxz1 + z)*x^3*sqrtxz1^-1
       - ln(z)*ln(1 + x*z^-1)
       - 2*ln(z)*ln(1 + x*z^-1)*x^-2*[1+x]^-1*z^-1*[1-z]^-1
       + 2*ln(z)*ln(1 + x*z^-1)*x^-2*[1+x]^-1*z^-1
       + ln(z)*ln(1 + x*z^-1)*x^-2*[1+x]^-1
       + ln(z)*ln(1 + x*z^-1)*x^-2*[1+x]^-1*z
       + 2*ln(z)*ln(1 + x*z^-1)*x^-2*z^-1*[1-z]^-1
       - 2*ln(z)*ln(1 + x*z^-1)*x^-2*z^-1
       - ln(z)*ln(1 + x*z^-1)*x^-2
       - ln(z)*ln(1 + x*z^-1)*x^-2*z
       - 6*ln(z)*ln(1 + x*z^-1)*x^-1*[1+x]^-1*z^-1*[1-z]^-1
       + 6*ln(z)*ln(1 + x*z^-1)*x^-1*[1+x]^-1*z^-1
       + 3*ln(z)*ln(1 + x*z^-1)*x^-1*[1+x]^-1
       + 3*ln(z)*ln(1 + x*z^-1)*x^-1*[1+x]^-1*z
       + 4*ln(z)*ln(1 + x*z^-1)*x^-1*z^-1*[1-z]^-1
       - 4*ln(z)*ln(1 + x*z^-1)*x^-1*z^-1
       - 2*ln(z)*ln(1 + x*z^-1)*x^-1
       - 2*ln(z)*ln(1 + x*z^-1)*x^-1*z
       - 6*ln(z)*ln(1 + x*z^-1)*[1+x]^-1*z^-1*[1-z]^-1
       + 6*ln(z)*ln(1 + x*z^-1)*[1+x]^-1*z^-1
       + 4*ln(z)*ln(1 + x*z^-1)*[1+x]^-1
       + 2*ln(z)*ln(1 + x*z^-1)*[1+x]^-1*z
       + ln(z)*ln(1 + x*z^-1)*z
       - 2*ln(z)*ln(1 + x*z^-1)*x*[1+x]^-1*z^-1*[1-z]^-1
       + 2*ln(z)*ln(1 + x*z^-1)*x*[1+x]^-1*z^-1
       + 2*ln(z)*ln(1 + x*z^-1)*x*[1+x]^-1
       - 4*ln(z)*ln(1 + x*z^-1)*x*z^-1*[1-z]^-1
       + 4*ln(z)*ln(1 + x*z^-1)*x*z^-1
       + 2*ln(z)*ln(1 + x*z^-1)*x
       + 2*ln(z)*ln(1 + x*z^-1)*x*z
       - 2*ln(z)*ln(1 + x*z^-1)*x^2*z^-1*[1-z]^-1
       + 2*ln(z)*ln(1 + x*z^-1)*x^2*z^-1
       + 2*ln(z)*ln(1 + x*z^-1)*x^2
       + 2*ln(z)*ln(1 + x*z)*x^-2*[1+x]^-1*z^-1*[1-z]^-1
       - 2*ln(z)*ln(1 + x*z)*x^-2*[1+x]^-1*z^-1
       - ln(z)*ln(1 + x*z)*x^-2*[1+x]^-1
       - ln(z)*ln(1 + x*z)*x^-2*[1+x]^-1*z
       - 2*ln(z)*ln(1 + x*z)*x^-2*z^-1*[1-z]^-1
       + 2*ln(z)*ln(1 + x*z)*x^-2*z^-1
       + ln(z)*ln(1 + x*z)*x^-2
       + ln(z)*ln(1 + x*z)*x^-2*z
       + 6*ln(z)*ln(1 + x*z)*x^-1*[1+x]^-1*z^-1*[1-z]^-1
       - 6*ln(z)*ln(1 + x*z)*x^-1*[1+x]^-1*z^-1
       - 3*ln(z)*ln(1 + x*z)*x^-1*[1+x]^-1
       - 3*ln(z)*ln(1 + x*z)*x^-1*[1+x]^-1*z
       - 4*ln(z)*ln(1 + x*z)*x^-1*z^-1*[1-z]^-1
       + 4*ln(z)*ln(1 + x*z)*x^-1*z^-1
       + 2*ln(z)*ln(1 + x*z)*x^-1
       + 2*ln(z)*ln(1 + x*z)*x^-1*z
       + 6*ln(z)*ln(1 + x*z)*[1+x]^-1*z^-1*[1-z]^-1
       - 6*ln(z)*ln(1 + x*z)*[1+x]^-1*z^-1
       - 4*ln(z)*ln(1 + x*z)*[1+x]^-1
       - 2*ln(z)*ln(1 + x*z)*[1+x]^-1*z
       - 2*ln(z)*ln(1 + x*z)*z^-1*[1+z]^-1
       + 2*ln(z)*ln(1 + x*z)*z^-1
       + 2*ln(z)*ln(1 + x*z)*x*[1+x]^-1*z^-1*[1-z]^-1
       - 2*ln(z)*ln(1 + x*z)*x*[1+x]^-1*z^-1
       - 2*ln(z)*ln(1 + x*z)*x*[1+x]^-1
       + 4*ln(z)*ln(1 + x*z)*x*z^-1*[1-z]^-1
       + 4*ln(z)*ln(1 + x*z)*x*z^-1*[1+z]^-1
       - 8*ln(z)*ln(1 + x*z)*x*z^-1
       - 4*ln(z)*ln(1 + x*z)*x*z
       + 2*ln(z)*ln(1 + x*z)*x^2*z^-1*[1-z]^-1
       - 2*ln(z)*ln(1 + x*z)*x^2*z^-1*[1+z]^-1
       - 4*ln(z)*ln(1 + x*z)*x^2
       - ln(z)*ln(z + x)
       - 2*ln(z)*ln(z + x)*z^-1*[1+z]^-1
       + 2*ln(z)*ln(z + x)*z^-1
       + ln(z)*ln(z + x)*z
       + 4*ln(z)*ln(z + x)*x*z^-1*[1+z]^-1
       - 4*ln(z)*ln(z + x)*x*z^-1
       + 2*ln(z)*ln(z + x)*x
       - 2*ln(z)*ln(z + x)*x*z
       - 2*ln(z)*ln(z + x)*x^2*z^-1*[1+z]^-1
       + 2*ln(z)*ln(z + x)*x^2*z^-1
       - 2*ln(z)*ln(z + x)*x^2
       - ln(z)^2
       - ln(z)^2*x^-2*[1+x]^-1*z^-1*[1-z]^-1
       + ln(z)^2*x^-2*[1+x]^-1*z^-1
       + 1/2*ln(z)^2*x^-2*[1+x]^-1
       + 1/2*ln(z)^2*x^-2*[1+x]^-1*z
       + ln(z)^2*x^-2*z^-1*[1-z]^-1
       - ln(z)^2*x^-2*z^-1
       - 1/2*ln(z)^2*x^-2
       - 1/2*ln(z)^2*x^-2*z
       - 3*ln(z)^2*x^-1*[1+x]^-1*z^-1*[1-z]^-1
       + 3*ln(z)^2*x^-1*[1+x]^-1*z^-1
       + 3/2*ln(z)^2*x^-1*[1+x]^-1
       + 3/2*ln(z)^2*x^-1*[1+x]^-1*z
       + 2*ln(z)^2*x^-1*z^-1*[1-z]^-1
       - 2*ln(z)^2*x^-1*z^-1
       - ln(z)^2*x^-1
       - ln(z)^2*x^-1*z
       - 3*ln(z)^2*[1+x]^-1*z^-1*[1-z]^-1
       + 3*ln(z)^2*[1+x]^-1*z^-1
       + 2*ln(z)^2*[1+x]^-1
       + ln(z)^2*[1+x]^-1*z
       - 2*ln(z)^2*z^-1*[1+z]^-1
       - 5*ln(z)^2*z^-1*[1+z]^-1*sqrtxz1
       - ln(z)^2*z^-1
       + 5*ln(z)^2*z^-1*sqrtxz1
       - ln(z)^2*[1-z]^-1
       - 15*ln(z)^2*[1+z]^-1*sqrtxz1^-1
       + 10*ln(z)^2*sqrtxz1^-1
       + 5*ln(z)^2*z*[1+z]^-1*sqrtxz1^-1
       - 5*ln(z)^2*z*sqrtxz1^-1
       - 1/2*ln(z)^2*z
       - ln(z)^2*x*[1+x]^-1*z^-1*[1-z]^-1
       + ln(z)^2*x*[1+x]^-1*z^-1
       + ln(z)^2*x*[1+x]^-1
       - 2*ln(z)^2*x*z^-1*[1-z]^-1
       + 4*ln(z)^2*x*z^-1*[1+z]^-1
       + 10*ln(z)^2*x*z^-1*[1+z]^-1*sqrtxz1
       + 4*ln(z)^2*x*z^-1
       - 10*ln(z)^2*x*z^-1*sqrtxz1
       + 2*ln(z)^2*x*[1-z]^-1
       + 50*ln(z)^2*x*[1+z]^-1*sqrtxz1^-1
       - 40*ln(z)^2*x*sqrtxz1^-1
       + 2*ln(z)^2*x
       - 10*ln(z)^2*x*z*[1+z]^-1*sqrtxz1^-1
       + 10*ln(z)^2*x*z*sqrtxz1^-1
       + 3*ln(z)^2*x*z
       - ln(z)^2*x^2*z^-1*[1-z]^-1
       - 2*ln(z)^2*x^2*z^-1*[1+z]^-1
       - 5*ln(z)^2*x^2*z^-1*[1+z]^-1*sqrtxz1
       + 3*ln(z)^2*x^2*z^-1
       + 5*ln(z)^2*x^2*z^-1*sqrtxz1
       - ln(z)^2*x^2*[1-z]^-1
       - 55*ln(z)^2*x^2*[1+z]^-1*sqrtxz1^-1
       + 50*ln(z)^2*x^2*sqrtxz1^-1
       + 2*ln(z)^2*x^2
       + 5*ln(z)^2*x^2*z*[1+z]^-1*sqrtxz1^-1
       - 5*ln(z)^2*x^2*z*sqrtxz1^-1
       + 20*ln(z)^2*x^3*[1+z]^-1*sqrtxz1^-1
       - 20*ln(z)^2*x^3*sqrtxz1^-1
       - 4*ln(z)*ln([1-z])*z^-1*[1+z]^-1*sqrtxz1
       + 4*ln(z)*ln([1-z])*z^-1*sqrtxz1
       - 12*ln(z)*ln([1-z])*[1+z]^-1*sqrtxz1^-1
       + 8*ln(z)*ln([1-z])*sqrtxz1^-1
       + 4*ln(z)*ln([1-z])*z*[1+z]^-1*sqrtxz1^-1
       - 4*ln(z)*ln([1-z])*z*sqrtxz1^-1
       + 8*ln(z)*ln([1-z])*x*z^-1*[1+z]^-1*sqrtxz1
       - 8*ln(z)*ln([1-z])*x*z^-1*sqrtxz1
       + 40*ln(z)*ln([1-z])*x*[1+z]^-1*sqrtxz1^-1
       - 32*ln(z)*ln([1-z])*x*sqrtxz1^-1
       - 8*ln(z)*ln([1-z])*x*z*[1+z]^-1*sqrtxz1^-1
       + 8*ln(z)*ln([1-z])*x*z*sqrtxz1^-1
       - 4*ln(z)*ln([1-z])*x^2*z^-1*[1+z]^-1*sqrtxz1
       + 4*ln(z)*ln([1-z])*x^2*z^-1*sqrtxz1
       - 44*ln(z)*ln([1-z])*x^2*[1+z]^-1*sqrtxz1^-1
       + 40*ln(z)*ln([1-z])*x^2*sqrtxz1^-1
       + 4*ln(z)*ln([1-z])*x^2*z*[1+z]^-1*sqrtxz1^-1
       - 4*ln(z)*ln([1-z])*x^2*z*sqrtxz1^-1
       + 16*ln(z)*ln([1-z])*x^3*[1+z]^-1*sqrtxz1^-1
       - 16*ln(z)*ln([1-z])*x^3*sqrtxz1^-1
       + 11/2*ln([1-z])
       - 8*ln([1-z])*z^-1*[1+z]^-1*sqrtxz1*rln2
       - 25/4*ln([1-z])*z^-1
       + 8*ln([1-z])*z^-1*sqrtxz1*rln2
       - 24*ln([1-z])*[1+z]^-1*sqrtxz1^-1*rln2
       + 16*ln([1-z])*sqrtxz1^-1*rln2
       + 8*ln([1-z])*z*[1+z]^-1*sqrtxz1^-1*rln2
       - 8*ln([1-z])*z*sqrtxz1^-1*rln2
       + 3/4*ln([1-z])*z
       + 16*ln([1-z])*x*z^-1*[1+z]^-1*sqrtxz1*rln2
       + 9*ln([1-z])*x*z^-1
       - 16*ln([1-z])*x*z^-1*sqrtxz1*rln2
       + 80*ln([1-z])*x*[1+z]^-1*sqrtxz1^-1*rln2
       - 64*ln([1-z])*x*sqrtxz1^-1*rln2
       - 8*ln([1-z])*x
       - 16*ln([1-z])*x*z*[1+z]^-1*sqrtxz1^-1*rln2
       + 16*ln([1-z])*x*z*sqrtxz1^-1*rln2
       - 3/2*ln([1-z])*x*z
       - 8*ln([1-z])*x^2*z^-1*[1+z]^-1*sqrtxz1*rln2
       + 8*ln([1-z])*x^2*z^-1*sqrtxz1*rln2
       - 88*ln([1-z])*x^2*[1+z]^-1*sqrtxz1^-1*rln2
       + 80*ln([1-z])*x^2*sqrtxz1^-1*rln2
       + 8*ln([1-z])*x^2*z*[1+z]^-1*sqrtxz1^-1*rln2
       - 8*ln([1-z])*x^2*z*sqrtxz1^-1*rln2
       + 32*ln([1-z])*x^3*[1+z]^-1*sqrtxz1^-1*rln2
       - 32*ln([1-z])*x^3*sqrtxz1^-1*rln2
       + 8*ln([1-z])*ln(1 + sqrtxz1 - z)*z^-1*[1+z]^-1*sqrtxz1
       - 8*ln([1-z])*ln(1 + sqrtxz1 - z)*z^-1*sqrtxz1
       + 24*ln([1-z])*ln(1 + sqrtxz1 - z)*[1+z]^-1*sqrtxz1^-1
       - 16*ln([1-z])*ln(1 + sqrtxz1 - z)*sqrtxz1^-1
       - 8*ln([1-z])*ln(1 + sqrtxz1 - z)*z*[1+z]^-1*sqrtxz1^-1
       + 8*ln([1-z])*ln(1 + sqrtxz1 - z)*z*sqrtxz1^-1
       - 16*ln([1-z])*ln(1 + sqrtxz1 - z)*x*z^-1*[1+z]^-1*sqrtxz1
       + 16*ln([1-z])*ln(1 + sqrtxz1 - z)*x*z^-1*sqrtxz1
       - 80*ln([1-z])*ln(1 + sqrtxz1 - z)*x*[1+z]^-1*sqrtxz1^-1
       + 64*ln([1-z])*ln(1 + sqrtxz1 - z)*x*sqrtxz1^-1
       + 16*ln([1-z])*ln(1 + sqrtxz1 - z)*x*z*[1+z]^-1*sqrtxz1^-1
       - 16*ln([1-z])*ln(1 + sqrtxz1 - z)*x*z*sqrtxz1^-1
       + 8*ln([1-z])*ln(1 + sqrtxz1 - z)*x^2*z^-1*[1+z]^-1*sqrtxz1
       - 8*ln([1-z])*ln(1 + sqrtxz1 - z)*x^2*z^-1*sqrtxz1
       + 88*ln([1-z])*ln(1 + sqrtxz1 - z)*x^2*[1+z]^-1*sqrtxz1^-1
       - 80*ln([1-z])*ln(1 + sqrtxz1 - z)*x^2*sqrtxz1^-1
       - 8*ln([1-z])*ln(1 + sqrtxz1 - z)*x^2*z*[1+z]^-1*sqrtxz1^-1
       + 8*ln([1-z])*ln(1 + sqrtxz1 - z)*x^2*z*sqrtxz1^-1
       - 32*ln([1-z])*ln(1 + sqrtxz1 - z)*x^3*[1+z]^-1*sqrtxz1^-1
       + 32*ln([1-z])*ln(1 + sqrtxz1 - z)*x^3*sqrtxz1^-1
       + 4*ln([1-z])*ln(1 - 2*z + z^2 + 4*x*z)*z^-1*[1+z]^-1*sqrtxz1
       - 4*ln([1-z])*ln(1 - 2*z + z^2 + 4*x*z)*z^-1*sqrtxz1
       + 12*ln([1-z])*ln(1 - 2*z + z^2 + 4*x*z)*[1+z]^-1*sqrtxz1^-1
       - 8*ln([1-z])*ln(1 - 2*z + z^2 + 4*x*z)*sqrtxz1^-1
       - 4*ln([1-z])*ln(1 - 2*z + z^2 + 4*x*z)*z*[1+z]^-1*sqrtxz1^-1
       + 4*ln([1-z])*ln(1 - 2*z + z^2 + 4*x*z)*z*sqrtxz1^-1
       - 8*ln([1-z])*ln(1 - 2*z + z^2 + 4*x*z)*x*z^-1*[1+z]^-1*sqrtxz1
       + 8*ln([1-z])*ln(1 - 2*z + z^2 + 4*x*z)*x*z^-1*sqrtxz1
       - 40*ln([1-z])*ln(1 - 2*z + z^2 + 4*x*z)*x*[1+z]^-1*sqrtxz1^-1
       + 32*ln([1-z])*ln(1 - 2*z + z^2 + 4*x*z)*x*sqrtxz1^-1
       + 8*ln([1-z])*ln(1 - 2*z + z^2 + 4*x*z)*x*z*[1+z]^-1*sqrtxz1^-1
       - 8*ln([1-z])*ln(1 - 2*z + z^2 + 4*x*z)*x*z*sqrtxz1^-1
       + 4*ln([1-z])*ln(1 - 2*z + z^2 + 4*x*z)*x^2*z^-1*[1+z]^-1*sqrtxz1
       - 4*ln([1-z])*ln(1 - 2*z + z^2 + 4*x*z)*x^2*z^-1*sqrtxz1
       + 44*ln([1-z])*ln(1 - 2*z + z^2 + 4*x*z)*x^2*[1+z]^-1*sqrtxz1^-1
       - 40*ln([1-z])*ln(1 - 2*z + z^2 + 4*x*z)*x^2*sqrtxz1^-1
       - 4*ln([1-z])*ln(1 - 2*z + z^2 + 4*x*z)*x^2*z*[1+z]^-1*sqrtxz1^-1
       + 4*ln([1-z])*ln(1 - 2*z + z^2 + 4*x*z)*x^2*z*sqrtxz1^-1
       - 16*ln([1-z])*ln(1 - 2*z + z^2 + 4*x*z)*x^3*[1+z]^-1*sqrtxz1^-1
       + 16*ln([1-z])*ln(1 - 2*z + z^2 + 4*x*z)*x^3*sqrtxz1^-1
       - ln([1-z])^2
       - 4*ln([1-z])^2*z^-1*[1+z]^-1*sqrtxz1
       + ln([1-z])^2*z^-1
       + 4*ln([1-z])^2*z^-1*sqrtxz1
       - 12*ln([1-z])^2*[1+z]^-1*sqrtxz1^-1
       + 8*ln([1-z])^2*sqrtxz1^-1
       + 4*ln([1-z])^2*z*[1+z]^-1*sqrtxz1^-1
       - 4*ln([1-z])^2*z*sqrtxz1^-1
       + 1/2*ln([1-z])^2*z
       + 8*ln([1-z])^2*x*z^-1*[1+z]^-1*sqrtxz1
       - 2*ln([1-z])^2*x*z^-1
       - 8*ln([1-z])^2*x*z^-1*sqrtxz1
       + 40*ln([1-z])^2*x*[1+z]^-1*sqrtxz1^-1
       - 32*ln([1-z])^2*x*sqrtxz1^-1
       + 2*ln([1-z])^2*x
       - 8*ln([1-z])^2*x*z*[1+z]^-1*sqrtxz1^-1
       + 8*ln([1-z])^2*x*z*sqrtxz1^-1
       - ln([1-z])^2*x*z
       - 4*ln([1-z])^2*x^2*z^-1*[1+z]^-1*sqrtxz1
       + 4*ln([1-z])^2*x^2*z^-1*sqrtxz1
       - 44*ln([1-z])^2*x^2*[1+z]^-1*sqrtxz1^-1
       + 40*ln([1-z])^2*x^2*sqrtxz1^-1
       + 4*ln([1-z])^2*x^2*z*[1+z]^-1*sqrtxz1^-1
       - 4*ln([1-z])^2*x^2*z*sqrtxz1^-1
       + 16*ln([1-z])^2*x^3*[1+z]^-1*sqrtxz1^-1
       - 16*ln([1-z])^2*x^3*sqrtxz1^-1
       + 8*ln(sqrtxz3)*ArcTan(sqrtxz3)*z*sqrtxz3
       - 4*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z^-1*[1+z]^-1
       - 4*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z^-1*[1+z]^-1*sqrtxz1
       + 2*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z^-1
       + 4*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z^-1*sqrtxz1
       - 12*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1+z]^-1*sqrtxz1^-1
       + 8*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*sqrtxz1^-1
       + 4*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z*[1+z]^-1*sqrtxz1^-1
       - 4*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z*sqrtxz1^-1
       + 8*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x*z^-1*[1+z]^-1
       + 8*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x*z^-1*[1+z]^-1*sqrtxz1
       - 4*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x*z^-1
       - 8*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x*z^-1*sqrtxz1
       + 40*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x*[1+z]^-1*sqrtxz1^-1
       - 32*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x*sqrtxz1^-1
       - 8*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x*z*[1+z]^-1*sqrtxz1^-1
       + 8*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x*z*sqrtxz1^-1
       - 4*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x^2*z^-1*[1+z]^-1
       - 4*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x^2*z^-1*[1+z]^-1*sqrtxz1
       + 4*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x^2*z^-1
       + 4*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x^2*z^-1*sqrtxz1
       - 44*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x^2*[1+z]^-1*sqrtxz1^-1
       + 40*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x^2*sqrtxz1^-1
       + 4*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x^2*z*[1+z]^-1*sqrtxz1^-1
       - 4*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x^2*z*sqrtxz1^-1
       + 16*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x^3*[1+z]^-1*sqrtxz1^-1
       - 16*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x^3*sqrtxz1^-1
       + 4*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z^-1*[1+z]^-1
       - 4*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z^-1*[1+z]^-1*sqrtxz1
       - 2*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z^-1
       + 4*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z^-1*sqrtxz1
       - 12*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1+z]^-1*sqrtxz1^-1
       + 8*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*sqrtxz1^-1
       + 4*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z*[1+z]^-1*sqrtxz1^-1
       - 4*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z*sqrtxz1^-1
       - 8*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x*z^-1*[1+z]^-1
       + 8*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x*z^-1*[1+z]^-1*sqrtxz1
       + 4*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x*z^-1
       - 8*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x*z^-1*sqrtxz1
       + 40*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x*[1+z]^-1*sqrtxz1^-1
       - 32*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x*sqrtxz1^-1
       - 8*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x*z*[1+z]^-1*sqrtxz1^-1
       + 8*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x*z*sqrtxz1^-1
       + 4*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x^2*z^-1*[1+z]^-1
       - 4*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x^2*z^-1*[1+z]^-1*sqrtxz1
       - 4*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x^2*z^-1
       + 4*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x^2*z^-1*sqrtxz1
       - 44*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x^2*[1+z]^-1*sqrtxz1^-1
       + 40*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x^2*sqrtxz1^-1
       + 4*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x^2*z*[1+z]^-1*sqrtxz1^-1
       - 4*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x^2*z*sqrtxz1^-1
       + 16*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x^3*[1+z]^-1*sqrtxz1^-1
       - 16*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x^3*sqrtxz1^-1
       + 4*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)
       + 4*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*z^-1*[1+z]^-1
       + 4*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*z^-1*[1+z]^-1*sqrtxz1
       - 4*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*z^-1
       - 4*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*z^-1*sqrtxz1
       + 12*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*[1+z]^-1*sqrtxz1^-1
       - 8*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*sqrtxz1^-1
       - 4*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*z*[1+z]^-1*sqrtxz1^-1
       + 4*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*z*sqrtxz1^-1
       - 2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*z
       - 8*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*x*z^-1*[1+z]^-1
       - 8*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*x*z^-1*[1+z]^-1*sqrtxz1
       + 8*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*x*z^-1
       + 8*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*x*z^-1*sqrtxz1
       - 40*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*x*[1+z]^-1*sqrtxz1^-1
       + 32*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*x*sqrtxz1^-1
       - 8*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*x
       + 8*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*x*z*[1+z]^-1*sqrtxz1^-1
       - 8*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*x*z*sqrtxz1^-1
       + 4*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*x*z
       + 4*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*x^2*z^-1*[1+z]^-1
       + 4*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*x^2*z^-1*[1+z]^-1*sqrtxz1
       - 4*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*x^2*z^-1
       - 4*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*x^2*z^-1*sqrtxz1
       + 44*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*x^2*[1+z]^-1*sqrtxz1^-1
       - 40*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*x^2*sqrtxz1^-1
       + 4*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*x^2
       - 4*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*x^2*z*[1+z]^-1*sqrtxz1^-1
       + 4*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*x^2*z*sqrtxz1^-1
       - 16*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*x^3*[1+z]^-1*sqrtxz1^-1
       + 16*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*x^3*sqrtxz1^-1
       - 4*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)
       - 4*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*z^-1*[1+z]^-1
       + 4*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*z^-1*[1+z]^-1*sqrtxz1
       + 4*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*z^-1
       - 4*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*z^-1*sqrtxz1
       + 12*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*[1+z]^-1*sqrtxz1^-1
       - 8*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*sqrtxz1^-1
       - 4*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*z*[1+z]^-1*sqrtxz1^-1
       + 4*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*z*sqrtxz1^-1
       + 2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*z
       + 8*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*x*z^-1*[1+z]^-1
       - 8*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*x*z^-1*[1+z]^-1*sqrtxz1
       - 8*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*x*z^-1
       + 8*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*x*z^-1*sqrtxz1
       - 40*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*x*[1+z]^-1*sqrtxz1^-1
       + 32*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*x*sqrtxz1^-1
       + 8*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*x
       + 8*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*x*z*[1+z]^-1*sqrtxz1^-1
       - 8*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*x*z*sqrtxz1^-1
       - 4*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*x*z
       - 4*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*x^2*z^-1*[1+z]^-1
       + 4*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*x^2*z^-1*[1+z]^-1*sqrtxz1
       + 4*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*x^2*z^-1
       - 4*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*x^2*z^-1*sqrtxz1
       + 44*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*x^2*[1+z]^-1*sqrtxz1^-1
       - 40*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*x^2*sqrtxz1^-1
       - 4*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*x^2
       - 4*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*x^2*z*[1+z]^-1*sqrtxz1^-1
       + 4*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*x^2*z*sqrtxz1^-1
       - 16*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*x^3*[1+z]^-1*sqrtxz1^-1
       + 16*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*x^3*sqrtxz1^-1
       - 8*Li2(sqrtxz1*[z-1]^-1)*z^-1*[1+z]^-1*sqrtxz1
       + 8*Li2(sqrtxz1*[z-1]^-1)*z^-1*sqrtxz1
       - 24*Li2(sqrtxz1*[z-1]^-1)*[1+z]^-1*sqrtxz1^-1
       + 16*Li2(sqrtxz1*[z-1]^-1)*sqrtxz1^-1
       + 8*Li2(sqrtxz1*[z-1]^-1)*z*[1+z]^-1*sqrtxz1^-1
       - 8*Li2(sqrtxz1*[z-1]^-1)*z*sqrtxz1^-1
       + 16*Li2(sqrtxz1*[z-1]^-1)*x*z^-1*[1+z]^-1*sqrtxz1
       - 16*Li2(sqrtxz1*[z-1]^-1)*x*z^-1*sqrtxz1
       + 80*Li2(sqrtxz1*[z-1]^-1)*x*[1+z]^-1*sqrtxz1^-1
       - 64*Li2(sqrtxz1*[z-1]^-1)*x*sqrtxz1^-1
       - 16*Li2(sqrtxz1*[z-1]^-1)*x*z*[1+z]^-1*sqrtxz1^-1
       + 16*Li2(sqrtxz1*[z-1]^-1)*x*z*sqrtxz1^-1
       - 8*Li2(sqrtxz1*[z-1]^-1)*x^2*z^-1*[1+z]^-1*sqrtxz1
       + 8*Li2(sqrtxz1*[z-1]^-1)*x^2*z^-1*sqrtxz1
       - 88*Li2(sqrtxz1*[z-1]^-1)*x^2*[1+z]^-1*sqrtxz1^-1
       + 80*Li2(sqrtxz1*[z-1]^-1)*x^2*sqrtxz1^-1
       + 8*Li2(sqrtxz1*[z-1]^-1)*x^2*z*[1+z]^-1*sqrtxz1^-1
       - 8*Li2(sqrtxz1*[z-1]^-1)*x^2*z*sqrtxz1^-1
       + 32*Li2(sqrtxz1*[z-1]^-1)*x^3*[1+z]^-1*sqrtxz1^-1
       - 32*Li2(sqrtxz1*[z-1]^-1)*x^3*sqrtxz1^-1
       - 8*Li2([1-z]*sqrtxz1^-1)*z^-1*[1+z]^-1*sqrtxz1
       + 8*Li2([1-z]*sqrtxz1^-1)*z^-1*sqrtxz1
       - 24*Li2([1-z]*sqrtxz1^-1)*[1+z]^-1*sqrtxz1^-1
       + 16*Li2([1-z]*sqrtxz1^-1)*sqrtxz1^-1
       + 8*Li2([1-z]*sqrtxz1^-1)*z*[1+z]^-1*sqrtxz1^-1
       - 8*Li2([1-z]*sqrtxz1^-1)*z*sqrtxz1^-1
       + 16*Li2([1-z]*sqrtxz1^-1)*x*z^-1*[1+z]^-1*sqrtxz1
       - 16*Li2([1-z]*sqrtxz1^-1)*x*z^-1*sqrtxz1
       + 80*Li2([1-z]*sqrtxz1^-1)*x*[1+z]^-1*sqrtxz1^-1
       - 64*Li2([1-z]*sqrtxz1^-1)*x*sqrtxz1^-1
       - 16*Li2([1-z]*sqrtxz1^-1)*x*z*[1+z]^-1*sqrtxz1^-1
       + 16*Li2([1-z]*sqrtxz1^-1)*x*z*sqrtxz1^-1
       - 8*Li2([1-z]*sqrtxz1^-1)*x^2*z^-1*[1+z]^-1*sqrtxz1
       + 8*Li2([1-z]*sqrtxz1^-1)*x^2*z^-1*sqrtxz1
       - 88*Li2([1-z]*sqrtxz1^-1)*x^2*[1+z]^-1*sqrtxz1^-1
       + 80*Li2([1-z]*sqrtxz1^-1)*x^2*sqrtxz1^-1
       + 8*Li2([1-z]*sqrtxz1^-1)*x^2*z*[1+z]^-1*sqrtxz1^-1
       - 8*Li2([1-z]*sqrtxz1^-1)*x^2*z*sqrtxz1^-1
       + 32*Li2([1-z]*sqrtxz1^-1)*x^3*[1+z]^-1*sqrtxz1^-1
       - 32*Li2([1-z]*sqrtxz1^-1)*x^3*sqrtxz1^-1
       + 2*Li2( - x*z^-1)
       + 2*Li2( - x*z^-1)*x^-2*[1+x]^-1*z^-1*[1-z]^-1
       - 2*Li2( - x*z^-1)*x^-2*[1+x]^-1*z^-1
       - Li2( - x*z^-1)*x^-2*[1+x]^-1
       - Li2( - x*z^-1)*x^-2*[1+x]^-1*z
       - 2*Li2( - x*z^-1)*x^-2*z^-1*[1-z]^-1
       + 2*Li2( - x*z^-1)*x^-2*z^-1
       + Li2( - x*z^-1)*x^-2
       + Li2( - x*z^-1)*x^-2*z
       + 6*Li2( - x*z^-1)*x^-1*[1+x]^-1*z^-1*[1-z]^-1
       - 6*Li2( - x*z^-1)*x^-1*[1+x]^-1*z^-1
       - 3*Li2( - x*z^-1)*x^-1*[1+x]^-1
       - 3*Li2( - x*z^-1)*x^-1*[1+x]^-1*z
       - 4*Li2( - x*z^-1)*x^-1*z^-1*[1-z]^-1
       + 4*Li2( - x*z^-1)*x^-1*z^-1
       + 2*Li2( - x*z^-1)*x^-1
       + 2*Li2( - x*z^-1)*x^-1*z
       + 6*Li2( - x*z^-1)*[1+x]^-1*z^-1*[1-z]^-1
       - 6*Li2( - x*z^-1)*[1+x]^-1*z^-1
       - 4*Li2( - x*z^-1)*[1+x]^-1
       - 2*Li2( - x*z^-1)*[1+x]^-1*z
       + 2*Li2( - x*z^-1)*z^-1*[1+z]^-1
       - 2*Li2( - x*z^-1)*z^-1
       - 2*Li2( - x*z^-1)*z
       + 2*Li2( - x*z^-1)*x*[1+x]^-1*z^-1*[1-z]^-1
       - 2*Li2( - x*z^-1)*x*[1+x]^-1*z^-1
       - 2*Li2( - x*z^-1)*x*[1+x]^-1
       + 4*Li2( - x*z^-1)*x*z^-1*[1-z]^-1
       - 4*Li2( - x*z^-1)*x*z^-1*[1+z]^-1
       - 4*Li2( - x*z^-1)*x
       + 2*Li2( - x*z^-1)*x^2*z^-1*[1-z]^-1
       + 2*Li2( - x*z^-1)*x^2*z^-1*[1+z]^-1
       - 4*Li2( - x*z^-1)*x^2*z^-1
       - 2*Li2( - x)
       - 4*Li2( - x)*x^-2*[1+x]^-1*z^-1*[1-z]^-1
       + 4*Li2( - x)*x^-2*[1+x]^-1*z^-1
       + 2*Li2( - x)*x^-2*[1+x]^-1
       + 2*Li2( - x)*x^-2*[1+x]^-1*z
       + 4*Li2( - x)*x^-2*z^-1*[1-z]^-1
       - 4*Li2( - x)*x^-2*z^-1
       - 2*Li2( - x)*x^-2
       - 2*Li2( - x)*x^-2*z
       - 12*Li2( - x)*x^-1*[1+x]^-1*z^-1*[1-z]^-1
       + 12*Li2( - x)*x^-1*[1+x]^-1*z^-1
       + 6*Li2( - x)*x^-1*[1+x]^-1
       + 6*Li2( - x)*x^-1*[1+x]^-1*z
       + 8*Li2( - x)*x^-1*z^-1*[1-z]^-1
       - 8*Li2( - x)*x^-1*z^-1
       - 4*Li2( - x)*x^-1
       - 4*Li2( - x)*x^-1*z
       - 12*Li2( - x)*[1+x]^-1*z^-1*[1-z]^-1
       + 12*Li2( - x)*[1+x]^-1*z^-1
       + 8*Li2( - x)*[1+x]^-1
       + 4*Li2( - x)*[1+x]^-1*z
       + 2*Li2( - x)*z
       - 4*Li2( - x)*x*[1+x]^-1*z^-1*[1-z]^-1
       + 4*Li2( - x)*x*[1+x]^-1*z^-1
       + 4*Li2( - x)*x*[1+x]^-1
       - 8*Li2( - x)*x*z^-1*[1-z]^-1
       + 8*Li2( - x)*x*z^-1
       + 4*Li2( - x)*x
       + 4*Li2( - x)*x*z
       - 4*Li2( - x)*x^2*z^-1*[1-z]^-1
       + 4*Li2( - x)*x^2*z^-1
       + 4*Li2( - x)*x^2
       + 2*Li2( - x*z)*x^-2*[1+x]^-1*z^-1*[1-z]^-1
       - 2*Li2( - x*z)*x^-2*[1+x]^-1*z^-1
       - Li2( - x*z)*x^-2*[1+x]^-1
       - Li2( - x*z)*x^-2*[1+x]^-1*z
       - 2*Li2( - x*z)*x^-2*z^-1*[1-z]^-1
       + 2*Li2( - x*z)*x^-2*z^-1
       + Li2( - x*z)*x^-2
       + Li2( - x*z)*x^-2*z
       + 6*Li2( - x*z)*x^-1*[1+x]^-1*z^-1*[1-z]^-1
       - 6*Li2( - x*z)*x^-1*[1+x]^-1*z^-1
       - 3*Li2( - x*z)*x^-1*[1+x]^-1
       - 3*Li2( - x*z)*x^-1*[1+x]^-1*z
       - 4*Li2( - x*z)*x^-1*z^-1*[1-z]^-1
       + 4*Li2( - x*z)*x^-1*z^-1
       + 2*Li2( - x*z)*x^-1
       + 2*Li2( - x*z)*x^-1*z
       + 6*Li2( - x*z)*[1+x]^-1*z^-1*[1-z]^-1
       - 6*Li2( - x*z)*[1+x]^-1*z^-1
       - 4*Li2( - x*z)*[1+x]^-1
       - 2*Li2( - x*z)*[1+x]^-1*z
       - 2*Li2( - x*z)*z^-1*[1+z]^-1
       + 2*Li2( - x*z)*z^-1
       + 2*Li2( - x*z)*x*[1+x]^-1*z^-1*[1-z]^-1
       - 2*Li2( - x*z)*x*[1+x]^-1*z^-1
       - 2*Li2( - x*z)*x*[1+x]^-1
       + 4*Li2( - x*z)*x*z^-1*[1-z]^-1
       + 4*Li2( - x*z)*x*z^-1*[1+z]^-1
       - 8*Li2( - x*z)*x*z^-1
       - 4*Li2( - x*z)*x*z
       + 2*Li2( - x*z)*x^2*z^-1*[1-z]^-1
       - 2*Li2( - x*z)*x^2*z^-1*[1+z]^-1
       - 4*Li2( - x*z)*x^2
       - 4*Li2( - 1/(1 + sqrtxz1 - z) + 1/(1 + sqrtxz1 - z)*sqrtxz1 + 1/(1 + sqrtxz1 - z)*z)*z^-1*[1+z]^-1*sqrtxz1
       + 4*Li2( - 1/(1 + sqrtxz1 - z) + 1/(1 + sqrtxz1 - z)*sqrtxz1 + 1/(1 + sqrtxz1 - z)*z)*z^-1*sqrtxz1
       - 12*Li2( - 1/(1 + sqrtxz1 - z) + 1/(1 + sqrtxz1 - z)*sqrtxz1 + 1/(1 + sqrtxz1 - z)*z)*[1+z]^-1*sqrtxz1^-1
       + 8*Li2( - 1/(1 + sqrtxz1 - z) + 1/(1 + sqrtxz1 - z)*sqrtxz1 + 1/(1 + sqrtxz1 - z)*z)*sqrtxz1^-1
       + 4*Li2( - 1/(1 + sqrtxz1 - z) + 1/(1 + sqrtxz1 - z)*sqrtxz1 + 1/(1 + sqrtxz1 - z)*z)*z*[1+z]^-1*sqrtxz1^-1
       - 4*Li2( - 1/(1 + sqrtxz1 - z) + 1/(1 + sqrtxz1 - z)*sqrtxz1 + 1/(1 + sqrtxz1 - z)*z)*z*sqrtxz1^-1
       + 8*Li2( - 1/(1 + sqrtxz1 - z) + 1/(1 + sqrtxz1 - z)*sqrtxz1 + 1/(1 + sqrtxz1 - z)*z)*x*z^-1*[1+z]^-1*sqrtxz1
       - 8*Li2( - 1/(1 + sqrtxz1 - z) + 1/(1 + sqrtxz1 - z)*sqrtxz1 + 1/(1 + sqrtxz1 - z)*z)*x*z^-1*sqrtxz1
       + 40*Li2( - 1/(1 + sqrtxz1 - z) + 1/(1 + sqrtxz1 - z)*sqrtxz1 + 1/(1 + sqrtxz1 - z)*z)*x*[1+z]^-1*sqrtxz1^-1
       - 32*Li2( - 1/(1 + sqrtxz1 - z) + 1/(1 + sqrtxz1 - z)*sqrtxz1 + 1/(1 + sqrtxz1 - z)*z)*x*sqrtxz1^-1
       - 8*Li2( - 1/(1 + sqrtxz1 - z) + 1/(1 + sqrtxz1 - z)*sqrtxz1 + 1/(1 + sqrtxz1 - z)*z)*x*z*[1+z]^-1*sqrtxz1^-1
       + 8*Li2( - 1/(1 + sqrtxz1 - z) + 1/(1 + sqrtxz1 - z)*sqrtxz1 + 1/(1 + sqrtxz1 - z)*z)*x*z*sqrtxz1^-1
       - 4*Li2( - 1/(1 + sqrtxz1 - z) + 1/(1 + sqrtxz1 - z)*sqrtxz1 + 1/(1 + sqrtxz1 - z)*z)*x^2*z^-1*[1+z]^-1*sqrtxz1
       + 4*Li2( - 1/(1 + sqrtxz1 - z) + 1/(1 + sqrtxz1 - z)*sqrtxz1 + 1/(1 + sqrtxz1 - z)*z)*x^2*z^-1*sqrtxz1
       - 44*Li2( - 1/(1 + sqrtxz1 - z) + 1/(1 + sqrtxz1 - z)*sqrtxz1 + 1/(1 + sqrtxz1 - z)*z)*x^2*[1+z]^-1*sqrtxz1^-1
       + 40*Li2( - 1/(1 + sqrtxz1 - z) + 1/(1 + sqrtxz1 - z)*sqrtxz1 + 1/(1 + sqrtxz1 - z)*z)*x^2*sqrtxz1^-1
       + 4*Li2( - 1/(1 + sqrtxz1 - z) + 1/(1 + sqrtxz1 - z)*sqrtxz1 + 1/(1 + sqrtxz1 - z)*z)*x^2*z*[1+z]^-1*sqrtxz1^-1
       - 4*Li2( - 1/(1 + sqrtxz1 - z) + 1/(1 + sqrtxz1 - z)*sqrtxz1 + 1/(1 + sqrtxz1 - z)*z)*x^2*z*sqrtxz1^-1
       + 16*Li2( - 1/(1 + sqrtxz1 - z) + 1/(1 + sqrtxz1 - z)*sqrtxz1 + 1/(1 + sqrtxz1 - z)*z)*x^3*[1+z]^-1*sqrtxz1^-1
       - 16*Li2( - 1/(1 + sqrtxz1 - z) + 1/(1 + sqrtxz1 - z)*sqrtxz1 + 1/(1 + sqrtxz1 - z)*z)*x^3*sqrtxz1^-1
       + 3*Li2(x)
       + 4*Li2(x)*x^-2*[1+x]^-1*z^-1*[1-z]^-1
       - 4*Li2(x)*x^-2*[1+x]^-1*z^-1
       - 2*Li2(x)*x^-2*[1+x]^-1
       - 2*Li2(x)*x^-2*[1+x]^-1*z
       - 4*Li2(x)*x^-2*z^-1*[1-z]^-1
       + 4*Li2(x)*x^-2*z^-1
       + 2*Li2(x)*x^-2
       + 2*Li2(x)*x^-2*z
       + 12*Li2(x)*x^-1*[1+x]^-1*z^-1*[1-z]^-1
       - 12*Li2(x)*x^-1*[1+x]^-1*z^-1
       - 6*Li2(x)*x^-1*[1+x]^-1
       - 6*Li2(x)*x^-1*[1+x]^-1*z
       - 8*Li2(x)*x^-1*z^-1*[1-z]^-1
       + 8*Li2(x)*x^-1*z^-1
       + 4*Li2(x)*x^-1
       + 4*Li2(x)*x^-1*z
       + 12*Li2(x)*[1+x]^-1*z^-1*[1-z]^-1
       - 12*Li2(x)*[1+x]^-1*z^-1
       - 8*Li2(x)*[1+x]^-1
       - 4*Li2(x)*[1+x]^-1*z
       + 1/2*Li2(x)*z^-1
       - 4*Li2(x)*[1-z]^-1
       + 4*Li2(x)*x*[1+x]^-1*z^-1*[1-z]^-1
       - 4*Li2(x)*x*[1+x]^-1*z^-1
       - 4*Li2(x)*x*[1+x]^-1
       + 8*Li2(x)*x*z^-1*[1-z]^-1
       - 9*Li2(x)*x*z^-1
       - 8*Li2(x)*x*[1-z]^-1
       + 2*Li2(x)*x
       - 2*Li2(x)*x*z
       + 4*Li2(x)*x^2*z^-1*[1-z]^-1
       - 4*Li2(x)*x^2*z^-1
       - 4*Li2(x)*x^2*[1-z]^-1
       - 2*Li2(x)*x^2
       + 2*Li2(z)
       - 4*Li2(z)*x
       + 4*InvTanInt( - sqrtxz3)*z*sqrtxz3
       + 8*InvTanInt(z*sqrtxz3)*z*sqrtxz3
       - 4*InvTanInt(sqrtxz3)*z*sqrtxz3
       )

       + NC * ( 65/4
       - 157/8*z^-1
       - 4*z^-1*rln2^2
       - 4*z^-1*sqrtxz1*rln2
       + 4*rln2^2
       + 43/8*z
       - 4*z*rln2^2
       + 41/2*x*z^-1
       + 8*x*z^-1*rln2^2
       + 8*x*z^-1*sqrtxz1*rln2
       - 17*x
       - 8*x*rln2^2
       - 11/2*x*z
       + 8*x*z*rln2^2
       + 8*x^2*rln2^2
       - 1/3*pi^2*x^-2*[1+x]^-1*z^-1
       - 1/3*pi^2*x^-2*[1+x]^-1*z
       + 1/3*pi^2*x^-2*z^-1
       + 1/3*pi^2*x^-2*z
       - pi^2*x^-1*[1+x]^-1*z^-1
       - pi^2*x^-1*[1+x]^-1*z
       + 2/3*pi^2*x^-1*z^-1
       + 2/3*pi^2*x^-1*z
       - 2/3*pi^2*[1+x]^-1*z^-1
       - 2/3*pi^2*[1+x]^-1
       - 2/3*pi^2*[1+x]^-1*z
       + 7/12*pi^2*z^-1
       + 1/6*pi^2
       + 1/3*pi^2*z
       - 2/3*pi^2*x*[1+x]^-1
       - 7/6*pi^2*x*z^-1
       + pi^2*x
       - pi^2*x*z
       - 2/3*pi^2*x^2
       + 6*ln(1 + sqrtxz1 - z)*z^-1*rln2
       + 4*ln(1 + sqrtxz1 - z)*z^-1*sqrtxz1
       - 8*ln(1 + sqrtxz1 - z)*rln2
       + 6*ln(1 + sqrtxz1 - z)*z*rln2
       - 12*ln(1 + sqrtxz1 - z)*x*z^-1*rln2
       - 8*ln(1 + sqrtxz1 - z)*x*z^-1*sqrtxz1
       + 16*ln(1 + sqrtxz1 - z)*x*rln2
       - 12*ln(1 + sqrtxz1 - z)*x*z*rln2
       - 12*ln(1 + sqrtxz1 - z)*x^2*rln2
       + 4*ln(1 + sqrtxz1 - z)^2
       - 2*ln(1 + sqrtxz1 - z)^2*z^-1
       - 2*ln(1 + sqrtxz1 - z)^2*z
       + 4*ln(1 + sqrtxz1 - z)^2*x*z^-1
       - 8*ln(1 + sqrtxz1 - z)^2*x
       + 4*ln(1 + sqrtxz1 - z)^2*x*z
       + 4*ln(1 + sqrtxz1 - z)^2*x^2
       - 2*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*z^-1
       - 2*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*z
       + 4*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*x*z^-1
       + 4*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*x*z
       + 4*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*x^2
       + 2*ln(1 + sqrtxz1 + z)*z^-1*rln2
       + 2*ln(1 + sqrtxz1 + z)*z*rln2
       - 4*ln(1 + sqrtxz1 + z)*x*z^-1*rln2
       - 4*ln(1 + sqrtxz1 + z)*x*z*rln2
       - 4*ln(1 + sqrtxz1 + z)*x^2*rln2
       + 2*ln(z*sqrtxz3)*ArcTan(z*sqrtxz3)*x^-1*sqrtxz3
       + 10*ln(z*sqrtxz3)*ArcTan(z*sqrtxz3)*z^-1*sqrtxz3
       + 18*ln(z*sqrtxz3)*ArcTan(z*sqrtxz3)*z*sqrtxz3
       + 10*ln(z*sqrtxz3)*ArcTan(z*sqrtxz3)*x*sqrtxz3
       + 13*ln(x)
       - 4*ln(x)*x^-2*[1+x]^-1*z^-1
       - 4*ln(x)*x^-2*[1+x]^-1*z
       + 4*ln(x)*x^-2*z^-1
       + 4*ln(x)*x^-2*z
       - 12*ln(x)*x^-1*[1+x]^-1*z^-1
       - 12*ln(x)*x^-1*[1+x]^-1*z
       + 8*ln(x)*x^-1*z^-1
       + 8*ln(x)*x^-1*z
       - 8*ln(x)*[1+x]^-1*z^-1
       - 8*ln(x)*[1+x]^-1
       - 8*ln(x)*[1+x]^-1*z
       - 15/2*ln(x)*z^-1
       - 4*ln(x)*z^-1*rln2
       - 2*ln(x)*z^-1*sqrtxz1
       + 6*ln(x)*rln2
       + 7/2*ln(x)*z
       - 4*ln(x)*z*rln2
       - 8*ln(x)*x*[1+x]^-1
       - ln(x)*x*z^-1
       + 8*ln(x)*x*z^-1*rln2
       + 4*ln(x)*x*z^-1*sqrtxz1
       + ln(x)*x
       - 12*ln(x)*x*rln2
       + 8*ln(x)*x*z*rln2
       + 8*ln(x)*x^2*rln2
       - 4*ln(x)*ln(1 + sqrtxz1 - z)
       + 2*ln(x)*ln(1 + sqrtxz1 - z)*z^-1
       + 2*ln(x)*ln(1 + sqrtxz1 - z)*z
       - 4*ln(x)*ln(1 + sqrtxz1 - z)*x*z^-1
       + 8*ln(x)*ln(1 + sqrtxz1 - z)*x
       - 4*ln(x)*ln(1 + sqrtxz1 - z)*x*z
       - 4*ln(x)*ln(1 + sqrtxz1 - z)*x^2
       - 2*ln(x)*ln(1 + sqrtxz1 + z)
       + 2*ln(x)*ln(1 + sqrtxz1 + z)*z^-1
       + 2*ln(x)*ln(1 + sqrtxz1 + z)*z
       - 4*ln(x)*ln(1 + sqrtxz1 + z)*x*z^-1
       + 4*ln(x)*ln(1 + sqrtxz1 + z)*x
       - 4*ln(x)*ln(1 + sqrtxz1 + z)*x*z
       - 4*ln(x)*ln(1 + sqrtxz1 + z)*x^2
       - 2*ln(x)*ln(1 + x*z^-1)
       + ln(x)*ln(1 + x*z^-1)*x^-2*[1+x]^-1*z^-1
       + ln(x)*ln(1 + x*z^-1)*x^-2*[1+x]^-1*z
       - ln(x)*ln(1 + x*z^-1)*x^-2*z^-1
       - ln(x)*ln(1 + x*z^-1)*x^-2*z
       + 3*ln(x)*ln(1 + x*z^-1)*x^-1*[1+x]^-1*z^-1
       + 3*ln(x)*ln(1 + x*z^-1)*x^-1*[1+x]^-1*z
       - 2*ln(x)*ln(1 + x*z^-1)*x^-1*z^-1
       - 2*ln(x)*ln(1 + x*z^-1)*x^-1*z
       + 2*ln(x)*ln(1 + x*z^-1)*[1+x]^-1*z^-1
       + 2*ln(x)*ln(1 + x*z^-1)*[1+x]^-1
       + 2*ln(x)*ln(1 + x*z^-1)*[1+x]^-1*z
       + ln(x)*ln(1 + x*z^-1)*z^-1
       + ln(x)*ln(1 + x*z^-1)*z
       + 2*ln(x)*ln(1 + x*z^-1)*x*[1+x]^-1
       + 2*ln(x)*ln(1 + x*z^-1)*x*z^-1
       + 2*ln(x)*ln(1 + x*z^-1)*x*z
       + 2*ln(x)*ln(1 + x*z^-1)*x^2
       - 2*ln(x)*ln(1 + x*z)
       + ln(x)*ln(1 + x*z)*x^-2*[1+x]^-1*z^-1
       + ln(x)*ln(1 + x*z)*x^-2*[1+x]^-1*z
       - ln(x)*ln(1 + x*z)*x^-2*z^-1
       - ln(x)*ln(1 + x*z)*x^-2*z
       + 3*ln(x)*ln(1 + x*z)*x^-1*[1+x]^-1*z^-1
       + 3*ln(x)*ln(1 + x*z)*x^-1*[1+x]^-1*z
       - 2*ln(x)*ln(1 + x*z)*x^-1*z^-1
       - 2*ln(x)*ln(1 + x*z)*x^-1*z
       + 2*ln(x)*ln(1 + x*z)*[1+x]^-1*z^-1
       + 2*ln(x)*ln(1 + x*z)*[1+x]^-1
       + 2*ln(x)*ln(1 + x*z)*[1+x]^-1*z
       + 2*ln(x)*ln(1 + x*z)*x*[1+x]^-1
       + 4*ln(x)*ln(1 + x*z)*x*z^-1
       + 4*ln(x)*ln(1 + x*z)*x*z
       + 4*ln(x)*ln(1 + x*z)*x^2
       + ln(x)*ln(z + x)*z^-1
       + ln(x)*ln(z + x)*z
       - 2*ln(x)*ln(z + x)*x*z^-1
       - 2*ln(x)*ln(z + x)*x*z
       - 2*ln(x)*ln(z + x)*x^2
       - 2*ln(x)^2
       + 3/2*ln(x)^2*x^-2*[1+x]^-1*z^-1
       + 3/2*ln(x)^2*x^-2*[1+x]^-1*z
       - 3/2*ln(x)^2*x^-2*z^-1
       - 3/2*ln(x)^2*x^-2*z
       + 9/2*ln(x)^2*x^-1*[1+x]^-1*z^-1
       + 9/2*ln(x)^2*x^-1*[1+x]^-1*z
       - 3*ln(x)^2*x^-1*z^-1
       - 3*ln(x)^2*x^-1*z
       + 3*ln(x)^2*[1+x]^-1*z^-1
       + 3*ln(x)^2*[1+x]^-1
       + 3*ln(x)^2*[1+x]^-1*z
       - ln(x)^2*z^-1
       + 3*ln(x)^2*x*[1+x]^-1
       + 2*ln(x)^2*x*z^-1
       - 2*ln(x)^2*x
       + 2*ln(x)^2*x*z
       + ln(x)^2*x^2
       - 4*ln(x)*ln([1-x])
       + 2*ln(x)*ln([1-x])*x^-2*[1+x]^-1*z^-1
       + 2*ln(x)*ln([1-x])*x^-2*[1+x]^-1*z
       - 2*ln(x)*ln([1-x])*x^-2*z^-1
       - 2*ln(x)*ln([1-x])*x^-2*z
       + 6*ln(x)*ln([1-x])*x^-1*[1+x]^-1*z^-1
       + 6*ln(x)*ln([1-x])*x^-1*[1+x]^-1*z
       - 4*ln(x)*ln([1-x])*x^-1*z^-1
       - 4*ln(x)*ln([1-x])*x^-1*z
       + 4*ln(x)*ln([1-x])*[1+x]^-1*z^-1
       + 4*ln(x)*ln([1-x])*[1+x]^-1
       + 4*ln(x)*ln([1-x])*[1+x]^-1*z
       + ln(x)*ln([1-x])*z^-1
       + 4*ln(x)*ln([1-x])*x*[1+x]^-1
       - 2*ln(x)*ln([1-x])*x*z^-1
       + 2*ln(x)*ln([1-x])*x^2
       + 4*ln(x)*ln([1+x])
       - 2*ln(x)*ln([1+x])*x^-2*[1+x]^-1*z^-1
       - 2*ln(x)*ln([1+x])*x^-2*[1+x]^-1*z
       + 2*ln(x)*ln([1+x])*x^-2*z^-1
       + 2*ln(x)*ln([1+x])*x^-2*z
       - 6*ln(x)*ln([1+x])*x^-1*[1+x]^-1*z^-1
       - 6*ln(x)*ln([1+x])*x^-1*[1+x]^-1*z
       + 4*ln(x)*ln([1+x])*x^-1*z^-1
       + 4*ln(x)*ln([1+x])*x^-1*z
       - 4*ln(x)*ln([1+x])*[1+x]^-1*z^-1
       - 4*ln(x)*ln([1+x])*[1+x]^-1
       - 4*ln(x)*ln([1+x])*[1+x]^-1*z
       - 2*ln(x)*ln([1+x])*z^-1
       - 2*ln(x)*ln([1+x])*z
       - 4*ln(x)*ln([1+x])*x*[1+x]^-1
       - 4*ln(x)*ln([1+x])*x*z^-1
       - 4*ln(x)*ln([1+x])*x*z
       - 4*ln(x)*ln([1+x])*x^2
       - ln(x)*ln(z)
       - ln(x)*ln(z)*x^-2*[1+x]^-1*z^-1
       - ln(x)*ln(z)*x^-2*[1+x]^-1*z
       + ln(x)*ln(z)*x^-2*z^-1
       + ln(x)*ln(z)*x^-2*z
       - 3*ln(x)*ln(z)*x^-1*[1+x]^-1*z^-1
       - 3*ln(x)*ln(z)*x^-1*[1+x]^-1*z
       + 2*ln(x)*ln(z)*x^-1*z^-1
       + 2*ln(x)*ln(z)*x^-1*z
       - 2*ln(x)*ln(z)*[1+x]^-1*z^-1
       - 2*ln(x)*ln(z)*[1+x]^-1
       - 2*ln(x)*ln(z)*[1+x]^-1*z
       - 7/2*ln(x)*ln(z)*z^-1
       - 2*ln(x)*ln(z)*z
       - 2*ln(x)*ln(z)*x*[1+x]^-1
       - ln(x)*ln(z)*x*z^-1
       - 6*ln(x)*ln(z)*x
       + 2*ln(x)*ln(z)*x*z
       + 2*ln(x)*ln(z)*x^2
       - ln(x)*ln([1-z])
       + 3/2*ln(x)*ln([1-z])*z^-1
       - 3*ln(x)*ln([1-z])*x*z^-1
       + 2*ln(x)*ln([1-z])*x
       - 2*ln(x)*ln([1-z])*x*z
       - 11/2*ln([1-x])
       + 25/4*ln([1-x])*z^-1
       - 3/4*ln([1-x])*z
       - 9*ln([1-x])*x*z^-1
       + 8*ln([1-x])*x
       + 3/2*ln([1-x])*x*z
       + ln([1-x])^2
       - ln([1-x])^2*z^-1
       - 1/2*ln([1-x])^2*z
       + 2*ln([1-x])^2*x*z^-1
       - 2*ln([1-x])^2*x
       + ln([1-x])^2*x*z
       + 2*ln([1-x])*ln(z)
       - 4*ln([1-x])*ln(z)*x
       + 2*ln([1-x])*ln([1-z])
       - 2*ln([1-x])*ln([1-z])*z^-1
       - ln([1-x])*ln([1-z])*z
       + 4*ln([1-x])*ln([1-z])*x*z^-1
       - 4*ln([1-x])*ln([1-z])*x
       + 2*ln([1-x])*ln([1-z])*x*z
       - 17/2*ln(z)
       + 5/4*ln(z)*z^-1
       - 4*ln(z)*z^-1*rln2
       - 2*ln(z)*z^-1*sqrtxz1
       + 2*ln(z)*rln2
       + 9/4*ln(z)*z
       - 4*ln(z)*z*rln2
       - 6*ln(z)*x*z^-1
       + 8*ln(z)*x*z^-1*rln2
       + 4*ln(z)*x*z^-1*sqrtxz1
       + 7*ln(z)*x
       - 4*ln(z)*x*rln2
       - 5/2*ln(z)*x*z
       + 8*ln(z)*x*z*rln2
       + 8*ln(z)*x^2*rln2
       - 2*ln(z)*ln(1 + sqrtxz1 - z)
       + 2*ln(z)*ln(1 + sqrtxz1 - z)*z^-1
       + 2*ln(z)*ln(1 + sqrtxz1 - z)*z
       - 4*ln(z)*ln(1 + sqrtxz1 - z)*x*z^-1
       + 4*ln(z)*ln(1 + sqrtxz1 - z)*x
       - 4*ln(z)*ln(1 + sqrtxz1 - z)*x*z
       - 4*ln(z)*ln(1 + sqrtxz1 - z)*x^2
       + 2*ln(z)*ln(1 + sqrtxz1 + z)*z^-1
       + 2*ln(z)*ln(1 + sqrtxz1 + z)*z
       - 4*ln(z)*ln(1 + sqrtxz1 + z)*x*z^-1
       - 4*ln(z)*ln(1 + sqrtxz1 + z)*x*z
       - 4*ln(z)*ln(1 + sqrtxz1 + z)*x^2
       + 2*ln(z)*ln(1 + x*z^-1)
       - ln(z)*ln(1 + x*z^-1)*x^-2*[1+x]^-1*z^-1
       - ln(z)*ln(1 + x*z^-1)*x^-2*[1+x]^-1*z
       + ln(z)*ln(1 + x*z^-1)*x^-2*z^-1
       + ln(z)*ln(1 + x*z^-1)*x^-2*z
       - 3*ln(z)*ln(1 + x*z^-1)*x^-1*[1+x]^-1*z^-1
       - 3*ln(z)*ln(1 + x*z^-1)*x^-1*[1+x]^-1*z
       + 2*ln(z)*ln(1 + x*z^-1)*x^-1*z^-1
       + 2*ln(z)*ln(1 + x*z^-1)*x^-1*z
       - 2*ln(z)*ln(1 + x*z^-1)*[1+x]^-1*z^-1
       - 2*ln(z)*ln(1 + x*z^-1)*[1+x]^-1
       - 2*ln(z)*ln(1 + x*z^-1)*[1+x]^-1*z
       - ln(z)*ln(1 + x*z^-1)*z^-1
       - ln(z)*ln(1 + x*z^-1)*z
       - 2*ln(z)*ln(1 + x*z^-1)*x*[1+x]^-1
       - 2*ln(z)*ln(1 + x*z^-1)*x*z^-1
       - 2*ln(z)*ln(1 + x*z^-1)*x*z
       - 2*ln(z)*ln(1 + x*z^-1)*x^2
       - 2*ln(z)*ln(1 + x*z)
       + ln(z)*ln(1 + x*z)*x^-2*[1+x]^-1*z^-1
       + ln(z)*ln(1 + x*z)*x^-2*[1+x]^-1*z
       - ln(z)*ln(1 + x*z)*x^-2*z^-1
       - ln(z)*ln(1 + x*z)*x^-2*z
       + 3*ln(z)*ln(1 + x*z)*x^-1*[1+x]^-1*z^-1
       + 3*ln(z)*ln(1 + x*z)*x^-1*[1+x]^-1*z
       - 2*ln(z)*ln(1 + x*z)*x^-1*z^-1
       - 2*ln(z)*ln(1 + x*z)*x^-1*z
       + 2*ln(z)*ln(1 + x*z)*[1+x]^-1*z^-1
       + 2*ln(z)*ln(1 + x*z)*[1+x]^-1
       + 2*ln(z)*ln(1 + x*z)*[1+x]^-1*z
       + 2*ln(z)*ln(1 + x*z)*x*[1+x]^-1
       + 4*ln(z)*ln(1 + x*z)*x*z^-1
       + 4*ln(z)*ln(1 + x*z)*x*z
       + 4*ln(z)*ln(1 + x*z)*x^2
       - ln(z)*ln(z + x)*z^-1
       - ln(z)*ln(z + x)*z
       + 2*ln(z)*ln(z + x)*x*z^-1
       + 2*ln(z)*ln(z + x)*x*z
       + 2*ln(z)*ln(z + x)*x^2
       + ln(z)^2
       - 1/2*ln(z)^2*x^-2*[1+x]^-1*z^-1
       - 1/2*ln(z)^2*x^-2*[1+x]^-1*z
       + 1/2*ln(z)^2*x^-2*z^-1
       + 1/2*ln(z)^2*x^-2*z
       - 3/2*ln(z)^2*x^-1*[1+x]^-1*z^-1
       - 3/2*ln(z)^2*x^-1*[1+x]^-1*z
       + ln(z)^2*x^-1*z^-1
       + ln(z)^2*x^-1*z
       - ln(z)^2*[1+x]^-1*z^-1
       - ln(z)^2*[1+x]^-1
       - ln(z)^2*[1+x]^-1*z
       + ln(z)^2*z^-1
       + 1/2*ln(z)^2*z
       - ln(z)^2*x*[1+x]^-1
       - 4*ln(z)^2*x*z^-1
       - 3*ln(z)^2*x*z
       - 2*ln(z)^2*x^2
       - 11/2*ln([1-z])
       + 25/4*ln([1-z])*z^-1
       - 3/4*ln([1-z])*z
       - 9*ln([1-z])*x*z^-1
       + 8*ln([1-z])*x
       + 3/2*ln([1-z])*x*z
       + ln([1-z])^2
       - ln([1-z])^2*z^-1
       - 1/2*ln([1-z])^2*z
       + 2*ln([1-z])^2*x*z^-1
       - 2*ln([1-z])^2*x
       + ln([1-z])^2*x*z
       - 2*ln(sqrtxz3)*ArcTan(sqrtxz3)*x^-1*sqrtxz3
       - 10*ln(sqrtxz3)*ArcTan(sqrtxz3)*z^-1*sqrtxz3
       - 18*ln(sqrtxz3)*ArcTan(sqrtxz3)*z*sqrtxz3
       - 10*ln(sqrtxz3)*ArcTan(sqrtxz3)*x*sqrtxz3
       - 2*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)
       + 4*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x
       + 2*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)
       - 4*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x
       - 2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)
       + 2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*z^-1
       + 2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*z
       - 4*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*x*z^-1
       + 4*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*x
       - 4*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*x*z
       - 4*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*x^2
       + 2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)
       - 2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*z^-1
       - 2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*z
       + 4*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*x*z^-1
       - 4*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*x
       + 4*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*x*z
       + 4*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*x^2
       - 2*Li2( - x*z^-1)
       + Li2( - x*z^-1)*x^-2*[1+x]^-1*z^-1
       + Li2( - x*z^-1)*x^-2*[1+x]^-1*z
       - Li2( - x*z^-1)*x^-2*z^-1
       - Li2( - x*z^-1)*x^-2*z
       + 3*Li2( - x*z^-1)*x^-1*[1+x]^-1*z^-1
       + 3*Li2( - x*z^-1)*x^-1*[1+x]^-1*z
       - 2*Li2( - x*z^-1)*x^-1*z^-1
       - 2*Li2( - x*z^-1)*x^-1*z
       + 2*Li2( - x*z^-1)*[1+x]^-1*z^-1
       + 2*Li2( - x*z^-1)*[1+x]^-1
       + 2*Li2( - x*z^-1)*[1+x]^-1*z
       + 2*Li2( - x*z^-1)*z^-1
       + 2*Li2( - x*z^-1)*z
       + 2*Li2( - x*z^-1)*x*[1+x]^-1
       + 4*Li2( - x)
       - 2*Li2( - x)*x^-2*[1+x]^-1*z^-1
       - 2*Li2( - x)*x^-2*[1+x]^-1*z
       + 2*Li2( - x)*x^-2*z^-1
       + 2*Li2( - x)*x^-2*z
       - 6*Li2( - x)*x^-1*[1+x]^-1*z^-1
       - 6*Li2( - x)*x^-1*[1+x]^-1*z
       + 4*Li2( - x)*x^-1*z^-1
       + 4*Li2( - x)*x^-1*z
       - 4*Li2( - x)*[1+x]^-1*z^-1
       - 4*Li2( - x)*[1+x]^-1
       - 4*Li2( - x)*[1+x]^-1*z
       - 2*Li2( - x)*z^-1
       - 2*Li2( - x)*z
       - 4*Li2( - x)*x*[1+x]^-1
       - 4*Li2( - x)*x*z^-1
       - 4*Li2( - x)*x*z
       - 4*Li2( - x)*x^2
       - 2*Li2( - x*z)
       + Li2( - x*z)*x^-2*[1+x]^-1*z^-1
       + Li2( - x*z)*x^-2*[1+x]^-1*z
       - Li2( - x*z)*x^-2*z^-1
       - Li2( - x*z)*x^-2*z
       + 3*Li2( - x*z)*x^-1*[1+x]^-1*z^-1
       + 3*Li2( - x*z)*x^-1*[1+x]^-1*z
       - 2*Li2( - x*z)*x^-1*z^-1
       - 2*Li2( - x*z)*x^-1*z
       + 2*Li2( - x*z)*[1+x]^-1*z^-1
       + 2*Li2( - x*z)*[1+x]^-1
       + 2*Li2( - x*z)*[1+x]^-1*z
       + 2*Li2( - x*z)*x*[1+x]^-1
       + 4*Li2( - x*z)*x*z^-1
       + 4*Li2( - x*z)*x*z
       + 4*Li2( - x*z)*x^2
       - 3*Li2(x)
       + 2*Li2(x)*x^-2*[1+x]^-1*z^-1
       + 2*Li2(x)*x^-2*[1+x]^-1*z
       - 2*Li2(x)*x^-2*z^-1
       - 2*Li2(x)*x^-2*z
       + 6*Li2(x)*x^-1*[1+x]^-1*z^-1
       + 6*Li2(x)*x^-1*[1+x]^-1*z
       - 4*Li2(x)*x^-1*z^-1
       - 4*Li2(x)*x^-1*z
       + 4*Li2(x)*[1+x]^-1*z^-1
       + 4*Li2(x)*[1+x]^-1
       + 4*Li2(x)*[1+x]^-1*z
       - 1/2*Li2(x)*z^-1
       + 4*Li2(x)*x*[1+x]^-1
       + Li2(x)*x*z^-1
       - 2*Li2(x)*x
       + 2*Li2(x)*x*z
       + 2*Li2(x)*x^2
       - 2*Li2(z)
       + 4*Li2(z)*x
       - InvTanInt( - sqrtxz3)*x^-1*sqrtxz3
       - 5*InvTanInt( - sqrtxz3)*z^-1*sqrtxz3
       - 9*InvTanInt( - sqrtxz3)*z*sqrtxz3
       - 5*InvTanInt( - sqrtxz3)*x*sqrtxz3
       - 2*InvTanInt(z*sqrtxz3)*x^-1*sqrtxz3
       - 10*InvTanInt(z*sqrtxz3)*z^-1*sqrtxz3
       - 18*InvTanInt(z*sqrtxz3)*z*sqrtxz3
       - 10*InvTanInt(z*sqrtxz3)*x*sqrtxz3
       + InvTanInt(sqrtxz3)*x^-1*sqrtxz3
       + 5*InvTanInt(sqrtxz3)*z^-1*sqrtxz3
       + 9*InvTanInt(sqrtxz3)*z*sqrtxz3
       + 5*InvTanInt(sqrtxz3)*x*sqrtxz3
       )

       + LMUF*NC^-1 * (  - 1
       + 3/2*z^-1
       - 3/2*z
       - 3/2*x*z^-1
       + x
       + 2*x*z
       + 1/2*ln(x)*z^-1
       - 1/2*ln(x)*z
       - ln(x)*x*z^-1
       - ln(x)*x*z
       + ln([1-x])
       - ln([1-x])*z^-1
       - 1/2*ln([1-x])*z
       + 2*ln([1-x])*x*z^-1
       - 2*ln([1-x])*x
       + ln([1-x])*x*z
       + ln(z)
       - ln(z)*z^-1
       - 1/2*ln(z)*z
       + 2*ln(z)*x*z^-1
       - 2*ln(z)*x
       + ln(z)*x*z
       + ln([1-z])
       - ln([1-z])*z^-1
       - 1/2*ln([1-z])*z
       + 2*ln([1-z])*x*z^-1
       - 2*ln([1-z])*x
       + ln([1-z])*x*z
       )

       + LMUF*NC * ( 1
       - 3/2*z^-1
       + 3/2*z
       + 3/2*x*z^-1
       - x
       - 2*x*z
       - 1/2*ln(x)*z^-1
       + 1/2*ln(x)*z
       + ln(x)*x*z^-1
       + ln(x)*x*z
       - ln([1-x])
       + ln([1-x])*z^-1
       + 1/2*ln([1-x])*z
       - 2*ln([1-x])*x*z^-1
       + 2*ln([1-x])*x
       - ln([1-x])*x*z
       - ln(z)
       + ln(z)*z^-1
       + 1/2*ln(z)*z
       - 2*ln(z)*x*z^-1
       + 2*ln(z)*x
       - ln(z)*x*z
       - ln([1-z])
       + ln([1-z])*z^-1
       + 1/2*ln([1-z])*z
       - 2*ln([1-z])*x*z^-1
       + 2*ln([1-z])*x
       - ln([1-z])*x*z
       )

       + LMUA*NC^-1 * (  - 9/2
       + 19/4*z^-1
       + 3/4*z
       - 15/2*x*z^-1
       + 7*x
       - 1/2*x*z
       - ln(x)
       + ln(x)*z^-1
       + 1/2*ln(x)*z
       - 2*ln(x)*x*z^-1
       + 2*ln(x)*x
       - ln(x)*x*z
       + ln([1-x])
       - ln([1-x])*z^-1
       - 1/2*ln([1-x])*z
       + 2*ln([1-x])*x*z^-1
       - 2*ln([1-x])*x
       + ln([1-x])*x*z
       + ln(z)
       + ln(z)*z^-1
       + 1/2*ln(z)*z
       - 2*ln(z)*x*z^-1
       - 2*ln(z)*x
       - ln(z)*x*z
       + ln([1-z])
       - ln([1-z])*z^-1
       - 1/2*ln([1-z])*z
       + 2*ln([1-z])*x*z^-1
       - 2*ln([1-z])*x
       + ln([1-z])*x*z
       )

       + LMUA*NC * ( 9/2
       - 19/4*z^-1
       - 3/4*z
       + 15/2*x*z^-1
       - 7*x
       + 1/2*x*z
       + ln(x)
       - ln(x)*z^-1
       - 1/2*ln(x)*z
       + 2*ln(x)*x*z^-1
       - 2*ln(x)*x
       + ln(x)*x*z
       - ln([1-x])
       + ln([1-x])*z^-1
       + 1/2*ln([1-x])*z
       - 2*ln([1-x])*x*z^-1
       + 2*ln([1-x])*x
       - ln([1-x])*x*z
       - ln(z)
       - ln(z)*z^-1
       - 1/2*ln(z)*z
       + 2*ln(z)*x*z^-1
       + 2*ln(z)*x
       + ln(z)*x*z
       - ln([1-z])
       + ln([1-z])*z^-1
       + 1/2*ln([1-z])*z
       - 2*ln([1-z])*x*z^-1
       + 2*ln([1-z])*x
       - ln([1-z])*x*z
       )

       + LMUA*LMUF*NC^-1 * (  - 1
       + z^-1
       + 1/2*z
       - 2*x*z^-1
       + 2*x
       - x*z
       )

       + LMUA*LMUF*NC * ( 1
       - z^-1
       - 1/2*z
       + 2*x*z^-1
       - 2*x
       + x*z ) );

id DC2G2Q      = (  + NC^-1 * (  - 2
       + 2*z^-1*[1-z]^-1
       - 47/16*z^-1
       - 4*[1-z]^-1
       - 1/8*poly2^-1
       - 1/4*z*[1-x-z]^-1
       + 51/16*z
       + 1/2*z^2*[1-x-z]^-1
       - 4*x*z^-1*[1-z]^-1
       + 11/4*x*z^-1
       + 3*x*[1-z]^-1
       + 1/8*x*poly2^-1
       + 37/8*x
       - 11/4*x*z
       + 2*x^2*z^-1*[1-z]^-1
       - 2*x^2*z^-1
       - 2*x^2*[1-z]^-1
       + 1/8*x^2*poly2^-1
       - 1/8*x^3*poly2^-1
       - 13/12*pi^2*z^-1
       - 3/4*pi^2*[1-z]^-1
       + 5/6*pi^2
       - 1/6*pi^2*z
       + 3/2*pi^2*x*z^-1
       + 3/2*pi^2*x*[1-z]^-1
       - 5/6*pi^2*x
       - 1/6*pi^2*x*z
       - 5/6*pi^2*x^2*z^-1
       - 1/2*pi^2*x^2*[1-z]^-1
       + 41/16*ln(x)
       - 8*ln(x)*z^-1*[1-z]^-1
       + 33/8*ln(x)*z^-1
       + 7/2*ln(x)*[1-z]^-1
       - 1/2*ln(x)*[1-x-z]^-2
       + ln(x)*[1-x-z]^-1
       - 3/16*ln(x)*poly2^-2
       - 3/8*ln(x)*poly2^-1
       - 1/4*ln(x)*z*[1-x-z]^-2
       + 2*ln(x)*z*[1-x-z]^-1
       - 3/4*ln(x)*z
       + 3/4*ln(x)*z^2*[1-x-z]^-2
       - 5/2*ln(x)*z^2*[1-x-z]^-1
       + 16*ln(x)*x*z^-1*[1-z]^-1
       - 17/4*ln(x)*x*z^-1
       - 4*ln(x)*x*[1-z]^-1
       + 3/2*ln(x)*x*[1-x-z]^-2
       - 3*ln(x)*x*[1-x-z]^-1
       - 1/2*ln(x)*x*[x-z]^-1
       + 3/16*ln(x)*x*poly2^-2
       + 1/2*ln(x)*x*poly2^-1
       - 251/16*ln(x)*x
       - 8*ln(x)*x^2*z^-1*[1-z]^-1
       + 8*ln(x)*x^2*z^-1
       + 8*ln(x)*x^2*[1-z]^-1
       - 3/2*ln(x)*x^2*[1-x-z]^-2
       + 3/2*ln(x)*x^2*[1-x-z]^-1
       + ln(x)*x^2*[x-z]^-1
       + 3/8*ln(x)*x^2*poly2^-2
       - 1/8*ln(x)*x^2*poly2^-1
       + ln(x)*x^2
       + 1/2*ln(x)*x^3*[1-x-z]^-2
       - 3/8*ln(x)*x^3*poly2^-2
       - 3/16*ln(x)*x^4*poly2^-2
       + 3/16*ln(x)*x^5*poly2^-2
       - ln(x)*ln(1 - sqrtxz2 + x)*z^-1*sqrtxz2^-1
       - ln(x)*ln(1 - sqrtxz2 + x)*[1-z]^-1*sqrtxz2^-1
       - 3/32*ln(x)*ln(1 - sqrtxz2 + x)*sqrtxz2^-1*poly2^-2
       - 5/32*ln(x)*ln(1 - sqrtxz2 + x)*sqrtxz2^-1*poly2^-1
       + 7/4*ln(x)*ln(1 - sqrtxz2 + x)*sqrtxz2^-1
       - 3*ln(x)*ln(1 - sqrtxz2 + x)*x*z^-1*sqrtxz2^-1
       + 3*ln(x)*ln(1 - sqrtxz2 + x)*x*[1-z]^-1*sqrtxz2^-1
       + 29/16*ln(x)*ln(1 - sqrtxz2 + x)*x*sqrtxz2^-1
       - 29/8*ln(x)*ln(1 - sqrtxz2 + x)*x*z*sqrtxz2^-1
       - 3*ln(x)*ln(1 - sqrtxz2 + x)*x^2*z^-1*sqrtxz2^-1
       - 3*ln(x)*ln(1 - sqrtxz2 + x)*x^2*[1-z]^-1*sqrtxz2^-1
       + 9/32*ln(x)*ln(1 - sqrtxz2 + x)*x^2*sqrtxz2^-1*poly2^-2
       + 3/16*ln(x)*ln(1 - sqrtxz2 + x)*x^2*sqrtxz2^-1*poly2^-1
       + 111/16*ln(x)*ln(1 - sqrtxz2 + x)*x^2*sqrtxz2^-1
       - ln(x)*ln(1 - sqrtxz2 + x)*x^3*z^-1*sqrtxz2^-1
       + ln(x)*ln(1 - sqrtxz2 + x)*x^3*[1-z]^-1*sqrtxz2^-1
       - 9/32*ln(x)*ln(1 - sqrtxz2 + x)*x^4*sqrtxz2^-1*poly2^-2
       - 1/32*ln(x)*ln(1 - sqrtxz2 + x)*x^4*sqrtxz2^-1*poly2^-1
       + 3/32*ln(x)*ln(1 - sqrtxz2 + x)*x^6*sqrtxz2^-1*poly2^-2
       + ln(x)*ln(1 + sqrtxz2 + x)*z^-1*sqrtxz2^-1
       + ln(x)*ln(1 + sqrtxz2 + x)*[1-z]^-1*sqrtxz2^-1
       + 3/32*ln(x)*ln(1 + sqrtxz2 + x)*sqrtxz2^-1*poly2^-2
       + 5/32*ln(x)*ln(1 + sqrtxz2 + x)*sqrtxz2^-1*poly2^-1
       - 7/4*ln(x)*ln(1 + sqrtxz2 + x)*sqrtxz2^-1
       + 3*ln(x)*ln(1 + sqrtxz2 + x)*x*z^-1*sqrtxz2^-1
       - 3*ln(x)*ln(1 + sqrtxz2 + x)*x*[1-z]^-1*sqrtxz2^-1
       - 29/16*ln(x)*ln(1 + sqrtxz2 + x)*x*sqrtxz2^-1
       + 29/8*ln(x)*ln(1 + sqrtxz2 + x)*x*z*sqrtxz2^-1
       + 3*ln(x)*ln(1 + sqrtxz2 + x)*x^2*z^-1*sqrtxz2^-1
       + 3*ln(x)*ln(1 + sqrtxz2 + x)*x^2*[1-z]^-1*sqrtxz2^-1
       - 9/32*ln(x)*ln(1 + sqrtxz2 + x)*x^2*sqrtxz2^-1*poly2^-2
       - 3/16*ln(x)*ln(1 + sqrtxz2 + x)*x^2*sqrtxz2^-1*poly2^-1
       - 111/16*ln(x)*ln(1 + sqrtxz2 + x)*x^2*sqrtxz2^-1
       + ln(x)*ln(1 + sqrtxz2 + x)*x^3*z^-1*sqrtxz2^-1
       - ln(x)*ln(1 + sqrtxz2 + x)*x^3*[1-z]^-1*sqrtxz2^-1
       + 9/32*ln(x)*ln(1 + sqrtxz2 + x)*x^4*sqrtxz2^-1*poly2^-2
       + 1/32*ln(x)*ln(1 + sqrtxz2 + x)*x^4*sqrtxz2^-1*poly2^-1
       - 3/32*ln(x)*ln(1 + sqrtxz2 + x)*x^6*sqrtxz2^-1*poly2^-2
       + 2*ln(x)*ln(1 + x)
       - 2*ln(x)*ln(1 + x)*z^-1
       - ln(x)*ln(1 + x)*z
       - 4*ln(x)*ln(1 + x)*x*z^-1
       + 4*ln(x)*ln(1 + x)*x
       - 2*ln(x)*ln(1 + x)*x*z
       - 2*ln(x)*ln(1 + x)*x^2*z^-1
       + 2*ln(x)*ln(1 + x)*x^2
       - 19/4*ln(x)^2
       + 3*ln(x)^2*z^-1*[1-z]^-1
       + 11/4*ln(x)^2*z^-1
       + 3*ln(x)^2*[1-z]^-1
       - 6*ln(x)^2*x*z^-1*[1-z]^-1
       - 7/2*ln(x)^2*x*z^-1
       - 6*ln(x)^2*x*[1-z]^-1
       + 13/2*ln(x)^2*x
       + 2*ln(x)^2*x*z
       + 3*ln(x)^2*x^2*z^-1*[1-z]^-1
       + ln(x)^2*x^2*z^-1
       + ln(x)^2*x^2*[1-z]^-1
       - 3/2*ln(x)^2*x^2
       + 13/2*ln(x)*ln([1-x])
       - 8*ln(x)*ln([1-x])*z^-1
       - 19/2*ln(x)*ln([1-x])*[1-z]^-1
       + ln(x)*ln([1-x])*z
       + 16*ln(x)*ln([1-x])*x*z^-1
       + 19*ln(x)*ln([1-x])*x*[1-z]^-1
       - 13*ln(x)*ln([1-x])*x
       - 2*ln(x)*ln([1-x])*x*z
       - 5*ln(x)*ln([1-x])*x^2*z^-1
       - 7*ln(x)*ln([1-x])*x^2*[1-z]^-1
       + ln(x)*ln([1-x])*x^2
       + 5/2*ln(x)*ln(z)
       - 2*ln(x)*ln(z)*z^-1*[1-z]^-1
       - 9/4*ln(x)*ln(z)*z^-1
       - 4*ln(x)*ln(z)*[1-z]^-1
       + 1/2*ln(x)*ln(z)*z
       + 4*ln(x)*ln(z)*x*z^-1*[1-z]^-1
       + 17/2*ln(x)*ln(z)*x*z^-1
       + 8*ln(x)*ln(z)*x*[1-z]^-1
       - 9*ln(x)*ln(z)*x
       - 2*ln(x)*ln(z)*x^2*z^-1*[1-z]^-1
       - 2*ln(x)*ln(z)*x^2*[1-z]^-1
       + 11/2*ln(x)*ln([1-z])
       - 7*ln(x)*ln([1-z])*z^-1
       - 6*ln(x)*ln([1-z])*[1-z]^-1
       + 14*ln(x)*ln([1-z])*x*z^-1
       + 12*ln(x)*ln([1-z])*x*[1-z]^-1
       - 10*ln(x)*ln([1-z])*x
       - ln(x)*ln([1-z])*x*z
       - 5*ln(x)*ln([1-z])*x^2*z^-1
       - 3*ln(x)*ln([1-z])*x^2*[1-z]^-1
       - 2*ln([1-x])
       + 9/4*ln([1-x])*z^-1
       + 3*ln([1-x])*[1-z]^-1
       - 3/8*ln([1-x])*z
       - 13/2*ln([1-x])*x*z^-1
       - 8*ln([1-x])*x*[1-z]^-1
       + 19/2*ln([1-x])*x
       + 3/4*ln([1-x])*x*z
       - ln([1-x])*x^2
       - 9/4*ln([1-x])^2
       + 3*ln([1-x])^2*z^-1
       + 3*ln([1-x])^2*[1-z]^-1
       - 1/4*ln([1-x])^2*z
       - 6*ln([1-x])^2*x*z^-1
       - 6*ln([1-x])^2*x*[1-z]^-1
       + 9/2*ln([1-x])^2*x
       + 1/2*ln([1-x])^2*x*z
       + 2*ln([1-x])^2*x^2*z^-1
       + 2*ln([1-x])^2*x^2*[1-z]^-1
       - 2*ln([1-x])*ln(z)
       + 7/2*ln([1-x])*ln(z)*z^-1
       + 4*ln([1-x])*ln(z)*[1-z]^-1
       - 7*ln([1-x])*ln(z)*x*z^-1
       - 8*ln([1-x])*ln(z)*x*[1-z]^-1
       + 4*ln([1-x])*ln(z)*x
       + 2*ln([1-x])*ln(z)*x^2*z^-1
       + 3*ln([1-x])*ln(z)*x^2*[1-z]^-1
       - 4*ln([1-x])*ln([1-z])
       + 9/2*ln([1-x])*ln([1-z])*z^-1
       + 4*ln([1-x])*ln([1-z])*[1-z]^-1
       - 1/2*ln([1-x])*ln([1-z])*z
       - 9*ln([1-x])*ln([1-z])*x*z^-1
       - 8*ln([1-x])*ln([1-z])*x*[1-z]^-1
       + 8*ln([1-x])*ln([1-z])*x
       + ln([1-x])*ln([1-z])*x*z
       + 3*ln([1-x])*ln([1-z])*x^2*z^-1
       + 2*ln([1-x])*ln([1-z])*x^2*[1-z]^-1
       - 83/16*ln(z)
       + 4*ln(z)*z^-1*[1-z]^-1
       - 13/8*ln(z)*z^-1
       - 5/2*ln(z)*[1-z]^-1
       - 3/16*ln(z)*poly2^-2
       - 3/8*ln(z)*poly2^-1
       + 9/8*ln(z)*z
       - 8*ln(z)*x*z^-1*[1-z]^-1
       + 2*ln(z)*x*z^-1
       + 3*ln(z)*x*[1-z]^-1
       + 1/2*ln(z)*x*[x-z]^-1
       - 3/16*ln(z)*x*poly2^-2
       - 1/2*ln(z)*x*poly2^-1
       + 187/16*ln(z)*x
       - 5/4*ln(z)*x*z
       + 4*ln(z)*x^2*z^-1*[1-z]^-1
       - 4*ln(z)*x^2*z^-1
       - 4*ln(z)*x^2*[1-z]^-1
       - ln(z)*x^2*[x-z]^-1
       + 3/8*ln(z)*x^2*poly2^-2
       - 1/8*ln(z)*x^2*poly2^-1
       + 3/8*ln(z)*x^3*poly2^-2
       - 3/16*ln(z)*x^4*poly2^-2
       - 3/16*ln(z)*x^5*poly2^-2
       - 1/4*ln(z)^2
       + 5/4*ln(z)^2*z^-1
       + ln(z)^2*[1-z]^-1
       + 1/4*ln(z)^2*z
       - 5/2*ln(z)^2*x*z^-1
       - 2*ln(z)^2*x*[1-z]^-1
       + 1/2*ln(z)^2*x
       - 1/2*ln(z)^2*x*z
       + 1/2*ln(z)^2*x^2*z^-1
       + ln(z)^2*x^2*[1-z]^-1
       - 7/2*ln(z)*ln([1-z])
       + 3*ln(z)*ln([1-z])*z^-1
       + 4*ln(z)*ln([1-z])*[1-z]^-1
       - 6*ln(z)*ln([1-z])*x*z^-1
       - 8*ln(z)*ln([1-z])*x*[1-z]^-1
       + 7*ln(z)*ln([1-z])*x
       + ln(z)*ln([1-z])*x^2*z^-1
       + 3*ln(z)*ln([1-z])*x^2*[1-z]^-1
       + ln([1-z])
       + 7/4*ln([1-z])*z^-1
       + 2*ln([1-z])*[1-z]^-1
       + 1/2*ln([1-z])*[1-x-z]^-2
       - ln([1-z])*[1-x-z]^-1
       + 1/4*ln([1-z])*z*[1-x-z]^-2
       - 2*ln([1-z])*z*[1-x-z]^-1
       + 17/8*ln([1-z])*z
       - 3/4*ln([1-z])*z^2*[1-x-z]^-2
       + 5/2*ln([1-z])*z^2*[1-x-z]^-1
       - 9/2*ln([1-z])*x*z^-1
       - 6*ln([1-z])*x*[1-z]^-1
       - 3/2*ln([1-z])*x*[1-x-z]^-2
       + 3*ln([1-z])*x*[1-x-z]^-1
       + 5/2*ln([1-z])*x
       + 3/4*ln([1-z])*x*z
       + 3/2*ln([1-z])*x^2*[1-x-z]^-2
       - 3/2*ln([1-z])*x^2*[1-x-z]^-1
       - 1/2*ln([1-z])*x^3*[1-x-z]^-2
       - 2*ln([1-z])^2
       + 5/2*ln([1-z])^2*z^-1
       + 3/2*ln([1-z])^2*[1-z]^-1
       - 1/4*ln([1-z])^2*z
       - 5*ln([1-z])^2*x*z^-1
       - 3*ln([1-z])^2*x*[1-z]^-1
       + 4*ln([1-z])^2*x
       + 1/2*ln([1-z])^2*x*z
       + 2*ln([1-z])^2*x^2*z^-1
       + 1/2*ln([1-z])^2*x^2*[1-z]^-1
       - Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*z^-1*sqrtxz2^-1
       - Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*[1-z]^-1*sqrtxz2^-1
       - 3/32*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*sqrtxz2^-1*poly2^-2
       - 5/32*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*sqrtxz2^-1*poly2^-1
       + 7/4*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*sqrtxz2^-1
       - 3*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x*z^-1*sqrtxz2^-1
       + 3*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x*[1-z]^-1*sqrtxz2^-1
       + 29/16*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x*sqrtxz2^-1
       - 29/8*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x*z*sqrtxz2^-1
       - 3*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x^2*z^-1*sqrtxz2^-1
       - 3*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x^2*[1-z]^-1*sqrtxz2^-1
       + 9/32*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x^2*sqrtxz2^-1*poly2^-2
       + 3/16*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x^2*sqrtxz2^-1*poly2^-1
       + 111/16*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x^2*sqrtxz2^-1
       - Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x^3*z^-1*sqrtxz2^-1
       + Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x^3*[1-z]^-1*sqrtxz2^-1
       - 9/32*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x^4*sqrtxz2^-1*poly2^-2
       - 1/32*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x^4*sqrtxz2^-1*poly2^-1
       + 3/32*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x^6*sqrtxz2^-1*poly2^-2
       + Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*z^-1*sqrtxz2^-1
       + Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*[1-z]^-1*sqrtxz2^-1
       + 3/32*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*sqrtxz2^-1*poly2^-2
       + 5/32*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*sqrtxz2^-1*poly2^-1
       - 7/4*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*sqrtxz2^-1
       + 3*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x*z^-1*sqrtxz2^-1
       - 3*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x*[1-z]^-1*sqrtxz2^-1
       - 29/16*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x*sqrtxz2^-1
       + 29/8*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x*z*sqrtxz2^-1
       + 3*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x^2*z^-1*sqrtxz2^-1
       + 3*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x^2*[1-z]^-1*sqrtxz2^-1
       - 9/32*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x^2*sqrtxz2^-1*poly2^-2
       - 3/16*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x^2*sqrtxz2^-1*poly2^-1
       - 111/16*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x^2*sqrtxz2^-1
       + Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x^3*z^-1*sqrtxz2^-1
       - Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x^3*[1-z]^-1*sqrtxz2^-1
       + 9/32*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x^4*sqrtxz2^-1*poly2^-2
       + 1/32*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x^4*sqrtxz2^-1*poly2^-1
       - 3/32*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x^6*sqrtxz2^-1*poly2^-2
       + Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*z^-1*sqrtxz2^-1
       + Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*[1-z]^-1*sqrtxz2^-1
       + 3/32*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*sqrtxz2^-1*poly2^-2
       + 5/32*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*sqrtxz2^-1*poly2^-1
       - 7/4*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*sqrtxz2^-1
       + 3*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x*z^-1*sqrtxz2^-1
       - 3*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x*[1-z]^-1*sqrtxz2^-1
       - 29/16*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x*sqrtxz2^-1
       + 29/8*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x*z*sqrtxz2^-1
       + 3*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x^2*z^-1*sqrtxz2^-1
       + 3*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x^2*[1-z]^-1*sqrtxz2^-1
       - 9/32*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x^2*sqrtxz2^-1*poly2^-2
       - 3/16*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x^2*sqrtxz2^-1*poly2^-1
       - 111/16*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x^2*sqrtxz2^-1
       + Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x^3*z^-1*sqrtxz2^-1
       - Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x^3*[1-z]^-1*sqrtxz2^-1
       + 9/32*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x^4*sqrtxz2^-1*poly2^-2
       + 1/32*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x^4*sqrtxz2^-1*poly2^-1
       - 3/32*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x^6*sqrtxz2^-1*poly2^-2
       - Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*z^-1*sqrtxz2^-1
       - Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*[1-z]^-1*sqrtxz2^-1
       - 3/32*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*sqrtxz2^-1*poly2^-2
       - 5/32*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*sqrtxz2^-1*poly2^-1
       + 7/4*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*sqrtxz2^-1
       - 3*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x*z^-1*sqrtxz2^-1
       + 3*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x*[1-z]^-1*sqrtxz2^-1
       + 29/16*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x*sqrtxz2^-1
       - 29/8*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x*z*sqrtxz2^-1
       - 3*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x^2*z^-1*sqrtxz2^-1
       - 3*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x^2*[1-z]^-1*sqrtxz2^-1
       + 9/32*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x^2*sqrtxz2^-1*poly2^-2
       + 3/16*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x^2*sqrtxz2^-1*poly2^-1
       + 111/16*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x^2*sqrtxz2^-1
       - Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x^3*z^-1*sqrtxz2^-1
       + Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x^3*[1-z]^-1*sqrtxz2^-1
       - 9/32*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x^4*sqrtxz2^-1*poly2^-2
       - 1/32*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x^4*sqrtxz2^-1*poly2^-1
       + 3/32*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x^6*sqrtxz2^-1*poly2^-2
       + 1/2*Li2(1 - x*z^-1)
       - 2*Li2(1 - x*z^-1)*z^-1*[1-z]^-1
       + 3/2*Li2(1 - x*z^-1)*z^-1
       + 1/2*Li2(1 - x*z^-1)*[1-z]^-1
       + 4*Li2(1 - x*z^-1)*x*z^-1*[1-z]^-1
       - 3*Li2(1 - x*z^-1)*x*z^-1
       - Li2(1 - x*z^-1)*x*[1-z]^-1
       - Li2(1 - x*z^-1)*x
       - 2*Li2(1 - x*z^-1)*x^2*z^-1*[1-z]^-1
       + 2*Li2(1 - x*z^-1)*x^2*z^-1
       + 2*Li2( - x)
       - 2*Li2( - x)*z^-1
       - Li2( - x)*z
       - 4*Li2( - x)*x*z^-1
       + 4*Li2( - x)*x
       - 2*Li2( - x)*x*z
       - 2*Li2( - x)*x^2*z^-1
       + 2*Li2( - x)*x^2
       + Li2(x)
       - 1/2*Li2(x)*z^-1
       - 2*Li2(x)*[1-z]^-1
       + Li2(x)*z
       + Li2(x)*x*z^-1
       + 4*Li2(x)*x*[1-z]^-1
       - 3*Li2(x)*x
       - Li2(x)*x*z
       - 2*Li2(x)*x^2*[1-z]^-1
       + Li2(x)*x^2
       - 2*Li2(z)
       + Li2(z)*[1-z]^-1
       - 2*Li2(z)*x*[1-z]^-1
       + 4*Li2(z)*x
       - Li2(z)*x^2*z^-1
       + Li2(z)*x^2*[1-z]^-1
       )

       + NF * (  - 1/6*ln(x^-1*[1-x]*z*[1-z])*z^-1
       - 1/6*ln(x^-1*[1-x]*z*[1-z])*[1-z]^-1
       - 1/6*ln(x)*z^-1
       - 1/6*ln(x)*[1-z]^-1
       + 1/6*ln([1-x])*z^-1
       + 1/6*ln([1-x])*[1-z]^-1
       + 1/6*ln(z)*z^-1
       + 1/6*ln(z)*[1-z]^-1
       + 1/6*ln([1-z])*z^-1
       + 1/6*ln([1-z])*[1-z]^-1
       )

       + NC * ( 14
       - 2*z^-1*[1-z]^-1
       - 33/16*z^-1
       + 3*[1-z]^-1
       + 1/8*poly2^-1
       + 1/4*z*[1-x-z]^-1
       - 51/16*z
       - 1/2*z^2*[1-x-z]^-1
       + 4*x*z^-1*[1-z]^-1
       + 3/4*x*z^-1
       - 4*x*[1-z]^-1
       - 1/8*x*poly2^-1
       - 117/8*x
       + 11/4*x*z
       - 1/8*x^2*poly2^-1
       + 1/8*x^3*poly2^-1
       + 1/3*pi^2*z^-1*[1-z]^-1
       + 1/3*pi^2*z^-1
       - 1/6*pi^2*[1-z]^-1
       - 3/2*pi^2
       + 1/6*pi^2*z
       - 2/3*pi^2*x*z^-1*[1-z]^-1
       + 1/3*pi^2*x*z^-1
       + 1/3*pi^2*x*[1-z]^-1
       + 1/2*pi^2*x
       + 1/6*pi^2*x*z
       + 207/16*ln(x)
       + 8*ln(x)*z^-1*[1-z]^-1
       - 101/8*ln(x)*z^-1
       - 13/2*ln(x)*[1-z]^-1
       + 1/2*ln(x)*[1-x-z]^-2
       - 1/2*ln(x)*[1-x-z]^-1
       + 3/16*ln(x)*poly2^-2
       + 7/8*ln(x)*poly2^-1
       + 1/4*ln(x)*z*[1-x-z]^-2
       - 7/2*ln(x)*z*[1-x-z]^-1
       + 7/4*ln(x)*z
       - 3/4*ln(x)*z^2*[1-x-z]^-2
       + 7/2*ln(x)*z^2*[1-x-z]^-1
       - 16*ln(x)*x*z^-1*[1-z]^-1
       + 73/4*ln(x)*x*z^-1
       + 10*ln(x)*x*[1-z]^-1
       - 3/2*ln(x)*x*[1-x-z]^-2
       + 3*ln(x)*x*[1-x-z]^-1
       - 3/16*ln(x)*x*poly2^-2
       - ln(x)*x*poly2^-1
       + 3/16*ln(x)*x
       + 3/2*ln(x)*x^2*[1-x-z]^-2
       - 3/2*ln(x)*x^2*[1-x-z]^-1
       - 3/8*ln(x)*x^2*poly2^-2
       + 5/8*ln(x)*x^2*poly2^-1
       - 1/2*ln(x)*x^3*[1-x-z]^-2
       + 3/8*ln(x)*x^3*poly2^-2
       - 1/2*ln(x)*x^3*poly2^-1
       + 3/16*ln(x)*x^4*poly2^-2
       - 3/16*ln(x)*x^5*poly2^-2
       + 1/2*ln(x)*ln(1 - sqrtxz2 + x)*z^-1*sqrtxz2^-1
       + 1/2*ln(x)*ln(1 - sqrtxz2 + x)*[1-z]^-1*sqrtxz2^-1
       + 3/32*ln(x)*ln(1 - sqrtxz2 + x)*sqrtxz2^-1*poly2^-2
       + 13/32*ln(x)*ln(1 - sqrtxz2 + x)*sqrtxz2^-1*poly2^-1
       - ln(x)*ln(1 - sqrtxz2 + x)*sqrtxz2^-1
       + 3/2*ln(x)*ln(1 - sqrtxz2 + x)*x*z^-1*sqrtxz2^-1
       - 3/2*ln(x)*ln(1 - sqrtxz2 + x)*x*[1-z]^-1*sqrtxz2^-1
       - 69/16*ln(x)*ln(1 - sqrtxz2 + x)*x*sqrtxz2^-1
       + 69/8*ln(x)*ln(1 - sqrtxz2 + x)*x*z*sqrtxz2^-1
       + ln(x)*ln(1 - sqrtxz2 + x)*x^2*z^-1*sqrtxz2^-1
       + ln(x)*ln(1 - sqrtxz2 + x)*x^2*[1-z]^-1*sqrtxz2^-1
       - 9/32*ln(x)*ln(1 - sqrtxz2 + x)*x^2*sqrtxz2^-1*poly2^-2
       - 3/16*ln(x)*ln(1 - sqrtxz2 + x)*x^2*sqrtxz2^-1*poly2^-1
       - 11/16*ln(x)*ln(1 - sqrtxz2 + x)*x^2*sqrtxz2^-1
       + 9/32*ln(x)*ln(1 - sqrtxz2 + x)*x^4*sqrtxz2^-1*poly2^-2
       - 7/32*ln(x)*ln(1 - sqrtxz2 + x)*x^4*sqrtxz2^-1*poly2^-1
       - 3/32*ln(x)*ln(1 - sqrtxz2 + x)*x^6*sqrtxz2^-1*poly2^-2
       - 1/2*ln(x)*ln(1 + sqrtxz2 + x)*z^-1*sqrtxz2^-1
       - 1/2*ln(x)*ln(1 + sqrtxz2 + x)*[1-z]^-1*sqrtxz2^-1
       - 3/32*ln(x)*ln(1 + sqrtxz2 + x)*sqrtxz2^-1*poly2^-2
       - 13/32*ln(x)*ln(1 + sqrtxz2 + x)*sqrtxz2^-1*poly2^-1
       + ln(x)*ln(1 + sqrtxz2 + x)*sqrtxz2^-1
       - 3/2*ln(x)*ln(1 + sqrtxz2 + x)*x*z^-1*sqrtxz2^-1
       + 3/2*ln(x)*ln(1 + sqrtxz2 + x)*x*[1-z]^-1*sqrtxz2^-1
       + 69/16*ln(x)*ln(1 + sqrtxz2 + x)*x*sqrtxz2^-1
       - 69/8*ln(x)*ln(1 + sqrtxz2 + x)*x*z*sqrtxz2^-1
       - ln(x)*ln(1 + sqrtxz2 + x)*x^2*z^-1*sqrtxz2^-1
       - ln(x)*ln(1 + sqrtxz2 + x)*x^2*[1-z]^-1*sqrtxz2^-1
       + 9/32*ln(x)*ln(1 + sqrtxz2 + x)*x^2*sqrtxz2^-1*poly2^-2
       + 3/16*ln(x)*ln(1 + sqrtxz2 + x)*x^2*sqrtxz2^-1*poly2^-1
       + 11/16*ln(x)*ln(1 + sqrtxz2 + x)*x^2*sqrtxz2^-1
       - 9/32*ln(x)*ln(1 + sqrtxz2 + x)*x^4*sqrtxz2^-1*poly2^-2
       + 7/32*ln(x)*ln(1 + sqrtxz2 + x)*x^4*sqrtxz2^-1*poly2^-1
       + 3/32*ln(x)*ln(1 + sqrtxz2 + x)*x^6*sqrtxz2^-1*poly2^-2
       - ln(x)*ln(1 + x)
       + ln(x)*ln(1 + x)*z
       - 2*ln(x)*ln(1 + x)*x
       + 2*ln(x)*ln(1 + x)*x*z
       - 2*ln(x)*ln(1 + x)*x^2
       + 33/4*ln(x)^2
       - 5/2*ln(x)^2*z^-1*[1-z]^-1
       - ln(x)^2*z^-1
       + 1/4*ln(x)^2*[1-z]^-1
       + 5*ln(x)^2*x*z^-1*[1-z]^-1
       - 4*ln(x)^2*x*z^-1
       - 1/2*ln(x)^2*x*[1-z]^-1
       - 5/2*ln(x)^2*x
       - 2*ln(x)^2*x*z
       + 3/2*ln(x)^2*x^2
       - 11/2*ln(x)*ln([1-x])
       - 2*ln(x)*ln([1-x])*z^-1*[1-z]^-1
       + 4*ln(x)*ln([1-x])*z^-1
       + 7/2*ln(x)*ln([1-x])*[1-z]^-1
       - ln(x)*ln([1-x])*z
       + 4*ln(x)*ln([1-x])*x*z^-1*[1-z]^-1
       - 8*ln(x)*ln([1-x])*x*z^-1
       - 7*ln(x)*ln([1-x])*x*[1-z]^-1
       + 11*ln(x)*ln([1-x])*x
       + 2*ln(x)*ln([1-x])*x*z
       - ln(x)*ln([1-x])*x^2
       - 11/2*ln(x)*ln(z)
       + 3*ln(x)*ln(z)*z^-1*[1-z]^-1
       + 3/4*ln(x)*ln(z)*z^-1
       + 1/2*ln(x)*ln(z)*[1-z]^-1
       - 1/2*ln(x)*ln(z)*z
       - 6*ln(x)*ln(z)*x*z^-1*[1-z]^-1
       + 5/2*ln(x)*ln(z)*x*z^-1
       + ln(x)*ln(z)*x*[1-z]^-1
       + 9*ln(x)*ln(z)*x
       - 21/2*ln(x)*ln([1-z])
       + 5*ln(x)*ln([1-z])*z^-1
       + 3*ln(x)*ln([1-z])*[1-z]^-1
       - 4*ln(x)*ln([1-z])*x*z^-1
       - 6*ln(x)*ln([1-z])*x*[1-z]^-1
       + 8*ln(x)*ln([1-z])*x
       + ln(x)*ln([1-z])*x*z
       - 14*ln([1-x])
       + 25/4*ln([1-x])*z^-1
       - 1/2*ln([1-x])*[1-z]^-1
       + 3/8*ln([1-x])*z
       - 11/2*ln([1-x])*x*z^-1
       + 2*ln([1-x])*x*[1-z]^-1
       + 21/2*ln([1-x])*x
       - 3/4*ln([1-x])*x*z
       + 7/4*ln([1-x])^2
       - 3/4*ln([1-x])^2*z^-1
       - 1/4*ln([1-x])^2*[1-z]^-1
       + 1/4*ln([1-x])^2*z
       + 3/2*ln([1-x])^2*x*z^-1
       + 1/2*ln([1-x])^2*x*[1-z]^-1
       - 7/2*ln([1-x])^2*x
       - 1/2*ln([1-x])^2*x*z
       + 5/2*ln([1-x])*ln(z)
       - 3/2*ln([1-x])*ln(z)*z^-1
       - 3/2*ln([1-x])*ln(z)*[1-z]^-1
       + 3*ln([1-x])*ln(z)*x*z^-1
       + 3*ln([1-x])*ln(z)*x*[1-z]^-1
       - 5*ln([1-x])*ln(z)*x
       + 9/2*ln([1-x])*ln([1-z])
       - 2*ln([1-x])*ln([1-z])*z^-1
       - ln([1-x])*ln([1-z])*[1-z]^-1
       + 1/2*ln([1-x])*ln([1-z])*z
       + 4*ln([1-x])*ln([1-z])*x*z^-1
       + 2*ln([1-x])*ln([1-z])*x*[1-z]^-1
       - 9*ln([1-x])*ln([1-z])*x
       - ln([1-x])*ln([1-z])*x*z
       - 77/16*ln(z)
       - 4*ln(z)*z^-1*[1-z]^-1
       + 69/8*ln(z)*z^-1
       + 7*ln(z)*[1-z]^-1
       + 3/16*ln(z)*poly2^-2
       + 7/8*ln(z)*poly2^-1
       - 9/8*ln(z)*z
       + 8*ln(z)*x*z^-1*[1-z]^-1
       - 10*ln(z)*x*z^-1
       - 8*ln(z)*x*[1-z]^-1
       + 3/16*ln(z)*x*poly2^-2
       + ln(z)*x*poly2^-1
       - 83/16*ln(z)*x
       + 5/4*ln(z)*x*z
       - 3/8*ln(z)*x^2*poly2^-2
       + 5/8*ln(z)*x^2*poly2^-1
       - 3/8*ln(z)*x^3*poly2^-2
       + 1/2*ln(z)*x^3*poly2^-1
       + 3/16*ln(z)*x^4*poly2^-2
       + 3/16*ln(z)*x^5*poly2^-2
       + 1/2*ln(z)^2
       - 1/2*ln(z)^2*z^-1*[1-z]^-1
       - 1/4*ln(z)^2*z^-1
       + 1/4*ln(z)^2*[1-z]^-1
       - 1/4*ln(z)^2*z
       + ln(z)^2*x*z^-1*[1-z]^-1
       + 1/2*ln(z)^2*x*z^-1
       - 1/2*ln(z)^2*x*[1-z]^-1
       - ln(z)^2*x
       + 1/2*ln(z)^2*x*z
       + 11/2*ln(z)*ln([1-z])
       - 5/2*ln(z)*ln([1-z])*z^-1
       - 5/2*ln(z)*ln([1-z])*[1-z]^-1
       + 5*ln(z)*ln([1-z])*x*z^-1
       + 5*ln(z)*ln([1-z])*x*[1-z]^-1
       - 11*ln(z)*ln([1-z])*x
       - 13*ln([1-z])
       + 23/4*ln([1-z])*z^-1
       - ln([1-z])*[1-z]^-1
       - 1/2*ln([1-z])*[1-x-z]^-2
       + 1/2*ln([1-z])*[1-x-z]^-1
       - 1/4*ln([1-z])*z*[1-x-z]^-2
       + 7/2*ln([1-z])*z*[1-x-z]^-1
       - 25/8*ln([1-z])*z
       + 3/4*ln([1-z])*z^2*[1-x-z]^-2
       - 7/2*ln([1-z])*z^2*[1-x-z]^-1
       - 7/2*ln([1-z])*x*z^-1
       + 4*ln([1-z])*x*[1-z]^-1
       + 3/2*ln([1-z])*x*[1-x-z]^-2
       - 3*ln([1-z])*x*[1-x-z]^-1
       + 13/2*ln([1-z])*x
       - 3/4*ln([1-z])*x*z
       - 3/2*ln([1-z])*x^2*[1-x-z]^-2
       + 3/2*ln([1-z])*x^2*[1-x-z]^-1
       + 1/2*ln([1-z])*x^3*[1-x-z]^-2
       + 11/4*ln([1-z])^2
       - 5/4*ln([1-z])^2*z^-1
       - ln([1-z])^2*[1-z]^-1
       + 1/4*ln([1-z])^2*z
       + 5/2*ln([1-z])^2*x*z^-1
       + 2*ln([1-z])^2*x*[1-z]^-1
       - 11/2*ln([1-z])^2*x
       - 1/2*ln([1-z])^2*x*z
       + 1/2*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*z^-1*sqrtxz2^-1
       + 1/2*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*[1-z]^-1*sqrtxz2^-1
       + 3/32*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*sqrtxz2^-1*poly2^-2
       + 13/32*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*sqrtxz2^-1*poly2^-1
       - Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*sqrtxz2^-1
       + 3/2*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x*z^-1*sqrtxz2^-1
       - 3/2*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x*[1-z]^-1*sqrtxz2^-1
       - 69/16*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x*sqrtxz2^-1
       + 69/8*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x*z*sqrtxz2^-1
       + Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x^2*z^-1*sqrtxz2^-1
       + Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x^2*[1-z]^-1*sqrtxz2^-1
       - 9/32*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x^2*sqrtxz2^-1*poly2^-2
       - 3/16*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x^2*sqrtxz2^-1*poly2^-1
       - 11/16*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x^2*sqrtxz2^-1
       + 9/32*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x^4*sqrtxz2^-1*poly2^-2
       - 7/32*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x^4*sqrtxz2^-1*poly2^-1
       - 3/32*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x^6*sqrtxz2^-1*poly2^-2
       - 1/2*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*z^-1*sqrtxz2^-1
       - 1/2*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*[1-z]^-1*sqrtxz2^-1
       - 3/32*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*sqrtxz2^-1*poly2^-2
       - 13/32*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*sqrtxz2^-1*poly2^-1
       + Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*sqrtxz2^-1
       - 3/2*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x*z^-1*sqrtxz2^-1
       + 3/2*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x*[1-z]^-1*sqrtxz2^-1
       + 69/16*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x*sqrtxz2^-1
       - 69/8*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x*z*sqrtxz2^-1
       - Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x^2*z^-1*sqrtxz2^-1
       - Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x^2*[1-z]^-1*sqrtxz2^-1
       + 9/32*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x^2*sqrtxz2^-1*poly2^-2
       + 3/16*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x^2*sqrtxz2^-1*poly2^-1
       + 11/16*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x^2*sqrtxz2^-1
       - 9/32*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x^4*sqrtxz2^-1*poly2^-2
       + 7/32*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x^4*sqrtxz2^-1*poly2^-1
       + 3/32*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x^6*sqrtxz2^-1*poly2^-2
       - 1/2*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*z^-1*sqrtxz2^-1
       - 1/2*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*[1-z]^-1*sqrtxz2^-1
       - 3/32*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*sqrtxz2^-1*poly2^-2
       - 13/32*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*sqrtxz2^-1*poly2^-1
       + Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*sqrtxz2^-1
       - 3/2*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x*z^-1*sqrtxz2^-1
       + 3/2*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x*[1-z]^-1*sqrtxz2^-1
       + 69/16*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x*sqrtxz2^-1
       - 69/8*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x*z*sqrtxz2^-1
       - Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x^2*z^-1*sqrtxz2^-1
       - Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x^2*[1-z]^-1*sqrtxz2^-1
       + 9/32*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x^2*sqrtxz2^-1*poly2^-2
       + 3/16*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x^2*sqrtxz2^-1*poly2^-1
       + 11/16*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x^2*sqrtxz2^-1
       - 9/32*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x^4*sqrtxz2^-1*poly2^-2
       + 7/32*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x^4*sqrtxz2^-1*poly2^-1
       + 3/32*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x^6*sqrtxz2^-1*poly2^-2
       + 1/2*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*z^-1*sqrtxz2^-1
       + 1/2*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*[1-z]^-1*sqrtxz2^-1
       + 3/32*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*sqrtxz2^-1*poly2^-2
       + 13/32*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*sqrtxz2^-1*poly2^-1
       - Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*sqrtxz2^-1
       + 3/2*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x*z^-1*sqrtxz2^-1
       - 3/2*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x*[1-z]^-1*sqrtxz2^-1
       - 69/16*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x*sqrtxz2^-1
       + 69/8*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x*z*sqrtxz2^-1
       + Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x^2*z^-1*sqrtxz2^-1
       + Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x^2*[1-z]^-1*sqrtxz2^-1
       - 9/32*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x^2*sqrtxz2^-1*poly2^-2
       - 3/16*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x^2*sqrtxz2^-1*poly2^-1
       - 11/16*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x^2*sqrtxz2^-1
       + 9/32*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x^4*sqrtxz2^-1*poly2^-2
       - 7/32*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x^4*sqrtxz2^-1*poly2^-1
       - 3/32*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x^6*sqrtxz2^-1*poly2^-2
       - 3/2*Li2(1 - x*z^-1)
       + Li2(1 - x*z^-1)*z^-1*[1-z]^-1
       - 1/2*Li2(1 - x*z^-1)*z^-1
       - 1/2*Li2(1 - x*z^-1)*[1-z]^-1
       - 2*Li2(1 - x*z^-1)*x*z^-1*[1-z]^-1
       + Li2(1 - x*z^-1)*x*z^-1
       + Li2(1 - x*z^-1)*x*[1-z]^-1
       + 3*Li2(1 - x*z^-1)*x
       - Li2( - x)
       + Li2( - x)*z
       - 2*Li2( - x)*x
       + 2*Li2( - x)*x*z
       - 2*Li2( - x)*x^2
       + 4*Li2(x)
       - 2*Li2(x)*z^-1*[1-z]^-1
       - 1/2*Li2(x)*z^-1
       + 2*Li2(x)*[1-z]^-1
       - Li2(x)*z
       + 4*Li2(x)*x*z^-1*[1-z]^-1
       - 5*Li2(x)*x*z^-1
       - 4*Li2(x)*x*[1-z]^-1
       + 5*Li2(x)*x
       + Li2(x)*x*z
       - Li2(x)*x^2
       + 2*Li2(z)
       - 1/2*Li2(z)*z^-1
       - 1/2*Li2(z)*[1-z]^-1
       + Li2(z)*x*z^-1
       + Li2(z)*x*[1-z]^-1
       - 4*Li2(z)*x
       )

       + LMUR*NF * (  - 1/3
       + 1/6*z^-1
       - 1/3*x*z^-1
       + 2/3*x
       )

       + LMUR*NC * ( 11/6
       - 11/12*z^-1
       + 11/6*x*z^-1
       - 11/3*x
       )

       + LMUF*NC^-1 * (  - 5/4
       + 3/4*z
       + 3/2*x
       - x*z
       - 1/4*ln(x)
       + 1/4*ln(x)*z
       - 1/2*ln(x)*x
       + 1/2*ln(x)*x*z
       + 1/4*ln([1-x])
       + 1/4*ln([1-x])*z
       - 1/2*ln([1-x])*x
       - 1/2*ln([1-x])*x*z
       + 1/4*ln(z)
       - 1/2*ln(z)*[1-z]^-1
       + 1/4*ln(z)*z
       + ln(z)*x*[1-z]^-1
       - 1/2*ln(z)*x
       - 1/2*ln(z)*x*z
       + 1/4*ln([1-z])
       + 1/4*ln([1-z])*z
       - 1/2*ln([1-z])*x
       - 1/2*ln([1-z])*x*z
       )

       + LMUF*NF * ( 1/3
       - 1/6*z^-1
       + 1/3*x*z^-1
       - 2/3*x
       )

       + LMUF*NC * ( 137/12
       - 61/12*z^-1
       - 3/4*z
       + 25/6*x*z^-1
       - 59/6*x
       + x*z
       + 17/4*ln(x)
       - 2*ln(x)*z^-1
       - 1/4*ln(x)*z
       - 2*ln(x)*x*z^-1
       + 9/2*ln(x)*x
       - 1/2*ln(x)*x*z
       - 9/4*ln([1-x])
       + ln([1-x])*z^-1
       - 1/4*ln([1-x])*z
       - 2*ln([1-x])*x*z^-1
       + 9/2*ln([1-x])*x
       + 1/2*ln([1-x])*x*z
       - 1/4*ln(z)
       + 1/2*ln(z)*[1-z]^-1
       - 1/4*ln(z)*z
       - ln(z)*x*[1-z]^-1
       + 1/2*ln(z)*x
       + 1/2*ln(z)*x*z
       - 1/4*ln([1-z])
       - 1/4*ln([1-z])*z
       + 1/2*ln([1-z])*x
       + 1/2*ln([1-z])*x*z
       )

       + LMUA*NC^-1 * (  - 1/4
       - 3/8*z
       + 1/4*x*z
       - 1/4*ln(x)
       - 1/4*ln(x)*z
       + 1/2*ln(x)*x
       + 1/2*ln(x)*x*z
       + 1/4*ln([1-x])
       + 1/4*ln([1-x])*z
       - 1/2*ln([1-x])*x
       - 1/2*ln([1-x])*x*z
       - 3/4*ln(z)
       + 1/2*ln(z)*[1-z]^-1
       - 1/4*ln(z)*z
       - ln(z)*x*[1-z]^-1
       + 3/2*ln(z)*x
       + 1/2*ln(z)*x*z
       + 5/4*ln([1-z])
       - 1/2*ln([1-z])*z^-1
       + 1/4*ln([1-z])*z
       + ln([1-z])*x*z^-1
       - 5/2*ln([1-z])*x
       - 1/2*ln([1-z])*x*z
       )

       + LMUA*NC * ( 1/4
       + 3/8*z
       - 1/4*x*z
       + 1/4*ln(x)
       + 1/4*ln(x)*z
       - 1/2*ln(x)*x
       - 1/2*ln(x)*x*z
       - 1/4*ln([1-x])
       - 1/4*ln([1-x])*z
       + 1/2*ln([1-x])*x
       + 1/2*ln([1-x])*x*z
       + 3/4*ln(z)
       - 1/2*ln(z)*[1-z]^-1
       + 1/4*ln(z)*z
       + ln(z)*x*[1-z]^-1
       - 3/2*ln(z)*x
       - 1/2*ln(z)*x*z
       - 5/4*ln([1-z])
       + 1/2*ln([1-z])*z^-1
       - 1/4*ln([1-z])*z
       - ln([1-z])*x*z^-1
       + 5/2*ln([1-z])*x
       + 1/2*ln([1-z])*x*z
       )

       + LMUA*LMUF*NC^-1 * (  - 1/4
       - 1/4*z
       + 1/2*x
       + 1/2*x*z
       )

       + LMUA*LMUF*NC * ( 1/4
       + 1/4*z
       - 1/2*x
       - 1/2*x*z
       )

       + Dd([1-z])*NC^-1 * (  - 35/16
       + 45/16*x
       + 4*zeta3
       - 8*zeta3*x
       + 5/24*pi^2
       - 1/4*pi^2*x
       - 19/8*ln(x)
       - 1/16*ln(x)*x
       + 1/4*ln(x)*pi^2
       - 1/2*ln(x)*pi^2*x
       - 37/32*ln(x)^2
       + 3/4*ln(x)^2*x
       - 1/48*ln(x)^3
       + 1/24*ln(x)^3*x
       + 3/8*ln(x)^2*ln([1-x])
       - 3/4*ln(x)^2*ln([1-x])*x
       + 21/8*ln(x)*ln([1-x])
       - 3*ln(x)*ln([1-x])*x
       - 5/4*ln(x)*ln([1-x])^2
       + 5/2*ln(x)*ln([1-x])^2*x
       + 1/2*ln(x)*Li2(x)
       - ln(x)*Li2(x)*x
       + 7/4*ln([1-x])
       - 1/2*ln([1-x])*x
       - 1/24*ln([1-x])*pi^2
       + 1/12*ln([1-x])*pi^2*x
       - 21/16*ln([1-x])^2
       + 3/2*ln([1-x])^2*x
       + 5/24*ln([1-x])^3
       - 5/12*ln([1-x])^3*x
       - 1/4*ln([1-x])*Li2(1 - x)
       + 1/2*ln([1-x])*Li2(1 - x)*x
       - 3/4*ln([1-x])*Li2(x)
       + 3/2*ln([1-x])*Li2(x)*x
       - 3/2*Li3(1 - x)
       + 3*Li3(1 - x)*x
       - 5/4*Li3(x)
       + 5/2*Li3(x)*x
       + 3/8*Li2(x)
       - 1/2*Li2(x)*x
       )

       + Dd([1-z])*NC * ( 263/16
       - 1/3*rln2^3
       - 281/16*x
       - 2/3*x*rln2^3
       - 13/4*zeta3
       + 9/2*zeta3*x
       - 23/24*pi^2
       + 1/6*pi^2*rln2
       + 19/12*pi^2*x
       + 1/3*pi^2*x*rln2
       + 1/2*ln(1 + x)*rln2^2
       + ln(1 + x)*x*rln2^2
       + 1/12*ln(1 + x)*pi^2
       + 1/6*ln(1 + x)*pi^2*x
       - 1/6*ln(1 + x)^3
       - 1/3*ln(1 + x)^3*x
       + 99/8*ln(x)
       + 1/16*ln(x)*x
       - 3/4*ln(x)*pi^2
       + 1/6*ln(x)*pi^2*x
       + ln(x)*ln(1 + x)
       + ln(x)*ln(1 + x)*x
       + 153/32*ln(x)^2
       - 25/4*ln(x)^2*x
       + ln(x)^2*ln(1 + x)
       + 2*ln(x)^2*ln(1 + x)*x
       + 29/48*ln(x)^3
       + 19/24*ln(x)^3*x
       - 11/8*ln(x)^2*ln([1-x])
       + 11/4*ln(x)^2*ln([1-x])*x
       - 69/8*ln(x)*ln([1-x])
       + 9*ln(x)*ln([1-x])*x
       - ln(x)*ln([1-x])*ln(1 + x)
       - 2*ln(x)*ln([1-x])*ln(1 + x)*x
       - 7/4*ln(x)*ln([1-x])^2
       - 11/2*ln(x)*ln([1-x])^2*x
       + ln(x)*Li2( - x)
       + 2*ln(x)*Li2( - x)*x
       + 1/2*ln(x)*Li2(x)
       + 5*ln(x)*Li2(x)*x
       - 19/4*ln([1-x])
       + 1/2*ln([1-x])*rln2^2
       + 15/4*ln([1-x])*x
       + ln([1-x])*x*rln2^2
       + 19/24*ln([1-x])*pi^2
       + 1/12*ln([1-x])*pi^2*x
       - ln([1-x])*ln(1 + x)*rln2
       - 2*ln([1-x])*ln(1 + x)*x*rln2
       + 1/2*ln([1-x])*ln(1 + x)^2
       + ln([1-x])*ln(1 + x)^2*x
       + 61/16*ln([1-x])^2
       - 4*ln([1-x])^2*x
       - 7/24*ln([1-x])^3
       + 7/12*ln([1-x])^3*x
       - 7/4*ln([1-x])*Li2(1 - x)
       - 5/2*ln([1-x])*Li2(1 - x)*x
       - ln([1-x])*Li2( - x)
       - 2*ln([1-x])*Li2( - x)*x
       - 15/4*ln([1-x])*Li2(x)
       - 9/2*ln([1-x])*Li2(x)*x
       - 1/4*ln([1+x])*pi^2
       - 1/2*ln([1+x])*pi^2*x
       + Li3(1/2 - 1/2*x)
       + 2*Li3(1/2 - 1/2*x)*x
       + Li3(1/2 + 1/2*x)
       + 2*Li3(1/2 + 1/2*x)*x
       - 7/2*Li3(1 - x)
       - 3*Li3(1 - x)*x
       + Li3(1/(1 + x) - 1/(1 + x)*x)
       + 2*Li3(1/(1 + x) - 1/(1 + x)*x)*x
       + 5/4*Li3(x)
       - 9/2*Li3(x)*x
       + Li2( - x)
       + Li2( - x)*x
       - 7/8*Li2(x)
       - 3/2*Li2(x)*x
       )

       + Dd([1-z])*LMUR*NF * (  - 1/3
       + 1/3*x
       - 1/6*ln(x)
       + 1/3*ln(x)*x
       + 1/6*ln([1-x])
       - 1/3*ln([1-x])*x
       )

       + Dd([1-z])*LMUR*NC * ( 11/6
       - 11/6*x
       + 11/12*ln(x)
       - 11/6*ln(x)*x
       - 11/12*ln([1-x])
       + 11/6*ln([1-x])*x
       )

       + Dd([1-z])*LMUF*NC^-1 * (  - 1/2
       - 7/8*x
       + 1/24*pi^2
       - 1/12*pi^2*x
       - 11/8*ln(x)
       + 5/4*ln(x)*x
       - 1/4*ln(x)^2
       + 1/2*ln(x)^2*x
       + ln(x)*ln([1-x])
       - 2*ln(x)*ln([1-x])*x
       + 7/4*ln([1-x])
       - 7/4*ln([1-x])*x
       - 1/2*ln([1-x])^2
       + ln([1-x])^2*x
       + 1/4*Li2(x)
       - 1/2*Li2(x)*x
       )

       + Dd([1-z])*LMUF*NF * ( 1/3
       - 1/3*x
       + 1/6*ln(x)
       - 1/3*ln(x)*x
       - 1/6*ln([1-x])
       + 1/3*ln([1-x])*x
       )

       + Dd([1-z])*LMUF*NC * ( 8/3
       - 43/24*x
       - 3/8*pi^2
       + 1/12*pi^2*x
       + 143/24*ln(x)
       - 89/12*ln(x)*x
       + ln(x)*ln(1 + x)
       + 2*ln(x)*ln(1 + x)*x
       + 7/4*ln(x)^2
       + 3/2*ln(x)^2*x
       - 2*ln(x)*ln([1-x])
       + 4*ln(x)*ln([1-x])*x
       - 41/6*ln([1-x])
       + 71/12*ln([1-x])*x
       + ln([1-x])^2
       - 2*ln([1-x])^2*x
       + Li2( - x)
       + 2*Li2( - x)*x
       + 7/4*Li2(x)
       + 5/2*Li2(x)*x
       )

       + Dd([1-z])*LMUF*LMUR*NF * (  - 1/6
       + 1/3*x
       )

       + Dd([1-z])*LMUF*LMUR*NC * ( 11/12
       - 11/6*x
       )

       + Dd([1-z])*LMUF^2*NC^-1 * (  - 3/16
       - 1/8*ln(x)
       + 1/4*ln(x)*x
       + 1/4*ln([1-x])
       - 1/2*ln([1-x])*x
       )

       + Dd([1-z])*LMUF^2*NF * ( 1/6
       - 1/3*x
       )

       + Dd([1-z])*LMUF^2*NC * ( 109/48
       - 7/6*x
       + 9/8*ln(x)
       + 3/4*ln(x)*x
       - 3/4*ln([1-x])
       + 3/2*ln([1-x])*x
       )

       + Dd([1-z])*LMUA*NC^-1 * ( 3/4
       - 3/4*x
       + 1/12*pi^2
       - 1/6*pi^2*x
       + 3/8*ln(x)
       - 3/4*ln(x)*x
       - 3/8*ln([1-x])
       + 3/4*ln([1-x])*x
       )

       + Dd([1-z])*LMUA*NC * (  - 3/4
       + 3/4*x
       - 1/12*pi^2
       + 1/6*pi^2*x
       - 3/8*ln(x)
       + 3/4*ln(x)*x
       + 3/8*ln([1-x])
       - 3/4*ln([1-x])*x
       )

       + Dd([1-z])*LMUA*LMUF*NC^-1 * ( 3/8
       - 3/4*x
       )

       + Dd([1-z])*LMUA*LMUF*NC * (  - 3/8
       + 3/4*x
       )

       + Dn(0,[1-z])*NC^-1 * ( 5/4
       + 1/8*x
       - 1/8*pi^2
       + 1/4*pi^2*x
       + 7/4*ln(x)
       - 2*ln(x)*x
       + 1/4*ln(x)^2
       - 1/2*ln(x)^2*x
       - ln(x)*ln([1-x])
       + 2*ln(x)*ln([1-x])*x
       - 17/8*ln([1-x])
       + 5/2*ln([1-x])*x
       + 1/2*ln([1-x])^2
       - ln([1-x])^2*x
       - 1/4*Li2(x)
       + 1/2*Li2(x)*x
       )

       + Dn(0,[1-z])*NC * (  - 21/4
       + 35/8*x
       + 11/24*pi^2
       - 1/4*pi^2*x
       - 29/4*ln(x)
       + 10*ln(x)*x
       - ln(x)*ln(1 + x)
       - 2*ln(x)*ln(1 + x)*x
       - 7/4*ln(x)^2
       - 3/2*ln(x)^2*x
       + 2*ln(x)*ln([1-x])
       - 4*ln(x)*ln([1-x])*x
       + 65/8*ln([1-x])
       - 17/2*ln([1-x])*x
       - ln([1-x])^2
       + 2*ln([1-x])^2*x
       - Li2( - x)
       - 2*Li2( - x)*x
       - 7/4*Li2(x)
       - 5/2*Li2(x)*x
       )

       + Dn(0,[1-z])*LMUR*NF * ( 1/6
       - 1/3*x
       )

       + Dn(0,[1-z])*LMUR*NC * (  - 11/12
       + 11/6*x
       )

       + Dn(0,[1-z])*LMUF*NC^-1 * ( 3/4
       - 3/4*x
       + 1/4*ln(x)
       - 1/2*ln(x)*x
       - 1/2*ln([1-x])
       + ln([1-x])*x
       )

       + Dn(0,[1-z])*LMUF*NF * (  - 1/6
       + 1/3*x
       )

       + Dn(0,[1-z])*LMUF*NC * (  - 35/6
       + 59/12*x
       - 9/4*ln(x)
       - 3/2*ln(x)*x
       + 3/2*ln([1-x])
       - 3*ln([1-x])*x
       )

       + Dn(0,[1-z])*LMUA*NC^-1 * ( 5/8
       - 1/4*x
       + 1/2*ln(x)
       - ln(x)*x
       - 1/2*ln([1-x])
       + ln([1-x])*x
       )

       + Dn(0,[1-z])*LMUA*NC * (  - 5/8
       + 1/4*x
       - 1/2*ln(x)
       + ln(x)*x
       + 1/2*ln([1-x])
       - ln([1-x])*x
       )

       + Dn(0,[1-z])*LMUA*LMUF*NC^-1 * ( 1/2
       - x
       )

       + Dn(0,[1-z])*LMUA*LMUF*NC * (  - 1/2
       + x
       )

       + Dn(1,[1-z])*NC^-1 * (  - 17/8
       + 5/2*x
       - 3/4*ln(x)
       + 3/2*ln(x)*x
       + ln([1-x])
       - 2*ln([1-x])*x
       )

       + Dn(1,[1-z])*NC * ( 65/8
       - 17/2*x
       + 11/4*ln(x)
       + 1/2*ln(x)*x
       - 2*ln([1-x])
       + 4*ln([1-x])*x
       )

       + Dn(1,[1-z])*LMUF*NC^-1 * (  - 1/2
       + x
       )

       + Dn(1,[1-z])*LMUF*NC * ( 1/2
       - x
       )

       + Dn(1,[1-z])*LMUA*NC^-1 * (  - 1
       + 2*x
       )

       + Dn(1,[1-z])*LMUA*NC * ( 1
       - 2*x
       )

       + Dn(2,[1-z])*NC^-1 * ( 3/4
       - 3/2*x
       )

       + Dn(2,[1-z])*NC * (  - 3/4
       + 3/2*x
       )

       + T(u1)*NC^-1 * (  - 1
       + 2*z^-1
       + 2*[1-z]^-1
       + x*z^-1
       + x*[1-z]^-1
       - 4*x
       + 11/12*pi^2*z^-1
       + 5/6*pi^2*[1-z]^-1
       - 7/12*pi^2
       - 11/6*pi^2*x*z^-1
       - 5/3*pi^2*x*[1-z]^-1
       + 7/6*pi^2*x
       + 2/3*pi^2*x^2*z^-1
       + 1/2*pi^2*x^2*[1-z]^-1
       - 6*ln(x)
       + 9/2*ln(x)*z^-1
       + 9/2*ln(x)*[1-z]^-1
       - 12*ln(x)*x*z^-1
       - 12*ln(x)*x*[1-z]^-1
       + 18*ln(x)*x
       + 13/4*ln(x)^2
       - 5*ln(x)^2*z^-1
       - 19/4*ln(x)^2*[1-z]^-1
       + 10*ln(x)^2*x*z^-1
       + 19/2*ln(x)^2*x*[1-z]^-1
       - 13/2*ln(x)^2*x
       - 7/2*ln(x)^2*x^2*z^-1
       - 3*ln(x)^2*x^2*[1-z]^-1
       - 11/2*ln(x)*ln([1-x])
       + 17/2*ln(x)*ln([1-x])*z^-1
       + 8*ln(x)*ln([1-x])*[1-z]^-1
       - 17*ln(x)*ln([1-x])*x*z^-1
       - 16*ln(x)*ln([1-x])*x*[1-z]^-1
       + 11*ln(x)*ln([1-x])*x
       + 6*ln(x)*ln([1-x])*x^2*z^-1
       + 5*ln(x)*ln([1-x])*x^2*[1-z]^-1
       - 4*ln(x)*ln(z)
       + 11/2*ln(x)*ln(z)*z^-1
       + 13/2*ln(x)*ln(z)*[1-z]^-1
       - 11*ln(x)*ln(z)*x*z^-1
       - 13*ln(x)*ln(z)*x*[1-z]^-1
       + 8*ln(x)*ln(z)*x
       + 3*ln(x)*ln(z)*x^2*z^-1
       + 5*ln(x)*ln(z)*x^2*[1-z]^-1
       - 5*ln(x)*ln([1-z])
       + 8*ln(x)*ln([1-z])*z^-1
       + 7*ln(x)*ln([1-z])*[1-z]^-1
       - 16*ln(x)*ln([1-z])*x*z^-1
       - 14*ln(x)*ln([1-z])*x*[1-z]^-1
       + 10*ln(x)*ln([1-z])*x
       + 6*ln(x)*ln([1-z])*x^2*z^-1
       + 4*ln(x)*ln([1-z])*x^2*[1-z]^-1
       + 1/2*ln(x)*ln([1-x-z])
       - ln(x)*ln([1-x-z])*z^-1
       - 1/2*ln(x)*ln([1-x-z])*[1-z]^-1
       + 2*ln(x)*ln([1-x-z])*x*z^-1
       + ln(x)*ln([1-x-z])*x*[1-z]^-1
       - ln(x)*ln([1-x-z])*x
       - ln(x)*ln([1-x-z])*x^2*z^-1
       + 1/2*ln(x)*ln([x-z])
       - 1/2*ln(x)*ln([x-z])*z^-1
       - ln(x)*ln([x-z])*[1-z]^-1
       + ln(x)*ln([x-z])*x*z^-1
       + 2*ln(x)*ln([x-z])*x*[1-z]^-1
       - ln(x)*ln([x-z])*x
       - ln(x)*ln([x-z])*x^2*[1-z]^-1
       + 4*ln([1-x])
       - 3*ln([1-x])*z^-1
       - 3*ln([1-x])*[1-z]^-1
       + 8*ln([1-x])*x*z^-1
       + 8*ln([1-x])*x*[1-z]^-1
       - 12*ln([1-x])*x
       + 2*ln([1-x])^2
       - 3*ln([1-x])^2*z^-1
       - 3*ln([1-x])^2*[1-z]^-1
       + 6*ln([1-x])^2*x*z^-1
       + 6*ln([1-x])^2*x*[1-z]^-1
       - 4*ln([1-x])^2*x
       - 2*ln([1-x])^2*x^2*z^-1
       - 2*ln([1-x])^2*x^2*[1-z]^-1
       + 3*ln([1-x])*ln(z)
       - 4*ln([1-x])*ln(z)*z^-1
       - 5*ln([1-x])*ln(z)*[1-z]^-1
       + 8*ln([1-x])*ln(z)*x*z^-1
       + 10*ln([1-x])*ln(z)*x*[1-z]^-1
       - 6*ln([1-x])*ln(z)*x
       - 2*ln([1-x])*ln(z)*x^2*z^-1
       - 4*ln([1-x])*ln(z)*x^2*[1-z]^-1
       + 7/2*ln([1-x])*ln([1-z])
       - 6*ln([1-x])*ln([1-z])*z^-1
       - 9/2*ln([1-x])*ln([1-z])*[1-z]^-1
       + 12*ln([1-x])*ln([1-z])*x*z^-1
       + 9*ln([1-x])*ln([1-z])*x*[1-z]^-1
       - 7*ln([1-x])*ln([1-z])*x
       - 5*ln([1-x])*ln([1-z])*x^2*z^-1
       - 2*ln([1-x])*ln([1-z])*x^2*[1-z]^-1
       + 3*ln(z)
       - 2*ln(z)*z^-1
       - 5/2*ln(z)*[1-z]^-1
       + 6*ln(z)*x*z^-1
       + 6*ln(z)*x*[1-z]^-1
       - 9*ln(z)*x
       + 5/4*ln(z)^2
       - 3/2*ln(z)^2*z^-1
       - 9/4*ln(z)^2*[1-z]^-1
       + 3*ln(z)^2*x*z^-1
       + 9/2*ln(z)^2*x*[1-z]^-1
       - 5/2*ln(z)^2*x
       - 1/2*ln(z)^2*x^2*z^-1
       - 2*ln(z)^2*x^2*[1-z]^-1
       + 2*ln(z)*ln([1-z])
       - 3*ln(z)*ln([1-z])*z^-1
       - 3*ln(z)*ln([1-z])*[1-z]^-1
       + 6*ln(z)*ln([1-z])*x*z^-1
       + 6*ln(z)*ln([1-z])*x*[1-z]^-1
       - 4*ln(z)*ln([1-z])*x
       - 2*ln(z)*ln([1-z])*x^2*z^-1
       - 2*ln(z)*ln([1-z])*x^2*[1-z]^-1
       - 1/2*ln(z)*ln([x-z])
       + 1/2*ln(z)*ln([x-z])*z^-1
       + ln(z)*ln([x-z])*[1-z]^-1
       - ln(z)*ln([x-z])*x*z^-1
       - 2*ln(z)*ln([x-z])*x*[1-z]^-1
       + ln(z)*ln([x-z])*x
       + ln(z)*ln([x-z])*x^2*[1-z]^-1
       + 3*ln([1-z])
       - 5/2*ln([1-z])*z^-1
       - 2*ln([1-z])*[1-z]^-1
       + 6*ln([1-z])*x*z^-1
       + 6*ln([1-z])*x*[1-z]^-1
       - 9*ln([1-z])*x
       + 3/2*ln([1-z])^2
       - 11/4*ln([1-z])^2*z^-1
       - 7/4*ln([1-z])^2*[1-z]^-1
       + 11/2*ln([1-z])^2*x*z^-1
       + 7/2*ln([1-z])^2*x*[1-z]^-1
       - 3*ln([1-z])^2*x
       - 5/2*ln([1-z])^2*x^2*z^-1
       - 1/2*ln([1-z])^2*x^2*[1-z]^-1
       - 1/2*ln([1-z])*ln([1-x-z])
       + ln([1-z])*ln([1-x-z])*z^-1
       + 1/2*ln([1-z])*ln([1-x-z])*[1-z]^-1
       - 2*ln([1-z])*ln([1-x-z])*x*z^-1
       - ln([1-z])*ln([1-x-z])*x*[1-z]^-1
       + ln([1-z])*ln([1-x-z])*x
       + ln([1-z])*ln([1-x-z])*x^2*z^-1
       - 1/2*Li2(x^-1*[1-x]*z*[1-z]^-1)
       + 1/2*Li2(x^-1*[1-x]*z*[1-z]^-1)*z^-1
       + Li2(x^-1*[1-x]*z*[1-z]^-1)*[1-z]^-1
       - Li2(x^-1*[1-x]*z*[1-z]^-1)*x*z^-1
       - 2*Li2(x^-1*[1-x]*z*[1-z]^-1)*x*[1-z]^-1
       + Li2(x^-1*[1-x]*z*[1-z]^-1)*x
       + Li2(x^-1*[1-x]*z*[1-z]^-1)*x^2*[1-z]^-1
       - 1/2*Li2([1-x]^-1*z)
       + Li2([1-x]^-1*z)*z^-1
       + 1/2*Li2([1-x]^-1*z)*[1-z]^-1
       - 2*Li2([1-x]^-1*z)*x*z^-1
       - Li2([1-x]^-1*z)*x*[1-z]^-1
       + Li2([1-x]^-1*z)*x
       + Li2([1-x]^-1*z)*x^2*z^-1
       + 1/2*Li2([1-x]*[1-z]^-1)
       - 1/2*Li2([1-x]*[1-z]^-1)*z^-1
       - Li2([1-x]*[1-z]^-1)*[1-z]^-1
       + Li2([1-x]*[1-z]^-1)*x*z^-1
       + 2*Li2([1-x]*[1-z]^-1)*x*[1-z]^-1
       - Li2([1-x]*[1-z]^-1)*x
       - Li2([1-x]*[1-z]^-1)*x^2*[1-z]^-1
       + 1/2*Li2(x*[1-x]^-1*z*[1-z]^-1)
       - Li2(x*[1-x]^-1*z*[1-z]^-1)*z^-1
       - 1/2*Li2(x*[1-x]^-1*z*[1-z]^-1)*[1-z]^-1
       + 2*Li2(x*[1-x]^-1*z*[1-z]^-1)*x*z^-1
       + Li2(x*[1-x]^-1*z*[1-z]^-1)*x*[1-z]^-1
       - Li2(x*[1-x]^-1*z*[1-z]^-1)*x
       - Li2(x*[1-x]^-1*z*[1-z]^-1)*x^2*z^-1
       - 1/2*Li2(z)*z^-1
       + 1/2*Li2(z)*[1-z]^-1
       + Li2(z)*x*z^-1
       - Li2(z)*x*[1-z]^-1
       - Li2(z)*x^2*z^-1
       + Li2(z)*x^2*[1-z]^-1
       )

       + T(u1)*NC * ( 2
       - z^-1
       - [1-z]^-1
       - 1/12*pi^2*z^-1
       - 1/12*pi^2*[1-z]^-1
       + 1/6*pi^2
       + 1/6*pi^2*x*z^-1
       + 1/6*pi^2*x*[1-z]^-1
       - 1/3*pi^2*x
       + 3*ln(x)
       - 3/2*ln(x)*z^-1
       - 3/2*ln(x)*[1-z]^-1
       + 6*ln(x)*x*z^-1
       + 6*ln(x)*x*[1-z]^-1
       - 18*ln(x)*x
       - 3*ln(x)^2
       + 3/2*ln(x)^2*z^-1
       + 3/2*ln(x)^2*[1-z]^-1
       - 3*ln(x)^2*x*z^-1
       - 3*ln(x)^2*x*[1-z]^-1
       + 6*ln(x)^2*x
       + 4*ln(x)*ln([1-x])
       - 2*ln(x)*ln([1-x])*z^-1
       - 2*ln(x)*ln([1-x])*[1-z]^-1
       + 4*ln(x)*ln([1-x])*x*z^-1
       + 4*ln(x)*ln([1-x])*x*[1-z]^-1
       - 8*ln(x)*ln([1-x])*x
       + 5*ln(x)*ln(z)
       - 5/2*ln(x)*ln(z)*z^-1
       - 5/2*ln(x)*ln(z)*[1-z]^-1
       + 5*ln(x)*ln(z)*x*z^-1
       + 5*ln(x)*ln(z)*x*[1-z]^-1
       - 10*ln(x)*ln(z)*x
       + 5*ln(x)*ln([1-z])
       - 5/2*ln(x)*ln([1-z])*z^-1
       - 5/2*ln(x)*ln([1-z])*[1-z]^-1
       + 5*ln(x)*ln([1-z])*x*z^-1
       + 5*ln(x)*ln([1-z])*x*[1-z]^-1
       - 10*ln(x)*ln([1-z])*x
       - ln(x)*ln([1-x-z])
       + 1/2*ln(x)*ln([1-x-z])*z^-1
       + 1/2*ln(x)*ln([1-x-z])*[1-z]^-1
       - ln(x)*ln([1-x-z])*x*z^-1
       - ln(x)*ln([1-x-z])*x*[1-z]^-1
       + 2*ln(x)*ln([1-x-z])*x
       - ln(x)*ln([x-z])
       + 1/2*ln(x)*ln([x-z])*z^-1
       + 1/2*ln(x)*ln([x-z])*[1-z]^-1
       - ln(x)*ln([x-z])*x*z^-1
       - ln(x)*ln([x-z])*x*[1-z]^-1
       + 2*ln(x)*ln([x-z])*x
       - ln([1-x])
       + 1/2*ln([1-x])*z^-1
       + 1/2*ln([1-x])*[1-z]^-1
       - 2*ln([1-x])*x*z^-1
       - 2*ln([1-x])*x*[1-z]^-1
       + 6*ln([1-x])*x
       - 1/2*ln([1-x])^2
       + 1/4*ln([1-x])^2*z^-1
       + 1/4*ln([1-x])^2*[1-z]^-1
       - 1/2*ln([1-x])^2*x*z^-1
       - 1/2*ln([1-x])^2*x*[1-z]^-1
       + ln([1-x])^2*x
       - 2*ln([1-x])*ln(z)
       + ln([1-x])*ln(z)*z^-1
       + ln([1-x])*ln(z)*[1-z]^-1
       - 2*ln([1-x])*ln(z)*x*z^-1
       - 2*ln([1-x])*ln(z)*x*[1-z]^-1
       + 4*ln([1-x])*ln(z)*x
       - 3*ln([1-x])*ln([1-z])
       + 3/2*ln([1-x])*ln([1-z])*z^-1
       + 3/2*ln([1-x])*ln([1-z])*[1-z]^-1
       - 3*ln([1-x])*ln([1-z])*x*z^-1
       - 3*ln([1-x])*ln([1-z])*x*[1-z]^-1
       + 6*ln([1-x])*ln([1-z])*x
       - 2*ln(z)
       + ln(z)*z^-1
       + ln(z)*[1-z]^-1
       - 4*ln(z)*x*z^-1
       - 4*ln(z)*x*[1-z]^-1
       + 12*ln(z)*x
       - 2*ln(z)^2
       + ln(z)^2*z^-1
       + ln(z)^2*[1-z]^-1
       - 2*ln(z)^2*x*z^-1
       - 2*ln(z)^2*x*[1-z]^-1
       + 4*ln(z)^2*x
       - 4*ln(z)*ln([1-z])
       + 2*ln(z)*ln([1-z])*z^-1
       + 2*ln(z)*ln([1-z])*[1-z]^-1
       - 4*ln(z)*ln([1-z])*x*z^-1
       - 4*ln(z)*ln([1-z])*x*[1-z]^-1
       + 8*ln(z)*ln([1-z])*x
       + ln(z)*ln([x-z])
       - 1/2*ln(z)*ln([x-z])*z^-1
       - 1/2*ln(z)*ln([x-z])*[1-z]^-1
       + ln(z)*ln([x-z])*x*z^-1
       + ln(z)*ln([x-z])*x*[1-z]^-1
       - 2*ln(z)*ln([x-z])*x
       - 2*ln([1-z])
       + ln([1-z])*z^-1
       + ln([1-z])*[1-z]^-1
       - 4*ln([1-z])*x*z^-1
       - 4*ln([1-z])*x*[1-z]^-1
       + 12*ln([1-z])*x
       - 3/2*ln([1-z])^2
       + 3/4*ln([1-z])^2*z^-1
       + 3/4*ln([1-z])^2*[1-z]^-1
       - 3/2*ln([1-z])^2*x*z^-1
       - 3/2*ln([1-z])^2*x*[1-z]^-1
       + 3*ln([1-z])^2*x
       + ln([1-z])*ln([1-x-z])
       - 1/2*ln([1-z])*ln([1-x-z])*z^-1
       - 1/2*ln([1-z])*ln([1-x-z])*[1-z]^-1
       + ln([1-z])*ln([1-x-z])*x*z^-1
       + ln([1-z])*ln([1-x-z])*x*[1-z]^-1
       - 2*ln([1-z])*ln([1-x-z])*x
       + Li2(x^-1*[1-x]^-1*z*[1-z])
       - 1/2*Li2(x^-1*[1-x]^-1*z*[1-z])*z^-1
       - 1/2*Li2(x^-1*[1-x]^-1*z*[1-z])*[1-z]^-1
       + Li2(x^-1*[1-x]^-1*z*[1-z])*x*z^-1
       + Li2(x^-1*[1-x]^-1*z*[1-z])*x*[1-z]^-1
       - 2*Li2(x^-1*[1-x]^-1*z*[1-z])*x
       - Li2([1-x]^-1*z)
       + 1/2*Li2([1-x]^-1*z)*z^-1
       + 1/2*Li2([1-x]^-1*z)*[1-z]^-1
       - Li2([1-x]^-1*z)*x*z^-1
       - Li2([1-x]^-1*z)*x*[1-z]^-1
       + 2*Li2([1-x]^-1*z)*x
       + Li2([1-x]*[1-z]^-1)
       - 1/2*Li2([1-x]*[1-z]^-1)*z^-1
       - 1/2*Li2([1-x]*[1-z]^-1)*[1-z]^-1
       + Li2([1-x]*[1-z]^-1)*x*z^-1
       + Li2([1-x]*[1-z]^-1)*x*[1-z]^-1
       - 2*Li2([1-x]*[1-z]^-1)*x
       )

       + T(u2)*NC^-1 * (  - 1
       + 2*z^-1
       + 2*[1-z]^-1
       + x*z^-1
       + x*[1-z]^-1
       - 4*x
       + 11/12*pi^2*z^-1
       + 5/6*pi^2*[1-z]^-1
       - 7/12*pi^2
       - 11/6*pi^2*x*z^-1
       - 5/3*pi^2*x*[1-z]^-1
       + 7/6*pi^2*x
       + 2/3*pi^2*x^2*z^-1
       + 1/2*pi^2*x^2*[1-z]^-1
       - 6*ln(x)
       + 9/2*ln(x)*z^-1
       + 9/2*ln(x)*[1-z]^-1
       - 12*ln(x)*x*z^-1
       - 12*ln(x)*x*[1-z]^-1
       + 18*ln(x)*x
       + 1/2*ln(x)*ln( - [1-x-z])
       - ln(x)*ln( - [1-x-z])*z^-1
       - 1/2*ln(x)*ln( - [1-x-z])*[1-z]^-1
       + 2*ln(x)*ln( - [1-x-z])*x*z^-1
       + ln(x)*ln( - [1-x-z])*x*[1-z]^-1
       - ln(x)*ln( - [1-x-z])*x
       - ln(x)*ln( - [1-x-z])*x^2*z^-1
       + 3*ln(x)^2
       - 9/2*ln(x)^2*z^-1
       - 9/2*ln(x)^2*[1-z]^-1
       + 9*ln(x)^2*x*z^-1
       + 9*ln(x)^2*x*[1-z]^-1
       - 6*ln(x)^2*x
       - 3*ln(x)^2*x^2*z^-1
       - 3*ln(x)^2*x^2*[1-z]^-1
       - 5*ln(x)*ln([1-x])
       + 15/2*ln(x)*ln([1-x])*z^-1
       + 15/2*ln(x)*ln([1-x])*[1-z]^-1
       - 15*ln(x)*ln([1-x])*x*z^-1
       - 15*ln(x)*ln([1-x])*x*[1-z]^-1
       + 10*ln(x)*ln([1-x])*x
       + 5*ln(x)*ln([1-x])*x^2*z^-1
       + 5*ln(x)*ln([1-x])*x^2*[1-z]^-1
       - 9/2*ln(x)*ln(z)
       + 13/2*ln(x)*ln(z)*z^-1
       + 7*ln(x)*ln(z)*[1-z]^-1
       - 13*ln(x)*ln(z)*x*z^-1
       - 14*ln(x)*ln(z)*x*[1-z]^-1
       + 9*ln(x)*ln(z)*x
       + 4*ln(x)*ln(z)*x^2*z^-1
       + 5*ln(x)*ln(z)*x^2*[1-z]^-1
       - 9/2*ln(x)*ln([1-z])
       + 7*ln(x)*ln([1-z])*z^-1
       + 13/2*ln(x)*ln([1-z])*[1-z]^-1
       - 14*ln(x)*ln([1-z])*x*z^-1
       - 13*ln(x)*ln([1-z])*x*[1-z]^-1
       + 9*ln(x)*ln([1-z])*x
       + 5*ln(x)*ln([1-z])*x^2*z^-1
       + 4*ln(x)*ln([1-z])*x^2*[1-z]^-1
       + 1/2*ln(x)*ln([x-z])
       - 1/2*ln(x)*ln([x-z])*z^-1
       - ln(x)*ln([x-z])*[1-z]^-1
       + ln(x)*ln([x-z])*x*z^-1
       + 2*ln(x)*ln([x-z])*x*[1-z]^-1
       - ln(x)*ln([x-z])*x
       - ln(x)*ln([x-z])*x^2*[1-z]^-1
       + 4*ln([1-x])
       - 3*ln([1-x])*z^-1
       - 3*ln([1-x])*[1-z]^-1
       + 8*ln([1-x])*x*z^-1
       + 8*ln([1-x])*x*[1-z]^-1
       - 12*ln([1-x])*x
       + 2*ln([1-x])^2
       - 3*ln([1-x])^2*z^-1
       - 3*ln([1-x])^2*[1-z]^-1
       + 6*ln([1-x])^2*x*z^-1
       + 6*ln([1-x])^2*x*[1-z]^-1
       - 4*ln([1-x])^2*x
       - 2*ln([1-x])^2*x^2*z^-1
       - 2*ln([1-x])^2*x^2*[1-z]^-1
       + 3*ln([1-x])*ln(z)
       - 4*ln([1-x])*ln(z)*z^-1
       - 5*ln([1-x])*ln(z)*[1-z]^-1
       + 8*ln([1-x])*ln(z)*x*z^-1
       + 10*ln([1-x])*ln(z)*x*[1-z]^-1
       - 6*ln([1-x])*ln(z)*x
       - 2*ln([1-x])*ln(z)*x^2*z^-1
       - 4*ln([1-x])*ln(z)*x^2*[1-z]^-1
       + 3*ln([1-x])*ln([1-z])
       - 5*ln([1-x])*ln([1-z])*z^-1
       - 4*ln([1-x])*ln([1-z])*[1-z]^-1
       + 10*ln([1-x])*ln([1-z])*x*z^-1
       + 8*ln([1-x])*ln([1-z])*x*[1-z]^-1
       - 6*ln([1-x])*ln([1-z])*x
       - 4*ln([1-x])*ln([1-z])*x^2*z^-1
       - 2*ln([1-x])*ln([1-z])*x^2*[1-z]^-1
       + 3*ln(z)
       - 2*ln(z)*z^-1
       - 5/2*ln(z)*[1-z]^-1
       + 6*ln(z)*x*z^-1
       + 6*ln(z)*x*[1-z]^-1
       - 9*ln(z)*x
       + 5/4*ln(z)^2
       - 3/2*ln(z)^2*z^-1
       - 9/4*ln(z)^2*[1-z]^-1
       + 3*ln(z)^2*x*z^-1
       + 9/2*ln(z)^2*x*[1-z]^-1
       - 5/2*ln(z)^2*x
       - 1/2*ln(z)^2*x^2*z^-1
       - 2*ln(z)^2*x^2*[1-z]^-1
       + 5/2*ln(z)*ln([1-z])
       - 4*ln(z)*ln([1-z])*z^-1
       - 7/2*ln(z)*ln([1-z])*[1-z]^-1
       + 8*ln(z)*ln([1-z])*x*z^-1
       + 7*ln(z)*ln([1-z])*x*[1-z]^-1
       - 5*ln(z)*ln([1-z])*x
       - 3*ln(z)*ln([1-z])*x^2*z^-1
       - 2*ln(z)*ln([1-z])*x^2*[1-z]^-1
       - 1/2*ln(z)*ln([x-z])
       + 1/2*ln(z)*ln([x-z])*z^-1
       + ln(z)*ln([x-z])*[1-z]^-1
       - ln(z)*ln([x-z])*x*z^-1
       - 2*ln(z)*ln([x-z])*x*[1-z]^-1
       + ln(z)*ln([x-z])*x
       + ln(z)*ln([x-z])*x^2*[1-z]^-1
       + 3*ln([1-z])
       - 5/2*ln([1-z])*z^-1
       - 2*ln([1-z])*[1-z]^-1
       + 6*ln([1-z])*x*z^-1
       + 6*ln([1-z])*x*[1-z]^-1
       - 9*ln([1-z])*x
       - 1/2*ln([1-z])*ln( - [1-x-z])
       + ln([1-z])*ln( - [1-x-z])*z^-1
       + 1/2*ln([1-z])*ln( - [1-x-z])*[1-z]^-1
       - 2*ln([1-z])*ln( - [1-x-z])*x*z^-1
       - ln([1-z])*ln( - [1-x-z])*x*[1-z]^-1
       + ln([1-z])*ln( - [1-x-z])*x
       + ln([1-z])*ln( - [1-x-z])*x^2*z^-1
       + 5/4*ln([1-z])^2
       - 9/4*ln([1-z])^2*z^-1
       - 3/2*ln([1-z])^2*[1-z]^-1
       + 9/2*ln([1-z])^2*x*z^-1
       + 3*ln([1-z])^2*x*[1-z]^-1
       - 5/2*ln([1-z])^2*x
       - 2*ln([1-z])^2*x^2*z^-1
       - 1/2*ln([1-z])^2*x^2*[1-z]^-1
       - 1/2*Li2(x^-1*[1-x]*z^-1*[1-z])
       + Li2(x^-1*[1-x]*z^-1*[1-z])*z^-1
       + 1/2*Li2(x^-1*[1-x]*z^-1*[1-z])*[1-z]^-1
       - 2*Li2(x^-1*[1-x]*z^-1*[1-z])*x*z^-1
       - Li2(x^-1*[1-x]*z^-1*[1-z])*x*[1-z]^-1
       + Li2(x^-1*[1-x]*z^-1*[1-z])*x
       + Li2(x^-1*[1-x]*z^-1*[1-z])*x^2*z^-1
       - 1/2*Li2(x^-1*[1-x]*z*[1-z]^-1)
       + 1/2*Li2(x^-1*[1-x]*z*[1-z]^-1)*z^-1
       + Li2(x^-1*[1-x]*z*[1-z]^-1)*[1-z]^-1
       - Li2(x^-1*[1-x]*z*[1-z]^-1)*x*z^-1
       - 2*Li2(x^-1*[1-x]*z*[1-z]^-1)*x*[1-z]^-1
       + Li2(x^-1*[1-x]*z*[1-z]^-1)*x
       + Li2(x^-1*[1-x]*z*[1-z]^-1)*x^2*[1-z]^-1
       + 1/2*Li2([1-x]*z^-1)
       - Li2([1-x]*z^-1)*z^-1
       - 1/2*Li2([1-x]*z^-1)*[1-z]^-1
       + 2*Li2([1-x]*z^-1)*x*z^-1
       + Li2([1-x]*z^-1)*x*[1-z]^-1
       - Li2([1-x]*z^-1)*x
       - Li2([1-x]*z^-1)*x^2*z^-1
       + 1/2*Li2([1-x]*[1-z]^-1)
       - 1/2*Li2([1-x]*[1-z]^-1)*z^-1
       - Li2([1-x]*[1-z]^-1)*[1-z]^-1
       + Li2([1-x]*[1-z]^-1)*x*z^-1
       + 2*Li2([1-x]*[1-z]^-1)*x*[1-z]^-1
       - Li2([1-x]*[1-z]^-1)*x
       - Li2([1-x]*[1-z]^-1)*x^2*[1-z]^-1
       - 1/2*Li2(z)*z^-1
       + 1/2*Li2(z)*[1-z]^-1
       + Li2(z)*x*z^-1
       - Li2(z)*x*[1-z]^-1
       - Li2(z)*x^2*z^-1
       + Li2(z)*x^2*[1-z]^-1
       )

       + T(u2)*NC * ( 2
       - z^-1
       - [1-z]^-1
       - 1/12*pi^2*z^-1
       - 1/12*pi^2*[1-z]^-1
       + 1/6*pi^2
       + 1/6*pi^2*x*z^-1
       + 1/6*pi^2*x*[1-z]^-1
       - 1/3*pi^2*x
       + 3*ln(x)
       - 3/2*ln(x)*z^-1
       - 3/2*ln(x)*[1-z]^-1
       + 6*ln(x)*x*z^-1
       + 6*ln(x)*x*[1-z]^-1
       - 18*ln(x)*x
       - ln(x)*ln( - [1-x-z])
       + 1/2*ln(x)*ln( - [1-x-z])*z^-1
       + 1/2*ln(x)*ln( - [1-x-z])*[1-z]^-1
       - ln(x)*ln( - [1-x-z])*x*z^-1
       - ln(x)*ln( - [1-x-z])*x*[1-z]^-1
       + 2*ln(x)*ln( - [1-x-z])*x
       - 7/2*ln(x)^2
       + 7/4*ln(x)^2*z^-1
       + 7/4*ln(x)^2*[1-z]^-1
       - 7/2*ln(x)^2*x*z^-1
       - 7/2*ln(x)^2*x*[1-z]^-1
       + 7*ln(x)^2*x
       + 3*ln(x)*ln([1-x])
       - 3/2*ln(x)*ln([1-x])*z^-1
       - 3/2*ln(x)*ln([1-x])*[1-z]^-1
       + 3*ln(x)*ln([1-x])*x*z^-1
       + 3*ln(x)*ln([1-x])*x*[1-z]^-1
       - 6*ln(x)*ln([1-x])*x
       + 6*ln(x)*ln(z)
       - 3*ln(x)*ln(z)*z^-1
       - 3*ln(x)*ln(z)*[1-z]^-1
       + 6*ln(x)*ln(z)*x*z^-1
       + 6*ln(x)*ln(z)*x*[1-z]^-1
       - 12*ln(x)*ln(z)*x
       + 6*ln(x)*ln([1-z])
       - 3*ln(x)*ln([1-z])*z^-1
       - 3*ln(x)*ln([1-z])*[1-z]^-1
       + 6*ln(x)*ln([1-z])*x*z^-1
       + 6*ln(x)*ln([1-z])*x*[1-z]^-1
       - 12*ln(x)*ln([1-z])*x
       - ln(x)*ln([x-z])
       + 1/2*ln(x)*ln([x-z])*z^-1
       + 1/2*ln(x)*ln([x-z])*[1-z]^-1
       - ln(x)*ln([x-z])*x*z^-1
       - ln(x)*ln([x-z])*x*[1-z]^-1
       + 2*ln(x)*ln([x-z])*x
       - ln([1-x])
       + 1/2*ln([1-x])*z^-1
       + 1/2*ln([1-x])*[1-z]^-1
       - 2*ln([1-x])*x*z^-1
       - 2*ln([1-x])*x*[1-z]^-1
       + 6*ln([1-x])*x
       - 1/2*ln([1-x])^2
       + 1/4*ln([1-x])^2*z^-1
       + 1/4*ln([1-x])^2*[1-z]^-1
       - 1/2*ln([1-x])^2*x*z^-1
       - 1/2*ln([1-x])^2*x*[1-z]^-1
       + ln([1-x])^2*x
       - 2*ln([1-x])*ln(z)
       + ln([1-x])*ln(z)*z^-1
       + ln([1-x])*ln(z)*[1-z]^-1
       - 2*ln([1-x])*ln(z)*x*z^-1
       - 2*ln([1-x])*ln(z)*x*[1-z]^-1
       + 4*ln([1-x])*ln(z)*x
       - 2*ln([1-x])*ln([1-z])
       + ln([1-x])*ln([1-z])*z^-1
       + ln([1-x])*ln([1-z])*[1-z]^-1
       - 2*ln([1-x])*ln([1-z])*x*z^-1
       - 2*ln([1-x])*ln([1-z])*x*[1-z]^-1
       + 4*ln([1-x])*ln([1-z])*x
       - 2*ln(z)
       + ln(z)*z^-1
       + ln(z)*[1-z]^-1
       - 4*ln(z)*x*z^-1
       - 4*ln(z)*x*[1-z]^-1
       + 12*ln(z)*x
       - 2*ln(z)^2
       + ln(z)^2*z^-1
       + ln(z)^2*[1-z]^-1
       - 2*ln(z)^2*x*z^-1
       - 2*ln(z)^2*x*[1-z]^-1
       + 4*ln(z)^2*x
       - 5*ln(z)*ln([1-z])
       + 5/2*ln(z)*ln([1-z])*z^-1
       + 5/2*ln(z)*ln([1-z])*[1-z]^-1
       - 5*ln(z)*ln([1-z])*x*z^-1
       - 5*ln(z)*ln([1-z])*x*[1-z]^-1
       + 10*ln(z)*ln([1-z])*x
       + ln(z)*ln([x-z])
       - 1/2*ln(z)*ln([x-z])*z^-1
       - 1/2*ln(z)*ln([x-z])*[1-z]^-1
       + ln(z)*ln([x-z])*x*z^-1
       + ln(z)*ln([x-z])*x*[1-z]^-1
       - 2*ln(z)*ln([x-z])*x
       - 2*ln([1-z])
       + ln([1-z])*z^-1
       + ln([1-z])*[1-z]^-1
       - 4*ln([1-z])*x*z^-1
       - 4*ln([1-z])*x*[1-z]^-1
       + 12*ln([1-z])*x
       + ln([1-z])*ln( - [1-x-z])
       - 1/2*ln([1-z])*ln( - [1-x-z])*z^-1
       - 1/2*ln([1-z])*ln( - [1-x-z])*[1-z]^-1
       + ln([1-z])*ln( - [1-x-z])*x*z^-1
       + ln([1-z])*ln( - [1-x-z])*x*[1-z]^-1
       - 2*ln([1-z])*ln( - [1-x-z])*x
       - 2*ln([1-z])^2
       + ln([1-z])^2*z^-1
       + ln([1-z])^2*[1-z]^-1
       - 2*ln([1-z])^2*x*z^-1
       - 2*ln([1-z])^2*x*[1-z]^-1
       + 4*ln([1-z])^2*x
       + Li2([1-x]*z^-1)
       - 1/2*Li2([1-x]*z^-1)*z^-1
       - 1/2*Li2([1-x]*z^-1)*[1-z]^-1
       + Li2([1-x]*z^-1)*x*z^-1
       + Li2([1-x]*z^-1)*x*[1-z]^-1
       - 2*Li2([1-x]*z^-1)*x
       + Li2([1-x]*[1-z]^-1)
       - 1/2*Li2([1-x]*[1-z]^-1)*z^-1
       - 1/2*Li2([1-x]*[1-z]^-1)*[1-z]^-1
       + Li2([1-x]*[1-z]^-1)*x*z^-1
       + Li2([1-x]*[1-z]^-1)*x*[1-z]^-1
       - 2*Li2([1-x]*[1-z]^-1)*x
       - Li2(x*[1-x]*z^-1*[1-z]^-1)
       + 1/2*Li2(x*[1-x]*z^-1*[1-z]^-1)*z^-1
       + 1/2*Li2(x*[1-x]*z^-1*[1-z]^-1)*[1-z]^-1
       - Li2(x*[1-x]*z^-1*[1-z]^-1)*x*z^-1
       - Li2(x*[1-x]*z^-1*[1-z]^-1)*x*[1-z]^-1
       + 2*Li2(x*[1-x]*z^-1*[1-z]^-1)*x
       )

       + T(u3)*NC^-1 * (  - 1
       + 2*z^-1
       + 2*[1-z]^-1
       + x*z^-1
       + x*[1-z]^-1
       - 4*x
       + 11/12*pi^2*z^-1
       + 5/6*pi^2*[1-z]^-1
       - 7/12*pi^2
       - 11/6*pi^2*x*z^-1
       - 5/3*pi^2*x*[1-z]^-1
       + 7/6*pi^2*x
       + 2/3*pi^2*x^2*z^-1
       + 1/2*pi^2*x^2*[1-z]^-1
       - 6*ln(x)
       + 9/2*ln(x)*z^-1
       + 9/2*ln(x)*[1-z]^-1
       - 12*ln(x)*x*z^-1
       - 12*ln(x)*x*[1-z]^-1
       + 18*ln(x)*x
       + 1/2*ln(x)*ln( - [x-z])
       - 1/2*ln(x)*ln( - [x-z])*z^-1
       - ln(x)*ln( - [x-z])*[1-z]^-1
       + ln(x)*ln( - [x-z])*x*z^-1
       + 2*ln(x)*ln( - [x-z])*x*[1-z]^-1
       - ln(x)*ln( - [x-z])*x
       - ln(x)*ln( - [x-z])*x^2*[1-z]^-1
       + 7/2*ln(x)^2
       - 21/4*ln(x)^2*z^-1
       - 21/4*ln(x)^2*[1-z]^-1
       + 21/2*ln(x)^2*x*z^-1
       + 21/2*ln(x)^2*x*[1-z]^-1
       - 7*ln(x)^2*x
       - 7/2*ln(x)^2*x^2*z^-1
       - 7/2*ln(x)^2*x^2*[1-z]^-1
       - 6*ln(x)*ln([1-x])
       + 9*ln(x)*ln([1-x])*z^-1
       + 9*ln(x)*ln([1-x])*[1-z]^-1
       - 18*ln(x)*ln([1-x])*x*z^-1
       - 18*ln(x)*ln([1-x])*x*[1-z]^-1
       + 12*ln(x)*ln([1-x])*x
       + 6*ln(x)*ln([1-x])*x^2*z^-1
       + 6*ln(x)*ln([1-x])*x^2*[1-z]^-1
       - 9/2*ln(x)*ln(z)
       + 6*ln(x)*ln(z)*z^-1
       + 15/2*ln(x)*ln(z)*[1-z]^-1
       - 12*ln(x)*ln(z)*x*z^-1
       - 15*ln(x)*ln(z)*x*[1-z]^-1
       + 9*ln(x)*ln(z)*x
       + 3*ln(x)*ln(z)*x^2*z^-1
       + 6*ln(x)*ln(z)*x^2*[1-z]^-1
       - 9/2*ln(x)*ln([1-z])
       + 15/2*ln(x)*ln([1-z])*z^-1
       + 6*ln(x)*ln([1-z])*[1-z]^-1
       - 15*ln(x)*ln([1-z])*x*z^-1
       - 12*ln(x)*ln([1-z])*x*[1-z]^-1
       + 9*ln(x)*ln([1-z])*x
       + 6*ln(x)*ln([1-z])*x^2*z^-1
       + 3*ln(x)*ln([1-z])*x^2*[1-z]^-1
       + 1/2*ln(x)*ln([1-x-z])
       - ln(x)*ln([1-x-z])*z^-1
       - 1/2*ln(x)*ln([1-x-z])*[1-z]^-1
       + 2*ln(x)*ln([1-x-z])*x*z^-1
       + ln(x)*ln([1-x-z])*x*[1-z]^-1
       - ln(x)*ln([1-x-z])*x
       - ln(x)*ln([1-x-z])*x^2*z^-1
       + 4*ln([1-x])
       - 3*ln([1-x])*z^-1
       - 3*ln([1-x])*[1-z]^-1
       + 8*ln([1-x])*x*z^-1
       + 8*ln([1-x])*x*[1-z]^-1
       - 12*ln([1-x])*x
       + 2*ln([1-x])^2
       - 3*ln([1-x])^2*z^-1
       - 3*ln([1-x])^2*[1-z]^-1
       + 6*ln([1-x])^2*x*z^-1
       + 6*ln([1-x])^2*x*[1-z]^-1
       - 4*ln([1-x])^2*x
       - 2*ln([1-x])^2*x^2*z^-1
       - 2*ln([1-x])^2*x^2*[1-z]^-1
       + 7/2*ln([1-x])*ln(z)
       - 9/2*ln([1-x])*ln(z)*z^-1
       - 6*ln([1-x])*ln(z)*[1-z]^-1
       + 9*ln([1-x])*ln(z)*x*z^-1
       + 12*ln([1-x])*ln(z)*x*[1-z]^-1
       - 7*ln([1-x])*ln(z)*x
       - 2*ln([1-x])*ln(z)*x^2*z^-1
       - 5*ln([1-x])*ln(z)*x^2*[1-z]^-1
       + 7/2*ln([1-x])*ln([1-z])
       - 6*ln([1-x])*ln([1-z])*z^-1
       - 9/2*ln([1-x])*ln([1-z])*[1-z]^-1
       + 12*ln([1-x])*ln([1-z])*x*z^-1
       + 9*ln([1-x])*ln([1-z])*x*[1-z]^-1
       - 7*ln([1-x])*ln([1-z])*x
       - 5*ln([1-x])*ln([1-z])*x^2*z^-1
       - 2*ln([1-x])*ln([1-z])*x^2*[1-z]^-1
       + 3*ln(z)
       - 2*ln(z)*z^-1
       - 5/2*ln(z)*[1-z]^-1
       + 6*ln(z)*x*z^-1
       + 6*ln(z)*x*[1-z]^-1
       - 9*ln(z)*x
       - 1/2*ln(z)*ln( - [x-z])
       + 1/2*ln(z)*ln( - [x-z])*z^-1
       + ln(z)*ln( - [x-z])*[1-z]^-1
       - ln(z)*ln( - [x-z])*x*z^-1
       - 2*ln(z)*ln( - [x-z])*x*[1-z]^-1
       + ln(z)*ln( - [x-z])*x
       + ln(z)*ln( - [x-z])*x^2*[1-z]^-1
       + 3/2*ln(z)^2
       - 7/4*ln(z)^2*z^-1
       - 11/4*ln(z)^2*[1-z]^-1
       + 7/2*ln(z)^2*x*z^-1
       + 11/2*ln(z)^2*x*[1-z]^-1
       - 3*ln(z)^2*x
       - 1/2*ln(z)^2*x^2*z^-1
       - 5/2*ln(z)^2*x^2*[1-z]^-1
       + 3/2*ln(z)*ln([1-z])
       - 5/2*ln(z)*ln([1-z])*z^-1
       - 2*ln(z)*ln([1-z])*[1-z]^-1
       + 5*ln(z)*ln([1-z])*x*z^-1
       + 4*ln(z)*ln([1-z])*x*[1-z]^-1
       - 3*ln(z)*ln([1-z])*x
       - 2*ln(z)*ln([1-z])*x^2*z^-1
       - ln(z)*ln([1-z])*x^2*[1-z]^-1
       + 3*ln([1-z])
       - 5/2*ln([1-z])*z^-1
       - 2*ln([1-z])*[1-z]^-1
       + 6*ln([1-z])*x*z^-1
       + 6*ln([1-z])*x*[1-z]^-1
       - 9*ln([1-z])*x
       + 3/2*ln([1-z])^2
       - 11/4*ln([1-z])^2*z^-1
       - 7/4*ln([1-z])^2*[1-z]^-1
       + 11/2*ln([1-z])^2*x*z^-1
       + 7/2*ln([1-z])^2*x*[1-z]^-1
       - 3*ln([1-z])^2*x
       - 5/2*ln([1-z])^2*x^2*z^-1
       - 1/2*ln([1-z])^2*x^2*[1-z]^-1
       - 1/2*ln([1-z])*ln([1-x-z])
       + ln([1-z])*ln([1-x-z])*z^-1
       + 1/2*ln([1-z])*ln([1-x-z])*[1-z]^-1
       - 2*ln([1-z])*ln([1-x-z])*x*z^-1
       - ln([1-z])*ln([1-x-z])*x*[1-z]^-1
       + ln([1-z])*ln([1-x-z])*x
       + ln([1-z])*ln([1-x-z])*x^2*z^-1
       - 1/2*Li2([1-x]^-1*[1-z])
       + 1/2*Li2([1-x]^-1*[1-z])*z^-1
       + Li2([1-x]^-1*[1-z])*[1-z]^-1
       - Li2([1-x]^-1*[1-z])*x*z^-1
       - 2*Li2([1-x]^-1*[1-z])*x*[1-z]^-1
       + Li2([1-x]^-1*[1-z])*x
       + Li2([1-x]^-1*[1-z])*x^2*[1-z]^-1
       - 1/2*Li2([1-x]^-1*z)
       + Li2([1-x]^-1*z)*z^-1
       + 1/2*Li2([1-x]^-1*z)*[1-z]^-1
       - 2*Li2([1-x]^-1*z)*x*z^-1
       - Li2([1-x]^-1*z)*x*[1-z]^-1
       + Li2([1-x]^-1*z)*x
       + Li2([1-x]^-1*z)*x^2*z^-1
       + 1/2*Li2(x*[1-x]^-1*z^-1*[1-z])
       - 1/2*Li2(x*[1-x]^-1*z^-1*[1-z])*z^-1
       - Li2(x*[1-x]^-1*z^-1*[1-z])*[1-z]^-1
       + Li2(x*[1-x]^-1*z^-1*[1-z])*x*z^-1
       + 2*Li2(x*[1-x]^-1*z^-1*[1-z])*x*[1-z]^-1
       - Li2(x*[1-x]^-1*z^-1*[1-z])*x
       - Li2(x*[1-x]^-1*z^-1*[1-z])*x^2*[1-z]^-1
       + 1/2*Li2(x*[1-x]^-1*z*[1-z]^-1)
       - Li2(x*[1-x]^-1*z*[1-z]^-1)*z^-1
       - 1/2*Li2(x*[1-x]^-1*z*[1-z]^-1)*[1-z]^-1
       + 2*Li2(x*[1-x]^-1*z*[1-z]^-1)*x*z^-1
       + Li2(x*[1-x]^-1*z*[1-z]^-1)*x*[1-z]^-1
       - Li2(x*[1-x]^-1*z*[1-z]^-1)*x
       - Li2(x*[1-x]^-1*z*[1-z]^-1)*x^2*z^-1
       - 1/2*Li2(z)*z^-1
       + 1/2*Li2(z)*[1-z]^-1
       + Li2(z)*x*z^-1
       - Li2(z)*x*[1-z]^-1
       - Li2(z)*x^2*z^-1
       + Li2(z)*x^2*[1-z]^-1
       )

       + T(u3)*NC * ( 2
       - z^-1
       - [1-z]^-1
       - 5/12*pi^2*z^-1
       - 5/12*pi^2*[1-z]^-1
       + 5/6*pi^2
       + 5/6*pi^2*x*z^-1
       + 5/6*pi^2*x*[1-z]^-1
       - 5/3*pi^2*x
       + 3*ln(x)
       - 3/2*ln(x)*z^-1
       - 3/2*ln(x)*[1-z]^-1
       + 6*ln(x)*x*z^-1
       + 6*ln(x)*x*[1-z]^-1
       - 18*ln(x)*x
       - ln(x)*ln( - [x-z])
       + 1/2*ln(x)*ln( - [x-z])*z^-1
       + 1/2*ln(x)*ln( - [x-z])*[1-z]^-1
       - ln(x)*ln( - [x-z])*x*z^-1
       - ln(x)*ln( - [x-z])*x*[1-z]^-1
       + 2*ln(x)*ln( - [x-z])*x
       - 7/2*ln(x)^2
       + 7/4*ln(x)^2*z^-1
       + 7/4*ln(x)^2*[1-z]^-1
       - 7/2*ln(x)^2*x*z^-1
       - 7/2*ln(x)^2*x*[1-z]^-1
       + 7*ln(x)^2*x
       + 3*ln(x)*ln([1-x])
       - 3/2*ln(x)*ln([1-x])*z^-1
       - 3/2*ln(x)*ln([1-x])*[1-z]^-1
       + 3*ln(x)*ln([1-x])*x*z^-1
       + 3*ln(x)*ln([1-x])*x*[1-z]^-1
       - 6*ln(x)*ln([1-x])*x
       + 6*ln(x)*ln(z)
       - 3*ln(x)*ln(z)*z^-1
       - 3*ln(x)*ln(z)*[1-z]^-1
       + 6*ln(x)*ln(z)*x*z^-1
       + 6*ln(x)*ln(z)*x*[1-z]^-1
       - 12*ln(x)*ln(z)*x
       + 6*ln(x)*ln([1-z])
       - 3*ln(x)*ln([1-z])*z^-1
       - 3*ln(x)*ln([1-z])*[1-z]^-1
       + 6*ln(x)*ln([1-z])*x*z^-1
       + 6*ln(x)*ln([1-z])*x*[1-z]^-1
       - 12*ln(x)*ln([1-z])*x
       - ln(x)*ln([1-x-z])
       + 1/2*ln(x)*ln([1-x-z])*z^-1
       + 1/2*ln(x)*ln([1-x-z])*[1-z]^-1
       - ln(x)*ln([1-x-z])*x*z^-1
       - ln(x)*ln([1-x-z])*x*[1-z]^-1
       + 2*ln(x)*ln([1-x-z])*x
       - ln([1-x])
       + 1/2*ln([1-x])*z^-1
       + 1/2*ln([1-x])*[1-z]^-1
       - 2*ln([1-x])*x*z^-1
       - 2*ln([1-x])*x*[1-z]^-1
       + 6*ln([1-x])*x
       - 3/2*ln([1-x])^2
       + 3/4*ln([1-x])^2*z^-1
       + 3/4*ln([1-x])^2*[1-z]^-1
       - 3/2*ln([1-x])^2*x*z^-1
       - 3/2*ln([1-x])^2*x*[1-z]^-1
       + 3*ln([1-x])^2*x
       - ln([1-x])*ln(z)
       + 1/2*ln([1-x])*ln(z)*z^-1
       + 1/2*ln([1-x])*ln(z)*[1-z]^-1
       - ln([1-x])*ln(z)*x*z^-1
       - ln([1-x])*ln(z)*x*[1-z]^-1
       + 2*ln([1-x])*ln(z)*x
       - ln([1-x])*ln([1-z])
       + 1/2*ln([1-x])*ln([1-z])*z^-1
       + 1/2*ln([1-x])*ln([1-z])*[1-z]^-1
       - ln([1-x])*ln([1-z])*x*z^-1
       - ln([1-x])*ln([1-z])*x*[1-z]^-1
       + 2*ln([1-x])*ln([1-z])*x
       - 2*ln(z)
       + ln(z)*z^-1
       + ln(z)*[1-z]^-1
       - 4*ln(z)*x*z^-1
       - 4*ln(z)*x*[1-z]^-1
       + 12*ln(z)*x
       + ln(z)*ln( - [x-z])
       - 1/2*ln(z)*ln( - [x-z])*z^-1
       - 1/2*ln(z)*ln( - [x-z])*[1-z]^-1
       + ln(z)*ln( - [x-z])*x*z^-1
       + ln(z)*ln( - [x-z])*x*[1-z]^-1
       - 2*ln(z)*ln( - [x-z])*x
       - 5/2*ln(z)^2
       + 5/4*ln(z)^2*z^-1
       + 5/4*ln(z)^2*[1-z]^-1
       - 5/2*ln(z)^2*x*z^-1
       - 5/2*ln(z)^2*x*[1-z]^-1
       + 5*ln(z)^2*x
       - 5*ln(z)*ln([1-z])
       + 5/2*ln(z)*ln([1-z])*z^-1
       + 5/2*ln(z)*ln([1-z])*[1-z]^-1
       - 5*ln(z)*ln([1-z])*x*z^-1
       - 5*ln(z)*ln([1-z])*x*[1-z]^-1
       + 10*ln(z)*ln([1-z])*x
       - 2*ln([1-z])
       + ln([1-z])*z^-1
       + ln([1-z])*[1-z]^-1
       - 4*ln([1-z])*x*z^-1
       - 4*ln([1-z])*x*[1-z]^-1
       + 12*ln([1-z])*x
       - 5/2*ln([1-z])^2
       + 5/4*ln([1-z])^2*z^-1
       + 5/4*ln([1-z])^2*[1-z]^-1
       - 5/2*ln([1-z])^2*x*z^-1
       - 5/2*ln([1-z])^2*x*[1-z]^-1
       + 5*ln([1-z])^2*x
       + ln([1-z])*ln([1-x-z])
       - 1/2*ln([1-z])*ln([1-x-z])*z^-1
       - 1/2*ln([1-z])*ln([1-x-z])*[1-z]^-1
       + ln([1-z])*ln([1-x-z])*x*z^-1
       + ln([1-z])*ln([1-x-z])*x*[1-z]^-1
       - 2*ln([1-z])*ln([1-x-z])*x
       - Li2([1-x]^-1*[1-z])
       + 1/2*Li2([1-x]^-1*[1-z])*z^-1
       + 1/2*Li2([1-x]^-1*[1-z])*[1-z]^-1
       - Li2([1-x]^-1*[1-z])*x*z^-1
       - Li2([1-x]^-1*[1-z])*x*[1-z]^-1
       + 2*Li2([1-x]^-1*[1-z])*x
       - Li2([1-x]^-1*z)
       + 1/2*Li2([1-x]^-1*z)*z^-1
       + 1/2*Li2([1-x]^-1*z)*[1-z]^-1
       - Li2([1-x]^-1*z)*x*z^-1
       - Li2([1-x]^-1*z)*x*[1-z]^-1
       + 2*Li2([1-x]^-1*z)*x
       - Li2(x*[1-x]*z^-1*[1-z]^-1)
       + 1/2*Li2(x*[1-x]*z^-1*[1-z]^-1)*z^-1
       + 1/2*Li2(x*[1-x]*z^-1*[1-z]^-1)*[1-z]^-1
       - Li2(x*[1-x]*z^-1*[1-z]^-1)*x*z^-1
       - Li2(x*[1-x]*z^-1*[1-z]^-1)*x*[1-z]^-1
       + 2*Li2(x*[1-x]*z^-1*[1-z]^-1)*x
       )

       + T(u4)*NC^-1 * (  - 1
       + 2*z^-1
       + 2*[1-z]^-1
       + x*z^-1
       + x*[1-z]^-1
       - 4*x
       + 11/12*pi^2*z^-1
       + 5/6*pi^2*[1-z]^-1
       - 7/12*pi^2
       - 11/6*pi^2*x*z^-1
       - 5/3*pi^2*x*[1-z]^-1
       + 7/6*pi^2*x
       + 2/3*pi^2*x^2*z^-1
       + 1/2*pi^2*x^2*[1-z]^-1
       - 6*ln(x)
       + 9/2*ln(x)*z^-1
       + 9/2*ln(x)*[1-z]^-1
       - 12*ln(x)*x*z^-1
       - 12*ln(x)*x*[1-z]^-1
       + 18*ln(x)*x
       + 1/2*ln(x)*ln( - [x-z])
       - 1/2*ln(x)*ln( - [x-z])*z^-1
       - ln(x)*ln( - [x-z])*[1-z]^-1
       + ln(x)*ln( - [x-z])*x*z^-1
       + 2*ln(x)*ln( - [x-z])*x*[1-z]^-1
       - ln(x)*ln( - [x-z])*x
       - ln(x)*ln( - [x-z])*x^2*[1-z]^-1
       + 1/2*ln(x)*ln( - [1-x-z])
       - ln(x)*ln( - [1-x-z])*z^-1
       - 1/2*ln(x)*ln( - [1-x-z])*[1-z]^-1
       + 2*ln(x)*ln( - [1-x-z])*x*z^-1
       + ln(x)*ln( - [1-x-z])*x*[1-z]^-1
       - ln(x)*ln( - [1-x-z])*x
       - ln(x)*ln( - [1-x-z])*x^2*z^-1
       + 13/4*ln(x)^2
       - 19/4*ln(x)^2*z^-1
       - 5*ln(x)^2*[1-z]^-1
       + 19/2*ln(x)^2*x*z^-1
       + 10*ln(x)^2*x*[1-z]^-1
       - 13/2*ln(x)^2*x
       - 3*ln(x)^2*x^2*z^-1
       - 7/2*ln(x)^2*x^2*[1-z]^-1
       - 11/2*ln(x)*ln([1-x])
       + 8*ln(x)*ln([1-x])*z^-1
       + 17/2*ln(x)*ln([1-x])*[1-z]^-1
       - 16*ln(x)*ln([1-x])*x*z^-1
       - 17*ln(x)*ln([1-x])*x*[1-z]^-1
       + 11*ln(x)*ln([1-x])*x
       + 5*ln(x)*ln([1-x])*x^2*z^-1
       + 6*ln(x)*ln([1-x])*x^2*[1-z]^-1
       - 5*ln(x)*ln(z)
       + 7*ln(x)*ln(z)*z^-1
       + 8*ln(x)*ln(z)*[1-z]^-1
       - 14*ln(x)*ln(z)*x*z^-1
       - 16*ln(x)*ln(z)*x*[1-z]^-1
       + 10*ln(x)*ln(z)*x
       + 4*ln(x)*ln(z)*x^2*z^-1
       + 6*ln(x)*ln(z)*x^2*[1-z]^-1
       - 4*ln(x)*ln([1-z])
       + 13/2*ln(x)*ln([1-z])*z^-1
       + 11/2*ln(x)*ln([1-z])*[1-z]^-1
       - 13*ln(x)*ln([1-z])*x*z^-1
       - 11*ln(x)*ln([1-z])*x*[1-z]^-1
       + 8*ln(x)*ln([1-z])*x
       + 5*ln(x)*ln([1-z])*x^2*z^-1
       + 3*ln(x)*ln([1-z])*x^2*[1-z]^-1
       + 4*ln([1-x])
       - 3*ln([1-x])*z^-1
       - 3*ln([1-x])*[1-z]^-1
       + 8*ln([1-x])*x*z^-1
       + 8*ln([1-x])*x*[1-z]^-1
       - 12*ln([1-x])*x
       + 2*ln([1-x])^2
       - 3*ln([1-x])^2*z^-1
       - 3*ln([1-x])^2*[1-z]^-1
       + 6*ln([1-x])^2*x*z^-1
       + 6*ln([1-x])^2*x*[1-z]^-1
       - 4*ln([1-x])^2*x
       - 2*ln([1-x])^2*x^2*z^-1
       - 2*ln([1-x])^2*x^2*[1-z]^-1
       + 7/2*ln([1-x])*ln(z)
       - 9/2*ln([1-x])*ln(z)*z^-1
       - 6*ln([1-x])*ln(z)*[1-z]^-1
       + 9*ln([1-x])*ln(z)*x*z^-1
       + 12*ln([1-x])*ln(z)*x*[1-z]^-1
       - 7*ln([1-x])*ln(z)*x
       - 2*ln([1-x])*ln(z)*x^2*z^-1
       - 5*ln([1-x])*ln(z)*x^2*[1-z]^-1
       + 3*ln([1-x])*ln([1-z])
       - 5*ln([1-x])*ln([1-z])*z^-1
       - 4*ln([1-x])*ln([1-z])*[1-z]^-1
       + 10*ln([1-x])*ln([1-z])*x*z^-1
       + 8*ln([1-x])*ln([1-z])*x*[1-z]^-1
       - 6*ln([1-x])*ln([1-z])*x
       - 4*ln([1-x])*ln([1-z])*x^2*z^-1
       - 2*ln([1-x])*ln([1-z])*x^2*[1-z]^-1
       + 3*ln(z)
       - 2*ln(z)*z^-1
       - 5/2*ln(z)*[1-z]^-1
       + 6*ln(z)*x*z^-1
       + 6*ln(z)*x*[1-z]^-1
       - 9*ln(z)*x
       - 1/2*ln(z)*ln( - [x-z])
       + 1/2*ln(z)*ln( - [x-z])*z^-1
       + ln(z)*ln( - [x-z])*[1-z]^-1
       - ln(z)*ln( - [x-z])*x*z^-1
       - 2*ln(z)*ln( - [x-z])*x*[1-z]^-1
       + ln(z)*ln( - [x-z])*x
       + ln(z)*ln( - [x-z])*x^2*[1-z]^-1
       + 3/2*ln(z)^2
       - 7/4*ln(z)^2*z^-1
       - 11/4*ln(z)^2*[1-z]^-1
       + 7/2*ln(z)^2*x*z^-1
       + 11/2*ln(z)^2*x*[1-z]^-1
       - 3*ln(z)^2*x
       - 1/2*ln(z)^2*x^2*z^-1
       - 5/2*ln(z)^2*x^2*[1-z]^-1
       + 2*ln(z)*ln([1-z])
       - 7/2*ln(z)*ln([1-z])*z^-1
       - 5/2*ln(z)*ln([1-z])*[1-z]^-1
       + 7*ln(z)*ln([1-z])*x*z^-1
       + 5*ln(z)*ln([1-z])*x*[1-z]^-1
       - 4*ln(z)*ln([1-z])*x
       - 3*ln(z)*ln([1-z])*x^2*z^-1
       - ln(z)*ln([1-z])*x^2*[1-z]^-1
       + 3*ln([1-z])
       - 5/2*ln([1-z])*z^-1
       - 2*ln([1-z])*[1-z]^-1
       + 6*ln([1-z])*x*z^-1
       + 6*ln([1-z])*x*[1-z]^-1
       - 9*ln([1-z])*x
       - 1/2*ln([1-z])*ln( - [1-x-z])
       + ln([1-z])*ln( - [1-x-z])*z^-1
       + 1/2*ln([1-z])*ln( - [1-x-z])*[1-z]^-1
       - 2*ln([1-z])*ln( - [1-x-z])*x*z^-1
       - ln([1-z])*ln( - [1-x-z])*x*[1-z]^-1
       + ln([1-z])*ln( - [1-x-z])*x
       + ln([1-z])*ln( - [1-x-z])*x^2*z^-1
       + 5/4*ln([1-z])^2
       - 9/4*ln([1-z])^2*z^-1
       - 3/2*ln([1-z])^2*[1-z]^-1
       + 9/2*ln([1-z])^2*x*z^-1
       + 3*ln([1-z])^2*x*[1-z]^-1
       - 5/2*ln([1-z])^2*x
       - 2*ln([1-z])^2*x^2*z^-1
       - 1/2*ln([1-z])^2*x^2*[1-z]^-1
       - 1/2*Li2(x^-1*[1-x]*z^-1*[1-z])
       + Li2(x^-1*[1-x]*z^-1*[1-z])*z^-1
       + 1/2*Li2(x^-1*[1-x]*z^-1*[1-z])*[1-z]^-1
       - 2*Li2(x^-1*[1-x]*z^-1*[1-z])*x*z^-1
       - Li2(x^-1*[1-x]*z^-1*[1-z])*x*[1-z]^-1
       + Li2(x^-1*[1-x]*z^-1*[1-z])*x
       + Li2(x^-1*[1-x]*z^-1*[1-z])*x^2*z^-1
       - 1/2*Li2([1-x]^-1*[1-z])
       + 1/2*Li2([1-x]^-1*[1-z])*z^-1
       + Li2([1-x]^-1*[1-z])*[1-z]^-1
       - Li2([1-x]^-1*[1-z])*x*z^-1
       - 2*Li2([1-x]^-1*[1-z])*x*[1-z]^-1
       + Li2([1-x]^-1*[1-z])*x
       + Li2([1-x]^-1*[1-z])*x^2*[1-z]^-1
       + 1/2*Li2([1-x]*z^-1)
       - Li2([1-x]*z^-1)*z^-1
       - 1/2*Li2([1-x]*z^-1)*[1-z]^-1
       + 2*Li2([1-x]*z^-1)*x*z^-1
       + Li2([1-x]*z^-1)*x*[1-z]^-1
       - Li2([1-x]*z^-1)*x
       - Li2([1-x]*z^-1)*x^2*z^-1
       + 1/2*Li2(x*[1-x]^-1*z^-1*[1-z])
       - 1/2*Li2(x*[1-x]^-1*z^-1*[1-z])*z^-1
       - Li2(x*[1-x]^-1*z^-1*[1-z])*[1-z]^-1
       + Li2(x*[1-x]^-1*z^-1*[1-z])*x*z^-1
       + 2*Li2(x*[1-x]^-1*z^-1*[1-z])*x*[1-z]^-1
       - Li2(x*[1-x]^-1*z^-1*[1-z])*x
       - Li2(x*[1-x]^-1*z^-1*[1-z])*x^2*[1-z]^-1
       - 1/2*Li2(z)*z^-1
       + 1/2*Li2(z)*[1-z]^-1
       + Li2(z)*x*z^-1
       - Li2(z)*x*[1-z]^-1
       - Li2(z)*x^2*z^-1
       + Li2(z)*x^2*[1-z]^-1
       )

       + T(u4)*NC * ( 2
       - z^-1
       - [1-z]^-1
       - 1/12*pi^2*z^-1
       - 1/12*pi^2*[1-z]^-1
       + 1/6*pi^2
       + 1/6*pi^2*x*z^-1
       + 1/6*pi^2*x*[1-z]^-1
       - 1/3*pi^2*x
       + 3*ln(x)
       - 3/2*ln(x)*z^-1
       - 3/2*ln(x)*[1-z]^-1
       + 6*ln(x)*x*z^-1
       + 6*ln(x)*x*[1-z]^-1
       - 18*ln(x)*x
       - ln(x)*ln( - [x-z])
       + 1/2*ln(x)*ln( - [x-z])*z^-1
       + 1/2*ln(x)*ln( - [x-z])*[1-z]^-1
       - ln(x)*ln( - [x-z])*x*z^-1
       - ln(x)*ln( - [x-z])*x*[1-z]^-1
       + 2*ln(x)*ln( - [x-z])*x
       - ln(x)*ln( - [1-x-z])
       + 1/2*ln(x)*ln( - [1-x-z])*z^-1
       + 1/2*ln(x)*ln( - [1-x-z])*[1-z]^-1
       - ln(x)*ln( - [1-x-z])*x*z^-1
       - ln(x)*ln( - [1-x-z])*x*[1-z]^-1
       + 2*ln(x)*ln( - [1-x-z])*x
       - 3*ln(x)^2
       + 3/2*ln(x)^2*z^-1
       + 3/2*ln(x)^2*[1-z]^-1
       - 3*ln(x)^2*x*z^-1
       - 3*ln(x)^2*x*[1-z]^-1
       + 6*ln(x)^2*x
       + 4*ln(x)*ln([1-x])
       - 2*ln(x)*ln([1-x])*z^-1
       - 2*ln(x)*ln([1-x])*[1-z]^-1
       + 4*ln(x)*ln([1-x])*x*z^-1
       + 4*ln(x)*ln([1-x])*x*[1-z]^-1
       - 8*ln(x)*ln([1-x])*x
       + 5*ln(x)*ln(z)
       - 5/2*ln(x)*ln(z)*z^-1
       - 5/2*ln(x)*ln(z)*[1-z]^-1
       + 5*ln(x)*ln(z)*x*z^-1
       + 5*ln(x)*ln(z)*x*[1-z]^-1
       - 10*ln(x)*ln(z)*x
       + 5*ln(x)*ln([1-z])
       - 5/2*ln(x)*ln([1-z])*z^-1
       - 5/2*ln(x)*ln([1-z])*[1-z]^-1
       + 5*ln(x)*ln([1-z])*x*z^-1
       + 5*ln(x)*ln([1-z])*x*[1-z]^-1
       - 10*ln(x)*ln([1-z])*x
       - ln([1-x])
       + 1/2*ln([1-x])*z^-1
       + 1/2*ln([1-x])*[1-z]^-1
       - 2*ln([1-x])*x*z^-1
       - 2*ln([1-x])*x*[1-z]^-1
       + 6*ln([1-x])*x
       - 1/2*ln([1-x])^2
       + 1/4*ln([1-x])^2*z^-1
       + 1/4*ln([1-x])^2*[1-z]^-1
       - 1/2*ln([1-x])^2*x*z^-1
       - 1/2*ln([1-x])^2*x*[1-z]^-1
       + ln([1-x])^2*x
       - 3*ln([1-x])*ln(z)
       + 3/2*ln([1-x])*ln(z)*z^-1
       + 3/2*ln([1-x])*ln(z)*[1-z]^-1
       - 3*ln([1-x])*ln(z)*x*z^-1
       - 3*ln([1-x])*ln(z)*x*[1-z]^-1
       + 6*ln([1-x])*ln(z)*x
       - 2*ln([1-x])*ln([1-z])
       + ln([1-x])*ln([1-z])*z^-1
       + ln([1-x])*ln([1-z])*[1-z]^-1
       - 2*ln([1-x])*ln([1-z])*x*z^-1
       - 2*ln([1-x])*ln([1-z])*x*[1-z]^-1
       + 4*ln([1-x])*ln([1-z])*x
       - 2*ln(z)
       + ln(z)*z^-1
       + ln(z)*[1-z]^-1
       - 4*ln(z)*x*z^-1
       - 4*ln(z)*x*[1-z]^-1
       + 12*ln(z)*x
       + ln(z)*ln( - [x-z])
       - 1/2*ln(z)*ln( - [x-z])*z^-1
       - 1/2*ln(z)*ln( - [x-z])*[1-z]^-1
       + ln(z)*ln( - [x-z])*x*z^-1
       + ln(z)*ln( - [x-z])*x*[1-z]^-1
       - 2*ln(z)*ln( - [x-z])*x
       - 3/2*ln(z)^2
       + 3/4*ln(z)^2*z^-1
       + 3/4*ln(z)^2*[1-z]^-1
       - 3/2*ln(z)^2*x*z^-1
       - 3/2*ln(z)^2*x*[1-z]^-1
       + 3*ln(z)^2*x
       - 4*ln(z)*ln([1-z])
       + 2*ln(z)*ln([1-z])*z^-1
       + 2*ln(z)*ln([1-z])*[1-z]^-1
       - 4*ln(z)*ln([1-z])*x*z^-1
       - 4*ln(z)*ln([1-z])*x*[1-z]^-1
       + 8*ln(z)*ln([1-z])*x
       - 2*ln([1-z])
       + ln([1-z])*z^-1
       + ln([1-z])*[1-z]^-1
       - 4*ln([1-z])*x*z^-1
       - 4*ln([1-z])*x*[1-z]^-1
       + 12*ln([1-z])*x
       + ln([1-z])*ln( - [1-x-z])
       - 1/2*ln([1-z])*ln( - [1-x-z])*z^-1
       - 1/2*ln([1-z])*ln( - [1-x-z])*[1-z]^-1
       + ln([1-z])*ln( - [1-x-z])*x*z^-1
       + ln([1-z])*ln( - [1-x-z])*x*[1-z]^-1
       - 2*ln([1-z])*ln( - [1-x-z])*x
       - 2*ln([1-z])^2
       + ln([1-z])^2*z^-1
       + ln([1-z])^2*[1-z]^-1
       - 2*ln([1-z])^2*x*z^-1
       - 2*ln([1-z])^2*x*[1-z]^-1
       + 4*ln([1-z])^2*x
       + Li2(x^-1*[1-x]^-1*z*[1-z])
       - 1/2*Li2(x^-1*[1-x]^-1*z*[1-z])*z^-1
       - 1/2*Li2(x^-1*[1-x]^-1*z*[1-z])*[1-z]^-1
       + Li2(x^-1*[1-x]^-1*z*[1-z])*x*z^-1
       + Li2(x^-1*[1-x]^-1*z*[1-z])*x*[1-z]^-1
       - 2*Li2(x^-1*[1-x]^-1*z*[1-z])*x
       - Li2([1-x]^-1*[1-z])
       + 1/2*Li2([1-x]^-1*[1-z])*z^-1
       + 1/2*Li2([1-x]^-1*[1-z])*[1-z]^-1
       - Li2([1-x]^-1*[1-z])*x*z^-1
       - Li2([1-x]^-1*[1-z])*x*[1-z]^-1
       + 2*Li2([1-x]^-1*[1-z])*x
       + Li2([1-x]*z^-1)
       - 1/2*Li2([1-x]*z^-1)*z^-1
       - 1/2*Li2([1-x]*z^-1)*[1-z]^-1
       + Li2([1-x]*z^-1)*x*z^-1
       + Li2([1-x]*z^-1)*x*[1-z]^-1
       - 2*Li2([1-x]*z^-1)*x ) );

id DC2Q2QNS    = (  + NC^-2 * ( 9
       + 6*[1-x]^-1*z^-1*[1-z]^-1
       - 6*[1-x]^-1*z^-1
       + 2*[1-x]^-1*z^-1*rln2^2
       + 3*[1-x]^-1*z^-1*sqrtxz1*rln2
       - 2*[1-x]^-1*[1-z]^-1
       + 2*[1-x]^-1
       + 4*[1-x]^-1*rln2^2
       + 3*[1-x]^-1*sqrtxz1*rln2
       - 6*[1-x]^-1*z
       + 2*[1-x]^-1*z*rln2^2
       - 4*z^-1*[1-z]^-1
       + z^-1
       - z^-1*rln2^2
       - z^-1*sqrtxz1*rln2
       - 5*[1-z]^-1
       + 1/2*poly2^-1
       - 2*rln2^2
       - 3*sqrtxz1*rln2
       + 1/4*z*[1-x-z]^-1
       + 7/2*z
       - 2*z*rln2^2
       - 2*x*z^-1*[1-z]^-1
       + 19/4*x*z^-1
       - x*z^-1*rln2^2
       - 2*x*z^-1*sqrtxz1*rln2
       + 7*x*[1-z]^-1
       - 1/4*x*[x-z]^-1
       - 1/2*x*poly2^-1
       - 15*x
       - 2*x*rln2^2
       - x*z
       + 1/2*x^2*[x-z]^-1
       + 1/2*x^2*poly2^-1
       - 1/2*x^3*poly2^-1
       + 1/6*pi^2*x^-2*[1+x]^-1*z^-1
       - 1/3*pi^2*x^-2*[1+x]^-1
       - 1/6*pi^2*x^-2*z^-1
       + 1/3*pi^2*x^-2
       - 1/3*pi^2*x^-1*[1+x]^-1*z
       + 1/6*pi^2*x^-1*z^-1
       - 1/3*pi^2*x^-1
       + 1/3*pi^2*x^-1*z
       - 2/3*pi^2*[1-x]^-1*z^-1*[1-z]^-1
       + 2/3*pi^2*[1-x]^-1*z^-1
       - 1/6*pi^2*[1-x]^-1*[1-z]^-1
       + 2/3*pi^2*[1-x]^-1
       + 2/3*pi^2*[1-x]^-1*z
       - 1/3*pi^2*[1+x]^-1*z^-1
       + 2/3*pi^2*[1+x]^-1
       - 1/2*pi^2*[1+x]^-1*z
       + 1/3*pi^2*z^-1*[1-z]^-1
       - 1/6*pi^2*z^-1
       + 1/3*pi^2*[1-z]^-1
       - 11/12*pi^2
       - 13/12*pi^2*z
       + 1/3*pi^2*x*z^-1*[1-z]^-1
       - 1/3*pi^2*x*z^-1
       - 1/6*pi^2*x*[1-z]^-1
       - 3/4*pi^2*x
       - 1/12*pi^2*x*z
       - 2*ln(1 + sqrtxz1 - z)*[1-x]^-1*z^-1*rln2
       - 3*ln(1 + sqrtxz1 - z)*[1-x]^-1*z^-1*sqrtxz1
       - 4*ln(1 + sqrtxz1 - z)*[1-x]^-1*rln2
       - 3*ln(1 + sqrtxz1 - z)*[1-x]^-1*sqrtxz1
       - 2*ln(1 + sqrtxz1 - z)*[1-x]^-1*z*rln2
       + ln(1 + sqrtxz1 - z)*z^-1*rln2
       + ln(1 + sqrtxz1 - z)*z^-1*sqrtxz1
       + 2*ln(1 + sqrtxz1 - z)*rln2
       + 3*ln(1 + sqrtxz1 - z)*sqrtxz1
       + 2*ln(1 + sqrtxz1 - z)*z*rln2
       + ln(1 + sqrtxz1 - z)*x*z^-1*rln2
       + 2*ln(1 + sqrtxz1 - z)*x*z^-1*sqrtxz1
       + 2*ln(1 + sqrtxz1 - z)*x*rln2
       - 2*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)
       + 2*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^-1
       + 4*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*[1-x]^-1
       + 2*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z
       - ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*z^-1
       - 2*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*z
       - ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*x*z^-1
       - 2*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*x
       - 2*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^-1*rln2
       - 4*ln(1 + sqrtxz1 + z)*[1-x]^-1*rln2
       - 2*ln(1 + sqrtxz1 + z)*[1-x]^-1*z*rln2
       + ln(1 + sqrtxz1 + z)*z^-1*rln2
       + 2*ln(1 + sqrtxz1 + z)*rln2
       + 2*ln(1 + sqrtxz1 + z)*z*rln2
       + ln(1 + sqrtxz1 + z)*x*z^-1*rln2
       + 2*ln(1 + sqrtxz1 + z)*x*rln2
       + 3*ln(z*sqrtxz3)*ArcTan(z*sqrtxz3)*x^-1*sqrtxz3
       - ln(z*sqrtxz3)*ArcTan(z*sqrtxz3)*z^-1*sqrtxz3
       + 3*ln(z*sqrtxz3)*ArcTan(z*sqrtxz3)*z*sqrtxz3
       + 3*ln(z*sqrtxz3)*ArcTan(z*sqrtxz3)*x*sqrtxz3
       - 15/2*ln(x)
       + 2*ln(x)*x^-2*[1+x]^-1*z^-1
       - 4*ln(x)*x^-2*[1+x]^-1
       - 2*ln(x)*x^-2*z^-1
       + 4*ln(x)*x^-2
       - 4*ln(x)*x^-1*[1+x]^-1*z
       + 2*ln(x)*x^-1*z^-1
       - 4*ln(x)*x^-1
       + 4*ln(x)*x^-1*z
       - 24*ln(x)*[1-x]^-1*z^-1*[1-z]^-1
       + 87/4*ln(x)*[1-x]^-1*z^-1
       + ln(x)*[1-x]^-1*z^-1*rln2
       + 3/2*ln(x)*[1-x]^-1*z^-1*sqrtxz1
       + 6*ln(x)*[1-x]^-1*[1-z]^-1
       - 3/4*ln(x)*[1-x]^-1*[x-z]^-1
       + 10*ln(x)*[1-x]^-1
       + 2*ln(x)*[1-x]^-1*rln2
       + 3/2*ln(x)*[1-x]^-1*sqrtxz1
       + 73/4*ln(x)*[1-x]^-1*z
       + ln(x)*[1-x]^-1*z*rln2
       - 2*ln(x)*[1+x]^-1*z^-1
       + 4*ln(x)*[1+x]^-1
       - 4*ln(x)*[1+x]^-1*z
       + 16*ln(x)*z^-1*[1-z]^-1
       - 33/2*ln(x)*z^-1
       - 1/2*ln(x)*z^-1*rln2
       - 1/2*ln(x)*z^-1*sqrtxz1
       + 2*ln(x)*[1-z]^-1
       - 1/4*ln(x)*[1-x-z]^-1
       + 3/4*ln(x)*[x-z]^-1
       + 3/4*ln(x)*poly2^-2
       - 1/2*ln(x)*poly2^-1
       - ln(x)*rln2
       - 3/2*ln(x)*sqrtxz1
       + 1/4*ln(x)*z*[1-x-z]^-2
       - 3/4*ln(x)*z*[1-x-z]^-1
       - 79/4*ln(x)*z
       - ln(x)*z*rln2
       - 1/4*ln(x)*z^2*[1-x-z]^-2
       + ln(x)*[1+x]
       - ln(x)*[1+x]*z
       + 8*ln(x)*x*z^-1*[1-z]^-1
       - 8*ln(x)*x*z^-1
       - 1/2*ln(x)*x*z^-1*rln2
       - ln(x)*x*z^-1*sqrtxz1
       - 8*ln(x)*x*[1-z]^-1
       + 3/4*ln(x)*x*[x-z]^-1
       - 3/4*ln(x)*x*poly2^-2
       - 11/2*ln(x)*x
       - ln(x)*x*rln2
       + 3/4*ln(x)*x*z
       + 1/4*ln(x)*x^2*[x-z]^-2
       + 3/2*ln(x)*x^2*[x-z]^-1
       - 1/2*ln(x)*x^3*[x-z]^-2
       - 1/2*ln(x)*x^3*poly2^-1
       - 3/4*ln(x)*x^4*poly2^-2
       + 3/4*ln(x)*x^5*poly2^-2
       - 1/2*ln(x)*ln(1 - sqrtxz2 + x)*z^-1*sqrtxz2^-1
       - 1/2*ln(x)*ln(1 - sqrtxz2 + x)*[1-z]^-1*sqrtxz2^-1
       + 3/8*ln(x)*ln(1 - sqrtxz2 + x)*sqrtxz2^-1*poly2^-2
       - 3/8*ln(x)*ln(1 - sqrtxz2 + x)*sqrtxz2^-1*poly2^-1
       + 2*ln(x)*ln(1 - sqrtxz2 + x)*sqrtxz2^-1
       + ln(x)*ln(1 - sqrtxz2 + x)*x*z^-1*sqrtxz2^-1
       - ln(x)*ln(1 - sqrtxz2 + x)*x*[1-z]^-1*sqrtxz2^-1
       + 3/4*ln(x)*ln(1 - sqrtxz2 + x)*x*sqrtxz2^-1
       - 3/2*ln(x)*ln(1 - sqrtxz2 + x)*x*z*sqrtxz2^-1
       - 1/2*ln(x)*ln(1 - sqrtxz2 + x)*x^2*z^-1*sqrtxz2^-1
       - 1/2*ln(x)*ln(1 - sqrtxz2 + x)*x^2*[1-z]^-1*sqrtxz2^-1
       - 3/8*ln(x)*ln(1 - sqrtxz2 + x)*x^2*sqrtxz2^-1*poly2^-2
       - 1/4*ln(x)*ln(1 - sqrtxz2 + x)*x^2*sqrtxz2^-1*poly2^-1
       + 2*ln(x)*ln(1 - sqrtxz2 + x)*x^2*sqrtxz2^-1
       - 3/8*ln(x)*ln(1 - sqrtxz2 + x)*x^4*sqrtxz2^-1*poly2^-2
       - 3/8*ln(x)*ln(1 - sqrtxz2 + x)*x^4*sqrtxz2^-1*poly2^-1
       + 3/8*ln(x)*ln(1 - sqrtxz2 + x)*x^6*sqrtxz2^-1*poly2^-2
       + 1/2*ln(x)*ln(1 + sqrtxz2 + x)*z^-1*sqrtxz2^-1
       + 1/2*ln(x)*ln(1 + sqrtxz2 + x)*[1-z]^-1*sqrtxz2^-1
       - 3/8*ln(x)*ln(1 + sqrtxz2 + x)*sqrtxz2^-1*poly2^-2
       + 3/8*ln(x)*ln(1 + sqrtxz2 + x)*sqrtxz2^-1*poly2^-1
       - 2*ln(x)*ln(1 + sqrtxz2 + x)*sqrtxz2^-1
       - ln(x)*ln(1 + sqrtxz2 + x)*x*z^-1*sqrtxz2^-1
       + ln(x)*ln(1 + sqrtxz2 + x)*x*[1-z]^-1*sqrtxz2^-1
       - 3/4*ln(x)*ln(1 + sqrtxz2 + x)*x*sqrtxz2^-1
       + 3/2*ln(x)*ln(1 + sqrtxz2 + x)*x*z*sqrtxz2^-1
       + 1/2*ln(x)*ln(1 + sqrtxz2 + x)*x^2*z^-1*sqrtxz2^-1
       + 1/2*ln(x)*ln(1 + sqrtxz2 + x)*x^2*[1-z]^-1*sqrtxz2^-1
       + 3/8*ln(x)*ln(1 + sqrtxz2 + x)*x^2*sqrtxz2^-1*poly2^-2
       + 1/4*ln(x)*ln(1 + sqrtxz2 + x)*x^2*sqrtxz2^-1*poly2^-1
       - 2*ln(x)*ln(1 + sqrtxz2 + x)*x^2*sqrtxz2^-1
       + 3/8*ln(x)*ln(1 + sqrtxz2 + x)*x^4*sqrtxz2^-1*poly2^-2
       + 3/8*ln(x)*ln(1 + sqrtxz2 + x)*x^4*sqrtxz2^-1*poly2^-1
       - 3/8*ln(x)*ln(1 + sqrtxz2 + x)*x^6*sqrtxz2^-1*poly2^-2
       + ln(x)*ln(1 + sqrtxz1 + z)
       - ln(x)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^-1
       - 2*ln(x)*ln(1 + sqrtxz1 + z)*[1-x]^-1
       - ln(x)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z
       + 1/2*ln(x)*ln(1 + sqrtxz1 + z)*z^-1
       + ln(x)*ln(1 + sqrtxz1 + z)*z
       + 1/2*ln(x)*ln(1 + sqrtxz1 + z)*x*z^-1
       + ln(x)*ln(1 + sqrtxz1 + z)*x
       - ln(x)*ln(1 + x*z^-1)
       - 1/2*ln(x)*ln(1 + x*z^-1)*x^-2*[1+x]^-1*z^-1
       + ln(x)*ln(1 + x*z^-1)*x^-2*[1+x]^-1
       + 1/2*ln(x)*ln(1 + x*z^-1)*x^-2*z^-1
       - ln(x)*ln(1 + x*z^-1)*x^-2
       + ln(x)*ln(1 + x*z^-1)*x^-1*[1+x]^-1*z
       - 1/2*ln(x)*ln(1 + x*z^-1)*x^-1*z^-1
       + ln(x)*ln(1 + x*z^-1)*x^-1
       - ln(x)*ln(1 + x*z^-1)*x^-1*z
       - 1/2*ln(x)*ln(1 + x*z^-1)*[1+x]^-1*z^-1
       + ln(x)*ln(1 + x*z^-1)*[1+x]^-1
       + 1/2*ln(x)*ln(1 + x*z^-1)*z^-1
       + ln(x)*ln(1 + x*z^-1)*z
       - 1/2*ln(x)*ln(1 + x*z^-1)*x*z^-1
       + ln(x)*ln(1 + x*z^-1)*x
       - 4*ln(x)*ln(1 + x)
       - 2*ln(x)*ln(1 + x)*[1+x]^-1*z^-1
       + 4*ln(x)*ln(1 + x)*[1+x]^-1
       - 2*ln(x)*ln(1 + x)*[1+x]^-1*z
       + 2*ln(x)*ln(1 + x)*z^-1
       + 2*ln(x)*ln(1 + x)*z
       - 2*ln(x)*ln(1 + x*z)
       - 1/2*ln(x)*ln(1 + x*z)*x^-2*[1+x]^-1*z^-1
       + ln(x)*ln(1 + x*z)*x^-2*[1+x]^-1
       + 1/2*ln(x)*ln(1 + x*z)*x^-2*z^-1
       - ln(x)*ln(1 + x*z)*x^-2
       + ln(x)*ln(1 + x*z)*x^-1*[1+x]^-1*z
       - 1/2*ln(x)*ln(1 + x*z)*x^-1*z^-1
       + ln(x)*ln(1 + x*z)*x^-1
       - ln(x)*ln(1 + x*z)*x^-1*z
       + ln(x)*ln(1 + x*z)*[1-x]^-1*z^-1
       + 2*ln(x)*ln(1 + x*z)*[1-x]^-1
       + ln(x)*ln(1 + x*z)*[1-x]^-1*z
       - 1/2*ln(x)*ln(1 + x*z)*[1+x]^-1*z^-1
       + ln(x)*ln(1 + x*z)*[1+x]^-1
       - ln(x)*ln(1 + x*z)*x*z^-1
       + ln(x)*ln(z + x)
       - ln(x)*ln(z + x)*[1-x]^-1*z^-1
       - 2*ln(x)*ln(z + x)*[1-x]^-1
       - ln(x)*ln(z + x)*[1-x]^-1*z
       + 1/2*ln(x)*ln(z + x)*z^-1
       + ln(x)*ln(z + x)*z
       + 1/2*ln(x)*ln(z + x)*x*z^-1
       + ln(x)*ln(z + x)*x
       + 21/4*ln(x)^2
       - 3/4*ln(x)^2*x^-2*[1+x]^-1*z^-1
       + 3/2*ln(x)^2*x^-2*[1+x]^-1
       + 3/4*ln(x)^2*x^-2*z^-1
       - 3/2*ln(x)^2*x^-2
       + 3/2*ln(x)^2*x^-1*[1+x]^-1*z
       - 3/4*ln(x)^2*x^-1*z^-1
       + 3/2*ln(x)^2*x^-1
       - 3/2*ln(x)^2*x^-1*z
       + 8*ln(x)^2*[1-x]^-1*z^-1*[1-z]^-1
       - 8*ln(x)^2*[1-x]^-1*z^-1
       + 1/2*ln(x)^2*[1-x]^-1*[1-z]^-1
       - 13/2*ln(x)^2*[1-x]^-1
       - 13/2*ln(x)^2*[1-x]^-1*z
       + 5/4*ln(x)^2*[1+x]^-1*z^-1
       - 5/2*ln(x)^2*[1+x]^-1
       + 2*ln(x)^2*[1+x]^-1*z
       - 11/2*ln(x)^2*z^-1*[1-z]^-1
       + 5*ln(x)^2*z^-1
       - 3/4*ln(x)^2*[1-z]^-1
       + 31/4*ln(x)^2*z
       - 5/2*ln(x)^2*x*z^-1*[1-z]^-1
       + 5/2*ln(x)^2*x*z^-1
       + 1/4*ln(x)^2*x*[1-z]^-1
       + 17/4*ln(x)^2*x
       + 1/4*ln(x)^2*x*z
       - 15/2*ln(x)*ln([1-x])
       - ln(x)*ln([1-x])*x^-2*[1+x]^-1*z^-1
       + 2*ln(x)*ln([1-x])*x^-2*[1+x]^-1
       + ln(x)*ln([1-x])*x^-2*z^-1
       - 2*ln(x)*ln([1-x])*x^-2
       + 2*ln(x)*ln([1-x])*x^-1*[1+x]^-1*z
       - ln(x)*ln([1-x])*x^-1*z^-1
       + 2*ln(x)*ln([1-x])*x^-1
       - 2*ln(x)*ln([1-x])*x^-1*z
       + 4*ln(x)*ln([1-x])*[1-x]^-1*z^-1*[1-z]^-1
       - 4*ln(x)*ln([1-x])*[1-x]^-1*z^-1
       - 14*ln(x)*ln([1-x])*[1-x]^-1*[1-z]^-1
       + 17/2*ln(x)*ln([1-x])*[1-x]^-1
       + 17/2*ln(x)*ln([1-x])*[1-x]^-1*z
       + ln(x)*ln([1-x])*[1+x]^-1*z^-1
       - 2*ln(x)*ln([1-x])*[1+x]^-1
       + 2*ln(x)*ln([1-x])*[1+x]^-1*z
       - 2*ln(x)*ln([1-x])*z^-1*[1-z]^-1
       + 2*ln(x)*ln([1-x])*z^-1
       + 21/2*ln(x)*ln([1-x])*[1-z]^-1
       - 23/2*ln(x)*ln([1-x])*z
       - 2*ln(x)*ln([1-x])*x*z^-1*[1-z]^-1
       + 2*ln(x)*ln([1-x])*x*z^-1
       + 7/2*ln(x)*ln([1-x])*x*[1-z]^-1
       - 9/2*ln(x)*ln([1-x])*x
       - 1/2*ln(x)*ln([1-x])*x*z
       + 2*ln(x)*ln([1+x])
       + ln(x)*ln([1+x])*x^-2*[1+x]^-1*z^-1
       - 2*ln(x)*ln([1+x])*x^-2*[1+x]^-1
       - ln(x)*ln([1+x])*x^-2*z^-1
       + 2*ln(x)*ln([1+x])*x^-2
       - 2*ln(x)*ln([1+x])*x^-1*[1+x]^-1*z
       + ln(x)*ln([1+x])*x^-1*z^-1
       - 2*ln(x)*ln([1+x])*x^-1
       + 2*ln(x)*ln([1+x])*x^-1*z
       + ln(x)*ln([1+x])*[1+x]^-1*z^-1
       - 2*ln(x)*ln([1+x])*[1+x]^-1
       - ln(x)*ln([1+x])*z^-1
       - 2*ln(x)*ln([1+x])*z
       + ln(x)*ln([1+x])*x*z^-1
       - 2*ln(x)*ln([1+x])*x
       - 2*ln(x)*ln(z)
       + 1/2*ln(x)*ln(z)*x^-2*[1+x]^-1*z^-1
       - ln(x)*ln(z)*x^-2*[1+x]^-1
       - 1/2*ln(x)*ln(z)*x^-2*z^-1
       + ln(x)*ln(z)*x^-2
       - ln(x)*ln(z)*x^-1*[1+x]^-1*z
       + 1/2*ln(x)*ln(z)*x^-1*z^-1
       - ln(x)*ln(z)*x^-1
       + ln(x)*ln(z)*x^-1*z
       - 8*ln(x)*ln(z)*[1-x]^-1*z^-1*[1-z]^-1
       + 8*ln(x)*ln(z)*[1-x]^-1*z^-1
       - 2*ln(x)*ln(z)*[1-x]^-1*[1-z]^-1
       + 23/4*ln(x)*ln(z)*[1-x]^-1
       + 27/4*ln(x)*ln(z)*[1-x]^-1*z
       + 1/2*ln(x)*ln(z)*[1+x]^-1*z^-1
       - ln(x)*ln(z)*[1+x]^-1
       + 5*ln(x)*ln(z)*z^-1*[1-z]^-1
       - 6*ln(x)*ln(z)*z^-1
       + 2*ln(x)*ln(z)*[1-z]^-1
       - 19/2*ln(x)*ln(z)*z
       + 3*ln(x)*ln(z)*x*z^-1*[1-z]^-1
       - 3*ln(x)*ln(z)*x*z^-1
       - 9/2*ln(x)*ln(z)*x
       - 7/2*ln(x)*ln([1-z])
       - 9*ln(x)*ln([1-z])*[1-x]^-1*[1-z]^-1
       + 8*ln(x)*ln([1-z])*[1-x]^-1
       + 8*ln(x)*ln([1-z])*[1-x]^-1*z
       + 6*ln(x)*ln([1-z])*[1-z]^-1
       - 10*ln(x)*ln([1-z])*z
       + 3*ln(x)*ln([1-z])*x*[1-z]^-1
       - 7*ln(x)*ln([1-z])*x
       - 1/2*ln(x)*ln([1-z])*x*z
       + 19/4*ln([1-x])
       + 9*ln([1-x])*[1-x]^-1*[1-z]^-1
       - 4*ln([1-x])*[1-x]^-1
       - 7*ln([1-x])*[1-x]^-1*z
       - 9*ln([1-x])*[1-z]^-1
       + 35/4*ln([1-x])*z
       + 9/4*ln([1-x])*x
       + 1/4*ln([1-x])*x*z
       + 9/4*ln([1-x])^2
       + 5/2*ln([1-x])^2*[1-x]^-1*[1-z]^-1
       - 3/2*ln([1-x])^2*[1-x]^-1
       - 3/2*ln([1-x])^2*[1-x]^-1*z
       - 9/4*ln([1-x])^2*[1-z]^-1
       + 13/4*ln([1-x])^2*z
       - 1/4*ln([1-x])^2*x*[1-z]^-1
       + 5/4*ln([1-x])^2*x
       + 1/4*ln([1-x])^2*x*z
       + 3*ln([1-x])*ln(z)
       + 6*ln([1-x])*ln(z)*[1-x]^-1*[1-z]^-1
       - 4*ln([1-x])*ln(z)*[1-x]^-1
       - 4*ln([1-x])*ln(z)*[1-x]^-1*z
       - 9/2*ln([1-x])*ln(z)*[1-z]^-1
       + 11/2*ln([1-x])*ln(z)*z
       - 3/2*ln([1-x])*ln(z)*x*[1-z]^-1
       + 5/2*ln([1-x])*ln(z)*x
       + 5/2*ln([1-x])*ln([1-z])
       + 4*ln([1-x])*ln([1-z])*[1-x]^-1*[1-z]^-1
       - 3*ln([1-x])*ln([1-z])*[1-x]^-1
       - 3*ln([1-x])*ln([1-z])*[1-x]^-1*z
       - 3*ln([1-x])*ln([1-z])*[1-z]^-1
       + 11/2*ln([1-x])*ln([1-z])*z
       - ln([1-x])*ln([1-z])*x*[1-z]^-1
       + 7/2*ln([1-x])*ln([1-z])*x
       + 1/2*ln([1-x])*ln([1-z])*x*z
       + 41/4*ln(z)
       + 12*ln(z)*[1-x]^-1*z^-1*[1-z]^-1
       - 27/2*ln(z)*[1-x]^-1*z^-1
       + 3*ln(z)*[1-x]^-1*z^-1*rln2
       + 3/2*ln(z)*[1-x]^-1*z^-1*sqrtxz1
       - 3/4*ln(z)*[1-x]^-1*[1-z]^-1
       + 3/4*ln(z)*[1-x]^-1*[x-z]^-1
       - 9*ln(z)*[1-x]^-1
       + 6*ln(z)*[1-x]^-1*rln2
       + 3/2*ln(z)*[1-x]^-1*sqrtxz1
       - 23/2*ln(z)*[1-x]^-1*z
       + 3*ln(z)*[1-x]^-1*z*rln2
       - 8*ln(z)*z^-1*[1-z]^-1
       + 8*ln(z)*z^-1
       - 3/2*ln(z)*z^-1*rln2
       - 1/2*ln(z)*z^-1*sqrtxz1
       - 4*ln(z)*[1-z]^-1
       - 3/4*ln(z)*[x-z]^-1
       + 3/4*ln(z)*poly2^-2
       - 1/2*ln(z)*poly2^-1
       - 3*ln(z)*rln2
       - 3/2*ln(z)*sqrtxz1
       + 15*ln(z)*z
       - 3*ln(z)*z*rln2
       - 4*ln(z)*x*z^-1*[1-z]^-1
       + 11/2*ln(z)*x*z^-1
       - 3/2*ln(z)*x*z^-1*rln2
       - ln(z)*x*z^-1*sqrtxz1
       + 4*ln(z)*x*[1-z]^-1
       - 3/4*ln(z)*x*[x-z]^-1
       + 3/4*ln(z)*x*poly2^-2
       + 15/4*ln(z)*x
       - 3*ln(z)*x*rln2
       - 1/4*ln(z)*x^2*[x-z]^-2
       - 3/2*ln(z)*x^2*[x-z]^-1
       + 1/2*ln(z)*x^3*[x-z]^-2
       + 1/2*ln(z)*x^3*poly2^-1
       - 3/4*ln(z)*x^4*poly2^-2
       - 3/4*ln(z)*x^5*poly2^-2
       + ln(z)*ln(1 + sqrtxz1 - z)
       - ln(z)*ln(1 + sqrtxz1 - z)*[1-x]^-1*z^-1
       - 2*ln(z)*ln(1 + sqrtxz1 - z)*[1-x]^-1
       - ln(z)*ln(1 + sqrtxz1 - z)*[1-x]^-1*z
       + 1/2*ln(z)*ln(1 + sqrtxz1 - z)*z^-1
       + ln(z)*ln(1 + sqrtxz1 - z)*z
       + 1/2*ln(z)*ln(1 + sqrtxz1 - z)*x*z^-1
       + ln(z)*ln(1 + sqrtxz1 - z)*x
       + 2*ln(z)*ln(1 + sqrtxz1 + z)
       - 2*ln(z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^-1
       - 4*ln(z)*ln(1 + sqrtxz1 + z)*[1-x]^-1
       - 2*ln(z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z
       + ln(z)*ln(1 + sqrtxz1 + z)*z^-1
       + 2*ln(z)*ln(1 + sqrtxz1 + z)*z
       + ln(z)*ln(1 + sqrtxz1 + z)*x*z^-1
       + 2*ln(z)*ln(1 + sqrtxz1 + z)*x
       + ln(z)*ln(1 + x*z^-1)
       + 1/2*ln(z)*ln(1 + x*z^-1)*x^-2*[1+x]^-1*z^-1
       - ln(z)*ln(1 + x*z^-1)*x^-2*[1+x]^-1
       - 1/2*ln(z)*ln(1 + x*z^-1)*x^-2*z^-1
       + ln(z)*ln(1 + x*z^-1)*x^-2
       - ln(z)*ln(1 + x*z^-1)*x^-1*[1+x]^-1*z
       + 1/2*ln(z)*ln(1 + x*z^-1)*x^-1*z^-1
       - ln(z)*ln(1 + x*z^-1)*x^-1
       + ln(z)*ln(1 + x*z^-1)*x^-1*z
       + 1/2*ln(z)*ln(1 + x*z^-1)*[1+x]^-1*z^-1
       - ln(z)*ln(1 + x*z^-1)*[1+x]^-1
       - 1/2*ln(z)*ln(1 + x*z^-1)*z^-1
       - ln(z)*ln(1 + x*z^-1)*z
       + 1/2*ln(z)*ln(1 + x*z^-1)*x*z^-1
       - ln(z)*ln(1 + x*z^-1)*x
       - 2*ln(z)*ln(1 + x*z)
       - 1/2*ln(z)*ln(1 + x*z)*x^-2*[1+x]^-1*z^-1
       + ln(z)*ln(1 + x*z)*x^-2*[1+x]^-1
       + 1/2*ln(z)*ln(1 + x*z)*x^-2*z^-1
       - ln(z)*ln(1 + x*z)*x^-2
       + ln(z)*ln(1 + x*z)*x^-1*[1+x]^-1*z
       - 1/2*ln(z)*ln(1 + x*z)*x^-1*z^-1
       + ln(z)*ln(1 + x*z)*x^-1
       - ln(z)*ln(1 + x*z)*x^-1*z
       + ln(z)*ln(1 + x*z)*[1-x]^-1*z^-1
       + 2*ln(z)*ln(1 + x*z)*[1-x]^-1
       + ln(z)*ln(1 + x*z)*[1-x]^-1*z
       - 1/2*ln(z)*ln(1 + x*z)*[1+x]^-1*z^-1
       + ln(z)*ln(1 + x*z)*[1+x]^-1
       - ln(z)*ln(1 + x*z)*x*z^-1
       - ln(z)*ln(z + x)
       + ln(z)*ln(z + x)*[1-x]^-1*z^-1
       + 2*ln(z)*ln(z + x)*[1-x]^-1
       + ln(z)*ln(z + x)*[1-x]^-1*z
       - 1/2*ln(z)*ln(z + x)*z^-1
       - ln(z)*ln(z + x)*z
       - 1/2*ln(z)*ln(z + x)*x*z^-1
       - ln(z)*ln(z + x)*x
       + 3/4*ln(z)^2
       + 1/4*ln(z)^2*x^-2*[1+x]^-1*z^-1
       - 1/2*ln(z)^2*x^-2*[1+x]^-1
       - 1/4*ln(z)^2*x^-2*z^-1
       + 1/2*ln(z)^2*x^-2
       - 1/2*ln(z)^2*x^-1*[1+x]^-1*z
       + 1/4*ln(z)^2*x^-1*z^-1
       - 1/2*ln(z)^2*x^-1
       + 1/2*ln(z)^2*x^-1*z
       + ln(z)^2*[1-x]^-1*z^-1*[1-z]^-1
       - ln(z)^2*[1-x]^-1*z^-1
       + 3/2*ln(z)^2*[1-x]^-1*[1-z]^-1
       - 3/2*ln(z)^2*[1-x]^-1
       - 3/2*ln(z)^2*[1-x]^-1*z
       + 1/4*ln(z)^2*[1+x]^-1*z^-1
       - 1/2*ln(z)^2*[1+x]^-1
       - 1/2*ln(z)^2*z^-1*[1-z]^-1
       + 1/4*ln(z)^2*z^-1
       - 1/4*ln(z)^2*[1-z]^-1
       + 1/4*ln(z)^2*z
       - 1/2*ln(z)^2*x*z^-1*[1-z]^-1
       + 3/4*ln(z)^2*x*z^-1
       + 1/4*ln(z)^2*x*[1-z]^-1
       - 1/4*ln(z)^2*x
       - 1/4*ln(z)^2*x*z
       + 3*ln(z)*ln([1-z])
       + 7*ln(z)*ln([1-z])*[1-x]^-1*[1-z]^-1
       - 13/2*ln(z)*ln([1-z])*[1-x]^-1
       - 13/2*ln(z)*ln([1-z])*[1-x]^-1*z
       - 11/2*ln(z)*ln([1-z])*[1-z]^-1
       + 9*ln(z)*ln([1-z])*z
       - 5/2*ln(z)*ln([1-z])*x*[1-z]^-1
       + 6*ln(z)*ln([1-z])*x
       + 6*ln([1-z])
       + 9*ln([1-z])*[1-x]^-1*[1-z]^-1
       - 5*ln([1-z])*[1-x]^-1
       - 8*ln([1-z])*[1-x]^-1*z
       - 9*ln([1-z])*[1-z]^-1
       + 1/4*ln([1-z])*[1-x-z]^-1
       - 1/4*ln([1-z])*z*[1-x-z]^-2
       + 3/4*ln([1-z])*z*[1-x-z]^-1
       + 8*ln([1-z])*z
       + 1/4*ln([1-z])*z^2*[1-x-z]^-2
       + 2*ln([1-z])*x
       + ln([1-z])*x*z
       + 3/4*ln([1-z])^2
       + 5/2*ln([1-z])^2*[1-x]^-1*[1-z]^-1
       - 9/4*ln([1-z])^2*[1-x]^-1
       - 9/4*ln([1-z])^2*[1-x]^-1*z
       - 3/2*ln([1-z])^2*[1-z]^-1
       + 13/4*ln([1-z])^2*z
       - ln([1-z])^2*x*[1-z]^-1
       + 11/4*ln([1-z])^2*x
       + 1/4*ln([1-z])^2*x*z
       - 3*ln(sqrtxz3)*ArcTan(sqrtxz3)*x^-1*sqrtxz3
       + ln(sqrtxz3)*ArcTan(sqrtxz3)*z^-1*sqrtxz3
       - 3*ln(sqrtxz3)*ArcTan(sqrtxz3)*z*sqrtxz3
       - 3*ln(sqrtxz3)*ArcTan(sqrtxz3)*x*sqrtxz3
       - 1/2*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*z^-1*sqrtxz2^-1
       - 1/2*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*[1-z]^-1*sqrtxz2^-1
       + 3/8*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*sqrtxz2^-1*poly2^-2
       - 3/8*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*sqrtxz2^-1*poly2^-1
       + 2*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*sqrtxz2^-1
       + Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x*z^-1*sqrtxz2^-1
       - Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x*[1-z]^-1*sqrtxz2^-1
       + 3/4*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x*sqrtxz2^-1
       - 3/2*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x*z*sqrtxz2^-1
       - 1/2*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x^2*z^-1*sqrtxz2^-1
       - 1/2*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x^2*[1-z]^-1*sqrtxz2^-1
       - 3/8*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x^2*sqrtxz2^-1*poly2^-2
       - 1/4*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x^2*sqrtxz2^-1*poly2^-1
       + 2*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x^2*sqrtxz2^-1
       - 3/8*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x^4*sqrtxz2^-1*poly2^-2
       - 3/8*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x^4*sqrtxz2^-1*poly2^-1
       + 3/8*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x^6*sqrtxz2^-1*poly2^-2
       + 1/2*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*z^-1*sqrtxz2^-1
       + 1/2*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*[1-z]^-1*sqrtxz2^-1
       - 3/8*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*sqrtxz2^-1*poly2^-2
       + 3/8*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*sqrtxz2^-1*poly2^-1
       - 2*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*sqrtxz2^-1
       - Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x*z^-1*sqrtxz2^-1
       + Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x*[1-z]^-1*sqrtxz2^-1
       - 3/4*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x*sqrtxz2^-1
       + 3/2*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x*z*sqrtxz2^-1
       + 1/2*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x^2*z^-1*sqrtxz2^-1
       + 1/2*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x^2*[1-z]^-1*sqrtxz2^-1
       + 3/8*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x^2*sqrtxz2^-1*poly2^-2
       + 1/4*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x^2*sqrtxz2^-1*poly2^-1
       - 2*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x^2*sqrtxz2^-1
       + 3/8*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x^4*sqrtxz2^-1*poly2^-2
       + 3/8*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x^4*sqrtxz2^-1*poly2^-1
       - 3/8*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x^6*sqrtxz2^-1*poly2^-2
       - Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)
       + Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*z^-1
       + 2*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1
       + Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*z
       - 1/2*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z^-1
       - Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z
       - 1/2*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x*z^-1
       - Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x
       + Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)
       - Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*z^-1
       - 2*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1
       - Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*z
       + 1/2*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z^-1
       + Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z
       + 1/2*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x*z^-1
       + Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x
       + 1/2*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*z^-1*sqrtxz2^-1
       + 1/2*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*[1-z]^-1*sqrtxz2^-1
       - 3/8*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*sqrtxz2^-1*poly2^-2
       + 3/8*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*sqrtxz2^-1*poly2^-1
       - 2*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*sqrtxz2^-1
       - Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x*z^-1*sqrtxz2^-1
       + Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x*[1-z]^-1*sqrtxz2^-1
       - 3/4*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x*sqrtxz2^-1
       + 3/2*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x*z*sqrtxz2^-1
       + 1/2*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x^2*z^-1*sqrtxz2^-1
       + 1/2*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x^2*[1-z]^-1*sqrtxz2^-1
       + 3/8*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x^2*sqrtxz2^-1*poly2^-2
       + 1/4*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x^2*sqrtxz2^-1*poly2^-1
       - 2*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x^2*sqrtxz2^-1
       + 3/8*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x^4*sqrtxz2^-1*poly2^-2
       + 3/8*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x^4*sqrtxz2^-1*poly2^-1
       - 3/8*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x^6*sqrtxz2^-1*poly2^-2
       - 1/2*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*z^-1*sqrtxz2^-1
       - 1/2*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*[1-z]^-1*sqrtxz2^-1
       + 3/8*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*sqrtxz2^-1*poly2^-2
       - 3/8*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*sqrtxz2^-1*poly2^-1
       + 2*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*sqrtxz2^-1
       + Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x*z^-1*sqrtxz2^-1
       - Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x*[1-z]^-1*sqrtxz2^-1
       + 3/4*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x*sqrtxz2^-1
       - 3/2*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x*z*sqrtxz2^-1
       - 1/2*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x^2*z^-1*sqrtxz2^-1
       - 1/2*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x^2*[1-z]^-1*sqrtxz2^-1
       - 3/8*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x^2*sqrtxz2^-1*poly2^-2
       - 1/4*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x^2*sqrtxz2^-1*poly2^-1
       + 2*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x^2*sqrtxz2^-1
       - 3/8*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x^4*sqrtxz2^-1*poly2^-2
       - 3/8*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x^4*sqrtxz2^-1*poly2^-1
       + 3/8*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x^6*sqrtxz2^-1*poly2^-2
       + Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)
       - Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*[1-x]^-1*z^-1
       - 2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*[1-x]^-1
       - Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*[1-x]^-1*z
       + 1/2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*z^-1
       + Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*z
       + 1/2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*x*z^-1
       + Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*x
       - Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)
       + Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*[1-x]^-1*z^-1
       + 2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*[1-x]^-1
       + Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*[1-x]^-1*z
       - 1/2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*z^-1
       - Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*z
       - 1/2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*x*z^-1
       - Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*x
       - 2*Li2(1 - x*z^-1)
       - 4*Li2(1 - x*z^-1)*[1-x]^-1*z^-1*[1-z]^-1
       + 4*Li2(1 - x*z^-1)*[1-x]^-1*z^-1
       + Li2(1 - x*z^-1)*[1-x]^-1*[1-z]^-1
       + 3*Li2(1 - x*z^-1)*[1-x]^-1
       + 3*Li2(1 - x*z^-1)*[1-x]^-1*z
       + 3*Li2(1 - x*z^-1)*z^-1*[1-z]^-1
       - 3*Li2(1 - x*z^-1)*z^-1
       - 1/2*Li2(1 - x*z^-1)*[1-z]^-1
       - 4*Li2(1 - x*z^-1)*z
       + Li2(1 - x*z^-1)*x*z^-1*[1-z]^-1
       - Li2(1 - x*z^-1)*x*z^-1
       - 1/2*Li2(1 - x*z^-1)*x*[1-z]^-1
       - 2*Li2(1 - x*z^-1)*x
       - 1/2*Li2( - x*z^-1)*x^-2*[1+x]^-1*z^-1
       + Li2( - x*z^-1)*x^-2*[1+x]^-1
       + 1/2*Li2( - x*z^-1)*x^-2*z^-1
       - Li2( - x*z^-1)*x^-2
       + Li2( - x*z^-1)*x^-1*[1+x]^-1*z
       - 1/2*Li2( - x*z^-1)*x^-1*z^-1
       + Li2( - x*z^-1)*x^-1
       - Li2( - x*z^-1)*x^-1*z
       - Li2( - x*z^-1)*[1-x]^-1*z^-1
       - 2*Li2( - x*z^-1)*[1-x]^-1
       - Li2( - x*z^-1)*[1-x]^-1*z
       - 1/2*Li2( - x*z^-1)*[1+x]^-1*z^-1
       + Li2( - x*z^-1)*[1+x]^-1
       + Li2( - x*z^-1)*z^-1
       + 2*Li2( - x*z^-1)*z
       + 2*Li2( - x*z^-1)*x
       - 2*Li2( - x)
       + Li2( - x)*x^-2*[1+x]^-1*z^-1
       - 2*Li2( - x)*x^-2*[1+x]^-1
       - Li2( - x)*x^-2*z^-1
       + 2*Li2( - x)*x^-2
       - 2*Li2( - x)*x^-1*[1+x]^-1*z
       + Li2( - x)*x^-1*z^-1
       - 2*Li2( - x)*x^-1
       + 2*Li2( - x)*x^-1*z
       - Li2( - x)*[1+x]^-1*z^-1
       + 2*Li2( - x)*[1+x]^-1
       - 2*Li2( - x)*[1+x]^-1*z
       + Li2( - x)*z^-1
       + Li2( - x)*x*z^-1
       - 2*Li2( - x)*x
       - 2*Li2( - x*z)
       - 1/2*Li2( - x*z)*x^-2*[1+x]^-1*z^-1
       + Li2( - x*z)*x^-2*[1+x]^-1
       + 1/2*Li2( - x*z)*x^-2*z^-1
       - Li2( - x*z)*x^-2
       + Li2( - x*z)*x^-1*[1+x]^-1*z
       - 1/2*Li2( - x*z)*x^-1*z^-1
       + Li2( - x*z)*x^-1
       - Li2( - x*z)*x^-1*z
       + Li2( - x*z)*[1-x]^-1*z^-1
       + 2*Li2( - x*z)*[1-x]^-1
       + Li2( - x*z)*[1-x]^-1*z
       - 1/2*Li2( - x*z)*[1+x]^-1*z^-1
       + Li2( - x*z)*[1+x]^-1
       - Li2( - x*z)*x*z^-1
       - 2*Li2(x)
       - Li2(x)*x^-2*[1+x]^-1*z^-1
       + 2*Li2(x)*x^-2*[1+x]^-1
       + Li2(x)*x^-2*z^-1
       - 2*Li2(x)*x^-2
       + 2*Li2(x)*x^-1*[1+x]^-1*z
       - Li2(x)*x^-1*z^-1
       + 2*Li2(x)*x^-1
       - 2*Li2(x)*x^-1*z
       + 4*Li2(x)*[1-x]^-1*z^-1*[1-z]^-1
       - 4*Li2(x)*[1-x]^-1*z^-1
       - 6*Li2(x)*[1-x]^-1*[1-z]^-1
       + 3/2*Li2(x)*[1-x]^-1
       + 3/2*Li2(x)*[1-x]^-1*z
       + Li2(x)*[1+x]^-1*z^-1
       - 2*Li2(x)*[1+x]^-1
       + 2*Li2(x)*[1+x]^-1*z
       - 2*Li2(x)*z^-1*[1-z]^-1
       + 2*Li2(x)*z^-1
       + 4*Li2(x)*[1-z]^-1
       - 3/2*Li2(x)*z
       - 2*Li2(x)*x*z^-1*[1-z]^-1
       + 2*Li2(x)*x*z^-1
       + 2*Li2(x)*x*[1-z]^-1
       + 1/2*Li2(x)*x
       + Li2(z)
       + 2*Li2(z)*[1-x]^-1*[1-z]^-1
       - 5/2*Li2(z)*[1-x]^-1
       - 5/2*Li2(z)*[1-x]^-1*z
       - 3/2*Li2(z)*[1-z]^-1
       + 7/2*Li2(z)*z
       - 1/2*Li2(z)*x*[1-z]^-1
       + 5/2*Li2(z)*x
       - 3/2*InvTanInt( - sqrtxz3)*x^-1*sqrtxz3
       + 1/2*InvTanInt( - sqrtxz3)*z^-1*sqrtxz3
       - 3/2*InvTanInt( - sqrtxz3)*z*sqrtxz3
       - 3/2*InvTanInt( - sqrtxz3)*x*sqrtxz3
       - 3*InvTanInt(z*sqrtxz3)*x^-1*sqrtxz3
       + InvTanInt(z*sqrtxz3)*z^-1*sqrtxz3
       - 3*InvTanInt(z*sqrtxz3)*z*sqrtxz3
       - 3*InvTanInt(z*sqrtxz3)*x*sqrtxz3
       + 3/2*InvTanInt(sqrtxz3)*x^-1*sqrtxz3
       - 1/2*InvTanInt(sqrtxz3)*z^-1*sqrtxz3
       + 3/2*InvTanInt(sqrtxz3)*z*sqrtxz3
       + 3/2*InvTanInt(sqrtxz3)*x*sqrtxz3
       )

       + NC^-1 * (  - 47/6
       - 3*[1-x]^-1*rln2^2
       - 3*[1-x]^-1*sqrtxz1*rln2
       - 2*[1-x]^-1*z*rln2^2
       - 8*[1-x]^-1*z^2*rln2^2
       + 113/72*z^-1
       + 2*rln2^2
       + 2*sqrtxz1*rln2
       + 11/24*z
       + 65/36*z^2
       + 8*z^2*rln2^2
       - 127/72*x*z^-1
       + 25/6*x
       + 47/24*x*z
       + 4*x*z*rln2^2
       - 31/36*x*z^2
       + 1/6*pi^2*x^-2*[1+x]^-1
       - 1/3*pi^2*x^-2*[1+x]^-1*z
       - 1/6*pi^2*x^-2
       + 1/3*pi^2*x^-2*z
       - 2/3*pi^2*x^-1*[1+x]^-1*z^2
       + 1/6*pi^2*x^-1
       - 1/3*pi^2*x^-1*z
       + 2/3*pi^2*x^-1*z^2
       + 5/12*pi^2*[1-x]^-1
       - 5/6*pi^2*[1-x]^-1*z
       + pi^2*[1-x]^-1*z^2
       - 1/3*pi^2*[1+x]^-1
       + 2/3*pi^2*[1+x]^-1*z
       - pi^2*[1+x]^-1*z^2
       - 1/12*pi^2*z^-1
       + 1/12*pi^2
       + 1/2*pi^2*z
       - 2/3*pi^2*z^2
       - 1/12*pi^2*x*z^-1
       + 5/12*pi^2*x
       - 1/3*pi^2*x*z
       + 4*ln(1 + sqrtxz1 - z)*[1-x]^-1*rln2
       + 3*ln(1 + sqrtxz1 - z)*[1-x]^-1*sqrtxz1
       + 12*ln(1 + sqrtxz1 - z)*[1-x]^-1*z^2*rln2
       - 3*ln(1 + sqrtxz1 - z)*rln2
       - 2*ln(1 + sqrtxz1 - z)*sqrtxz1
       + 2*ln(1 + sqrtxz1 - z)*z*rln2
       - 12*ln(1 + sqrtxz1 - z)*z^2*rln2
       + ln(1 + sqrtxz1 - z)*x*rln2
       - 6*ln(1 + sqrtxz1 - z)*x*z*rln2
       + ln(1 + sqrtxz1 - z)^2
       - ln(1 + sqrtxz1 - z)^2*[1-x]^-1
       + 2*ln(1 + sqrtxz1 - z)^2*[1-x]^-1*z
       - 4*ln(1 + sqrtxz1 - z)^2*[1-x]^-1*z^2
       - 2*ln(1 + sqrtxz1 - z)^2*z
       + 4*ln(1 + sqrtxz1 - z)^2*z^2
       - ln(1 + sqrtxz1 - z)^2*x
       + 2*ln(1 + sqrtxz1 - z)^2*x*z
       + ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)
       - 2*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*[1-x]^-1
       - 4*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z
       - 4*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^2
       + 2*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*z
       + 4*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*z^2
       + ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*x
       + 2*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*x*z
       + 2*ln(1 + sqrtxz1 + z)*[1-x]^-1*rln2
       + 4*ln(1 + sqrtxz1 + z)*[1-x]^-1*z*rln2
       + 4*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^2*rln2
       - ln(1 + sqrtxz1 + z)*rln2
       - 2*ln(1 + sqrtxz1 + z)*z*rln2
       - 4*ln(1 + sqrtxz1 + z)*z^2*rln2
       - ln(1 + sqrtxz1 + z)*x*rln2
       - 2*ln(1 + sqrtxz1 + z)*x*z*rln2
       + 4*ln(z*sqrtxz3)*ArcTan(z*sqrtxz3)*z*sqrtxz3
       - 25/4*ln(x)
       + 2*ln(x)*x^-2*[1+x]^-1
       - 4*ln(x)*x^-2*[1+x]^-1*z
       - 2*ln(x)*x^-2
       + 4*ln(x)*x^-2*z
       - 8*ln(x)*x^-1*[1+x]^-1*z^2
       + 2*ln(x)*x^-1
       - 4*ln(x)*x^-1*z
       + 8*ln(x)*x^-1*z^2
       + 2/3*ln(x)*[1-x]^-1*z^-1
       + 9/4*ln(x)*[1-x]^-1
       - 5/2*ln(x)*[1-x]^-1*rln2
       - 3/2*ln(x)*[1-x]^-1*sqrtxz1
       - 11/4*ln(x)*[1-x]^-1*z
       + ln(x)*[1-x]^-1*z*rln2
       - 2/3*ln(x)*[1-x]^-1*z^2
       - 8*ln(x)*[1-x]^-1*z^2*rln2
       - 2*ln(x)*[1+x]^-1
       + 4*ln(x)*[1+x]^-1*z
       - 8*ln(x)*[1+x]^-1*z^2
       + 7/6*ln(x)*z^-1
       + 2*ln(x)*rln2
       + ln(x)*sqrtxz1
       + 3/2*ln(x)*z
       - 2*ln(x)*z*rln2
       + 4/3*ln(x)*z^2
       + 8*ln(x)*z^2*rln2
       - 11/6*ln(x)*x*z^-1
       + 1/2*ln(x)*x*[x-z]^-1
       + 9/4*ln(x)*x
       - ln(x)*x*rln2
       + 3/2*ln(x)*x*z
       + 4*ln(x)*x*z*rln2
       - 2/3*ln(x)*x*z^2
       - ln(x)*x^2*[x-z]^-1
       - ln(x)*ln(1 - sqrtxz2 + x)*sqrtxz2^-1
       + ln(x)*ln(1 - sqrtxz2 + x)*x*sqrtxz2^-1
       - 2*ln(x)*ln(1 - sqrtxz2 + x)*x*z*sqrtxz2^-1
       - ln(x)*ln(1 - sqrtxz2 + x)*x^2*sqrtxz2^-1
       + ln(x)*ln(1 + sqrtxz2 + x)*sqrtxz2^-1
       - ln(x)*ln(1 + sqrtxz2 + x)*x*sqrtxz2^-1
       + 2*ln(x)*ln(1 + sqrtxz2 + x)*x*z*sqrtxz2^-1
       + ln(x)*ln(1 + sqrtxz2 + x)*x^2*sqrtxz2^-1
       - ln(x)*ln(1 + sqrtxz1 - z)
       + ln(x)*ln(1 + sqrtxz1 - z)*[1-x]^-1
       - 2*ln(x)*ln(1 + sqrtxz1 - z)*[1-x]^-1*z
       + 4*ln(x)*ln(1 + sqrtxz1 - z)*[1-x]^-1*z^2
       + 2*ln(x)*ln(1 + sqrtxz1 - z)*z
       - 4*ln(x)*ln(1 + sqrtxz1 - z)*z^2
       + ln(x)*ln(1 + sqrtxz1 - z)*x
       - 2*ln(x)*ln(1 + sqrtxz1 - z)*x*z
       - ln(x)*ln(1 + sqrtxz1 + z)
       + 3/2*ln(x)*ln(1 + sqrtxz1 + z)*[1-x]^-1
       + ln(x)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z
       + 4*ln(x)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^2
       - 4*ln(x)*ln(1 + sqrtxz1 + z)*z^2
       - 2*ln(x)*ln(1 + sqrtxz1 + z)*x*z
       + 1/2*ln(x)*ln(1 + x*z^-1)
       - 1/2*ln(x)*ln(1 + x*z^-1)*x^-2*[1+x]^-1
       + ln(x)*ln(1 + x*z^-1)*x^-2*[1+x]^-1*z
       + 1/2*ln(x)*ln(1 + x*z^-1)*x^-2
       - ln(x)*ln(1 + x*z^-1)*x^-2*z
       + 2*ln(x)*ln(1 + x*z^-1)*x^-1*[1+x]^-1*z^2
       - 1/2*ln(x)*ln(1 + x*z^-1)*x^-1
       + ln(x)*ln(1 + x*z^-1)*x^-1*z
       - 2*ln(x)*ln(1 + x*z^-1)*x^-1*z^2
       - 1/2*ln(x)*ln(1 + x*z^-1)*[1+x]^-1
       + ln(x)*ln(1 + x*z^-1)*[1+x]^-1*z
       - ln(x)*ln(1 + x*z^-1)*z
       + 2*ln(x)*ln(1 + x*z^-1)*z^2
       - 1/2*ln(x)*ln(1 + x*z^-1)*x
       + ln(x)*ln(1 + x*z^-1)*x*z
       + 3*ln(x)*ln(1 + x)
       - 2*ln(x)*ln(1 + x)*[1+x]^-1
       + 4*ln(x)*ln(1 + x)*[1+x]^-1*z
       - 4*ln(x)*ln(1 + x)*[1+x]^-1*z^2
       - 6*ln(x)*ln(1 + x)*z
       + 4*ln(x)*ln(1 + x)*z^2
       + ln(x)*ln(1 + x)*x
       - 2*ln(x)*ln(1 + x)*x*z
       + ln(x)*ln(1 + x*z)
       - 1/2*ln(x)*ln(1 + x*z)*x^-2*[1+x]^-1
       + ln(x)*ln(1 + x*z)*x^-2*[1+x]^-1*z
       + 1/2*ln(x)*ln(1 + x*z)*x^-2
       - ln(x)*ln(1 + x*z)*x^-2*z
       + 2*ln(x)*ln(1 + x*z)*x^-1*[1+x]^-1*z^2
       - 1/2*ln(x)*ln(1 + x*z)*x^-1
       + ln(x)*ln(1 + x*z)*x^-1*z
       - 2*ln(x)*ln(1 + x*z)*x^-1*z^2
       - ln(x)*ln(1 + x*z)*[1-x]^-1
       - 2*ln(x)*ln(1 + x*z)*[1-x]^-1*z
       - 2*ln(x)*ln(1 + x*z)*[1-x]^-1*z^2
       - 1/2*ln(x)*ln(1 + x*z)*[1+x]^-1
       + ln(x)*ln(1 + x*z)*[1+x]^-1*z
       + 4*ln(x)*ln(1 + x*z)*z^2
       + 2*ln(x)*ln(1 + x*z)*x*z
       - 1/2*ln(x)*ln(z + x)
       + ln(x)*ln(z + x)*[1-x]^-1
       + 2*ln(x)*ln(z + x)*[1-x]^-1*z
       + 2*ln(x)*ln(z + x)*[1-x]^-1*z^2
       - ln(x)*ln(z + x)*z
       - 2*ln(x)*ln(z + x)*z^2
       - 1/2*ln(x)*ln(z + x)*x
       - ln(x)*ln(z + x)*x*z
       - 5/2*ln(x)^2
       - 3/4*ln(x)^2*x^-2*[1+x]^-1
       + 3/2*ln(x)^2*x^-2*[1+x]^-1*z
       + 3/4*ln(x)^2*x^-2
       - 3/2*ln(x)^2*x^-2*z
       + 3*ln(x)^2*x^-1*[1+x]^-1*z^2
       - 3/4*ln(x)^2*x^-1
       + 3/2*ln(x)^2*x^-1*z
       - 3*ln(x)^2*x^-1*z^2
       + ln(x)^2*[1-x]^-1
       - 2*ln(x)^2*[1-x]^-1*z
       + ln(x)^2*[1-x]^-1*z^2
       + 5/4*ln(x)^2*[1+x]^-1
       - 5/2*ln(x)^2*[1+x]^-1*z
       + 4*ln(x)^2*[1+x]^-1*z^2
       + 1/2*ln(x)^2*z^-1
       + 3*ln(x)^2*z
       - 2*ln(x)^2*z^2
       + 1/2*ln(x)^2*x*z^-1
       - 2*ln(x)^2*x
       + 2*ln(x)^2*x*z
       + ln(x)*ln([1-x])
       - ln(x)*ln([1-x])*x^-2*[1+x]^-1
       + 2*ln(x)*ln([1-x])*x^-2*[1+x]^-1*z
       + ln(x)*ln([1-x])*x^-2
       - 2*ln(x)*ln([1-x])*x^-2*z
       + 4*ln(x)*ln([1-x])*x^-1*[1+x]^-1*z^2
       - ln(x)*ln([1-x])*x^-1
       + 2*ln(x)*ln([1-x])*x^-1*z
       - 4*ln(x)*ln([1-x])*x^-1*z^2
       - 3/2*ln(x)*ln([1-x])*[1-x]^-1
       + 3*ln(x)*ln([1-x])*[1-x]^-1*z
       - 4*ln(x)*ln([1-x])*[1-x]^-1*z^2
       + ln(x)*ln([1-x])*[1+x]^-1
       - 2*ln(x)*ln([1-x])*[1+x]^-1*z
       + 4*ln(x)*ln([1-x])*[1+x]^-1*z^2
       - 2*ln(x)*ln([1-x])*z
       + 4*ln(x)*ln([1-x])*z^2
       - ln(x)*ln([1+x])
       + ln(x)*ln([1+x])*x^-2*[1+x]^-1
       - 2*ln(x)*ln([1+x])*x^-2*[1+x]^-1*z
       - ln(x)*ln([1+x])*x^-2
       + 2*ln(x)*ln([1+x])*x^-2*z
       - 4*ln(x)*ln([1+x])*x^-1*[1+x]^-1*z^2
       + ln(x)*ln([1+x])*x^-1
       - 2*ln(x)*ln([1+x])*x^-1*z
       + 4*ln(x)*ln([1+x])*x^-1*z^2
       + ln(x)*ln([1+x])*[1+x]^-1
       - 2*ln(x)*ln([1+x])*[1+x]^-1*z
       + 2*ln(x)*ln([1+x])*z
       - 4*ln(x)*ln([1+x])*z^2
       + ln(x)*ln([1+x])*x
       - 2*ln(x)*ln([1+x])*x*z
       + 1/2*ln(x)*ln(z)*x^-2*[1+x]^-1
       - ln(x)*ln(z)*x^-2*[1+x]^-1*z
       - 1/2*ln(x)*ln(z)*x^-2
       + ln(x)*ln(z)*x^-2*z
       - 2*ln(x)*ln(z)*x^-1*[1+x]^-1*z^2
       + 1/2*ln(x)*ln(z)*x^-1
       - ln(x)*ln(z)*x^-1*z
       + 2*ln(x)*ln(z)*x^-1*z^2
       - 1/2*ln(x)*ln(z)*[1-x]^-1
       + 7/2*ln(x)*ln(z)*[1-x]^-1*z
       - 4*ln(x)*ln(z)*[1-x]^-1*z^2
       + 1/2*ln(x)*ln(z)*[1+x]^-1
       - ln(x)*ln(z)*[1+x]^-1*z
       - 1/2*ln(x)*ln(z)*z^-1
       - 3*ln(x)*ln(z)*z
       + 2*ln(x)*ln(z)*z^2
       - 1/2*ln(x)*ln(z)*x*z^-1
       - ln(x)*ln(z)*x
       + ln(x)*ln(z)*x*z
       + ln(x)*ln([1-z])
       - 1/2*ln(x)*ln([1-z])*z^-1
       - 1/2*ln(x)*ln([1-z])*x*z^-1
       + ln(x)*ln([1-z])*x
       + 13/4*ln([1-x])
       - 13/12*ln([1-x])*z^-1
       - 1/4*ln([1-x])*z
       - 2/3*ln([1-x])*z^2
       + 17/12*ln([1-x])*x*z^-1
       - 9/4*ln([1-x])*x
       - 3/4*ln([1-x])*x*z
       + 1/3*ln([1-x])*x*z^2
       + 1/2*ln([1-x])*ln(z)
       + ln([1-x])*ln(z)*z
       + 1/2*ln([1-x])*ln(z)*x
       + 3/2*ln(z)*[1-x]^-1
       - 7/2*ln(z)*[1-x]^-1*rln2
       - 3/2*ln(z)*[1-x]^-1*sqrtxz1
       + 3/2*ln(z)*[1-x]^-1*z
       - 5*ln(z)*[1-x]^-1*z*rln2
       - 8*ln(z)*[1-x]^-1*z^2*rln2
       - 13/12*ln(z)*z^-1
       - 1/2*ln(z)*[1-z]^-1
       + 2*ln(z)*rln2
       + ln(z)*sqrtxz1
       - 7/4*ln(z)*z
       + 2*ln(z)*z*rln2
       - 2/3*ln(z)*z^2
       + 8*ln(z)*z^2*rln2
       + 17/12*ln(z)*x*z^-1
       + 1/2*ln(z)*x*[1-z]^-1
       - 1/2*ln(z)*x*[x-z]^-1
       - 2*ln(z)*x
       + ln(z)*x*rln2
       - 3/4*ln(z)*x*z
       + 4*ln(z)*x*z*rln2
       + 1/3*ln(z)*x*z^2
       + ln(z)*x^2*[x-z]^-1
       - ln(z)*ln(1 + sqrtxz1 - z)
       + 3/2*ln(z)*ln(1 + sqrtxz1 - z)*[1-x]^-1
       + ln(z)*ln(1 + sqrtxz1 - z)*[1-x]^-1*z
       + 4*ln(z)*ln(1 + sqrtxz1 - z)*[1-x]^-1*z^2
       - 4*ln(z)*ln(1 + sqrtxz1 - z)*z^2
       - 2*ln(z)*ln(1 + sqrtxz1 - z)*x*z
       - ln(z)*ln(1 + sqrtxz1 + z)
       + 2*ln(z)*ln(1 + sqrtxz1 + z)*[1-x]^-1
       + 4*ln(z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z
       + 4*ln(z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^2
       - 2*ln(z)*ln(1 + sqrtxz1 + z)*z
       - 4*ln(z)*ln(1 + sqrtxz1 + z)*z^2
       - ln(z)*ln(1 + sqrtxz1 + z)*x
       - 2*ln(z)*ln(1 + sqrtxz1 + z)*x*z
       - 1/2*ln(z)*ln(1 + x*z^-1)
       + 1/2*ln(z)*ln(1 + x*z^-1)*x^-2*[1+x]^-1
       - ln(z)*ln(1 + x*z^-1)*x^-2*[1+x]^-1*z
       - 1/2*ln(z)*ln(1 + x*z^-1)*x^-2
       + ln(z)*ln(1 + x*z^-1)*x^-2*z
       - 2*ln(z)*ln(1 + x*z^-1)*x^-1*[1+x]^-1*z^2
       + 1/2*ln(z)*ln(1 + x*z^-1)*x^-1
       - ln(z)*ln(1 + x*z^-1)*x^-1*z
       + 2*ln(z)*ln(1 + x*z^-1)*x^-1*z^2
       + 1/2*ln(z)*ln(1 + x*z^-1)*[1+x]^-1
       - ln(z)*ln(1 + x*z^-1)*[1+x]^-1*z
       + ln(z)*ln(1 + x*z^-1)*z
       - 2*ln(z)*ln(1 + x*z^-1)*z^2
       + 1/2*ln(z)*ln(1 + x*z^-1)*x
       - ln(z)*ln(1 + x*z^-1)*x*z
       + ln(z)*ln(1 + x*z)
       - 1/2*ln(z)*ln(1 + x*z)*x^-2*[1+x]^-1
       + ln(z)*ln(1 + x*z)*x^-2*[1+x]^-1*z
       + 1/2*ln(z)*ln(1 + x*z)*x^-2
       - ln(z)*ln(1 + x*z)*x^-2*z
       + 2*ln(z)*ln(1 + x*z)*x^-1*[1+x]^-1*z^2
       - 1/2*ln(z)*ln(1 + x*z)*x^-1
       + ln(z)*ln(1 + x*z)*x^-1*z
       - 2*ln(z)*ln(1 + x*z)*x^-1*z^2
       - ln(z)*ln(1 + x*z)*[1-x]^-1
       - 2*ln(z)*ln(1 + x*z)*[1-x]^-1*z
       - 2*ln(z)*ln(1 + x*z)*[1-x]^-1*z^2
       - 1/2*ln(z)*ln(1 + x*z)*[1+x]^-1
       + ln(z)*ln(1 + x*z)*[1+x]^-1*z
       + 4*ln(z)*ln(1 + x*z)*z^2
       + 2*ln(z)*ln(1 + x*z)*x*z
       + 1/2*ln(z)*ln(z + x)
       - ln(z)*ln(z + x)*[1-x]^-1
       - 2*ln(z)*ln(z + x)*[1-x]^-1*z
       - 2*ln(z)*ln(z + x)*[1-x]^-1*z^2
       + ln(z)*ln(z + x)*z
       + 2*ln(z)*ln(z + x)*z^2
       + 1/2*ln(z)*ln(z + x)*x
       + ln(z)*ln(z + x)*x*z
       - 1/2*ln(z)^2
       + 1/4*ln(z)^2*x^-2*[1+x]^-1
       - 1/2*ln(z)^2*x^-2*[1+x]^-1*z
       - 1/4*ln(z)^2*x^-2
       + 1/2*ln(z)^2*x^-2*z
       - ln(z)^2*x^-1*[1+x]^-1*z^2
       + 1/4*ln(z)^2*x^-1
       - 1/2*ln(z)^2*x^-1*z
       + ln(z)^2*x^-1*z^2
       + 1/2*ln(z)^2*[1-x]^-1
       - ln(z)^2*[1-x]^-1*z
       + ln(z)^2*[1-x]^-1*z^2
       + 1/4*ln(z)^2*[1+x]^-1
       - 1/2*ln(z)^2*[1+x]^-1*z
       + 3*ln(z)^2*z
       - 2*ln(z)^2*z^2
       + ln(z)^2*x
       - ln(z)^2*x*z
       + ln(z)*ln([1-z])
       - 1/2*ln(z)*ln([1-z])*[1-x]^-1
       + ln(z)*ln([1-z])*[1-x]^-1*z
       - 2*ln(z)*ln([1-z])*z
       + 13/4*ln([1-z])
       - 13/12*ln([1-z])*z^-1
       - 1/4*ln([1-z])*z
       - 2/3*ln([1-z])*z^2
       + 17/12*ln([1-z])*x*z^-1
       - 9/4*ln([1-z])*x
       - 3/4*ln([1-z])*x*z
       + 1/3*ln([1-z])*x*z^2
       - 4*ln(sqrtxz3)*ArcTan(sqrtxz3)*z*sqrtxz3
       - Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*sqrtxz2^-1
       + Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x*sqrtxz2^-1
       - 2*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x*z*sqrtxz2^-1
       - Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x^2*sqrtxz2^-1
       + Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*sqrtxz2^-1
       - Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x*sqrtxz2^-1
       + 2*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x*z*sqrtxz2^-1
       + Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x^2*sqrtxz2^-1
       - 1/2*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1
       - 3*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*z
       + 2*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z
       + Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x
       + 1/2*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1
       + 3*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*z
       - 2*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z
       - Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x
       + Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*sqrtxz2^-1
       - Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x*sqrtxz2^-1
       + 2*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x*z*sqrtxz2^-1
       + Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x^2*sqrtxz2^-1
       - Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*sqrtxz2^-1
       + Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x*sqrtxz2^-1
       - 2*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x*z*sqrtxz2^-1
       - Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x^2*sqrtxz2^-1
       - Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)
       + 3/2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*[1-x]^-1
       + Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*[1-x]^-1*z
       + 4*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*[1-x]^-1*z^2
       - 4*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*z^2
       - 2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*x*z
       + Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)
       - 3/2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*[1-x]^-1
       - Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*[1-x]^-1*z
       - 4*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*[1-x]^-1*z^2
       + 4*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*z^2
       + 2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*x*z
       - Li2(1 - x*z^-1)
       + 1/2*Li2(1 - x*z^-1)*[1-x]^-1
       - Li2(1 - x*z^-1)*[1-x]^-1*z
       + 2*Li2(1 - x*z^-1)*z
       - 1/2*Li2( - x*z^-1)*x^-2*[1+x]^-1
       + Li2( - x*z^-1)*x^-2*[1+x]^-1*z
       + 1/2*Li2( - x*z^-1)*x^-2
       - Li2( - x*z^-1)*x^-2*z
       + 2*Li2( - x*z^-1)*x^-1*[1+x]^-1*z^2
       - 1/2*Li2( - x*z^-1)*x^-1
       + Li2( - x*z^-1)*x^-1*z
       - 2*Li2( - x*z^-1)*x^-1*z^2
       + Li2( - x*z^-1)*[1-x]^-1
       + 2*Li2( - x*z^-1)*[1-x]^-1*z
       + 2*Li2( - x*z^-1)*[1-x]^-1*z^2
       - 1/2*Li2( - x*z^-1)*[1+x]^-1
       + Li2( - x*z^-1)*[1+x]^-1*z
       - 2*Li2( - x*z^-1)*z
       - Li2( - x*z^-1)*x
       + 2*Li2( - x)
       + Li2( - x)*x^-2*[1+x]^-1
       - 2*Li2( - x)*x^-2*[1+x]^-1*z
       - Li2( - x)*x^-2
       + 2*Li2( - x)*x^-2*z
       - 4*Li2( - x)*x^-1*[1+x]^-1*z^2
       + Li2( - x)*x^-1
       - 2*Li2( - x)*x^-1*z
       + 4*Li2( - x)*x^-1*z^2
       - Li2( - x)*[1+x]^-1
       + 2*Li2( - x)*[1+x]^-1*z
       - 4*Li2( - x)*[1+x]^-1*z^2
       - 4*Li2( - x)*z
       + 2*Li2( - x)*x
       - 4*Li2( - x)*x*z
       + Li2( - x*z)
       - 1/2*Li2( - x*z)*x^-2*[1+x]^-1
       + Li2( - x*z)*x^-2*[1+x]^-1*z
       + 1/2*Li2( - x*z)*x^-2
       - Li2( - x*z)*x^-2*z
       + 2*Li2( - x*z)*x^-1*[1+x]^-1*z^2
       - 1/2*Li2( - x*z)*x^-1
       + Li2( - x*z)*x^-1*z
       - 2*Li2( - x*z)*x^-1*z^2
       - Li2( - x*z)*[1-x]^-1
       - 2*Li2( - x*z)*[1-x]^-1*z
       - 2*Li2( - x*z)*[1-x]^-1*z^2
       - 1/2*Li2( - x*z)*[1+x]^-1
       + Li2( - x*z)*[1+x]^-1*z
       + 4*Li2( - x*z)*z^2
       + 2*Li2( - x*z)*x*z
       - Li2(x)*x^-2*[1+x]^-1
       + 2*Li2(x)*x^-2*[1+x]^-1*z
       + Li2(x)*x^-2
       - 2*Li2(x)*x^-2*z
       + 4*Li2(x)*x^-1*[1+x]^-1*z^2
       - Li2(x)*x^-1
       + 2*Li2(x)*x^-1*z
       - 4*Li2(x)*x^-1*z^2
       - 3/2*Li2(x)*[1-x]^-1
       + 3*Li2(x)*[1-x]^-1*z
       - 4*Li2(x)*[1-x]^-1*z^2
       + Li2(x)*[1+x]^-1
       - 2*Li2(x)*[1+x]^-1*z
       + 4*Li2(x)*[1+x]^-1*z^2
       + 1/2*Li2(x)*z^-1
       - 2*Li2(x)*z
       + 4*Li2(x)*z^2
       + 1/2*Li2(x)*x*z^-1
       - Li2(x)*x
       + 1/2*Li2(z)
       - 1/2*Li2(z)*[1-x]^-1
       + Li2(z)*[1-x]^-1*z
       - 3*Li2(z)*z
       - 1/2*Li2(z)*x
       - 2*InvTanInt( - sqrtxz3)*z*sqrtxz3
       - 4*InvTanInt(z*sqrtxz3)*z*sqrtxz3
       + 2*InvTanInt(sqrtxz3)*z*sqrtxz3
       )

       + NC^-1*NF * (  - 1/3
       + 5/9*z
       + 8/9*x
       - 1/2*ln(x)*[1-x]^-1
       - 1/2*ln(x)*[1-x]^-1*z
       + 2/3*ln(x)*z
       + 2/3*ln(x)*x
       - 1/3*ln([1-x])*z
       - 1/3*ln([1-x])*x
       - 1/3*ln([1-z])*z
       - 1/3*ln([1-z])*x
       )

       + NC * ( 47/6
       + 3*[1-x]^-1*rln2^2
       + 3*[1-x]^-1*sqrtxz1*rln2
       + 2*[1-x]^-1*z*rln2^2
       + 8*[1-x]^-1*z^2*rln2^2
       - 113/72*z^-1
       - 2*rln2^2
       - 2*sqrtxz1*rln2
       - 11/24*z
       - 65/36*z^2
       - 8*z^2*rln2^2
       + 127/72*x*z^-1
       - 25/6*x
       - 47/24*x*z
       - 4*x*z*rln2^2
       + 31/36*x*z^2
       - 1/6*pi^2*x^-2*[1+x]^-1
       + 1/3*pi^2*x^-2*[1+x]^-1*z
       + 1/6*pi^2*x^-2
       - 1/3*pi^2*x^-2*z
       + 2/3*pi^2*x^-1*[1+x]^-1*z^2
       - 1/6*pi^2*x^-1
       + 1/3*pi^2*x^-1*z
       - 2/3*pi^2*x^-1*z^2
       - 5/12*pi^2*[1-x]^-1
       + 5/6*pi^2*[1-x]^-1*z
       - pi^2*[1-x]^-1*z^2
       + 1/3*pi^2*[1+x]^-1
       - 2/3*pi^2*[1+x]^-1*z
       + pi^2*[1+x]^-1*z^2
       + 1/12*pi^2*z^-1
       - 1/12*pi^2
       - 1/2*pi^2*z
       + 2/3*pi^2*z^2
       + 1/12*pi^2*x*z^-1
       - 5/12*pi^2*x
       + 1/3*pi^2*x*z
       - 4*ln(1 + sqrtxz1 - z)*[1-x]^-1*rln2
       - 3*ln(1 + sqrtxz1 - z)*[1-x]^-1*sqrtxz1
       - 12*ln(1 + sqrtxz1 - z)*[1-x]^-1*z^2*rln2
       + 3*ln(1 + sqrtxz1 - z)*rln2
       + 2*ln(1 + sqrtxz1 - z)*sqrtxz1
       - 2*ln(1 + sqrtxz1 - z)*z*rln2
       + 12*ln(1 + sqrtxz1 - z)*z^2*rln2
       - ln(1 + sqrtxz1 - z)*x*rln2
       + 6*ln(1 + sqrtxz1 - z)*x*z*rln2
       - ln(1 + sqrtxz1 - z)^2
       + ln(1 + sqrtxz1 - z)^2*[1-x]^-1
       - 2*ln(1 + sqrtxz1 - z)^2*[1-x]^-1*z
       + 4*ln(1 + sqrtxz1 - z)^2*[1-x]^-1*z^2
       + 2*ln(1 + sqrtxz1 - z)^2*z
       - 4*ln(1 + sqrtxz1 - z)^2*z^2
       + ln(1 + sqrtxz1 - z)^2*x
       - 2*ln(1 + sqrtxz1 - z)^2*x*z
       - ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)
       + 2*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*[1-x]^-1
       + 4*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z
       + 4*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^2
       - 2*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*z
       - 4*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*z^2
       - ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*x
       - 2*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*x*z
       - 2*ln(1 + sqrtxz1 + z)*[1-x]^-1*rln2
       - 4*ln(1 + sqrtxz1 + z)*[1-x]^-1*z*rln2
       - 4*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^2*rln2
       + ln(1 + sqrtxz1 + z)*rln2
       + 2*ln(1 + sqrtxz1 + z)*z*rln2
       + 4*ln(1 + sqrtxz1 + z)*z^2*rln2
       + ln(1 + sqrtxz1 + z)*x*rln2
       + 2*ln(1 + sqrtxz1 + z)*x*z*rln2
       - 4*ln(z*sqrtxz3)*ArcTan(z*sqrtxz3)*z*sqrtxz3
       + 25/4*ln(x)
       - 2*ln(x)*x^-2*[1+x]^-1
       + 4*ln(x)*x^-2*[1+x]^-1*z
       + 2*ln(x)*x^-2
       - 4*ln(x)*x^-2*z
       + 8*ln(x)*x^-1*[1+x]^-1*z^2
       - 2*ln(x)*x^-1
       + 4*ln(x)*x^-1*z
       - 8*ln(x)*x^-1*z^2
       - 2/3*ln(x)*[1-x]^-1*z^-1
       - 9/4*ln(x)*[1-x]^-1
       + 5/2*ln(x)*[1-x]^-1*rln2
       + 3/2*ln(x)*[1-x]^-1*sqrtxz1
       + 11/4*ln(x)*[1-x]^-1*z
       - ln(x)*[1-x]^-1*z*rln2
       + 2/3*ln(x)*[1-x]^-1*z^2
       + 8*ln(x)*[1-x]^-1*z^2*rln2
       + 2*ln(x)*[1+x]^-1
       - 4*ln(x)*[1+x]^-1*z
       + 8*ln(x)*[1+x]^-1*z^2
       - 7/6*ln(x)*z^-1
       - 2*ln(x)*rln2
       - ln(x)*sqrtxz1
       - 3/2*ln(x)*z
       + 2*ln(x)*z*rln2
       - 4/3*ln(x)*z^2
       - 8*ln(x)*z^2*rln2
       + 11/6*ln(x)*x*z^-1
       - 1/2*ln(x)*x*[x-z]^-1
       - 9/4*ln(x)*x
       + ln(x)*x*rln2
       - 3/2*ln(x)*x*z
       - 4*ln(x)*x*z*rln2
       + 2/3*ln(x)*x*z^2
       + ln(x)*x^2*[x-z]^-1
       + ln(x)*ln(1 - sqrtxz2 + x)*sqrtxz2^-1
       - ln(x)*ln(1 - sqrtxz2 + x)*x*sqrtxz2^-1
       + 2*ln(x)*ln(1 - sqrtxz2 + x)*x*z*sqrtxz2^-1
       + ln(x)*ln(1 - sqrtxz2 + x)*x^2*sqrtxz2^-1
       - ln(x)*ln(1 + sqrtxz2 + x)*sqrtxz2^-1
       + ln(x)*ln(1 + sqrtxz2 + x)*x*sqrtxz2^-1
       - 2*ln(x)*ln(1 + sqrtxz2 + x)*x*z*sqrtxz2^-1
       - ln(x)*ln(1 + sqrtxz2 + x)*x^2*sqrtxz2^-1
       + ln(x)*ln(1 + sqrtxz1 - z)
       - ln(x)*ln(1 + sqrtxz1 - z)*[1-x]^-1
       + 2*ln(x)*ln(1 + sqrtxz1 - z)*[1-x]^-1*z
       - 4*ln(x)*ln(1 + sqrtxz1 - z)*[1-x]^-1*z^2
       - 2*ln(x)*ln(1 + sqrtxz1 - z)*z
       + 4*ln(x)*ln(1 + sqrtxz1 - z)*z^2
       - ln(x)*ln(1 + sqrtxz1 - z)*x
       + 2*ln(x)*ln(1 + sqrtxz1 - z)*x*z
       + ln(x)*ln(1 + sqrtxz1 + z)
       - 3/2*ln(x)*ln(1 + sqrtxz1 + z)*[1-x]^-1
       - ln(x)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z
       - 4*ln(x)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^2
       + 4*ln(x)*ln(1 + sqrtxz1 + z)*z^2
       + 2*ln(x)*ln(1 + sqrtxz1 + z)*x*z
       - 1/2*ln(x)*ln(1 + x*z^-1)
       + 1/2*ln(x)*ln(1 + x*z^-1)*x^-2*[1+x]^-1
       - ln(x)*ln(1 + x*z^-1)*x^-2*[1+x]^-1*z
       - 1/2*ln(x)*ln(1 + x*z^-1)*x^-2
       + ln(x)*ln(1 + x*z^-1)*x^-2*z
       - 2*ln(x)*ln(1 + x*z^-1)*x^-1*[1+x]^-1*z^2
       + 1/2*ln(x)*ln(1 + x*z^-1)*x^-1
       - ln(x)*ln(1 + x*z^-1)*x^-1*z
       + 2*ln(x)*ln(1 + x*z^-1)*x^-1*z^2
       + 1/2*ln(x)*ln(1 + x*z^-1)*[1+x]^-1
       - ln(x)*ln(1 + x*z^-1)*[1+x]^-1*z
       + ln(x)*ln(1 + x*z^-1)*z
       - 2*ln(x)*ln(1 + x*z^-1)*z^2
       + 1/2*ln(x)*ln(1 + x*z^-1)*x
       - ln(x)*ln(1 + x*z^-1)*x*z
       - 3*ln(x)*ln(1 + x)
       + 2*ln(x)*ln(1 + x)*[1+x]^-1
       - 4*ln(x)*ln(1 + x)*[1+x]^-1*z
       + 4*ln(x)*ln(1 + x)*[1+x]^-1*z^2
       + 6*ln(x)*ln(1 + x)*z
       - 4*ln(x)*ln(1 + x)*z^2
       - ln(x)*ln(1 + x)*x
       + 2*ln(x)*ln(1 + x)*x*z
       - ln(x)*ln(1 + x*z)
       + 1/2*ln(x)*ln(1 + x*z)*x^-2*[1+x]^-1
       - ln(x)*ln(1 + x*z)*x^-2*[1+x]^-1*z
       - 1/2*ln(x)*ln(1 + x*z)*x^-2
       + ln(x)*ln(1 + x*z)*x^-2*z
       - 2*ln(x)*ln(1 + x*z)*x^-1*[1+x]^-1*z^2
       + 1/2*ln(x)*ln(1 + x*z)*x^-1
       - ln(x)*ln(1 + x*z)*x^-1*z
       + 2*ln(x)*ln(1 + x*z)*x^-1*z^2
       + ln(x)*ln(1 + x*z)*[1-x]^-1
       + 2*ln(x)*ln(1 + x*z)*[1-x]^-1*z
       + 2*ln(x)*ln(1 + x*z)*[1-x]^-1*z^2
       + 1/2*ln(x)*ln(1 + x*z)*[1+x]^-1
       - ln(x)*ln(1 + x*z)*[1+x]^-1*z
       - 4*ln(x)*ln(1 + x*z)*z^2
       - 2*ln(x)*ln(1 + x*z)*x*z
       + 1/2*ln(x)*ln(z + x)
       - ln(x)*ln(z + x)*[1-x]^-1
       - 2*ln(x)*ln(z + x)*[1-x]^-1*z
       - 2*ln(x)*ln(z + x)*[1-x]^-1*z^2
       + ln(x)*ln(z + x)*z
       + 2*ln(x)*ln(z + x)*z^2
       + 1/2*ln(x)*ln(z + x)*x
       + ln(x)*ln(z + x)*x*z
       + 5/2*ln(x)^2
       + 3/4*ln(x)^2*x^-2*[1+x]^-1
       - 3/2*ln(x)^2*x^-2*[1+x]^-1*z
       - 3/4*ln(x)^2*x^-2
       + 3/2*ln(x)^2*x^-2*z
       - 3*ln(x)^2*x^-1*[1+x]^-1*z^2
       + 3/4*ln(x)^2*x^-1
       - 3/2*ln(x)^2*x^-1*z
       + 3*ln(x)^2*x^-1*z^2
       - ln(x)^2*[1-x]^-1
       + 2*ln(x)^2*[1-x]^-1*z
       - ln(x)^2*[1-x]^-1*z^2
       - 5/4*ln(x)^2*[1+x]^-1
       + 5/2*ln(x)^2*[1+x]^-1*z
       - 4*ln(x)^2*[1+x]^-1*z^2
       - 1/2*ln(x)^2*z^-1
       - 3*ln(x)^2*z
       + 2*ln(x)^2*z^2
       - 1/2*ln(x)^2*x*z^-1
       + 2*ln(x)^2*x
       - 2*ln(x)^2*x*z
       - ln(x)*ln([1-x])
       + ln(x)*ln([1-x])*x^-2*[1+x]^-1
       - 2*ln(x)*ln([1-x])*x^-2*[1+x]^-1*z
       - ln(x)*ln([1-x])*x^-2
       + 2*ln(x)*ln([1-x])*x^-2*z
       - 4*ln(x)*ln([1-x])*x^-1*[1+x]^-1*z^2
       + ln(x)*ln([1-x])*x^-1
       - 2*ln(x)*ln([1-x])*x^-1*z
       + 4*ln(x)*ln([1-x])*x^-1*z^2
       + 3/2*ln(x)*ln([1-x])*[1-x]^-1
       - 3*ln(x)*ln([1-x])*[1-x]^-1*z
       + 4*ln(x)*ln([1-x])*[1-x]^-1*z^2
       - ln(x)*ln([1-x])*[1+x]^-1
       + 2*ln(x)*ln([1-x])*[1+x]^-1*z
       - 4*ln(x)*ln([1-x])*[1+x]^-1*z^2
       + 2*ln(x)*ln([1-x])*z
       - 4*ln(x)*ln([1-x])*z^2
       + ln(x)*ln([1+x])
       - ln(x)*ln([1+x])*x^-2*[1+x]^-1
       + 2*ln(x)*ln([1+x])*x^-2*[1+x]^-1*z
       + ln(x)*ln([1+x])*x^-2
       - 2*ln(x)*ln([1+x])*x^-2*z
       + 4*ln(x)*ln([1+x])*x^-1*[1+x]^-1*z^2
       - ln(x)*ln([1+x])*x^-1
       + 2*ln(x)*ln([1+x])*x^-1*z
       - 4*ln(x)*ln([1+x])*x^-1*z^2
       - ln(x)*ln([1+x])*[1+x]^-1
       + 2*ln(x)*ln([1+x])*[1+x]^-1*z
       - 2*ln(x)*ln([1+x])*z
       + 4*ln(x)*ln([1+x])*z^2
       - ln(x)*ln([1+x])*x
       + 2*ln(x)*ln([1+x])*x*z
       - 1/2*ln(x)*ln(z)*x^-2*[1+x]^-1
       + ln(x)*ln(z)*x^-2*[1+x]^-1*z
       + 1/2*ln(x)*ln(z)*x^-2
       - ln(x)*ln(z)*x^-2*z
       + 2*ln(x)*ln(z)*x^-1*[1+x]^-1*z^2
       - 1/2*ln(x)*ln(z)*x^-1
       + ln(x)*ln(z)*x^-1*z
       - 2*ln(x)*ln(z)*x^-1*z^2
       + 1/2*ln(x)*ln(z)*[1-x]^-1
       - 7/2*ln(x)*ln(z)*[1-x]^-1*z
       + 4*ln(x)*ln(z)*[1-x]^-1*z^2
       - 1/2*ln(x)*ln(z)*[1+x]^-1
       + ln(x)*ln(z)*[1+x]^-1*z
       + 1/2*ln(x)*ln(z)*z^-1
       + 3*ln(x)*ln(z)*z
       - 2*ln(x)*ln(z)*z^2
       + 1/2*ln(x)*ln(z)*x*z^-1
       + ln(x)*ln(z)*x
       - ln(x)*ln(z)*x*z
       - ln(x)*ln([1-z])
       + 1/2*ln(x)*ln([1-z])*z^-1
       + 1/2*ln(x)*ln([1-z])*x*z^-1
       - ln(x)*ln([1-z])*x
       - 13/4*ln([1-x])
       + 13/12*ln([1-x])*z^-1
       + 1/4*ln([1-x])*z
       + 2/3*ln([1-x])*z^2
       - 17/12*ln([1-x])*x*z^-1
       + 9/4*ln([1-x])*x
       + 3/4*ln([1-x])*x*z
       - 1/3*ln([1-x])*x*z^2
       - 1/2*ln([1-x])*ln(z)
       - ln([1-x])*ln(z)*z
       - 1/2*ln([1-x])*ln(z)*x
       - 3/2*ln(z)*[1-x]^-1
       + 7/2*ln(z)*[1-x]^-1*rln2
       + 3/2*ln(z)*[1-x]^-1*sqrtxz1
       - 3/2*ln(z)*[1-x]^-1*z
       + 5*ln(z)*[1-x]^-1*z*rln2
       + 8*ln(z)*[1-x]^-1*z^2*rln2
       + 13/12*ln(z)*z^-1
       + 1/2*ln(z)*[1-z]^-1
       - 2*ln(z)*rln2
       - ln(z)*sqrtxz1
       + 7/4*ln(z)*z
       - 2*ln(z)*z*rln2
       + 2/3*ln(z)*z^2
       - 8*ln(z)*z^2*rln2
       - 17/12*ln(z)*x*z^-1
       - 1/2*ln(z)*x*[1-z]^-1
       + 1/2*ln(z)*x*[x-z]^-1
       + 2*ln(z)*x
       - ln(z)*x*rln2
       + 3/4*ln(z)*x*z
       - 4*ln(z)*x*z*rln2
       - 1/3*ln(z)*x*z^2
       - ln(z)*x^2*[x-z]^-1
       + ln(z)*ln(1 + sqrtxz1 - z)
       - 3/2*ln(z)*ln(1 + sqrtxz1 - z)*[1-x]^-1
       - ln(z)*ln(1 + sqrtxz1 - z)*[1-x]^-1*z
       - 4*ln(z)*ln(1 + sqrtxz1 - z)*[1-x]^-1*z^2
       + 4*ln(z)*ln(1 + sqrtxz1 - z)*z^2
       + 2*ln(z)*ln(1 + sqrtxz1 - z)*x*z
       + ln(z)*ln(1 + sqrtxz1 + z)
       - 2*ln(z)*ln(1 + sqrtxz1 + z)*[1-x]^-1
       - 4*ln(z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z
       - 4*ln(z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^2
       + 2*ln(z)*ln(1 + sqrtxz1 + z)*z
       + 4*ln(z)*ln(1 + sqrtxz1 + z)*z^2
       + ln(z)*ln(1 + sqrtxz1 + z)*x
       + 2*ln(z)*ln(1 + sqrtxz1 + z)*x*z
       + 1/2*ln(z)*ln(1 + x*z^-1)
       - 1/2*ln(z)*ln(1 + x*z^-1)*x^-2*[1+x]^-1
       + ln(z)*ln(1 + x*z^-1)*x^-2*[1+x]^-1*z
       + 1/2*ln(z)*ln(1 + x*z^-1)*x^-2
       - ln(z)*ln(1 + x*z^-1)*x^-2*z
       + 2*ln(z)*ln(1 + x*z^-1)*x^-1*[1+x]^-1*z^2
       - 1/2*ln(z)*ln(1 + x*z^-1)*x^-1
       + ln(z)*ln(1 + x*z^-1)*x^-1*z
       - 2*ln(z)*ln(1 + x*z^-1)*x^-1*z^2
       - 1/2*ln(z)*ln(1 + x*z^-1)*[1+x]^-1
       + ln(z)*ln(1 + x*z^-1)*[1+x]^-1*z
       - ln(z)*ln(1 + x*z^-1)*z
       + 2*ln(z)*ln(1 + x*z^-1)*z^2
       - 1/2*ln(z)*ln(1 + x*z^-1)*x
       + ln(z)*ln(1 + x*z^-1)*x*z
       - ln(z)*ln(1 + x*z)
       + 1/2*ln(z)*ln(1 + x*z)*x^-2*[1+x]^-1
       - ln(z)*ln(1 + x*z)*x^-2*[1+x]^-1*z
       - 1/2*ln(z)*ln(1 + x*z)*x^-2
       + ln(z)*ln(1 + x*z)*x^-2*z
       - 2*ln(z)*ln(1 + x*z)*x^-1*[1+x]^-1*z^2
       + 1/2*ln(z)*ln(1 + x*z)*x^-1
       - ln(z)*ln(1 + x*z)*x^-1*z
       + 2*ln(z)*ln(1 + x*z)*x^-1*z^2
       + ln(z)*ln(1 + x*z)*[1-x]^-1
       + 2*ln(z)*ln(1 + x*z)*[1-x]^-1*z
       + 2*ln(z)*ln(1 + x*z)*[1-x]^-1*z^2
       + 1/2*ln(z)*ln(1 + x*z)*[1+x]^-1
       - ln(z)*ln(1 + x*z)*[1+x]^-1*z
       - 4*ln(z)*ln(1 + x*z)*z^2
       - 2*ln(z)*ln(1 + x*z)*x*z
       - 1/2*ln(z)*ln(z + x)
       + ln(z)*ln(z + x)*[1-x]^-1
       + 2*ln(z)*ln(z + x)*[1-x]^-1*z
       + 2*ln(z)*ln(z + x)*[1-x]^-1*z^2
       - ln(z)*ln(z + x)*z
       - 2*ln(z)*ln(z + x)*z^2
       - 1/2*ln(z)*ln(z + x)*x
       - ln(z)*ln(z + x)*x*z
       + 1/2*ln(z)^2
       - 1/4*ln(z)^2*x^-2*[1+x]^-1
       + 1/2*ln(z)^2*x^-2*[1+x]^-1*z
       + 1/4*ln(z)^2*x^-2
       - 1/2*ln(z)^2*x^-2*z
       + ln(z)^2*x^-1*[1+x]^-1*z^2
       - 1/4*ln(z)^2*x^-1
       + 1/2*ln(z)^2*x^-1*z
       - ln(z)^2*x^-1*z^2
       - 1/2*ln(z)^2*[1-x]^-1
       + ln(z)^2*[1-x]^-1*z
       - ln(z)^2*[1-x]^-1*z^2
       - 1/4*ln(z)^2*[1+x]^-1
       + 1/2*ln(z)^2*[1+x]^-1*z
       - 3*ln(z)^2*z
       + 2*ln(z)^2*z^2
       - ln(z)^2*x
       + ln(z)^2*x*z
       - ln(z)*ln([1-z])
       + 1/2*ln(z)*ln([1-z])*[1-x]^-1
       - ln(z)*ln([1-z])*[1-x]^-1*z
       + 2*ln(z)*ln([1-z])*z
       - 13/4*ln([1-z])
       + 13/12*ln([1-z])*z^-1
       + 1/4*ln([1-z])*z
       + 2/3*ln([1-z])*z^2
       - 17/12*ln([1-z])*x*z^-1
       + 9/4*ln([1-z])*x
       + 3/4*ln([1-z])*x*z
       - 1/3*ln([1-z])*x*z^2
       + 4*ln(sqrtxz3)*ArcTan(sqrtxz3)*z*sqrtxz3
       + Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*sqrtxz2^-1
       - Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x*sqrtxz2^-1
       + 2*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x*z*sqrtxz2^-1
       + Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x^2*sqrtxz2^-1
       - Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*sqrtxz2^-1
       + Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x*sqrtxz2^-1
       - 2*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x*z*sqrtxz2^-1
       - Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x^2*sqrtxz2^-1
       + 1/2*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1
       + 3*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*z
       - 2*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z
       - Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x
       - 1/2*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1
       - 3*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*z
       + 2*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z
       + Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x
       - Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*sqrtxz2^-1
       + Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x*sqrtxz2^-1
       - 2*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x*z*sqrtxz2^-1
       - Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x^2*sqrtxz2^-1
       + Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*sqrtxz2^-1
       - Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x*sqrtxz2^-1
       + 2*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x*z*sqrtxz2^-1
       + Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x^2*sqrtxz2^-1
       + Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)
       - 3/2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*[1-x]^-1
       - Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*[1-x]^-1*z
       - 4*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*[1-x]^-1*z^2
       + 4*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*z^2
       + 2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*x*z
       - Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)
       + 3/2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*[1-x]^-1
       + Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*[1-x]^-1*z
       + 4*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*[1-x]^-1*z^2
       - 4*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*z^2
       - 2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*x*z
       + Li2(1 - x*z^-1)
       - 1/2*Li2(1 - x*z^-1)*[1-x]^-1
       + Li2(1 - x*z^-1)*[1-x]^-1*z
       - 2*Li2(1 - x*z^-1)*z
       + 1/2*Li2( - x*z^-1)*x^-2*[1+x]^-1
       - Li2( - x*z^-1)*x^-2*[1+x]^-1*z
       - 1/2*Li2( - x*z^-1)*x^-2
       + Li2( - x*z^-1)*x^-2*z
       - 2*Li2( - x*z^-1)*x^-1*[1+x]^-1*z^2
       + 1/2*Li2( - x*z^-1)*x^-1
       - Li2( - x*z^-1)*x^-1*z
       + 2*Li2( - x*z^-1)*x^-1*z^2
       - Li2( - x*z^-1)*[1-x]^-1
       - 2*Li2( - x*z^-1)*[1-x]^-1*z
       - 2*Li2( - x*z^-1)*[1-x]^-1*z^2
       + 1/2*Li2( - x*z^-1)*[1+x]^-1
       - Li2( - x*z^-1)*[1+x]^-1*z
       + 2*Li2( - x*z^-1)*z
       + Li2( - x*z^-1)*x
       - 2*Li2( - x)
       - Li2( - x)*x^-2*[1+x]^-1
       + 2*Li2( - x)*x^-2*[1+x]^-1*z
       + Li2( - x)*x^-2
       - 2*Li2( - x)*x^-2*z
       + 4*Li2( - x)*x^-1*[1+x]^-1*z^2
       - Li2( - x)*x^-1
       + 2*Li2( - x)*x^-1*z
       - 4*Li2( - x)*x^-1*z^2
       + Li2( - x)*[1+x]^-1
       - 2*Li2( - x)*[1+x]^-1*z
       + 4*Li2( - x)*[1+x]^-1*z^2
       + 4*Li2( - x)*z
       - 2*Li2( - x)*x
       + 4*Li2( - x)*x*z
       - Li2( - x*z)
       + 1/2*Li2( - x*z)*x^-2*[1+x]^-1
       - Li2( - x*z)*x^-2*[1+x]^-1*z
       - 1/2*Li2( - x*z)*x^-2
       + Li2( - x*z)*x^-2*z
       - 2*Li2( - x*z)*x^-1*[1+x]^-1*z^2
       + 1/2*Li2( - x*z)*x^-1
       - Li2( - x*z)*x^-1*z
       + 2*Li2( - x*z)*x^-1*z^2
       + Li2( - x*z)*[1-x]^-1
       + 2*Li2( - x*z)*[1-x]^-1*z
       + 2*Li2( - x*z)*[1-x]^-1*z^2
       + 1/2*Li2( - x*z)*[1+x]^-1
       - Li2( - x*z)*[1+x]^-1*z
       - 4*Li2( - x*z)*z^2
       - 2*Li2( - x*z)*x*z
       + Li2(x)*x^-2*[1+x]^-1
       - 2*Li2(x)*x^-2*[1+x]^-1*z
       - Li2(x)*x^-2
       + 2*Li2(x)*x^-2*z
       - 4*Li2(x)*x^-1*[1+x]^-1*z^2
       + Li2(x)*x^-1
       - 2*Li2(x)*x^-1*z
       + 4*Li2(x)*x^-1*z^2
       + 3/2*Li2(x)*[1-x]^-1
       - 3*Li2(x)*[1-x]^-1*z
       + 4*Li2(x)*[1-x]^-1*z^2
       - Li2(x)*[1+x]^-1
       + 2*Li2(x)*[1+x]^-1*z
       - 4*Li2(x)*[1+x]^-1*z^2
       - 1/2*Li2(x)*z^-1
       + 2*Li2(x)*z
       - 4*Li2(x)*z^2
       - 1/2*Li2(x)*x*z^-1
       + Li2(x)*x
       - 1/2*Li2(z)
       + 1/2*Li2(z)*[1-x]^-1
       - Li2(z)*[1-x]^-1*z
       + 3*Li2(z)*z
       + 1/2*Li2(z)*x
       + 2*InvTanInt( - sqrtxz3)*z*sqrtxz3
       + 4*InvTanInt(z*sqrtxz3)*z*sqrtxz3
       - 2*InvTanInt(sqrtxz3)*z*sqrtxz3
       )

       + NC*NF * ( 1/3
       - 5/9*z
       - 8/9*x
       + 1/2*ln(x)*[1-x]^-1
       + 1/2*ln(x)*[1-x]^-1*z
       - 2/3*ln(x)*z
       - 2/3*ln(x)*x
       + 1/3*ln([1-x])*z
       + 1/3*ln([1-x])*x
       + 1/3*ln([1-z])*z
       + 1/3*ln([1-z])*x
       )

       + NC^2 * ( 5/3
       + 2*[1-x]^-1*[1-z]^-1
       + 1/2*[1-x]^-1
       - 5/2*[1-x]^-1*z
       - 7/2*[1-z]^-1
       + 1/4*z*[1-x-z]^-1
       + 20/9*z
       + 3/2*x*[1-z]^-1
       + 1/4*x*[x-z]^-1
       - 13/9*x
       - 1/4*x*z
       - 1/2*x^2*[x-z]^-1
       - 1/2*pi^2*[1-x]^-1*[1-z]^-1
       + 1/4*pi^2*[1-x]^-1
       + 1/4*pi^2*[1-x]^-1*z
       + 1/4*pi^2*[1-z]^-1
       - 1/12*pi^2
       - 13/12*pi^2*z
       + 1/4*pi^2*x*[1-z]^-1
       - 13/12*pi^2*x
       - 1/12*pi^2*x*z
       - 1/4*ln(x)
       - 9*ln(x)*[1-x]^-1*[1-z]^-1
       - 3/4*ln(x)*[1-x]^-1*[x-z]^-1
       + 7/4*ln(x)*[1-x]^-1
       + 17/4*ln(x)*[1-x]^-1*z
       + 9*ln(x)*[1-z]^-1
       - 3/4*ln(x)*[1-x-z]^-1
       + 3/4*ln(x)*[x-z]^-1
       + 1/4*ln(x)*z*[1-x-z]^-2
       - 1/4*ln(x)*z*[1-x-z]^-1
       - 67/12*ln(x)*z
       - 1/4*ln(x)*z^2*[1-x-z]^-2
       + ln(x)*[1+x]
       - ln(x)*[1+x]*z
       + 3/4*ln(x)*x*[x-z]^-1
       - 79/12*ln(x)*x
       + 3/4*ln(x)*x*z
       - 1/4*ln(x)*x^2*[x-z]^-2
       + 1/2*ln(x)*x^2*[x-z]^-1
       + 1/2*ln(x)*x^3*[x-z]^-2
       + 1/4*ln(x)^2
       + 4*ln(x)^2*[1-x]^-1*[1-z]^-1
       - 5/2*ln(x)^2*[1-x]^-1
       - 5/2*ln(x)^2*[1-x]^-1*z
       - 2*ln(x)^2*[1-z]^-1
       + 17/4*ln(x)^2*z
       - 2*ln(x)^2*x*[1-z]^-1
       + 17/4*ln(x)^2*x
       + 1/4*ln(x)^2*x*z
       - 1/2*ln(x)*ln([1-x])
       - 5*ln(x)*ln([1-x])*[1-x]^-1*[1-z]^-1
       + 4*ln(x)*ln([1-x])*[1-x]^-1
       + 4*ln(x)*ln([1-x])*[1-x]^-1*z
       + 5/2*ln(x)*ln([1-x])*[1-z]^-1
       - 13/2*ln(x)*ln([1-x])*z
       + 5/2*ln(x)*ln([1-x])*x*[1-z]^-1
       - 13/2*ln(x)*ln([1-x])*x
       - 1/2*ln(x)*ln([1-x])*x*z
       - 4*ln(x)*ln(z)*[1-x]^-1*[1-z]^-1
       + 7/4*ln(x)*ln(z)*[1-x]^-1
       + 7/4*ln(x)*ln(z)*[1-x]^-1*z
       + 2*ln(x)*ln(z)*[1-z]^-1
       - 5/2*ln(x)*ln(z)*z
       + 2*ln(x)*ln(z)*x*[1-z]^-1
       - 5/2*ln(x)*ln(z)*x
       - 1/2*ln(x)*ln([1-z])
       - 6*ln(x)*ln([1-z])*[1-x]^-1*[1-z]^-1
       + 4*ln(x)*ln([1-z])*[1-x]^-1
       + 4*ln(x)*ln([1-z])*[1-x]^-1*z
       + 3*ln(x)*ln([1-z])*[1-z]^-1
       - 7*ln(x)*ln([1-z])*z
       + 3*ln(x)*ln([1-z])*x*[1-z]^-1
       - 7*ln(x)*ln([1-z])*x
       - 1/2*ln(x)*ln([1-z])*x*z
       - 3/4*ln([1-x])
       + 6*ln([1-x])*[1-x]^-1*[1-z]^-1
       - 2*ln([1-x])*[1-x]^-1
       - 4*ln([1-x])*[1-x]^-1*z
       - 6*ln([1-x])*[1-z]^-1
       + 47/12*ln([1-x])*z
       + 47/12*ln([1-x])*x
       + 1/4*ln([1-x])*x*z
       + 1/4*ln([1-x])^2
       + 2*ln([1-x])^2*[1-x]^-1*[1-z]^-1
       - ln([1-x])^2*[1-x]^-1
       - ln([1-x])^2*[1-x]^-1*z
       - ln([1-x])^2*[1-z]^-1
       + 11/4*ln([1-x])^2*z
       - ln([1-x])^2*x*[1-z]^-1
       + 11/4*ln([1-x])^2*x
       + 1/4*ln([1-x])^2*x*z
       + 2*ln([1-x])*ln(z)*[1-x]^-1*[1-z]^-1
       - ln([1-x])*ln(z)*[1-x]^-1
       - ln([1-x])*ln(z)*[1-x]^-1*z
       - ln([1-x])*ln(z)*[1-z]^-1
       + 3/2*ln([1-x])*ln(z)*z
       - ln([1-x])*ln(z)*x*[1-z]^-1
       + 3/2*ln([1-x])*ln(z)*x
       + 1/2*ln([1-x])*ln([1-z])
       + 4*ln([1-x])*ln([1-z])*[1-x]^-1*[1-z]^-1
       - 2*ln([1-x])*ln([1-z])*[1-x]^-1
       - 2*ln([1-x])*ln([1-z])*[1-x]^-1*z
       - 2*ln([1-x])*ln([1-z])*[1-z]^-1
       + 11/2*ln([1-x])*ln([1-z])*z
       - 2*ln([1-x])*ln([1-z])*x*[1-z]^-1
       + 11/2*ln([1-x])*ln([1-z])*x
       + 1/2*ln([1-x])*ln([1-z])*x*z
       + 1/2*ln(z)
       + 9/4*ln(z)*[1-x]^-1*[1-z]^-1
       + 3/4*ln(z)*[1-x]^-1*[x-z]^-1
       - ln(z)*[1-x]^-1
       - 2*ln(z)*[1-x]^-1*z
       - 5/2*ln(z)*[1-z]^-1
       - 3/4*ln(z)*[x-z]^-1
       + 9/2*ln(z)*z
       - 1/2*ln(z)*x*[1-z]^-1
       - 3/4*ln(z)*x*[x-z]^-1
       + 4*ln(z)*x
       + 1/4*ln(z)*x^2*[x-z]^-2
       - 1/2*ln(z)*x^2*[x-z]^-1
       - 1/2*ln(z)*x^3*[x-z]^-2
       - 1/4*ln(z)^2
       + ln(z)^2*[1-x]^-1*[1-z]^-1
       - 1/2*ln(z)^2*[1-x]^-1
       - 1/2*ln(z)^2*[1-x]^-1*z
       - 1/4*ln(z)^2*z
       - 1/4*ln(z)^2*x
       - 1/4*ln(z)^2*x*z
       + ln(z)*ln([1-z])*[1-x]^-1*[1-z]^-1
       - 1/2*ln(z)*ln([1-z])*[1-x]^-1
       - 1/2*ln(z)*ln([1-z])*[1-x]^-1*z
       - ln(z)*ln([1-z])*[1-z]^-1
       + 3*ln(z)*ln([1-z])*z
       - ln(z)*ln([1-z])*x*[1-z]^-1
       + 3*ln(z)*ln([1-z])*x
       + 6*ln([1-z])*[1-x]^-1*[1-z]^-1
       - 2*ln([1-z])*[1-x]^-1
       - 4*ln([1-z])*[1-x]^-1*z
       - 6*ln([1-z])*[1-z]^-1
       + 3/4*ln([1-z])*[1-x-z]^-1
       - 1/4*ln([1-z])*z*[1-x-z]^-2
       + 1/4*ln([1-z])*z*[1-x-z]^-1
       + 19/6*ln([1-z])*z
       + 1/4*ln([1-z])*z^2*[1-x-z]^-2
       + 19/6*ln([1-z])*x
       + ln([1-z])*x*z
       + 1/4*ln([1-z])^2
       + 2*ln([1-z])^2*[1-x]^-1*[1-z]^-1
       - ln([1-z])^2*[1-x]^-1
       - ln([1-z])^2*[1-x]^-1*z
       - ln([1-z])^2*[1-z]^-1
       + 11/4*ln([1-z])^2*z
       - ln([1-z])^2*x*[1-z]^-1
       + 11/4*ln([1-z])^2*x
       + 1/4*ln([1-z])^2*x*z
       + Li2(1 - x*z^-1)*[1-x]^-1*[1-z]^-1
       - 1/2*Li2(1 - x*z^-1)*[1-x]^-1
       - 1/2*Li2(1 - x*z^-1)*[1-x]^-1*z
       - 1/2*Li2(1 - x*z^-1)*[1-z]^-1
       - 1/2*Li2(1 - x*z^-1)*x*[1-z]^-1
       + 1/2*Li2(x)*z
       + 1/2*Li2(x)*x
       + 3/2*Li2(z)*z
       + 3/2*Li2(z)*x
       )

       + LMUR*NC^-1*NF * ( 1/3*z
       + 1/3*x
       )

       + LMUR * (  - 11/6*z
       - 11/6*x
       )

       + LMUR*NC*NF * (  - 1/3*z
       - 1/3*x
       )

       + LMUR*NC^2 * ( 11/6*z
       + 11/6*x
       )

       + LMUF*NC^-2 * (  - 1/2
       - 3/4*z
       + 1/4*x
       - 1/2*x*z
       + 1/4*ln(x)
       - 1/2*ln(x)*[1-x]^-1
       - 1/2*ln(x)*[1-x]^-1*z
       + 3/4*ln(x)*z
       + 3/4*ln(x)*x
       + 1/4*ln(x)*x*z
       - 1/4*ln([1-x])
       - 5/4*ln([1-x])*z
       - 5/4*ln([1-x])*x
       - 1/4*ln([1-x])*x*z
       - 1/4*ln(z)
       + 1/2*ln(z)*[1-z]^-1
       - 1/4*ln(z)*z
       + 1/2*ln(z)*x*[1-z]^-1
       - 1/4*ln(z)*x
       - 1/4*ln(z)*x*z
       - 1/4*ln([1-z])
       - 1/4*ln([1-z])*z
       - 1/4*ln([1-z])*x
       - 1/4*ln([1-z])*x*z
       )

       + LMUF*NC^-1 * (  - 5/2
       + 5/4*z^-1
       - 5/4*x*z^-1
       + 5/2*x
       - ln(x)
       + 1/2*ln(x)*z^-1
       + 1/2*ln(x)*x*z^-1
       - ln(x)*x
       )

       + LMUF * ( 1
       + 3/2*z
       - 1/2*x
       + x*z
       - 1/2*ln(x)
       + ln(x)*[1-x]^-1
       + ln(x)*[1-x]^-1*z
       - 3/2*ln(x)*z
       - 3/2*ln(x)*x
       - 1/2*ln(x)*x*z
       + 1/2*ln([1-x])
       + 5/2*ln([1-x])*z
       + 5/2*ln([1-x])*x
       + 1/2*ln([1-x])*x*z
       + 1/2*ln(z)
       - ln(z)*[1-z]^-1
       + 1/2*ln(z)*z
       - ln(z)*x*[1-z]^-1
       + 1/2*ln(z)*x
       + 1/2*ln(z)*x*z
       + 1/2*ln([1-z])
       + 1/2*ln([1-z])*z
       + 1/2*ln([1-z])*x
       + 1/2*ln([1-z])*x*z
       )

       + LMUF*NC * ( 5/2
       - 5/4*z^-1
       + 5/4*x*z^-1
       - 5/2*x
       + ln(x)
       - 1/2*ln(x)*z^-1
       - 1/2*ln(x)*x*z^-1
       + ln(x)*x
       )

       + LMUF*NC^2 * (  - 1/2
       - 3/4*z
       + 1/4*x
       - 1/2*x*z
       + 1/4*ln(x)
       - 1/2*ln(x)*[1-x]^-1
       - 1/2*ln(x)*[1-x]^-1*z
       + 3/4*ln(x)*z
       + 3/4*ln(x)*x
       + 1/4*ln(x)*x*z
       - 1/4*ln([1-x])
       - 5/4*ln([1-x])*z
       - 5/4*ln([1-x])*x
       - 1/4*ln([1-x])*x*z
       - 1/4*ln(z)
       + 1/2*ln(z)*[1-z]^-1
       - 1/4*ln(z)*z
       + 1/2*ln(z)*x*[1-z]^-1
       - 1/4*ln(z)*x
       - 1/4*ln(z)*x*z
       - 1/4*ln([1-z])
       - 1/4*ln([1-z])*z
       - 1/4*ln([1-z])*x
       - 1/4*ln([1-z])*x*z
       )

       + LMUA*NC^-2 * (  - 1/2
       + 1/4*z
       - 3/4*x
       - 1/2*x*z
       + 1/4*ln(x)
       - 1/2*ln(x)*[1-x]^-1
       - 1/2*ln(x)*[1-x]^-1*z
       + 1/4*ln(x)*z
       + 1/4*ln(x)*x
       + 1/4*ln(x)*x*z
       - 1/4*ln([1-x])
       - 1/4*ln([1-x])*z
       - 1/4*ln([1-x])*x
       - 1/4*ln([1-x])*x*z
       + 1/4*ln(z)
       - 1/2*ln(z)*[1-z]^-1
       + 3/4*ln(z)*z
       - 1/2*ln(z)*x*[1-z]^-1
       + 3/4*ln(z)*x
       + 1/4*ln(z)*x*z
       - 1/4*ln([1-z])
       - 5/4*ln([1-z])*z
       - 5/4*ln([1-z])*x
       - 1/4*ln([1-z])*x*z
       )

       + LMUA*NC^-1 * (  - 3/4
       - 1/6*z^-1
       + 1/4*z
       + 2/3*z^2
       - 1/6*x*z^-1
       - 1/4*x
       + 3/4*x*z
       - 1/3*x*z^2
       - 1/2*ln(z)
       - ln(z)*z
       - 1/2*ln(z)*x
       )

       + LMUA * ( 1
       - 1/2*z
       + 3/2*x
       + x*z
       - 1/2*ln(x)
       + ln(x)*[1-x]^-1
       + ln(x)*[1-x]^-1*z
       - 1/2*ln(x)*z
       - 1/2*ln(x)*x
       - 1/2*ln(x)*x*z
       + 1/2*ln([1-x])
       + 1/2*ln([1-x])*z
       + 1/2*ln([1-x])*x
       + 1/2*ln([1-x])*x*z
       - 1/2*ln(z)
       + ln(z)*[1-z]^-1
       - 3/2*ln(z)*z
       + ln(z)*x*[1-z]^-1
       - 3/2*ln(z)*x
       - 1/2*ln(z)*x*z
       + 1/2*ln([1-z])
       + 5/2*ln([1-z])*z
       + 5/2*ln([1-z])*x
       + 1/2*ln([1-z])*x*z
       )

       + LMUA*NC * ( 3/4
       + 1/6*z^-1
       - 1/4*z
       - 2/3*z^2
       + 1/6*x*z^-1
       + 1/4*x
       - 3/4*x*z
       + 1/3*x*z^2
       + 1/2*ln(z)
       + ln(z)*z
       + 1/2*ln(z)*x
       )

       + LMUA*NC^2 * (  - 1/2
       + 1/4*z
       - 3/4*x
       - 1/2*x*z
       + 1/4*ln(x)
       - 1/2*ln(x)*[1-x]^-1
       - 1/2*ln(x)*[1-x]^-1*z
       + 1/4*ln(x)*z
       + 1/4*ln(x)*x
       + 1/4*ln(x)*x*z
       - 1/4*ln([1-x])
       - 1/4*ln([1-x])*z
       - 1/4*ln([1-x])*x
       - 1/4*ln([1-x])*x*z
       + 1/4*ln(z)
       - 1/2*ln(z)*[1-z]^-1
       + 3/4*ln(z)*z
       - 1/2*ln(z)*x*[1-z]^-1
       + 3/4*ln(z)*x
       + 1/4*ln(z)*x*z
       - 1/4*ln([1-z])
       - 5/4*ln([1-z])*z
       - 5/4*ln([1-z])*x
       - 1/4*ln([1-z])*x*z
       )

       + LMUA*LMUF*NC^-2 * ( 1/4
       + 1/4*z
       + 1/4*x
       + 1/4*x*z
       )

       + LMUA*LMUF * (  - 1/2
       - 1/2*z
       - 1/2*x
       - 1/2*x*z
       )

       + LMUA*LMUF*NC^2 * ( 1/4
       + 1/4*z
       + 1/4*x
       + 1/4*x*z
       )

       + Dd([1-x])*NC^-2 * ( 1/8
       - 1/8*z
       - 5*zeta3*[1-z]^-1
       + 2*zeta3
       + 2*zeta3*z
       + 1/8*pi^2*[1-z]^-1
       - 7/24*pi^2
       - 1/12*pi^2*z
       + 15/8*ln(z)
       - 19/4*ln(z)*[1-z]^-1
       + 29/8*ln(z)*z
       - 1/4*ln(z)*pi^2*[1-z]^-1
       + 1/6*ln(z)*pi^2
       + 1/6*ln(z)*pi^2*z
       - 21/16*ln(z)^2
       + 9/16*ln(z)^2*[1-z]^-1
       - 7/16*ln(z)^2*z
       + 25/48*ln(z)^3
       - 5/6*ln(z)^3*[1-z]^-1
       + 25/48*ln(z)^3*z
       + 1/8*ln(z)^2*ln([1-z])
       - 1/4*ln(z)^2*ln([1-z])*[1-z]^-1
       + 1/8*ln(z)^2*ln([1-z])*z
       + 3/4*ln(z)*ln([1-z])
       - 3/4*ln(z)*ln([1-z])*z
       - 1/2*ln(z)*ln([1-z])^2
       + 1/4*ln(z)*ln([1-z])^2*[1-z]^-1
       - 1/2*ln(z)*ln([1-z])^2*z
       + ln(z)*Li2(z)
       - 2*ln(z)*Li2(z)*[1-z]^-1
       + ln(z)*Li2(z)*z
       + 5/8*ln([1-z])
       + 7/2*ln([1-z])*z
       + 1/12*ln([1-z])*pi^2*[1-z]^-1
       + 5/24*ln([1-z])*pi^2
       + 5/24*ln([1-z])*pi^2*z
       - 1/4*ln([1-z])^3
       - 1/4*ln([1-z])^3*z
       - 1/4*ln([1-z])*Li2(1 - z)
       - 1/4*ln([1-z])*Li2(1 - z)*z
       - 1/4*ln([1-z])*Li2(z)
       - 1/2*ln([1-z])*Li2(z)*[1-z]^-1
       - 1/4*ln([1-z])*Li2(z)*z
       + 1/2*Li3(1 - z)
       - 3/2*Li3(1 - z)*[1-z]^-1
       + 1/2*Li3(1 - z)*z
       - 3*Li3(z)
       + 5*Li3(z)*[1-z]^-1
       - 3*Li3(z)*z
       + 9/4*Li2(z)
       - 3/4*Li2(z)*[1-z]^-1
       )

       + Dd([1-x])*NC^-1 * ( 11/18
       + 107/216*z^-1
       + 2/9*z
       - 287/216*z^2
       - zeta3
       - zeta3*z
       + 1/12*pi^2
       + 1/8*pi^2*z
       + 5/2*ln(z)
       + 7/36*ln(z)*z^-1
       + 9/4*ln(z)*z
       + 31/36*ln(z)*z^2
       - 1/12*ln(z)*pi^2
       - 1/12*ln(z)*pi^2*z
       + 1/16*ln(z)^2
       - 1/6*ln(z)^2*z^-1
       + 5/16*ln(z)^2*z
       - 5/24*ln(z)^3
       - 5/24*ln(z)^3*z
       + 3/4*ln(z)*ln([1-z])^2
       + 3/4*ln(z)*ln([1-z])^2*z
       + 5/3*ln([1-z])
       + 7/36*ln([1-z])*z^-1
       - 7/6*ln([1-z])*z
       - 25/36*ln([1-z])*z^2
       - 1/6*ln([1-z])*pi^2
       - 1/6*ln([1-z])*pi^2*z
       - 1/8*ln([1-z])^2
       - 1/6*ln([1-z])^2*z^-1
       + 1/8*ln([1-z])^2*z
       + 1/6*ln([1-z])^2*z^2
       + 1/2*ln([1-z])*Li2(1 - z)
       + 1/2*ln([1-z])*Li2(1 - z)*z
       + ln([1-z])*Li2(z)
       + ln([1-z])*Li2(z)*z
       + 1/2*Li3(1 - z)
       + 1/2*Li3(1 - z)*z
       + Li3(z)
       + Li3(z)*z
       - 1/4*Li2(z)
       + 1/3*Li2(z)*z^-1
       - Li2(z)*z
       - 1/3*Li2(z)*z^2
       )

       + Dd([1-x])*NC^-1*NF * ( 19/108
       + 37/108*z
       - 1/36*pi^2
       - 1/36*pi^2*z
       + 1/36*ln(z)
       + 5/18*ln(z)*[1-z]^-1
       - 11/36*ln(z)*z
       - 1/24*ln(z)^2
       + 1/12*ln(z)^2*[1-z]^-1
       - 1/24*ln(z)^2*z
       - 1/9*ln([1-z])
       - 4/9*ln([1-z])*z
       + 1/12*ln([1-z])^2
       + 1/12*ln([1-z])^2*z
       )

       + Dd([1-x])*NC * (  - 11/18
       - 107/216*z^-1
       - 2/9*z
       + 287/216*z^2
       + zeta3
       + zeta3*z
       - 1/12*pi^2
       - 1/8*pi^2*z
       - 5/2*ln(z)
       - 7/36*ln(z)*z^-1
       - 9/4*ln(z)*z
       - 31/36*ln(z)*z^2
       + 1/12*ln(z)*pi^2
       + 1/12*ln(z)*pi^2*z
       - 1/16*ln(z)^2
       + 1/6*ln(z)^2*z^-1
       - 5/16*ln(z)^2*z
       + 5/24*ln(z)^3
       + 5/24*ln(z)^3*z
       - 3/4*ln(z)*ln([1-z])^2
       - 3/4*ln(z)*ln([1-z])^2*z
       - 5/3*ln([1-z])
       - 7/36*ln([1-z])*z^-1
       + 7/6*ln([1-z])*z
       + 25/36*ln([1-z])*z^2
       + 1/6*ln([1-z])*pi^2
       + 1/6*ln([1-z])*pi^2*z
       + 1/8*ln([1-z])^2
       + 1/6*ln([1-z])^2*z^-1
       - 1/8*ln([1-z])^2*z
       - 1/6*ln([1-z])^2*z^2
       - 1/2*ln([1-z])*Li2(1 - z)
       - 1/2*ln([1-z])*Li2(1 - z)*z
       - ln([1-z])*Li2(z)
       - ln([1-z])*Li2(z)*z
       - 1/2*Li3(1 - z)
       - 1/2*Li3(1 - z)*z
       - Li3(z)
       - Li3(z)*z
       + 1/4*Li2(z)
       - 1/3*Li2(z)*z^-1
       + Li2(z)*z
       + 1/3*Li2(z)*z^2
       )

       + Dd([1-x])*NC*NF * (  - 19/108
       - 37/108*z
       + 1/36*pi^2
       + 1/36*pi^2*z
       - 1/36*ln(z)
       - 5/18*ln(z)*[1-z]^-1
       + 11/36*ln(z)*z
       + 1/24*ln(z)^2
       - 1/12*ln(z)^2*[1-z]^-1
       + 1/24*ln(z)^2*z
       + 1/9*ln([1-z])
       + 4/9*ln([1-z])*z
       - 1/12*ln([1-z])^2
       - 1/12*ln([1-z])^2*z
       )

       + Dd([1-x])*NC^2 * ( 197/216
       + 611/216*z
       - 7*zeta3*[1-z]^-1
       + 5/4*zeta3
       + 5/4*zeta3*z
       + 1/8*pi^2*[1-z]^-1
       - 29/72*pi^2
       - 1/9*pi^2*z
       + 233/72*ln(z)
       - 113/36*ln(z)*[1-z]^-1
       + 101/72*ln(z)*z
       - 5/12*ln(z)*pi^2*[1-z]^-1
       + 1/4*ln(z)*pi^2
       + 1/4*ln(z)*pi^2*z
       - 25/24*ln(z)^2
       + 49/48*ln(z)^2*[1-z]^-1
       - 5/12*ln(z)^2*z
       + 5/16*ln(z)^3
       - 5/12*ln(z)^3*[1-z]^-1
       + 5/16*ln(z)^3*z
       + 1/8*ln(z)^2*ln([1-z])
       - 1/4*ln(z)^2*ln([1-z])*[1-z]^-1
       + 1/8*ln(z)^2*ln([1-z])*z
       + 3/4*ln(z)*ln([1-z])
       - 3/4*ln(z)*ln([1-z])*z
       - ln(z)*ln([1-z])^2
       + 5/4*ln(z)*ln([1-z])^2*[1-z]^-1
       - ln(z)*ln([1-z])^2*z
       + 3/2*ln(z)*Li2(z)
       - 3*ln(z)*Li2(z)*[1-z]^-1
       + 3/2*ln(z)*Li2(z)*z
       + 67/72*ln([1-z])
       - 7/9*ln([1-z])*z
       - 1/12*ln([1-z])*pi^2*[1-z]^-1
       + 3/8*ln([1-z])*pi^2
       + 3/8*ln([1-z])*pi^2*z
       + 11/24*ln([1-z])^2
       + 11/24*ln([1-z])^2*z
       - 1/4*ln([1-z])^3
       - 1/4*ln([1-z])^3*z
       - 1/4*ln([1-z])*Li2(1 - z)
       - 1/4*ln([1-z])*Li2(1 - z)*z
       - 3/4*ln([1-z])*Li2(z)
       + 1/2*ln([1-z])*Li2(z)*[1-z]^-1
       - 3/4*ln([1-z])*Li2(z)*z
       - Li3(1 - z)
       + 3/2*Li3(1 - z)*[1-z]^-1
       - Li3(1 - z)*z
       - 4*Li3(z)
       + 7*Li3(z)*[1-z]^-1
       - 4*Li3(z)*z
       + 7/4*Li2(z)
       - 3/4*Li2(z)*[1-z]^-1
       - 1/2*Li2(z)*z
       )

       + Dd([1-x])*LMUR*NC^-1*NF * ( 1/6
       - 1/6*z
       - 1/6*ln(z)
       + 1/3*ln(z)*[1-z]^-1
       - 1/6*ln(z)*z
       - 1/6*ln([1-z])
       - 1/6*ln([1-z])*z
       )

       + Dd([1-x])*LMUR * (  - 11/12
       + 11/12*z
       + 11/12*ln(z)
       - 11/6*ln(z)*[1-z]^-1
       + 11/12*ln(z)*z
       + 11/12*ln([1-z])
       + 11/12*ln([1-z])*z
       )

       + Dd([1-x])*LMUR*NC*NF * (  - 1/6
       + 1/6*z
       + 1/6*ln(z)
       - 1/3*ln(z)*[1-z]^-1
       + 1/6*ln(z)*z
       + 1/6*ln([1-z])
       + 1/6*ln([1-z])*z
       )

       + Dd([1-x])*LMUR*NC^2 * ( 11/12
       - 11/12*z
       - 11/12*ln(z)
       + 11/6*ln(z)*[1-z]^-1
       - 11/12*ln(z)*z
       - 11/12*ln([1-z])
       - 11/12*ln([1-z])*z
       )

       + Dd([1-x])*LMUF*NC^-2 * (  - 3/8
       + 3/8*z
       - 1/12*pi^2
       - 1/12*pi^2*z
       + 3/8*ln(z)
       - 3/4*ln(z)*[1-z]^-1
       + 3/8*ln(z)*z
       + 3/8*ln([1-z])
       + 3/8*ln([1-z])*z
       )

       + Dd([1-x])*LMUF * ( 3/4
       - 3/4*z
       + 1/6*pi^2
       + 1/6*pi^2*z
       - 3/4*ln(z)
       + 3/2*ln(z)*[1-z]^-1
       - 3/4*ln(z)*z
       - 3/4*ln([1-z])
       - 3/4*ln([1-z])*z
       )

       + Dd([1-x])*LMUF*NC^2 * (  - 3/8
       + 3/8*z
       - 1/12*pi^2
       - 1/12*pi^2*z
       + 3/8*ln(z)
       - 3/4*ln(z)*[1-z]^-1
       + 3/8*ln(z)*z
       + 3/8*ln([1-z])
       + 3/8*ln([1-z])*z
       )

       + Dd([1-x])*LMUA*NC^-2 * (  - 7/8
       - 25/8*z
       - 1/8*pi^2
       - 1/8*pi^2*z
       + 17/8*ln(z)
       - 3/2*ln(z)*[1-z]^-1
       + 7/8*ln(z)*z
       - ln(z)^2
       + 3/2*ln(z)^2*[1-z]^-1
       - ln(z)^2*z
       + 1/2*ln(z)*ln([1-z])
       - ln(z)*ln([1-z])*[1-z]^-1
       + 1/2*ln(z)*ln([1-z])*z
       + 3/8*ln([1-z])
       + 3/8*ln([1-z])*z
       + 3/4*ln([1-z])^2
       + 3/4*ln([1-z])^2*z
       + 1/4*Li2(z)
       + 1/4*Li2(z)*z
       )

       + Dd([1-x])*LMUA*NC^-1 * (  - 5/3
       - 7/36*z^-1
       + 7/6*z
       + 25/36*z^2
       + 1/12*pi^2
       + 1/12*pi^2*z
       - 1/4*ln(z)
       + 1/3*ln(z)*z^-1
       - ln(z)*z
       - 1/3*ln(z)*z^2
       + 1/2*ln(z)^2
       + 1/2*ln(z)^2*z
       + 1/4*ln([1-z])
       + 1/3*ln([1-z])*z^-1
       - 1/4*ln([1-z])*z
       - 1/3*ln([1-z])*z^2
       - 1/2*Li2(z)
       - 1/2*Li2(z)*z
       )

       + Dd([1-x])*LMUA*NC^-1*NF * (  - 1/18
       + 11/18*z
       + 1/6*ln(z)
       - 1/3*ln(z)*[1-z]^-1
       + 1/6*ln(z)*z
       )

       + Dd([1-x])*LMUA * ( 29/9
       + 19/18*z
       + 1/3*pi^2
       + 1/3*pi^2*z
       - 14/3*ln(z)
       + 29/6*ln(z)*[1-z]^-1
       - 13/6*ln(z)*z
       + 7/4*ln(z)^2
       - 5/2*ln(z)^2*[1-z]^-1
       + 7/4*ln(z)^2*z
       - ln(z)*ln([1-z])
       + 2*ln(z)*ln([1-z])*[1-z]^-1
       - ln(z)*ln([1-z])*z
       - 3/4*ln([1-z])
       - 3/4*ln([1-z])*z
       - 3/2*ln([1-z])^2
       - 3/2*ln([1-z])^2*z
       - 1/2*Li2(z)
       - 1/2*Li2(z)*z
       )

       + Dd([1-x])*LMUA*NC * ( 5/3
       + 7/36*z^-1
       - 7/6*z
       - 25/36*z^2
       - 1/12*pi^2
       - 1/12*pi^2*z
       + 1/4*ln(z)
       - 1/3*ln(z)*z^-1
       + ln(z)*z
       + 1/3*ln(z)*z^2
       - 1/2*ln(z)^2
       - 1/2*ln(z)^2*z
       - 1/4*ln([1-z])
       - 1/3*ln([1-z])*z^-1
       + 1/4*ln([1-z])*z
       + 1/3*ln([1-z])*z^2
       + 1/2*Li2(z)
       + 1/2*Li2(z)*z
       )

       + Dd([1-x])*LMUA*NC*NF * ( 1/18
       - 11/18*z
       - 1/6*ln(z)
       + 1/3*ln(z)*[1-z]^-1
       - 1/6*ln(z)*z
       )

       + Dd([1-x])*LMUA*NC^2 * (  - 169/72
       + 149/72*z
       - 5/24*pi^2
       - 5/24*pi^2*z
       + 61/24*ln(z)
       - 10/3*ln(z)*[1-z]^-1
       + 31/24*ln(z)*z
       - 3/4*ln(z)^2
       + ln(z)^2*[1-z]^-1
       - 3/4*ln(z)^2*z
       + 1/2*ln(z)*ln([1-z])
       - ln(z)*ln([1-z])*[1-z]^-1
       + 1/2*ln(z)*ln([1-z])*z
       + 3/8*ln([1-z])
       + 3/8*ln([1-z])*z
       + 3/4*ln([1-z])^2
       + 3/4*ln([1-z])^2*z
       + 1/4*Li2(z)
       + 1/4*Li2(z)*z
       )

       + Dd([1-x])*LMUA*LMUR*NC^-1*NF * ( 1/6
       + 1/6*z
       )

       + Dd([1-x])*LMUA*LMUR * (  - 11/12
       - 11/12*z
       )

       + Dd([1-x])*LMUA*LMUR*NC*NF * (  - 1/6
       - 1/6*z
       )

       + Dd([1-x])*LMUA*LMUR*NC^2 * ( 11/12
       + 11/12*z
       )

       + Dd([1-x])*LMUA*LMUF*NC^-2 * (  - 3/8
       - 3/8*z
       )

       + Dd([1-x])*LMUA*LMUF * ( 3/4
       + 3/4*z
       )

       + Dd([1-x])*LMUA*LMUF*NC^2 * (  - 3/8
       - 3/8*z
       )

       + Dd([1-x])*LMUA^2*NC^-2 * (  - 5/8
       - 1/8*z
       + 3/8*ln(z)
       - 1/2*ln(z)*[1-z]^-1
       + 3/8*ln(z)*z
       - 1/2*ln([1-z])
       - 1/2*ln([1-z])*z
       )

       + Dd([1-x])*LMUA^2*NC^-1 * (  - 1/8
       - 1/6*z^-1
       + 1/8*z
       + 1/6*z^2
       - 1/4*ln(z)
       - 1/4*ln(z)*z
       )

       + Dd([1-x])*LMUA^2*NC^-1*NF * (  - 1/12
       - 1/12*z
       )

       + Dd([1-x])*LMUA^2 * ( 41/24
       + 17/24*z
       - 3/4*ln(z)
       + ln(z)*[1-z]^-1
       - 3/4*ln(z)*z
       + ln([1-z])
       + ln([1-z])*z
       )

       + Dd([1-x])*LMUA^2*NC * ( 1/8
       + 1/6*z^-1
       - 1/8*z
       - 1/6*z^2
       + 1/4*ln(z)
       + 1/4*ln(z)*z
       )

       + Dd([1-x])*LMUA^2*NC*NF * ( 1/12
       + 1/12*z
       )

       + Dd([1-x])*LMUA^2*NC^2 * (  - 13/12
       - 7/12*z
       + 3/8*ln(z)
       - 1/2*ln(z)*[1-z]^-1
       + 3/8*ln(z)*z
       - 1/2*ln([1-z])
       - 1/2*ln([1-z])*z
       )

       + Dd([1-x])*Dd([1-z])*NC^-2 * ( 511/64
       - 15/4*zeta3
       + 29/48*pi^2
       - 7/360*pi^4
       )

       + Dd([1-x])*Dd([1-z])*NC^-1*NF * (  - 127/48
       - 1/3*zeta3
       - 19/108*pi^2
       )

       + Dd([1-x])*Dd([1-z])*NC*NF * ( 127/48
       + 1/3*zeta3
       + 19/108*pi^2
       )

       + Dd([1-x])*Dd([1-z])*NC^2 * (  - 1537/192
       + 41/12*zeta3
       - 277/432*pi^2
       + 1/36*pi^4
       )

       + Dd([1-x])*Dd([1-z])*LMUR*NC^-1*NF * (  - 4/3
       )

       + Dd([1-x])*Dd([1-z])*LMUR * ( 22/3
       )

       + Dd([1-x])*Dd([1-z])*LMUR*NC*NF * ( 4/3
       )

       + Dd([1-x])*Dd([1-z])*LMUR*NC^2 * (  - 22/3
       )

       + Dd([1-x])*Dd([1-z])*LMUF*NC^-2 * ( 93/32
       - 5/2*zeta3
       + 1/8*pi^2
       )

       + Dd([1-x])*Dd([1-z])*LMUF*NC^-1*NF * (  - 1/24
       - 1/18*pi^2
       )

       + Dd([1-x])*Dd([1-z])*LMUF * (  - 131/24
       + 7/2*zeta3
       + 1/18*pi^2
       )

       + Dd([1-x])*Dd([1-z])*LMUF*NC*NF * ( 1/24
       + 1/18*pi^2
       )

       + Dd([1-x])*Dd([1-z])*LMUF*NC^2 * ( 245/96
       - zeta3
       - 13/72*pi^2
       )

       + Dd([1-x])*Dd([1-z])*LMUF*LMUR*NC^-1*NF * (  - 1/4
       )

       + Dd([1-x])*Dd([1-z])*LMUF*LMUR * ( 11/8
       )

       + Dd([1-x])*Dd([1-z])*LMUF*LMUR*NC*NF * ( 1/4
       )

       + Dd([1-x])*Dd([1-z])*LMUF*LMUR*NC^2 * (  - 11/8
       )

       + Dd([1-x])*Dd([1-z])*LMUF^2*NC^-2 * ( 9/32
       - 1/12*pi^2
       )

       + Dd([1-x])*Dd([1-z])*LMUF^2*NC^-1*NF * ( 1/8
       )

       + Dd([1-x])*Dd([1-z])*LMUF^2 * (  - 5/4
       + 1/6*pi^2
       )

       + Dd([1-x])*Dd([1-z])*LMUF^2*NC*NF * (  - 1/8
       )

       + Dd([1-x])*Dd([1-z])*LMUF^2*NC^2 * ( 31/32
       - 1/12*pi^2
       )

       + Dd([1-x])*Dd([1-z])*LMUA*NC^-2 * ( 93/32
       - 5/2*zeta3
       + 1/8*pi^2
       )

       + Dd([1-x])*Dd([1-z])*LMUA*NC^-1*NF * (  - 1/24
       - 1/18*pi^2
       )

       + Dd([1-x])*Dd([1-z])*LMUA * (  - 131/24
       + 7/2*zeta3
       + 1/18*pi^2
       )

       + Dd([1-x])*Dd([1-z])*LMUA*NC*NF * ( 1/24
       + 1/18*pi^2
       )

       + Dd([1-x])*Dd([1-z])*LMUA*NC^2 * ( 245/96
       - zeta3
       - 13/72*pi^2
       )

       + Dd([1-x])*Dd([1-z])*LMUA*LMUR*NC^-1*NF * (  - 1/4
       )

       + Dd([1-x])*Dd([1-z])*LMUA*LMUR * ( 11/8
       )

       + Dd([1-x])*Dd([1-z])*LMUA*LMUR*NC*NF * ( 1/4
       )

       + Dd([1-x])*Dd([1-z])*LMUA*LMUR*NC^2 * (  - 11/8
       )

       + Dd([1-x])*Dd([1-z])*LMUA*LMUF*NC^-2 * ( 9/16
       )

       + Dd([1-x])*Dd([1-z])*LMUA*LMUF * (  - 9/8
       )

       + Dd([1-x])*Dd([1-z])*LMUA*LMUF*NC^2 * ( 9/16
       )

       + Dd([1-x])*Dd([1-z])*LMUA^2*NC^-2 * ( 9/32
       - 1/12*pi^2
       )

       + Dd([1-x])*Dd([1-z])*LMUA^2*NC^-1*NF * ( 1/8
       )

       + Dd([1-x])*Dd([1-z])*LMUA^2 * (  - 5/4
       + 1/6*pi^2
       )

       + Dd([1-x])*Dd([1-z])*LMUA^2*NC*NF * (  - 1/8
       )

       + Dd([1-x])*Dd([1-z])*LMUA^2*NC^2 * ( 31/32
       - 1/12*pi^2
       )

       + Dd([1-x])*Dd([1-z]) * ( 1/48
       + 1/3*zeta3
       + 1/27*pi^2
       - 1/120*pi^4
       )

       + Dd([1-x])*Dn(0,[1-z])*NC^-2 * ( 2*zeta3
       )

       + Dd([1-x])*Dn(0,[1-z])*NC^-1*NF * (  - 14/27
       + 1/18*pi^2
       )

       + Dd([1-x])*Dn(0,[1-z])*NC*NF * ( 14/27
       - 1/18*pi^2
       )

       + Dd([1-x])*Dn(0,[1-z])*NC^2 * (  - 101/27
       + 11/2*zeta3
       + 11/36*pi^2
       )

       + Dd([1-x])*Dn(0,[1-z])*LMUF*NC^-2 * ( 1/6*pi^2
       )

       + Dd([1-x])*Dn(0,[1-z])*LMUF * (  - 1/3*pi^2
       )

       + Dd([1-x])*Dn(0,[1-z])*LMUF*NC^2 * ( 1/6*pi^2
       )

       + Dd([1-x])*Dn(0,[1-z])*LMUA*NC^-2 * ( 4
       + 1/6*pi^2
       )

       + Dd([1-x])*Dn(0,[1-z])*LMUA*NC^-1*NF * (  - 5/9
       )

       + Dd([1-x])*Dn(0,[1-z])*LMUA * (  - 77/18
       - 1/2*pi^2
       )

       + Dd([1-x])*Dn(0,[1-z])*LMUA*NC*NF * ( 5/9
       )

       + Dd([1-x])*Dn(0,[1-z])*LMUA*NC^2 * ( 5/18
       + 1/3*pi^2
       )

       + Dd([1-x])*Dn(0,[1-z])*LMUA*LMUR*NC^-1*NF * (  - 1/3
       )

       + Dd([1-x])*Dn(0,[1-z])*LMUA*LMUR * ( 11/6
       )

       + Dd([1-x])*Dn(0,[1-z])*LMUA*LMUR*NC*NF * ( 1/3
       )

       + Dd([1-x])*Dn(0,[1-z])*LMUA*LMUR*NC^2 * (  - 11/6
       )

       + Dd([1-x])*Dn(0,[1-z])*LMUA*LMUF*NC^-2 * ( 3/4
       )

       + Dd([1-x])*Dn(0,[1-z])*LMUA*LMUF * (  - 3/2
       )

       + Dd([1-x])*Dn(0,[1-z])*LMUA*LMUF*NC^2 * ( 3/4
       )

       + Dd([1-x])*Dn(0,[1-z])*LMUA^2*NC^-2 * ( 3/4
       )

       + Dd([1-x])*Dn(0,[1-z])*LMUA^2*NC^-1*NF * ( 1/6
       )

       + Dd([1-x])*Dn(0,[1-z])*LMUA^2 * (  - 29/12
       )

       + Dd([1-x])*Dn(0,[1-z])*LMUA^2*NC*NF * (  - 1/6
       )

       + Dd([1-x])*Dn(0,[1-z])*LMUA^2*NC^2 * ( 5/3
       )

       + Dd([1-x])*Dn(0,[1-z]) * ( 101/27
       - 15/2*zeta3
       - 11/36*pi^2
       )

       + Dd([1-x])*Dn(1,[1-z])*NC^-2 * (  - 4
       - 1/3*pi^2
       )

       + Dd([1-x])*Dn(1,[1-z])*NC^-1*NF * ( 5/9
       )

       + Dd([1-x])*Dn(1,[1-z])*NC*NF * (  - 5/9
       )

       + Dd([1-x])*Dn(1,[1-z])*NC^2 * (  - 5/18
       - 1/2*pi^2
       )

       + Dd([1-x])*Dn(1,[1-z])*LMUR*NC^-1*NF * ( 1/3
       )

       + Dd([1-x])*Dn(1,[1-z])*LMUR * (  - 11/6
       )

       + Dd([1-x])*Dn(1,[1-z])*LMUR*NC*NF * (  - 1/3
       )

       + Dd([1-x])*Dn(1,[1-z])*LMUR*NC^2 * ( 11/6
       )

       + Dd([1-x])*Dn(1,[1-z])*LMUF*NC^-2 * (  - 3/4
       )

       + Dd([1-x])*Dn(1,[1-z])*LMUF * ( 3/2
       )

       + Dd([1-x])*Dn(1,[1-z])*LMUF*NC^2 * (  - 3/4
       )

       + Dd([1-x])*Dn(1,[1-z])*LMUA*NC^-2 * (  - 3/4
       )

       + Dd([1-x])*Dn(1,[1-z])*LMUA * ( 3/2
       )

       + Dd([1-x])*Dn(1,[1-z])*LMUA*NC^2 * (  - 3/4
       )

       + Dd([1-x])*Dn(1,[1-z])*LMUA^2*NC^-2 * ( 1
       )

       + Dd([1-x])*Dn(1,[1-z])*LMUA^2 * (  - 2
       )

       + Dd([1-x])*Dn(1,[1-z])*LMUA^2*NC^2 * ( 1
       )

       + Dd([1-x])*Dn(1,[1-z]) * ( 77/18
       + 5/6*pi^2
       )

       + Dd([1-x])*Dn(2,[1-z])*NC^-1*NF * (  - 1/6
       )

       + Dd([1-x])*Dn(2,[1-z])*NC*NF * ( 1/6
       )

       + Dd([1-x])*Dn(2,[1-z])*NC^2 * (  - 11/12
       )

       + Dd([1-x])*Dn(2,[1-z])*LMUA*NC^-2 * (  - 3/2
       )

       + Dd([1-x])*Dn(2,[1-z])*LMUA * ( 3
       )

       + Dd([1-x])*Dn(2,[1-z])*LMUA*NC^2 * (  - 3/2
       )

       + Dd([1-x])*Dn(2,[1-z]) * ( 11/12
       )

       + Dd([1-x])*Dn(3,[1-z])*NC^-2 * ( 1/2
       )

       + Dd([1-x])*Dn(3,[1-z])*NC^2 * ( 1/2
       )

       + Dd([1-x])*Dn(3,[1-z]) * (  - 1
       )

       + Dd([1-x]) * (  - 28/27
       - 73/27*z
       + 12*zeta3*[1-z]^-1
       - 13/4*zeta3
       - 13/4*zeta3*z
       - 1/4*pi^2*[1-z]^-1
       + 25/36*pi^2
       + 7/36*pi^2*z
       - 46/9*ln(z)
       + 71/9*ln(z)*[1-z]^-1
       - 181/36*ln(z)*z
       + 2/3*ln(z)*pi^2*[1-z]^-1
       - 5/12*ln(z)*pi^2
       - 5/12*ln(z)*pi^2*z
       + 113/48*ln(z)^2
       - 19/12*ln(z)^2*[1-z]^-1
       + 41/48*ln(z)^2*z
       - 5/6*ln(z)^3
       + 5/4*ln(z)^3*[1-z]^-1
       - 5/6*ln(z)^3*z
       - 1/4*ln(z)^2*ln([1-z])
       + 1/2*ln(z)^2*ln([1-z])*[1-z]^-1
       - 1/4*ln(z)^2*ln([1-z])*z
       - 3/2*ln(z)*ln([1-z])
       + 3/2*ln(z)*ln([1-z])*z
       + 3/2*ln(z)*ln([1-z])^2
       - 3/2*ln(z)*ln([1-z])^2*[1-z]^-1
       + 3/2*ln(z)*ln([1-z])^2*z
       - 5/2*ln(z)*Li2(z)
       + 5*ln(z)*Li2(z)*[1-z]^-1
       - 5/2*ln(z)*Li2(z)*z
       - 14/9*ln([1-z])
       - 49/18*ln([1-z])*z
       - 7/12*ln([1-z])*pi^2
       - 7/12*ln([1-z])*pi^2*z
       - 11/24*ln([1-z])^2
       - 11/24*ln([1-z])^2*z
       + 1/2*ln([1-z])^3
       + 1/2*ln([1-z])^3*z
       + 1/2*ln([1-z])*Li2(1 - z)
       + 1/2*ln([1-z])*Li2(1 - z)*z
       + ln([1-z])*Li2(z)
       + ln([1-z])*Li2(z)*z
       + 1/2*Li3(1 - z)
       + 1/2*Li3(1 - z)*z
       + 7*Li3(z)
       - 12*Li3(z)*[1-z]^-1
       + 7*Li3(z)*z
       - 4*Li2(z)
       + 3/2*Li2(z)*[1-z]^-1
       + 1/2*Li2(z)*z
       )

       + Dd([1-z])*NC^-2 * (  - 57/8
       + 59/8*x
       + 7/2*zeta3*[1-x]^-1
       - 11/4*zeta3
       - 11/4*zeta3*x
       - 1/8*pi^2*[1-x]^-1
       + 1/4*pi^2
       - 1/24*pi^2*x
       - 1/8*ln(x)
       + 5*ln(x)*[1-x]^-1
       - 4*ln(x)*[1+x]
       - 3*ln(x)*x
       + 1/2*ln(x)*pi^2*[1-x]^-1
       - 1/3*ln(x)*pi^2
       - 1/3*ln(x)*pi^2*x
       + 2*ln(x)*ln(1 + x)
       + 2*ln(x)*ln(1 + x)*x
       - 5/8*ln(x)^2
       + 9/16*ln(x)^2*[1-x]^-1
       + 1/4*ln(x)^2*x
       + 5/48*ln(x)^3
       + 5/48*ln(x)^3*x
       - 9/8*ln(x)^2*ln([1-x])
       + 9/4*ln(x)^2*ln([1-x])*[1-x]^-1
       - 9/8*ln(x)^2*ln([1-x])*x
       - 1/4*ln(x)*ln([1-x])
       + 1/4*ln(x)*ln([1-x])*x
       + 1/4*ln(x)*ln([1-x])^2
       - 5/4*ln(x)*ln([1-x])^2*[1-x]^-1
       + 1/4*ln(x)*ln([1-x])^2*x
       - ln(x)*Li2(x)
       + 5/2*ln(x)*Li2(x)*[1-x]^-1
       - ln(x)*Li2(x)*x
       + 1/2*ln([1-x])
       + 27/8*ln([1-x])*x
       - 1/12*ln([1-x])*pi^2*[1-x]^-1
       + 7/24*ln([1-x])*pi^2
       + 7/24*ln([1-x])*pi^2*x
       - 1/4*ln([1-x])^3
       - 1/4*ln([1-x])^3*x
       - 1/4*ln([1-x])*Li2(1 - x)
       - 1/4*ln([1-x])*Li2(1 - x)*x
       - 3/4*ln([1-x])*Li2(x)
       + 1/2*ln([1-x])*Li2(x)*[1-x]^-1
       - 3/4*ln([1-x])*Li2(x)*x
       - Li3(1 - x)
       + 3/2*Li3(1 - x)*[1-x]^-1
       - Li3(1 - x)*x
       + 7/4*Li3(x)
       - 7/2*Li3(x)*[1-x]^-1
       + 7/4*Li3(x)*x
       + 2*Li2( - x)
       + 2*Li2( - x)*x
       + 3/4*Li2(x)*[1-x]^-1
       + 3/4*Li2(x)*x
       )

       + Dd([1-z])*NC^-1 * (  - 9/2
       + 9/2*x
       + 1/4*pi^2
       - 1/4*pi^2*x
       - 25/8*ln(x)
       + 3/8*ln(x)*x
       + 1/6*ln(x)*pi^2
       + 1/6*ln(x)*pi^2*x
       - 17/16*ln(x)^2
       + 15/16*ln(x)^2*x
       - 5/24*ln(x)^3
       - 5/24*ln(x)^3*x
       + 5/4*ln(x)*ln([1-x])
       - 5/4*ln(x)*ln([1-x])*x
       + 3/4*ln(x)*ln([1-x])^2
       + 3/4*ln(x)*ln([1-x])^2*x
       - 1/2*ln(x)*Li2(x)
       - 1/2*ln(x)*Li2(x)*x
       + 3/2*ln([1-x])
       - 3/2*ln([1-x])*x
       - 1/6*ln([1-x])*pi^2
       - 1/6*ln([1-x])*pi^2*x
       - 5/8*ln([1-x])^2
       + 5/8*ln([1-x])^2*x
       + 1/2*ln([1-x])*Li2(1 - x)
       + 1/2*ln([1-x])*Li2(1 - x)*x
       + ln([1-x])*Li2(x)
       + ln([1-x])*Li2(x)*x
       + 1/2*Li3(1 - x)
       + 1/2*Li3(1 - x)*x
       - 1/4*Li2(x)
       + 1/4*Li2(x)*x
       )

       + Dd([1-z])*NC^-1*NF * ( 37/108
       + 19/108*x
       + 1/18*pi^2*[1-x]^-1
       - 1/18*pi^2
       - 1/18*pi^2*x
       + 5/12*ln(x)
       - 5/6*ln(x)*[1-x]^-1
       + 7/12*ln(x)*x
       + 5/24*ln(x)^2
       - 5/12*ln(x)^2*[1-x]^-1
       + 5/24*ln(x)^2*x
       - 1/6*ln(x)*ln([1-x])
       + 1/3*ln(x)*ln([1-x])*[1-x]^-1
       - 1/6*ln(x)*ln([1-x])*x
       - 1/9*ln([1-x])
       - 4/9*ln([1-x])*x
       + 1/12*ln([1-x])^2
       + 1/12*ln([1-x])^2*x
       + 1/6*Li2(x)
       - 1/3*Li2(x)*[1-x]^-1
       + 1/6*Li2(x)*x
       )

       + Dd([1-z])*NC * ( 9/2
       - 9/2*x
       - 1/4*pi^2
       + 1/4*pi^2*x
       + 25/8*ln(x)
       - 3/8*ln(x)*x
       - 1/6*ln(x)*pi^2
       - 1/6*ln(x)*pi^2*x
       + 17/16*ln(x)^2
       - 15/16*ln(x)^2*x
       + 5/24*ln(x)^3
       + 5/24*ln(x)^3*x
       - 5/4*ln(x)*ln([1-x])
       + 5/4*ln(x)*ln([1-x])*x
       - 3/4*ln(x)*ln([1-x])^2
       - 3/4*ln(x)*ln([1-x])^2*x
       + 1/2*ln(x)*Li2(x)
       + 1/2*ln(x)*Li2(x)*x
       - 3/2*ln([1-x])
       + 3/2*ln([1-x])*x
       + 1/6*ln([1-x])*pi^2
       + 1/6*ln([1-x])*pi^2*x
       + 5/8*ln([1-x])^2
       - 5/8*ln([1-x])^2*x
       - 1/2*ln([1-x])*Li2(1 - x)
       - 1/2*ln([1-x])*Li2(1 - x)*x
       - ln([1-x])*Li2(x)
       - ln([1-x])*Li2(x)*x
       - 1/2*Li3(1 - x)
       - 1/2*Li3(1 - x)*x
       + 1/4*Li2(x)
       - 1/4*Li2(x)*x
       )

       + Dd([1-z])*NC*NF * (  - 37/108
       - 19/108*x
       - 1/18*pi^2*[1-x]^-1
       + 1/18*pi^2
       + 1/18*pi^2*x
       - 5/12*ln(x)
       + 5/6*ln(x)*[1-x]^-1
       - 7/12*ln(x)*x
       - 5/24*ln(x)^2
       + 5/12*ln(x)^2*[1-x]^-1
       - 5/24*ln(x)^2*x
       + 1/6*ln(x)*ln([1-x])
       - 1/3*ln(x)*ln([1-x])*[1-x]^-1
       + 1/6*ln(x)*ln([1-x])*x
       + 1/9*ln([1-x])
       + 4/9*ln([1-x])*x
       - 1/12*ln([1-x])^2
       - 1/12*ln([1-x])^2*x
       - 1/6*Li2(x)
       + 1/3*Li2(x)*[1-x]^-1
       - 1/6*Li2(x)*x
       )

       + Dd([1-z])*NC^2 * (  - 91/216
       + 845/216*x
       + 5/2*zeta3*[1-x]^-1
       - 4*zeta3
       - 4*zeta3*x
       + 13/72*pi^2*[1-x]^-1
       - 13/72*pi^2
       - 7/18*pi^2*x
       + 97/24*ln(x)
       - 5/6*ln(x)*[1-x]^-1
       - 4*ln(x)*[1+x]
       + 61/12*ln(x)*x
       + 5/6*ln(x)*pi^2*[1-x]^-1
       - 1/2*ln(x)*pi^2
       - 1/2*ln(x)*pi^2*x
       + 2*ln(x)*ln(1 + x)
       - 2*ln(x)*ln(1 + x)*[1+x]
       + 2*ln(x)*ln(1 + x)*x
       + 49/48*ln(x)^2
       - 83/48*ln(x)^2*[1-x]^-1
       + 55/48*ln(x)^2*x
       + 11/48*ln(x)^3
       - 1/4*ln(x)^3*[1-x]^-1
       + 11/48*ln(x)^3*x
       - 5/8*ln(x)^2*ln([1-x])
       + 5/4*ln(x)^2*ln([1-x])*[1-x]^-1
       - 5/8*ln(x)^2*ln([1-x])*x
       - 7/6*ln(x)*ln([1-x])
       + 11/6*ln(x)*ln([1-x])*[1-x]^-1
       - 2/3*ln(x)*ln([1-x])*x
       + 3/4*ln(x)*ln([1-x])^2
       - 9/4*ln(x)*ln([1-x])^2*[1-x]^-1
       + 3/4*ln(x)*ln([1-x])^2*x
       + 1/2*ln(x)*Li2(x)*[1-x]^-1
       + 19/18*ln([1-x])
       - 47/72*ln([1-x])*x
       + 1/12*ln([1-x])*pi^2*[1-x]^-1
       + 7/24*ln([1-x])*pi^2
       + 7/24*ln([1-x])*pi^2*x
       + 11/24*ln([1-x])^2
       + 11/24*ln([1-x])^2*x
       - 1/4*ln([1-x])^3
       - 1/4*ln([1-x])^3*x
       - 1/4*ln([1-x])*Li2(1 - x)
       - 1/4*ln([1-x])*Li2(1 - x)*x
       - 1/4*ln([1-x])*Li2(x)
       - 1/2*ln([1-x])*Li2(x)*[1-x]^-1
       - 1/4*ln([1-x])*Li2(x)*x
       + Li3(1 - x)
       - 5/2*Li3(1 - x)*[1-x]^-1
       + Li3(1 - x)*x
       + 5/4*Li3(x)
       - 5/2*Li3(x)*[1-x]^-1
       + 5/4*Li3(x)*x
       + 2*Li2( - x)
       - 2*Li2( - x)*[1+x]
       + 2*Li2( - x)*x
       + 5/12*Li2(x)
       - 13/12*Li2(x)*[1-x]^-1
       + 7/6*Li2(x)*x
       )

       + Dd([1-z])*LMUR*NC^-1*NF * ( 1/6
       - 1/6*x
       + 1/6*ln(x)
       - 1/3*ln(x)*[1-x]^-1
       + 1/6*ln(x)*x
       - 1/6*ln([1-x])
       - 1/6*ln([1-x])*x
       )

       + Dd([1-z])*LMUR * (  - 11/12
       + 11/12*x
       - 11/12*ln(x)
       + 11/6*ln(x)*[1-x]^-1
       - 11/12*ln(x)*x
       + 11/12*ln([1-x])
       + 11/12*ln([1-x])*x
       )

       + Dd([1-z])*LMUR*NC*NF * (  - 1/6
       + 1/6*x
       - 1/6*ln(x)
       + 1/3*ln(x)*[1-x]^-1
       - 1/6*ln(x)*x
       + 1/6*ln([1-x])
       + 1/6*ln([1-x])*x
       )

       + Dd([1-z])*LMUR*NC^2 * ( 11/12
       - 11/12*x
       + 11/12*ln(x)
       - 11/6*ln(x)*[1-x]^-1
       + 11/12*ln(x)*x
       - 11/12*ln([1-x])
       - 11/12*ln([1-x])*x
       )

       + Dd([1-z])*LMUF*NC^-2 * (  - 7/8
       - 25/8*x
       - 1/8*pi^2
       - 1/8*pi^2*x
       - 3/8*ln(x)
       + 3/2*ln(x)*[1-x]^-1
       + 3/8*ln(x)*x
       + 1/2*ln(x)^2
       - 1/2*ln(x)^2*[1-x]^-1
       + 1/2*ln(x)^2*x
       - 3/2*ln(x)*ln([1-x])
       + 3*ln(x)*ln([1-x])*[1-x]^-1
       - 3/2*ln(x)*ln([1-x])*x
       + 3/8*ln([1-x])
       + 3/8*ln([1-x])*x
       + 3/4*ln([1-x])^2
       + 3/4*ln([1-x])^2*x
       + 1/4*Li2(x)
       + 1/4*Li2(x)*x
       )

       + Dd([1-z])*LMUF*NC^-1 * (  - 3/2
       + 3/2*x
       + 1/12*pi^2
       + 1/12*pi^2*x
       - 3/2*ln(x)
       + 3/2*ln(x)*x
       - 1/2*ln(x)^2
       - 1/2*ln(x)^2*x
       + 5/4*ln([1-x])
       - 5/4*ln([1-x])*x
       - 1/2*Li2(x)
       - 1/2*Li2(x)*x
       )

       + Dd([1-z])*LMUF*NC^-1*NF * (  - 1/18
       + 11/18*x
       + 1/6*ln(x)
       - 1/3*ln(x)*[1-x]^-1
       + 1/6*ln(x)*x
       )

       + Dd([1-z])*LMUF * ( 29/9
       + 19/18*x
       + 1/3*pi^2
       + 1/3*pi^2*x
       + 1/3*ln(x)
       - 7/6*ln(x)*[1-x]^-1
       - 7/6*ln(x)*x
       - 5/4*ln(x)^2
       + 3/2*ln(x)^2*[1-x]^-1
       - 5/4*ln(x)^2*x
       + 3*ln(x)*ln([1-x])
       - 6*ln(x)*ln([1-x])*[1-x]^-1
       + 3*ln(x)*ln([1-x])*x
       - 3/4*ln([1-x])
       - 3/4*ln([1-x])*x
       - 3/2*ln([1-x])^2
       - 3/2*ln([1-x])^2*x
       - 1/2*Li2(x)
       - 1/2*Li2(x)*x
       )

       + Dd([1-z])*LMUF*NC * ( 3/2
       - 3/2*x
       - 1/12*pi^2
       - 1/12*pi^2*x
       + 3/2*ln(x)
       - 3/2*ln(x)*x
       + 1/2*ln(x)^2
       + 1/2*ln(x)^2*x
       - 5/4*ln([1-x])
       + 5/4*ln([1-x])*x
       + 1/2*Li2(x)
       + 1/2*Li2(x)*x
       )

       + Dd([1-z])*LMUF*NC*NF * ( 1/18
       - 11/18*x
       - 1/6*ln(x)
       + 1/3*ln(x)*[1-x]^-1
       - 1/6*ln(x)*x
       )

       + Dd([1-z])*LMUF*NC^2 * (  - 169/72
       + 149/72*x
       - 5/24*pi^2
       - 5/24*pi^2*x
       + 1/24*ln(x)
       - 1/3*ln(x)*[1-x]^-1
       + 19/24*ln(x)*x
       + 3/4*ln(x)^2
       - ln(x)^2*[1-x]^-1
       + 3/4*ln(x)^2*x
       - 3/2*ln(x)*ln([1-x])
       + 3*ln(x)*ln([1-x])*[1-x]^-1
       - 3/2*ln(x)*ln([1-x])*x
       + 3/8*ln([1-x])
       + 3/8*ln([1-x])*x
       + 3/4*ln([1-x])^2
       + 3/4*ln([1-x])^2*x
       + 1/4*Li2(x)
       + 1/4*Li2(x)*x
       )

       + Dd([1-z])*LMUF*LMUR*NC^-1*NF * ( 1/6
       + 1/6*x
       )

       + Dd([1-z])*LMUF*LMUR * (  - 11/12
       - 11/12*x
       )

       + Dd([1-z])*LMUF*LMUR*NC*NF * (  - 1/6
       - 1/6*x
       )

       + Dd([1-z])*LMUF*LMUR*NC^2 * ( 11/12
       + 11/12*x
       )

       + Dd([1-z])*LMUF^2*NC^-2 * (  - 5/8
       - 1/8*x
       + 3/8*ln(x)
       - 1/2*ln(x)*[1-x]^-1
       + 3/8*ln(x)*x
       - 1/2*ln([1-x])
       - 1/2*ln([1-x])*x
       )

       + Dd([1-z])*LMUF^2*NC^-1 * (  - 5/8
       + 5/8*x
       - 1/4*ln(x)
       - 1/4*ln(x)*x
       )

       + Dd([1-z])*LMUF^2*NC^-1*NF * (  - 1/12
       - 1/12*x
       )

       + Dd([1-z])*LMUF^2 * ( 41/24
       + 17/24*x
       - 3/4*ln(x)
       + ln(x)*[1-x]^-1
       - 3/4*ln(x)*x
       + ln([1-x])
       + ln([1-x])*x
       )

       + Dd([1-z])*LMUF^2*NC * ( 5/8
       - 5/8*x
       + 1/4*ln(x)
       + 1/4*ln(x)*x
       )

       + Dd([1-z])*LMUF^2*NC*NF * ( 1/12
       + 1/12*x
       )

       + Dd([1-z])*LMUF^2*NC^2 * (  - 13/12
       - 7/12*x
       + 3/8*ln(x)
       - 1/2*ln(x)*[1-x]^-1
       + 3/8*ln(x)*x
       - 1/2*ln([1-x])
       - 1/2*ln([1-x])*x
       )

       + Dd([1-z])*LMUA*NC^-2 * (  - 3/8
       + 3/8*x
       - 1/12*pi^2
       - 1/12*pi^2*x
       - 3/8*ln(x)
       + 3/4*ln(x)*[1-x]^-1
       - 3/8*ln(x)*x
       + 3/8*ln([1-x])
       + 3/8*ln([1-x])*x
       )

       + Dd([1-z])*LMUA * ( 3/4
       - 3/4*x
       + 1/6*pi^2
       + 1/6*pi^2*x
       + 3/4*ln(x)
       - 3/2*ln(x)*[1-x]^-1
       + 3/4*ln(x)*x
       - 3/4*ln([1-x])
       - 3/4*ln([1-x])*x
       )

       + Dd([1-z])*LMUA*NC^2 * (  - 3/8
       + 3/8*x
       - 1/12*pi^2
       - 1/12*pi^2*x
       - 3/8*ln(x)
       + 3/4*ln(x)*[1-x]^-1
       - 3/8*ln(x)*x
       + 3/8*ln([1-x])
       + 3/8*ln([1-x])*x
       )

       + Dd([1-z])*LMUA*LMUF*NC^-2 * (  - 3/8
       - 3/8*x
       )

       + Dd([1-z])*LMUA*LMUF * ( 3/4
       + 3/4*x
       )

       + Dd([1-z])*LMUA*LMUF*NC^2 * (  - 3/8
       - 3/8*x
       )

       + Dd([1-z])*Dn(0,[1-x])*NC^-2 * ( 2*zeta3
       )

       + Dd([1-z])*Dn(0,[1-x])*NC^-1*NF * (  - 14/27
       + 1/18*pi^2
       )

       + Dd([1-z])*Dn(0,[1-x])*NC*NF * ( 14/27
       - 1/18*pi^2
       )

       + Dd([1-z])*Dn(0,[1-x])*NC^2 * (  - 101/27
       + 11/2*zeta3
       + 11/36*pi^2
       )

       + Dd([1-z])*Dn(0,[1-x])*LMUF*NC^-2 * ( 4
       + 1/6*pi^2
       )

       + Dd([1-z])*Dn(0,[1-x])*LMUF*NC^-1*NF * (  - 5/9
       )

       + Dd([1-z])*Dn(0,[1-x])*LMUF * (  - 77/18
       - 1/2*pi^2
       )

       + Dd([1-z])*Dn(0,[1-x])*LMUF*NC*NF * ( 5/9
       )

       + Dd([1-z])*Dn(0,[1-x])*LMUF*NC^2 * ( 5/18
       + 1/3*pi^2
       )

       + Dd([1-z])*Dn(0,[1-x])*LMUF*LMUR*NC^-1*NF * (  - 1/3
       )

       + Dd([1-z])*Dn(0,[1-x])*LMUF*LMUR * ( 11/6
       )

       + Dd([1-z])*Dn(0,[1-x])*LMUF*LMUR*NC*NF * ( 1/3
       )

       + Dd([1-z])*Dn(0,[1-x])*LMUF*LMUR*NC^2 * (  - 11/6
       )

       + Dd([1-z])*Dn(0,[1-x])*LMUF^2*NC^-2 * ( 3/4
       )

       + Dd([1-z])*Dn(0,[1-x])*LMUF^2*NC^-1*NF * ( 1/6
       )

       + Dd([1-z])*Dn(0,[1-x])*LMUF^2 * (  - 29/12
       )

       + Dd([1-z])*Dn(0,[1-x])*LMUF^2*NC*NF * (  - 1/6
       )

       + Dd([1-z])*Dn(0,[1-x])*LMUF^2*NC^2 * ( 5/3
       )

       + Dd([1-z])*Dn(0,[1-x])*LMUA*NC^-2 * ( 1/6*pi^2
       )

       + Dd([1-z])*Dn(0,[1-x])*LMUA * (  - 1/3*pi^2
       )

       + Dd([1-z])*Dn(0,[1-x])*LMUA*NC^2 * ( 1/6*pi^2
       )

       + Dd([1-z])*Dn(0,[1-x])*LMUA*LMUF*NC^-2 * ( 3/4
       )

       + Dd([1-z])*Dn(0,[1-x])*LMUA*LMUF * (  - 3/2
       )

       + Dd([1-z])*Dn(0,[1-x])*LMUA*LMUF*NC^2 * ( 3/4
       )

       + Dd([1-z])*Dn(0,[1-x]) * ( 101/27
       - 15/2*zeta3
       - 11/36*pi^2
       )

       + Dd([1-z])*Dn(1,[1-x])*NC^-2 * (  - 4
       - 1/3*pi^2
       )

       + Dd([1-z])*Dn(1,[1-x])*NC^-1*NF * ( 5/9
       )

       + Dd([1-z])*Dn(1,[1-x])*NC*NF * (  - 5/9
       )

       + Dd([1-z])*Dn(1,[1-x])*NC^2 * (  - 5/18
       - 1/2*pi^2
       )

       + Dd([1-z])*Dn(1,[1-x])*LMUR*NC^-1*NF * ( 1/3
       )

       + Dd([1-z])*Dn(1,[1-x])*LMUR * (  - 11/6
       )

       + Dd([1-z])*Dn(1,[1-x])*LMUR*NC*NF * (  - 1/3
       )

       + Dd([1-z])*Dn(1,[1-x])*LMUR*NC^2 * ( 11/6
       )

       + Dd([1-z])*Dn(1,[1-x])*LMUF*NC^-2 * (  - 3/4
       )

       + Dd([1-z])*Dn(1,[1-x])*LMUF * ( 3/2
       )

       + Dd([1-z])*Dn(1,[1-x])*LMUF*NC^2 * (  - 3/4
       )

       + Dd([1-z])*Dn(1,[1-x])*LMUF^2*NC^-2 * ( 1
       )

       + Dd([1-z])*Dn(1,[1-x])*LMUF^2 * (  - 2
       )

       + Dd([1-z])*Dn(1,[1-x])*LMUF^2*NC^2 * ( 1
       )

       + Dd([1-z])*Dn(1,[1-x])*LMUA*NC^-2 * (  - 3/4
       )

       + Dd([1-z])*Dn(1,[1-x])*LMUA * ( 3/2
       )

       + Dd([1-z])*Dn(1,[1-x])*LMUA*NC^2 * (  - 3/4
       )

       + Dd([1-z])*Dn(1,[1-x]) * ( 77/18
       + 5/6*pi^2
       )

       + Dd([1-z])*Dn(2,[1-x])*NC^-1*NF * (  - 1/6
       )

       + Dd([1-z])*Dn(2,[1-x])*NC*NF * ( 1/6
       )

       + Dd([1-z])*Dn(2,[1-x])*NC^2 * (  - 11/12
       )

       + Dd([1-z])*Dn(2,[1-x])*LMUF*NC^-2 * (  - 3/2
       )

       + Dd([1-z])*Dn(2,[1-x])*LMUF * ( 3
       )

       + Dd([1-z])*Dn(2,[1-x])*LMUF*NC^2 * (  - 3/2
       )

       + Dd([1-z])*Dn(2,[1-x]) * ( 11/12
       )

       + Dd([1-z])*Dn(3,[1-x])*NC^-2 * ( 1/2
       )

       + Dd([1-z])*Dn(3,[1-x])*NC^2 * ( 1/2
       )

       + Dd([1-z])*Dn(3,[1-x]) * (  - 1
       )

       + Dd([1-z]) * ( 815/108
       - 1219/108*x
       - 6*zeta3*[1-x]^-1
       + 27/4*zeta3
       + 27/4*zeta3*x
       - 1/18*pi^2*[1-x]^-1
       - 5/72*pi^2
       + 31/72*pi^2*x
       - 47/12*ln(x)
       - 25/6*ln(x)*[1-x]^-1
       + 8*ln(x)*[1+x]
       - 25/12*ln(x)*x
       - 4/3*ln(x)*pi^2*[1-x]^-1
       + 5/6*ln(x)*pi^2
       + 5/6*ln(x)*pi^2*x
       - 4*ln(x)*ln(1 + x)
       + 2*ln(x)*ln(1 + x)*[1+x]
       - 4*ln(x)*ln(1 + x)*x
       - 19/48*ln(x)^2
       + 7/6*ln(x)^2*[1-x]^-1
       - 67/48*ln(x)^2*x
       - 1/3*ln(x)^3
       + 1/4*ln(x)^3*[1-x]^-1
       - 1/3*ln(x)^3*x
       + 7/4*ln(x)^2*ln([1-x])
       - 7/2*ln(x)^2*ln([1-x])*[1-x]^-1
       + 7/4*ln(x)^2*ln([1-x])*x
       + 17/12*ln(x)*ln([1-x])
       - 11/6*ln(x)*ln([1-x])*[1-x]^-1
       + 5/12*ln(x)*ln([1-x])*x
       - ln(x)*ln([1-x])^2
       + 7/2*ln(x)*ln([1-x])^2*[1-x]^-1
       - ln(x)*ln([1-x])^2*x
       + ln(x)*Li2(x)
       - 3*ln(x)*Li2(x)*[1-x]^-1
       + ln(x)*Li2(x)*x
       - 14/9*ln([1-x])
       - 49/18*ln([1-x])*x
       - 7/12*ln([1-x])*pi^2
       - 7/12*ln([1-x])*pi^2*x
       - 11/24*ln([1-x])^2
       - 11/24*ln([1-x])^2*x
       + 1/2*ln([1-x])^3
       + 1/2*ln([1-x])^3*x
       + 1/2*ln([1-x])*Li2(1 - x)
       + 1/2*ln([1-x])*Li2(1 - x)*x
       + ln([1-x])*Li2(x)
       + ln([1-x])*Li2(x)*x
       + Li3(1 - x)*[1-x]^-1
       - 3*Li3(x)
       + 6*Li3(x)*[1-x]^-1
       - 3*Li3(x)*x
       - 4*Li2( - x)
       + 2*Li2( - x)*[1+x]
       - 4*Li2( - x)*x
       - 5/12*Li2(x)
       + 1/3*Li2(x)*[1-x]^-1
       - 23/12*Li2(x)*x
       )

       + Dn(0,[1-x])*NC^-2 * ( 1/2
       + 7/2*z
       + 5/24*pi^2
       + 5/24*pi^2*z
       - 7/4*ln(z)
       + 3/4*ln(z)*[1-z]^-1
       - 1/2*ln(z)*z
       + ln(z)^2
       - 3/2*ln(z)^2*[1-z]^-1
       + ln(z)^2*z
       - 1/2*ln(z)*ln([1-z])
       + ln(z)*ln([1-z])*[1-z]^-1
       - 1/2*ln(z)*ln([1-z])*z
       - 3/4*ln([1-z])^2
       - 3/4*ln([1-z])^2*z
       - 1/4*Li2(z)
       - 1/4*Li2(z)*z
       )

       + Dn(0,[1-x])*NC^-1 * ( 5/3
       + 7/36*z^-1
       - 7/6*z
       - 25/36*z^2
       - 1/12*pi^2
       - 1/12*pi^2*z
       + 1/4*ln(z)
       - 1/3*ln(z)*z^-1
       + ln(z)*z
       + 1/3*ln(z)*z^2
       - 1/2*ln(z)^2
       - 1/2*ln(z)^2*z
       - 1/4*ln([1-z])
       - 1/3*ln([1-z])*z^-1
       + 1/4*ln([1-z])*z
       + 1/3*ln([1-z])*z^2
       + 1/2*Li2(z)
       + 1/2*Li2(z)*z
       )

       + Dn(0,[1-x])*NC^-1*NF * (  - 1/9
       - 4/9*z
       + 1/6*ln([1-z])
       + 1/6*ln([1-z])*z
       )

       + Dn(0,[1-x])*NC * (  - 5/3
       - 7/36*z^-1
       + 7/6*z
       + 25/36*z^2
       + 1/12*pi^2
       + 1/12*pi^2*z
       - 1/4*ln(z)
       + 1/3*ln(z)*z^-1
       - ln(z)*z
       - 1/3*ln(z)*z^2
       + 1/2*ln(z)^2
       + 1/2*ln(z)^2*z
       + 1/4*ln([1-z])
       + 1/3*ln([1-z])*z^-1
       - 1/4*ln([1-z])*z
       - 1/3*ln([1-z])*z^2
       - 1/2*Li2(z)
       - 1/2*Li2(z)*z
       )

       + Dn(0,[1-x])*NC*NF * ( 1/9
       + 4/9*z
       - 1/6*ln([1-z])
       - 1/6*ln([1-z])*z
       )

       + Dn(0,[1-x])*NC^2 * ( 19/18
       - 7/9*z
       + 7/24*pi^2
       + 7/24*pi^2*z
       - 5/4*ln(z)
       + 3/4*ln(z)*[1-z]^-1
       + 3/4*ln(z)^2
       - ln(z)^2*[1-z]^-1
       + 3/4*ln(z)^2*z
       - 1/2*ln(z)*ln([1-z])
       + ln(z)*ln([1-z])*[1-z]^-1
       - 1/2*ln(z)*ln([1-z])*z
       + 11/12*ln([1-z])
       + 11/12*ln([1-z])*z
       - 3/4*ln([1-z])^2
       - 3/4*ln([1-z])^2*z
       - 1/4*Li2(z)
       - 1/4*Li2(z)*z
       )

       + Dn(0,[1-x])*LMUR*NC^-1*NF * (  - 1/6
       - 1/6*z
       )

       + Dn(0,[1-x])*LMUR * ( 11/12
       + 11/12*z
       )

       + Dn(0,[1-x])*LMUR*NC*NF * ( 1/6
       + 1/6*z
       )

       + Dn(0,[1-x])*LMUR*NC^2 * (  - 11/12
       - 11/12*z
       )

       + Dn(0,[1-x])*LMUF*NC^-2 * (  - 1/8
       + 7/8*z
       + 1/2*ln(z)
       - ln(z)*[1-z]^-1
       + 1/2*ln(z)*z
       + 1/2*ln([1-z])
       + 1/2*ln([1-z])*z
       )

       + Dn(0,[1-x])*LMUF * ( 1/4
       - 7/4*z
       - ln(z)
       + 2*ln(z)*[1-z]^-1
       - ln(z)*z
       - ln([1-z])
       - ln([1-z])*z
       )

       + Dn(0,[1-x])*LMUF*NC^2 * (  - 1/8
       + 7/8*z
       + 1/2*ln(z)
       - ln(z)*[1-z]^-1
       + 1/2*ln(z)*z
       + 1/2*ln([1-z])
       + 1/2*ln([1-z])*z
       )

       + Dn(0,[1-x])*LMUA*NC^-2 * ( 7/8
       - 1/8*z
       - 3/4*ln(z)
       + ln(z)*[1-z]^-1
       - 3/4*ln(z)*z
       + ln([1-z])
       + ln([1-z])*z
       )

       + Dn(0,[1-x])*LMUA*NC^-1 * ( 1/4
       + 1/3*z^-1
       - 1/4*z
       - 1/3*z^2
       + 1/2*ln(z)
       + 1/2*ln(z)*z
       )

       + Dn(0,[1-x])*LMUA * (  - 7/4
       + 1/4*z
       + 3/2*ln(z)
       - 2*ln(z)*[1-z]^-1
       + 3/2*ln(z)*z
       - 2*ln([1-z])
       - 2*ln([1-z])*z
       )

       + Dn(0,[1-x])*LMUA*NC * (  - 1/4
       - 1/3*z^-1
       + 1/4*z
       + 1/3*z^2
       - 1/2*ln(z)
       - 1/2*ln(z)*z
       )

       + Dn(0,[1-x])*LMUA*NC^2 * ( 7/8
       - 1/8*z
       - 3/4*ln(z)
       + ln(z)*[1-z]^-1
       - 3/4*ln(z)*z
       + ln([1-z])
       + ln([1-z])*z
       )

       + Dn(0,[1-x])*LMUA*LMUF*NC^-2 * (  - 1/2
       - 1/2*z
       )

       + Dn(0,[1-x])*LMUA*LMUF * ( 1
       + z
       )

       + Dn(0,[1-x])*LMUA*LMUF*NC^2 * (  - 1/2
       - 1/2*z
       )

       + Dn(0,[1-x])*Dn(0,[1-z])*NC^-2 * (  - 4
       - 1/3*pi^2
       )

       + Dn(0,[1-x])*Dn(0,[1-z])*NC^-1*NF * ( 5/9
       )

       + Dn(0,[1-x])*Dn(0,[1-z])*NC*NF * (  - 5/9
       )

       + Dn(0,[1-x])*Dn(0,[1-z])*NC^2 * (  - 5/18
       - 1/2*pi^2
       )

       + Dn(0,[1-x])*Dn(0,[1-z])*LMUR*NC^-1*NF * ( 1/3
       )

       + Dn(0,[1-x])*Dn(0,[1-z])*LMUR * (  - 11/6
       )

       + Dn(0,[1-x])*Dn(0,[1-z])*LMUR*NC*NF * (  - 1/3
       )

       + Dn(0,[1-x])*Dn(0,[1-z])*LMUR*NC^2 * ( 11/6
       )

       + Dn(0,[1-x])*Dn(0,[1-z])*LMUF*NC^-2 * (  - 3/4
       )

       + Dn(0,[1-x])*Dn(0,[1-z])*LMUF * ( 3/2
       )

       + Dn(0,[1-x])*Dn(0,[1-z])*LMUF*NC^2 * (  - 3/4
       )

       + Dn(0,[1-x])*Dn(0,[1-z])*LMUA*NC^-2 * (  - 3/4
       )

       + Dn(0,[1-x])*Dn(0,[1-z])*LMUA * ( 3/2
       )

       + Dn(0,[1-x])*Dn(0,[1-z])*LMUA*NC^2 * (  - 3/4
       )

       + Dn(0,[1-x])*Dn(0,[1-z])*LMUA*LMUF*NC^-2 * ( 1
       )

       + Dn(0,[1-x])*Dn(0,[1-z])*LMUA*LMUF * (  - 2
       )

       + Dn(0,[1-x])*Dn(0,[1-z])*LMUA*LMUF*NC^2 * ( 1
       )

       + Dn(0,[1-x])*Dn(0,[1-z]) * ( 77/18
       + 5/6*pi^2
       )

       + Dn(0,[1-x])*Dn(1,[1-z])*NC^-1*NF * (  - 1/3
       )

       + Dn(0,[1-x])*Dn(1,[1-z])*NC*NF * ( 1/3
       )

       + Dn(0,[1-x])*Dn(1,[1-z])*NC^2 * (  - 11/6
       )

       + Dn(0,[1-x])*Dn(1,[1-z])*LMUF*NC^-2 * (  - 1
       )

       + Dn(0,[1-x])*Dn(1,[1-z])*LMUF * ( 2
       )

       + Dn(0,[1-x])*Dn(1,[1-z])*LMUF*NC^2 * (  - 1
       )

       + Dn(0,[1-x])*Dn(1,[1-z])*LMUA*NC^-2 * (  - 2
       )

       + Dn(0,[1-x])*Dn(1,[1-z])*LMUA * ( 4
       )

       + Dn(0,[1-x])*Dn(1,[1-z])*LMUA*NC^2 * (  - 2
       )

       + Dn(0,[1-x])*Dn(1,[1-z]) * ( 11/6
       )

       + Dn(0,[1-x])*Dn(2,[1-z])*NC^-2 * ( 3/2
       )

       + Dn(0,[1-x])*Dn(2,[1-z])*NC^2 * ( 3/2
       )

       + Dn(0,[1-x])*Dn(2,[1-z]) * (  - 3
       )

       + Dn(0,[1-x]) * (  - 14/9
       - 49/18*z
       - 1/2*pi^2
       - 1/2*pi^2*z
       + 3*ln(z)
       - 3/2*ln(z)*[1-z]^-1
       + 1/2*ln(z)*z
       - 7/4*ln(z)^2
       + 5/2*ln(z)^2*[1-z]^-1
       - 7/4*ln(z)^2*z
       + ln(z)*ln([1-z])
       - 2*ln(z)*ln([1-z])*[1-z]^-1
       + ln(z)*ln([1-z])*z
       - 11/12*ln([1-z])
       - 11/12*ln([1-z])*z
       + 3/2*ln([1-z])^2
       + 3/2*ln([1-z])^2*z
       + 1/2*Li2(z)
       + 1/2*Li2(z)*z
       )

       + Dn(0,[1-z])*NC^-2 * ( 1/2
       + 7/2*x
       + 5/24*pi^2
       + 5/24*pi^2*x
       - 3/4*ln(x)*[1-x]^-1
       - 3/4*ln(x)*x
       - 1/2*ln(x)^2
       + 1/2*ln(x)^2*[1-x]^-1
       - 1/2*ln(x)^2*x
       + 3/2*ln(x)*ln([1-x])
       - 3*ln(x)*ln([1-x])*[1-x]^-1
       + 3/2*ln(x)*ln([1-x])*x
       - 3/4*ln([1-x])^2
       - 3/4*ln([1-x])^2*x
       - 1/4*Li2(x)
       - 1/4*Li2(x)*x
       )

       + Dn(0,[1-z])*NC^-1 * ( 3/2
       - 3/2*x
       - 1/12*pi^2
       - 1/12*pi^2*x
       + 3/2*ln(x)
       - 3/2*ln(x)*x
       + 1/2*ln(x)^2
       + 1/2*ln(x)^2*x
       - 5/4*ln([1-x])
       + 5/4*ln([1-x])*x
       + 1/2*Li2(x)
       + 1/2*Li2(x)*x
       )

       + Dn(0,[1-z])*NC^-1*NF * (  - 1/9
       - 4/9*x
       - 1/3*ln(x)
       + 2/3*ln(x)*[1-x]^-1
       - 1/3*ln(x)*x
       + 1/6*ln([1-x])
       + 1/6*ln([1-x])*x
       )

       + Dn(0,[1-z])*NC * (  - 3/2
       + 3/2*x
       + 1/12*pi^2
       + 1/12*pi^2*x
       - 3/2*ln(x)
       + 3/2*ln(x)*x
       - 1/2*ln(x)^2
       - 1/2*ln(x)^2*x
       + 5/4*ln([1-x])
       - 5/4*ln([1-x])*x
       - 1/2*Li2(x)
       - 1/2*Li2(x)*x
       )

       + Dn(0,[1-z])*NC*NF * ( 1/9
       + 4/9*x
       + 1/3*ln(x)
       - 2/3*ln(x)*[1-x]^-1
       + 1/3*ln(x)*x
       - 1/6*ln([1-x])
       - 1/6*ln([1-x])*x
       )

       + Dn(0,[1-z])*NC^2 * ( 19/18
       - 7/9*x
       + 7/24*pi^2
       + 7/24*pi^2*x
       - 4/3*ln(x)
       + 35/12*ln(x)*[1-x]^-1
       - 25/12*ln(x)*x
       - 3/4*ln(x)^2
       + ln(x)^2*[1-x]^-1
       - 3/4*ln(x)^2*x
       + 3/2*ln(x)*ln([1-x])
       - 3*ln(x)*ln([1-x])*[1-x]^-1
       + 3/2*ln(x)*ln([1-x])*x
       + 11/12*ln([1-x])
       + 11/12*ln([1-x])*x
       - 3/4*ln([1-x])^2
       - 3/4*ln([1-x])^2*x
       - 1/4*Li2(x)
       - 1/4*Li2(x)*x
       )

       + Dn(0,[1-z])*LMUR*NC^-1*NF * (  - 1/6
       - 1/6*x
       )

       + Dn(0,[1-z])*LMUR * ( 11/12
       + 11/12*x
       )

       + Dn(0,[1-z])*LMUR*NC*NF * ( 1/6
       + 1/6*x
       )

       + Dn(0,[1-z])*LMUR*NC^2 * (  - 11/12
       - 11/12*x
       )

       + Dn(0,[1-z])*LMUF*NC^-2 * ( 7/8
       - 1/8*x
       - 3/4*ln(x)
       + ln(x)*[1-x]^-1
       - 3/4*ln(x)*x
       + ln([1-x])
       + ln([1-x])*x
       )

       + Dn(0,[1-z])*LMUF*NC^-1 * ( 5/4
       - 5/4*x
       + 1/2*ln(x)
       + 1/2*ln(x)*x
       )

       + Dn(0,[1-z])*LMUF * (  - 7/4
       + 1/4*x
       + 3/2*ln(x)
       - 2*ln(x)*[1-x]^-1
       + 3/2*ln(x)*x
       - 2*ln([1-x])
       - 2*ln([1-x])*x
       )

       + Dn(0,[1-z])*LMUF*NC * (  - 5/4
       + 5/4*x
       - 1/2*ln(x)
       - 1/2*ln(x)*x
       )

       + Dn(0,[1-z])*LMUF*NC^2 * ( 7/8
       - 1/8*x
       - 3/4*ln(x)
       + ln(x)*[1-x]^-1
       - 3/4*ln(x)*x
       + ln([1-x])
       + ln([1-x])*x
       )

       + Dn(0,[1-z])*LMUA*NC^-2 * (  - 1/8
       + 7/8*x
       - 1/2*ln(x)
       + ln(x)*[1-x]^-1
       - 1/2*ln(x)*x
       + 1/2*ln([1-x])
       + 1/2*ln([1-x])*x
       )

       + Dn(0,[1-z])*LMUA * ( 1/4
       - 7/4*x
       + ln(x)
       - 2*ln(x)*[1-x]^-1
       + ln(x)*x
       - ln([1-x])
       - ln([1-x])*x
       )

       + Dn(0,[1-z])*LMUA*NC^2 * (  - 1/8
       + 7/8*x
       - 1/2*ln(x)
       + ln(x)*[1-x]^-1
       - 1/2*ln(x)*x
       + 1/2*ln([1-x])
       + 1/2*ln([1-x])*x
       )

       + Dn(0,[1-z])*LMUA*LMUF*NC^-2 * (  - 1/2
       - 1/2*x
       )

       + Dn(0,[1-z])*LMUA*LMUF * ( 1
       + x
       )

       + Dn(0,[1-z])*LMUA*LMUF*NC^2 * (  - 1/2
       - 1/2*x
       )

       + Dn(0,[1-z])*Dn(1,[1-x])*NC^-1*NF * (  - 1/3
       )

       + Dn(0,[1-z])*Dn(1,[1-x])*NC*NF * ( 1/3
       )

       + Dn(0,[1-z])*Dn(1,[1-x])*NC^2 * (  - 11/6
       )

       + Dn(0,[1-z])*Dn(1,[1-x])*LMUF*NC^-2 * (  - 2
       )

       + Dn(0,[1-z])*Dn(1,[1-x])*LMUF * ( 4
       )

       + Dn(0,[1-z])*Dn(1,[1-x])*LMUF*NC^2 * (  - 2
       )

       + Dn(0,[1-z])*Dn(1,[1-x])*LMUA*NC^-2 * (  - 1
       )

       + Dn(0,[1-z])*Dn(1,[1-x])*LMUA * ( 2
       )

       + Dn(0,[1-z])*Dn(1,[1-x])*LMUA*NC^2 * (  - 1
       )

       + Dn(0,[1-z])*Dn(1,[1-x]) * ( 11/6
       )

       + Dn(0,[1-z])*Dn(2,[1-x])*NC^-2 * ( 3/2
       )

       + Dn(0,[1-z])*Dn(2,[1-x])*NC^2 * ( 3/2
       )

       + Dn(0,[1-z])*Dn(2,[1-x]) * (  - 3
       )

       + Dn(0,[1-z]) * (  - 14/9
       - 49/18*x
       - 1/2*pi^2
       - 1/2*pi^2*x
       + 4/3*ln(x)
       - 13/6*ln(x)*[1-x]^-1
       + 17/6*ln(x)*x
       + 5/4*ln(x)^2
       - 3/2*ln(x)^2*[1-x]^-1
       + 5/4*ln(x)^2*x
       - 3*ln(x)*ln([1-x])
       + 6*ln(x)*ln([1-x])*[1-x]^-1
       - 3*ln(x)*ln([1-x])*x
       - 11/12*ln([1-x])
       - 11/12*ln([1-x])*x
       + 3/2*ln([1-x])^2
       + 3/2*ln([1-x])^2*x
       + 1/2*Li2(x)
       + 1/2*Li2(x)*x
       )

       + Dn(1,[1-x])*NC^-2 * ( 1/4*ln(z)
       + 1/4*ln(z)*z
       - 3/2*ln([1-z])
       - 3/2*ln([1-z])*z
       )

       + Dn(1,[1-x])*NC^-1 * (  - 1/4
       - 1/3*z^-1
       + 1/4*z
       + 1/3*z^2
       - 1/2*ln(z)
       - 1/2*ln(z)*z
       )

       + Dn(1,[1-x])*NC^-1*NF * ( 1/6
       + 1/6*z
       )

       + Dn(1,[1-x])*NC * ( 1/4
       + 1/3*z^-1
       - 1/4*z
       - 1/3*z^2
       + 1/2*ln(z)
       + 1/2*ln(z)*z
       )

       + Dn(1,[1-x])*NC*NF * (  - 1/6
       - 1/6*z
       )

       + Dn(1,[1-x])*NC^2 * ( 11/12
       + 11/12*z
       + 1/4*ln(z)
       + 1/4*ln(z)*z
       - 3/2*ln([1-z])
       - 3/2*ln([1-z])*z
       )

       + Dn(1,[1-x])*LMUF*NC^-2 * ( 1
       + z
       )

       + Dn(1,[1-x])*LMUF * (  - 2
       - 2*z
       )

       + Dn(1,[1-x])*LMUF*NC^2 * ( 1
       + z
       )

       + Dn(1,[1-x])*LMUA*NC^-2 * ( 1/2
       + 1/2*z
       )

       + Dn(1,[1-x])*LMUA * (  - 1
       - z
       )

       + Dn(1,[1-x])*LMUA*NC^2 * ( 1/2
       + 1/2*z
       )

       + Dn(1,[1-x])*Dn(1,[1-z])*NC^-2 * ( 3
       )

       + Dn(1,[1-x])*Dn(1,[1-z])*NC^2 * ( 3
       )

       + Dn(1,[1-x])*Dn(1,[1-z]) * (  - 6
       )

       + Dn(1,[1-x]) * (  - 11/12
       - 11/12*z
       - 1/2*ln(z)
       - 1/2*ln(z)*z
       + 3*ln([1-z])
       + 3*ln([1-z])*z
       )

       + Dn(1,[1-z])*NC^-2 * ( 5/4*ln(x)
       - 2*ln(x)*[1-x]^-1
       + 5/4*ln(x)*x
       - 3/2*ln([1-x])
       - 3/2*ln([1-x])*x
       )

       + Dn(1,[1-z])*NC^-1 * (  - 5/4
       + 5/4*x
       - 1/2*ln(x)
       - 1/2*ln(x)*x
       )

       + Dn(1,[1-z])*NC^-1*NF * ( 1/6
       + 1/6*x
       )

       + Dn(1,[1-z])*NC * ( 5/4
       - 5/4*x
       + 1/2*ln(x)
       + 1/2*ln(x)*x
       )

       + Dn(1,[1-z])*NC*NF * (  - 1/6
       - 1/6*x
       )

       + Dn(1,[1-z])*NC^2 * ( 11/12
       + 11/12*x
       + 5/4*ln(x)
       - 2*ln(x)*[1-x]^-1
       + 5/4*ln(x)*x
       - 3/2*ln([1-x])
       - 3/2*ln([1-x])*x
       )

       + Dn(1,[1-z])*LMUF*NC^-2 * ( 1/2
       + 1/2*x
       )

       + Dn(1,[1-z])*LMUF * (  - 1
       - x
       )

       + Dn(1,[1-z])*LMUF*NC^2 * ( 1/2
       + 1/2*x
       )

       + Dn(1,[1-z])*LMUA*NC^-2 * ( 1
       + x
       )

       + Dn(1,[1-z])*LMUA * (  - 2
       - 2*x
       )

       + Dn(1,[1-z])*LMUA*NC^2 * ( 1
       + x
       )

       + Dn(1,[1-z]) * (  - 11/12
       - 11/12*x
       - 5/2*ln(x)
       + 4*ln(x)*[1-x]^-1
       - 5/2*ln(x)*x
       + 3*ln([1-x])
       + 3*ln([1-x])*x
       )

       + Dn(2,[1-x])*NC^-2 * (  - 3/4
       - 3/4*z
       )

       + Dn(2,[1-x])*NC^2 * (  - 3/4
       - 3/4*z
       )

       + Dn(2,[1-x]) * ( 3/2
       + 3/2*z
       )

       + Dn(2,[1-z])*NC^-2 * (  - 3/4
       - 3/4*x
       )

       + Dn(2,[1-z])*NC^2 * (  - 3/4
       - 3/4*x
       )

       + Dn(2,[1-z]) * ( 3/2
       + 3/2*x
       )

       + T(u1)*NC^-2 * (  - 3
       - 4*[1-x]^-1*[1-z]^-1
       - 2*[1-x]^-1
       + 6*[1-x]^-1*z
       + 9*[1-z]^-1
       - 7*z
       - 5*x*[1-z]^-1
       + 5*x
       + 2/3*pi^2*[1-x]^-1*[1-z]^-1
       - 5/12*pi^2*[1-x]^-1
       - 5/12*pi^2*[1-x]^-1*z
       - 7/12*pi^2*[1-z]^-1
       + 1/2*pi^2
       + 2/3*pi^2*z
       - 1/12*pi^2*x*[1-z]^-1
       + 1/6*pi^2*x
       + 9*ln(x)
       + 18*ln(x)*[1-x]^-1*[1-z]^-1
       - 9*ln(x)*[1-x]^-1
       - 15*ln(x)*[1-x]^-1*z
       - 18*ln(x)*[1-z]^-1
       + 18*ln(x)*z
       + 6*ln(x)*x
       - 3*ln(x)^2
       - 6*ln(x)^2*[1-x]^-1*[1-z]^-1
       + 9/2*ln(x)^2*[1-x]^-1
       + 9/2*ln(x)^2*[1-x]^-1*z
       + 9/2*ln(x)^2*[1-z]^-1
       - 6*ln(x)^2*z
       + 3/2*ln(x)^2*x*[1-z]^-1
       - 3*ln(x)^2*x
       + 5*ln(x)*ln([1-x])
       + 9*ln(x)*ln([1-x])*[1-x]^-1*[1-z]^-1
       - 13/2*ln(x)*ln([1-x])*[1-x]^-1
       - 13/2*ln(x)*ln([1-x])*[1-x]^-1*z
       - 7*ln(x)*ln([1-x])*[1-z]^-1
       + 9*ln(x)*ln([1-x])*z
       - 2*ln(x)*ln([1-x])*x*[1-z]^-1
       + 4*ln(x)*ln([1-x])*x
       + 5*ln(x)*ln(z)
       + 10*ln(x)*ln(z)*[1-x]^-1*[1-z]^-1
       - 15/2*ln(x)*ln(z)*[1-x]^-1
       - 15/2*ln(x)*ln(z)*[1-x]^-1*z
       - 15/2*ln(x)*ln(z)*[1-z]^-1
       + 10*ln(x)*ln(z)*z
       - 5/2*ln(x)*ln(z)*x*[1-z]^-1
       + 5*ln(x)*ln(z)*x
       + 4*ln(x)*ln([1-z])
       + 9*ln(x)*ln([1-z])*[1-x]^-1*[1-z]^-1
       - 7*ln(x)*ln([1-z])*[1-x]^-1
       - 7*ln(x)*ln([1-z])*[1-x]^-1*z
       - 13/2*ln(x)*ln([1-z])*[1-z]^-1
       + 9*ln(x)*ln([1-z])*z
       - 5/2*ln(x)*ln([1-z])*x*[1-z]^-1
       + 5*ln(x)*ln([1-z])*x
       - ln(x)*ln([1-x-z])*[1-x]^-1*[1-z]^-1
       + ln(x)*ln([1-x-z])*[1-x]^-1
       + ln(x)*ln([1-x-z])*[1-x]^-1*z
       + 1/2*ln(x)*ln([1-x-z])*[1-z]^-1
       - ln(x)*ln([1-x-z])*z
       + 1/2*ln(x)*ln([1-x-z])*x*[1-z]^-1
       - ln(x)*ln([1-x-z])*x
       - ln(x)*ln([x-z])
       - 2*ln(x)*ln([x-z])*[1-x]^-1*[1-z]^-1
       + 3/2*ln(x)*ln([x-z])*[1-x]^-1
       + 3/2*ln(x)*ln([x-z])*[1-x]^-1*z
       + 3/2*ln(x)*ln([x-z])*[1-z]^-1
       - 2*ln(x)*ln([x-z])*z
       + 1/2*ln(x)*ln([x-z])*x*[1-z]^-1
       - ln(x)*ln([x-z])*x
       - 4*ln([1-x])
       - 9*ln([1-x])*[1-x]^-1*[1-z]^-1
       + 4*ln([1-x])*[1-x]^-1
       + 7*ln([1-x])*[1-x]^-1*z
       + 9*ln([1-x])*[1-z]^-1
       - 9*ln([1-x])*z
       - 3*ln([1-x])*x
       - 2*ln([1-x])^2
       - 5/2*ln([1-x])^2*[1-x]^-1*[1-z]^-1
       + 3/2*ln([1-x])^2*[1-x]^-1
       + 3/2*ln([1-x])^2*[1-x]^-1*z
       + 9/4*ln([1-x])^2*[1-z]^-1
       - 5/2*ln([1-x])^2*z
       + 1/4*ln([1-x])^2*x*[1-z]^-1
       - 1/2*ln([1-x])^2*x
       - 4*ln([1-x])*ln(z)
       - 6*ln([1-x])*ln(z)*[1-x]^-1*[1-z]^-1
       + 4*ln([1-x])*ln(z)*[1-x]^-1
       + 4*ln([1-x])*ln(z)*[1-x]^-1*z
       + 5*ln([1-x])*ln(z)*[1-z]^-1
       - 6*ln([1-x])*ln(z)*z
       + ln([1-x])*ln(z)*x*[1-z]^-1
       - 2*ln([1-x])*ln(z)*x
       - 2*ln([1-x])*ln([1-z])
       - 5*ln([1-x])*ln([1-z])*[1-x]^-1*[1-z]^-1
       + 4*ln([1-x])*ln([1-z])*[1-x]^-1
       + 4*ln([1-x])*ln([1-z])*[1-x]^-1*z
       + 7/2*ln([1-x])*ln([1-z])*[1-z]^-1
       - 5*ln([1-x])*ln([1-z])*z
       + 3/2*ln([1-x])*ln([1-z])*x*[1-z]^-1
       - 3*ln([1-x])*ln([1-z])*x
       - 6*ln(z)
       - 12*ln(z)*[1-x]^-1*[1-z]^-1
       + 6*ln(z)*[1-x]^-1
       + 10*ln(z)*[1-x]^-1*z
       + 12*ln(z)*[1-z]^-1
       - 12*ln(z)*z
       - 4*ln(z)*x
       - 2*ln(z)^2
       - 4*ln(z)^2*[1-x]^-1*[1-z]^-1
       + 3*ln(z)^2*[1-x]^-1
       + 3*ln(z)^2*[1-x]^-1*z
       + 3*ln(z)^2*[1-z]^-1
       - 4*ln(z)^2*z
       + ln(z)^2*x*[1-z]^-1
       - 2*ln(z)^2*x
       - 2*ln(z)*ln([1-z])
       - 6*ln(z)*ln([1-z])*[1-x]^-1*[1-z]^-1
       + 5*ln(z)*ln([1-z])*[1-x]^-1
       + 5*ln(z)*ln([1-z])*[1-x]^-1*z
       + 4*ln(z)*ln([1-z])*[1-z]^-1
       - 6*ln(z)*ln([1-z])*z
       + 2*ln(z)*ln([1-z])*x*[1-z]^-1
       - 4*ln(z)*ln([1-z])*x
       + ln(z)*ln([x-z])
       + 2*ln(z)*ln([x-z])*[1-x]^-1*[1-z]^-1
       - 3/2*ln(z)*ln([x-z])*[1-x]^-1
       - 3/2*ln(z)*ln([x-z])*[1-x]^-1*z
       - 3/2*ln(z)*ln([x-z])*[1-z]^-1
       + 2*ln(z)*ln([x-z])*z
       - 1/2*ln(z)*ln([x-z])*x*[1-z]^-1
       + ln(z)*ln([x-z])*x
       - 5*ln([1-z])
       - 9*ln([1-z])*[1-x]^-1*[1-z]^-1
       + 5*ln([1-z])*[1-x]^-1
       + 8*ln([1-z])*[1-x]^-1*z
       + 9*ln([1-z])*[1-z]^-1
       - 9*ln([1-z])*z
       - 3*ln([1-z])*x
       - 1/2*ln([1-z])^2
       - 2*ln([1-z])^2*[1-x]^-1*[1-z]^-1
       + 7/4*ln([1-z])^2*[1-x]^-1
       + 7/4*ln([1-z])^2*[1-x]^-1*z
       + 5/4*ln([1-z])^2*[1-z]^-1
       - 2*ln([1-z])^2*z
       + 3/4*ln([1-z])^2*x*[1-z]^-1
       - 3/2*ln([1-z])^2*x
       + ln([1-z])*ln([1-x-z])*[1-x]^-1*[1-z]^-1
       - ln([1-z])*ln([1-x-z])*[1-x]^-1
       - ln([1-z])*ln([1-x-z])*[1-x]^-1*z
       - 1/2*ln([1-z])*ln([1-x-z])*[1-z]^-1
       + ln([1-z])*ln([1-x-z])*z
       - 1/2*ln([1-z])*ln([1-x-z])*x*[1-z]^-1
       + ln([1-z])*ln([1-x-z])*x
       + Li2(x^-1*[1-x]^-1*z*[1-z])*[1-x]^-1*[1-z]^-1
       - Li2(x^-1*[1-x]^-1*z*[1-z])*[1-x]^-1
       - Li2(x^-1*[1-x]^-1*z*[1-z])*[1-x]^-1*z
       - 1/2*Li2(x^-1*[1-x]^-1*z*[1-z])*[1-z]^-1
       + Li2(x^-1*[1-x]^-1*z*[1-z])*z
       - 1/2*Li2(x^-1*[1-x]^-1*z*[1-z])*x*[1-z]^-1
       + Li2(x^-1*[1-x]^-1*z*[1-z])*x
       + Li2(x^-1*[1-x]*z*[1-z]^-1)
       + Li2(x^-1*[1-x]*z*[1-z]^-1)*[1-x]^-1*[1-z]^-1
       - 1/2*Li2(x^-1*[1-x]*z*[1-z]^-1)*[1-x]^-1
       - 1/2*Li2(x^-1*[1-x]*z*[1-z]^-1)*[1-x]^-1*z
       - Li2(x^-1*[1-x]*z*[1-z]^-1)*[1-z]^-1
       + Li2(x^-1*[1-x]*z*[1-z]^-1)*z
       - Li2([1-x]^-1*z)*[1-x]^-1*[1-z]^-1
       + Li2([1-x]^-1*z)*[1-x]^-1
       + Li2([1-x]^-1*z)*[1-x]^-1*z
       + 1/2*Li2([1-x]^-1*z)*[1-z]^-1
       - Li2([1-x]^-1*z)*z
       + 1/2*Li2([1-x]^-1*z)*x*[1-z]^-1
       - Li2([1-x]^-1*z)*x
       - Li2([1-x]*[1-z]^-1)
       - 1/2*Li2([1-x]*[1-z]^-1)*[1-x]^-1
       - 1/2*Li2([1-x]*[1-z]^-1)*[1-x]^-1*z
       + 1/2*Li2([1-x]*[1-z]^-1)*[1-z]^-1
       - 1/2*Li2([1-x]*[1-z]^-1)*x*[1-z]^-1
       + Li2([1-x]*[1-z]^-1)*x
       + Li2(z)
       + Li2(z)*[1-x]^-1*[1-z]^-1
       - 1/2*Li2(z)*[1-x]^-1
       - 1/2*Li2(z)*[1-x]^-1*z
       - Li2(z)*[1-z]^-1
       + Li2(z)*z
       )

       + T(u1) * ( 3
       + 4*[1-x]^-1*[1-z]^-1
       + 2*[1-x]^-1
       - 6*[1-x]^-1*z
       - 9*[1-z]^-1
       + 7*z
       + 5*x*[1-z]^-1
       - 5*x
       - 2/3*pi^2*[1-x]^-1*[1-z]^-1
       + 5/12*pi^2*[1-x]^-1
       + 5/12*pi^2*[1-x]^-1*z
       + 7/12*pi^2*[1-z]^-1
       - 1/2*pi^2
       - 2/3*pi^2*z
       + 1/12*pi^2*x*[1-z]^-1
       - 1/6*pi^2*x
       - 9*ln(x)
       - 18*ln(x)*[1-x]^-1*[1-z]^-1
       + 9*ln(x)*[1-x]^-1
       + 15*ln(x)*[1-x]^-1*z
       + 18*ln(x)*[1-z]^-1
       - 18*ln(x)*z
       - 6*ln(x)*x
       + 3*ln(x)^2
       + 6*ln(x)^2*[1-x]^-1*[1-z]^-1
       - 9/2*ln(x)^2*[1-x]^-1
       - 9/2*ln(x)^2*[1-x]^-1*z
       - 9/2*ln(x)^2*[1-z]^-1
       + 6*ln(x)^2*z
       - 3/2*ln(x)^2*x*[1-z]^-1
       + 3*ln(x)^2*x
       - 5*ln(x)*ln([1-x])
       - 9*ln(x)*ln([1-x])*[1-x]^-1*[1-z]^-1
       + 13/2*ln(x)*ln([1-x])*[1-x]^-1
       + 13/2*ln(x)*ln([1-x])*[1-x]^-1*z
       + 7*ln(x)*ln([1-x])*[1-z]^-1
       - 9*ln(x)*ln([1-x])*z
       + 2*ln(x)*ln([1-x])*x*[1-z]^-1
       - 4*ln(x)*ln([1-x])*x
       - 5*ln(x)*ln(z)
       - 10*ln(x)*ln(z)*[1-x]^-1*[1-z]^-1
       + 15/2*ln(x)*ln(z)*[1-x]^-1
       + 15/2*ln(x)*ln(z)*[1-x]^-1*z
       + 15/2*ln(x)*ln(z)*[1-z]^-1
       - 10*ln(x)*ln(z)*z
       + 5/2*ln(x)*ln(z)*x*[1-z]^-1
       - 5*ln(x)*ln(z)*x
       - 4*ln(x)*ln([1-z])
       - 9*ln(x)*ln([1-z])*[1-x]^-1*[1-z]^-1
       + 7*ln(x)*ln([1-z])*[1-x]^-1
       + 7*ln(x)*ln([1-z])*[1-x]^-1*z
       + 13/2*ln(x)*ln([1-z])*[1-z]^-1
       - 9*ln(x)*ln([1-z])*z
       + 5/2*ln(x)*ln([1-z])*x*[1-z]^-1
       - 5*ln(x)*ln([1-z])*x
       + ln(x)*ln([1-x-z])*[1-x]^-1*[1-z]^-1
       - ln(x)*ln([1-x-z])*[1-x]^-1
       - ln(x)*ln([1-x-z])*[1-x]^-1*z
       - 1/2*ln(x)*ln([1-x-z])*[1-z]^-1
       + ln(x)*ln([1-x-z])*z
       - 1/2*ln(x)*ln([1-x-z])*x*[1-z]^-1
       + ln(x)*ln([1-x-z])*x
       + ln(x)*ln([x-z])
       + 2*ln(x)*ln([x-z])*[1-x]^-1*[1-z]^-1
       - 3/2*ln(x)*ln([x-z])*[1-x]^-1
       - 3/2*ln(x)*ln([x-z])*[1-x]^-1*z
       - 3/2*ln(x)*ln([x-z])*[1-z]^-1
       + 2*ln(x)*ln([x-z])*z
       - 1/2*ln(x)*ln([x-z])*x*[1-z]^-1
       + ln(x)*ln([x-z])*x
       + 4*ln([1-x])
       + 9*ln([1-x])*[1-x]^-1*[1-z]^-1
       - 4*ln([1-x])*[1-x]^-1
       - 7*ln([1-x])*[1-x]^-1*z
       - 9*ln([1-x])*[1-z]^-1
       + 9*ln([1-x])*z
       + 3*ln([1-x])*x
       + 2*ln([1-x])^2
       + 5/2*ln([1-x])^2*[1-x]^-1*[1-z]^-1
       - 3/2*ln([1-x])^2*[1-x]^-1
       - 3/2*ln([1-x])^2*[1-x]^-1*z
       - 9/4*ln([1-x])^2*[1-z]^-1
       + 5/2*ln([1-x])^2*z
       - 1/4*ln([1-x])^2*x*[1-z]^-1
       + 1/2*ln([1-x])^2*x
       + 4*ln([1-x])*ln(z)
       + 6*ln([1-x])*ln(z)*[1-x]^-1*[1-z]^-1
       - 4*ln([1-x])*ln(z)*[1-x]^-1
       - 4*ln([1-x])*ln(z)*[1-x]^-1*z
       - 5*ln([1-x])*ln(z)*[1-z]^-1
       + 6*ln([1-x])*ln(z)*z
       - ln([1-x])*ln(z)*x*[1-z]^-1
       + 2*ln([1-x])*ln(z)*x
       + 2*ln([1-x])*ln([1-z])
       + 5*ln([1-x])*ln([1-z])*[1-x]^-1*[1-z]^-1
       - 4*ln([1-x])*ln([1-z])*[1-x]^-1
       - 4*ln([1-x])*ln([1-z])*[1-x]^-1*z
       - 7/2*ln([1-x])*ln([1-z])*[1-z]^-1
       + 5*ln([1-x])*ln([1-z])*z
       - 3/2*ln([1-x])*ln([1-z])*x*[1-z]^-1
       + 3*ln([1-x])*ln([1-z])*x
       + 6*ln(z)
       + 12*ln(z)*[1-x]^-1*[1-z]^-1
       - 6*ln(z)*[1-x]^-1
       - 10*ln(z)*[1-x]^-1*z
       - 12*ln(z)*[1-z]^-1
       + 12*ln(z)*z
       + 4*ln(z)*x
       + 2*ln(z)^2
       + 4*ln(z)^2*[1-x]^-1*[1-z]^-1
       - 3*ln(z)^2*[1-x]^-1
       - 3*ln(z)^2*[1-x]^-1*z
       - 3*ln(z)^2*[1-z]^-1
       + 4*ln(z)^2*z
       - ln(z)^2*x*[1-z]^-1
       + 2*ln(z)^2*x
       + 2*ln(z)*ln([1-z])
       + 6*ln(z)*ln([1-z])*[1-x]^-1*[1-z]^-1
       - 5*ln(z)*ln([1-z])*[1-x]^-1
       - 5*ln(z)*ln([1-z])*[1-x]^-1*z
       - 4*ln(z)*ln([1-z])*[1-z]^-1
       + 6*ln(z)*ln([1-z])*z
       - 2*ln(z)*ln([1-z])*x*[1-z]^-1
       + 4*ln(z)*ln([1-z])*x
       - ln(z)*ln([x-z])
       - 2*ln(z)*ln([x-z])*[1-x]^-1*[1-z]^-1
       + 3/2*ln(z)*ln([x-z])*[1-x]^-1
       + 3/2*ln(z)*ln([x-z])*[1-x]^-1*z
       + 3/2*ln(z)*ln([x-z])*[1-z]^-1
       - 2*ln(z)*ln([x-z])*z
       + 1/2*ln(z)*ln([x-z])*x*[1-z]^-1
       - ln(z)*ln([x-z])*x
       + 5*ln([1-z])
       + 9*ln([1-z])*[1-x]^-1*[1-z]^-1
       - 5*ln([1-z])*[1-x]^-1
       - 8*ln([1-z])*[1-x]^-1*z
       - 9*ln([1-z])*[1-z]^-1
       + 9*ln([1-z])*z
       + 3*ln([1-z])*x
       + 1/2*ln([1-z])^2
       + 2*ln([1-z])^2*[1-x]^-1*[1-z]^-1
       - 7/4*ln([1-z])^2*[1-x]^-1
       - 7/4*ln([1-z])^2*[1-x]^-1*z
       - 5/4*ln([1-z])^2*[1-z]^-1
       + 2*ln([1-z])^2*z
       - 3/4*ln([1-z])^2*x*[1-z]^-1
       + 3/2*ln([1-z])^2*x
       - ln([1-z])*ln([1-x-z])*[1-x]^-1*[1-z]^-1
       + ln([1-z])*ln([1-x-z])*[1-x]^-1
       + ln([1-z])*ln([1-x-z])*[1-x]^-1*z
       + 1/2*ln([1-z])*ln([1-x-z])*[1-z]^-1
       - ln([1-z])*ln([1-x-z])*z
       + 1/2*ln([1-z])*ln([1-x-z])*x*[1-z]^-1
       - ln([1-z])*ln([1-x-z])*x
       - Li2(x^-1*[1-x]^-1*z*[1-z])*[1-x]^-1*[1-z]^-1
       + Li2(x^-1*[1-x]^-1*z*[1-z])*[1-x]^-1
       + Li2(x^-1*[1-x]^-1*z*[1-z])*[1-x]^-1*z
       + 1/2*Li2(x^-1*[1-x]^-1*z*[1-z])*[1-z]^-1
       - Li2(x^-1*[1-x]^-1*z*[1-z])*z
       + 1/2*Li2(x^-1*[1-x]^-1*z*[1-z])*x*[1-z]^-1
       - Li2(x^-1*[1-x]^-1*z*[1-z])*x
       - Li2(x^-1*[1-x]*z*[1-z]^-1)
       - Li2(x^-1*[1-x]*z*[1-z]^-1)*[1-x]^-1*[1-z]^-1
       + 1/2*Li2(x^-1*[1-x]*z*[1-z]^-1)*[1-x]^-1
       + 1/2*Li2(x^-1*[1-x]*z*[1-z]^-1)*[1-x]^-1*z
       + Li2(x^-1*[1-x]*z*[1-z]^-1)*[1-z]^-1
       - Li2(x^-1*[1-x]*z*[1-z]^-1)*z
       + Li2([1-x]^-1*z)*[1-x]^-1*[1-z]^-1
       - Li2([1-x]^-1*z)*[1-x]^-1
       - Li2([1-x]^-1*z)*[1-x]^-1*z
       - 1/2*Li2([1-x]^-1*z)*[1-z]^-1
       + Li2([1-x]^-1*z)*z
       - 1/2*Li2([1-x]^-1*z)*x*[1-z]^-1
       + Li2([1-x]^-1*z)*x
       + Li2([1-x]*[1-z]^-1)
       + 1/2*Li2([1-x]*[1-z]^-1)*[1-x]^-1
       + 1/2*Li2([1-x]*[1-z]^-1)*[1-x]^-1*z
       - 1/2*Li2([1-x]*[1-z]^-1)*[1-z]^-1
       + 1/2*Li2([1-x]*[1-z]^-1)*x*[1-z]^-1
       - Li2([1-x]*[1-z]^-1)*x
       - Li2(z)
       - Li2(z)*[1-x]^-1*[1-z]^-1
       + 1/2*Li2(z)*[1-x]^-1
       + 1/2*Li2(z)*[1-x]^-1*z
       + Li2(z)*[1-z]^-1
       - Li2(z)*z
       )

       + T(u2)*NC^-2 * (  - 3
       - 4*[1-x]^-1*[1-z]^-1
       - 2*[1-x]^-1
       + 6*[1-x]^-1*z
       + 9*[1-z]^-1
       - 7*z
       - 5*x*[1-z]^-1
       + 5*x
       + 2/3*pi^2*[1-x]^-1*[1-z]^-1
       - 5/12*pi^2*[1-x]^-1
       - 5/12*pi^2*[1-x]^-1*z
       - 7/12*pi^2*[1-z]^-1
       + 1/2*pi^2
       + 2/3*pi^2*z
       - 1/12*pi^2*x*[1-z]^-1
       + 1/6*pi^2*x
       + 9*ln(x)
       + 18*ln(x)*[1-x]^-1*[1-z]^-1
       - 9*ln(x)*[1-x]^-1
       - 15*ln(x)*[1-x]^-1*z
       - 18*ln(x)*[1-z]^-1
       + 18*ln(x)*z
       + 6*ln(x)*x
       - ln(x)*ln( - [1-x-z])*[1-x]^-1*[1-z]^-1
       + ln(x)*ln( - [1-x-z])*[1-x]^-1
       + ln(x)*ln( - [1-x-z])*[1-x]^-1*z
       + 1/2*ln(x)*ln( - [1-x-z])*[1-z]^-1
       - ln(x)*ln( - [1-x-z])*z
       + 1/2*ln(x)*ln( - [1-x-z])*x*[1-z]^-1
       - ln(x)*ln( - [1-x-z])*x
       - 3*ln(x)^2
       - 13/2*ln(x)^2*[1-x]^-1*[1-z]^-1
       + 5*ln(x)^2*[1-x]^-1
       + 5*ln(x)^2*[1-x]^-1*z
       + 19/4*ln(x)^2*[1-z]^-1
       - 13/2*ln(x)^2*z
       + 7/4*ln(x)^2*x*[1-z]^-1
       - 7/2*ln(x)^2*x
       + 5*ln(x)*ln([1-x])
       + 8*ln(x)*ln([1-x])*[1-x]^-1*[1-z]^-1
       - 11/2*ln(x)*ln([1-x])*[1-x]^-1
       - 11/2*ln(x)*ln([1-x])*[1-x]^-1*z
       - 13/2*ln(x)*ln([1-x])*[1-z]^-1
       + 8*ln(x)*ln([1-x])*z
       - 3/2*ln(x)*ln([1-x])*x*[1-z]^-1
       + 3*ln(x)*ln([1-x])*x
       + 5*ln(x)*ln(z)
       + 11*ln(x)*ln(z)*[1-x]^-1*[1-z]^-1
       - 17/2*ln(x)*ln(z)*[1-x]^-1
       - 17/2*ln(x)*ln(z)*[1-x]^-1*z
       - 8*ln(x)*ln(z)*[1-z]^-1
       + 11*ln(x)*ln(z)*z
       - 3*ln(x)*ln(z)*x*[1-z]^-1
       + 6*ln(x)*ln(z)*x
       + 4*ln(x)*ln([1-z])
       + 10*ln(x)*ln([1-z])*[1-x]^-1*[1-z]^-1
       - 8*ln(x)*ln([1-z])*[1-x]^-1
       - 8*ln(x)*ln([1-z])*[1-x]^-1*z
       - 7*ln(x)*ln([1-z])*[1-z]^-1
       + 10*ln(x)*ln([1-z])*z
       - 3*ln(x)*ln([1-z])*x*[1-z]^-1
       + 6*ln(x)*ln([1-z])*x
       - ln(x)*ln([x-z])
       - 2*ln(x)*ln([x-z])*[1-x]^-1*[1-z]^-1
       + 3/2*ln(x)*ln([x-z])*[1-x]^-1
       + 3/2*ln(x)*ln([x-z])*[1-x]^-1*z
       + 3/2*ln(x)*ln([x-z])*[1-z]^-1
       - 2*ln(x)*ln([x-z])*z
       + 1/2*ln(x)*ln([x-z])*x*[1-z]^-1
       - ln(x)*ln([x-z])*x
       - 4*ln([1-x])
       - 9*ln([1-x])*[1-x]^-1*[1-z]^-1
       + 4*ln([1-x])*[1-x]^-1
       + 7*ln([1-x])*[1-x]^-1*z
       + 9*ln([1-x])*[1-z]^-1
       - 9*ln([1-x])*z
       - 3*ln([1-x])*x
       - 2*ln([1-x])^2
       - 5/2*ln([1-x])^2*[1-x]^-1*[1-z]^-1
       + 3/2*ln([1-x])^2*[1-x]^-1
       + 3/2*ln([1-x])^2*[1-x]^-1*z
       + 9/4*ln([1-x])^2*[1-z]^-1
       - 5/2*ln([1-x])^2*z
       + 1/4*ln([1-x])^2*x*[1-z]^-1
       - 1/2*ln([1-x])^2*x
       - 4*ln([1-x])*ln(z)
       - 6*ln([1-x])*ln(z)*[1-x]^-1*[1-z]^-1
       + 4*ln([1-x])*ln(z)*[1-x]^-1
       + 4*ln([1-x])*ln(z)*[1-x]^-1*z
       + 5*ln([1-x])*ln(z)*[1-z]^-1
       - 6*ln([1-x])*ln(z)*z
       + ln([1-x])*ln(z)*x*[1-z]^-1
       - 2*ln([1-x])*ln(z)*x
       - 2*ln([1-x])*ln([1-z])
       - 4*ln([1-x])*ln([1-z])*[1-x]^-1*[1-z]^-1
       + 3*ln([1-x])*ln([1-z])*[1-x]^-1
       + 3*ln([1-x])*ln([1-z])*[1-x]^-1*z
       + 3*ln([1-x])*ln([1-z])*[1-z]^-1
       - 4*ln([1-x])*ln([1-z])*z
       + ln([1-x])*ln([1-z])*x*[1-z]^-1
       - 2*ln([1-x])*ln([1-z])*x
       - 6*ln(z)
       - 12*ln(z)*[1-x]^-1*[1-z]^-1
       + 6*ln(z)*[1-x]^-1
       + 10*ln(z)*[1-x]^-1*z
       + 12*ln(z)*[1-z]^-1
       - 12*ln(z)*z
       - 4*ln(z)*x
       - 2*ln(z)^2
       - 4*ln(z)^2*[1-x]^-1*[1-z]^-1
       + 3*ln(z)^2*[1-x]^-1
       + 3*ln(z)^2*[1-x]^-1*z
       + 3*ln(z)^2*[1-z]^-1
       - 4*ln(z)^2*z
       + ln(z)^2*x*[1-z]^-1
       - 2*ln(z)^2*x
       - 2*ln(z)*ln([1-z])
       - 7*ln(z)*ln([1-z])*[1-x]^-1*[1-z]^-1
       + 6*ln(z)*ln([1-z])*[1-x]^-1
       + 6*ln(z)*ln([1-z])*[1-x]^-1*z
       + 9/2*ln(z)*ln([1-z])*[1-z]^-1
       - 7*ln(z)*ln([1-z])*z
       + 5/2*ln(z)*ln([1-z])*x*[1-z]^-1
       - 5*ln(z)*ln([1-z])*x
       + ln(z)*ln([x-z])
       + 2*ln(z)*ln([x-z])*[1-x]^-1*[1-z]^-1
       - 3/2*ln(z)*ln([x-z])*[1-x]^-1
       - 3/2*ln(z)*ln([x-z])*[1-x]^-1*z
       - 3/2*ln(z)*ln([x-z])*[1-z]^-1
       + 2*ln(z)*ln([x-z])*z
       - 1/2*ln(z)*ln([x-z])*x*[1-z]^-1
       + ln(z)*ln([x-z])*x
       - 5*ln([1-z])
       - 9*ln([1-z])*[1-x]^-1*[1-z]^-1
       + 5*ln([1-z])*[1-x]^-1
       + 8*ln([1-z])*[1-x]^-1*z
       + 9*ln([1-z])*[1-z]^-1
       - 9*ln([1-z])*z
       - 3*ln([1-z])*x
       + ln([1-z])*ln( - [1-x-z])*[1-x]^-1*[1-z]^-1
       - ln([1-z])*ln( - [1-x-z])*[1-x]^-1
       - ln([1-z])*ln( - [1-x-z])*[1-x]^-1*z
       - 1/2*ln([1-z])*ln( - [1-x-z])*[1-z]^-1
       + ln([1-z])*ln( - [1-x-z])*z
       - 1/2*ln([1-z])*ln( - [1-x-z])*x*[1-z]^-1
       + ln([1-z])*ln( - [1-x-z])*x
       - 1/2*ln([1-z])^2
       - 5/2*ln([1-z])^2*[1-x]^-1*[1-z]^-1
       + 9/4*ln([1-z])^2*[1-x]^-1
       + 9/4*ln([1-z])^2*[1-x]^-1*z
       + 3/2*ln([1-z])^2*[1-z]^-1
       - 5/2*ln([1-z])^2*z
       + ln([1-z])^2*x*[1-z]^-1
       - 2*ln([1-z])^2*x
       + Li2(x^-1*[1-x]*z*[1-z]^-1)
       + Li2(x^-1*[1-x]*z*[1-z]^-1)*[1-x]^-1*[1-z]^-1
       - 1/2*Li2(x^-1*[1-x]*z*[1-z]^-1)*[1-x]^-1
       - 1/2*Li2(x^-1*[1-x]*z*[1-z]^-1)*[1-x]^-1*z
       - Li2(x^-1*[1-x]*z*[1-z]^-1)*[1-z]^-1
       + Li2(x^-1*[1-x]*z*[1-z]^-1)*z
       + Li2([1-x]*z^-1)*[1-x]^-1*[1-z]^-1
       - Li2([1-x]*z^-1)*[1-x]^-1
       - Li2([1-x]*z^-1)*[1-x]^-1*z
       - 1/2*Li2([1-x]*z^-1)*[1-z]^-1
       + Li2([1-x]*z^-1)*z
       - 1/2*Li2([1-x]*z^-1)*x*[1-z]^-1
       + Li2([1-x]*z^-1)*x
       - Li2([1-x]*[1-z]^-1)
       - 1/2*Li2([1-x]*[1-z]^-1)*[1-x]^-1
       - 1/2*Li2([1-x]*[1-z]^-1)*[1-x]^-1*z
       + 1/2*Li2([1-x]*[1-z]^-1)*[1-z]^-1
       - 1/2*Li2([1-x]*[1-z]^-1)*x*[1-z]^-1
       + Li2([1-x]*[1-z]^-1)*x
       - Li2(x*[1-x]*z^-1*[1-z]^-1)*[1-x]^-1*[1-z]^-1
       + Li2(x*[1-x]*z^-1*[1-z]^-1)*[1-x]^-1
       + Li2(x*[1-x]*z^-1*[1-z]^-1)*[1-x]^-1*z
       + 1/2*Li2(x*[1-x]*z^-1*[1-z]^-1)*[1-z]^-1
       - Li2(x*[1-x]*z^-1*[1-z]^-1)*z
       + 1/2*Li2(x*[1-x]*z^-1*[1-z]^-1)*x*[1-z]^-1
       - Li2(x*[1-x]*z^-1*[1-z]^-1)*x
       + Li2(z)
       + Li2(z)*[1-x]^-1*[1-z]^-1
       - 1/2*Li2(z)*[1-x]^-1
       - 1/2*Li2(z)*[1-x]^-1*z
       - Li2(z)*[1-z]^-1
       + Li2(z)*z
       )

       + T(u2) * ( 3
       + 4*[1-x]^-1*[1-z]^-1
       + 2*[1-x]^-1
       - 6*[1-x]^-1*z
       - 9*[1-z]^-1
       + 7*z
       + 5*x*[1-z]^-1
       - 5*x
       - 2/3*pi^2*[1-x]^-1*[1-z]^-1
       + 5/12*pi^2*[1-x]^-1
       + 5/12*pi^2*[1-x]^-1*z
       + 7/12*pi^2*[1-z]^-1
       - 1/2*pi^2
       - 2/3*pi^2*z
       + 1/12*pi^2*x*[1-z]^-1
       - 1/6*pi^2*x
       - 9*ln(x)
       - 18*ln(x)*[1-x]^-1*[1-z]^-1
       + 9*ln(x)*[1-x]^-1
       + 15*ln(x)*[1-x]^-1*z
       + 18*ln(x)*[1-z]^-1
       - 18*ln(x)*z
       - 6*ln(x)*x
       + ln(x)*ln( - [1-x-z])*[1-x]^-1*[1-z]^-1
       - ln(x)*ln( - [1-x-z])*[1-x]^-1
       - ln(x)*ln( - [1-x-z])*[1-x]^-1*z
       - 1/2*ln(x)*ln( - [1-x-z])*[1-z]^-1
       + ln(x)*ln( - [1-x-z])*z
       - 1/2*ln(x)*ln( - [1-x-z])*x*[1-z]^-1
       + ln(x)*ln( - [1-x-z])*x
       + 3*ln(x)^2
       + 13/2*ln(x)^2*[1-x]^-1*[1-z]^-1
       - 5*ln(x)^2*[1-x]^-1
       - 5*ln(x)^2*[1-x]^-1*z
       - 19/4*ln(x)^2*[1-z]^-1
       + 13/2*ln(x)^2*z
       - 7/4*ln(x)^2*x*[1-z]^-1
       + 7/2*ln(x)^2*x
       - 5*ln(x)*ln([1-x])
       - 8*ln(x)*ln([1-x])*[1-x]^-1*[1-z]^-1
       + 11/2*ln(x)*ln([1-x])*[1-x]^-1
       + 11/2*ln(x)*ln([1-x])*[1-x]^-1*z
       + 13/2*ln(x)*ln([1-x])*[1-z]^-1
       - 8*ln(x)*ln([1-x])*z
       + 3/2*ln(x)*ln([1-x])*x*[1-z]^-1
       - 3*ln(x)*ln([1-x])*x
       - 5*ln(x)*ln(z)
       - 11*ln(x)*ln(z)*[1-x]^-1*[1-z]^-1
       + 17/2*ln(x)*ln(z)*[1-x]^-1
       + 17/2*ln(x)*ln(z)*[1-x]^-1*z
       + 8*ln(x)*ln(z)*[1-z]^-1
       - 11*ln(x)*ln(z)*z
       + 3*ln(x)*ln(z)*x*[1-z]^-1
       - 6*ln(x)*ln(z)*x
       - 4*ln(x)*ln([1-z])
       - 10*ln(x)*ln([1-z])*[1-x]^-1*[1-z]^-1
       + 8*ln(x)*ln([1-z])*[1-x]^-1
       + 8*ln(x)*ln([1-z])*[1-x]^-1*z
       + 7*ln(x)*ln([1-z])*[1-z]^-1
       - 10*ln(x)*ln([1-z])*z
       + 3*ln(x)*ln([1-z])*x*[1-z]^-1
       - 6*ln(x)*ln([1-z])*x
       + ln(x)*ln([x-z])
       + 2*ln(x)*ln([x-z])*[1-x]^-1*[1-z]^-1
       - 3/2*ln(x)*ln([x-z])*[1-x]^-1
       - 3/2*ln(x)*ln([x-z])*[1-x]^-1*z
       - 3/2*ln(x)*ln([x-z])*[1-z]^-1
       + 2*ln(x)*ln([x-z])*z
       - 1/2*ln(x)*ln([x-z])*x*[1-z]^-1
       + ln(x)*ln([x-z])*x
       + 4*ln([1-x])
       + 9*ln([1-x])*[1-x]^-1*[1-z]^-1
       - 4*ln([1-x])*[1-x]^-1
       - 7*ln([1-x])*[1-x]^-1*z
       - 9*ln([1-x])*[1-z]^-1
       + 9*ln([1-x])*z
       + 3*ln([1-x])*x
       + 2*ln([1-x])^2
       + 5/2*ln([1-x])^2*[1-x]^-1*[1-z]^-1
       - 3/2*ln([1-x])^2*[1-x]^-1
       - 3/2*ln([1-x])^2*[1-x]^-1*z
       - 9/4*ln([1-x])^2*[1-z]^-1
       + 5/2*ln([1-x])^2*z
       - 1/4*ln([1-x])^2*x*[1-z]^-1
       + 1/2*ln([1-x])^2*x
       + 4*ln([1-x])*ln(z)
       + 6*ln([1-x])*ln(z)*[1-x]^-1*[1-z]^-1
       - 4*ln([1-x])*ln(z)*[1-x]^-1
       - 4*ln([1-x])*ln(z)*[1-x]^-1*z
       - 5*ln([1-x])*ln(z)*[1-z]^-1
       + 6*ln([1-x])*ln(z)*z
       - ln([1-x])*ln(z)*x*[1-z]^-1
       + 2*ln([1-x])*ln(z)*x
       + 2*ln([1-x])*ln([1-z])
       + 4*ln([1-x])*ln([1-z])*[1-x]^-1*[1-z]^-1
       - 3*ln([1-x])*ln([1-z])*[1-x]^-1
       - 3*ln([1-x])*ln([1-z])*[1-x]^-1*z
       - 3*ln([1-x])*ln([1-z])*[1-z]^-1
       + 4*ln([1-x])*ln([1-z])*z
       - ln([1-x])*ln([1-z])*x*[1-z]^-1
       + 2*ln([1-x])*ln([1-z])*x
       + 6*ln(z)
       + 12*ln(z)*[1-x]^-1*[1-z]^-1
       - 6*ln(z)*[1-x]^-1
       - 10*ln(z)*[1-x]^-1*z
       - 12*ln(z)*[1-z]^-1
       + 12*ln(z)*z
       + 4*ln(z)*x
       + 2*ln(z)^2
       + 4*ln(z)^2*[1-x]^-1*[1-z]^-1
       - 3*ln(z)^2*[1-x]^-1
       - 3*ln(z)^2*[1-x]^-1*z
       - 3*ln(z)^2*[1-z]^-1
       + 4*ln(z)^2*z
       - ln(z)^2*x*[1-z]^-1
       + 2*ln(z)^2*x
       + 2*ln(z)*ln([1-z])
       + 7*ln(z)*ln([1-z])*[1-x]^-1*[1-z]^-1
       - 6*ln(z)*ln([1-z])*[1-x]^-1
       - 6*ln(z)*ln([1-z])*[1-x]^-1*z
       - 9/2*ln(z)*ln([1-z])*[1-z]^-1
       + 7*ln(z)*ln([1-z])*z
       - 5/2*ln(z)*ln([1-z])*x*[1-z]^-1
       + 5*ln(z)*ln([1-z])*x
       - ln(z)*ln([x-z])
       - 2*ln(z)*ln([x-z])*[1-x]^-1*[1-z]^-1
       + 3/2*ln(z)*ln([x-z])*[1-x]^-1
       + 3/2*ln(z)*ln([x-z])*[1-x]^-1*z
       + 3/2*ln(z)*ln([x-z])*[1-z]^-1
       - 2*ln(z)*ln([x-z])*z
       + 1/2*ln(z)*ln([x-z])*x*[1-z]^-1
       - ln(z)*ln([x-z])*x
       + 5*ln([1-z])
       + 9*ln([1-z])*[1-x]^-1*[1-z]^-1
       - 5*ln([1-z])*[1-x]^-1
       - 8*ln([1-z])*[1-x]^-1*z
       - 9*ln([1-z])*[1-z]^-1
       + 9*ln([1-z])*z
       + 3*ln([1-z])*x
       - ln([1-z])*ln( - [1-x-z])*[1-x]^-1*[1-z]^-1
       + ln([1-z])*ln( - [1-x-z])*[1-x]^-1
       + ln([1-z])*ln( - [1-x-z])*[1-x]^-1*z
       + 1/2*ln([1-z])*ln( - [1-x-z])*[1-z]^-1
       - ln([1-z])*ln( - [1-x-z])*z
       + 1/2*ln([1-z])*ln( - [1-x-z])*x*[1-z]^-1
       - ln([1-z])*ln( - [1-x-z])*x
       + 1/2*ln([1-z])^2
       + 5/2*ln([1-z])^2*[1-x]^-1*[1-z]^-1
       - 9/4*ln([1-z])^2*[1-x]^-1
       - 9/4*ln([1-z])^2*[1-x]^-1*z
       - 3/2*ln([1-z])^2*[1-z]^-1
       + 5/2*ln([1-z])^2*z
       - ln([1-z])^2*x*[1-z]^-1
       + 2*ln([1-z])^2*x
       - Li2(x^-1*[1-x]*z*[1-z]^-1)
       - Li2(x^-1*[1-x]*z*[1-z]^-1)*[1-x]^-1*[1-z]^-1
       + 1/2*Li2(x^-1*[1-x]*z*[1-z]^-1)*[1-x]^-1
       + 1/2*Li2(x^-1*[1-x]*z*[1-z]^-1)*[1-x]^-1*z
       + Li2(x^-1*[1-x]*z*[1-z]^-1)*[1-z]^-1
       - Li2(x^-1*[1-x]*z*[1-z]^-1)*z
       - Li2([1-x]*z^-1)*[1-x]^-1*[1-z]^-1
       + Li2([1-x]*z^-1)*[1-x]^-1
       + Li2([1-x]*z^-1)*[1-x]^-1*z
       + 1/2*Li2([1-x]*z^-1)*[1-z]^-1
       - Li2([1-x]*z^-1)*z
       + 1/2*Li2([1-x]*z^-1)*x*[1-z]^-1
       - Li2([1-x]*z^-1)*x
       + Li2([1-x]*[1-z]^-1)
       + 1/2*Li2([1-x]*[1-z]^-1)*[1-x]^-1
       + 1/2*Li2([1-x]*[1-z]^-1)*[1-x]^-1*z
       - 1/2*Li2([1-x]*[1-z]^-1)*[1-z]^-1
       + 1/2*Li2([1-x]*[1-z]^-1)*x*[1-z]^-1
       - Li2([1-x]*[1-z]^-1)*x
       + Li2(x*[1-x]*z^-1*[1-z]^-1)*[1-x]^-1*[1-z]^-1
       - Li2(x*[1-x]*z^-1*[1-z]^-1)*[1-x]^-1
       - Li2(x*[1-x]*z^-1*[1-z]^-1)*[1-x]^-1*z
       - 1/2*Li2(x*[1-x]*z^-1*[1-z]^-1)*[1-z]^-1
       + Li2(x*[1-x]*z^-1*[1-z]^-1)*z
       - 1/2*Li2(x*[1-x]*z^-1*[1-z]^-1)*x*[1-z]^-1
       + Li2(x*[1-x]*z^-1*[1-z]^-1)*x
       - Li2(z)
       - Li2(z)*[1-x]^-1*[1-z]^-1
       + 1/2*Li2(z)*[1-x]^-1
       + 1/2*Li2(z)*[1-x]^-1*z
       + Li2(z)*[1-z]^-1
       - Li2(z)*z
       )

       + T(u3)*NC^-2 * (  - 3
       - 4*[1-x]^-1*[1-z]^-1
       - 2*[1-x]^-1
       + 6*[1-x]^-1*z
       + 9*[1-z]^-1
       - 7*z
       - 5*x*[1-z]^-1
       + 5*x
       + 4/3*pi^2*[1-x]^-1*[1-z]^-1
       - 13/12*pi^2*[1-x]^-1
       - 13/12*pi^2*[1-x]^-1*z
       - 11/12*pi^2*[1-z]^-1
       + 1/2*pi^2
       + 4/3*pi^2*z
       - 5/12*pi^2*x*[1-z]^-1
       + 5/6*pi^2*x
       + 9*ln(x)
       + 18*ln(x)*[1-x]^-1*[1-z]^-1
       - 9*ln(x)*[1-x]^-1
       - 15*ln(x)*[1-x]^-1*z
       - 18*ln(x)*[1-z]^-1
       + 18*ln(x)*z
       + 6*ln(x)*x
       - ln(x)*ln( - [x-z])
       - 2*ln(x)*ln( - [x-z])*[1-x]^-1*[1-z]^-1
       + 3/2*ln(x)*ln( - [x-z])*[1-x]^-1
       + 3/2*ln(x)*ln( - [x-z])*[1-x]^-1*z
       + 3/2*ln(x)*ln( - [x-z])*[1-z]^-1
       - 2*ln(x)*ln( - [x-z])*z
       + 1/2*ln(x)*ln( - [x-z])*x*[1-z]^-1
       - ln(x)*ln( - [x-z])*x
       - 7/2*ln(x)^2
       - 7*ln(x)^2*[1-x]^-1*[1-z]^-1
       + 21/4*ln(x)^2*[1-x]^-1
       + 21/4*ln(x)^2*[1-x]^-1*z
       + 21/4*ln(x)^2*[1-z]^-1
       - 7*ln(x)^2*z
       + 7/4*ln(x)^2*x*[1-z]^-1
       - 7/2*ln(x)^2*x
       + 6*ln(x)*ln([1-x])
       + 9*ln(x)*ln([1-x])*[1-x]^-1*[1-z]^-1
       - 6*ln(x)*ln([1-x])*[1-x]^-1
       - 6*ln(x)*ln([1-x])*[1-x]^-1*z
       - 15/2*ln(x)*ln([1-x])*[1-z]^-1
       + 9*ln(x)*ln([1-x])*z
       - 3/2*ln(x)*ln([1-x])*x*[1-z]^-1
       + 3*ln(x)*ln([1-x])*x
       + 6*ln(x)*ln(z)
       + 12*ln(x)*ln(z)*[1-x]^-1*[1-z]^-1
       - 9*ln(x)*ln(z)*[1-x]^-1
       - 9*ln(x)*ln(z)*[1-x]^-1*z
       - 9*ln(x)*ln(z)*[1-z]^-1
       + 12*ln(x)*ln(z)*z
       - 3*ln(x)*ln(z)*x*[1-z]^-1
       + 6*ln(x)*ln(z)*x
       + 3*ln(x)*ln([1-z])
       + 9*ln(x)*ln([1-z])*[1-x]^-1*[1-z]^-1
       - 15/2*ln(x)*ln([1-z])*[1-x]^-1
       - 15/2*ln(x)*ln([1-z])*[1-x]^-1*z
       - 6*ln(x)*ln([1-z])*[1-z]^-1
       + 9*ln(x)*ln([1-z])*z
       - 3*ln(x)*ln([1-z])*x*[1-z]^-1
       + 6*ln(x)*ln([1-z])*x
       - ln(x)*ln([1-x-z])*[1-x]^-1*[1-z]^-1
       + ln(x)*ln([1-x-z])*[1-x]^-1
       + ln(x)*ln([1-x-z])*[1-x]^-1*z
       + 1/2*ln(x)*ln([1-x-z])*[1-z]^-1
       - ln(x)*ln([1-x-z])*z
       + 1/2*ln(x)*ln([1-x-z])*x*[1-z]^-1
       - ln(x)*ln([1-x-z])*x
       - 4*ln([1-x])
       - 9*ln([1-x])*[1-x]^-1*[1-z]^-1
       + 4*ln([1-x])*[1-x]^-1
       + 7*ln([1-x])*[1-x]^-1*z
       + 9*ln([1-x])*[1-z]^-1
       - 9*ln([1-x])*z
       - 3*ln([1-x])*x
       - 2*ln([1-x])^2
       - 7/2*ln([1-x])^2*[1-x]^-1*[1-z]^-1
       + 5/2*ln([1-x])^2*[1-x]^-1
       + 5/2*ln([1-x])^2*[1-x]^-1*z
       + 11/4*ln([1-x])^2*[1-z]^-1
       - 7/2*ln([1-x])^2*z
       + 3/4*ln([1-x])^2*x*[1-z]^-1
       - 3/2*ln([1-x])^2*x
       - 5*ln([1-x])*ln(z)
       - 6*ln([1-x])*ln(z)*[1-x]^-1*[1-z]^-1
       + 7/2*ln([1-x])*ln(z)*[1-x]^-1
       + 7/2*ln([1-x])*ln(z)*[1-x]^-1*z
       + 11/2*ln([1-x])*ln(z)*[1-z]^-1
       - 6*ln([1-x])*ln(z)*z
       + 1/2*ln([1-x])*ln(z)*x*[1-z]^-1
       - ln([1-x])*ln(z)*x
       - 2*ln([1-x])*ln([1-z])
       - 3*ln([1-x])*ln([1-z])*[1-x]^-1*[1-z]^-1
       + 2*ln([1-x])*ln([1-z])*[1-x]^-1
       + 2*ln([1-x])*ln([1-z])*[1-x]^-1*z
       + 5/2*ln([1-x])*ln([1-z])*[1-z]^-1
       - 3*ln([1-x])*ln([1-z])*z
       + 1/2*ln([1-x])*ln([1-z])*x*[1-z]^-1
       - ln([1-x])*ln([1-z])*x
       - 6*ln(z)
       - 12*ln(z)*[1-x]^-1*[1-z]^-1
       + 6*ln(z)*[1-x]^-1
       + 10*ln(z)*[1-x]^-1*z
       + 12*ln(z)*[1-z]^-1
       - 12*ln(z)*z
       - 4*ln(z)*x
       + ln(z)*ln( - [x-z])
       + 2*ln(z)*ln( - [x-z])*[1-x]^-1*[1-z]^-1
       - 3/2*ln(z)*ln( - [x-z])*[1-x]^-1
       - 3/2*ln(z)*ln( - [x-z])*[1-x]^-1*z
       - 3/2*ln(z)*ln( - [x-z])*[1-z]^-1
       + 2*ln(z)*ln( - [x-z])*z
       - 1/2*ln(z)*ln( - [x-z])*x*[1-z]^-1
       + ln(z)*ln( - [x-z])*x
       - 5/2*ln(z)^2
       - 5*ln(z)^2*[1-x]^-1*[1-z]^-1
       + 15/4*ln(z)^2*[1-x]^-1
       + 15/4*ln(z)^2*[1-x]^-1*z
       + 15/4*ln(z)^2*[1-z]^-1
       - 5*ln(z)^2*z
       + 5/4*ln(z)^2*x*[1-z]^-1
       - 5/2*ln(z)^2*x
       - ln(z)*ln([1-z])
       - 6*ln(z)*ln([1-z])*[1-x]^-1*[1-z]^-1
       + 11/2*ln(z)*ln([1-z])*[1-x]^-1
       + 11/2*ln(z)*ln([1-z])*[1-x]^-1*z
       + 7/2*ln(z)*ln([1-z])*[1-z]^-1
       - 6*ln(z)*ln([1-z])*z
       + 5/2*ln(z)*ln([1-z])*x*[1-z]^-1
       - 5*ln(z)*ln([1-z])*x
       - 5*ln([1-z])
       - 9*ln([1-z])*[1-x]^-1*[1-z]^-1
       + 5*ln([1-z])*[1-x]^-1
       + 8*ln([1-z])*[1-x]^-1*z
       + 9*ln([1-z])*[1-z]^-1
       - 9*ln([1-z])*z
       - 3*ln([1-z])*x
       - 1/2*ln([1-z])^2
       - 3*ln([1-z])^2*[1-x]^-1*[1-z]^-1
       + 11/4*ln([1-z])^2*[1-x]^-1
       + 11/4*ln([1-z])^2*[1-x]^-1*z
       + 7/4*ln([1-z])^2*[1-z]^-1
       - 3*ln([1-z])^2*z
       + 5/4*ln([1-z])^2*x*[1-z]^-1
       - 5/2*ln([1-z])^2*x
       + ln([1-z])*ln([1-x-z])*[1-x]^-1*[1-z]^-1
       - ln([1-z])*ln([1-x-z])*[1-x]^-1
       - ln([1-z])*ln([1-x-z])*[1-x]^-1*z
       - 1/2*ln([1-z])*ln([1-x-z])*[1-z]^-1
       + ln([1-z])*ln([1-x-z])*z
       - 1/2*ln([1-z])*ln([1-x-z])*x*[1-z]^-1
       + ln([1-z])*ln([1-x-z])*x
       + Li2([1-x]^-1*[1-z])
       + 1/2*Li2([1-x]^-1*[1-z])*[1-x]^-1
       + 1/2*Li2([1-x]^-1*[1-z])*[1-x]^-1*z
       - 1/2*Li2([1-x]^-1*[1-z])*[1-z]^-1
       + 1/2*Li2([1-x]^-1*[1-z])*x*[1-z]^-1
       - Li2([1-x]^-1*[1-z])*x
       - Li2([1-x]^-1*z)*[1-x]^-1*[1-z]^-1
       + Li2([1-x]^-1*z)*[1-x]^-1
       + Li2([1-x]^-1*z)*[1-x]^-1*z
       + 1/2*Li2([1-x]^-1*z)*[1-z]^-1
       - Li2([1-x]^-1*z)*z
       + 1/2*Li2([1-x]^-1*z)*x*[1-z]^-1
       - Li2([1-x]^-1*z)*x
       - Li2(x*[1-x]^-1*z^-1*[1-z])
       - Li2(x*[1-x]^-1*z^-1*[1-z])*[1-x]^-1*[1-z]^-1
       + 1/2*Li2(x*[1-x]^-1*z^-1*[1-z])*[1-x]^-1
       + 1/2*Li2(x*[1-x]^-1*z^-1*[1-z])*[1-x]^-1*z
       + Li2(x*[1-x]^-1*z^-1*[1-z])*[1-z]^-1
       - Li2(x*[1-x]^-1*z^-1*[1-z])*z
       - Li2(x*[1-x]*z^-1*[1-z]^-1)*[1-x]^-1*[1-z]^-1
       + Li2(x*[1-x]*z^-1*[1-z]^-1)*[1-x]^-1
       + Li2(x*[1-x]*z^-1*[1-z]^-1)*[1-x]^-1*z
       + 1/2*Li2(x*[1-x]*z^-1*[1-z]^-1)*[1-z]^-1
       - Li2(x*[1-x]*z^-1*[1-z]^-1)*z
       + 1/2*Li2(x*[1-x]*z^-1*[1-z]^-1)*x*[1-z]^-1
       - Li2(x*[1-x]*z^-1*[1-z]^-1)*x
       + Li2(z)
       + Li2(z)*[1-x]^-1*[1-z]^-1
       - 1/2*Li2(z)*[1-x]^-1
       - 1/2*Li2(z)*[1-x]^-1*z
       - Li2(z)*[1-z]^-1
       + Li2(z)*z
       )

       + T(u3) * ( 3
       + 4*[1-x]^-1*[1-z]^-1
       + 2*[1-x]^-1
       - 6*[1-x]^-1*z
       - 9*[1-z]^-1
       + 7*z
       + 5*x*[1-z]^-1
       - 5*x
       - 4/3*pi^2*[1-x]^-1*[1-z]^-1
       + 13/12*pi^2*[1-x]^-1
       + 13/12*pi^2*[1-x]^-1*z
       + 11/12*pi^2*[1-z]^-1
       - 1/2*pi^2
       - 4/3*pi^2*z
       + 5/12*pi^2*x*[1-z]^-1
       - 5/6*pi^2*x
       - 9*ln(x)
       - 18*ln(x)*[1-x]^-1*[1-z]^-1
       + 9*ln(x)*[1-x]^-1
       + 15*ln(x)*[1-x]^-1*z
       + 18*ln(x)*[1-z]^-1
       - 18*ln(x)*z
       - 6*ln(x)*x
       + ln(x)*ln( - [x-z])
       + 2*ln(x)*ln( - [x-z])*[1-x]^-1*[1-z]^-1
       - 3/2*ln(x)*ln( - [x-z])*[1-x]^-1
       - 3/2*ln(x)*ln( - [x-z])*[1-x]^-1*z
       - 3/2*ln(x)*ln( - [x-z])*[1-z]^-1
       + 2*ln(x)*ln( - [x-z])*z
       - 1/2*ln(x)*ln( - [x-z])*x*[1-z]^-1
       + ln(x)*ln( - [x-z])*x
       + 7/2*ln(x)^2
       + 7*ln(x)^2*[1-x]^-1*[1-z]^-1
       - 21/4*ln(x)^2*[1-x]^-1
       - 21/4*ln(x)^2*[1-x]^-1*z
       - 21/4*ln(x)^2*[1-z]^-1
       + 7*ln(x)^2*z
       - 7/4*ln(x)^2*x*[1-z]^-1
       + 7/2*ln(x)^2*x
       - 6*ln(x)*ln([1-x])
       - 9*ln(x)*ln([1-x])*[1-x]^-1*[1-z]^-1
       + 6*ln(x)*ln([1-x])*[1-x]^-1
       + 6*ln(x)*ln([1-x])*[1-x]^-1*z
       + 15/2*ln(x)*ln([1-x])*[1-z]^-1
       - 9*ln(x)*ln([1-x])*z
       + 3/2*ln(x)*ln([1-x])*x*[1-z]^-1
       - 3*ln(x)*ln([1-x])*x
       - 6*ln(x)*ln(z)
       - 12*ln(x)*ln(z)*[1-x]^-1*[1-z]^-1
       + 9*ln(x)*ln(z)*[1-x]^-1
       + 9*ln(x)*ln(z)*[1-x]^-1*z
       + 9*ln(x)*ln(z)*[1-z]^-1
       - 12*ln(x)*ln(z)*z
       + 3*ln(x)*ln(z)*x*[1-z]^-1
       - 6*ln(x)*ln(z)*x
       - 3*ln(x)*ln([1-z])
       - 9*ln(x)*ln([1-z])*[1-x]^-1*[1-z]^-1
       + 15/2*ln(x)*ln([1-z])*[1-x]^-1
       + 15/2*ln(x)*ln([1-z])*[1-x]^-1*z
       + 6*ln(x)*ln([1-z])*[1-z]^-1
       - 9*ln(x)*ln([1-z])*z
       + 3*ln(x)*ln([1-z])*x*[1-z]^-1
       - 6*ln(x)*ln([1-z])*x
       + ln(x)*ln([1-x-z])*[1-x]^-1*[1-z]^-1
       - ln(x)*ln([1-x-z])*[1-x]^-1
       - ln(x)*ln([1-x-z])*[1-x]^-1*z
       - 1/2*ln(x)*ln([1-x-z])*[1-z]^-1
       + ln(x)*ln([1-x-z])*z
       - 1/2*ln(x)*ln([1-x-z])*x*[1-z]^-1
       + ln(x)*ln([1-x-z])*x
       + 4*ln([1-x])
       + 9*ln([1-x])*[1-x]^-1*[1-z]^-1
       - 4*ln([1-x])*[1-x]^-1
       - 7*ln([1-x])*[1-x]^-1*z
       - 9*ln([1-x])*[1-z]^-1
       + 9*ln([1-x])*z
       + 3*ln([1-x])*x
       + 2*ln([1-x])^2
       + 7/2*ln([1-x])^2*[1-x]^-1*[1-z]^-1
       - 5/2*ln([1-x])^2*[1-x]^-1
       - 5/2*ln([1-x])^2*[1-x]^-1*z
       - 11/4*ln([1-x])^2*[1-z]^-1
       + 7/2*ln([1-x])^2*z
       - 3/4*ln([1-x])^2*x*[1-z]^-1
       + 3/2*ln([1-x])^2*x
       + 5*ln([1-x])*ln(z)
       + 6*ln([1-x])*ln(z)*[1-x]^-1*[1-z]^-1
       - 7/2*ln([1-x])*ln(z)*[1-x]^-1
       - 7/2*ln([1-x])*ln(z)*[1-x]^-1*z
       - 11/2*ln([1-x])*ln(z)*[1-z]^-1
       + 6*ln([1-x])*ln(z)*z
       - 1/2*ln([1-x])*ln(z)*x*[1-z]^-1
       + ln([1-x])*ln(z)*x
       + 2*ln([1-x])*ln([1-z])
       + 3*ln([1-x])*ln([1-z])*[1-x]^-1*[1-z]^-1
       - 2*ln([1-x])*ln([1-z])*[1-x]^-1
       - 2*ln([1-x])*ln([1-z])*[1-x]^-1*z
       - 5/2*ln([1-x])*ln([1-z])*[1-z]^-1
       + 3*ln([1-x])*ln([1-z])*z
       - 1/2*ln([1-x])*ln([1-z])*x*[1-z]^-1
       + ln([1-x])*ln([1-z])*x
       + 6*ln(z)
       + 12*ln(z)*[1-x]^-1*[1-z]^-1
       - 6*ln(z)*[1-x]^-1
       - 10*ln(z)*[1-x]^-1*z
       - 12*ln(z)*[1-z]^-1
       + 12*ln(z)*z
       + 4*ln(z)*x
       - ln(z)*ln( - [x-z])
       - 2*ln(z)*ln( - [x-z])*[1-x]^-1*[1-z]^-1
       + 3/2*ln(z)*ln( - [x-z])*[1-x]^-1
       + 3/2*ln(z)*ln( - [x-z])*[1-x]^-1*z
       + 3/2*ln(z)*ln( - [x-z])*[1-z]^-1
       - 2*ln(z)*ln( - [x-z])*z
       + 1/2*ln(z)*ln( - [x-z])*x*[1-z]^-1
       - ln(z)*ln( - [x-z])*x
       + 5/2*ln(z)^2
       + 5*ln(z)^2*[1-x]^-1*[1-z]^-1
       - 15/4*ln(z)^2*[1-x]^-1
       - 15/4*ln(z)^2*[1-x]^-1*z
       - 15/4*ln(z)^2*[1-z]^-1
       + 5*ln(z)^2*z
       - 5/4*ln(z)^2*x*[1-z]^-1
       + 5/2*ln(z)^2*x
       + ln(z)*ln([1-z])
       + 6*ln(z)*ln([1-z])*[1-x]^-1*[1-z]^-1
       - 11/2*ln(z)*ln([1-z])*[1-x]^-1
       - 11/2*ln(z)*ln([1-z])*[1-x]^-1*z
       - 7/2*ln(z)*ln([1-z])*[1-z]^-1
       + 6*ln(z)*ln([1-z])*z
       - 5/2*ln(z)*ln([1-z])*x*[1-z]^-1
       + 5*ln(z)*ln([1-z])*x
       + 5*ln([1-z])
       + 9*ln([1-z])*[1-x]^-1*[1-z]^-1
       - 5*ln([1-z])*[1-x]^-1
       - 8*ln([1-z])*[1-x]^-1*z
       - 9*ln([1-z])*[1-z]^-1
       + 9*ln([1-z])*z
       + 3*ln([1-z])*x
       + 1/2*ln([1-z])^2
       + 3*ln([1-z])^2*[1-x]^-1*[1-z]^-1
       - 11/4*ln([1-z])^2*[1-x]^-1
       - 11/4*ln([1-z])^2*[1-x]^-1*z
       - 7/4*ln([1-z])^2*[1-z]^-1
       + 3*ln([1-z])^2*z
       - 5/4*ln([1-z])^2*x*[1-z]^-1
       + 5/2*ln([1-z])^2*x
       - ln([1-z])*ln([1-x-z])*[1-x]^-1*[1-z]^-1
       + ln([1-z])*ln([1-x-z])*[1-x]^-1
       + ln([1-z])*ln([1-x-z])*[1-x]^-1*z
       + 1/2*ln([1-z])*ln([1-x-z])*[1-z]^-1
       - ln([1-z])*ln([1-x-z])*z
       + 1/2*ln([1-z])*ln([1-x-z])*x*[1-z]^-1
       - ln([1-z])*ln([1-x-z])*x
       - Li2([1-x]^-1*[1-z])
       - 1/2*Li2([1-x]^-1*[1-z])*[1-x]^-1
       - 1/2*Li2([1-x]^-1*[1-z])*[1-x]^-1*z
       + 1/2*Li2([1-x]^-1*[1-z])*[1-z]^-1
       - 1/2*Li2([1-x]^-1*[1-z])*x*[1-z]^-1
       + Li2([1-x]^-1*[1-z])*x
       + Li2([1-x]^-1*z)*[1-x]^-1*[1-z]^-1
       - Li2([1-x]^-1*z)*[1-x]^-1
       - Li2([1-x]^-1*z)*[1-x]^-1*z
       - 1/2*Li2([1-x]^-1*z)*[1-z]^-1
       + Li2([1-x]^-1*z)*z
       - 1/2*Li2([1-x]^-1*z)*x*[1-z]^-1
       + Li2([1-x]^-1*z)*x
       + Li2(x*[1-x]^-1*z^-1*[1-z])
       + Li2(x*[1-x]^-1*z^-1*[1-z])*[1-x]^-1*[1-z]^-1
       - 1/2*Li2(x*[1-x]^-1*z^-1*[1-z])*[1-x]^-1
       - 1/2*Li2(x*[1-x]^-1*z^-1*[1-z])*[1-x]^-1*z
       - Li2(x*[1-x]^-1*z^-1*[1-z])*[1-z]^-1
       + Li2(x*[1-x]^-1*z^-1*[1-z])*z
       + Li2(x*[1-x]*z^-1*[1-z]^-1)*[1-x]^-1*[1-z]^-1
       - Li2(x*[1-x]*z^-1*[1-z]^-1)*[1-x]^-1
       - Li2(x*[1-x]*z^-1*[1-z]^-1)*[1-x]^-1*z
       - 1/2*Li2(x*[1-x]*z^-1*[1-z]^-1)*[1-z]^-1
       + Li2(x*[1-x]*z^-1*[1-z]^-1)*z
       - 1/2*Li2(x*[1-x]*z^-1*[1-z]^-1)*x*[1-z]^-1
       + Li2(x*[1-x]*z^-1*[1-z]^-1)*x
       - Li2(z)
       - Li2(z)*[1-x]^-1*[1-z]^-1
       + 1/2*Li2(z)*[1-x]^-1
       + 1/2*Li2(z)*[1-x]^-1*z
       + Li2(z)*[1-z]^-1
       - Li2(z)*z
       )

       + T(u4)*NC^-2 * (  - 3
       - 4*[1-x]^-1*[1-z]^-1
       - 2*[1-x]^-1
       + 6*[1-x]^-1*z
       + 9*[1-z]^-1
       - 7*z
       - 5*x*[1-z]^-1
       + 5*x
       + 2/3*pi^2*[1-x]^-1*[1-z]^-1
       - 5/12*pi^2*[1-x]^-1
       - 5/12*pi^2*[1-x]^-1*z
       - 7/12*pi^2*[1-z]^-1
       + 1/2*pi^2
       + 2/3*pi^2*z
       - 1/12*pi^2*x*[1-z]^-1
       + 1/6*pi^2*x
       + 9*ln(x)
       + 18*ln(x)*[1-x]^-1*[1-z]^-1
       - 9*ln(x)*[1-x]^-1
       - 15*ln(x)*[1-x]^-1*z
       - 18*ln(x)*[1-z]^-1
       + 18*ln(x)*z
       + 6*ln(x)*x
       - ln(x)*ln( - [x-z])
       - 2*ln(x)*ln( - [x-z])*[1-x]^-1*[1-z]^-1
       + 3/2*ln(x)*ln( - [x-z])*[1-x]^-1
       + 3/2*ln(x)*ln( - [x-z])*[1-x]^-1*z
       + 3/2*ln(x)*ln( - [x-z])*[1-z]^-1
       - 2*ln(x)*ln( - [x-z])*z
       + 1/2*ln(x)*ln( - [x-z])*x*[1-z]^-1
       - ln(x)*ln( - [x-z])*x
       - ln(x)*ln( - [1-x-z])*[1-x]^-1*[1-z]^-1
       + ln(x)*ln( - [1-x-z])*[1-x]^-1
       + ln(x)*ln( - [1-x-z])*[1-x]^-1*z
       + 1/2*ln(x)*ln( - [1-x-z])*[1-z]^-1
       - ln(x)*ln( - [1-x-z])*z
       + 1/2*ln(x)*ln( - [1-x-z])*x*[1-z]^-1
       - ln(x)*ln( - [1-x-z])*x
       - 7/2*ln(x)^2
       - 13/2*ln(x)^2*[1-x]^-1*[1-z]^-1
       + 19/4*ln(x)^2*[1-x]^-1
       + 19/4*ln(x)^2*[1-x]^-1*z
       + 5*ln(x)^2*[1-z]^-1
       - 13/2*ln(x)^2*z
       + 3/2*ln(x)^2*x*[1-z]^-1
       - 3*ln(x)^2*x
       + 6*ln(x)*ln([1-x])
       + 10*ln(x)*ln([1-x])*[1-x]^-1*[1-z]^-1
       - 7*ln(x)*ln([1-x])*[1-x]^-1
       - 7*ln(x)*ln([1-x])*[1-x]^-1*z
       - 8*ln(x)*ln([1-x])*[1-z]^-1
       + 10*ln(x)*ln([1-x])*z
       - 2*ln(x)*ln([1-x])*x*[1-z]^-1
       + 4*ln(x)*ln([1-x])*x
       + 6*ln(x)*ln(z)
       + 11*ln(x)*ln(z)*[1-x]^-1*[1-z]^-1
       - 8*ln(x)*ln(z)*[1-x]^-1
       - 8*ln(x)*ln(z)*[1-x]^-1*z
       - 17/2*ln(x)*ln(z)*[1-z]^-1
       + 11*ln(x)*ln(z)*z
       - 5/2*ln(x)*ln(z)*x*[1-z]^-1
       + 5*ln(x)*ln(z)*x
       + 3*ln(x)*ln([1-z])
       + 8*ln(x)*ln([1-z])*[1-x]^-1*[1-z]^-1
       - 13/2*ln(x)*ln([1-z])*[1-x]^-1
       - 13/2*ln(x)*ln([1-z])*[1-x]^-1*z
       - 11/2*ln(x)*ln([1-z])*[1-z]^-1
       + 8*ln(x)*ln([1-z])*z
       - 5/2*ln(x)*ln([1-z])*x*[1-z]^-1
       + 5*ln(x)*ln([1-z])*x
       - 4*ln([1-x])
       - 9*ln([1-x])*[1-x]^-1*[1-z]^-1
       + 4*ln([1-x])*[1-x]^-1
       + 7*ln([1-x])*[1-x]^-1*z
       + 9*ln([1-x])*[1-z]^-1
       - 9*ln([1-x])*z
       - 3*ln([1-x])*x
       - 2*ln([1-x])^2
       - 5/2*ln([1-x])^2*[1-x]^-1*[1-z]^-1
       + 3/2*ln([1-x])^2*[1-x]^-1
       + 3/2*ln([1-x])^2*[1-x]^-1*z
       + 9/4*ln([1-x])^2*[1-z]^-1
       - 5/2*ln([1-x])^2*z
       + 1/4*ln([1-x])^2*x*[1-z]^-1
       - 1/2*ln([1-x])^2*x
       - 5*ln([1-x])*ln(z)
       - 8*ln([1-x])*ln(z)*[1-x]^-1*[1-z]^-1
       + 11/2*ln([1-x])*ln(z)*[1-x]^-1
       + 11/2*ln([1-x])*ln(z)*[1-x]^-1*z
       + 13/2*ln([1-x])*ln(z)*[1-z]^-1
       - 8*ln([1-x])*ln(z)*z
       + 3/2*ln([1-x])*ln(z)*x*[1-z]^-1
       - 3*ln([1-x])*ln(z)*x
       - 2*ln([1-x])*ln([1-z])
       - 4*ln([1-x])*ln([1-z])*[1-x]^-1*[1-z]^-1
       + 3*ln([1-x])*ln([1-z])*[1-x]^-1
       + 3*ln([1-x])*ln([1-z])*[1-x]^-1*z
       + 3*ln([1-x])*ln([1-z])*[1-z]^-1
       - 4*ln([1-x])*ln([1-z])*z
       + ln([1-x])*ln([1-z])*x*[1-z]^-1
       - 2*ln([1-x])*ln([1-z])*x
       - 6*ln(z)
       - 12*ln(z)*[1-x]^-1*[1-z]^-1
       + 6*ln(z)*[1-x]^-1
       + 10*ln(z)*[1-x]^-1*z
       + 12*ln(z)*[1-z]^-1
       - 12*ln(z)*z
       - 4*ln(z)*x
       + ln(z)*ln( - [x-z])
       + 2*ln(z)*ln( - [x-z])*[1-x]^-1*[1-z]^-1
       - 3/2*ln(z)*ln( - [x-z])*[1-x]^-1
       - 3/2*ln(z)*ln( - [x-z])*[1-x]^-1*z
       - 3/2*ln(z)*ln( - [x-z])*[1-z]^-1
       + 2*ln(z)*ln( - [x-z])*z
       - 1/2*ln(z)*ln( - [x-z])*x*[1-z]^-1
       + ln(z)*ln( - [x-z])*x
       - 5/2*ln(z)^2
       - 4*ln(z)^2*[1-x]^-1*[1-z]^-1
       + 11/4*ln(z)^2*[1-x]^-1
       + 11/4*ln(z)^2*[1-x]^-1*z
       + 13/4*ln(z)^2*[1-z]^-1
       - 4*ln(z)^2*z
       + 3/4*ln(z)^2*x*[1-z]^-1
       - 3/2*ln(z)^2*x
       - ln(z)*ln([1-z])
       - 5*ln(z)*ln([1-z])*[1-x]^-1*[1-z]^-1
       + 9/2*ln(z)*ln([1-z])*[1-x]^-1
       + 9/2*ln(z)*ln([1-z])*[1-x]^-1*z
       + 3*ln(z)*ln([1-z])*[1-z]^-1
       - 5*ln(z)*ln([1-z])*z
       + 2*ln(z)*ln([1-z])*x*[1-z]^-1
       - 4*ln(z)*ln([1-z])*x
       - 5*ln([1-z])
       - 9*ln([1-z])*[1-x]^-1*[1-z]^-1
       + 5*ln([1-z])*[1-x]^-1
       + 8*ln([1-z])*[1-x]^-1*z
       + 9*ln([1-z])*[1-z]^-1
       - 9*ln([1-z])*z
       - 3*ln([1-z])*x
       + ln([1-z])*ln( - [1-x-z])*[1-x]^-1*[1-z]^-1
       - ln([1-z])*ln( - [1-x-z])*[1-x]^-1
       - ln([1-z])*ln( - [1-x-z])*[1-x]^-1*z
       - 1/2*ln([1-z])*ln( - [1-x-z])*[1-z]^-1
       + ln([1-z])*ln( - [1-x-z])*z
       - 1/2*ln([1-z])*ln( - [1-x-z])*x*[1-z]^-1
       + ln([1-z])*ln( - [1-x-z])*x
       - 1/2*ln([1-z])^2
       - 5/2*ln([1-z])^2*[1-x]^-1*[1-z]^-1
       + 9/4*ln([1-z])^2*[1-x]^-1
       + 9/4*ln([1-z])^2*[1-x]^-1*z
       + 3/2*ln([1-z])^2*[1-z]^-1
       - 5/2*ln([1-z])^2*z
       + ln([1-z])^2*x*[1-z]^-1
       - 2*ln([1-z])^2*x
       + Li2(x^-1*[1-x]^-1*z*[1-z])*[1-x]^-1*[1-z]^-1
       - Li2(x^-1*[1-x]^-1*z*[1-z])*[1-x]^-1
       - Li2(x^-1*[1-x]^-1*z*[1-z])*[1-x]^-1*z
       - 1/2*Li2(x^-1*[1-x]^-1*z*[1-z])*[1-z]^-1
       + Li2(x^-1*[1-x]^-1*z*[1-z])*z
       - 1/2*Li2(x^-1*[1-x]^-1*z*[1-z])*x*[1-z]^-1
       + Li2(x^-1*[1-x]^-1*z*[1-z])*x
       + Li2([1-x]^-1*[1-z])
       + 1/2*Li2([1-x]^-1*[1-z])*[1-x]^-1
       + 1/2*Li2([1-x]^-1*[1-z])*[1-x]^-1*z
       - 1/2*Li2([1-x]^-1*[1-z])*[1-z]^-1
       + 1/2*Li2([1-x]^-1*[1-z])*x*[1-z]^-1
       - Li2([1-x]^-1*[1-z])*x
       + Li2([1-x]*z^-1)*[1-x]^-1*[1-z]^-1
       - Li2([1-x]*z^-1)*[1-x]^-1
       - Li2([1-x]*z^-1)*[1-x]^-1*z
       - 1/2*Li2([1-x]*z^-1)*[1-z]^-1
       + Li2([1-x]*z^-1)*z
       - 1/2*Li2([1-x]*z^-1)*x*[1-z]^-1
       + Li2([1-x]*z^-1)*x
       - Li2(x*[1-x]^-1*z^-1*[1-z])
       - Li2(x*[1-x]^-1*z^-1*[1-z])*[1-x]^-1*[1-z]^-1
       + 1/2*Li2(x*[1-x]^-1*z^-1*[1-z])*[1-x]^-1
       + 1/2*Li2(x*[1-x]^-1*z^-1*[1-z])*[1-x]^-1*z
       + Li2(x*[1-x]^-1*z^-1*[1-z])*[1-z]^-1
       - Li2(x*[1-x]^-1*z^-1*[1-z])*z
       + Li2(z)
       + Li2(z)*[1-x]^-1*[1-z]^-1
       - 1/2*Li2(z)*[1-x]^-1
       - 1/2*Li2(z)*[1-x]^-1*z
       - Li2(z)*[1-z]^-1
       + Li2(z)*z
       )

       + T(u4) * ( 3
       + 4*[1-x]^-1*[1-z]^-1
       + 2*[1-x]^-1
       - 6*[1-x]^-1*z
       - 9*[1-z]^-1
       + 7*z
       + 5*x*[1-z]^-1
       - 5*x
       - 2/3*pi^2*[1-x]^-1*[1-z]^-1
       + 5/12*pi^2*[1-x]^-1
       + 5/12*pi^2*[1-x]^-1*z
       + 7/12*pi^2*[1-z]^-1
       - 1/2*pi^2
       - 2/3*pi^2*z
       + 1/12*pi^2*x*[1-z]^-1
       - 1/6*pi^2*x
       - 9*ln(x)
       - 18*ln(x)*[1-x]^-1*[1-z]^-1
       + 9*ln(x)*[1-x]^-1
       + 15*ln(x)*[1-x]^-1*z
       + 18*ln(x)*[1-z]^-1
       - 18*ln(x)*z
       - 6*ln(x)*x
       + ln(x)*ln( - [x-z])
       + 2*ln(x)*ln( - [x-z])*[1-x]^-1*[1-z]^-1
       - 3/2*ln(x)*ln( - [x-z])*[1-x]^-1
       - 3/2*ln(x)*ln( - [x-z])*[1-x]^-1*z
       - 3/2*ln(x)*ln( - [x-z])*[1-z]^-1
       + 2*ln(x)*ln( - [x-z])*z
       - 1/2*ln(x)*ln( - [x-z])*x*[1-z]^-1
       + ln(x)*ln( - [x-z])*x
       + ln(x)*ln( - [1-x-z])*[1-x]^-1*[1-z]^-1
       - ln(x)*ln( - [1-x-z])*[1-x]^-1
       - ln(x)*ln( - [1-x-z])*[1-x]^-1*z
       - 1/2*ln(x)*ln( - [1-x-z])*[1-z]^-1
       + ln(x)*ln( - [1-x-z])*z
       - 1/2*ln(x)*ln( - [1-x-z])*x*[1-z]^-1
       + ln(x)*ln( - [1-x-z])*x
       + 7/2*ln(x)^2
       + 13/2*ln(x)^2*[1-x]^-1*[1-z]^-1
       - 19/4*ln(x)^2*[1-x]^-1
       - 19/4*ln(x)^2*[1-x]^-1*z
       - 5*ln(x)^2*[1-z]^-1
       + 13/2*ln(x)^2*z
       - 3/2*ln(x)^2*x*[1-z]^-1
       + 3*ln(x)^2*x
       - 6*ln(x)*ln([1-x])
       - 10*ln(x)*ln([1-x])*[1-x]^-1*[1-z]^-1
       + 7*ln(x)*ln([1-x])*[1-x]^-1
       + 7*ln(x)*ln([1-x])*[1-x]^-1*z
       + 8*ln(x)*ln([1-x])*[1-z]^-1
       - 10*ln(x)*ln([1-x])*z
       + 2*ln(x)*ln([1-x])*x*[1-z]^-1
       - 4*ln(x)*ln([1-x])*x
       - 6*ln(x)*ln(z)
       - 11*ln(x)*ln(z)*[1-x]^-1*[1-z]^-1
       + 8*ln(x)*ln(z)*[1-x]^-1
       + 8*ln(x)*ln(z)*[1-x]^-1*z
       + 17/2*ln(x)*ln(z)*[1-z]^-1
       - 11*ln(x)*ln(z)*z
       + 5/2*ln(x)*ln(z)*x*[1-z]^-1
       - 5*ln(x)*ln(z)*x
       - 3*ln(x)*ln([1-z])
       - 8*ln(x)*ln([1-z])*[1-x]^-1*[1-z]^-1
       + 13/2*ln(x)*ln([1-z])*[1-x]^-1
       + 13/2*ln(x)*ln([1-z])*[1-x]^-1*z
       + 11/2*ln(x)*ln([1-z])*[1-z]^-1
       - 8*ln(x)*ln([1-z])*z
       + 5/2*ln(x)*ln([1-z])*x*[1-z]^-1
       - 5*ln(x)*ln([1-z])*x
       + 4*ln([1-x])
       + 9*ln([1-x])*[1-x]^-1*[1-z]^-1
       - 4*ln([1-x])*[1-x]^-1
       - 7*ln([1-x])*[1-x]^-1*z
       - 9*ln([1-x])*[1-z]^-1
       + 9*ln([1-x])*z
       + 3*ln([1-x])*x
       + 2*ln([1-x])^2
       + 5/2*ln([1-x])^2*[1-x]^-1*[1-z]^-1
       - 3/2*ln([1-x])^2*[1-x]^-1
       - 3/2*ln([1-x])^2*[1-x]^-1*z
       - 9/4*ln([1-x])^2*[1-z]^-1
       + 5/2*ln([1-x])^2*z
       - 1/4*ln([1-x])^2*x*[1-z]^-1
       + 1/2*ln([1-x])^2*x
       + 5*ln([1-x])*ln(z)
       + 8*ln([1-x])*ln(z)*[1-x]^-1*[1-z]^-1
       - 11/2*ln([1-x])*ln(z)*[1-x]^-1
       - 11/2*ln([1-x])*ln(z)*[1-x]^-1*z
       - 13/2*ln([1-x])*ln(z)*[1-z]^-1
       + 8*ln([1-x])*ln(z)*z
       - 3/2*ln([1-x])*ln(z)*x*[1-z]^-1
       + 3*ln([1-x])*ln(z)*x
       + 2*ln([1-x])*ln([1-z])
       + 4*ln([1-x])*ln([1-z])*[1-x]^-1*[1-z]^-1
       - 3*ln([1-x])*ln([1-z])*[1-x]^-1
       - 3*ln([1-x])*ln([1-z])*[1-x]^-1*z
       - 3*ln([1-x])*ln([1-z])*[1-z]^-1
       + 4*ln([1-x])*ln([1-z])*z
       - ln([1-x])*ln([1-z])*x*[1-z]^-1
       + 2*ln([1-x])*ln([1-z])*x
       + 6*ln(z)
       + 12*ln(z)*[1-x]^-1*[1-z]^-1
       - 6*ln(z)*[1-x]^-1
       - 10*ln(z)*[1-x]^-1*z
       - 12*ln(z)*[1-z]^-1
       + 12*ln(z)*z
       + 4*ln(z)*x
       - ln(z)*ln( - [x-z])
       - 2*ln(z)*ln( - [x-z])*[1-x]^-1*[1-z]^-1
       + 3/2*ln(z)*ln( - [x-z])*[1-x]^-1
       + 3/2*ln(z)*ln( - [x-z])*[1-x]^-1*z
       + 3/2*ln(z)*ln( - [x-z])*[1-z]^-1
       - 2*ln(z)*ln( - [x-z])*z
       + 1/2*ln(z)*ln( - [x-z])*x*[1-z]^-1
       - ln(z)*ln( - [x-z])*x
       + 5/2*ln(z)^2
       + 4*ln(z)^2*[1-x]^-1*[1-z]^-1
       - 11/4*ln(z)^2*[1-x]^-1
       - 11/4*ln(z)^2*[1-x]^-1*z
       - 13/4*ln(z)^2*[1-z]^-1
       + 4*ln(z)^2*z
       - 3/4*ln(z)^2*x*[1-z]^-1
       + 3/2*ln(z)^2*x
       + ln(z)*ln([1-z])
       + 5*ln(z)*ln([1-z])*[1-x]^-1*[1-z]^-1
       - 9/2*ln(z)*ln([1-z])*[1-x]^-1
       - 9/2*ln(z)*ln([1-z])*[1-x]^-1*z
       - 3*ln(z)*ln([1-z])*[1-z]^-1
       + 5*ln(z)*ln([1-z])*z
       - 2*ln(z)*ln([1-z])*x*[1-z]^-1
       + 4*ln(z)*ln([1-z])*x
       + 5*ln([1-z])
       + 9*ln([1-z])*[1-x]^-1*[1-z]^-1
       - 5*ln([1-z])*[1-x]^-1
       - 8*ln([1-z])*[1-x]^-1*z
       - 9*ln([1-z])*[1-z]^-1
       + 9*ln([1-z])*z
       + 3*ln([1-z])*x
       - ln([1-z])*ln( - [1-x-z])*[1-x]^-1*[1-z]^-1
       + ln([1-z])*ln( - [1-x-z])*[1-x]^-1
       + ln([1-z])*ln( - [1-x-z])*[1-x]^-1*z
       + 1/2*ln([1-z])*ln( - [1-x-z])*[1-z]^-1
       - ln([1-z])*ln( - [1-x-z])*z
       + 1/2*ln([1-z])*ln( - [1-x-z])*x*[1-z]^-1
       - ln([1-z])*ln( - [1-x-z])*x
       + 1/2*ln([1-z])^2
       + 5/2*ln([1-z])^2*[1-x]^-1*[1-z]^-1
       - 9/4*ln([1-z])^2*[1-x]^-1
       - 9/4*ln([1-z])^2*[1-x]^-1*z
       - 3/2*ln([1-z])^2*[1-z]^-1
       + 5/2*ln([1-z])^2*z
       - ln([1-z])^2*x*[1-z]^-1
       + 2*ln([1-z])^2*x
       - Li2(x^-1*[1-x]^-1*z*[1-z])*[1-x]^-1*[1-z]^-1
       + Li2(x^-1*[1-x]^-1*z*[1-z])*[1-x]^-1
       + Li2(x^-1*[1-x]^-1*z*[1-z])*[1-x]^-1*z
       + 1/2*Li2(x^-1*[1-x]^-1*z*[1-z])*[1-z]^-1
       - Li2(x^-1*[1-x]^-1*z*[1-z])*z
       + 1/2*Li2(x^-1*[1-x]^-1*z*[1-z])*x*[1-z]^-1
       - Li2(x^-1*[1-x]^-1*z*[1-z])*x
       - Li2([1-x]^-1*[1-z])
       - 1/2*Li2([1-x]^-1*[1-z])*[1-x]^-1
       - 1/2*Li2([1-x]^-1*[1-z])*[1-x]^-1*z
       + 1/2*Li2([1-x]^-1*[1-z])*[1-z]^-1
       - 1/2*Li2([1-x]^-1*[1-z])*x*[1-z]^-1
       + Li2([1-x]^-1*[1-z])*x
       - Li2([1-x]*z^-1)*[1-x]^-1*[1-z]^-1
       + Li2([1-x]*z^-1)*[1-x]^-1
       + Li2([1-x]*z^-1)*[1-x]^-1*z
       + 1/2*Li2([1-x]*z^-1)*[1-z]^-1
       - Li2([1-x]*z^-1)*z
       + 1/2*Li2([1-x]*z^-1)*x*[1-z]^-1
       - Li2([1-x]*z^-1)*x
       + Li2(x*[1-x]^-1*z^-1*[1-z])
       + Li2(x*[1-x]^-1*z^-1*[1-z])*[1-x]^-1*[1-z]^-1
       - 1/2*Li2(x*[1-x]^-1*z^-1*[1-z])*[1-x]^-1
       - 1/2*Li2(x*[1-x]^-1*z^-1*[1-z])*[1-x]^-1*z
       - Li2(x*[1-x]^-1*z^-1*[1-z])*[1-z]^-1
       + Li2(x*[1-x]^-1*z^-1*[1-z])*z
       - Li2(z)
       - Li2(z)*[1-x]^-1*[1-z]^-1
       + 1/2*Li2(z)*[1-x]^-1
       + 1/2*Li2(z)*[1-x]^-1*z
       + Li2(z)*[1-z]^-1
       - Li2(z)*z
       )

       + T(r1)*NC^2 * (  - 5/2
       - 2*[1-x]^-1*[1-z]^-1
       - 1/2*[1-x]^-1
       + 5/2*[1-x]^-1*z
       + 7/2*[1-z]^-1
       - 2*z
       - 3/2*x*[1-z]^-1
       + 3*x
       + 2/3*pi^2*[1-x]^-1*[1-z]^-1
       - 1/3*pi^2*[1-x]^-1
       - 1/3*pi^2*[1-x]^-1*z
       - 1/3*pi^2*[1-z]^-1
       + 2/3*pi^2*z
       - 1/3*pi^2*x*[1-z]^-1
       + 2/3*pi^2*x
       - 3/2*ln(x)
       + 9*ln(x)*[1-x]^-1*[1-z]^-1
       - 3*ln(x)*[1-x]^-1
       - 6*ln(x)*[1-x]^-1*z
       - 9*ln(x)*[1-z]^-1
       + 9*ln(x)*z
       + 9*ln(x)*x
       - 7/2*ln(x)^2*[1-x]^-1*[1-z]^-1
       + 7/4*ln(x)^2*[1-x]^-1
       + 7/4*ln(x)^2*[1-x]^-1*z
       + 7/4*ln(x)^2*[1-z]^-1
       - 7/2*ln(x)^2*z
       + 7/4*ln(x)^2*x*[1-z]^-1
       - 7/2*ln(x)^2*x
       + 6*ln(x)*ln([1-x])*[1-x]^-1*[1-z]^-1
       - 3*ln(x)*ln([1-x])*[1-x]^-1
       - 3*ln(x)*ln([1-x])*[1-x]^-1*z
       - 3*ln(x)*ln([1-x])*[1-z]^-1
       + 6*ln(x)*ln([1-x])*z
       - 3*ln(x)*ln([1-x])*x*[1-z]^-1
       + 6*ln(x)*ln([1-x])*x
       + 3*ln(x)*ln(z)*[1-x]^-1*[1-z]^-1
       - 3/2*ln(x)*ln(z)*[1-x]^-1
       - 3/2*ln(x)*ln(z)*[1-x]^-1*z
       - 3/2*ln(x)*ln(z)*[1-z]^-1
       + 3*ln(x)*ln(z)*z
       - 3/2*ln(x)*ln(z)*x*[1-z]^-1
       + 3*ln(x)*ln(z)*x
       + 6*ln(x)*ln([1-z])*[1-x]^-1*[1-z]^-1
       - 3*ln(x)*ln([1-z])*[1-x]^-1
       - 3*ln(x)*ln([1-z])*[1-x]^-1*z
       - 3*ln(x)*ln([1-z])*[1-z]^-1
       + 6*ln(x)*ln([1-z])*z
       - 3*ln(x)*ln([1-z])*x*[1-z]^-1
       + 6*ln(x)*ln([1-z])*x
       - ln(x)*ln([1-x-z])*[1-x]^-1*[1-z]^-1
       + 1/2*ln(x)*ln([1-x-z])*[1-x]^-1
       + 1/2*ln(x)*ln([1-x-z])*[1-x]^-1*z
       + 1/2*ln(x)*ln([1-x-z])*[1-z]^-1
       - ln(x)*ln([1-x-z])*z
       + 1/2*ln(x)*ln([1-x-z])*x*[1-z]^-1
       - ln(x)*ln([1-x-z])*x
       + ln([1-x])
       - 6*ln([1-x])*[1-x]^-1*[1-z]^-1
       + 2*ln([1-x])*[1-x]^-1
       + 4*ln([1-x])*[1-x]^-1*z
       + 6*ln([1-x])*[1-z]^-1
       - 6*ln([1-x])*z
       - 6*ln([1-x])*x
       - 2*ln([1-x])^2*[1-x]^-1*[1-z]^-1
       + ln([1-x])^2*[1-x]^-1
       + ln([1-x])^2*[1-x]^-1*z
       + ln([1-x])^2*[1-z]^-1
       - 2*ln([1-x])^2*z
       + ln([1-x])^2*x*[1-z]^-1
       - 2*ln([1-x])^2*x
       - 2*ln([1-x])*ln(z)*[1-x]^-1*[1-z]^-1
       + ln([1-x])*ln(z)*[1-x]^-1
       + ln([1-x])*ln(z)*[1-x]^-1*z
       + ln([1-x])*ln(z)*[1-z]^-1
       - 2*ln([1-x])*ln(z)*z
       + ln([1-x])*ln(z)*x*[1-z]^-1
       - 2*ln([1-x])*ln(z)*x
       - 5*ln([1-x])*ln([1-z])*[1-x]^-1*[1-z]^-1
       + 5/2*ln([1-x])*ln([1-z])*[1-x]^-1
       + 5/2*ln([1-x])*ln([1-z])*[1-x]^-1*z
       + 5/2*ln([1-x])*ln([1-z])*[1-z]^-1
       - 5*ln([1-x])*ln([1-z])*z
       + 5/2*ln([1-x])*ln([1-z])*x*[1-z]^-1
       - 5*ln([1-x])*ln([1-z])*x
       + 1/2*ln(z)
       - 3*ln(z)*[1-x]^-1*[1-z]^-1
       + ln(z)*[1-x]^-1
       + 2*ln(z)*[1-x]^-1*z
       + 3*ln(z)*[1-z]^-1
       - 3*ln(z)*z
       - 3*ln(z)*x
       - 1/2*ln(z)^2*[1-x]^-1*[1-z]^-1
       + 1/4*ln(z)^2*[1-x]^-1
       + 1/4*ln(z)^2*[1-x]^-1*z
       + 1/4*ln(z)^2*[1-z]^-1
       - 1/2*ln(z)^2*z
       + 1/4*ln(z)^2*x*[1-z]^-1
       - 1/2*ln(z)^2*x
       - 2*ln(z)*ln([1-z])*[1-x]^-1*[1-z]^-1
       + ln(z)*ln([1-z])*[1-x]^-1
       + ln(z)*ln([1-z])*[1-x]^-1*z
       + ln(z)*ln([1-z])*[1-z]^-1
       - 2*ln(z)*ln([1-z])*z
       + ln(z)*ln([1-z])*x*[1-z]^-1
       - 2*ln(z)*ln([1-z])*x
       + ln([1-z])
       - 6*ln([1-z])*[1-x]^-1*[1-z]^-1
       + 2*ln([1-z])*[1-x]^-1
       + 4*ln([1-z])*[1-x]^-1*z
       + 6*ln([1-z])*[1-z]^-1
       - 6*ln([1-z])*z
       - 6*ln([1-z])*x
       - 5/2*ln([1-z])^2*[1-x]^-1*[1-z]^-1
       + 5/4*ln([1-z])^2*[1-x]^-1
       + 5/4*ln([1-z])^2*[1-x]^-1*z
       + 5/4*ln([1-z])^2*[1-z]^-1
       - 5/2*ln([1-z])^2*z
       + 5/4*ln([1-z])^2*x*[1-z]^-1
       - 5/2*ln([1-z])^2*x
       + ln([1-z])*ln([1-x-z])*[1-x]^-1*[1-z]^-1
       - 1/2*ln([1-z])*ln([1-x-z])*[1-x]^-1
       - 1/2*ln([1-z])*ln([1-x-z])*[1-x]^-1*z
       - 1/2*ln([1-z])*ln([1-x-z])*[1-z]^-1
       + ln([1-z])*ln([1-x-z])*z
       - 1/2*ln([1-z])*ln([1-x-z])*x*[1-z]^-1
       + ln([1-z])*ln([1-x-z])*x
       + Li2([1-x]^-1*z)*[1-x]^-1*[1-z]^-1
       - 1/2*Li2([1-x]^-1*z)*[1-x]^-1
       - 1/2*Li2([1-x]^-1*z)*[1-x]^-1*z
       - 1/2*Li2([1-x]^-1*z)*[1-z]^-1
       + Li2([1-x]^-1*z)*z
       - 1/2*Li2([1-x]^-1*z)*x*[1-z]^-1
       + Li2([1-x]^-1*z)*x
       - Li2(x*[1-x]^-1*z*[1-z]^-1)*[1-x]^-1*[1-z]^-1
       + 1/2*Li2(x*[1-x]^-1*z*[1-z]^-1)*[1-x]^-1
       + 1/2*Li2(x*[1-x]^-1*z*[1-z]^-1)*[1-x]^-1*z
       + 1/2*Li2(x*[1-x]^-1*z*[1-z]^-1)*[1-z]^-1
       - Li2(x*[1-x]^-1*z*[1-z]^-1)*z
       + 1/2*Li2(x*[1-x]^-1*z*[1-z]^-1)*x*[1-z]^-1
       - Li2(x*[1-x]^-1*z*[1-z]^-1)*x
       - Li2(z)*[1-x]^-1*[1-z]^-1
       + 1/2*Li2(z)*[1-x]^-1
       + 1/2*Li2(z)*[1-x]^-1*z
       + 1/2*Li2(z)*[1-z]^-1
       - Li2(z)*z
       + 1/2*Li2(z)*x*[1-z]^-1
       - Li2(z)*x
       )

       + T(r1) * ( 5/2
       + 2*[1-x]^-1*[1-z]^-1
       + 1/2*[1-x]^-1
       - 5/2*[1-x]^-1*z
       - 7/2*[1-z]^-1
       + 2*z
       + 3/2*x*[1-z]^-1
       - 3*x
       - 2/3*pi^2*[1-x]^-1*[1-z]^-1
       + 1/3*pi^2*[1-x]^-1
       + 1/3*pi^2*[1-x]^-1*z
       + 1/3*pi^2*[1-z]^-1
       - 2/3*pi^2*z
       + 1/3*pi^2*x*[1-z]^-1
       - 2/3*pi^2*x
       + 3/2*ln(x)
       - 9*ln(x)*[1-x]^-1*[1-z]^-1
       + 3*ln(x)*[1-x]^-1
       + 6*ln(x)*[1-x]^-1*z
       + 9*ln(x)*[1-z]^-1
       - 9*ln(x)*z
       - 9*ln(x)*x
       + 7/2*ln(x)^2*[1-x]^-1*[1-z]^-1
       - 7/4*ln(x)^2*[1-x]^-1
       - 7/4*ln(x)^2*[1-x]^-1*z
       - 7/4*ln(x)^2*[1-z]^-1
       + 7/2*ln(x)^2*z
       - 7/4*ln(x)^2*x*[1-z]^-1
       + 7/2*ln(x)^2*x
       - 6*ln(x)*ln([1-x])*[1-x]^-1*[1-z]^-1
       + 3*ln(x)*ln([1-x])*[1-x]^-1
       + 3*ln(x)*ln([1-x])*[1-x]^-1*z
       + 3*ln(x)*ln([1-x])*[1-z]^-1
       - 6*ln(x)*ln([1-x])*z
       + 3*ln(x)*ln([1-x])*x*[1-z]^-1
       - 6*ln(x)*ln([1-x])*x
       - 3*ln(x)*ln(z)*[1-x]^-1*[1-z]^-1
       + 3/2*ln(x)*ln(z)*[1-x]^-1
       + 3/2*ln(x)*ln(z)*[1-x]^-1*z
       + 3/2*ln(x)*ln(z)*[1-z]^-1
       - 3*ln(x)*ln(z)*z
       + 3/2*ln(x)*ln(z)*x*[1-z]^-1
       - 3*ln(x)*ln(z)*x
       - 6*ln(x)*ln([1-z])*[1-x]^-1*[1-z]^-1
       + 3*ln(x)*ln([1-z])*[1-x]^-1
       + 3*ln(x)*ln([1-z])*[1-x]^-1*z
       + 3*ln(x)*ln([1-z])*[1-z]^-1
       - 6*ln(x)*ln([1-z])*z
       + 3*ln(x)*ln([1-z])*x*[1-z]^-1
       - 6*ln(x)*ln([1-z])*x
       + ln(x)*ln([1-x-z])*[1-x]^-1*[1-z]^-1
       - 1/2*ln(x)*ln([1-x-z])*[1-x]^-1
       - 1/2*ln(x)*ln([1-x-z])*[1-x]^-1*z
       - 1/2*ln(x)*ln([1-x-z])*[1-z]^-1
       + ln(x)*ln([1-x-z])*z
       - 1/2*ln(x)*ln([1-x-z])*x*[1-z]^-1
       + ln(x)*ln([1-x-z])*x
       - ln([1-x])
       + 6*ln([1-x])*[1-x]^-1*[1-z]^-1
       - 2*ln([1-x])*[1-x]^-1
       - 4*ln([1-x])*[1-x]^-1*z
       - 6*ln([1-x])*[1-z]^-1
       + 6*ln([1-x])*z
       + 6*ln([1-x])*x
       + 2*ln([1-x])^2*[1-x]^-1*[1-z]^-1
       - ln([1-x])^2*[1-x]^-1
       - ln([1-x])^2*[1-x]^-1*z
       - ln([1-x])^2*[1-z]^-1
       + 2*ln([1-x])^2*z
       - ln([1-x])^2*x*[1-z]^-1
       + 2*ln([1-x])^2*x
       + 2*ln([1-x])*ln(z)*[1-x]^-1*[1-z]^-1
       - ln([1-x])*ln(z)*[1-x]^-1
       - ln([1-x])*ln(z)*[1-x]^-1*z
       - ln([1-x])*ln(z)*[1-z]^-1
       + 2*ln([1-x])*ln(z)*z
       - ln([1-x])*ln(z)*x*[1-z]^-1
       + 2*ln([1-x])*ln(z)*x
       + 5*ln([1-x])*ln([1-z])*[1-x]^-1*[1-z]^-1
       - 5/2*ln([1-x])*ln([1-z])*[1-x]^-1
       - 5/2*ln([1-x])*ln([1-z])*[1-x]^-1*z
       - 5/2*ln([1-x])*ln([1-z])*[1-z]^-1
       + 5*ln([1-x])*ln([1-z])*z
       - 5/2*ln([1-x])*ln([1-z])*x*[1-z]^-1
       + 5*ln([1-x])*ln([1-z])*x
       - 1/2*ln(z)
       + 3*ln(z)*[1-x]^-1*[1-z]^-1
       - ln(z)*[1-x]^-1
       - 2*ln(z)*[1-x]^-1*z
       - 3*ln(z)*[1-z]^-1
       + 3*ln(z)*z
       + 3*ln(z)*x
       + 1/2*ln(z)^2*[1-x]^-1*[1-z]^-1
       - 1/4*ln(z)^2*[1-x]^-1
       - 1/4*ln(z)^2*[1-x]^-1*z
       - 1/4*ln(z)^2*[1-z]^-1
       + 1/2*ln(z)^2*z
       - 1/4*ln(z)^2*x*[1-z]^-1
       + 1/2*ln(z)^2*x
       + 2*ln(z)*ln([1-z])*[1-x]^-1*[1-z]^-1
       - ln(z)*ln([1-z])*[1-x]^-1
       - ln(z)*ln([1-z])*[1-x]^-1*z
       - ln(z)*ln([1-z])*[1-z]^-1
       + 2*ln(z)*ln([1-z])*z
       - ln(z)*ln([1-z])*x*[1-z]^-1
       + 2*ln(z)*ln([1-z])*x
       - ln([1-z])
       + 6*ln([1-z])*[1-x]^-1*[1-z]^-1
       - 2*ln([1-z])*[1-x]^-1
       - 4*ln([1-z])*[1-x]^-1*z
       - 6*ln([1-z])*[1-z]^-1
       + 6*ln([1-z])*z
       + 6*ln([1-z])*x
       + 5/2*ln([1-z])^2*[1-x]^-1*[1-z]^-1
       - 5/4*ln([1-z])^2*[1-x]^-1
       - 5/4*ln([1-z])^2*[1-x]^-1*z
       - 5/4*ln([1-z])^2*[1-z]^-1
       + 5/2*ln([1-z])^2*z
       - 5/4*ln([1-z])^2*x*[1-z]^-1
       + 5/2*ln([1-z])^2*x
       - ln([1-z])*ln([1-x-z])*[1-x]^-1*[1-z]^-1
       + 1/2*ln([1-z])*ln([1-x-z])*[1-x]^-1
       + 1/2*ln([1-z])*ln([1-x-z])*[1-x]^-1*z
       + 1/2*ln([1-z])*ln([1-x-z])*[1-z]^-1
       - ln([1-z])*ln([1-x-z])*z
       + 1/2*ln([1-z])*ln([1-x-z])*x*[1-z]^-1
       - ln([1-z])*ln([1-x-z])*x
       - Li2([1-x]^-1*z)*[1-x]^-1*[1-z]^-1
       + 1/2*Li2([1-x]^-1*z)*[1-x]^-1
       + 1/2*Li2([1-x]^-1*z)*[1-x]^-1*z
       + 1/2*Li2([1-x]^-1*z)*[1-z]^-1
       - Li2([1-x]^-1*z)*z
       + 1/2*Li2([1-x]^-1*z)*x*[1-z]^-1
       - Li2([1-x]^-1*z)*x
       + Li2(x*[1-x]^-1*z*[1-z]^-1)*[1-x]^-1*[1-z]^-1
       - 1/2*Li2(x*[1-x]^-1*z*[1-z]^-1)*[1-x]^-1
       - 1/2*Li2(x*[1-x]^-1*z*[1-z]^-1)*[1-x]^-1*z
       - 1/2*Li2(x*[1-x]^-1*z*[1-z]^-1)*[1-z]^-1
       + Li2(x*[1-x]^-1*z*[1-z]^-1)*z
       - 1/2*Li2(x*[1-x]^-1*z*[1-z]^-1)*x*[1-z]^-1
       + Li2(x*[1-x]^-1*z*[1-z]^-1)*x
       + Li2(z)*[1-x]^-1*[1-z]^-1
       - 1/2*Li2(z)*[1-x]^-1
       - 1/2*Li2(z)*[1-x]^-1*z
       - 1/2*Li2(z)*[1-z]^-1
       + Li2(z)*z
       - 1/2*Li2(z)*x*[1-z]^-1
       + Li2(z)*x
       )

       + T(r2)*NC^2 * (  - 5/2
       - 2*[1-x]^-1*[1-z]^-1
       - 1/2*[1-x]^-1
       + 5/2*[1-x]^-1*z
       + 7/2*[1-z]^-1
       - 2*z
       - 3/2*x*[1-z]^-1
       + 3*x
       + 2/3*pi^2*[1-x]^-1*[1-z]^-1
       - 1/3*pi^2*[1-x]^-1
       - 1/3*pi^2*[1-x]^-1*z
       - 1/3*pi^2*[1-z]^-1
       + 2/3*pi^2*z
       - 1/3*pi^2*x*[1-z]^-1
       + 2/3*pi^2*x
       - 3/2*ln(x)
       + 9*ln(x)*[1-x]^-1*[1-z]^-1
       - 3*ln(x)*[1-x]^-1
       - 6*ln(x)*[1-x]^-1*z
       - 9*ln(x)*[1-z]^-1
       + 9*ln(x)*z
       + 9*ln(x)*x
       - ln(x)*ln( - [1-x-z])*[1-x]^-1*[1-z]^-1
       + 1/2*ln(x)*ln( - [1-x-z])*[1-x]^-1
       + 1/2*ln(x)*ln( - [1-x-z])*[1-x]^-1*z
       + 1/2*ln(x)*ln( - [1-x-z])*[1-z]^-1
       - ln(x)*ln( - [1-x-z])*z
       + 1/2*ln(x)*ln( - [1-x-z])*x*[1-z]^-1
       - ln(x)*ln( - [1-x-z])*x
       - 3*ln(x)^2*[1-x]^-1*[1-z]^-1
       + 3/2*ln(x)^2*[1-x]^-1
       + 3/2*ln(x)^2*[1-x]^-1*z
       + 3/2*ln(x)^2*[1-z]^-1
       - 3*ln(x)^2*z
       + 3/2*ln(x)^2*x*[1-z]^-1
       - 3*ln(x)^2*x
       + 5*ln(x)*ln([1-x])*[1-x]^-1*[1-z]^-1
       - 5/2*ln(x)*ln([1-x])*[1-x]^-1
       - 5/2*ln(x)*ln([1-x])*[1-x]^-1*z
       - 5/2*ln(x)*ln([1-x])*[1-z]^-1
       + 5*ln(x)*ln([1-x])*z
       - 5/2*ln(x)*ln([1-x])*x*[1-z]^-1
       + 5*ln(x)*ln([1-x])*x
       + 4*ln(x)*ln(z)*[1-x]^-1*[1-z]^-1
       - 2*ln(x)*ln(z)*[1-x]^-1
       - 2*ln(x)*ln(z)*[1-x]^-1*z
       - 2*ln(x)*ln(z)*[1-z]^-1
       + 4*ln(x)*ln(z)*z
       - 2*ln(x)*ln(z)*x*[1-z]^-1
       + 4*ln(x)*ln(z)*x
       + 5*ln(x)*ln([1-z])*[1-x]^-1*[1-z]^-1
       - 5/2*ln(x)*ln([1-z])*[1-x]^-1
       - 5/2*ln(x)*ln([1-z])*[1-x]^-1*z
       - 5/2*ln(x)*ln([1-z])*[1-z]^-1
       + 5*ln(x)*ln([1-z])*z
       - 5/2*ln(x)*ln([1-z])*x*[1-z]^-1
       + 5*ln(x)*ln([1-z])*x
       + ln([1-x])
       - 6*ln([1-x])*[1-x]^-1*[1-z]^-1
       + 2*ln([1-x])*[1-x]^-1
       + 4*ln([1-x])*[1-x]^-1*z
       + 6*ln([1-x])*[1-z]^-1
       - 6*ln([1-x])*z
       - 6*ln([1-x])*x
       - 2*ln([1-x])^2*[1-x]^-1*[1-z]^-1
       + ln([1-x])^2*[1-x]^-1
       + ln([1-x])^2*[1-x]^-1*z
       + ln([1-x])^2*[1-z]^-1
       - 2*ln([1-x])^2*z
       + ln([1-x])^2*x*[1-z]^-1
       - 2*ln([1-x])^2*x
       - 2*ln([1-x])*ln(z)*[1-x]^-1*[1-z]^-1
       + ln([1-x])*ln(z)*[1-x]^-1
       + ln([1-x])*ln(z)*[1-x]^-1*z
       + ln([1-x])*ln(z)*[1-z]^-1
       - 2*ln([1-x])*ln(z)*z
       + ln([1-x])*ln(z)*x*[1-z]^-1
       - 2*ln([1-x])*ln(z)*x
       - 4*ln([1-x])*ln([1-z])*[1-x]^-1*[1-z]^-1
       + 2*ln([1-x])*ln([1-z])*[1-x]^-1
       + 2*ln([1-x])*ln([1-z])*[1-x]^-1*z
       + 2*ln([1-x])*ln([1-z])*[1-z]^-1
       - 4*ln([1-x])*ln([1-z])*z
       + 2*ln([1-x])*ln([1-z])*x*[1-z]^-1
       - 4*ln([1-x])*ln([1-z])*x
       + 1/2*ln(z)
       - 3*ln(z)*[1-x]^-1*[1-z]^-1
       + ln(z)*[1-x]^-1
       + 2*ln(z)*[1-x]^-1*z
       + 3*ln(z)*[1-z]^-1
       - 3*ln(z)*z
       - 3*ln(z)*x
       - 1/2*ln(z)^2*[1-x]^-1*[1-z]^-1
       + 1/4*ln(z)^2*[1-x]^-1
       + 1/4*ln(z)^2*[1-x]^-1*z
       + 1/4*ln(z)^2*[1-z]^-1
       - 1/2*ln(z)^2*z
       + 1/4*ln(z)^2*x*[1-z]^-1
       - 1/2*ln(z)^2*x
       - 3*ln(z)*ln([1-z])*[1-x]^-1*[1-z]^-1
       + 3/2*ln(z)*ln([1-z])*[1-x]^-1
       + 3/2*ln(z)*ln([1-z])*[1-x]^-1*z
       + 3/2*ln(z)*ln([1-z])*[1-z]^-1
       - 3*ln(z)*ln([1-z])*z
       + 3/2*ln(z)*ln([1-z])*x*[1-z]^-1
       - 3*ln(z)*ln([1-z])*x
       + ln([1-z])
       - 6*ln([1-z])*[1-x]^-1*[1-z]^-1
       + 2*ln([1-z])*[1-x]^-1
       + 4*ln([1-z])*[1-x]^-1*z
       + 6*ln([1-z])*[1-z]^-1
       - 6*ln([1-z])*z
       - 6*ln([1-z])*x
       + ln([1-z])*ln( - [1-x-z])*[1-x]^-1*[1-z]^-1
       - 1/2*ln([1-z])*ln( - [1-x-z])*[1-x]^-1
       - 1/2*ln([1-z])*ln( - [1-x-z])*[1-x]^-1*z
       - 1/2*ln([1-z])*ln( - [1-x-z])*[1-z]^-1
       + ln([1-z])*ln( - [1-x-z])*z
       - 1/2*ln([1-z])*ln( - [1-x-z])*x*[1-z]^-1
       + ln([1-z])*ln( - [1-x-z])*x
       - 2*ln([1-z])^2*[1-x]^-1*[1-z]^-1
       + ln([1-z])^2*[1-x]^-1
       + ln([1-z])^2*[1-x]^-1*z
       + ln([1-z])^2*[1-z]^-1
       - 2*ln([1-z])^2*z
       + ln([1-z])^2*x*[1-z]^-1
       - 2*ln([1-z])^2*x
       + Li2(x^-1*[1-x]*z^-1*[1-z])*[1-x]^-1*[1-z]^-1
       - 1/2*Li2(x^-1*[1-x]*z^-1*[1-z])*[1-x]^-1
       - 1/2*Li2(x^-1*[1-x]*z^-1*[1-z])*[1-x]^-1*z
       - 1/2*Li2(x^-1*[1-x]*z^-1*[1-z])*[1-z]^-1
       + Li2(x^-1*[1-x]*z^-1*[1-z])*z
       - 1/2*Li2(x^-1*[1-x]*z^-1*[1-z])*x*[1-z]^-1
       + Li2(x^-1*[1-x]*z^-1*[1-z])*x
       - Li2([1-x]*z^-1)*[1-x]^-1*[1-z]^-1
       + 1/2*Li2([1-x]*z^-1)*[1-x]^-1
       + 1/2*Li2([1-x]*z^-1)*[1-x]^-1*z
       + 1/2*Li2([1-x]*z^-1)*[1-z]^-1
       - Li2([1-x]*z^-1)*z
       + 1/2*Li2([1-x]*z^-1)*x*[1-z]^-1
       - Li2([1-x]*z^-1)*x
       - Li2(z)*[1-x]^-1*[1-z]^-1
       + 1/2*Li2(z)*[1-x]^-1
       + 1/2*Li2(z)*[1-x]^-1*z
       + 1/2*Li2(z)*[1-z]^-1
       - Li2(z)*z
       + 1/2*Li2(z)*x*[1-z]^-1
       - Li2(z)*x
       )

       + T(r2) * ( 5/2
       + 2*[1-x]^-1*[1-z]^-1
       + 1/2*[1-x]^-1
       - 5/2*[1-x]^-1*z
       - 7/2*[1-z]^-1
       + 2*z
       + 3/2*x*[1-z]^-1
       - 3*x
       - 2/3*pi^2*[1-x]^-1*[1-z]^-1
       + 1/3*pi^2*[1-x]^-1
       + 1/3*pi^2*[1-x]^-1*z
       + 1/3*pi^2*[1-z]^-1
       - 2/3*pi^2*z
       + 1/3*pi^2*x*[1-z]^-1
       - 2/3*pi^2*x
       + 3/2*ln(x)
       - 9*ln(x)*[1-x]^-1*[1-z]^-1
       + 3*ln(x)*[1-x]^-1
       + 6*ln(x)*[1-x]^-1*z
       + 9*ln(x)*[1-z]^-1
       - 9*ln(x)*z
       - 9*ln(x)*x
       + ln(x)*ln( - [1-x-z])*[1-x]^-1*[1-z]^-1
       - 1/2*ln(x)*ln( - [1-x-z])*[1-x]^-1
       - 1/2*ln(x)*ln( - [1-x-z])*[1-x]^-1*z
       - 1/2*ln(x)*ln( - [1-x-z])*[1-z]^-1
       + ln(x)*ln( - [1-x-z])*z
       - 1/2*ln(x)*ln( - [1-x-z])*x*[1-z]^-1
       + ln(x)*ln( - [1-x-z])*x
       + 3*ln(x)^2*[1-x]^-1*[1-z]^-1
       - 3/2*ln(x)^2*[1-x]^-1
       - 3/2*ln(x)^2*[1-x]^-1*z
       - 3/2*ln(x)^2*[1-z]^-1
       + 3*ln(x)^2*z
       - 3/2*ln(x)^2*x*[1-z]^-1
       + 3*ln(x)^2*x
       - 5*ln(x)*ln([1-x])*[1-x]^-1*[1-z]^-1
       + 5/2*ln(x)*ln([1-x])*[1-x]^-1
       + 5/2*ln(x)*ln([1-x])*[1-x]^-1*z
       + 5/2*ln(x)*ln([1-x])*[1-z]^-1
       - 5*ln(x)*ln([1-x])*z
       + 5/2*ln(x)*ln([1-x])*x*[1-z]^-1
       - 5*ln(x)*ln([1-x])*x
       - 4*ln(x)*ln(z)*[1-x]^-1*[1-z]^-1
       + 2*ln(x)*ln(z)*[1-x]^-1
       + 2*ln(x)*ln(z)*[1-x]^-1*z
       + 2*ln(x)*ln(z)*[1-z]^-1
       - 4*ln(x)*ln(z)*z
       + 2*ln(x)*ln(z)*x*[1-z]^-1
       - 4*ln(x)*ln(z)*x
       - 5*ln(x)*ln([1-z])*[1-x]^-1*[1-z]^-1
       + 5/2*ln(x)*ln([1-z])*[1-x]^-1
       + 5/2*ln(x)*ln([1-z])*[1-x]^-1*z
       + 5/2*ln(x)*ln([1-z])*[1-z]^-1
       - 5*ln(x)*ln([1-z])*z
       + 5/2*ln(x)*ln([1-z])*x*[1-z]^-1
       - 5*ln(x)*ln([1-z])*x
       - ln([1-x])
       + 6*ln([1-x])*[1-x]^-1*[1-z]^-1
       - 2*ln([1-x])*[1-x]^-1
       - 4*ln([1-x])*[1-x]^-1*z
       - 6*ln([1-x])*[1-z]^-1
       + 6*ln([1-x])*z
       + 6*ln([1-x])*x
       + 2*ln([1-x])^2*[1-x]^-1*[1-z]^-1
       - ln([1-x])^2*[1-x]^-1
       - ln([1-x])^2*[1-x]^-1*z
       - ln([1-x])^2*[1-z]^-1
       + 2*ln([1-x])^2*z
       - ln([1-x])^2*x*[1-z]^-1
       + 2*ln([1-x])^2*x
       + 2*ln([1-x])*ln(z)*[1-x]^-1*[1-z]^-1
       - ln([1-x])*ln(z)*[1-x]^-1
       - ln([1-x])*ln(z)*[1-x]^-1*z
       - ln([1-x])*ln(z)*[1-z]^-1
       + 2*ln([1-x])*ln(z)*z
       - ln([1-x])*ln(z)*x*[1-z]^-1
       + 2*ln([1-x])*ln(z)*x
       + 4*ln([1-x])*ln([1-z])*[1-x]^-1*[1-z]^-1
       - 2*ln([1-x])*ln([1-z])*[1-x]^-1
       - 2*ln([1-x])*ln([1-z])*[1-x]^-1*z
       - 2*ln([1-x])*ln([1-z])*[1-z]^-1
       + 4*ln([1-x])*ln([1-z])*z
       - 2*ln([1-x])*ln([1-z])*x*[1-z]^-1
       + 4*ln([1-x])*ln([1-z])*x
       - 1/2*ln(z)
       + 3*ln(z)*[1-x]^-1*[1-z]^-1
       - ln(z)*[1-x]^-1
       - 2*ln(z)*[1-x]^-1*z
       - 3*ln(z)*[1-z]^-1
       + 3*ln(z)*z
       + 3*ln(z)*x
       + 1/2*ln(z)^2*[1-x]^-1*[1-z]^-1
       - 1/4*ln(z)^2*[1-x]^-1
       - 1/4*ln(z)^2*[1-x]^-1*z
       - 1/4*ln(z)^2*[1-z]^-1
       + 1/2*ln(z)^2*z
       - 1/4*ln(z)^2*x*[1-z]^-1
       + 1/2*ln(z)^2*x
       + 3*ln(z)*ln([1-z])*[1-x]^-1*[1-z]^-1
       - 3/2*ln(z)*ln([1-z])*[1-x]^-1
       - 3/2*ln(z)*ln([1-z])*[1-x]^-1*z
       - 3/2*ln(z)*ln([1-z])*[1-z]^-1
       + 3*ln(z)*ln([1-z])*z
       - 3/2*ln(z)*ln([1-z])*x*[1-z]^-1
       + 3*ln(z)*ln([1-z])*x
       - ln([1-z])
       + 6*ln([1-z])*[1-x]^-1*[1-z]^-1
       - 2*ln([1-z])*[1-x]^-1
       - 4*ln([1-z])*[1-x]^-1*z
       - 6*ln([1-z])*[1-z]^-1
       + 6*ln([1-z])*z
       + 6*ln([1-z])*x
       - ln([1-z])*ln( - [1-x-z])*[1-x]^-1*[1-z]^-1
       + 1/2*ln([1-z])*ln( - [1-x-z])*[1-x]^-1
       + 1/2*ln([1-z])*ln( - [1-x-z])*[1-x]^-1*z
       + 1/2*ln([1-z])*ln( - [1-x-z])*[1-z]^-1
       - ln([1-z])*ln( - [1-x-z])*z
       + 1/2*ln([1-z])*ln( - [1-x-z])*x*[1-z]^-1
       - ln([1-z])*ln( - [1-x-z])*x
       + 2*ln([1-z])^2*[1-x]^-1*[1-z]^-1
       - ln([1-z])^2*[1-x]^-1
       - ln([1-z])^2*[1-x]^-1*z
       - ln([1-z])^2*[1-z]^-1
       + 2*ln([1-z])^2*z
       - ln([1-z])^2*x*[1-z]^-1
       + 2*ln([1-z])^2*x
       - Li2(x^-1*[1-x]*z^-1*[1-z])*[1-x]^-1*[1-z]^-1
       + 1/2*Li2(x^-1*[1-x]*z^-1*[1-z])*[1-x]^-1
       + 1/2*Li2(x^-1*[1-x]*z^-1*[1-z])*[1-x]^-1*z
       + 1/2*Li2(x^-1*[1-x]*z^-1*[1-z])*[1-z]^-1
       - Li2(x^-1*[1-x]*z^-1*[1-z])*z
       + 1/2*Li2(x^-1*[1-x]*z^-1*[1-z])*x*[1-z]^-1
       - Li2(x^-1*[1-x]*z^-1*[1-z])*x
       + Li2([1-x]*z^-1)*[1-x]^-1*[1-z]^-1
       - 1/2*Li2([1-x]*z^-1)*[1-x]^-1
       - 1/2*Li2([1-x]*z^-1)*[1-x]^-1*z
       - 1/2*Li2([1-x]*z^-1)*[1-z]^-1
       + Li2([1-x]*z^-1)*z
       - 1/2*Li2([1-x]*z^-1)*x*[1-z]^-1
       + Li2([1-x]*z^-1)*x
       + Li2(z)*[1-x]^-1*[1-z]^-1
       - 1/2*Li2(z)*[1-x]^-1
       - 1/2*Li2(z)*[1-x]^-1*z
       - 1/2*Li2(z)*[1-z]^-1
       + Li2(z)*z
       - 1/2*Li2(z)*x*[1-z]^-1
       + Li2(z)*x
       )

       - 32/3
       - 6*[1-x]^-1*z^-1*[1-z]^-1
       + 6*[1-x]^-1*z^-1
       - 2*[1-x]^-1*z^-1*rln2^2
       - 3*[1-x]^-1*z^-1*sqrtxz1*rln2
       - 5/2*[1-x]^-1
       - 4*[1-x]^-1*rln2^2
       - 3*[1-x]^-1*sqrtxz1*rln2
       + 17/2*[1-x]^-1*z
       - 2*[1-x]^-1*z*rln2^2
       + 4*z^-1*[1-z]^-1
       - z^-1
       + z^-1*rln2^2
       + z^-1*sqrtxz1*rln2
       + 17/2*[1-z]^-1
       - 1/2*poly2^-1
       + 2*rln2^2
       + 3*sqrtxz1*rln2
       - 1/2*z*[1-x-z]^-1
       - 103/18*z
       + 2*z*rln2^2
       + 2*x*z^-1*[1-z]^-1
       - 19/4*x*z^-1
       + x*z^-1*rln2^2
       + 2*x*z^-1*sqrtxz1*rln2
       - 17/2*x*[1-z]^-1
       + 1/2*x*poly2^-1
       + 148/9*x
       + 2*x*rln2^2
       + 5/4*x*z
       - 1/2*x^2*poly2^-1
       + 1/2*x^3*poly2^-1
       - 1/6*pi^2*x^-2*[1+x]^-1*z^-1
       + 1/3*pi^2*x^-2*[1+x]^-1
       + 1/6*pi^2*x^-2*z^-1
       - 1/3*pi^2*x^-2
       + 1/3*pi^2*x^-1*[1+x]^-1*z
       - 1/6*pi^2*x^-1*z^-1
       + 1/3*pi^2*x^-1
       - 1/3*pi^2*x^-1*z
       + 2/3*pi^2*[1-x]^-1*z^-1*[1-z]^-1
       - 2/3*pi^2*[1-x]^-1*z^-1
       + 2/3*pi^2*[1-x]^-1*[1-z]^-1
       - 11/12*pi^2*[1-x]^-1
       - 11/12*pi^2*[1-x]^-1*z
       + 1/3*pi^2*[1+x]^-1*z^-1
       - 2/3*pi^2*[1+x]^-1
       + 1/2*pi^2*[1+x]^-1*z
       - 1/3*pi^2*z^-1*[1-z]^-1
       + 1/6*pi^2*z^-1
       - 7/12*pi^2*[1-z]^-1
       + pi^2
       + 13/6*pi^2*z
       - 1/3*pi^2*x*z^-1*[1-z]^-1
       + 1/3*pi^2*x*z^-1
       - 1/12*pi^2*x*[1-z]^-1
       + 11/6*pi^2*x
       + 1/6*pi^2*x*z
       + 2*ln(1 + sqrtxz1 - z)*[1-x]^-1*z^-1*rln2
       + 3*ln(1 + sqrtxz1 - z)*[1-x]^-1*z^-1*sqrtxz1
       + 4*ln(1 + sqrtxz1 - z)*[1-x]^-1*rln2
       + 3*ln(1 + sqrtxz1 - z)*[1-x]^-1*sqrtxz1
       + 2*ln(1 + sqrtxz1 - z)*[1-x]^-1*z*rln2
       - ln(1 + sqrtxz1 - z)*z^-1*rln2
       - ln(1 + sqrtxz1 - z)*z^-1*sqrtxz1
       - 2*ln(1 + sqrtxz1 - z)*rln2
       - 3*ln(1 + sqrtxz1 - z)*sqrtxz1
       - 2*ln(1 + sqrtxz1 - z)*z*rln2
       - ln(1 + sqrtxz1 - z)*x*z^-1*rln2
       - 2*ln(1 + sqrtxz1 - z)*x*z^-1*sqrtxz1
       - 2*ln(1 + sqrtxz1 - z)*x*rln2
       + 2*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)
       - 2*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^-1
       - 4*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*[1-x]^-1
       - 2*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z
       + ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*z^-1
       + 2*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*z
       + ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*x*z^-1
       + 2*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*x
       + 2*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^-1*rln2
       + 4*ln(1 + sqrtxz1 + z)*[1-x]^-1*rln2
       + 2*ln(1 + sqrtxz1 + z)*[1-x]^-1*z*rln2
       - ln(1 + sqrtxz1 + z)*z^-1*rln2
       - 2*ln(1 + sqrtxz1 + z)*rln2
       - 2*ln(1 + sqrtxz1 + z)*z*rln2
       - ln(1 + sqrtxz1 + z)*x*z^-1*rln2
       - 2*ln(1 + sqrtxz1 + z)*x*rln2
       - 3*ln(z*sqrtxz3)*ArcTan(z*sqrtxz3)*x^-1*sqrtxz3
       + ln(z*sqrtxz3)*ArcTan(z*sqrtxz3)*z^-1*sqrtxz3
       - 3*ln(z*sqrtxz3)*ArcTan(z*sqrtxz3)*z*sqrtxz3
       - 3*ln(z*sqrtxz3)*ArcTan(z*sqrtxz3)*x*sqrtxz3
       + 31/4*ln(x)
       - 2*ln(x)*x^-2*[1+x]^-1*z^-1
       + 4*ln(x)*x^-2*[1+x]^-1
       + 2*ln(x)*x^-2*z^-1
       - 4*ln(x)*x^-2
       + 4*ln(x)*x^-1*[1+x]^-1*z
       - 2*ln(x)*x^-1*z^-1
       + 4*ln(x)*x^-1
       - 4*ln(x)*x^-1*z
       + 24*ln(x)*[1-x]^-1*z^-1*[1-z]^-1
       - 87/4*ln(x)*[1-x]^-1*z^-1
       - ln(x)*[1-x]^-1*z^-1*rln2
       - 3/2*ln(x)*[1-x]^-1*z^-1*sqrtxz1
       + 3*ln(x)*[1-x]^-1*[1-z]^-1
       + 3/2*ln(x)*[1-x]^-1*[x-z]^-1
       - 47/4*ln(x)*[1-x]^-1
       - 2*ln(x)*[1-x]^-1*rln2
       - 3/2*ln(x)*[1-x]^-1*sqrtxz1
       - 45/2*ln(x)*[1-x]^-1*z
       - ln(x)*[1-x]^-1*z*rln2
       + 2*ln(x)*[1+x]^-1*z^-1
       - 4*ln(x)*[1+x]^-1
       + 4*ln(x)*[1+x]^-1*z
       - 16*ln(x)*z^-1*[1-z]^-1
       + 33/2*ln(x)*z^-1
       + 1/2*ln(x)*z^-1*rln2
       + 1/2*ln(x)*z^-1*sqrtxz1
       - 11*ln(x)*[1-z]^-1
       + ln(x)*[1-x-z]^-1
       - 3/2*ln(x)*[x-z]^-1
       - 3/4*ln(x)*poly2^-2
       + 1/2*ln(x)*poly2^-1
       + ln(x)*rln2
       + 3/2*ln(x)*sqrtxz1
       - 1/2*ln(x)*z*[1-x-z]^-2
       + ln(x)*z*[1-x-z]^-1
       + 76/3*ln(x)*z
       + ln(x)*z*rln2
       + 1/2*ln(x)*z^2*[1-x-z]^-2
       - 2*ln(x)*[1+x]
       + 2*ln(x)*[1+x]*z
       - 8*ln(x)*x*z^-1*[1-z]^-1
       + 8*ln(x)*x*z^-1
       + 1/2*ln(x)*x*z^-1*rln2
       + ln(x)*x*z^-1*sqrtxz1
       + 8*ln(x)*x*[1-z]^-1
       - 3/2*ln(x)*x*[x-z]^-1
       + 3/4*ln(x)*x*poly2^-2
       + 145/12*ln(x)*x
       + ln(x)*x*rln2
       - 3/2*ln(x)*x*z
       - 2*ln(x)*x^2*[x-z]^-1
       + 1/2*ln(x)*x^3*poly2^-1
       + 3/4*ln(x)*x^4*poly2^-2
       - 3/4*ln(x)*x^5*poly2^-2
       + 1/2*ln(x)*ln(1 - sqrtxz2 + x)*z^-1*sqrtxz2^-1
       + 1/2*ln(x)*ln(1 - sqrtxz2 + x)*[1-z]^-1*sqrtxz2^-1
       - 3/8*ln(x)*ln(1 - sqrtxz2 + x)*sqrtxz2^-1*poly2^-2
       + 3/8*ln(x)*ln(1 - sqrtxz2 + x)*sqrtxz2^-1*poly2^-1
       - 2*ln(x)*ln(1 - sqrtxz2 + x)*sqrtxz2^-1
       - ln(x)*ln(1 - sqrtxz2 + x)*x*z^-1*sqrtxz2^-1
       + ln(x)*ln(1 - sqrtxz2 + x)*x*[1-z]^-1*sqrtxz2^-1
       - 3/4*ln(x)*ln(1 - sqrtxz2 + x)*x*sqrtxz2^-1
       + 3/2*ln(x)*ln(1 - sqrtxz2 + x)*x*z*sqrtxz2^-1
       + 1/2*ln(x)*ln(1 - sqrtxz2 + x)*x^2*z^-1*sqrtxz2^-1
       + 1/2*ln(x)*ln(1 - sqrtxz2 + x)*x^2*[1-z]^-1*sqrtxz2^-1
       + 3/8*ln(x)*ln(1 - sqrtxz2 + x)*x^2*sqrtxz2^-1*poly2^-2
       + 1/4*ln(x)*ln(1 - sqrtxz2 + x)*x^2*sqrtxz2^-1*poly2^-1
       - 2*ln(x)*ln(1 - sqrtxz2 + x)*x^2*sqrtxz2^-1
       + 3/8*ln(x)*ln(1 - sqrtxz2 + x)*x^4*sqrtxz2^-1*poly2^-2
       + 3/8*ln(x)*ln(1 - sqrtxz2 + x)*x^4*sqrtxz2^-1*poly2^-1
       - 3/8*ln(x)*ln(1 - sqrtxz2 + x)*x^6*sqrtxz2^-1*poly2^-2
       - 1/2*ln(x)*ln(1 + sqrtxz2 + x)*z^-1*sqrtxz2^-1
       - 1/2*ln(x)*ln(1 + sqrtxz2 + x)*[1-z]^-1*sqrtxz2^-1
       + 3/8*ln(x)*ln(1 + sqrtxz2 + x)*sqrtxz2^-1*poly2^-2
       - 3/8*ln(x)*ln(1 + sqrtxz2 + x)*sqrtxz2^-1*poly2^-1
       + 2*ln(x)*ln(1 + sqrtxz2 + x)*sqrtxz2^-1
       + ln(x)*ln(1 + sqrtxz2 + x)*x*z^-1*sqrtxz2^-1
       - ln(x)*ln(1 + sqrtxz2 + x)*x*[1-z]^-1*sqrtxz2^-1
       + 3/4*ln(x)*ln(1 + sqrtxz2 + x)*x*sqrtxz2^-1
       - 3/2*ln(x)*ln(1 + sqrtxz2 + x)*x*z*sqrtxz2^-1
       - 1/2*ln(x)*ln(1 + sqrtxz2 + x)*x^2*z^-1*sqrtxz2^-1
       - 1/2*ln(x)*ln(1 + sqrtxz2 + x)*x^2*[1-z]^-1*sqrtxz2^-1
       - 3/8*ln(x)*ln(1 + sqrtxz2 + x)*x^2*sqrtxz2^-1*poly2^-2
       - 1/4*ln(x)*ln(1 + sqrtxz2 + x)*x^2*sqrtxz2^-1*poly2^-1
       + 2*ln(x)*ln(1 + sqrtxz2 + x)*x^2*sqrtxz2^-1
       - 3/8*ln(x)*ln(1 + sqrtxz2 + x)*x^4*sqrtxz2^-1*poly2^-2
       - 3/8*ln(x)*ln(1 + sqrtxz2 + x)*x^4*sqrtxz2^-1*poly2^-1
       + 3/8*ln(x)*ln(1 + sqrtxz2 + x)*x^6*sqrtxz2^-1*poly2^-2
       - ln(x)*ln(1 + sqrtxz1 + z)
       + ln(x)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^-1
       + 2*ln(x)*ln(1 + sqrtxz1 + z)*[1-x]^-1
       + ln(x)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z
       - 1/2*ln(x)*ln(1 + sqrtxz1 + z)*z^-1
       - ln(x)*ln(1 + sqrtxz1 + z)*z
       - 1/2*ln(x)*ln(1 + sqrtxz1 + z)*x*z^-1
       - ln(x)*ln(1 + sqrtxz1 + z)*x
       + ln(x)*ln(1 + x*z^-1)
       + 1/2*ln(x)*ln(1 + x*z^-1)*x^-2*[1+x]^-1*z^-1
       - ln(x)*ln(1 + x*z^-1)*x^-2*[1+x]^-1
       - 1/2*ln(x)*ln(1 + x*z^-1)*x^-2*z^-1
       + ln(x)*ln(1 + x*z^-1)*x^-2
       - ln(x)*ln(1 + x*z^-1)*x^-1*[1+x]^-1*z
       + 1/2*ln(x)*ln(1 + x*z^-1)*x^-1*z^-1
       - ln(x)*ln(1 + x*z^-1)*x^-1
       + ln(x)*ln(1 + x*z^-1)*x^-1*z
       + 1/2*ln(x)*ln(1 + x*z^-1)*[1+x]^-1*z^-1
       - ln(x)*ln(1 + x*z^-1)*[1+x]^-1
       - 1/2*ln(x)*ln(1 + x*z^-1)*z^-1
       - ln(x)*ln(1 + x*z^-1)*z
       + 1/2*ln(x)*ln(1 + x*z^-1)*x*z^-1
       - ln(x)*ln(1 + x*z^-1)*x
       + 4*ln(x)*ln(1 + x)
       + 2*ln(x)*ln(1 + x)*[1+x]^-1*z^-1
       - 4*ln(x)*ln(1 + x)*[1+x]^-1
       + 2*ln(x)*ln(1 + x)*[1+x]^-1*z
       - 2*ln(x)*ln(1 + x)*z^-1
       - 2*ln(x)*ln(1 + x)*z
       + 2*ln(x)*ln(1 + x*z)
       + 1/2*ln(x)*ln(1 + x*z)*x^-2*[1+x]^-1*z^-1
       - ln(x)*ln(1 + x*z)*x^-2*[1+x]^-1
       - 1/2*ln(x)*ln(1 + x*z)*x^-2*z^-1
       + ln(x)*ln(1 + x*z)*x^-2
       - ln(x)*ln(1 + x*z)*x^-1*[1+x]^-1*z
       + 1/2*ln(x)*ln(1 + x*z)*x^-1*z^-1
       - ln(x)*ln(1 + x*z)*x^-1
       + ln(x)*ln(1 + x*z)*x^-1*z
       - ln(x)*ln(1 + x*z)*[1-x]^-1*z^-1
       - 2*ln(x)*ln(1 + x*z)*[1-x]^-1
       - ln(x)*ln(1 + x*z)*[1-x]^-1*z
       + 1/2*ln(x)*ln(1 + x*z)*[1+x]^-1*z^-1
       - ln(x)*ln(1 + x*z)*[1+x]^-1
       + ln(x)*ln(1 + x*z)*x*z^-1
       - ln(x)*ln(z + x)
       + ln(x)*ln(z + x)*[1-x]^-1*z^-1
       + 2*ln(x)*ln(z + x)*[1-x]^-1
       + ln(x)*ln(z + x)*[1-x]^-1*z
       - 1/2*ln(x)*ln(z + x)*z^-1
       - ln(x)*ln(z + x)*z
       - 1/2*ln(x)*ln(z + x)*x*z^-1
       - ln(x)*ln(z + x)*x
       - 11/2*ln(x)^2
       + 3/4*ln(x)^2*x^-2*[1+x]^-1*z^-1
       - 3/2*ln(x)^2*x^-2*[1+x]^-1
       - 3/4*ln(x)^2*x^-2*z^-1
       + 3/2*ln(x)^2*x^-2
       - 3/2*ln(x)^2*x^-1*[1+x]^-1*z
       + 3/4*ln(x)^2*x^-1*z^-1
       - 3/2*ln(x)^2*x^-1
       + 3/2*ln(x)^2*x^-1*z
       - 8*ln(x)^2*[1-x]^-1*z^-1*[1-z]^-1
       + 8*ln(x)^2*[1-x]^-1*z^-1
       - 9/2*ln(x)^2*[1-x]^-1*[1-z]^-1
       + 9*ln(x)^2*[1-x]^-1
       + 9*ln(x)^2*[1-x]^-1*z
       - 5/4*ln(x)^2*[1+x]^-1*z^-1
       + 5/2*ln(x)^2*[1+x]^-1
       - 2*ln(x)^2*[1+x]^-1*z
       + 11/2*ln(x)^2*z^-1*[1-z]^-1
       - 5*ln(x)^2*z^-1
       + 11/4*ln(x)^2*[1-z]^-1
       - 12*ln(x)^2*z
       + 5/2*ln(x)^2*x*z^-1*[1-z]^-1
       - 5/2*ln(x)^2*x*z^-1
       + 7/4*ln(x)^2*x*[1-z]^-1
       - 17/2*ln(x)^2*x
       - 1/2*ln(x)^2*x*z
       + 8*ln(x)*ln([1-x])
       + ln(x)*ln([1-x])*x^-2*[1+x]^-1*z^-1
       - 2*ln(x)*ln([1-x])*x^-2*[1+x]^-1
       - ln(x)*ln([1-x])*x^-2*z^-1
       + 2*ln(x)*ln([1-x])*x^-2
       - 2*ln(x)*ln([1-x])*x^-1*[1+x]^-1*z
       + ln(x)*ln([1-x])*x^-1*z^-1
       - 2*ln(x)*ln([1-x])*x^-1
       + 2*ln(x)*ln([1-x])*x^-1*z
       - 4*ln(x)*ln([1-x])*[1-x]^-1*z^-1*[1-z]^-1
       + 4*ln(x)*ln([1-x])*[1-x]^-1*z^-1
       + 19*ln(x)*ln([1-x])*[1-x]^-1*[1-z]^-1
       - 25/2*ln(x)*ln([1-x])*[1-x]^-1
       - 25/2*ln(x)*ln([1-x])*[1-x]^-1*z
       - ln(x)*ln([1-x])*[1+x]^-1*z^-1
       + 2*ln(x)*ln([1-x])*[1+x]^-1
       - 2*ln(x)*ln([1-x])*[1+x]^-1*z
       + 2*ln(x)*ln([1-x])*z^-1*[1-z]^-1
       - 2*ln(x)*ln([1-x])*z^-1
       - 13*ln(x)*ln([1-x])*[1-z]^-1
       + 18*ln(x)*ln([1-x])*z
       + 2*ln(x)*ln([1-x])*x*z^-1*[1-z]^-1
       - 2*ln(x)*ln([1-x])*x*z^-1
       - 6*ln(x)*ln([1-x])*x*[1-z]^-1
       + 11*ln(x)*ln([1-x])*x
       + ln(x)*ln([1-x])*x*z
       - 2*ln(x)*ln([1+x])
       - ln(x)*ln([1+x])*x^-2*[1+x]^-1*z^-1
       + 2*ln(x)*ln([1+x])*x^-2*[1+x]^-1
       + ln(x)*ln([1+x])*x^-2*z^-1
       - 2*ln(x)*ln([1+x])*x^-2
       + 2*ln(x)*ln([1+x])*x^-1*[1+x]^-1*z
       - ln(x)*ln([1+x])*x^-1*z^-1
       + 2*ln(x)*ln([1+x])*x^-1
       - 2*ln(x)*ln([1+x])*x^-1*z
       - ln(x)*ln([1+x])*[1+x]^-1*z^-1
       + 2*ln(x)*ln([1+x])*[1+x]^-1
       + ln(x)*ln([1+x])*z^-1
       + 2*ln(x)*ln([1+x])*z
       - ln(x)*ln([1+x])*x*z^-1
       + 2*ln(x)*ln([1+x])*x
       + 2*ln(x)*ln(z)
       - 1/2*ln(x)*ln(z)*x^-2*[1+x]^-1*z^-1
       + ln(x)*ln(z)*x^-2*[1+x]^-1
       + 1/2*ln(x)*ln(z)*x^-2*z^-1
       - ln(x)*ln(z)*x^-2
       + ln(x)*ln(z)*x^-1*[1+x]^-1*z
       - 1/2*ln(x)*ln(z)*x^-1*z^-1
       + ln(x)*ln(z)*x^-1
       - ln(x)*ln(z)*x^-1*z
       + 8*ln(x)*ln(z)*[1-x]^-1*z^-1*[1-z]^-1
       - 8*ln(x)*ln(z)*[1-x]^-1*z^-1
       + 6*ln(x)*ln(z)*[1-x]^-1*[1-z]^-1
       - 15/2*ln(x)*ln(z)*[1-x]^-1
       - 17/2*ln(x)*ln(z)*[1-x]^-1*z
       - 1/2*ln(x)*ln(z)*[1+x]^-1*z^-1
       + ln(x)*ln(z)*[1+x]^-1
       - 5*ln(x)*ln(z)*z^-1*[1-z]^-1
       + 6*ln(x)*ln(z)*z^-1
       - 4*ln(x)*ln(z)*[1-z]^-1
       + 12*ln(x)*ln(z)*z
       - 3*ln(x)*ln(z)*x*z^-1*[1-z]^-1
       + 3*ln(x)*ln(z)*x*z^-1
       - 2*ln(x)*ln(z)*x*[1-z]^-1
       + 7*ln(x)*ln(z)*x
       + 4*ln(x)*ln([1-z])
       + 15*ln(x)*ln([1-z])*[1-x]^-1*[1-z]^-1
       - 12*ln(x)*ln([1-z])*[1-x]^-1
       - 12*ln(x)*ln([1-z])*[1-x]^-1*z
       - 9*ln(x)*ln([1-z])*[1-z]^-1
       + 17*ln(x)*ln([1-z])*z
       - 6*ln(x)*ln([1-z])*x*[1-z]^-1
       + 14*ln(x)*ln([1-z])*x
       + ln(x)*ln([1-z])*x*z
       - 4*ln([1-x])
       - 15*ln([1-x])*[1-x]^-1*[1-z]^-1
       + 6*ln([1-x])*[1-x]^-1
       + 11*ln([1-x])*[1-x]^-1*z
       + 15*ln([1-x])*[1-z]^-1
       - 38/3*ln([1-x])*z
       - 37/6*ln([1-x])*x
       - 1/2*ln([1-x])*x*z
       - 5/2*ln([1-x])^2
       - 9/2*ln([1-x])^2*[1-x]^-1*[1-z]^-1
       + 5/2*ln([1-x])^2*[1-x]^-1
       + 5/2*ln([1-x])^2*[1-x]^-1*z
       + 13/4*ln([1-x])^2*[1-z]^-1
       - 6*ln([1-x])^2*z
       + 5/4*ln([1-x])^2*x*[1-z]^-1
       - 4*ln([1-x])^2*x
       - 1/2*ln([1-x])^2*x*z
       - 3*ln([1-x])*ln(z)
       - 8*ln([1-x])*ln(z)*[1-x]^-1*[1-z]^-1
       + 5*ln([1-x])*ln(z)*[1-x]^-1
       + 5*ln([1-x])*ln(z)*[1-x]^-1*z
       + 11/2*ln([1-x])*ln(z)*[1-z]^-1
       - 7*ln([1-x])*ln(z)*z
       + 5/2*ln([1-x])*ln(z)*x*[1-z]^-1
       - 4*ln([1-x])*ln(z)*x
       - 3*ln([1-x])*ln([1-z])
       - 8*ln([1-x])*ln([1-z])*[1-x]^-1*[1-z]^-1
       + 5*ln([1-x])*ln([1-z])*[1-x]^-1
       + 5*ln([1-x])*ln([1-z])*[1-x]^-1*z
       + 5*ln([1-x])*ln([1-z])*[1-z]^-1
       - 11*ln([1-x])*ln([1-z])*z
       + 3*ln([1-x])*ln([1-z])*x*[1-z]^-1
       - 9*ln([1-x])*ln([1-z])*x
       - ln([1-x])*ln([1-z])*x*z
       - 43/4*ln(z)
       - 12*ln(z)*[1-x]^-1*z^-1*[1-z]^-1
       + 27/2*ln(z)*[1-x]^-1*z^-1
       - 3*ln(z)*[1-x]^-1*z^-1*rln2
       - 3/2*ln(z)*[1-x]^-1*z^-1*sqrtxz1
       - 3/2*ln(z)*[1-x]^-1*[1-z]^-1
       - 3/2*ln(z)*[1-x]^-1*[x-z]^-1
       + 10*ln(z)*[1-x]^-1
       - 6*ln(z)*[1-x]^-1*rln2
       - 3/2*ln(z)*[1-x]^-1*sqrtxz1
       + 27/2*ln(z)*[1-x]^-1*z
       - 3*ln(z)*[1-x]^-1*z*rln2
       + 8*ln(z)*z^-1*[1-z]^-1
       - 8*ln(z)*z^-1
       + 3/2*ln(z)*z^-1*rln2
       + 1/2*ln(z)*z^-1*sqrtxz1
       + 13/2*ln(z)*[1-z]^-1
       + 3/2*ln(z)*[x-z]^-1
       - 3/4*ln(z)*poly2^-2
       + 1/2*ln(z)*poly2^-1
       + 3*ln(z)*rln2
       + 3/2*ln(z)*sqrtxz1
       - 39/2*ln(z)*z
       + 3*ln(z)*z*rln2
       + 4*ln(z)*x*z^-1*[1-z]^-1
       - 11/2*ln(z)*x*z^-1
       + 3/2*ln(z)*x*z^-1*rln2
       + ln(z)*x*z^-1*sqrtxz1
       - 7/2*ln(z)*x*[1-z]^-1
       + 3/2*ln(z)*x*[x-z]^-1
       - 3/4*ln(z)*x*poly2^-2
       - 31/4*ln(z)*x
       + 3*ln(z)*x*rln2
       + 2*ln(z)*x^2*[x-z]^-1
       - 1/2*ln(z)*x^3*poly2^-1
       + 3/4*ln(z)*x^4*poly2^-2
       + 3/4*ln(z)*x^5*poly2^-2
       - ln(z)*ln(1 + sqrtxz1 - z)
       + ln(z)*ln(1 + sqrtxz1 - z)*[1-x]^-1*z^-1
       + 2*ln(z)*ln(1 + sqrtxz1 - z)*[1-x]^-1
       + ln(z)*ln(1 + sqrtxz1 - z)*[1-x]^-1*z
       - 1/2*ln(z)*ln(1 + sqrtxz1 - z)*z^-1
       - ln(z)*ln(1 + sqrtxz1 - z)*z
       - 1/2*ln(z)*ln(1 + sqrtxz1 - z)*x*z^-1
       - ln(z)*ln(1 + sqrtxz1 - z)*x
       - 2*ln(z)*ln(1 + sqrtxz1 + z)
       + 2*ln(z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^-1
       + 4*ln(z)*ln(1 + sqrtxz1 + z)*[1-x]^-1
       + 2*ln(z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z
       - ln(z)*ln(1 + sqrtxz1 + z)*z^-1
       - 2*ln(z)*ln(1 + sqrtxz1 + z)*z
       - ln(z)*ln(1 + sqrtxz1 + z)*x*z^-1
       - 2*ln(z)*ln(1 + sqrtxz1 + z)*x
       - ln(z)*ln(1 + x*z^-1)
       - 1/2*ln(z)*ln(1 + x*z^-1)*x^-2*[1+x]^-1*z^-1
       + ln(z)*ln(1 + x*z^-1)*x^-2*[1+x]^-1
       + 1/2*ln(z)*ln(1 + x*z^-1)*x^-2*z^-1
       - ln(z)*ln(1 + x*z^-1)*x^-2
       + ln(z)*ln(1 + x*z^-1)*x^-1*[1+x]^-1*z
       - 1/2*ln(z)*ln(1 + x*z^-1)*x^-1*z^-1
       + ln(z)*ln(1 + x*z^-1)*x^-1
       - ln(z)*ln(1 + x*z^-1)*x^-1*z
       - 1/2*ln(z)*ln(1 + x*z^-1)*[1+x]^-1*z^-1
       + ln(z)*ln(1 + x*z^-1)*[1+x]^-1
       + 1/2*ln(z)*ln(1 + x*z^-1)*z^-1
       + ln(z)*ln(1 + x*z^-1)*z
       - 1/2*ln(z)*ln(1 + x*z^-1)*x*z^-1
       + ln(z)*ln(1 + x*z^-1)*x
       + 2*ln(z)*ln(1 + x*z)
       + 1/2*ln(z)*ln(1 + x*z)*x^-2*[1+x]^-1*z^-1
       - ln(z)*ln(1 + x*z)*x^-2*[1+x]^-1
       - 1/2*ln(z)*ln(1 + x*z)*x^-2*z^-1
       + ln(z)*ln(1 + x*z)*x^-2
       - ln(z)*ln(1 + x*z)*x^-1*[1+x]^-1*z
       + 1/2*ln(z)*ln(1 + x*z)*x^-1*z^-1
       - ln(z)*ln(1 + x*z)*x^-1
       + ln(z)*ln(1 + x*z)*x^-1*z
       - ln(z)*ln(1 + x*z)*[1-x]^-1*z^-1
       - 2*ln(z)*ln(1 + x*z)*[1-x]^-1
       - ln(z)*ln(1 + x*z)*[1-x]^-1*z
       + 1/2*ln(z)*ln(1 + x*z)*[1+x]^-1*z^-1
       - ln(z)*ln(1 + x*z)*[1+x]^-1
       + ln(z)*ln(1 + x*z)*x*z^-1
       + ln(z)*ln(z + x)
       - ln(z)*ln(z + x)*[1-x]^-1*z^-1
       - 2*ln(z)*ln(z + x)*[1-x]^-1
       - ln(z)*ln(z + x)*[1-x]^-1*z
       + 1/2*ln(z)*ln(z + x)*z^-1
       + ln(z)*ln(z + x)*z
       + 1/2*ln(z)*ln(z + x)*x*z^-1
       + ln(z)*ln(z + x)*x
       - 1/2*ln(z)^2
       - 1/4*ln(z)^2*x^-2*[1+x]^-1*z^-1
       + 1/2*ln(z)^2*x^-2*[1+x]^-1
       + 1/4*ln(z)^2*x^-2*z^-1
       - 1/2*ln(z)^2*x^-2
       + 1/2*ln(z)^2*x^-1*[1+x]^-1*z
       - 1/4*ln(z)^2*x^-1*z^-1
       + 1/2*ln(z)^2*x^-1
       - 1/2*ln(z)^2*x^-1*z
       - ln(z)^2*[1-x]^-1*z^-1*[1-z]^-1
       + ln(z)^2*[1-x]^-1*z^-1
       - 5/2*ln(z)^2*[1-x]^-1*[1-z]^-1
       + 2*ln(z)^2*[1-x]^-1
       + 2*ln(z)^2*[1-x]^-1*z
       - 1/4*ln(z)^2*[1+x]^-1*z^-1
       + 1/2*ln(z)^2*[1+x]^-1
       + 1/2*ln(z)^2*z^-1*[1-z]^-1
       - 1/4*ln(z)^2*z^-1
       + 1/4*ln(z)^2*[1-z]^-1
       + 1/2*ln(z)^2*x*z^-1*[1-z]^-1
       - 3/4*ln(z)^2*x*z^-1
       - 1/4*ln(z)^2*x*[1-z]^-1
       + 1/2*ln(z)^2*x
       + 1/2*ln(z)^2*x*z
       - 3*ln(z)*ln([1-z])
       - 8*ln(z)*ln([1-z])*[1-x]^-1*[1-z]^-1
       + 7*ln(z)*ln([1-z])*[1-x]^-1
       + 7*ln(z)*ln([1-z])*[1-x]^-1*z
       + 13/2*ln(z)*ln([1-z])*[1-z]^-1
       - 12*ln(z)*ln([1-z])*z
       + 7/2*ln(z)*ln([1-z])*x*[1-z]^-1
       - 9*ln(z)*ln([1-z])*x
       - 6*ln([1-z])
       - 15*ln([1-z])*[1-x]^-1*[1-z]^-1
       + 7*ln([1-z])*[1-x]^-1
       + 12*ln([1-z])*[1-x]^-1*z
       + 15*ln([1-z])*[1-z]^-1
       - ln([1-z])*[1-x-z]^-1
       + 1/2*ln([1-z])*z*[1-x-z]^-2
       - ln([1-z])*z*[1-x-z]^-1
       - 67/6*ln([1-z])*z
       - 1/2*ln([1-z])*z^2*[1-x-z]^-2
       - 31/6*ln([1-z])*x
       - 2*ln([1-z])*x*z
       - ln([1-z])^2
       - 9/2*ln([1-z])^2*[1-x]^-1*[1-z]^-1
       + 13/4*ln([1-z])^2*[1-x]^-1
       + 13/4*ln([1-z])^2*[1-x]^-1*z
       + 5/2*ln([1-z])^2*[1-z]^-1
       - 6*ln([1-z])^2*z
       + 2*ln([1-z])^2*x*[1-z]^-1
       - 11/2*ln([1-z])^2*x
       - 1/2*ln([1-z])^2*x*z
       + 3*ln(sqrtxz3)*ArcTan(sqrtxz3)*x^-1*sqrtxz3
       - ln(sqrtxz3)*ArcTan(sqrtxz3)*z^-1*sqrtxz3
       + 3*ln(sqrtxz3)*ArcTan(sqrtxz3)*z*sqrtxz3
       + 3*ln(sqrtxz3)*ArcTan(sqrtxz3)*x*sqrtxz3
       + 1/2*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*z^-1*sqrtxz2^-1
       + 1/2*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*[1-z]^-1*sqrtxz2^-1
       - 3/8*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*sqrtxz2^-1*poly2^-2
       + 3/8*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*sqrtxz2^-1*poly2^-1
       - 2*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*sqrtxz2^-1
       - Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x*z^-1*sqrtxz2^-1
       + Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x*[1-z]^-1*sqrtxz2^-1
       - 3/4*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x*sqrtxz2^-1
       + 3/2*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x*z*sqrtxz2^-1
       + 1/2*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x^2*z^-1*sqrtxz2^-1
       + 1/2*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x^2*[1-z]^-1*sqrtxz2^-1
       + 3/8*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x^2*sqrtxz2^-1*poly2^-2
       + 1/4*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x^2*sqrtxz2^-1*poly2^-1
       - 2*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x^2*sqrtxz2^-1
       + 3/8*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x^4*sqrtxz2^-1*poly2^-2
       + 3/8*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x^4*sqrtxz2^-1*poly2^-1
       - 3/8*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x^6*sqrtxz2^-1*poly2^-2
       - 1/2*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*z^-1*sqrtxz2^-1
       - 1/2*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*[1-z]^-1*sqrtxz2^-1
       + 3/8*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*sqrtxz2^-1*poly2^-2
       - 3/8*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*sqrtxz2^-1*poly2^-1
       + 2*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*sqrtxz2^-1
       + Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x*z^-1*sqrtxz2^-1
       - Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x*[1-z]^-1*sqrtxz2^-1
       + 3/4*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x*sqrtxz2^-1
       - 3/2*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x*z*sqrtxz2^-1
       - 1/2*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x^2*z^-1*sqrtxz2^-1
       - 1/2*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x^2*[1-z]^-1*sqrtxz2^-1
       - 3/8*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x^2*sqrtxz2^-1*poly2^-2
       - 1/4*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x^2*sqrtxz2^-1*poly2^-1
       + 2*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x^2*sqrtxz2^-1
       - 3/8*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x^4*sqrtxz2^-1*poly2^-2
       - 3/8*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x^4*sqrtxz2^-1*poly2^-1
       + 3/8*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x^6*sqrtxz2^-1*poly2^-2
       + Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)
       - Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*z^-1
       - 2*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1
       - Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*z
       + 1/2*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z^-1
       + Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z
       + 1/2*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x*z^-1
       + Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x
       - Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)
       + Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*z^-1
       + 2*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1
       + Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*z
       - 1/2*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z^-1
       - Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z
       - 1/2*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x*z^-1
       - Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x
       - 1/2*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*z^-1*sqrtxz2^-1
       - 1/2*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*[1-z]^-1*sqrtxz2^-1
       + 3/8*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*sqrtxz2^-1*poly2^-2
       - 3/8*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*sqrtxz2^-1*poly2^-1
       + 2*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*sqrtxz2^-1
       + Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x*z^-1*sqrtxz2^-1
       - Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x*[1-z]^-1*sqrtxz2^-1
       + 3/4*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x*sqrtxz2^-1
       - 3/2*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x*z*sqrtxz2^-1
       - 1/2*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x^2*z^-1*sqrtxz2^-1
       - 1/2*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x^2*[1-z]^-1*sqrtxz2^-1
       - 3/8*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x^2*sqrtxz2^-1*poly2^-2
       - 1/4*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x^2*sqrtxz2^-1*poly2^-1
       + 2*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x^2*sqrtxz2^-1
       - 3/8*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x^4*sqrtxz2^-1*poly2^-2
       - 3/8*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x^4*sqrtxz2^-1*poly2^-1
       + 3/8*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x^6*sqrtxz2^-1*poly2^-2
       + 1/2*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*z^-1*sqrtxz2^-1
       + 1/2*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*[1-z]^-1*sqrtxz2^-1
       - 3/8*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*sqrtxz2^-1*poly2^-2
       + 3/8*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*sqrtxz2^-1*poly2^-1
       - 2*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*sqrtxz2^-1
       - Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x*z^-1*sqrtxz2^-1
       + Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x*[1-z]^-1*sqrtxz2^-1
       - 3/4*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x*sqrtxz2^-1
       + 3/2*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x*z*sqrtxz2^-1
       + 1/2*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x^2*z^-1*sqrtxz2^-1
       + 1/2*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x^2*[1-z]^-1*sqrtxz2^-1
       + 3/8*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x^2*sqrtxz2^-1*poly2^-2
       + 1/4*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x^2*sqrtxz2^-1*poly2^-1
       - 2*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x^2*sqrtxz2^-1
       + 3/8*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x^4*sqrtxz2^-1*poly2^-2
       + 3/8*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x^4*sqrtxz2^-1*poly2^-1
       - 3/8*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x^6*sqrtxz2^-1*poly2^-2
       - Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)
       + Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*[1-x]^-1*z^-1
       + 2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*[1-x]^-1
       + Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*[1-x]^-1*z
       - 1/2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*z^-1
       - Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*z
       - 1/2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*x*z^-1
       - Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*x
       + Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)
       - Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*[1-x]^-1*z^-1
       - 2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*[1-x]^-1
       - Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*[1-x]^-1*z
       + 1/2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*z^-1
       + Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*z
       + 1/2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*x*z^-1
       + Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*x
       + 2*Li2(1 - x*z^-1)
       + 4*Li2(1 - x*z^-1)*[1-x]^-1*z^-1*[1-z]^-1
       - 4*Li2(1 - x*z^-1)*[1-x]^-1*z^-1
       - 2*Li2(1 - x*z^-1)*[1-x]^-1*[1-z]^-1
       - 5/2*Li2(1 - x*z^-1)*[1-x]^-1
       - 5/2*Li2(1 - x*z^-1)*[1-x]^-1*z
       - 3*Li2(1 - x*z^-1)*z^-1*[1-z]^-1
       + 3*Li2(1 - x*z^-1)*z^-1
       + Li2(1 - x*z^-1)*[1-z]^-1
       + 4*Li2(1 - x*z^-1)*z
       - Li2(1 - x*z^-1)*x*z^-1*[1-z]^-1
       + Li2(1 - x*z^-1)*x*z^-1
       + Li2(1 - x*z^-1)*x*[1-z]^-1
       + 2*Li2(1 - x*z^-1)*x
       + 1/2*Li2( - x*z^-1)*x^-2*[1+x]^-1*z^-1
       - Li2( - x*z^-1)*x^-2*[1+x]^-1
       - 1/2*Li2( - x*z^-1)*x^-2*z^-1
       + Li2( - x*z^-1)*x^-2
       - Li2( - x*z^-1)*x^-1*[1+x]^-1*z
       + 1/2*Li2( - x*z^-1)*x^-1*z^-1
       - Li2( - x*z^-1)*x^-1
       + Li2( - x*z^-1)*x^-1*z
       + Li2( - x*z^-1)*[1-x]^-1*z^-1
       + 2*Li2( - x*z^-1)*[1-x]^-1
       + Li2( - x*z^-1)*[1-x]^-1*z
       + 1/2*Li2( - x*z^-1)*[1+x]^-1*z^-1
       - Li2( - x*z^-1)*[1+x]^-1
       - Li2( - x*z^-1)*z^-1
       - 2*Li2( - x*z^-1)*z
       - 2*Li2( - x*z^-1)*x
       + 2*Li2( - x)
       - Li2( - x)*x^-2*[1+x]^-1*z^-1
       + 2*Li2( - x)*x^-2*[1+x]^-1
       + Li2( - x)*x^-2*z^-1
       - 2*Li2( - x)*x^-2
       + 2*Li2( - x)*x^-1*[1+x]^-1*z
       - Li2( - x)*x^-1*z^-1
       + 2*Li2( - x)*x^-1
       - 2*Li2( - x)*x^-1*z
       + Li2( - x)*[1+x]^-1*z^-1
       - 2*Li2( - x)*[1+x]^-1
       + 2*Li2( - x)*[1+x]^-1*z
       - Li2( - x)*z^-1
       - Li2( - x)*x*z^-1
       + 2*Li2( - x)*x
       + 2*Li2( - x*z)
       + 1/2*Li2( - x*z)*x^-2*[1+x]^-1*z^-1
       - Li2( - x*z)*x^-2*[1+x]^-1
       - 1/2*Li2( - x*z)*x^-2*z^-1
       + Li2( - x*z)*x^-2
       - Li2( - x*z)*x^-1*[1+x]^-1*z
       + 1/2*Li2( - x*z)*x^-1*z^-1
       - Li2( - x*z)*x^-1
       + Li2( - x*z)*x^-1*z
       - Li2( - x*z)*[1-x]^-1*z^-1
       - 2*Li2( - x*z)*[1-x]^-1
       - Li2( - x*z)*[1-x]^-1*z
       + 1/2*Li2( - x*z)*[1+x]^-1*z^-1
       - Li2( - x*z)*[1+x]^-1
       + Li2( - x*z)*x*z^-1
       + 2*Li2(x)
       + Li2(x)*x^-2*[1+x]^-1*z^-1
       - 2*Li2(x)*x^-2*[1+x]^-1
       - Li2(x)*x^-2*z^-1
       + 2*Li2(x)*x^-2
       - 2*Li2(x)*x^-1*[1+x]^-1*z
       + Li2(x)*x^-1*z^-1
       - 2*Li2(x)*x^-1
       + 2*Li2(x)*x^-1*z
       - 4*Li2(x)*[1-x]^-1*z^-1*[1-z]^-1
       + 4*Li2(x)*[1-x]^-1*z^-1
       + 6*Li2(x)*[1-x]^-1*[1-z]^-1
       - 3/2*Li2(x)*[1-x]^-1
       - 3/2*Li2(x)*[1-x]^-1*z
       - Li2(x)*[1+x]^-1*z^-1
       + 2*Li2(x)*[1+x]^-1
       - 2*Li2(x)*[1+x]^-1*z
       + 2*Li2(x)*z^-1*[1-z]^-1
       - 2*Li2(x)*z^-1
       - 4*Li2(x)*[1-z]^-1
       + Li2(x)*z
       + 2*Li2(x)*x*z^-1*[1-z]^-1
       - 2*Li2(x)*x*z^-1
       - 2*Li2(x)*x*[1-z]^-1
       - Li2(x)*x
       - Li2(z)
       - 2*Li2(z)*[1-x]^-1*[1-z]^-1
       + 5/2*Li2(z)*[1-x]^-1
       + 5/2*Li2(z)*[1-x]^-1*z
       + 3/2*Li2(z)*[1-z]^-1
       - 5*Li2(z)*z
       + 1/2*Li2(z)*x*[1-z]^-1
       - 4*Li2(z)*x
       + 3/2*InvTanInt( - sqrtxz3)*x^-1*sqrtxz3
       - 1/2*InvTanInt( - sqrtxz3)*z^-1*sqrtxz3
       + 3/2*InvTanInt( - sqrtxz3)*z*sqrtxz3
       + 3/2*InvTanInt( - sqrtxz3)*x*sqrtxz3
       + 3*InvTanInt(z*sqrtxz3)*x^-1*sqrtxz3
       - InvTanInt(z*sqrtxz3)*z^-1*sqrtxz3
       + 3*InvTanInt(z*sqrtxz3)*z*sqrtxz3
       + 3*InvTanInt(z*sqrtxz3)*x*sqrtxz3
       - 3/2*InvTanInt(sqrtxz3)*x^-1*sqrtxz3
       + 1/2*InvTanInt(sqrtxz3)*z^-1*sqrtxz3
       - 3/2*InvTanInt(sqrtxz3)*z*sqrtxz3
       - 3/2*InvTanInt(sqrtxz3)*x*sqrtxz3 );

id DC2Q2QPS    = (  + NC^-1 * (  - 4
       + 4*z^-1
       - 4*x*z^-1
       + 4*x
       - 2*ln(z*sqrtxz3)*ArcTan(z*sqrtxz3)*x^-1*sqrtxz3
       - 2*ln(z*sqrtxz3)*ArcTan(z*sqrtxz3)*z^-1*sqrtxz3
       - 2*ln(z*sqrtxz3)*ArcTan(z*sqrtxz3)*z*sqrtxz3
       - 2*ln(z*sqrtxz3)*ArcTan(z*sqrtxz3)*x*sqrtxz3
       - 3/2*ln(x)
       + 3/2*ln(x)*z^-1
       + 3/2*ln(x)*x*z^-1
       - 3/2*ln(x)*x
       + ln(x)*ln(z)
       + ln(x)*ln(z)*z^-1
       + ln(x)*ln(z)*x*z^-1
       + ln(x)*ln(z)*x
       + 3/2*ln(z)
       + 3/2*ln(z)*z^-1
       - 3/2*ln(z)*x*z^-1
       - 3/2*ln(z)*x
       + 2*ln(sqrtxz3)*ArcTan(sqrtxz3)*x^-1*sqrtxz3
       + 2*ln(sqrtxz3)*ArcTan(sqrtxz3)*z^-1*sqrtxz3
       + 2*ln(sqrtxz3)*ArcTan(sqrtxz3)*z*sqrtxz3
       + 2*ln(sqrtxz3)*ArcTan(sqrtxz3)*x*sqrtxz3
       + InvTanInt( - sqrtxz3)*x^-1*sqrtxz3
       + InvTanInt( - sqrtxz3)*z^-1*sqrtxz3
       + InvTanInt( - sqrtxz3)*z*sqrtxz3
       + InvTanInt( - sqrtxz3)*x*sqrtxz3
       + 2*InvTanInt(z*sqrtxz3)*x^-1*sqrtxz3
       + 2*InvTanInt(z*sqrtxz3)*z^-1*sqrtxz3
       + 2*InvTanInt(z*sqrtxz3)*z*sqrtxz3
       + 2*InvTanInt(z*sqrtxz3)*x*sqrtxz3
       - InvTanInt(sqrtxz3)*x^-1*sqrtxz3
       - InvTanInt(sqrtxz3)*z^-1*sqrtxz3
       - InvTanInt(sqrtxz3)*z*sqrtxz3
       - InvTanInt(sqrtxz3)*x*sqrtxz3
       )

       + NC * ( 4
       - 4*z^-1
       + 4*x*z^-1
       - 4*x
       + 2*ln(z*sqrtxz3)*ArcTan(z*sqrtxz3)*x^-1*sqrtxz3
       + 2*ln(z*sqrtxz3)*ArcTan(z*sqrtxz3)*z^-1*sqrtxz3
       + 2*ln(z*sqrtxz3)*ArcTan(z*sqrtxz3)*z*sqrtxz3
       + 2*ln(z*sqrtxz3)*ArcTan(z*sqrtxz3)*x*sqrtxz3
       + 3/2*ln(x)
       - 3/2*ln(x)*z^-1
       - 3/2*ln(x)*x*z^-1
       + 3/2*ln(x)*x
       - ln(x)*ln(z)
       - ln(x)*ln(z)*z^-1
       - ln(x)*ln(z)*x*z^-1
       - ln(x)*ln(z)*x
       - 3/2*ln(z)
       - 3/2*ln(z)*z^-1
       + 3/2*ln(z)*x*z^-1
       + 3/2*ln(z)*x
       - 2*ln(sqrtxz3)*ArcTan(sqrtxz3)*x^-1*sqrtxz3
       - 2*ln(sqrtxz3)*ArcTan(sqrtxz3)*z^-1*sqrtxz3
       - 2*ln(sqrtxz3)*ArcTan(sqrtxz3)*z*sqrtxz3
       - 2*ln(sqrtxz3)*ArcTan(sqrtxz3)*x*sqrtxz3
       - InvTanInt( - sqrtxz3)*x^-1*sqrtxz3
       - InvTanInt( - sqrtxz3)*z^-1*sqrtxz3
       - InvTanInt( - sqrtxz3)*z*sqrtxz3
       - InvTanInt( - sqrtxz3)*x*sqrtxz3
       - 2*InvTanInt(z*sqrtxz3)*x^-1*sqrtxz3
       - 2*InvTanInt(z*sqrtxz3)*z^-1*sqrtxz3
       - 2*InvTanInt(z*sqrtxz3)*z*sqrtxz3
       - 2*InvTanInt(z*sqrtxz3)*x*sqrtxz3
       + InvTanInt(sqrtxz3)*x^-1*sqrtxz3
       + InvTanInt(sqrtxz3)*z^-1*sqrtxz3
       + InvTanInt(sqrtxz3)*z*sqrtxz3
       + InvTanInt(sqrtxz3)*x*sqrtxz3
       )

       + Dd([1-z])*NC^-1 * ( 1
       - x
       + 3/2*ln(x)
       - 1/2*ln(x)*x
       + 1/2*ln(x)^2
       + 1/4*ln(x)^2*x
       )

       + Dd([1-z])*NC * (  - 1
       + x
       - 3/2*ln(x)
       + 1/2*ln(x)*x
       - 1/2*ln(x)^2
       - 1/4*ln(x)^2*x ) );

id DC2Q2QB     = (  + NC^-2 * ( 6
       - 4*[1-x]^-1*z^-1*[1+z]^-1*rln2^2
       - 8*[1-x]^-1*z^-1*[1+z]^-1*sqrtxz1*rln2
       - 6*[1-x]^-1*z^-1*[1+z]^-1*sqrtxz1*rln2^2
       + 4*[1-x]^-1*z^-1*rln2^2
       + 8*[1-x]^-1*z^-1*sqrtxz1*rln2
       + 6*[1-x]^-1*z^-1*sqrtxz1*rln2^2
       + 8*[1-x]^-1*[1+z]^-1*sqrtxz1^-1*rln2
       + 6*[1-x]^-1*[1+z]^-1*sqrtxz1^-1*rln2^2
       - 16*[1-x]^-1*sqrtxz1^-1*rln2
       - 12*[1-x]^-1*sqrtxz1^-1*rln2^2
       - 3*[1-x]^-1*rln2^2
       + [1-x]^-1*sqrtxz1*rln2
       + 8*[1-x]^-1*z*[1+z]^-1*sqrtxz1^-1*rln2
       + 6*[1-x]^-1*z*[1+z]^-1*sqrtxz1^-1*rln2^2
       - 8*[1-x]^-1*z*sqrtxz1^-1*rln2
       - 6*[1-x]^-1*z*sqrtxz1^-1*rln2^2
       + 3*[1-x]^-1*z*rln2^2
       + 4*z^-1*[1+z]^-1*rln2^2
       + 8*z^-1*[1+z]^-1*sqrtxz1*rln2
       + 6*z^-1*[1+z]^-1*sqrtxz1*rln2^2
       - 2*z^-1
       - 4*z^-1*rln2^2
       - 8*z^-1*sqrtxz1*rln2
       - 6*z^-1*sqrtxz1*rln2^2
       - 8*[1+z]^-1*sqrtxz1^-1*rln2
       - 6*[1+z]^-1*sqrtxz1^-1*rln2^2
       + 16*sqrtxz1^-1*rln2
       + 12*sqrtxz1^-1*rln2^2
       + 4*rln2^2
       - 8*z*[1+z]^-1*sqrtxz1^-1*rln2
       - 6*z*[1+z]^-1*sqrtxz1^-1*rln2^2
       + 8*z*sqrtxz1^-1*rln2
       + 6*z*sqrtxz1^-1*rln2^2
       - z
       - 4*z*rln2^2
       + 7/4*x*z^-1
       - 32*x*[1+z]^-1*sqrtxz1^-1*rln2
       - 24*x*[1+z]^-1*sqrtxz1^-1*rln2^2
       + 32*x*sqrtxz1^-1*rln2
       + 24*x*sqrtxz1^-1*rln2^2
       - 11/2*x
       + 3/4*x*z
       + 1/3*pi^2*x^-1*[1+x]^-1*z^-1*[1-z]^-1
       - 1/3*pi^2*x^-1*[1+x]^-1*z^-1
       - 1/3*pi^2*x^-1*[1+x]^-1
       - 1/3*pi^2*x^-1*[1+x]^-1*z
       - 1/3*pi^2*x^-1*z^-1*[1-z]^-1
       + 1/3*pi^2*x^-1*z^-1
       + 1/3*pi^2*x^-1
       + 1/3*pi^2*x^-1*z
       - 1/3*pi^2*[1-x]^-1*z^-1*[1+z]^-1*sqrtxz1
       - 1/6*pi^2*[1-x]^-1*z^-1
       + 1/3*pi^2*[1-x]^-1*z^-1*sqrtxz1
       + 1/3*pi^2*[1-x]^-1*[1+z]^-1*sqrtxz1^-1
       + 1/6*pi^2*[1-x]^-1*[1+z]^-1
       - 2/3*pi^2*[1-x]^-1*sqrtxz1^-1
       + 1/4*pi^2*[1-x]^-1
       + 1/3*pi^2*[1-x]^-1*z*[1+z]^-1*sqrtxz1^-1
       - 1/3*pi^2*[1-x]^-1*z*sqrtxz1^-1
       - 1/12*pi^2*[1-x]^-1*z
       + 1/3*pi^2*[1+x]^-1*[1-z]^-1
       - 1/2*pi^2*[1+x]^-1
       - 1/2*pi^2*[1+x]^-1*z
       + 1/3*pi^2*z^-1*[1-z]^-1
       + 1/3*pi^2*z^-1*[1+z]^-1*sqrtxz1
       - 1/4*pi^2*z^-1
       - 1/3*pi^2*z^-1*sqrtxz1
       - 1/3*pi^2*[1-z]^-1
       - 1/3*pi^2*[1+z]^-1*sqrtxz1^-1
       + 2/3*pi^2*sqrtxz1^-1
       - 1/6*pi^2
       - 1/3*pi^2*z*[1+z]^-1*sqrtxz1^-1
       + 1/3*pi^2*z*sqrtxz1^-1
       + 1/3*pi^2*z
       + 1/12*pi^2*x*z^-1
       - 4/3*pi^2*x*[1+z]^-1*sqrtxz1^-1
       + 4/3*pi^2*x*sqrtxz1^-1
       - 1/3*pi^2*x
       + 4*ln(1 + sqrtxz1 - z)*[1-x]^-1*z^-1*[1+z]^-1*rln2
       + 8*ln(1 + sqrtxz1 - z)*[1-x]^-1*z^-1*[1+z]^-1*sqrtxz1
       + 8*ln(1 + sqrtxz1 - z)*[1-x]^-1*z^-1*[1+z]^-1*sqrtxz1*rln2
       - 4*ln(1 + sqrtxz1 - z)*[1-x]^-1*z^-1*rln2
       - 8*ln(1 + sqrtxz1 - z)*[1-x]^-1*z^-1*sqrtxz1
       - 8*ln(1 + sqrtxz1 - z)*[1-x]^-1*z^-1*sqrtxz1*rln2
       - 8*ln(1 + sqrtxz1 - z)*[1-x]^-1*[1+z]^-1*sqrtxz1^-1
       - 8*ln(1 + sqrtxz1 - z)*[1-x]^-1*[1+z]^-1*sqrtxz1^-1*rln2
       + 16*ln(1 + sqrtxz1 - z)*[1-x]^-1*sqrtxz1^-1
       + 16*ln(1 + sqrtxz1 - z)*[1-x]^-1*sqrtxz1^-1*rln2
       + 3*ln(1 + sqrtxz1 - z)*[1-x]^-1*rln2
       - ln(1 + sqrtxz1 - z)*[1-x]^-1*sqrtxz1
       - 8*ln(1 + sqrtxz1 - z)*[1-x]^-1*z*[1+z]^-1*sqrtxz1^-1
       - 8*ln(1 + sqrtxz1 - z)*[1-x]^-1*z*[1+z]^-1*sqrtxz1^-1*rln2
       + 8*ln(1 + sqrtxz1 - z)*[1-x]^-1*z*sqrtxz1^-1
       + 8*ln(1 + sqrtxz1 - z)*[1-x]^-1*z*sqrtxz1^-1*rln2
       - 3*ln(1 + sqrtxz1 - z)*[1-x]^-1*z*rln2
       - 4*ln(1 + sqrtxz1 - z)*z^-1*[1+z]^-1*rln2
       - 8*ln(1 + sqrtxz1 - z)*z^-1*[1+z]^-1*sqrtxz1
       - 8*ln(1 + sqrtxz1 - z)*z^-1*[1+z]^-1*sqrtxz1*rln2
       + 4*ln(1 + sqrtxz1 - z)*z^-1*rln2
       + 8*ln(1 + sqrtxz1 - z)*z^-1*sqrtxz1
       + 8*ln(1 + sqrtxz1 - z)*z^-1*sqrtxz1*rln2
       + 8*ln(1 + sqrtxz1 - z)*[1+z]^-1*sqrtxz1^-1
       + 8*ln(1 + sqrtxz1 - z)*[1+z]^-1*sqrtxz1^-1*rln2
       - 16*ln(1 + sqrtxz1 - z)*sqrtxz1^-1
       - 16*ln(1 + sqrtxz1 - z)*sqrtxz1^-1*rln2
       - 4*ln(1 + sqrtxz1 - z)*rln2
       + 8*ln(1 + sqrtxz1 - z)*z*[1+z]^-1*sqrtxz1^-1
       + 8*ln(1 + sqrtxz1 - z)*z*[1+z]^-1*sqrtxz1^-1*rln2
       - 8*ln(1 + sqrtxz1 - z)*z*sqrtxz1^-1
       - 8*ln(1 + sqrtxz1 - z)*z*sqrtxz1^-1*rln2
       + 4*ln(1 + sqrtxz1 - z)*z*rln2
       + 32*ln(1 + sqrtxz1 - z)*x*[1+z]^-1*sqrtxz1^-1
       + 32*ln(1 + sqrtxz1 - z)*x*[1+z]^-1*sqrtxz1^-1*rln2
       - 32*ln(1 + sqrtxz1 - z)*x*sqrtxz1^-1
       - 32*ln(1 + sqrtxz1 - z)*x*sqrtxz1^-1*rln2
       - 2*ln(1 + sqrtxz1 - z)^2*[1-x]^-1*z^-1*[1+z]^-1*sqrtxz1
       + 2*ln(1 + sqrtxz1 - z)^2*[1-x]^-1*z^-1*sqrtxz1
       + 2*ln(1 + sqrtxz1 - z)^2*[1-x]^-1*[1+z]^-1*sqrtxz1^-1
       - 4*ln(1 + sqrtxz1 - z)^2*[1-x]^-1*sqrtxz1^-1
       + 2*ln(1 + sqrtxz1 - z)^2*[1-x]^-1*z*[1+z]^-1*sqrtxz1^-1
       - 2*ln(1 + sqrtxz1 - z)^2*[1-x]^-1*z*sqrtxz1^-1
       + 2*ln(1 + sqrtxz1 - z)^2*z^-1*[1+z]^-1*sqrtxz1
       - 2*ln(1 + sqrtxz1 - z)^2*z^-1*sqrtxz1
       - 2*ln(1 + sqrtxz1 - z)^2*[1+z]^-1*sqrtxz1^-1
       + 4*ln(1 + sqrtxz1 - z)^2*sqrtxz1^-1
       - 2*ln(1 + sqrtxz1 - z)^2*z*[1+z]^-1*sqrtxz1^-1
       + 2*ln(1 + sqrtxz1 - z)^2*z*sqrtxz1^-1
       - 8*ln(1 + sqrtxz1 - z)^2*x*[1+z]^-1*sqrtxz1^-1
       + 8*ln(1 + sqrtxz1 - z)^2*x*sqrtxz1^-1
       + 4*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)
       - 4*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^-1*[1+z]^-1
       - 4*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^-1*[1+z]^-1*sqrtxz1
       + 4*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^-1
       + 4*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^-1*sqrtxz1
       + 4*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*[1+z]^-1*sqrtxz1^-1
       - 8*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*sqrtxz1^-1
       - 3*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*[1-x]^-1
       + 4*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z*[1+z]^-1*sqrtxz1^-1
       - 4*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z*sqrtxz1^-1
       + 3*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z
       + 4*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*z^-1*[1+z]^-1
       + 4*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*z^-1*[1+z]^-1*sqrtxz1
       - 4*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*z^-1
       - 4*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*z^-1*sqrtxz1
       - 4*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*[1+z]^-1*sqrtxz1^-1
       + 8*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*sqrtxz1^-1
       - 4*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*z*[1+z]^-1*sqrtxz1^-1
       + 4*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*z*sqrtxz1^-1
       - 4*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*z
       - 16*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*x*[1+z]^-1*sqrtxz1^-1
       + 16*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*x*sqrtxz1^-1
       + 4*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^-1*[1+z]^-1*rln2
       + 4*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^-1*[1+z]^-1*sqrtxz1*rln2
       - 4*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^-1*rln2
       - 4*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^-1*sqrtxz1*rln2
       - 4*ln(1 + sqrtxz1 + z)*[1-x]^-1*[1+z]^-1*sqrtxz1^-1*rln2
       + 8*ln(1 + sqrtxz1 + z)*[1-x]^-1*sqrtxz1^-1*rln2
       + 3*ln(1 + sqrtxz1 + z)*[1-x]^-1*rln2
       - 4*ln(1 + sqrtxz1 + z)*[1-x]^-1*z*[1+z]^-1*sqrtxz1^-1*rln2
       + 4*ln(1 + sqrtxz1 + z)*[1-x]^-1*z*sqrtxz1^-1*rln2
       - 3*ln(1 + sqrtxz1 + z)*[1-x]^-1*z*rln2
       - 4*ln(1 + sqrtxz1 + z)*z^-1*[1+z]^-1*rln2
       - 4*ln(1 + sqrtxz1 + z)*z^-1*[1+z]^-1*sqrtxz1*rln2
       + 4*ln(1 + sqrtxz1 + z)*z^-1*rln2
       + 4*ln(1 + sqrtxz1 + z)*z^-1*sqrtxz1*rln2
       + 4*ln(1 + sqrtxz1 + z)*[1+z]^-1*sqrtxz1^-1*rln2
       - 8*ln(1 + sqrtxz1 + z)*sqrtxz1^-1*rln2
       - 4*ln(1 + sqrtxz1 + z)*rln2
       + 4*ln(1 + sqrtxz1 + z)*z*[1+z]^-1*sqrtxz1^-1*rln2
       - 4*ln(1 + sqrtxz1 + z)*z*sqrtxz1^-1*rln2
       + 4*ln(1 + sqrtxz1 + z)*z*rln2
       + 16*ln(1 + sqrtxz1 + z)*x*[1+z]^-1*sqrtxz1^-1*rln2
       - 16*ln(1 + sqrtxz1 + z)*x*sqrtxz1^-1*rln2
       - 1/2*ln(1 - 2*z + z^2 + 4*x*z)^2*[1-x]^-1*z^-1*[1+z]^-1*sqrtxz1
       + 1/2*ln(1 - 2*z + z^2 + 4*x*z)^2*[1-x]^-1*z^-1*sqrtxz1
       + 1/2*ln(1 - 2*z + z^2 + 4*x*z)^2*[1-x]^-1*[1+z]^-1*sqrtxz1^-1
       - ln(1 - 2*z + z^2 + 4*x*z)^2*[1-x]^-1*sqrtxz1^-1
       + 1/2*ln(1 - 2*z + z^2 + 4*x*z)^2*[1-x]^-1*z*[1+z]^-1*sqrtxz1^-1
       - 1/2*ln(1 - 2*z + z^2 + 4*x*z)^2*[1-x]^-1*z*sqrtxz1^-1
       + 1/2*ln(1 - 2*z + z^2 + 4*x*z)^2*z^-1*[1+z]^-1*sqrtxz1
       - 1/2*ln(1 - 2*z + z^2 + 4*x*z)^2*z^-1*sqrtxz1
       - 1/2*ln(1 - 2*z + z^2 + 4*x*z)^2*[1+z]^-1*sqrtxz1^-1
       + ln(1 - 2*z + z^2 + 4*x*z)^2*sqrtxz1^-1
       - 1/2*ln(1 - 2*z + z^2 + 4*x*z)^2*z*[1+z]^-1*sqrtxz1^-1
       + 1/2*ln(1 - 2*z + z^2 + 4*x*z)^2*z*sqrtxz1^-1
       - 2*ln(1 - 2*z + z^2 + 4*x*z)^2*x*[1+z]^-1*sqrtxz1^-1
       + 2*ln(1 - 2*z + z^2 + 4*x*z)^2*x*sqrtxz1^-1
       + 2*ln(x)
       + 4*ln(x)*x^-1*[1+x]^-1*z^-1*[1-z]^-1
       - 4*ln(x)*x^-1*[1+x]^-1*z^-1
       - 4*ln(x)*x^-1*[1+x]^-1
       - 4*ln(x)*x^-1*[1+x]^-1*z
       - 4*ln(x)*x^-1*z^-1*[1-z]^-1
       + 4*ln(x)*x^-1*z^-1
       + 4*ln(x)*x^-1
       + 4*ln(x)*x^-1*z
       - 2*ln(x)*[1-x]^-1*z^-1*[1+z]^-1*rln2
       - 4*ln(x)*[1-x]^-1*z^-1*[1+z]^-1*sqrtxz1
       - 3/4*ln(x)*[1-x]^-1*z^-1
       + 2*ln(x)*[1-x]^-1*z^-1*rln2
       + 4*ln(x)*[1-x]^-1*z^-1*sqrtxz1
       + 4*ln(x)*[1-x]^-1*[1+z]^-1*sqrtxz1^-1
       - 8*ln(x)*[1-x]^-1*sqrtxz1^-1
       + ln(x)*[1-x]^-1
       - 3/2*ln(x)*[1-x]^-1*rln2
       + 1/2*ln(x)*[1-x]^-1*sqrtxz1
       + 4*ln(x)*[1-x]^-1*z*[1+z]^-1*sqrtxz1^-1
       - 4*ln(x)*[1-x]^-1*z*sqrtxz1^-1
       - 1/4*ln(x)*[1-x]^-1*z
       + 3/2*ln(x)*[1-x]^-1*z*rln2
       + 4*ln(x)*[1+x]^-1*[1-z]^-1
       - 4*ln(x)*[1+x]^-1
       - 4*ln(x)*[1+x]^-1*z
       + 4*ln(x)*z^-1*[1-z]^-1
       + 2*ln(x)*z^-1*[1+z]^-1*rln2
       + 4*ln(x)*z^-1*[1+z]^-1*sqrtxz1
       - 9/2*ln(x)*z^-1
       - 2*ln(x)*z^-1*rln2
       - 4*ln(x)*z^-1*sqrtxz1
       - 4*ln(x)*[1-z]^-1
       - 4*ln(x)*[1+z]^-1*sqrtxz1^-1
       + 8*ln(x)*sqrtxz1^-1
       + 2*ln(x)*rln2
       - 4*ln(x)*z*[1+z]^-1*sqrtxz1^-1
       + 4*ln(x)*z*sqrtxz1^-1
       - 2*ln(x)*z*rln2
       - 1/2*ln(x)*x*z^-1
       - 16*ln(x)*x*[1+z]^-1*sqrtxz1^-1
       + 16*ln(x)*x*sqrtxz1^-1
       + 2*ln(x)*x
       - 2*ln(x)*ln(1 + sqrtxz1 - z)*[1-x]^-1*z^-1*[1+z]^-1*sqrtxz1
       + 2*ln(x)*ln(1 + sqrtxz1 - z)*[1-x]^-1*z^-1*sqrtxz1
       + 2*ln(x)*ln(1 + sqrtxz1 - z)*[1-x]^-1*[1+z]^-1*sqrtxz1^-1
       - 4*ln(x)*ln(1 + sqrtxz1 - z)*[1-x]^-1*sqrtxz1^-1
       + 2*ln(x)*ln(1 + sqrtxz1 - z)*[1-x]^-1*z*[1+z]^-1*sqrtxz1^-1
       - 2*ln(x)*ln(1 + sqrtxz1 - z)*[1-x]^-1*z*sqrtxz1^-1
       + 2*ln(x)*ln(1 + sqrtxz1 - z)*z^-1*[1+z]^-1*sqrtxz1
       - 2*ln(x)*ln(1 + sqrtxz1 - z)*z^-1*sqrtxz1
       - 2*ln(x)*ln(1 + sqrtxz1 - z)*[1+z]^-1*sqrtxz1^-1
       + 4*ln(x)*ln(1 + sqrtxz1 - z)*sqrtxz1^-1
       - 2*ln(x)*ln(1 + sqrtxz1 - z)*z*[1+z]^-1*sqrtxz1^-1
       + 2*ln(x)*ln(1 + sqrtxz1 - z)*z*sqrtxz1^-1
       - 8*ln(x)*ln(1 + sqrtxz1 - z)*x*[1+z]^-1*sqrtxz1^-1
       + 8*ln(x)*ln(1 + sqrtxz1 - z)*x*sqrtxz1^-1
       - 2*ln(x)*ln(1 + sqrtxz1 + z)
       + 2*ln(x)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^-1*[1+z]^-1
       + 2*ln(x)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^-1*[1+z]^-1*sqrtxz1
       - 2*ln(x)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^-1
       - 2*ln(x)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^-1*sqrtxz1
       - 2*ln(x)*ln(1 + sqrtxz1 + z)*[1-x]^-1*[1+z]^-1*sqrtxz1^-1
       + 4*ln(x)*ln(1 + sqrtxz1 + z)*[1-x]^-1*sqrtxz1^-1
       + 3/2*ln(x)*ln(1 + sqrtxz1 + z)*[1-x]^-1
       - 2*ln(x)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z*[1+z]^-1*sqrtxz1^-1
       + 2*ln(x)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z*sqrtxz1^-1
       - 3/2*ln(x)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z
       - 2*ln(x)*ln(1 + sqrtxz1 + z)*z^-1*[1+z]^-1
       - 2*ln(x)*ln(1 + sqrtxz1 + z)*z^-1*[1+z]^-1*sqrtxz1
       + 2*ln(x)*ln(1 + sqrtxz1 + z)*z^-1
       + 2*ln(x)*ln(1 + sqrtxz1 + z)*z^-1*sqrtxz1
       + 2*ln(x)*ln(1 + sqrtxz1 + z)*[1+z]^-1*sqrtxz1^-1
       - 4*ln(x)*ln(1 + sqrtxz1 + z)*sqrtxz1^-1
       + 2*ln(x)*ln(1 + sqrtxz1 + z)*z*[1+z]^-1*sqrtxz1^-1
       - 2*ln(x)*ln(1 + sqrtxz1 + z)*z*sqrtxz1^-1
       + 2*ln(x)*ln(1 + sqrtxz1 + z)*z
       + 8*ln(x)*ln(1 + sqrtxz1 + z)*x*[1+z]^-1*sqrtxz1^-1
       - 8*ln(x)*ln(1 + sqrtxz1 + z)*x*sqrtxz1^-1
       + ln(x)*ln(1 + x*z^-1)
       - ln(x)*ln(1 + x*z^-1)*x^-1*[1+x]^-1*z^-1*[1-z]^-1
       + ln(x)*ln(1 + x*z^-1)*x^-1*[1+x]^-1*z^-1
       + ln(x)*ln(1 + x*z^-1)*x^-1*[1+x]^-1
       + ln(x)*ln(1 + x*z^-1)*x^-1*[1+x]^-1*z
       + ln(x)*ln(1 + x*z^-1)*x^-1*z^-1*[1-z]^-1
       - ln(x)*ln(1 + x*z^-1)*x^-1*z^-1
       - ln(x)*ln(1 + x*z^-1)*x^-1
       - ln(x)*ln(1 + x*z^-1)*x^-1*z
       - ln(x)*ln(1 + x*z^-1)*z^-1*[1-z]^-1
       + ln(x)*ln(1 + x*z^-1)*z^-1
       + ln(x)*ln(1 + x*z^-1)*z
       - 2*ln(x)*ln(1 + x)*[1+x]^-1
       - 2*ln(x)*ln(1 + x)*[1+x]^-1*z
       + 2*ln(x)*ln(1 + x)*z
       - 2*ln(x)*ln(1 + x)*x
       + 2*ln(x)*ln(1 + x*z)
       - ln(x)*ln(1 + x*z)*x^-1*[1+x]^-1*z^-1*[1-z]^-1
       + ln(x)*ln(1 + x*z)*x^-1*[1+x]^-1*z^-1
       + ln(x)*ln(1 + x*z)*x^-1*[1+x]^-1
       + ln(x)*ln(1 + x*z)*x^-1*[1+x]^-1*z
       + ln(x)*ln(1 + x*z)*x^-1*z^-1*[1-z]^-1
       - ln(x)*ln(1 + x*z)*x^-1*z^-1
       - ln(x)*ln(1 + x*z)*x^-1
       - ln(x)*ln(1 + x*z)*x^-1*z
       - ln(x)*ln(1 + x*z)*[1-x]^-1*z^-1*[1+z]^-1
       + ln(x)*ln(1 + x*z)*[1-x]^-1*z^-1
       - ln(x)*ln(1 + x*z)*[1-x]^-1
       + ln(x)*ln(1 + x*z)*[1-x]^-1*z
       - ln(x)*ln(1 + x*z)*z^-1*[1-z]^-1
       + ln(x)*ln(1 + x*z)*z^-1*[1+z]^-1
       - ln(x)*ln(z + x)
       + ln(x)*ln(z + x)*[1-x]^-1*z^-1*[1+z]^-1
       - ln(x)*ln(z + x)*[1-x]^-1*z^-1
       + ln(x)*ln(z + x)*[1-x]^-1
       - ln(x)*ln(z + x)*[1-x]^-1*z
       - ln(x)*ln(z + x)*z^-1*[1+z]^-1
       + ln(x)*ln(z + x)*z^-1
       + ln(x)*ln(z + x)*z
       - 1/2*ln(x)^2
       - 3/2*ln(x)^2*x^-1*[1+x]^-1*z^-1*[1-z]^-1
       + 3/2*ln(x)^2*x^-1*[1+x]^-1*z^-1
       + 3/2*ln(x)^2*x^-1*[1+x]^-1
       + 3/2*ln(x)^2*x^-1*[1+x]^-1*z
       + 3/2*ln(x)^2*x^-1*z^-1*[1-z]^-1
       - 3/2*ln(x)^2*x^-1*z^-1
       - 3/2*ln(x)^2*x^-1
       - 3/2*ln(x)^2*x^-1*z
       + 3/2*ln(x)^2*[1-x]^-1*z^-1*[1+z]^-1*sqrtxz1
       - 1/2*ln(x)^2*[1-x]^-1*z^-1
       - 3/2*ln(x)^2*[1-x]^-1*z^-1*sqrtxz1
       - 3/2*ln(x)^2*[1-x]^-1*[1+z]^-1*sqrtxz1^-1
       + 3*ln(x)^2*[1-x]^-1*sqrtxz1^-1
       + ln(x)^2*[1-x]^-1
       - 3/2*ln(x)^2*[1-x]^-1*z*[1+z]^-1*sqrtxz1^-1
       + 3/2*ln(x)^2*[1-x]^-1*z*sqrtxz1^-1
       - 1/2*ln(x)^2*[1-x]^-1*z
       - 3/2*ln(x)^2*[1+x]^-1*[1-z]^-1
       + 2*ln(x)^2*[1+x]^-1
       + 2*ln(x)^2*[1+x]^-1*z
       - 3/2*ln(x)^2*z^-1*[1-z]^-1
       - 3/2*ln(x)^2*z^-1*[1+z]^-1*sqrtxz1
       + 7/4*ln(x)^2*z^-1
       + 3/2*ln(x)^2*z^-1*sqrtxz1
       + 3/2*ln(x)^2*[1-z]^-1
       + 3/2*ln(x)^2*[1+z]^-1*sqrtxz1^-1
       - 3*ln(x)^2*sqrtxz1^-1
       + 3/2*ln(x)^2*z*[1+z]^-1*sqrtxz1^-1
       - 3/2*ln(x)^2*z*sqrtxz1^-1
       + 1/4*ln(x)^2*x*z^-1
       + 6*ln(x)^2*x*[1+z]^-1*sqrtxz1^-1
       - 6*ln(x)^2*x*sqrtxz1^-1
       + ln(x)*ln([1-x])
       - 2*ln(x)*ln([1-x])*x^-1*[1+x]^-1*z^-1*[1-z]^-1
       + 2*ln(x)*ln([1-x])*x^-1*[1+x]^-1*z^-1
       + 2*ln(x)*ln([1-x])*x^-1*[1+x]^-1
       + 2*ln(x)*ln([1-x])*x^-1*[1+x]^-1*z
       + 2*ln(x)*ln([1-x])*x^-1*z^-1*[1-z]^-1
       - 2*ln(x)*ln([1-x])*x^-1*z^-1
       - 2*ln(x)*ln([1-x])*x^-1
       - 2*ln(x)*ln([1-x])*x^-1*z
       + ln(x)*ln([1-x])*[1-x]^-1*z^-1
       - 2*ln(x)*ln([1-x])*[1-x]^-1
       + ln(x)*ln([1-x])*[1-x]^-1*z
       - 2*ln(x)*ln([1-x])*[1+x]^-1*[1-z]^-1
       + 2*ln(x)*ln([1-x])*[1+x]^-1
       + 2*ln(x)*ln([1-x])*[1+x]^-1*z
       - 2*ln(x)*ln([1-x])*z^-1*[1-z]^-1
       + 3/2*ln(x)*ln([1-x])*z^-1
       + 2*ln(x)*ln([1-x])*[1-z]^-1
       - ln(x)*ln([1-x])*z
       - 1/2*ln(x)*ln([1-x])*x*z^-1
       + ln(x)*ln([1-x])*x
       - 2*ln(x)*ln([1+x])
       + 2*ln(x)*ln([1+x])*x^-1*[1+x]^-1*z^-1*[1-z]^-1
       - 2*ln(x)*ln([1+x])*x^-1*[1+x]^-1*z^-1
       - 2*ln(x)*ln([1+x])*x^-1*[1+x]^-1
       - 2*ln(x)*ln([1+x])*x^-1*[1+x]^-1*z
       - 2*ln(x)*ln([1+x])*x^-1*z^-1*[1-z]^-1
       + 2*ln(x)*ln([1+x])*x^-1*z^-1
       + 2*ln(x)*ln([1+x])*x^-1
       + 2*ln(x)*ln([1+x])*x^-1*z
       + 2*ln(x)*ln([1+x])*z^-1*[1-z]^-1
       - 2*ln(x)*ln([1+x])*z^-1
       - 2*ln(x)*ln([1+x])*z
       + ln(x)*ln(z)
       + ln(x)*ln(z)*x^-1*[1+x]^-1*z^-1*[1-z]^-1
       - ln(x)*ln(z)*x^-1*[1+x]^-1*z^-1
       - ln(x)*ln(z)*x^-1*[1+x]^-1
       - ln(x)*ln(z)*x^-1*[1+x]^-1*z
       - ln(x)*ln(z)*x^-1*z^-1*[1-z]^-1
       + ln(x)*ln(z)*x^-1*z^-1
       + ln(x)*ln(z)*x^-1
       + ln(x)*ln(z)*x^-1*z
       - ln(x)*ln(z)*[1-x]^-1*z^-1*[1+z]^-1
       + ln(x)*ln(z)*[1-x]^-1*z^-1*[1+z]^-1*sqrtxz1
       + ln(x)*ln(z)*[1-x]^-1*z^-1
       - ln(x)*ln(z)*[1-x]^-1*z^-1*sqrtxz1
       - ln(x)*ln(z)*[1-x]^-1*[1+z]^-1*sqrtxz1^-1
       + 2*ln(x)*ln(z)*[1-x]^-1*sqrtxz1^-1
       - ln(x)*ln(z)*[1-x]^-1
       - ln(x)*ln(z)*[1-x]^-1*z*[1+z]^-1*sqrtxz1^-1
       + ln(x)*ln(z)*[1-x]^-1*z*sqrtxz1^-1
       + ln(x)*ln(z)*[1-x]^-1*z
       + ln(x)*ln(z)*z^-1*[1-z]^-1
       + ln(x)*ln(z)*z^-1*[1+z]^-1
       - ln(x)*ln(z)*z^-1*[1+z]^-1*sqrtxz1
       - 2*ln(x)*ln(z)*z^-1
       + ln(x)*ln(z)*z^-1*sqrtxz1
       - 1/2*ln(x)*ln(z)*[1-z]^-1
       + ln(x)*ln(z)*[1+z]^-1*sqrtxz1^-1
       - 2*ln(x)*ln(z)*sqrtxz1^-1
       + ln(x)*ln(z)*z*[1+z]^-1*sqrtxz1^-1
       - ln(x)*ln(z)*z*sqrtxz1^-1
       - 2*ln(x)*ln(z)*z
       - 1/2*ln(x)*ln(z)*x*[1-z]^-1
       + 4*ln(x)*ln(z)*x*[1+z]^-1*sqrtxz1^-1
       - 4*ln(x)*ln(z)*x*sqrtxz1^-1
       + ln(x)*ln(z)*x
       - 2*ln(x)*ln([1-z])*[1-x]^-1*z^-1*[1+z]^-1*sqrtxz1
       + 2*ln(x)*ln([1-z])*[1-x]^-1*z^-1*sqrtxz1
       + 2*ln(x)*ln([1-z])*[1-x]^-1*[1+z]^-1*sqrtxz1^-1
       - 4*ln(x)*ln([1-z])*[1-x]^-1*sqrtxz1^-1
       + 2*ln(x)*ln([1-z])*[1-x]^-1*z*[1+z]^-1*sqrtxz1^-1
       - 2*ln(x)*ln([1-z])*[1-x]^-1*z*sqrtxz1^-1
       + 2*ln(x)*ln([1-z])*z^-1*[1+z]^-1*sqrtxz1
       - 2*ln(x)*ln([1-z])*z^-1*sqrtxz1
       - 2*ln(x)*ln([1-z])*[1+z]^-1*sqrtxz1^-1
       + 4*ln(x)*ln([1-z])*sqrtxz1^-1
       - 2*ln(x)*ln([1-z])*z*[1+z]^-1*sqrtxz1^-1
       + 2*ln(x)*ln([1-z])*z*sqrtxz1^-1
       - 8*ln(x)*ln([1-z])*x*[1+z]^-1*sqrtxz1^-1
       + 8*ln(x)*ln([1-z])*x*sqrtxz1^-1
       + 3/2*ln(z)
       - 6*ln(z)*[1-x]^-1*z^-1*[1+z]^-1*rln2
       - 4*ln(z)*[1-x]^-1*z^-1*[1+z]^-1*sqrtxz1
       - 8*ln(z)*[1-x]^-1*z^-1*[1+z]^-1*sqrtxz1*rln2
       + 6*ln(z)*[1-x]^-1*z^-1*rln2
       + 4*ln(z)*[1-x]^-1*z^-1*sqrtxz1
       + 8*ln(z)*[1-x]^-1*z^-1*sqrtxz1*rln2
       + 4*ln(z)*[1-x]^-1*[1+z]^-1*sqrtxz1^-1
       + 8*ln(z)*[1-x]^-1*[1+z]^-1*sqrtxz1^-1*rln2
       - 8*ln(z)*[1-x]^-1*sqrtxz1^-1
       - 16*ln(z)*[1-x]^-1*sqrtxz1^-1*rln2
       - 1/2*ln(z)*[1-x]^-1
       - 9/2*ln(z)*[1-x]^-1*rln2
       + 1/2*ln(z)*[1-x]^-1*sqrtxz1
       + 4*ln(z)*[1-x]^-1*z*[1+z]^-1*sqrtxz1^-1
       + 8*ln(z)*[1-x]^-1*z*[1+z]^-1*sqrtxz1^-1*rln2
       - 4*ln(z)*[1-x]^-1*z*sqrtxz1^-1
       - 8*ln(z)*[1-x]^-1*z*sqrtxz1^-1*rln2
       - 1/2*ln(z)*[1-x]^-1*z
       + 9/2*ln(z)*[1-x]^-1*z*rln2
       + 6*ln(z)*z^-1*[1+z]^-1*rln2
       + 4*ln(z)*z^-1*[1+z]^-1*sqrtxz1
       + 8*ln(z)*z^-1*[1+z]^-1*sqrtxz1*rln2
       - 6*ln(z)*z^-1*rln2
       - 4*ln(z)*z^-1*sqrtxz1
       - 8*ln(z)*z^-1*sqrtxz1*rln2
       - 1/2*ln(z)*[1-z]^-1
       - 4*ln(z)*[1+z]^-1*sqrtxz1^-1
       - 8*ln(z)*[1+z]^-1*sqrtxz1^-1*rln2
       + 8*ln(z)*sqrtxz1^-1
       + 16*ln(z)*sqrtxz1^-1*rln2
       + 6*ln(z)*rln2
       - 4*ln(z)*z*[1+z]^-1*sqrtxz1^-1
       - 8*ln(z)*z*[1+z]^-1*sqrtxz1^-1*rln2
       + 4*ln(z)*z*sqrtxz1^-1
       + 8*ln(z)*z*sqrtxz1^-1*rln2
       - 6*ln(z)*z*rln2
       + 1/2*ln(z)*x*[1-z]^-1
       - 16*ln(z)*x*[1+z]^-1*sqrtxz1^-1
       - 32*ln(z)*x*[1+z]^-1*sqrtxz1^-1*rln2
       + 16*ln(z)*x*sqrtxz1^-1
       + 32*ln(z)*x*sqrtxz1^-1*rln2
       - 3/2*ln(z)*x
       - 2*ln(z)*ln(1 + sqrtxz1 - z)
       + 2*ln(z)*ln(1 + sqrtxz1 - z)*[1-x]^-1*z^-1*[1+z]^-1
       + 4*ln(z)*ln(1 + sqrtxz1 - z)*[1-x]^-1*z^-1*[1+z]^-1*sqrtxz1
       - 2*ln(z)*ln(1 + sqrtxz1 - z)*[1-x]^-1*z^-1
       - 4*ln(z)*ln(1 + sqrtxz1 - z)*[1-x]^-1*z^-1*sqrtxz1
       - 4*ln(z)*ln(1 + sqrtxz1 - z)*[1-x]^-1*[1+z]^-1*sqrtxz1^-1
       + 8*ln(z)*ln(1 + sqrtxz1 - z)*[1-x]^-1*sqrtxz1^-1
       + 3/2*ln(z)*ln(1 + sqrtxz1 - z)*[1-x]^-1
       - 4*ln(z)*ln(1 + sqrtxz1 - z)*[1-x]^-1*z*[1+z]^-1*sqrtxz1^-1
       + 4*ln(z)*ln(1 + sqrtxz1 - z)*[1-x]^-1*z*sqrtxz1^-1
       - 3/2*ln(z)*ln(1 + sqrtxz1 - z)*[1-x]^-1*z
       - 2*ln(z)*ln(1 + sqrtxz1 - z)*z^-1*[1+z]^-1
       - 4*ln(z)*ln(1 + sqrtxz1 - z)*z^-1*[1+z]^-1*sqrtxz1
       + 2*ln(z)*ln(1 + sqrtxz1 - z)*z^-1
       + 4*ln(z)*ln(1 + sqrtxz1 - z)*z^-1*sqrtxz1
       + 4*ln(z)*ln(1 + sqrtxz1 - z)*[1+z]^-1*sqrtxz1^-1
       - 8*ln(z)*ln(1 + sqrtxz1 - z)*sqrtxz1^-1
       + 4*ln(z)*ln(1 + sqrtxz1 - z)*z*[1+z]^-1*sqrtxz1^-1
       - 4*ln(z)*ln(1 + sqrtxz1 - z)*z*sqrtxz1^-1
       + 2*ln(z)*ln(1 + sqrtxz1 - z)*z
       + 16*ln(z)*ln(1 + sqrtxz1 - z)*x*[1+z]^-1*sqrtxz1^-1
       - 16*ln(z)*ln(1 + sqrtxz1 - z)*x*sqrtxz1^-1
       - 4*ln(z)*ln(1 + sqrtxz1 + z)
       + 4*ln(z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^-1*[1+z]^-1
       + 4*ln(z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^-1*[1+z]^-1*sqrtxz1
       - 4*ln(z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^-1
       - 4*ln(z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^-1*sqrtxz1
       - 4*ln(z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*[1+z]^-1*sqrtxz1^-1
       + 8*ln(z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*sqrtxz1^-1
       + 3*ln(z)*ln(1 + sqrtxz1 + z)*[1-x]^-1
       - 4*ln(z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z*[1+z]^-1*sqrtxz1^-1
       + 4*ln(z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z*sqrtxz1^-1
       - 3*ln(z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z
       - 4*ln(z)*ln(1 + sqrtxz1 + z)*z^-1*[1+z]^-1
       - 4*ln(z)*ln(1 + sqrtxz1 + z)*z^-1*[1+z]^-1*sqrtxz1
       + 4*ln(z)*ln(1 + sqrtxz1 + z)*z^-1
       + 4*ln(z)*ln(1 + sqrtxz1 + z)*z^-1*sqrtxz1
       + 4*ln(z)*ln(1 + sqrtxz1 + z)*[1+z]^-1*sqrtxz1^-1
       - 8*ln(z)*ln(1 + sqrtxz1 + z)*sqrtxz1^-1
       + 4*ln(z)*ln(1 + sqrtxz1 + z)*z*[1+z]^-1*sqrtxz1^-1
       - 4*ln(z)*ln(1 + sqrtxz1 + z)*z*sqrtxz1^-1
       + 4*ln(z)*ln(1 + sqrtxz1 + z)*z
       + 16*ln(z)*ln(1 + sqrtxz1 + z)*x*[1+z]^-1*sqrtxz1^-1
       - 16*ln(z)*ln(1 + sqrtxz1 + z)*x*sqrtxz1^-1
       + 2*ln(z)*ln(1 + z)*[1-x]^-1*[1+z]^-1
       - ln(z)*ln(1 + z)*[1-x]^-1
       + ln(z)*ln(1 + z)*[1-x]^-1*z
       - ln(z)*ln(1 + x*z^-1)
       + ln(z)*ln(1 + x*z^-1)*x^-1*[1+x]^-1*z^-1*[1-z]^-1
       - ln(z)*ln(1 + x*z^-1)*x^-1*[1+x]^-1*z^-1
       - ln(z)*ln(1 + x*z^-1)*x^-1*[1+x]^-1
       - ln(z)*ln(1 + x*z^-1)*x^-1*[1+x]^-1*z
       - ln(z)*ln(1 + x*z^-1)*x^-1*z^-1*[1-z]^-1
       + ln(z)*ln(1 + x*z^-1)*x^-1*z^-1
       + ln(z)*ln(1 + x*z^-1)*x^-1
       + ln(z)*ln(1 + x*z^-1)*x^-1*z
       + ln(z)*ln(1 + x*z^-1)*z^-1*[1-z]^-1
       - ln(z)*ln(1 + x*z^-1)*z^-1
       - ln(z)*ln(1 + x*z^-1)*z
       + 2*ln(z)*ln(1 + x*z)
       - ln(z)*ln(1 + x*z)*x^-1*[1+x]^-1*z^-1*[1-z]^-1
       + ln(z)*ln(1 + x*z)*x^-1*[1+x]^-1*z^-1
       + ln(z)*ln(1 + x*z)*x^-1*[1+x]^-1
       + ln(z)*ln(1 + x*z)*x^-1*[1+x]^-1*z
       + ln(z)*ln(1 + x*z)*x^-1*z^-1*[1-z]^-1
       - ln(z)*ln(1 + x*z)*x^-1*z^-1
       - ln(z)*ln(1 + x*z)*x^-1
       - ln(z)*ln(1 + x*z)*x^-1*z
       - ln(z)*ln(1 + x*z)*[1-x]^-1*z^-1*[1+z]^-1
       + ln(z)*ln(1 + x*z)*[1-x]^-1*z^-1
       - ln(z)*ln(1 + x*z)*[1-x]^-1
       + ln(z)*ln(1 + x*z)*[1-x]^-1*z
       - ln(z)*ln(1 + x*z)*z^-1*[1-z]^-1
       + ln(z)*ln(1 + x*z)*z^-1*[1+z]^-1
       + ln(z)*ln(z + x)
       - ln(z)*ln(z + x)*[1-x]^-1*z^-1*[1+z]^-1
       + ln(z)*ln(z + x)*[1-x]^-1*z^-1
       - ln(z)*ln(z + x)*[1-x]^-1
       + ln(z)*ln(z + x)*[1-x]^-1*z
       + ln(z)*ln(z + x)*z^-1*[1+z]^-1
       - ln(z)*ln(z + x)*z^-1
       - ln(z)*ln(z + x)*z
       + 1/2*ln(z)^2
       + 1/2*ln(z)^2*x^-1*[1+x]^-1*z^-1*[1-z]^-1
       - 1/2*ln(z)^2*x^-1*[1+x]^-1*z^-1
       - 1/2*ln(z)^2*x^-1*[1+x]^-1
       - 1/2*ln(z)^2*x^-1*[1+x]^-1*z
       - 1/2*ln(z)^2*x^-1*z^-1*[1-z]^-1
       + 1/2*ln(z)^2*x^-1*z^-1
       + 1/2*ln(z)^2*x^-1
       + 1/2*ln(z)^2*x^-1*z
       - ln(z)^2*[1-x]^-1*z^-1*[1+z]^-1
       - 5/2*ln(z)^2*[1-x]^-1*z^-1*[1+z]^-1*sqrtxz1
       + ln(z)^2*[1-x]^-1*z^-1
       + 5/2*ln(z)^2*[1-x]^-1*z^-1*sqrtxz1
       + 5/2*ln(z)^2*[1-x]^-1*[1+z]^-1*sqrtxz1^-1
       - 1/2*ln(z)^2*[1-x]^-1*[1+z]^-1
       - 5*ln(z)^2*[1-x]^-1*sqrtxz1^-1
       - 1/4*ln(z)^2*[1-x]^-1
       + 5/2*ln(z)^2*[1-x]^-1*z*[1+z]^-1*sqrtxz1^-1
       - 5/2*ln(z)^2*[1-x]^-1*z*sqrtxz1^-1
       + 1/4*ln(z)^2*[1-x]^-1*z
       + 1/2*ln(z)^2*z^-1*[1-z]^-1
       + ln(z)^2*z^-1*[1+z]^-1
       + 5/2*ln(z)^2*z^-1*[1+z]^-1*sqrtxz1
       - 3/2*ln(z)^2*z^-1
       - 5/2*ln(z)^2*z^-1*sqrtxz1
       - 5/2*ln(z)^2*[1+z]^-1*sqrtxz1^-1
       + 5*ln(z)^2*sqrtxz1^-1
       - 5/2*ln(z)^2*z*[1+z]^-1*sqrtxz1^-1
       + 5/2*ln(z)^2*z*sqrtxz1^-1
       - 3/2*ln(z)^2*z
       - 10*ln(z)^2*x*[1+z]^-1*sqrtxz1^-1
       + 10*ln(z)^2*x*sqrtxz1^-1
       - 2*ln(z)*ln([1-z])*[1-x]^-1*z^-1*[1+z]^-1*sqrtxz1
       + 2*ln(z)*ln([1-z])*[1-x]^-1*z^-1*sqrtxz1
       + 2*ln(z)*ln([1-z])*[1-x]^-1*[1+z]^-1*sqrtxz1^-1
       - 4*ln(z)*ln([1-z])*[1-x]^-1*sqrtxz1^-1
       + 2*ln(z)*ln([1-z])*[1-x]^-1*z*[1+z]^-1*sqrtxz1^-1
       - 2*ln(z)*ln([1-z])*[1-x]^-1*z*sqrtxz1^-1
       + 2*ln(z)*ln([1-z])*z^-1*[1+z]^-1*sqrtxz1
       - 2*ln(z)*ln([1-z])*z^-1*sqrtxz1
       - 2*ln(z)*ln([1-z])*[1+z]^-1*sqrtxz1^-1
       + 4*ln(z)*ln([1-z])*sqrtxz1^-1
       - 2*ln(z)*ln([1-z])*z*[1+z]^-1*sqrtxz1^-1
       + 2*ln(z)*ln([1-z])*z*sqrtxz1^-1
       - 8*ln(z)*ln([1-z])*x*[1+z]^-1*sqrtxz1^-1
       + 8*ln(z)*ln([1-z])*x*sqrtxz1^-1
       - 4*ln([1-z])*[1-x]^-1*z^-1*[1+z]^-1*sqrtxz1*rln2
       + 4*ln([1-z])*[1-x]^-1*z^-1*sqrtxz1*rln2
       + 4*ln([1-z])*[1-x]^-1*[1+z]^-1*sqrtxz1^-1*rln2
       - 8*ln([1-z])*[1-x]^-1*sqrtxz1^-1*rln2
       + 4*ln([1-z])*[1-x]^-1*z*[1+z]^-1*sqrtxz1^-1*rln2
       - 4*ln([1-z])*[1-x]^-1*z*sqrtxz1^-1*rln2
       + 4*ln([1-z])*z^-1*[1+z]^-1*sqrtxz1*rln2
       - 4*ln([1-z])*z^-1*sqrtxz1*rln2
       - 4*ln([1-z])*[1+z]^-1*sqrtxz1^-1*rln2
       + 8*ln([1-z])*sqrtxz1^-1*rln2
       - 4*ln([1-z])*z*[1+z]^-1*sqrtxz1^-1*rln2
       + 4*ln([1-z])*z*sqrtxz1^-1*rln2
       - 16*ln([1-z])*x*[1+z]^-1*sqrtxz1^-1*rln2
       + 16*ln([1-z])*x*sqrtxz1^-1*rln2
       + 4*ln([1-z])*ln(1 + sqrtxz1 - z)*[1-x]^-1*z^-1*[1+z]^-1*sqrtxz1
       - 4*ln([1-z])*ln(1 + sqrtxz1 - z)*[1-x]^-1*z^-1*sqrtxz1
       - 4*ln([1-z])*ln(1 + sqrtxz1 - z)*[1-x]^-1*[1+z]^-1*sqrtxz1^-1
       + 8*ln([1-z])*ln(1 + sqrtxz1 - z)*[1-x]^-1*sqrtxz1^-1
       - 4*ln([1-z])*ln(1 + sqrtxz1 - z)*[1-x]^-1*z*[1+z]^-1*sqrtxz1^-1
       + 4*ln([1-z])*ln(1 + sqrtxz1 - z)*[1-x]^-1*z*sqrtxz1^-1
       - 4*ln([1-z])*ln(1 + sqrtxz1 - z)*z^-1*[1+z]^-1*sqrtxz1
       + 4*ln([1-z])*ln(1 + sqrtxz1 - z)*z^-1*sqrtxz1
       + 4*ln([1-z])*ln(1 + sqrtxz1 - z)*[1+z]^-1*sqrtxz1^-1
       - 8*ln([1-z])*ln(1 + sqrtxz1 - z)*sqrtxz1^-1
       + 4*ln([1-z])*ln(1 + sqrtxz1 - z)*z*[1+z]^-1*sqrtxz1^-1
       - 4*ln([1-z])*ln(1 + sqrtxz1 - z)*z*sqrtxz1^-1
       + 16*ln([1-z])*ln(1 + sqrtxz1 - z)*x*[1+z]^-1*sqrtxz1^-1
       - 16*ln([1-z])*ln(1 + sqrtxz1 - z)*x*sqrtxz1^-1
       + 2*ln([1-z])*ln(1 - 2*z + z^2 + 4*x*z)*[1-x]^-1*z^-1*[1+z]^-1*sqrtxz1
       - 2*ln([1-z])*ln(1 - 2*z + z^2 + 4*x*z)*[1-x]^-1*z^-1*sqrtxz1
       - 2*ln([1-z])*ln(1 - 2*z + z^2 + 4*x*z)*[1-x]^-1*[1+z]^-1*sqrtxz1^-1
       + 4*ln([1-z])*ln(1 - 2*z + z^2 + 4*x*z)*[1-x]^-1*sqrtxz1^-1
       - 2*ln([1-z])*ln(1 - 2*z + z^2 + 4*x*z)*[1-x]^-1*z*[1+z]^-1*sqrtxz1^-1
       + 2*ln([1-z])*ln(1 - 2*z + z^2 + 4*x*z)*[1-x]^-1*z*sqrtxz1^-1
       - 2*ln([1-z])*ln(1 - 2*z + z^2 + 4*x*z)*z^-1*[1+z]^-1*sqrtxz1
       + 2*ln([1-z])*ln(1 - 2*z + z^2 + 4*x*z)*z^-1*sqrtxz1
       + 2*ln([1-z])*ln(1 - 2*z + z^2 + 4*x*z)*[1+z]^-1*sqrtxz1^-1
       - 4*ln([1-z])*ln(1 - 2*z + z^2 + 4*x*z)*sqrtxz1^-1
       + 2*ln([1-z])*ln(1 - 2*z + z^2 + 4*x*z)*z*[1+z]^-1*sqrtxz1^-1
       - 2*ln([1-z])*ln(1 - 2*z + z^2 + 4*x*z)*z*sqrtxz1^-1
       + 8*ln([1-z])*ln(1 - 2*z + z^2 + 4*x*z)*x*[1+z]^-1*sqrtxz1^-1
       - 8*ln([1-z])*ln(1 - 2*z + z^2 + 4*x*z)*x*sqrtxz1^-1
       - 2*ln([1-z])^2*[1-x]^-1*z^-1*[1+z]^-1*sqrtxz1
       + 2*ln([1-z])^2*[1-x]^-1*z^-1*sqrtxz1
       + 2*ln([1-z])^2*[1-x]^-1*[1+z]^-1*sqrtxz1^-1
       - 4*ln([1-z])^2*[1-x]^-1*sqrtxz1^-1
       + 2*ln([1-z])^2*[1-x]^-1*z*[1+z]^-1*sqrtxz1^-1
       - 2*ln([1-z])^2*[1-x]^-1*z*sqrtxz1^-1
       + 2*ln([1-z])^2*z^-1*[1+z]^-1*sqrtxz1
       - 2*ln([1-z])^2*z^-1*sqrtxz1
       - 2*ln([1-z])^2*[1+z]^-1*sqrtxz1^-1
       + 4*ln([1-z])^2*sqrtxz1^-1
       - 2*ln([1-z])^2*z*[1+z]^-1*sqrtxz1^-1
       + 2*ln([1-z])^2*z*sqrtxz1^-1
       - 8*ln([1-z])^2*x*[1+z]^-1*sqrtxz1^-1
       + 8*ln([1-z])^2*x*sqrtxz1^-1
       + 2*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)
       - 2*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*z^-1*[1+z]^-1
       - 2*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*z^-1*[1+z]^-1*sqrtxz1
       + 2*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*z^-1
       + 2*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*z^-1*sqrtxz1
       + 2*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*[1+z]^-1*sqrtxz1^-1
       - 4*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*sqrtxz1^-1
       - 3/2*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1
       + 2*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*z*[1+z]^-1*sqrtxz1^-1
       - 2*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*z*sqrtxz1^-1
       + 3/2*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*z
       + 2*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z^-1*[1+z]^-1
       + 2*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z^-1*[1+z]^-1*sqrtxz1
       - 2*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z^-1
       - 2*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z^-1*sqrtxz1
       - 2*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1+z]^-1*sqrtxz1^-1
       + 4*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*sqrtxz1^-1
       - 2*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z*[1+z]^-1*sqrtxz1^-1
       + 2*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z*sqrtxz1^-1
       - 2*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z
       - 8*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x*[1+z]^-1*sqrtxz1^-1
       + 8*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x*sqrtxz1^-1
       - 2*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)
       + 2*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*z^-1*[1+z]^-1
       - 2*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*z^-1*[1+z]^-1*sqrtxz1
       - 2*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*z^-1
       + 2*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*z^-1*sqrtxz1
       + 2*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*[1+z]^-1*sqrtxz1^-1
       - 4*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*sqrtxz1^-1
       + 3/2*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1
       + 2*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*z*[1+z]^-1*sqrtxz1^-1
       - 2*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*z*sqrtxz1^-1
       - 3/2*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*z
       - 2*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z^-1*[1+z]^-1
       + 2*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z^-1*[1+z]^-1*sqrtxz1
       + 2*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z^-1
       - 2*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z^-1*sqrtxz1
       - 2*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1+z]^-1*sqrtxz1^-1
       + 4*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*sqrtxz1^-1
       - 2*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z*[1+z]^-1*sqrtxz1^-1
       + 2*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z*sqrtxz1^-1
       + 2*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z
       - 8*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x*[1+z]^-1*sqrtxz1^-1
       + 8*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x*sqrtxz1^-1
       - 2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)
       + 2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*[1-x]^-1*z^-1*[1+z]^-1
       + 2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*[1-x]^-1*z^-1*[1+z]^-1*sqrtxz1
       - 2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*[1-x]^-1*z^-1
       - 2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*[1-x]^-1*z^-1*sqrtxz1
       - 2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*[1-x]^-1*[1+z]^-1*sqrtxz1^-1
       + 4*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*[1-x]^-1*sqrtxz1^-1
       + 3/2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*[1-x]^-1
       - 2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*[1-x]^-1*z*[1+z]^-1*sqrtxz1^-1
       + 2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*[1-x]^-1*z*sqrtxz1^-1
       - 3/2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*[1-x]^-1*z
       - 2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*z^-1*[1+z]^-1
       - 2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*z^-1*[1+z]^-1*sqrtxz1
       + 2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*z^-1
       + 2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*z^-1*sqrtxz1
       + 2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*[1+z]^-1*sqrtxz1^-1
       - 4*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*sqrtxz1^-1
       + 2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*z*[1+z]^-1*sqrtxz1^-1
       - 2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*z*sqrtxz1^-1
       + 2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*z
       + 8*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*x*[1+z]^-1*sqrtxz1^-1
       - 8*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*x*sqrtxz1^-1
       + 2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)
       - 2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*[1-x]^-1*z^-1*[1+z]^-1
       + 2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*[1-x]^-1*z^-1*[1+z]^-1*sqrtxz1
       + 2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*[1-x]^-1*z^-1
       - 2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*[1-x]^-1*z^-1*sqrtxz1
       - 2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*[1-x]^-1*[1+z]^-1*sqrtxz1^-1
       + 4*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*[1-x]^-1*sqrtxz1^-1
       - 3/2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*[1-x]^-1
       - 2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*[1-x]^-1*z*[1+z]^-1*sqrtxz1^-1
       + 2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*[1-x]^-1*z*sqrtxz1^-1
       + 3/2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*[1-x]^-1*z
       + 2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*z^-1*[1+z]^-1
       - 2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*z^-1*[1+z]^-1*sqrtxz1
       - 2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*z^-1
       + 2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*z^-1*sqrtxz1
       + 2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*[1+z]^-1*sqrtxz1^-1
       - 4*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*sqrtxz1^-1
       + 2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*z*[1+z]^-1*sqrtxz1^-1
       - 2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*z*sqrtxz1^-1
       - 2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*z
       + 8*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*x*[1+z]^-1*sqrtxz1^-1
       - 8*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*x*sqrtxz1^-1
       - 4*Li2(sqrtxz1*[z-1]^-1)*[1-x]^-1*z^-1*[1+z]^-1*sqrtxz1
       + 4*Li2(sqrtxz1*[z-1]^-1)*[1-x]^-1*z^-1*sqrtxz1
       + 4*Li2(sqrtxz1*[z-1]^-1)*[1-x]^-1*[1+z]^-1*sqrtxz1^-1
       - 8*Li2(sqrtxz1*[z-1]^-1)*[1-x]^-1*sqrtxz1^-1
       + 4*Li2(sqrtxz1*[z-1]^-1)*[1-x]^-1*z*[1+z]^-1*sqrtxz1^-1
       - 4*Li2(sqrtxz1*[z-1]^-1)*[1-x]^-1*z*sqrtxz1^-1
       + 4*Li2(sqrtxz1*[z-1]^-1)*z^-1*[1+z]^-1*sqrtxz1
       - 4*Li2(sqrtxz1*[z-1]^-1)*z^-1*sqrtxz1
       - 4*Li2(sqrtxz1*[z-1]^-1)*[1+z]^-1*sqrtxz1^-1
       + 8*Li2(sqrtxz1*[z-1]^-1)*sqrtxz1^-1
       - 4*Li2(sqrtxz1*[z-1]^-1)*z*[1+z]^-1*sqrtxz1^-1
       + 4*Li2(sqrtxz1*[z-1]^-1)*z*sqrtxz1^-1
       - 16*Li2(sqrtxz1*[z-1]^-1)*x*[1+z]^-1*sqrtxz1^-1
       + 16*Li2(sqrtxz1*[z-1]^-1)*x*sqrtxz1^-1
       - 4*Li2([1-z]*sqrtxz1^-1)*[1-x]^-1*z^-1*[1+z]^-1*sqrtxz1
       + 4*Li2([1-z]*sqrtxz1^-1)*[1-x]^-1*z^-1*sqrtxz1
       + 4*Li2([1-z]*sqrtxz1^-1)*[1-x]^-1*[1+z]^-1*sqrtxz1^-1
       - 8*Li2([1-z]*sqrtxz1^-1)*[1-x]^-1*sqrtxz1^-1
       + 4*Li2([1-z]*sqrtxz1^-1)*[1-x]^-1*z*[1+z]^-1*sqrtxz1^-1
       - 4*Li2([1-z]*sqrtxz1^-1)*[1-x]^-1*z*sqrtxz1^-1
       + 4*Li2([1-z]*sqrtxz1^-1)*z^-1*[1+z]^-1*sqrtxz1
       - 4*Li2([1-z]*sqrtxz1^-1)*z^-1*sqrtxz1
       - 4*Li2([1-z]*sqrtxz1^-1)*[1+z]^-1*sqrtxz1^-1
       + 8*Li2([1-z]*sqrtxz1^-1)*sqrtxz1^-1
       - 4*Li2([1-z]*sqrtxz1^-1)*z*[1+z]^-1*sqrtxz1^-1
       + 4*Li2([1-z]*sqrtxz1^-1)*z*sqrtxz1^-1
       - 16*Li2([1-z]*sqrtxz1^-1)*x*[1+z]^-1*sqrtxz1^-1
       + 16*Li2([1-z]*sqrtxz1^-1)*x*sqrtxz1^-1
       + 2*Li2( - z)*[1-x]^-1*[1+z]^-1
       - Li2( - z)*[1-x]^-1
       + Li2( - z)*[1-x]^-1*z
       - Li2( - x*z^-1)*x^-1*[1+x]^-1*z^-1*[1-z]^-1
       + Li2( - x*z^-1)*x^-1*[1+x]^-1*z^-1
       + Li2( - x*z^-1)*x^-1*[1+x]^-1
       + Li2( - x*z^-1)*x^-1*[1+x]^-1*z
       + Li2( - x*z^-1)*x^-1*z^-1*[1-z]^-1
       - Li2( - x*z^-1)*x^-1*z^-1
       - Li2( - x*z^-1)*x^-1
       - Li2( - x*z^-1)*x^-1*z
       + Li2( - x*z^-1)*[1-x]^-1*z^-1*[1+z]^-1
       - Li2( - x*z^-1)*[1-x]^-1*z^-1
       + Li2( - x*z^-1)*[1-x]^-1
       - Li2( - x*z^-1)*[1-x]^-1*z
       - Li2( - x*z^-1)*z^-1*[1-z]^-1
       - Li2( - x*z^-1)*z^-1*[1+z]^-1
       + 2*Li2( - x*z^-1)*z^-1
       + 2*Li2( - x*z^-1)*z
       - 2*Li2( - x)
       + 2*Li2( - x)*x^-1*[1+x]^-1*z^-1*[1-z]^-1
       - 2*Li2( - x)*x^-1*[1+x]^-1*z^-1
       - 2*Li2( - x)*x^-1*[1+x]^-1
       - 2*Li2( - x)*x^-1*[1+x]^-1*z
       - 2*Li2( - x)*x^-1*z^-1*[1-z]^-1
       + 2*Li2( - x)*x^-1*z^-1
       + 2*Li2( - x)*x^-1
       + 2*Li2( - x)*x^-1*z
       - 2*Li2( - x)*[1+x]^-1
       - 2*Li2( - x)*[1+x]^-1*z
       + 2*Li2( - x)*z^-1*[1-z]^-1
       - 2*Li2( - x)*z^-1
       - 2*Li2( - x)*x
       + 2*Li2( - x*z)
       - Li2( - x*z)*x^-1*[1+x]^-1*z^-1*[1-z]^-1
       + Li2( - x*z)*x^-1*[1+x]^-1*z^-1
       + Li2( - x*z)*x^-1*[1+x]^-1
       + Li2( - x*z)*x^-1*[1+x]^-1*z
       + Li2( - x*z)*x^-1*z^-1*[1-z]^-1
       - Li2( - x*z)*x^-1*z^-1
       - Li2( - x*z)*x^-1
       - Li2( - x*z)*x^-1*z
       - Li2( - x*z)*[1-x]^-1*z^-1*[1+z]^-1
       + Li2( - x*z)*[1-x]^-1*z^-1
       - Li2( - x*z)*[1-x]^-1
       + Li2( - x*z)*[1-x]^-1*z
       - Li2( - x*z)*z^-1*[1-z]^-1
       + Li2( - x*z)*z^-1*[1+z]^-1
       - 2*Li2( - 1/(1 + sqrtxz1 - z) + 1/(1 + sqrtxz1 - z)*sqrtxz1 + 1/(1 + sqrtxz1 - z)*z)*[1-x]^-1*z^-1*[1+z]^-1*sqrtxz1
       + 2*Li2( - 1/(1 + sqrtxz1 - z) + 1/(1 + sqrtxz1 - z)*sqrtxz1 + 1/(1 + sqrtxz1 - z)*z)*[1-x]^-1*z^-1*sqrtxz1
       + 2*Li2( - 1/(1 + sqrtxz1 - z) + 1/(1 + sqrtxz1 - z)*sqrtxz1 + 1/(1 + sqrtxz1 - z)*z)*[1-x]^-1*[1+z]^-1*sqrtxz1^-1
       - 4*Li2( - 1/(1 + sqrtxz1 - z) + 1/(1 + sqrtxz1 - z)*sqrtxz1 + 1/(1 + sqrtxz1 - z)*z)*[1-x]^-1*sqrtxz1^-1
       + 2*Li2( - 1/(1 + sqrtxz1 - z) + 1/(1 + sqrtxz1 - z)*sqrtxz1 + 1/(1 + sqrtxz1 - z)*z)*[1-x]^-1*z*[1+z]^-1*sqrtxz1^-1
       - 2*Li2( - 1/(1 + sqrtxz1 - z) + 1/(1 + sqrtxz1 - z)*sqrtxz1 + 1/(1 + sqrtxz1 - z)*z)*[1-x]^-1*z*sqrtxz1^-1
       + 2*Li2( - 1/(1 + sqrtxz1 - z) + 1/(1 + sqrtxz1 - z)*sqrtxz1 + 1/(1 + sqrtxz1 - z)*z)*z^-1*[1+z]^-1*sqrtxz1
       - 2*Li2( - 1/(1 + sqrtxz1 - z) + 1/(1 + sqrtxz1 - z)*sqrtxz1 + 1/(1 + sqrtxz1 - z)*z)*z^-1*sqrtxz1
       - 2*Li2( - 1/(1 + sqrtxz1 - z) + 1/(1 + sqrtxz1 - z)*sqrtxz1 + 1/(1 + sqrtxz1 - z)*z)*[1+z]^-1*sqrtxz1^-1
       + 4*Li2( - 1/(1 + sqrtxz1 - z) + 1/(1 + sqrtxz1 - z)*sqrtxz1 + 1/(1 + sqrtxz1 - z)*z)*sqrtxz1^-1
       - 2*Li2( - 1/(1 + sqrtxz1 - z) + 1/(1 + sqrtxz1 - z)*sqrtxz1 + 1/(1 + sqrtxz1 - z)*z)*z*[1+z]^-1*sqrtxz1^-1
       + 2*Li2( - 1/(1 + sqrtxz1 - z) + 1/(1 + sqrtxz1 - z)*sqrtxz1 + 1/(1 + sqrtxz1 - z)*z)*z*sqrtxz1^-1
       - 8*Li2( - 1/(1 + sqrtxz1 - z) + 1/(1 + sqrtxz1 - z)*sqrtxz1 + 1/(1 + sqrtxz1 - z)*z)*x*[1+z]^-1*sqrtxz1^-1
       + 8*Li2( - 1/(1 + sqrtxz1 - z) + 1/(1 + sqrtxz1 - z)*sqrtxz1 + 1/(1 + sqrtxz1 - z)*z)*x*sqrtxz1^-1
       + Li2(x)
       - 2*Li2(x)*x^-1*[1+x]^-1*z^-1*[1-z]^-1
       + 2*Li2(x)*x^-1*[1+x]^-1*z^-1
       + 2*Li2(x)*x^-1*[1+x]^-1
       + 2*Li2(x)*x^-1*[1+x]^-1*z
       + 2*Li2(x)*x^-1*z^-1*[1-z]^-1
       - 2*Li2(x)*x^-1*z^-1
       - 2*Li2(x)*x^-1
       - 2*Li2(x)*x^-1*z
       + Li2(x)*[1-x]^-1*z^-1
       - 2*Li2(x)*[1-x]^-1
       + Li2(x)*[1-x]^-1*z
       - 2*Li2(x)*[1+x]^-1*[1-z]^-1
       + 2*Li2(x)*[1+x]^-1
       + 2*Li2(x)*[1+x]^-1*z
       - 2*Li2(x)*z^-1*[1-z]^-1
       + 3/2*Li2(x)*z^-1
       + 2*Li2(x)*[1-z]^-1
       - Li2(x)*z
       - 1/2*Li2(x)*x*z^-1
       + Li2(x)*x
       )

       + NC^-1 * (  - 17/6
       + 3*[1-x]^-1*rln2^2
       + 3*[1-x]^-1*sqrtxz1*rln2
       + 2*[1-x]^-1*z*rln2^2
       + 8*[1-x]^-1*z^2*rln2^2
       + 113/72*z^-1
       - 2*rln2^2
       - 2*sqrtxz1*rln2
       - 85/24*z
       + 65/36*z^2
       - 8*z^2*rln2^2
       - 127/72*x*z^-1
       + 25/6*x
       + 47/24*x*z
       - 4*x*z*rln2^2
       - 31/36*x*z^2
       - 1/6*pi^2*x^-2*[1+x]^-1
       + 1/3*pi^2*x^-2*[1+x]^-1*z
       + 1/6*pi^2*x^-2
       - 1/3*pi^2*x^-2*z
       + 2/3*pi^2*x^-1*[1+x]^-1*z^2
       - 1/6*pi^2*x^-1
       + 1/3*pi^2*x^-1*z
       - 2/3*pi^2*x^-1*z^2
       - 5/12*pi^2*[1-x]^-1
       + 5/6*pi^2*[1-x]^-1*z
       - pi^2*[1-x]^-1*z^2
       + 1/3*pi^2*[1+x]^-1
       - 2/3*pi^2*[1+x]^-1*z
       + pi^2*[1+x]^-1*z^2
       - 1/12*pi^2*z^-1
       + 5/12*pi^2
       - 1/6*pi^2*z
       + 2/3*pi^2*z^2
       - 1/12*pi^2*x*z^-1
       + 1/12*pi^2*x
       + 1/3*pi^2*x*z
       - 4*ln(1 + sqrtxz1 - z)*[1-x]^-1*rln2
       - 3*ln(1 + sqrtxz1 - z)*[1-x]^-1*sqrtxz1
       - 12*ln(1 + sqrtxz1 - z)*[1-x]^-1*z^2*rln2
       + 3*ln(1 + sqrtxz1 - z)*rln2
       + 2*ln(1 + sqrtxz1 - z)*sqrtxz1
       - 2*ln(1 + sqrtxz1 - z)*z*rln2
       + 12*ln(1 + sqrtxz1 - z)*z^2*rln2
       - ln(1 + sqrtxz1 - z)*x*rln2
       + 6*ln(1 + sqrtxz1 - z)*x*z*rln2
       - ln(1 + sqrtxz1 - z)^2
       + ln(1 + sqrtxz1 - z)^2*[1-x]^-1
       - 2*ln(1 + sqrtxz1 - z)^2*[1-x]^-1*z
       + 4*ln(1 + sqrtxz1 - z)^2*[1-x]^-1*z^2
       + 2*ln(1 + sqrtxz1 - z)^2*z
       - 4*ln(1 + sqrtxz1 - z)^2*z^2
       + ln(1 + sqrtxz1 - z)^2*x
       - 2*ln(1 + sqrtxz1 - z)^2*x*z
       - ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)
       + 2*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*[1-x]^-1
       + 4*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z
       + 4*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^2
       - 2*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*z
       - 4*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*z^2
       - ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*x
       - 2*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*x*z
       - 2*ln(1 + sqrtxz1 + z)*[1-x]^-1*rln2
       - 4*ln(1 + sqrtxz1 + z)*[1-x]^-1*z*rln2
       - 4*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^2*rln2
       + ln(1 + sqrtxz1 + z)*rln2
       + 2*ln(1 + sqrtxz1 + z)*z*rln2
       + 4*ln(1 + sqrtxz1 + z)*z^2*rln2
       + ln(1 + sqrtxz1 + z)*x*rln2
       + 2*ln(1 + sqrtxz1 + z)*x*z*rln2
       - 4*ln(z*sqrtxz3)*ArcTan(z*sqrtxz3)*z*sqrtxz3
       - 19/4*ln(x)
       - 2*ln(x)*x^-2*[1+x]^-1
       + 4*ln(x)*x^-2*[1+x]^-1*z
       + 2*ln(x)*x^-2
       - 4*ln(x)*x^-2*z
       + 8*ln(x)*x^-1*[1+x]^-1*z^2
       - 2*ln(x)*x^-1
       + 4*ln(x)*x^-1*z
       - 8*ln(x)*x^-1*z^2
       + 2/3*ln(x)*[1-x]^-1*z^-1
       + 1/4*ln(x)*[1-x]^-1
       + 5/2*ln(x)*[1-x]^-1*rln2
       + 3/2*ln(x)*[1-x]^-1*sqrtxz1
       + 1/4*ln(x)*[1-x]^-1*z
       - ln(x)*[1-x]^-1*z*rln2
       - 2/3*ln(x)*[1-x]^-1*z^2
       + 8*ln(x)*[1-x]^-1*z^2*rln2
       + 2*ln(x)*[1+x]^-1
       - 4*ln(x)*[1+x]^-1*z
       + 8*ln(x)*[1+x]^-1*z^2
       + 7/6*ln(x)*z^-1
       - 1/2*ln(x)*poly2^-1
       - 2*ln(x)*rln2
       - ln(x)*sqrtxz1
       - 1/2*ln(x)*z
       + 2*ln(x)*z*rln2
       + 4/3*ln(x)*z^2
       - 8*ln(x)*z^2*rln2
       - 11/6*ln(x)*x*z^-1
       - 1/2*ln(x)*x*[x-z]^-1
       + 1/2*ln(x)*x*poly2^-1
       + 3/4*ln(x)*x
       + ln(x)*x*rln2
       + 3/2*ln(x)*x*z
       - 4*ln(x)*x*z*rln2
       - 2/3*ln(x)*x*z^2
       + ln(x)*x^2*[x-z]^-1
       + 1/2*ln(x)*x^2*poly2^-1
       - 1/2*ln(x)*x^3*poly2^-1
       - 1/4*ln(x)*ln(1 - sqrtxz2 + x)*sqrtxz2^-1*poly2^-1
       - 3/4*ln(x)*ln(1 - sqrtxz2 + x)*sqrtxz2^-1
       - 3/2*ln(x)*ln(1 - sqrtxz2 + x)*x*sqrtxz2^-1
       + 3*ln(x)*ln(1 - sqrtxz2 + x)*x*z*sqrtxz2^-1
       + 1/2*ln(x)*ln(1 - sqrtxz2 + x)*x^2*sqrtxz2^-1*poly2^-1
       - 3/4*ln(x)*ln(1 - sqrtxz2 + x)*x^2*sqrtxz2^-1
       - 1/4*ln(x)*ln(1 - sqrtxz2 + x)*x^4*sqrtxz2^-1*poly2^-1
       + 1/4*ln(x)*ln(1 + sqrtxz2 + x)*sqrtxz2^-1*poly2^-1
       + 3/4*ln(x)*ln(1 + sqrtxz2 + x)*sqrtxz2^-1
       + 3/2*ln(x)*ln(1 + sqrtxz2 + x)*x*sqrtxz2^-1
       - 3*ln(x)*ln(1 + sqrtxz2 + x)*x*z*sqrtxz2^-1
       - 1/2*ln(x)*ln(1 + sqrtxz2 + x)*x^2*sqrtxz2^-1*poly2^-1
       + 3/4*ln(x)*ln(1 + sqrtxz2 + x)*x^2*sqrtxz2^-1
       + 1/4*ln(x)*ln(1 + sqrtxz2 + x)*x^4*sqrtxz2^-1*poly2^-1
       + ln(x)*ln(1 + sqrtxz1 - z)
       - ln(x)*ln(1 + sqrtxz1 - z)*[1-x]^-1
       + 2*ln(x)*ln(1 + sqrtxz1 - z)*[1-x]^-1*z
       - 4*ln(x)*ln(1 + sqrtxz1 - z)*[1-x]^-1*z^2
       - 2*ln(x)*ln(1 + sqrtxz1 - z)*z
       + 4*ln(x)*ln(1 + sqrtxz1 - z)*z^2
       - ln(x)*ln(1 + sqrtxz1 - z)*x
       + 2*ln(x)*ln(1 + sqrtxz1 - z)*x*z
       + ln(x)*ln(1 + sqrtxz1 + z)
       - 3/2*ln(x)*ln(1 + sqrtxz1 + z)*[1-x]^-1
       - ln(x)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z
       - 4*ln(x)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^2
       + 4*ln(x)*ln(1 + sqrtxz1 + z)*z^2
       + 2*ln(x)*ln(1 + sqrtxz1 + z)*x*z
       - 1/2*ln(x)*ln(1 + x*z^-1)
       + 1/2*ln(x)*ln(1 + x*z^-1)*x^-2*[1+x]^-1
       - ln(x)*ln(1 + x*z^-1)*x^-2*[1+x]^-1*z
       - 1/2*ln(x)*ln(1 + x*z^-1)*x^-2
       + ln(x)*ln(1 + x*z^-1)*x^-2*z
       - 2*ln(x)*ln(1 + x*z^-1)*x^-1*[1+x]^-1*z^2
       + 1/2*ln(x)*ln(1 + x*z^-1)*x^-1
       - ln(x)*ln(1 + x*z^-1)*x^-1*z
       + 2*ln(x)*ln(1 + x*z^-1)*x^-1*z^2
       + 1/2*ln(x)*ln(1 + x*z^-1)*[1+x]^-1
       - ln(x)*ln(1 + x*z^-1)*[1+x]^-1*z
       + ln(x)*ln(1 + x*z^-1)*z
       - 2*ln(x)*ln(1 + x*z^-1)*z^2
       + 1/2*ln(x)*ln(1 + x*z^-1)*x
       - ln(x)*ln(1 + x*z^-1)*x*z
       - 3*ln(x)*ln(1 + x)
       + 2*ln(x)*ln(1 + x)*[1+x]^-1
       - 4*ln(x)*ln(1 + x)*[1+x]^-1*z
       + 4*ln(x)*ln(1 + x)*[1+x]^-1*z^2
       + 6*ln(x)*ln(1 + x)*z
       - 4*ln(x)*ln(1 + x)*z^2
       - ln(x)*ln(1 + x)*x
       + 2*ln(x)*ln(1 + x)*x*z
       - ln(x)*ln(1 + x*z)
       + 1/2*ln(x)*ln(1 + x*z)*x^-2*[1+x]^-1
       - ln(x)*ln(1 + x*z)*x^-2*[1+x]^-1*z
       - 1/2*ln(x)*ln(1 + x*z)*x^-2
       + ln(x)*ln(1 + x*z)*x^-2*z
       - 2*ln(x)*ln(1 + x*z)*x^-1*[1+x]^-1*z^2
       + 1/2*ln(x)*ln(1 + x*z)*x^-1
       - ln(x)*ln(1 + x*z)*x^-1*z
       + 2*ln(x)*ln(1 + x*z)*x^-1*z^2
       + ln(x)*ln(1 + x*z)*[1-x]^-1
       + 2*ln(x)*ln(1 + x*z)*[1-x]^-1*z
       + 2*ln(x)*ln(1 + x*z)*[1-x]^-1*z^2
       + 1/2*ln(x)*ln(1 + x*z)*[1+x]^-1
       - ln(x)*ln(1 + x*z)*[1+x]^-1*z
       - 4*ln(x)*ln(1 + x*z)*z^2
       - 2*ln(x)*ln(1 + x*z)*x*z
       + 1/2*ln(x)*ln(z + x)
       - ln(x)*ln(z + x)*[1-x]^-1
       - 2*ln(x)*ln(z + x)*[1-x]^-1*z
       - 2*ln(x)*ln(z + x)*[1-x]^-1*z^2
       + ln(x)*ln(z + x)*z
       + 2*ln(x)*ln(z + x)*z^2
       + 1/2*ln(x)*ln(z + x)*x
       + ln(x)*ln(z + x)*x*z
       + 1/2*ln(x)^2
       + 3/4*ln(x)^2*x^-2*[1+x]^-1
       - 3/2*ln(x)^2*x^-2*[1+x]^-1*z
       - 3/4*ln(x)^2*x^-2
       + 3/2*ln(x)^2*x^-2*z
       - 3*ln(x)^2*x^-1*[1+x]^-1*z^2
       + 3/4*ln(x)^2*x^-1
       - 3/2*ln(x)^2*x^-1*z
       + 3*ln(x)^2*x^-1*z^2
       - ln(x)^2*[1-x]^-1
       + 2*ln(x)^2*[1-x]^-1*z
       - ln(x)^2*[1-x]^-1*z^2
       - 5/4*ln(x)^2*[1+x]^-1
       + 5/2*ln(x)^2*[1+x]^-1*z
       - 4*ln(x)^2*[1+x]^-1*z^2
       + 1/2*ln(x)^2*z^-1
       - 3*ln(x)^2*z
       + 2*ln(x)^2*z^2
       + 1/2*ln(x)^2*x*z^-1
       - 2*ln(x)^2*x*z
       - ln(x)*ln([1-x])
       + ln(x)*ln([1-x])*x^-2*[1+x]^-1
       - 2*ln(x)*ln([1-x])*x^-2*[1+x]^-1*z
       - ln(x)*ln([1-x])*x^-2
       + 2*ln(x)*ln([1-x])*x^-2*z
       - 4*ln(x)*ln([1-x])*x^-1*[1+x]^-1*z^2
       + ln(x)*ln([1-x])*x^-1
       - 2*ln(x)*ln([1-x])*x^-1*z
       + 4*ln(x)*ln([1-x])*x^-1*z^2
       + 3/2*ln(x)*ln([1-x])*[1-x]^-1
       - 3*ln(x)*ln([1-x])*[1-x]^-1*z
       + 4*ln(x)*ln([1-x])*[1-x]^-1*z^2
       - ln(x)*ln([1-x])*[1+x]^-1
       + 2*ln(x)*ln([1-x])*[1+x]^-1*z
       - 4*ln(x)*ln([1-x])*[1+x]^-1*z^2
       + 2*ln(x)*ln([1-x])*z
       - 4*ln(x)*ln([1-x])*z^2
       + ln(x)*ln([1+x])
       - ln(x)*ln([1+x])*x^-2*[1+x]^-1
       + 2*ln(x)*ln([1+x])*x^-2*[1+x]^-1*z
       + ln(x)*ln([1+x])*x^-2
       - 2*ln(x)*ln([1+x])*x^-2*z
       + 4*ln(x)*ln([1+x])*x^-1*[1+x]^-1*z^2
       - ln(x)*ln([1+x])*x^-1
       + 2*ln(x)*ln([1+x])*x^-1*z
       - 4*ln(x)*ln([1+x])*x^-1*z^2
       - ln(x)*ln([1+x])*[1+x]^-1
       + 2*ln(x)*ln([1+x])*[1+x]^-1*z
       - 2*ln(x)*ln([1+x])*z
       + 4*ln(x)*ln([1+x])*z^2
       - ln(x)*ln([1+x])*x
       + 2*ln(x)*ln([1+x])*x*z
       - ln(x)*ln(z)
       - 1/2*ln(x)*ln(z)*x^-2*[1+x]^-1
       + ln(x)*ln(z)*x^-2*[1+x]^-1*z
       + 1/2*ln(x)*ln(z)*x^-2
       - ln(x)*ln(z)*x^-2*z
       + 2*ln(x)*ln(z)*x^-1*[1+x]^-1*z^2
       - 1/2*ln(x)*ln(z)*x^-1
       + ln(x)*ln(z)*x^-1*z
       - 2*ln(x)*ln(z)*x^-1*z^2
       + 7/2*ln(x)*ln(z)*[1-x]^-1
       - 1/2*ln(x)*ln(z)*[1-x]^-1*z
       + 4*ln(x)*ln(z)*[1-x]^-1*z^2
       - 1/2*ln(x)*ln(z)*[1+x]^-1
       + ln(x)*ln(z)*[1+x]^-1*z
       - 1/2*ln(x)*ln(z)*z^-1
       - ln(x)*ln(z)*z
       - 2*ln(x)*ln(z)*z^2
       - 1/2*ln(x)*ln(z)*x*z^-1
       - ln(x)*ln(z)*x*z
       + ln(x)*ln([1-z])
       - 1/2*ln(x)*ln([1-z])*z^-1
       - 1/2*ln(x)*ln([1-z])*x*z^-1
       + ln(x)*ln([1-z])*x
       + 13/4*ln([1-x])
       - 13/12*ln([1-x])*z^-1
       - 1/4*ln([1-x])*z
       - 2/3*ln([1-x])*z^2
       + 17/12*ln([1-x])*x*z^-1
       - 9/4*ln([1-x])*x
       - 3/4*ln([1-x])*x*z
       + 1/3*ln([1-x])*x*z^2
       + 1/2*ln([1-x])*ln(z)
       + ln([1-x])*ln(z)*z
       + 1/2*ln([1-x])*ln(z)*x
       + 9/2*ln(z)
       - 3/2*ln(z)*[1-x]^-1
       + 7/2*ln(z)*[1-x]^-1*rln2
       + 3/2*ln(z)*[1-x]^-1*sqrtxz1
       - 3/2*ln(z)*[1-x]^-1*z
       + 5*ln(z)*[1-x]^-1*z*rln2
       + 8*ln(z)*[1-x]^-1*z^2*rln2
       - 13/12*ln(z)*z^-1
       - 1/2*ln(z)*[1-z]^-1
       - 1/2*ln(z)*poly2^-1
       - 2*ln(z)*rln2
       - ln(z)*sqrtxz1
       + 1/4*ln(z)*z
       - 2*ln(z)*z*rln2
       - 2/3*ln(z)*z^2
       - 8*ln(z)*z^2*rln2
       + 17/12*ln(z)*x*z^-1
       + 1/2*ln(z)*x*[1-z]^-1
       + 1/2*ln(z)*x*[x-z]^-1
       - 1/2*ln(z)*x*poly2^-1
       - 1/2*ln(z)*x
       - ln(z)*x*rln2
       - 3/4*ln(z)*x*z
       - 4*ln(z)*x*z*rln2
       + 1/3*ln(z)*x*z^2
       - ln(z)*x^2*[x-z]^-1
       + 1/2*ln(z)*x^2*poly2^-1
       + 1/2*ln(z)*x^3*poly2^-1
       + ln(z)*ln(1 + sqrtxz1 - z)
       - 3/2*ln(z)*ln(1 + sqrtxz1 - z)*[1-x]^-1
       - ln(z)*ln(1 + sqrtxz1 - z)*[1-x]^-1*z
       - 4*ln(z)*ln(1 + sqrtxz1 - z)*[1-x]^-1*z^2
       + 4*ln(z)*ln(1 + sqrtxz1 - z)*z^2
       + 2*ln(z)*ln(1 + sqrtxz1 - z)*x*z
       + ln(z)*ln(1 + sqrtxz1 + z)
       - 2*ln(z)*ln(1 + sqrtxz1 + z)*[1-x]^-1
       - 4*ln(z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z
       - 4*ln(z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^2
       + 2*ln(z)*ln(1 + sqrtxz1 + z)*z
       + 4*ln(z)*ln(1 + sqrtxz1 + z)*z^2
       + ln(z)*ln(1 + sqrtxz1 + z)*x
       + 2*ln(z)*ln(1 + sqrtxz1 + z)*x*z
       + 1/2*ln(z)*ln(1 + x*z^-1)
       - 1/2*ln(z)*ln(1 + x*z^-1)*x^-2*[1+x]^-1
       + ln(z)*ln(1 + x*z^-1)*x^-2*[1+x]^-1*z
       + 1/2*ln(z)*ln(1 + x*z^-1)*x^-2
       - ln(z)*ln(1 + x*z^-1)*x^-2*z
       + 2*ln(z)*ln(1 + x*z^-1)*x^-1*[1+x]^-1*z^2
       - 1/2*ln(z)*ln(1 + x*z^-1)*x^-1
       + ln(z)*ln(1 + x*z^-1)*x^-1*z
       - 2*ln(z)*ln(1 + x*z^-1)*x^-1*z^2
       - 1/2*ln(z)*ln(1 + x*z^-1)*[1+x]^-1
       + ln(z)*ln(1 + x*z^-1)*[1+x]^-1*z
       - ln(z)*ln(1 + x*z^-1)*z
       + 2*ln(z)*ln(1 + x*z^-1)*z^2
       - 1/2*ln(z)*ln(1 + x*z^-1)*x
       + ln(z)*ln(1 + x*z^-1)*x*z
       - ln(z)*ln(1 + x*z)
       + 1/2*ln(z)*ln(1 + x*z)*x^-2*[1+x]^-1
       - ln(z)*ln(1 + x*z)*x^-2*[1+x]^-1*z
       - 1/2*ln(z)*ln(1 + x*z)*x^-2
       + ln(z)*ln(1 + x*z)*x^-2*z
       - 2*ln(z)*ln(1 + x*z)*x^-1*[1+x]^-1*z^2
       + 1/2*ln(z)*ln(1 + x*z)*x^-1
       - ln(z)*ln(1 + x*z)*x^-1*z
       + 2*ln(z)*ln(1 + x*z)*x^-1*z^2
       + ln(z)*ln(1 + x*z)*[1-x]^-1
       + 2*ln(z)*ln(1 + x*z)*[1-x]^-1*z
       + 2*ln(z)*ln(1 + x*z)*[1-x]^-1*z^2
       + 1/2*ln(z)*ln(1 + x*z)*[1+x]^-1
       - ln(z)*ln(1 + x*z)*[1+x]^-1*z
       - 4*ln(z)*ln(1 + x*z)*z^2
       - 2*ln(z)*ln(1 + x*z)*x*z
       - 1/2*ln(z)*ln(z + x)
       + ln(z)*ln(z + x)*[1-x]^-1
       + 2*ln(z)*ln(z + x)*[1-x]^-1*z
       + 2*ln(z)*ln(z + x)*[1-x]^-1*z^2
       - ln(z)*ln(z + x)*z
       - 2*ln(z)*ln(z + x)*z^2
       - 1/2*ln(z)*ln(z + x)*x
       - ln(z)*ln(z + x)*x*z
       + 3/2*ln(z)^2
       - 1/4*ln(z)^2*x^-2*[1+x]^-1
       + 1/2*ln(z)^2*x^-2*[1+x]^-1*z
       + 1/4*ln(z)^2*x^-2
       - 1/2*ln(z)^2*x^-2*z
       + ln(z)^2*x^-1*[1+x]^-1*z^2
       - 1/4*ln(z)^2*x^-1
       + 1/2*ln(z)^2*x^-1*z
       - ln(z)^2*x^-1*z^2
       - 1/2*ln(z)^2*[1-x]^-1
       + ln(z)^2*[1-x]^-1*z
       - ln(z)^2*[1-x]^-1*z^2
       - 1/4*ln(z)^2*[1+x]^-1
       + 1/2*ln(z)^2*[1+x]^-1*z
       - ln(z)^2*z
       + 2*ln(z)^2*z^2
       + ln(z)^2*x*z
       - ln(z)*ln([1-z])
       + 1/2*ln(z)*ln([1-z])*[1-x]^-1
       - ln(z)*ln([1-z])*[1-x]^-1*z
       + 2*ln(z)*ln([1-z])*z
       + 13/4*ln([1-z])
       - 13/12*ln([1-z])*z^-1
       - 1/4*ln([1-z])*z
       - 2/3*ln([1-z])*z^2
       + 17/12*ln([1-z])*x*z^-1
       - 9/4*ln([1-z])*x
       - 3/4*ln([1-z])*x*z
       + 1/3*ln([1-z])*x*z^2
       + 4*ln(sqrtxz3)*ArcTan(sqrtxz3)*z*sqrtxz3
       - 1/4*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*sqrtxz2^-1*poly2^-1
       - 3/4*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*sqrtxz2^-1
       - 3/2*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x*sqrtxz2^-1
       + 3*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x*z*sqrtxz2^-1
       + 1/2*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x^2*sqrtxz2^-1*poly2^-1
       - 3/4*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x^2*sqrtxz2^-1
       - 1/4*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x^4*sqrtxz2^-1*poly2^-1
       + 1/4*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*sqrtxz2^-1*poly2^-1
       + 3/4*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*sqrtxz2^-1
       + 3/2*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x*sqrtxz2^-1
       - 3*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x*z*sqrtxz2^-1
       - 1/2*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x^2*sqrtxz2^-1*poly2^-1
       + 3/4*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x^2*sqrtxz2^-1
       + 1/4*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x^4*sqrtxz2^-1*poly2^-1
       + 1/2*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1
       + 3*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*z
       - 2*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z
       - Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x
       - 1/2*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1
       - 3*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*z
       + 2*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z
       + Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x
       + 1/4*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*sqrtxz2^-1*poly2^-1
       + 3/4*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*sqrtxz2^-1
       + 3/2*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x*sqrtxz2^-1
       - 3*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x*z*sqrtxz2^-1
       - 1/2*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x^2*sqrtxz2^-1*poly2^-1
       + 3/4*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x^2*sqrtxz2^-1
       + 1/4*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x^4*sqrtxz2^-1*poly2^-1
       - 1/4*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*sqrtxz2^-1*poly2^-1
       - 3/4*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*sqrtxz2^-1
       - 3/2*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x*sqrtxz2^-1
       + 3*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x*z*sqrtxz2^-1
       + 1/2*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x^2*sqrtxz2^-1*poly2^-1
       - 3/4*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x^2*sqrtxz2^-1
       - 1/4*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x^4*sqrtxz2^-1*poly2^-1
       + Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)
       - 3/2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*[1-x]^-1
       - Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*[1-x]^-1*z
       - 4*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*[1-x]^-1*z^2
       + 4*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*z^2
       + 2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*x*z
       - Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)
       + 3/2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*[1-x]^-1
       + Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*[1-x]^-1*z
       + 4*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*[1-x]^-1*z^2
       - 4*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*z^2
       - 2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*x*z
       + Li2(1 - x*z^-1)
       - 1/2*Li2(1 - x*z^-1)*[1-x]^-1
       + Li2(1 - x*z^-1)*[1-x]^-1*z
       - 2*Li2(1 - x*z^-1)*z
       + 1/2*Li2( - x*z^-1)*x^-2*[1+x]^-1
       - Li2( - x*z^-1)*x^-2*[1+x]^-1*z
       - 1/2*Li2( - x*z^-1)*x^-2
       + Li2( - x*z^-1)*x^-2*z
       - 2*Li2( - x*z^-1)*x^-1*[1+x]^-1*z^2
       + 1/2*Li2( - x*z^-1)*x^-1
       - Li2( - x*z^-1)*x^-1*z
       + 2*Li2( - x*z^-1)*x^-1*z^2
       - Li2( - x*z^-1)*[1-x]^-1
       - 2*Li2( - x*z^-1)*[1-x]^-1*z
       - 2*Li2( - x*z^-1)*[1-x]^-1*z^2
       + 1/2*Li2( - x*z^-1)*[1+x]^-1
       - Li2( - x*z^-1)*[1+x]^-1*z
       + 2*Li2( - x*z^-1)*z
       + Li2( - x*z^-1)*x
       - 2*Li2( - x)
       - Li2( - x)*x^-2*[1+x]^-1
       + 2*Li2( - x)*x^-2*[1+x]^-1*z
       + Li2( - x)*x^-2
       - 2*Li2( - x)*x^-2*z
       + 4*Li2( - x)*x^-1*[1+x]^-1*z^2
       - Li2( - x)*x^-1
       + 2*Li2( - x)*x^-1*z
       - 4*Li2( - x)*x^-1*z^2
       + Li2( - x)*[1+x]^-1
       - 2*Li2( - x)*[1+x]^-1*z
       + 4*Li2( - x)*[1+x]^-1*z^2
       + 4*Li2( - x)*z
       - 2*Li2( - x)*x
       + 4*Li2( - x)*x*z
       - Li2( - x*z)
       + 1/2*Li2( - x*z)*x^-2*[1+x]^-1
       - Li2( - x*z)*x^-2*[1+x]^-1*z
       - 1/2*Li2( - x*z)*x^-2
       + Li2( - x*z)*x^-2*z
       - 2*Li2( - x*z)*x^-1*[1+x]^-1*z^2
       + 1/2*Li2( - x*z)*x^-1
       - Li2( - x*z)*x^-1*z
       + 2*Li2( - x*z)*x^-1*z^2
       + Li2( - x*z)*[1-x]^-1
       + 2*Li2( - x*z)*[1-x]^-1*z
       + 2*Li2( - x*z)*[1-x]^-1*z^2
       + 1/2*Li2( - x*z)*[1+x]^-1
       - Li2( - x*z)*[1+x]^-1*z
       - 4*Li2( - x*z)*z^2
       - 2*Li2( - x*z)*x*z
       - 2*Li2(x)
       + Li2(x)*x^-2*[1+x]^-1
       - 2*Li2(x)*x^-2*[1+x]^-1*z
       - Li2(x)*x^-2
       + 2*Li2(x)*x^-2*z
       - 4*Li2(x)*x^-1*[1+x]^-1*z^2
       + Li2(x)*x^-1
       - 2*Li2(x)*x^-1*z
       + 4*Li2(x)*x^-1*z^2
       + 3/2*Li2(x)*[1-x]^-1
       - 3*Li2(x)*[1-x]^-1*z
       + 4*Li2(x)*[1-x]^-1*z^2
       - Li2(x)*[1+x]^-1
       + 2*Li2(x)*[1+x]^-1*z
       - 4*Li2(x)*[1+x]^-1*z^2
       + 1/2*Li2(x)*z^-1
       + 2*Li2(x)*z
       - 4*Li2(x)*z^2
       + 1/2*Li2(x)*x*z^-1
       - Li2(x)*x
       - 3/2*Li2(z)
       + 1/2*Li2(z)*[1-x]^-1
       - Li2(z)*[1-x]^-1*z
       + Li2(z)*z
       - 1/2*Li2(z)*x
       + 2*InvTanInt( - sqrtxz3)*z*sqrtxz3
       + 4*InvTanInt(z*sqrtxz3)*z*sqrtxz3
       - 2*InvTanInt(sqrtxz3)*z*sqrtxz3
       )

       + NC * ( 17/6
       - 3*[1-x]^-1*rln2^2
       - 3*[1-x]^-1*sqrtxz1*rln2
       - 2*[1-x]^-1*z*rln2^2
       - 8*[1-x]^-1*z^2*rln2^2
       - 113/72*z^-1
       + 2*rln2^2
       + 2*sqrtxz1*rln2
       + 85/24*z
       - 65/36*z^2
       + 8*z^2*rln2^2
       + 127/72*x*z^-1
       - 25/6*x
       - 47/24*x*z
       + 4*x*z*rln2^2
       + 31/36*x*z^2
       + 1/6*pi^2*x^-2*[1+x]^-1
       - 1/3*pi^2*x^-2*[1+x]^-1*z
       - 1/6*pi^2*x^-2
       + 1/3*pi^2*x^-2*z
       - 2/3*pi^2*x^-1*[1+x]^-1*z^2
       + 1/6*pi^2*x^-1
       - 1/3*pi^2*x^-1*z
       + 2/3*pi^2*x^-1*z^2
       + 5/12*pi^2*[1-x]^-1
       - 5/6*pi^2*[1-x]^-1*z
       + pi^2*[1-x]^-1*z^2
       - 1/3*pi^2*[1+x]^-1
       + 2/3*pi^2*[1+x]^-1*z
       - pi^2*[1+x]^-1*z^2
       + 1/12*pi^2*z^-1
       - 5/12*pi^2
       + 1/6*pi^2*z
       - 2/3*pi^2*z^2
       + 1/12*pi^2*x*z^-1
       - 1/12*pi^2*x
       - 1/3*pi^2*x*z
       + 4*ln(1 + sqrtxz1 - z)*[1-x]^-1*rln2
       + 3*ln(1 + sqrtxz1 - z)*[1-x]^-1*sqrtxz1
       + 12*ln(1 + sqrtxz1 - z)*[1-x]^-1*z^2*rln2
       - 3*ln(1 + sqrtxz1 - z)*rln2
       - 2*ln(1 + sqrtxz1 - z)*sqrtxz1
       + 2*ln(1 + sqrtxz1 - z)*z*rln2
       - 12*ln(1 + sqrtxz1 - z)*z^2*rln2
       + ln(1 + sqrtxz1 - z)*x*rln2
       - 6*ln(1 + sqrtxz1 - z)*x*z*rln2
       + ln(1 + sqrtxz1 - z)^2
       - ln(1 + sqrtxz1 - z)^2*[1-x]^-1
       + 2*ln(1 + sqrtxz1 - z)^2*[1-x]^-1*z
       - 4*ln(1 + sqrtxz1 - z)^2*[1-x]^-1*z^2
       - 2*ln(1 + sqrtxz1 - z)^2*z
       + 4*ln(1 + sqrtxz1 - z)^2*z^2
       - ln(1 + sqrtxz1 - z)^2*x
       + 2*ln(1 + sqrtxz1 - z)^2*x*z
       + ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)
       - 2*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*[1-x]^-1
       - 4*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z
       - 4*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^2
       + 2*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*z
       + 4*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*z^2
       + ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*x
       + 2*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*x*z
       + 2*ln(1 + sqrtxz1 + z)*[1-x]^-1*rln2
       + 4*ln(1 + sqrtxz1 + z)*[1-x]^-1*z*rln2
       + 4*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^2*rln2
       - ln(1 + sqrtxz1 + z)*rln2
       - 2*ln(1 + sqrtxz1 + z)*z*rln2
       - 4*ln(1 + sqrtxz1 + z)*z^2*rln2
       - ln(1 + sqrtxz1 + z)*x*rln2
       - 2*ln(1 + sqrtxz1 + z)*x*z*rln2
       + 4*ln(z*sqrtxz3)*ArcTan(z*sqrtxz3)*z*sqrtxz3
       + 19/4*ln(x)
       + 2*ln(x)*x^-2*[1+x]^-1
       - 4*ln(x)*x^-2*[1+x]^-1*z
       - 2*ln(x)*x^-2
       + 4*ln(x)*x^-2*z
       - 8*ln(x)*x^-1*[1+x]^-1*z^2
       + 2*ln(x)*x^-1
       - 4*ln(x)*x^-1*z
       + 8*ln(x)*x^-1*z^2
       - 2/3*ln(x)*[1-x]^-1*z^-1
       - 1/4*ln(x)*[1-x]^-1
       - 5/2*ln(x)*[1-x]^-1*rln2
       - 3/2*ln(x)*[1-x]^-1*sqrtxz1
       - 1/4*ln(x)*[1-x]^-1*z
       + ln(x)*[1-x]^-1*z*rln2
       + 2/3*ln(x)*[1-x]^-1*z^2
       - 8*ln(x)*[1-x]^-1*z^2*rln2
       - 2*ln(x)*[1+x]^-1
       + 4*ln(x)*[1+x]^-1*z
       - 8*ln(x)*[1+x]^-1*z^2
       - 7/6*ln(x)*z^-1
       + 1/2*ln(x)*poly2^-1
       + 2*ln(x)*rln2
       + ln(x)*sqrtxz1
       + 1/2*ln(x)*z
       - 2*ln(x)*z*rln2
       - 4/3*ln(x)*z^2
       + 8*ln(x)*z^2*rln2
       + 11/6*ln(x)*x*z^-1
       + 1/2*ln(x)*x*[x-z]^-1
       - 1/2*ln(x)*x*poly2^-1
       - 3/4*ln(x)*x
       - ln(x)*x*rln2
       - 3/2*ln(x)*x*z
       + 4*ln(x)*x*z*rln2
       + 2/3*ln(x)*x*z^2
       - ln(x)*x^2*[x-z]^-1
       - 1/2*ln(x)*x^2*poly2^-1
       + 1/2*ln(x)*x^3*poly2^-1
       + 1/4*ln(x)*ln(1 - sqrtxz2 + x)*sqrtxz2^-1*poly2^-1
       + 3/4*ln(x)*ln(1 - sqrtxz2 + x)*sqrtxz2^-1
       + 3/2*ln(x)*ln(1 - sqrtxz2 + x)*x*sqrtxz2^-1
       - 3*ln(x)*ln(1 - sqrtxz2 + x)*x*z*sqrtxz2^-1
       - 1/2*ln(x)*ln(1 - sqrtxz2 + x)*x^2*sqrtxz2^-1*poly2^-1
       + 3/4*ln(x)*ln(1 - sqrtxz2 + x)*x^2*sqrtxz2^-1
       + 1/4*ln(x)*ln(1 - sqrtxz2 + x)*x^4*sqrtxz2^-1*poly2^-1
       - 1/4*ln(x)*ln(1 + sqrtxz2 + x)*sqrtxz2^-1*poly2^-1
       - 3/4*ln(x)*ln(1 + sqrtxz2 + x)*sqrtxz2^-1
       - 3/2*ln(x)*ln(1 + sqrtxz2 + x)*x*sqrtxz2^-1
       + 3*ln(x)*ln(1 + sqrtxz2 + x)*x*z*sqrtxz2^-1
       + 1/2*ln(x)*ln(1 + sqrtxz2 + x)*x^2*sqrtxz2^-1*poly2^-1
       - 3/4*ln(x)*ln(1 + sqrtxz2 + x)*x^2*sqrtxz2^-1
       - 1/4*ln(x)*ln(1 + sqrtxz2 + x)*x^4*sqrtxz2^-1*poly2^-1
       - ln(x)*ln(1 + sqrtxz1 - z)
       + ln(x)*ln(1 + sqrtxz1 - z)*[1-x]^-1
       - 2*ln(x)*ln(1 + sqrtxz1 - z)*[1-x]^-1*z
       + 4*ln(x)*ln(1 + sqrtxz1 - z)*[1-x]^-1*z^2
       + 2*ln(x)*ln(1 + sqrtxz1 - z)*z
       - 4*ln(x)*ln(1 + sqrtxz1 - z)*z^2
       + ln(x)*ln(1 + sqrtxz1 - z)*x
       - 2*ln(x)*ln(1 + sqrtxz1 - z)*x*z
       - ln(x)*ln(1 + sqrtxz1 + z)
       + 3/2*ln(x)*ln(1 + sqrtxz1 + z)*[1-x]^-1
       + ln(x)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z
       + 4*ln(x)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^2
       - 4*ln(x)*ln(1 + sqrtxz1 + z)*z^2
       - 2*ln(x)*ln(1 + sqrtxz1 + z)*x*z
       + 1/2*ln(x)*ln(1 + x*z^-1)
       - 1/2*ln(x)*ln(1 + x*z^-1)*x^-2*[1+x]^-1
       + ln(x)*ln(1 + x*z^-1)*x^-2*[1+x]^-1*z
       + 1/2*ln(x)*ln(1 + x*z^-1)*x^-2
       - ln(x)*ln(1 + x*z^-1)*x^-2*z
       + 2*ln(x)*ln(1 + x*z^-1)*x^-1*[1+x]^-1*z^2
       - 1/2*ln(x)*ln(1 + x*z^-1)*x^-1
       + ln(x)*ln(1 + x*z^-1)*x^-1*z
       - 2*ln(x)*ln(1 + x*z^-1)*x^-1*z^2
       - 1/2*ln(x)*ln(1 + x*z^-1)*[1+x]^-1
       + ln(x)*ln(1 + x*z^-1)*[1+x]^-1*z
       - ln(x)*ln(1 + x*z^-1)*z
       + 2*ln(x)*ln(1 + x*z^-1)*z^2
       - 1/2*ln(x)*ln(1 + x*z^-1)*x
       + ln(x)*ln(1 + x*z^-1)*x*z
       + 3*ln(x)*ln(1 + x)
       - 2*ln(x)*ln(1 + x)*[1+x]^-1
       + 4*ln(x)*ln(1 + x)*[1+x]^-1*z
       - 4*ln(x)*ln(1 + x)*[1+x]^-1*z^2
       - 6*ln(x)*ln(1 + x)*z
       + 4*ln(x)*ln(1 + x)*z^2
       + ln(x)*ln(1 + x)*x
       - 2*ln(x)*ln(1 + x)*x*z
       + ln(x)*ln(1 + x*z)
       - 1/2*ln(x)*ln(1 + x*z)*x^-2*[1+x]^-1
       + ln(x)*ln(1 + x*z)*x^-2*[1+x]^-1*z
       + 1/2*ln(x)*ln(1 + x*z)*x^-2
       - ln(x)*ln(1 + x*z)*x^-2*z
       + 2*ln(x)*ln(1 + x*z)*x^-1*[1+x]^-1*z^2
       - 1/2*ln(x)*ln(1 + x*z)*x^-1
       + ln(x)*ln(1 + x*z)*x^-1*z
       - 2*ln(x)*ln(1 + x*z)*x^-1*z^2
       - ln(x)*ln(1 + x*z)*[1-x]^-1
       - 2*ln(x)*ln(1 + x*z)*[1-x]^-1*z
       - 2*ln(x)*ln(1 + x*z)*[1-x]^-1*z^2
       - 1/2*ln(x)*ln(1 + x*z)*[1+x]^-1
       + ln(x)*ln(1 + x*z)*[1+x]^-1*z
       + 4*ln(x)*ln(1 + x*z)*z^2
       + 2*ln(x)*ln(1 + x*z)*x*z
       - 1/2*ln(x)*ln(z + x)
       + ln(x)*ln(z + x)*[1-x]^-1
       + 2*ln(x)*ln(z + x)*[1-x]^-1*z
       + 2*ln(x)*ln(z + x)*[1-x]^-1*z^2
       - ln(x)*ln(z + x)*z
       - 2*ln(x)*ln(z + x)*z^2
       - 1/2*ln(x)*ln(z + x)*x
       - ln(x)*ln(z + x)*x*z
       - 1/2*ln(x)^2
       - 3/4*ln(x)^2*x^-2*[1+x]^-1
       + 3/2*ln(x)^2*x^-2*[1+x]^-1*z
       + 3/4*ln(x)^2*x^-2
       - 3/2*ln(x)^2*x^-2*z
       + 3*ln(x)^2*x^-1*[1+x]^-1*z^2
       - 3/4*ln(x)^2*x^-1
       + 3/2*ln(x)^2*x^-1*z
       - 3*ln(x)^2*x^-1*z^2
       + ln(x)^2*[1-x]^-1
       - 2*ln(x)^2*[1-x]^-1*z
       + ln(x)^2*[1-x]^-1*z^2
       + 5/4*ln(x)^2*[1+x]^-1
       - 5/2*ln(x)^2*[1+x]^-1*z
       + 4*ln(x)^2*[1+x]^-1*z^2
       - 1/2*ln(x)^2*z^-1
       + 3*ln(x)^2*z
       - 2*ln(x)^2*z^2
       - 1/2*ln(x)^2*x*z^-1
       + 2*ln(x)^2*x*z
       + ln(x)*ln([1-x])
       - ln(x)*ln([1-x])*x^-2*[1+x]^-1
       + 2*ln(x)*ln([1-x])*x^-2*[1+x]^-1*z
       + ln(x)*ln([1-x])*x^-2
       - 2*ln(x)*ln([1-x])*x^-2*z
       + 4*ln(x)*ln([1-x])*x^-1*[1+x]^-1*z^2
       - ln(x)*ln([1-x])*x^-1
       + 2*ln(x)*ln([1-x])*x^-1*z
       - 4*ln(x)*ln([1-x])*x^-1*z^2
       - 3/2*ln(x)*ln([1-x])*[1-x]^-1
       + 3*ln(x)*ln([1-x])*[1-x]^-1*z
       - 4*ln(x)*ln([1-x])*[1-x]^-1*z^2
       + ln(x)*ln([1-x])*[1+x]^-1
       - 2*ln(x)*ln([1-x])*[1+x]^-1*z
       + 4*ln(x)*ln([1-x])*[1+x]^-1*z^2
       - 2*ln(x)*ln([1-x])*z
       + 4*ln(x)*ln([1-x])*z^2
       - ln(x)*ln([1+x])
       + ln(x)*ln([1+x])*x^-2*[1+x]^-1
       - 2*ln(x)*ln([1+x])*x^-2*[1+x]^-1*z
       - ln(x)*ln([1+x])*x^-2
       + 2*ln(x)*ln([1+x])*x^-2*z
       - 4*ln(x)*ln([1+x])*x^-1*[1+x]^-1*z^2
       + ln(x)*ln([1+x])*x^-1
       - 2*ln(x)*ln([1+x])*x^-1*z
       + 4*ln(x)*ln([1+x])*x^-1*z^2
       + ln(x)*ln([1+x])*[1+x]^-1
       - 2*ln(x)*ln([1+x])*[1+x]^-1*z
       + 2*ln(x)*ln([1+x])*z
       - 4*ln(x)*ln([1+x])*z^2
       + ln(x)*ln([1+x])*x
       - 2*ln(x)*ln([1+x])*x*z
       + ln(x)*ln(z)
       + 1/2*ln(x)*ln(z)*x^-2*[1+x]^-1
       - ln(x)*ln(z)*x^-2*[1+x]^-1*z
       - 1/2*ln(x)*ln(z)*x^-2
       + ln(x)*ln(z)*x^-2*z
       - 2*ln(x)*ln(z)*x^-1*[1+x]^-1*z^2
       + 1/2*ln(x)*ln(z)*x^-1
       - ln(x)*ln(z)*x^-1*z
       + 2*ln(x)*ln(z)*x^-1*z^2
       - 7/2*ln(x)*ln(z)*[1-x]^-1
       + 1/2*ln(x)*ln(z)*[1-x]^-1*z
       - 4*ln(x)*ln(z)*[1-x]^-1*z^2
       + 1/2*ln(x)*ln(z)*[1+x]^-1
       - ln(x)*ln(z)*[1+x]^-1*z
       + 1/2*ln(x)*ln(z)*z^-1
       + ln(x)*ln(z)*z
       + 2*ln(x)*ln(z)*z^2
       + 1/2*ln(x)*ln(z)*x*z^-1
       + ln(x)*ln(z)*x*z
       - ln(x)*ln([1-z])
       + 1/2*ln(x)*ln([1-z])*z^-1
       + 1/2*ln(x)*ln([1-z])*x*z^-1
       - ln(x)*ln([1-z])*x
       - 13/4*ln([1-x])
       + 13/12*ln([1-x])*z^-1
       + 1/4*ln([1-x])*z
       + 2/3*ln([1-x])*z^2
       - 17/12*ln([1-x])*x*z^-1
       + 9/4*ln([1-x])*x
       + 3/4*ln([1-x])*x*z
       - 1/3*ln([1-x])*x*z^2
       - 1/2*ln([1-x])*ln(z)
       - ln([1-x])*ln(z)*z
       - 1/2*ln([1-x])*ln(z)*x
       - 9/2*ln(z)
       + 3/2*ln(z)*[1-x]^-1
       - 7/2*ln(z)*[1-x]^-1*rln2
       - 3/2*ln(z)*[1-x]^-1*sqrtxz1
       + 3/2*ln(z)*[1-x]^-1*z
       - 5*ln(z)*[1-x]^-1*z*rln2
       - 8*ln(z)*[1-x]^-1*z^2*rln2
       + 13/12*ln(z)*z^-1
       + 1/2*ln(z)*[1-z]^-1
       + 1/2*ln(z)*poly2^-1
       + 2*ln(z)*rln2
       + ln(z)*sqrtxz1
       - 1/4*ln(z)*z
       + 2*ln(z)*z*rln2
       + 2/3*ln(z)*z^2
       + 8*ln(z)*z^2*rln2
       - 17/12*ln(z)*x*z^-1
       - 1/2*ln(z)*x*[1-z]^-1
       - 1/2*ln(z)*x*[x-z]^-1
       + 1/2*ln(z)*x*poly2^-1
       + 1/2*ln(z)*x
       + ln(z)*x*rln2
       + 3/4*ln(z)*x*z
       + 4*ln(z)*x*z*rln2
       - 1/3*ln(z)*x*z^2
       + ln(z)*x^2*[x-z]^-1
       - 1/2*ln(z)*x^2*poly2^-1
       - 1/2*ln(z)*x^3*poly2^-1
       - ln(z)*ln(1 + sqrtxz1 - z)
       + 3/2*ln(z)*ln(1 + sqrtxz1 - z)*[1-x]^-1
       + ln(z)*ln(1 + sqrtxz1 - z)*[1-x]^-1*z
       + 4*ln(z)*ln(1 + sqrtxz1 - z)*[1-x]^-1*z^2
       - 4*ln(z)*ln(1 + sqrtxz1 - z)*z^2
       - 2*ln(z)*ln(1 + sqrtxz1 - z)*x*z
       - ln(z)*ln(1 + sqrtxz1 + z)
       + 2*ln(z)*ln(1 + sqrtxz1 + z)*[1-x]^-1
       + 4*ln(z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z
       + 4*ln(z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^2
       - 2*ln(z)*ln(1 + sqrtxz1 + z)*z
       - 4*ln(z)*ln(1 + sqrtxz1 + z)*z^2
       - ln(z)*ln(1 + sqrtxz1 + z)*x
       - 2*ln(z)*ln(1 + sqrtxz1 + z)*x*z
       - 1/2*ln(z)*ln(1 + x*z^-1)
       + 1/2*ln(z)*ln(1 + x*z^-1)*x^-2*[1+x]^-1
       - ln(z)*ln(1 + x*z^-1)*x^-2*[1+x]^-1*z
       - 1/2*ln(z)*ln(1 + x*z^-1)*x^-2
       + ln(z)*ln(1 + x*z^-1)*x^-2*z
       - 2*ln(z)*ln(1 + x*z^-1)*x^-1*[1+x]^-1*z^2
       + 1/2*ln(z)*ln(1 + x*z^-1)*x^-1
       - ln(z)*ln(1 + x*z^-1)*x^-1*z
       + 2*ln(z)*ln(1 + x*z^-1)*x^-1*z^2
       + 1/2*ln(z)*ln(1 + x*z^-1)*[1+x]^-1
       - ln(z)*ln(1 + x*z^-1)*[1+x]^-1*z
       + ln(z)*ln(1 + x*z^-1)*z
       - 2*ln(z)*ln(1 + x*z^-1)*z^2
       + 1/2*ln(z)*ln(1 + x*z^-1)*x
       - ln(z)*ln(1 + x*z^-1)*x*z
       + ln(z)*ln(1 + x*z)
       - 1/2*ln(z)*ln(1 + x*z)*x^-2*[1+x]^-1
       + ln(z)*ln(1 + x*z)*x^-2*[1+x]^-1*z
       + 1/2*ln(z)*ln(1 + x*z)*x^-2
       - ln(z)*ln(1 + x*z)*x^-2*z
       + 2*ln(z)*ln(1 + x*z)*x^-1*[1+x]^-1*z^2
       - 1/2*ln(z)*ln(1 + x*z)*x^-1
       + ln(z)*ln(1 + x*z)*x^-1*z
       - 2*ln(z)*ln(1 + x*z)*x^-1*z^2
       - ln(z)*ln(1 + x*z)*[1-x]^-1
       - 2*ln(z)*ln(1 + x*z)*[1-x]^-1*z
       - 2*ln(z)*ln(1 + x*z)*[1-x]^-1*z^2
       - 1/2*ln(z)*ln(1 + x*z)*[1+x]^-1
       + ln(z)*ln(1 + x*z)*[1+x]^-1*z
       + 4*ln(z)*ln(1 + x*z)*z^2
       + 2*ln(z)*ln(1 + x*z)*x*z
       + 1/2*ln(z)*ln(z + x)
       - ln(z)*ln(z + x)*[1-x]^-1
       - 2*ln(z)*ln(z + x)*[1-x]^-1*z
       - 2*ln(z)*ln(z + x)*[1-x]^-1*z^2
       + ln(z)*ln(z + x)*z
       + 2*ln(z)*ln(z + x)*z^2
       + 1/2*ln(z)*ln(z + x)*x
       + ln(z)*ln(z + x)*x*z
       - 3/2*ln(z)^2
       + 1/4*ln(z)^2*x^-2*[1+x]^-1
       - 1/2*ln(z)^2*x^-2*[1+x]^-1*z
       - 1/4*ln(z)^2*x^-2
       + 1/2*ln(z)^2*x^-2*z
       - ln(z)^2*x^-1*[1+x]^-1*z^2
       + 1/4*ln(z)^2*x^-1
       - 1/2*ln(z)^2*x^-1*z
       + ln(z)^2*x^-1*z^2
       + 1/2*ln(z)^2*[1-x]^-1
       - ln(z)^2*[1-x]^-1*z
       + ln(z)^2*[1-x]^-1*z^2
       + 1/4*ln(z)^2*[1+x]^-1
       - 1/2*ln(z)^2*[1+x]^-1*z
       + ln(z)^2*z
       - 2*ln(z)^2*z^2
       - ln(z)^2*x*z
       + ln(z)*ln([1-z])
       - 1/2*ln(z)*ln([1-z])*[1-x]^-1
       + ln(z)*ln([1-z])*[1-x]^-1*z
       - 2*ln(z)*ln([1-z])*z
       - 13/4*ln([1-z])
       + 13/12*ln([1-z])*z^-1
       + 1/4*ln([1-z])*z
       + 2/3*ln([1-z])*z^2
       - 17/12*ln([1-z])*x*z^-1
       + 9/4*ln([1-z])*x
       + 3/4*ln([1-z])*x*z
       - 1/3*ln([1-z])*x*z^2
       - 4*ln(sqrtxz3)*ArcTan(sqrtxz3)*z*sqrtxz3
       + 1/4*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*sqrtxz2^-1*poly2^-1
       + 3/4*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*sqrtxz2^-1
       + 3/2*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x*sqrtxz2^-1
       - 3*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x*z*sqrtxz2^-1
       - 1/2*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x^2*sqrtxz2^-1*poly2^-1
       + 3/4*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x^2*sqrtxz2^-1
       + 1/4*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x^4*sqrtxz2^-1*poly2^-1
       - 1/4*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*sqrtxz2^-1*poly2^-1
       - 3/4*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*sqrtxz2^-1
       - 3/2*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x*sqrtxz2^-1
       + 3*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x*z*sqrtxz2^-1
       + 1/2*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x^2*sqrtxz2^-1*poly2^-1
       - 3/4*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x^2*sqrtxz2^-1
       - 1/4*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x^4*sqrtxz2^-1*poly2^-1
       - 1/2*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1
       - 3*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*z
       + 2*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z
       + Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x
       + 1/2*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1
       + 3*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*z
       - 2*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z
       - Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x
       - 1/4*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*sqrtxz2^-1*poly2^-1
       - 3/4*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*sqrtxz2^-1
       - 3/2*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x*sqrtxz2^-1
       + 3*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x*z*sqrtxz2^-1
       + 1/2*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x^2*sqrtxz2^-1*poly2^-1
       - 3/4*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x^2*sqrtxz2^-1
       - 1/4*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x^4*sqrtxz2^-1*poly2^-1
       + 1/4*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*sqrtxz2^-1*poly2^-1
       + 3/4*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*sqrtxz2^-1
       + 3/2*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x*sqrtxz2^-1
       - 3*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x*z*sqrtxz2^-1
       - 1/2*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x^2*sqrtxz2^-1*poly2^-1
       + 3/4*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x^2*sqrtxz2^-1
       + 1/4*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x^4*sqrtxz2^-1*poly2^-1
       - Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)
       + 3/2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*[1-x]^-1
       + Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*[1-x]^-1*z
       + 4*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*[1-x]^-1*z^2
       - 4*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*z^2
       - 2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*x*z
       + Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)
       - 3/2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*[1-x]^-1
       - Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*[1-x]^-1*z
       - 4*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*[1-x]^-1*z^2
       + 4*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*z^2
       + 2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*x*z
       - Li2(1 - x*z^-1)
       + 1/2*Li2(1 - x*z^-1)*[1-x]^-1
       - Li2(1 - x*z^-1)*[1-x]^-1*z
       + 2*Li2(1 - x*z^-1)*z
       - 1/2*Li2( - x*z^-1)*x^-2*[1+x]^-1
       + Li2( - x*z^-1)*x^-2*[1+x]^-1*z
       + 1/2*Li2( - x*z^-1)*x^-2
       - Li2( - x*z^-1)*x^-2*z
       + 2*Li2( - x*z^-1)*x^-1*[1+x]^-1*z^2
       - 1/2*Li2( - x*z^-1)*x^-1
       + Li2( - x*z^-1)*x^-1*z
       - 2*Li2( - x*z^-1)*x^-1*z^2
       + Li2( - x*z^-1)*[1-x]^-1
       + 2*Li2( - x*z^-1)*[1-x]^-1*z
       + 2*Li2( - x*z^-1)*[1-x]^-1*z^2
       - 1/2*Li2( - x*z^-1)*[1+x]^-1
       + Li2( - x*z^-1)*[1+x]^-1*z
       - 2*Li2( - x*z^-1)*z
       - Li2( - x*z^-1)*x
       + 2*Li2( - x)
       + Li2( - x)*x^-2*[1+x]^-1
       - 2*Li2( - x)*x^-2*[1+x]^-1*z
       - Li2( - x)*x^-2
       + 2*Li2( - x)*x^-2*z
       - 4*Li2( - x)*x^-1*[1+x]^-1*z^2
       + Li2( - x)*x^-1
       - 2*Li2( - x)*x^-1*z
       + 4*Li2( - x)*x^-1*z^2
       - Li2( - x)*[1+x]^-1
       + 2*Li2( - x)*[1+x]^-1*z
       - 4*Li2( - x)*[1+x]^-1*z^2
       - 4*Li2( - x)*z
       + 2*Li2( - x)*x
       - 4*Li2( - x)*x*z
       + Li2( - x*z)
       - 1/2*Li2( - x*z)*x^-2*[1+x]^-1
       + Li2( - x*z)*x^-2*[1+x]^-1*z
       + 1/2*Li2( - x*z)*x^-2
       - Li2( - x*z)*x^-2*z
       + 2*Li2( - x*z)*x^-1*[1+x]^-1*z^2
       - 1/2*Li2( - x*z)*x^-1
       + Li2( - x*z)*x^-1*z
       - 2*Li2( - x*z)*x^-1*z^2
       - Li2( - x*z)*[1-x]^-1
       - 2*Li2( - x*z)*[1-x]^-1*z
       - 2*Li2( - x*z)*[1-x]^-1*z^2
       - 1/2*Li2( - x*z)*[1+x]^-1
       + Li2( - x*z)*[1+x]^-1*z
       + 4*Li2( - x*z)*z^2
       + 2*Li2( - x*z)*x*z
       + 2*Li2(x)
       - Li2(x)*x^-2*[1+x]^-1
       + 2*Li2(x)*x^-2*[1+x]^-1*z
       + Li2(x)*x^-2
       - 2*Li2(x)*x^-2*z
       + 4*Li2(x)*x^-1*[1+x]^-1*z^2
       - Li2(x)*x^-1
       + 2*Li2(x)*x^-1*z
       - 4*Li2(x)*x^-1*z^2
       - 3/2*Li2(x)*[1-x]^-1
       + 3*Li2(x)*[1-x]^-1*z
       - 4*Li2(x)*[1-x]^-1*z^2
       + Li2(x)*[1+x]^-1
       - 2*Li2(x)*[1+x]^-1*z
       + 4*Li2(x)*[1+x]^-1*z^2
       - 1/2*Li2(x)*z^-1
       - 2*Li2(x)*z
       + 4*Li2(x)*z^2
       - 1/2*Li2(x)*x*z^-1
       + Li2(x)*x
       + 3/2*Li2(z)
       - 1/2*Li2(z)*[1-x]^-1
       + Li2(z)*[1-x]^-1*z
       - Li2(z)*z
       + 1/2*Li2(z)*x
       - 2*InvTanInt( - sqrtxz3)*z*sqrtxz3
       - 4*InvTanInt(z*sqrtxz3)*z*sqrtxz3
       + 2*InvTanInt(sqrtxz3)*z*sqrtxz3
       )

       + LMUF*NC^-1 * (  - 5/2
       + 5/4*z^-1
       - 5/4*x*z^-1
       + 5/2*x
       - ln(x)
       + 1/2*ln(x)*z^-1
       + 1/2*ln(x)*x*z^-1
       - ln(x)*x
       )

       + LMUF*NC * ( 5/2
       - 5/4*z^-1
       + 5/4*x*z^-1
       - 5/2*x
       + ln(x)
       - 1/2*ln(x)*z^-1
       - 1/2*ln(x)*x*z^-1
       + ln(x)*x
       )

       + LMUA*NC^-1 * (  - 3/4
       - 1/6*z^-1
       + 1/4*z
       + 2/3*z^2
       - 1/6*x*z^-1
       - 1/4*x
       + 3/4*x*z
       - 1/3*x*z^2
       - 1/2*ln(z)
       - ln(z)*z
       - 1/2*ln(z)*x
       )

       + LMUA*NC * ( 3/4
       + 1/6*z^-1
       - 1/4*z
       - 2/3*z^2
       + 1/6*x*z^-1
       + 1/4*x
       - 3/4*x*z
       + 1/3*x*z^2
       + 1/2*ln(z)
       + ln(z)*z
       + 1/2*ln(z)*x
       )

       + Dd([1-x])*NC^-2 * (  - 1/8
       - 2/3*[1+z]^-1*rln2^3
       + 1/3*rln2^3
       + 1/8*z
       - 1/3*z*rln2^3
       + 1/3*pi^2*[1+z]^-1*rln2
       + 1/24*pi^2
       - 1/6*pi^2*rln2
       + 1/24*pi^2*z
       + 1/6*pi^2*z*rln2
       + ln(1 + z)*[1+z]^-1*rln2^2
       - 1/2*ln(1 + z)*rln2^2
       + 1/2*ln(1 + z)*z*rln2^2
       - 1/6*ln(1 + z)*pi^2*[1+z]^-1
       + 1/12*ln(1 + z)*pi^2
       - 1/12*ln(1 + z)*pi^2*z
       + 7/8*ln(z)
       - 1/8*ln(z)*z
       - 1/2*ln(z)*ln(1 + z)
       - 1/2*ln(z)*ln(1 + z)*z
       + 1/2*ln(z)^2
       + 1/2*ln(z)^2*z
       + 3/4*ln(z)^2*ln(1 + z)
       - 3/2*ln(z)^2*ln(1 + z)*[1+z]^-1
       - 3/4*ln(z)^2*ln(1 + z)*z
       - 5/24*ln(z)^3
       + 5/12*ln(z)^3*[1+z]^-1
       + 5/24*ln(z)^3*z
       + ln(z)*ln([1-z])*ln(1 + z)
       - 2*ln(z)*ln([1-z])*ln(1 + z)*[1+z]^-1
       - ln(z)*ln([1-z])*ln(1 + z)*z
       + 1/2*ln(z)*Li2( - z)
       - ln(z)*Li2( - z)*[1+z]^-1
       - 1/2*ln(z)*Li2( - z)*z
       + ln([1-z])
       + ln([1-z])*[1+z]^-1*rln2^2
       - 1/2*ln([1-z])*rln2^2
       - ln([1-z])*z
       + 1/2*ln([1-z])*z*rln2^2
       - 1/6*ln([1-z])*pi^2*[1+z]^-1
       + 1/12*ln([1-z])*pi^2
       - 1/12*ln([1-z])*pi^2*z
       - 2*ln([1-z])*ln(1 + z)*[1+z]^-1*rln2
       + ln([1-z])*ln(1 + z)*rln2
       - ln([1-z])*ln(1 + z)*z*rln2
       - 1/2*ln([1-z])*ln(1 + z)^2
       + ln([1-z])*ln(1 + z)^2*[1+z]^-1
       + 1/2*ln([1-z])*ln(1 + z)^2*z
       + ln([1-z])*Li2( - z)
       - 2*ln([1-z])*Li2( - z)*[1+z]^-1
       - ln([1-z])*Li2( - z)*z
       - 1/3*ln([1+z])*pi^2*[1+z]^-1
       + 1/6*ln([1+z])*pi^2
       - 1/6*ln([1+z])*pi^2*z
       - Li3(1/2 - 1/2*z)
       + 2*Li3(1/2 - 1/2*z)*[1+z]^-1
       + Li3(1/2 - 1/2*z)*z
       - Li3(1/2 + 1/2*z)
       + 2*Li3(1/2 + 1/2*z)*[1+z]^-1
       + Li3(1/2 + 1/2*z)*z
       + Li3(1 - z)
       - 2*Li3(1 - z)*[1+z]^-1
       - Li3(1 - z)*z
       + 1/2*Li3( - z)
       - Li3( - z)*[1+z]^-1
       - 1/2*Li3( - z)*z
       + Li3(1/(1 + z))
       - 2*Li3(1/(1 + z))*[1+z]^-1
       - Li3(1/(1 + z))*z
       - Li3(1/(1 + z) - 1/(1 + z)*z)
       + 2*Li3(1/(1 + z) - 1/(1 + z)*z)*[1+z]^-1
       + Li3(1/(1 + z) - 1/(1 + z)*z)*z
       + 1/2*Li3(z)
       - Li3(z)*[1+z]^-1
       - 1/2*Li3(z)*z
       - 1/2*Li2( - z)
       - 1/2*Li2( - z)*z
       - 1/2*Li2(z)
       - 1/2*Li2(z)*z
       )

       + Dd([1-x])*NC^-1 * ( 11/18
       + 107/216*z^-1
       + 2/9*z
       - 287/216*z^2
       - zeta3
       - zeta3*z
       + 1/12*pi^2
       + 1/8*pi^2*z
       + 5/2*ln(z)
       + 7/36*ln(z)*z^-1
       + 9/4*ln(z)*z
       + 31/36*ln(z)*z^2
       - 1/12*ln(z)*pi^2
       - 1/12*ln(z)*pi^2*z
       + 1/16*ln(z)^2
       - 1/6*ln(z)^2*z^-1
       + 5/16*ln(z)^2*z
       - 5/24*ln(z)^3
       - 5/24*ln(z)^3*z
       + 3/4*ln(z)*ln([1-z])^2
       + 3/4*ln(z)*ln([1-z])^2*z
       + 5/3*ln([1-z])
       + 7/36*ln([1-z])*z^-1
       - 7/6*ln([1-z])*z
       - 25/36*ln([1-z])*z^2
       - 1/6*ln([1-z])*pi^2
       - 1/6*ln([1-z])*pi^2*z
       - 1/8*ln([1-z])^2
       - 1/6*ln([1-z])^2*z^-1
       + 1/8*ln([1-z])^2*z
       + 1/6*ln([1-z])^2*z^2
       + 1/2*ln([1-z])*Li2(1 - z)
       + 1/2*ln([1-z])*Li2(1 - z)*z
       + ln([1-z])*Li2(z)
       + ln([1-z])*Li2(z)*z
       + 1/2*Li3(1 - z)
       + 1/2*Li3(1 - z)*z
       + Li3(z)
       + Li3(z)*z
       - 1/4*Li2(z)
       + 1/3*Li2(z)*z^-1
       - Li2(z)*z
       - 1/3*Li2(z)*z^2
       )

       + Dd([1-x])*NC * (  - 11/18
       - 107/216*z^-1
       - 2/9*z
       + 287/216*z^2
       + zeta3
       + zeta3*z
       - 1/12*pi^2
       - 1/8*pi^2*z
       - 5/2*ln(z)
       - 7/36*ln(z)*z^-1
       - 9/4*ln(z)*z
       - 31/36*ln(z)*z^2
       + 1/12*ln(z)*pi^2
       + 1/12*ln(z)*pi^2*z
       - 1/16*ln(z)^2
       + 1/6*ln(z)^2*z^-1
       - 5/16*ln(z)^2*z
       + 5/24*ln(z)^3
       + 5/24*ln(z)^3*z
       - 3/4*ln(z)*ln([1-z])^2
       - 3/4*ln(z)*ln([1-z])^2*z
       - 5/3*ln([1-z])
       - 7/36*ln([1-z])*z^-1
       + 7/6*ln([1-z])*z
       + 25/36*ln([1-z])*z^2
       + 1/6*ln([1-z])*pi^2
       + 1/6*ln([1-z])*pi^2*z
       + 1/8*ln([1-z])^2
       + 1/6*ln([1-z])^2*z^-1
       - 1/8*ln([1-z])^2*z
       - 1/6*ln([1-z])^2*z^2
       - 1/2*ln([1-z])*Li2(1 - z)
       - 1/2*ln([1-z])*Li2(1 - z)*z
       - ln([1-z])*Li2(z)
       - ln([1-z])*Li2(z)*z
       - 1/2*Li3(1 - z)
       - 1/2*Li3(1 - z)*z
       - Li3(z)
       - Li3(z)*z
       + 1/4*Li2(z)
       - 1/3*Li2(z)*z^-1
       + Li2(z)*z
       + 1/3*Li2(z)*z^2
       )

       + Dd([1-x])*LMUA*NC^-2 * (  - 1
       + z
       + 1/6*pi^2*[1+z]^-1
       - 1/12*pi^2
       + 1/12*pi^2*z
       - 1/2*ln(z)
       - 1/2*ln(z)*z
       - ln(z)*ln(1 + z)
       + 2*ln(z)*ln(1 + z)*[1+z]^-1
       + ln(z)*ln(1 + z)*z
       + 1/4*ln(z)^2
       - 1/2*ln(z)^2*[1+z]^-1
       - 1/4*ln(z)^2*z
       - Li2( - z)
       + 2*Li2( - z)*[1+z]^-1
       + Li2( - z)*z
       )

       + Dd([1-x])*LMUA*NC^-1 * (  - 5/3
       - 7/36*z^-1
       + 7/6*z
       + 25/36*z^2
       + 1/12*pi^2
       + 1/12*pi^2*z
       - 1/4*ln(z)
       + 1/3*ln(z)*z^-1
       - ln(z)*z
       - 1/3*ln(z)*z^2
       + 1/2*ln(z)^2
       + 1/2*ln(z)^2*z
       + 1/4*ln([1-z])
       + 1/3*ln([1-z])*z^-1
       - 1/4*ln([1-z])*z
       - 1/3*ln([1-z])*z^2
       - 1/2*Li2(z)
       - 1/2*Li2(z)*z
       )

       + Dd([1-x])*LMUA * ( 1
       - z
       - 1/6*pi^2*[1+z]^-1
       + 1/12*pi^2
       - 1/12*pi^2*z
       + 1/2*ln(z)
       + 1/2*ln(z)*z
       + ln(z)*ln(1 + z)
       - 2*ln(z)*ln(1 + z)*[1+z]^-1
       - ln(z)*ln(1 + z)*z
       - 1/4*ln(z)^2
       + 1/2*ln(z)^2*[1+z]^-1
       + 1/4*ln(z)^2*z
       + Li2( - z)
       - 2*Li2( - z)*[1+z]^-1
       - Li2( - z)*z
       )

       + Dd([1-x])*LMUA*NC * ( 5/3
       + 7/36*z^-1
       - 7/6*z
       - 25/36*z^2
       - 1/12*pi^2
       - 1/12*pi^2*z
       + 1/4*ln(z)
       - 1/3*ln(z)*z^-1
       + ln(z)*z
       + 1/3*ln(z)*z^2
       - 1/2*ln(z)^2
       - 1/2*ln(z)^2*z
       - 1/4*ln([1-z])
       - 1/3*ln([1-z])*z^-1
       + 1/4*ln([1-z])*z
       + 1/3*ln([1-z])*z^2
       + 1/2*Li2(z)
       + 1/2*Li2(z)*z
       )

       + Dd([1-x])*LMUA^2*NC^-1 * (  - 1/8
       - 1/6*z^-1
       + 1/8*z
       + 1/6*z^2
       - 1/4*ln(z)
       - 1/4*ln(z)*z
       )

       + Dd([1-x])*LMUA^2*NC * ( 1/8
       + 1/6*z^-1
       - 1/8*z
       - 1/6*z^2
       + 1/4*ln(z)
       + 1/4*ln(z)*z
       )

       + Dd([1-x]) * ( 1/8
       + 2/3*[1+z]^-1*rln2^3
       - 1/3*rln2^3
       - 1/8*z
       + 1/3*z*rln2^3
       - 1/3*pi^2*[1+z]^-1*rln2
       - 1/24*pi^2
       + 1/6*pi^2*rln2
       - 1/24*pi^2*z
       - 1/6*pi^2*z*rln2
       - ln(1 + z)*[1+z]^-1*rln2^2
       + 1/2*ln(1 + z)*rln2^2
       - 1/2*ln(1 + z)*z*rln2^2
       + 1/6*ln(1 + z)*pi^2*[1+z]^-1
       - 1/12*ln(1 + z)*pi^2
       + 1/12*ln(1 + z)*pi^2*z
       - 7/8*ln(z)
       + 1/8*ln(z)*z
       + 1/2*ln(z)*ln(1 + z)
       + 1/2*ln(z)*ln(1 + z)*z
       - 1/2*ln(z)^2
       - 1/2*ln(z)^2*z
       - 3/4*ln(z)^2*ln(1 + z)
       + 3/2*ln(z)^2*ln(1 + z)*[1+z]^-1
       + 3/4*ln(z)^2*ln(1 + z)*z
       + 5/24*ln(z)^3
       - 5/12*ln(z)^3*[1+z]^-1
       - 5/24*ln(z)^3*z
       - ln(z)*ln([1-z])*ln(1 + z)
       + 2*ln(z)*ln([1-z])*ln(1 + z)*[1+z]^-1
       + ln(z)*ln([1-z])*ln(1 + z)*z
       - 1/2*ln(z)*Li2( - z)
       + ln(z)*Li2( - z)*[1+z]^-1
       + 1/2*ln(z)*Li2( - z)*z
       - ln([1-z])
       - ln([1-z])*[1+z]^-1*rln2^2
       + 1/2*ln([1-z])*rln2^2
       + ln([1-z])*z
       - 1/2*ln([1-z])*z*rln2^2
       + 1/6*ln([1-z])*pi^2*[1+z]^-1
       - 1/12*ln([1-z])*pi^2
       + 1/12*ln([1-z])*pi^2*z
       + 2*ln([1-z])*ln(1 + z)*[1+z]^-1*rln2
       - ln([1-z])*ln(1 + z)*rln2
       + ln([1-z])*ln(1 + z)*z*rln2
       + 1/2*ln([1-z])*ln(1 + z)^2
       - ln([1-z])*ln(1 + z)^2*[1+z]^-1
       - 1/2*ln([1-z])*ln(1 + z)^2*z
       - ln([1-z])*Li2( - z)
       + 2*ln([1-z])*Li2( - z)*[1+z]^-1
       + ln([1-z])*Li2( - z)*z
       + 1/3*ln([1+z])*pi^2*[1+z]^-1
       - 1/6*ln([1+z])*pi^2
       + 1/6*ln([1+z])*pi^2*z
       + Li3(1/2 - 1/2*z)
       - 2*Li3(1/2 - 1/2*z)*[1+z]^-1
       - Li3(1/2 - 1/2*z)*z
       + Li3(1/2 + 1/2*z)
       - 2*Li3(1/2 + 1/2*z)*[1+z]^-1
       - Li3(1/2 + 1/2*z)*z
       - Li3(1 - z)
       + 2*Li3(1 - z)*[1+z]^-1
       + Li3(1 - z)*z
       - 1/2*Li3( - z)
       + Li3( - z)*[1+z]^-1
       + 1/2*Li3( - z)*z
       - Li3(1/(1 + z))
       + 2*Li3(1/(1 + z))*[1+z]^-1
       + Li3(1/(1 + z))*z
       + Li3(1/(1 + z) - 1/(1 + z)*z)
       - 2*Li3(1/(1 + z) - 1/(1 + z)*z)*[1+z]^-1
       - Li3(1/(1 + z) - 1/(1 + z)*z)*z
       - 1/2*Li3(z)
       + Li3(z)*[1+z]^-1
       + 1/2*Li3(z)*z
       + 1/2*Li2( - z)
       + 1/2*Li2( - z)*z
       + 1/2*Li2(z)
       + 1/2*Li2(z)*z
       )

       + Dd([1-z])*NC^-2 * ( 13/8
       + 2/3*[1+x]^-1*rln2^3
       - 1/3*rln2^3
       - 13/8*x
       + 1/3*x*rln2^3
       + zeta3*[1+x]^-1
       - 1/2*zeta3
       + 1/2*zeta3*x
       - 1/3*pi^2*[1+x]^-1*rln2
       - 5/24*pi^2
       + 1/6*pi^2*rln2
       - 5/24*pi^2*x
       - 1/6*pi^2*x*rln2
       - ln(1 + x)*[1+x]^-1*rln2^2
       + 1/2*ln(1 + x)*rln2^2
       - 1/2*ln(1 + x)*x*rln2^2
       - 1/6*ln(1 + x)*pi^2*[1+x]^-1
       + 1/12*ln(1 + x)*pi^2
       - 1/12*ln(1 + x)*pi^2*x
       - 1/6*ln(1 + x)^3
       + 1/3*ln(1 + x)^3*[1+x]^-1
       + 1/6*ln(1 + x)^3*x
       + 9/8*ln(x)
       - 7/8*ln(x)*x
       - 1/3*ln(x)*pi^2*[1+x]^-1
       + 1/6*ln(x)*pi^2
       - 1/6*ln(x)*pi^2*x
       - 3/2*ln(x)*ln(1 + x)
       - 3/2*ln(x)*ln(1 + x)*x
       + 1/2*ln(x)^2
       + ln(x)^2*x
       + ln(x)^2*ln(1 + x)
       - 2*ln(x)^2*ln(1 + x)*[1+x]^-1
       - ln(x)^2*ln(1 + x)*x
       - 1/8*ln(x)^3
       + 1/4*ln(x)^3*[1+x]^-1
       + 1/8*ln(x)^3*x
       - ln(x)*ln([1-x])*ln(1 + x)
       + 2*ln(x)*ln([1-x])*ln(1 + x)*[1+x]^-1
       + ln(x)*ln([1-x])*ln(1 + x)*x
       + ln(x)*Li2( - x)
       - 2*ln(x)*Li2( - x)*[1+x]^-1
       - ln(x)*Li2( - x)*x
       - ln([1-x])
       - ln([1-x])*[1+x]^-1*rln2^2
       + 1/2*ln([1-x])*rln2^2
       + ln([1-x])*x
       - 1/2*ln([1-x])*x*rln2^2
       + 1/6*ln([1-x])*pi^2*[1+x]^-1
       - 1/12*ln([1-x])*pi^2
       + 1/12*ln([1-x])*pi^2*x
       + 2*ln([1-x])*ln(1 + x)*[1+x]^-1*rln2
       - ln([1-x])*ln(1 + x)*rln2
       + ln([1-x])*ln(1 + x)*x*rln2
       + 1/2*ln([1-x])*ln(1 + x)^2
       - ln([1-x])*ln(1 + x)^2*[1+x]^-1
       - 1/2*ln([1-x])*ln(1 + x)^2*x
       - ln([1-x])*Li2( - x)
       + 2*ln([1-x])*Li2( - x)*[1+x]^-1
       + ln([1-x])*Li2( - x)*x
       + 1/2*ln([1+x])*pi^2*[1+x]^-1
       - 1/4*ln([1+x])*pi^2
       + 1/4*ln([1+x])*pi^2*x
       + Li3(1/2 - 1/2*x)
       - 2*Li3(1/2 - 1/2*x)*[1+x]^-1
       - Li3(1/2 - 1/2*x)*x
       + Li3(1/2 + 1/2*x)
       - 2*Li3(1/2 + 1/2*x)*[1+x]^-1
       - Li3(1/2 + 1/2*x)*x
       - Li3(1 - x)
       + 2*Li3(1 - x)*[1+x]^-1
       + Li3(1 - x)*x
       + Li3(1/(1 + x) - 1/(1 + x)*x)
       - 2*Li3(1/(1 + x) - 1/(1 + x)*x)*[1+x]^-1
       - Li3(1/(1 + x) - 1/(1 + x)*x)*x
       - 1/2*Li3(x)
       + Li3(x)*[1+x]^-1
       + 1/2*Li3(x)*x
       - 3/2*Li2( - x)
       - 3/2*Li2( - x)*x
       + 1/2*Li2(x)
       + 1/2*Li2(x)*x
       )

       + Dd([1-z])*NC^-1 * (  - 9/2
       + 9/2*x
       + 1/4*pi^2
       - 1/4*pi^2*x
       - 25/8*ln(x)
       + 3/8*ln(x)*x
       + 1/6*ln(x)*pi^2
       + 1/6*ln(x)*pi^2*x
       - 17/16*ln(x)^2
       + 15/16*ln(x)^2*x
       - 5/24*ln(x)^3
       - 5/24*ln(x)^3*x
       + 5/4*ln(x)*ln([1-x])
       - 5/4*ln(x)*ln([1-x])*x
       + 3/4*ln(x)*ln([1-x])^2
       + 3/4*ln(x)*ln([1-x])^2*x
       - 1/2*ln(x)*Li2(x)
       - 1/2*ln(x)*Li2(x)*x
       + 3/2*ln([1-x])
       - 3/2*ln([1-x])*x
       - 1/6*ln([1-x])*pi^2
       - 1/6*ln([1-x])*pi^2*x
       - 5/8*ln([1-x])^2
       + 5/8*ln([1-x])^2*x
       + 1/2*ln([1-x])*Li2(1 - x)
       + 1/2*ln([1-x])*Li2(1 - x)*x
       + ln([1-x])*Li2(x)
       + ln([1-x])*Li2(x)*x
       + 1/2*Li3(1 - x)
       + 1/2*Li3(1 - x)*x
       - 1/4*Li2(x)
       + 1/4*Li2(x)*x
       )

       + Dd([1-z])*NC * ( 9/2
       - 9/2*x
       - 1/4*pi^2
       + 1/4*pi^2*x
       + 25/8*ln(x)
       - 3/8*ln(x)*x
       - 1/6*ln(x)*pi^2
       - 1/6*ln(x)*pi^2*x
       + 17/16*ln(x)^2
       - 15/16*ln(x)^2*x
       + 5/24*ln(x)^3
       + 5/24*ln(x)^3*x
       - 5/4*ln(x)*ln([1-x])
       + 5/4*ln(x)*ln([1-x])*x
       - 3/4*ln(x)*ln([1-x])^2
       - 3/4*ln(x)*ln([1-x])^2*x
       + 1/2*ln(x)*Li2(x)
       + 1/2*ln(x)*Li2(x)*x
       - 3/2*ln([1-x])
       + 3/2*ln([1-x])*x
       + 1/6*ln([1-x])*pi^2
       + 1/6*ln([1-x])*pi^2*x
       + 5/8*ln([1-x])^2
       - 5/8*ln([1-x])^2*x
       - 1/2*ln([1-x])*Li2(1 - x)
       - 1/2*ln([1-x])*Li2(1 - x)*x
       - ln([1-x])*Li2(x)
       - ln([1-x])*Li2(x)*x
       - 1/2*Li3(1 - x)
       - 1/2*Li3(1 - x)*x
       + 1/4*Li2(x)
       - 1/4*Li2(x)*x
       )

       + Dd([1-z])*LMUF*NC^-2 * ( 1
       - x
       - 1/6*pi^2*[1+x]^-1
       + 1/12*pi^2
       - 1/12*pi^2*x
       + 1/2*ln(x)
       + 1/2*ln(x)*x
       + ln(x)*ln(1 + x)
       - 2*ln(x)*ln(1 + x)*[1+x]^-1
       - ln(x)*ln(1 + x)*x
       - 1/4*ln(x)^2
       + 1/2*ln(x)^2*[1+x]^-1
       + 1/4*ln(x)^2*x
       + Li2( - x)
       - 2*Li2( - x)*[1+x]^-1
       - Li2( - x)*x
       )

       + Dd([1-z])*LMUF*NC^-1 * (  - 3/2
       + 3/2*x
       + 1/12*pi^2
       + 1/12*pi^2*x
       - 3/2*ln(x)
       + 3/2*ln(x)*x
       - 1/2*ln(x)^2
       - 1/2*ln(x)^2*x
       + 5/4*ln([1-x])
       - 5/4*ln([1-x])*x
       - 1/2*Li2(x)
       - 1/2*Li2(x)*x
       )

       + Dd([1-z])*LMUF * (  - 1
       + x
       + 1/6*pi^2*[1+x]^-1
       - 1/12*pi^2
       + 1/12*pi^2*x
       - 1/2*ln(x)
       - 1/2*ln(x)*x
       - ln(x)*ln(1 + x)
       + 2*ln(x)*ln(1 + x)*[1+x]^-1
       + ln(x)*ln(1 + x)*x
       + 1/4*ln(x)^2
       - 1/2*ln(x)^2*[1+x]^-1
       - 1/4*ln(x)^2*x
       - Li2( - x)
       + 2*Li2( - x)*[1+x]^-1
       + Li2( - x)*x
       )

       + Dd([1-z])*LMUF*NC * ( 3/2
       - 3/2*x
       - 1/12*pi^2
       - 1/12*pi^2*x
       + 3/2*ln(x)
       - 3/2*ln(x)*x
       + 1/2*ln(x)^2
       + 1/2*ln(x)^2*x
       - 5/4*ln([1-x])
       + 5/4*ln([1-x])*x
       + 1/2*Li2(x)
       + 1/2*Li2(x)*x
       )

       + Dd([1-z])*LMUF^2*NC^-1 * (  - 5/8
       + 5/8*x
       - 1/4*ln(x)
       - 1/4*ln(x)*x
       )

       + Dd([1-z])*LMUF^2*NC * ( 5/8
       - 5/8*x
       + 1/4*ln(x)
       + 1/4*ln(x)*x
       )

       + Dd([1-z]) * (  - 13/8
       - 2/3*[1+x]^-1*rln2^3
       + 1/3*rln2^3
       + 13/8*x
       - 1/3*x*rln2^3
       - zeta3*[1+x]^-1
       + 1/2*zeta3
       - 1/2*zeta3*x
       + 1/3*pi^2*[1+x]^-1*rln2
       + 5/24*pi^2
       - 1/6*pi^2*rln2
       + 5/24*pi^2*x
       + 1/6*pi^2*x*rln2
       + ln(1 + x)*[1+x]^-1*rln2^2
       - 1/2*ln(1 + x)*rln2^2
       + 1/2*ln(1 + x)*x*rln2^2
       + 1/6*ln(1 + x)*pi^2*[1+x]^-1
       - 1/12*ln(1 + x)*pi^2
       + 1/12*ln(1 + x)*pi^2*x
       + 1/6*ln(1 + x)^3
       - 1/3*ln(1 + x)^3*[1+x]^-1
       - 1/6*ln(1 + x)^3*x
       - 9/8*ln(x)
       + 7/8*ln(x)*x
       + 1/3*ln(x)*pi^2*[1+x]^-1
       - 1/6*ln(x)*pi^2
       + 1/6*ln(x)*pi^2*x
       + 3/2*ln(x)*ln(1 + x)
       + 3/2*ln(x)*ln(1 + x)*x
       - 1/2*ln(x)^2
       - ln(x)^2*x
       - ln(x)^2*ln(1 + x)
       + 2*ln(x)^2*ln(1 + x)*[1+x]^-1
       + ln(x)^2*ln(1 + x)*x
       + 1/8*ln(x)^3
       - 1/4*ln(x)^3*[1+x]^-1
       - 1/8*ln(x)^3*x
       + ln(x)*ln([1-x])*ln(1 + x)
       - 2*ln(x)*ln([1-x])*ln(1 + x)*[1+x]^-1
       - ln(x)*ln([1-x])*ln(1 + x)*x
       - ln(x)*Li2( - x)
       + 2*ln(x)*Li2( - x)*[1+x]^-1
       + ln(x)*Li2( - x)*x
       + ln([1-x])
       + ln([1-x])*[1+x]^-1*rln2^2
       - 1/2*ln([1-x])*rln2^2
       - ln([1-x])*x
       + 1/2*ln([1-x])*x*rln2^2
       - 1/6*ln([1-x])*pi^2*[1+x]^-1
       + 1/12*ln([1-x])*pi^2
       - 1/12*ln([1-x])*pi^2*x
       - 2*ln([1-x])*ln(1 + x)*[1+x]^-1*rln2
       + ln([1-x])*ln(1 + x)*rln2
       - ln([1-x])*ln(1 + x)*x*rln2
       - 1/2*ln([1-x])*ln(1 + x)^2
       + ln([1-x])*ln(1 + x)^2*[1+x]^-1
       + 1/2*ln([1-x])*ln(1 + x)^2*x
       + ln([1-x])*Li2( - x)
       - 2*ln([1-x])*Li2( - x)*[1+x]^-1
       - ln([1-x])*Li2( - x)*x
       - 1/2*ln([1+x])*pi^2*[1+x]^-1
       + 1/4*ln([1+x])*pi^2
       - 1/4*ln([1+x])*pi^2*x
       - Li3(1/2 - 1/2*x)
       + 2*Li3(1/2 - 1/2*x)*[1+x]^-1
       + Li3(1/2 - 1/2*x)*x
       - Li3(1/2 + 1/2*x)
       + 2*Li3(1/2 + 1/2*x)*[1+x]^-1
       + Li3(1/2 + 1/2*x)*x
       + Li3(1 - x)
       - 2*Li3(1 - x)*[1+x]^-1
       - Li3(1 - x)*x
       - Li3(1/(1 + x) - 1/(1 + x)*x)
       + 2*Li3(1/(1 + x) - 1/(1 + x)*x)*[1+x]^-1
       + Li3(1/(1 + x) - 1/(1 + x)*x)*x
       + 1/2*Li3(x)
       - Li3(x)*[1+x]^-1
       - 1/2*Li3(x)*x
       + 3/2*Li2( - x)
       + 3/2*Li2( - x)*x
       - 1/2*Li2(x)
       - 1/2*Li2(x)*x
       )

       + Dn(0,[1-x])*NC^-2 * ( 1
       - z
       - 1/6*pi^2*[1+z]^-1
       + 1/12*pi^2
       - 1/12*pi^2*z
       + 1/2*ln(z)
       + 1/2*ln(z)*z
       + ln(z)*ln(1 + z)
       - 2*ln(z)*ln(1 + z)*[1+z]^-1
       - ln(z)*ln(1 + z)*z
       - 1/4*ln(z)^2
       + 1/2*ln(z)^2*[1+z]^-1
       + 1/4*ln(z)^2*z
       + Li2( - z)
       - 2*Li2( - z)*[1+z]^-1
       - Li2( - z)*z
       )

       + Dn(0,[1-x])*NC^-1 * ( 5/3
       + 7/36*z^-1
       - 7/6*z
       - 25/36*z^2
       - 1/12*pi^2
       - 1/12*pi^2*z
       + 1/4*ln(z)
       - 1/3*ln(z)*z^-1
       + ln(z)*z
       + 1/3*ln(z)*z^2
       - 1/2*ln(z)^2
       - 1/2*ln(z)^2*z
       - 1/4*ln([1-z])
       - 1/3*ln([1-z])*z^-1
       + 1/4*ln([1-z])*z
       + 1/3*ln([1-z])*z^2
       + 1/2*Li2(z)
       + 1/2*Li2(z)*z
       )

       + Dn(0,[1-x])*NC * (  - 5/3
       - 7/36*z^-1
       + 7/6*z
       + 25/36*z^2
       + 1/12*pi^2
       + 1/12*pi^2*z
       - 1/4*ln(z)
       + 1/3*ln(z)*z^-1
       - ln(z)*z
       - 1/3*ln(z)*z^2
       + 1/2*ln(z)^2
       + 1/2*ln(z)^2*z
       + 1/4*ln([1-z])
       + 1/3*ln([1-z])*z^-1
       - 1/4*ln([1-z])*z
       - 1/3*ln([1-z])*z^2
       - 1/2*Li2(z)
       - 1/2*Li2(z)*z
       )

       + Dn(0,[1-x])*LMUA*NC^-1 * ( 1/4
       + 1/3*z^-1
       - 1/4*z
       - 1/3*z^2
       + 1/2*ln(z)
       + 1/2*ln(z)*z
       )

       + Dn(0,[1-x])*LMUA*NC * (  - 1/4
       - 1/3*z^-1
       + 1/4*z
       + 1/3*z^2
       - 1/2*ln(z)
       - 1/2*ln(z)*z
       )

       + Dn(0,[1-x]) * (  - 1
       + z
       + 1/6*pi^2*[1+z]^-1
       - 1/12*pi^2
       + 1/12*pi^2*z
       - 1/2*ln(z)
       - 1/2*ln(z)*z
       - ln(z)*ln(1 + z)
       + 2*ln(z)*ln(1 + z)*[1+z]^-1
       + ln(z)*ln(1 + z)*z
       + 1/4*ln(z)^2
       - 1/2*ln(z)^2*[1+z]^-1
       - 1/4*ln(z)^2*z
       - Li2( - z)
       + 2*Li2( - z)*[1+z]^-1
       + Li2( - z)*z
       )

       + Dn(0,[1-z])*NC^-2 * (  - 1
       + x
       + 1/6*pi^2*[1+x]^-1
       - 1/12*pi^2
       + 1/12*pi^2*x
       - 1/2*ln(x)
       - 1/2*ln(x)*x
       - ln(x)*ln(1 + x)
       + 2*ln(x)*ln(1 + x)*[1+x]^-1
       + ln(x)*ln(1 + x)*x
       + 1/4*ln(x)^2
       - 1/2*ln(x)^2*[1+x]^-1
       - 1/4*ln(x)^2*x
       - Li2( - x)
       + 2*Li2( - x)*[1+x]^-1
       + Li2( - x)*x
       )

       + Dn(0,[1-z])*NC^-1 * ( 3/2
       - 3/2*x
       - 1/12*pi^2
       - 1/12*pi^2*x
       + 3/2*ln(x)
       - 3/2*ln(x)*x
       + 1/2*ln(x)^2
       + 1/2*ln(x)^2*x
       - 5/4*ln([1-x])
       + 5/4*ln([1-x])*x
       + 1/2*Li2(x)
       + 1/2*Li2(x)*x
       )

       + Dn(0,[1-z])*NC * (  - 3/2
       + 3/2*x
       + 1/12*pi^2
       + 1/12*pi^2*x
       - 3/2*ln(x)
       + 3/2*ln(x)*x
       - 1/2*ln(x)^2
       - 1/2*ln(x)^2*x
       + 5/4*ln([1-x])
       - 5/4*ln([1-x])*x
       - 1/2*Li2(x)
       - 1/2*Li2(x)*x
       )

       + Dn(0,[1-z])*LMUF*NC^-1 * ( 5/4
       - 5/4*x
       + 1/2*ln(x)
       + 1/2*ln(x)*x
       )

       + Dn(0,[1-z])*LMUF*NC * (  - 5/4
       + 5/4*x
       - 1/2*ln(x)
       - 1/2*ln(x)*x
       )

       + Dn(0,[1-z]) * ( 1
       - x
       - 1/6*pi^2*[1+x]^-1
       + 1/12*pi^2
       - 1/12*pi^2*x
       + 1/2*ln(x)
       + 1/2*ln(x)*x
       + ln(x)*ln(1 + x)
       - 2*ln(x)*ln(1 + x)*[1+x]^-1
       - ln(x)*ln(1 + x)*x
       - 1/4*ln(x)^2
       + 1/2*ln(x)^2*[1+x]^-1
       + 1/4*ln(x)^2*x
       + Li2( - x)
       - 2*Li2( - x)*[1+x]^-1
       - Li2( - x)*x
       )

       + Dn(1,[1-x])*NC^-1 * (  - 1/4
       - 1/3*z^-1
       + 1/4*z
       + 1/3*z^2
       - 1/2*ln(z)
       - 1/2*ln(z)*z
       )

       + Dn(1,[1-x])*NC * ( 1/4
       + 1/3*z^-1
       - 1/4*z
       - 1/3*z^2
       + 1/2*ln(z)
       + 1/2*ln(z)*z
       )

       + Dn(1,[1-z])*NC^-1 * (  - 5/4
       + 5/4*x
       - 1/2*ln(x)
       - 1/2*ln(x)*x
       )

       + Dn(1,[1-z])*NC * ( 5/4
       - 5/4*x
       + 1/2*ln(x)
       + 1/2*ln(x)*x
       )

       - 6
       + 4*[1-x]^-1*z^-1*[1+z]^-1*rln2^2
       + 8*[1-x]^-1*z^-1*[1+z]^-1*sqrtxz1*rln2
       + 6*[1-x]^-1*z^-1*[1+z]^-1*sqrtxz1*rln2^2
       - 4*[1-x]^-1*z^-1*rln2^2
       - 8*[1-x]^-1*z^-1*sqrtxz1*rln2
       - 6*[1-x]^-1*z^-1*sqrtxz1*rln2^2
       - 8*[1-x]^-1*[1+z]^-1*sqrtxz1^-1*rln2
       - 6*[1-x]^-1*[1+z]^-1*sqrtxz1^-1*rln2^2
       + 16*[1-x]^-1*sqrtxz1^-1*rln2
       + 12*[1-x]^-1*sqrtxz1^-1*rln2^2
       + 3*[1-x]^-1*rln2^2
       - [1-x]^-1*sqrtxz1*rln2
       - 8*[1-x]^-1*z*[1+z]^-1*sqrtxz1^-1*rln2
       - 6*[1-x]^-1*z*[1+z]^-1*sqrtxz1^-1*rln2^2
       + 8*[1-x]^-1*z*sqrtxz1^-1*rln2
       + 6*[1-x]^-1*z*sqrtxz1^-1*rln2^2
       - 3*[1-x]^-1*z*rln2^2
       - 4*z^-1*[1+z]^-1*rln2^2
       - 8*z^-1*[1+z]^-1*sqrtxz1*rln2
       - 6*z^-1*[1+z]^-1*sqrtxz1*rln2^2
       + 2*z^-1
       + 4*z^-1*rln2^2
       + 8*z^-1*sqrtxz1*rln2
       + 6*z^-1*sqrtxz1*rln2^2
       + 8*[1+z]^-1*sqrtxz1^-1*rln2
       + 6*[1+z]^-1*sqrtxz1^-1*rln2^2
       - 16*sqrtxz1^-1*rln2
       - 12*sqrtxz1^-1*rln2^2
       - 4*rln2^2
       + 8*z*[1+z]^-1*sqrtxz1^-1*rln2
       + 6*z*[1+z]^-1*sqrtxz1^-1*rln2^2
       - 8*z*sqrtxz1^-1*rln2
       - 6*z*sqrtxz1^-1*rln2^2
       + z
       + 4*z*rln2^2
       - 7/4*x*z^-1
       + 32*x*[1+z]^-1*sqrtxz1^-1*rln2
       + 24*x*[1+z]^-1*sqrtxz1^-1*rln2^2
       - 32*x*sqrtxz1^-1*rln2
       - 24*x*sqrtxz1^-1*rln2^2
       + 11/2*x
       - 3/4*x*z
       - 1/3*pi^2*x^-1*[1+x]^-1*z^-1*[1-z]^-1
       + 1/3*pi^2*x^-1*[1+x]^-1*z^-1
       + 1/3*pi^2*x^-1*[1+x]^-1
       + 1/3*pi^2*x^-1*[1+x]^-1*z
       + 1/3*pi^2*x^-1*z^-1*[1-z]^-1
       - 1/3*pi^2*x^-1*z^-1
       - 1/3*pi^2*x^-1
       - 1/3*pi^2*x^-1*z
       + 1/3*pi^2*[1-x]^-1*z^-1*[1+z]^-1*sqrtxz1
       + 1/6*pi^2*[1-x]^-1*z^-1
       - 1/3*pi^2*[1-x]^-1*z^-1*sqrtxz1
       - 1/3*pi^2*[1-x]^-1*[1+z]^-1*sqrtxz1^-1
       - 1/6*pi^2*[1-x]^-1*[1+z]^-1
       + 2/3*pi^2*[1-x]^-1*sqrtxz1^-1
       - 1/4*pi^2*[1-x]^-1
       - 1/3*pi^2*[1-x]^-1*z*[1+z]^-1*sqrtxz1^-1
       + 1/3*pi^2*[1-x]^-1*z*sqrtxz1^-1
       + 1/12*pi^2*[1-x]^-1*z
       - 1/3*pi^2*[1+x]^-1*[1-z]^-1
       + 1/2*pi^2*[1+x]^-1
       + 1/2*pi^2*[1+x]^-1*z
       - 1/3*pi^2*z^-1*[1-z]^-1
       - 1/3*pi^2*z^-1*[1+z]^-1*sqrtxz1
       + 1/4*pi^2*z^-1
       + 1/3*pi^2*z^-1*sqrtxz1
       + 1/3*pi^2*[1-z]^-1
       + 1/3*pi^2*[1+z]^-1*sqrtxz1^-1
       - 2/3*pi^2*sqrtxz1^-1
       + 1/6*pi^2
       + 1/3*pi^2*z*[1+z]^-1*sqrtxz1^-1
       - 1/3*pi^2*z*sqrtxz1^-1
       - 1/3*pi^2*z
       - 1/12*pi^2*x*z^-1
       + 4/3*pi^2*x*[1+z]^-1*sqrtxz1^-1
       - 4/3*pi^2*x*sqrtxz1^-1
       + 1/3*pi^2*x
       - 4*ln(1 + sqrtxz1 - z)*[1-x]^-1*z^-1*[1+z]^-1*rln2
       - 8*ln(1 + sqrtxz1 - z)*[1-x]^-1*z^-1*[1+z]^-1*sqrtxz1
       - 8*ln(1 + sqrtxz1 - z)*[1-x]^-1*z^-1*[1+z]^-1*sqrtxz1*rln2
       + 4*ln(1 + sqrtxz1 - z)*[1-x]^-1*z^-1*rln2
       + 8*ln(1 + sqrtxz1 - z)*[1-x]^-1*z^-1*sqrtxz1
       + 8*ln(1 + sqrtxz1 - z)*[1-x]^-1*z^-1*sqrtxz1*rln2
       + 8*ln(1 + sqrtxz1 - z)*[1-x]^-1*[1+z]^-1*sqrtxz1^-1
       + 8*ln(1 + sqrtxz1 - z)*[1-x]^-1*[1+z]^-1*sqrtxz1^-1*rln2
       - 16*ln(1 + sqrtxz1 - z)*[1-x]^-1*sqrtxz1^-1
       - 16*ln(1 + sqrtxz1 - z)*[1-x]^-1*sqrtxz1^-1*rln2
       - 3*ln(1 + sqrtxz1 - z)*[1-x]^-1*rln2
       + ln(1 + sqrtxz1 - z)*[1-x]^-1*sqrtxz1
       + 8*ln(1 + sqrtxz1 - z)*[1-x]^-1*z*[1+z]^-1*sqrtxz1^-1
       + 8*ln(1 + sqrtxz1 - z)*[1-x]^-1*z*[1+z]^-1*sqrtxz1^-1*rln2
       - 8*ln(1 + sqrtxz1 - z)*[1-x]^-1*z*sqrtxz1^-1
       - 8*ln(1 + sqrtxz1 - z)*[1-x]^-1*z*sqrtxz1^-1*rln2
       + 3*ln(1 + sqrtxz1 - z)*[1-x]^-1*z*rln2
       + 4*ln(1 + sqrtxz1 - z)*z^-1*[1+z]^-1*rln2
       + 8*ln(1 + sqrtxz1 - z)*z^-1*[1+z]^-1*sqrtxz1
       + 8*ln(1 + sqrtxz1 - z)*z^-1*[1+z]^-1*sqrtxz1*rln2
       - 4*ln(1 + sqrtxz1 - z)*z^-1*rln2
       - 8*ln(1 + sqrtxz1 - z)*z^-1*sqrtxz1
       - 8*ln(1 + sqrtxz1 - z)*z^-1*sqrtxz1*rln2
       - 8*ln(1 + sqrtxz1 - z)*[1+z]^-1*sqrtxz1^-1
       - 8*ln(1 + sqrtxz1 - z)*[1+z]^-1*sqrtxz1^-1*rln2
       + 16*ln(1 + sqrtxz1 - z)*sqrtxz1^-1
       + 16*ln(1 + sqrtxz1 - z)*sqrtxz1^-1*rln2
       + 4*ln(1 + sqrtxz1 - z)*rln2
       - 8*ln(1 + sqrtxz1 - z)*z*[1+z]^-1*sqrtxz1^-1
       - 8*ln(1 + sqrtxz1 - z)*z*[1+z]^-1*sqrtxz1^-1*rln2
       + 8*ln(1 + sqrtxz1 - z)*z*sqrtxz1^-1
       + 8*ln(1 + sqrtxz1 - z)*z*sqrtxz1^-1*rln2
       - 4*ln(1 + sqrtxz1 - z)*z*rln2
       - 32*ln(1 + sqrtxz1 - z)*x*[1+z]^-1*sqrtxz1^-1
       - 32*ln(1 + sqrtxz1 - z)*x*[1+z]^-1*sqrtxz1^-1*rln2
       + 32*ln(1 + sqrtxz1 - z)*x*sqrtxz1^-1
       + 32*ln(1 + sqrtxz1 - z)*x*sqrtxz1^-1*rln2
       + 2*ln(1 + sqrtxz1 - z)^2*[1-x]^-1*z^-1*[1+z]^-1*sqrtxz1
       - 2*ln(1 + sqrtxz1 - z)^2*[1-x]^-1*z^-1*sqrtxz1
       - 2*ln(1 + sqrtxz1 - z)^2*[1-x]^-1*[1+z]^-1*sqrtxz1^-1
       + 4*ln(1 + sqrtxz1 - z)^2*[1-x]^-1*sqrtxz1^-1
       - 2*ln(1 + sqrtxz1 - z)^2*[1-x]^-1*z*[1+z]^-1*sqrtxz1^-1
       + 2*ln(1 + sqrtxz1 - z)^2*[1-x]^-1*z*sqrtxz1^-1
       - 2*ln(1 + sqrtxz1 - z)^2*z^-1*[1+z]^-1*sqrtxz1
       + 2*ln(1 + sqrtxz1 - z)^2*z^-1*sqrtxz1
       + 2*ln(1 + sqrtxz1 - z)^2*[1+z]^-1*sqrtxz1^-1
       - 4*ln(1 + sqrtxz1 - z)^2*sqrtxz1^-1
       + 2*ln(1 + sqrtxz1 - z)^2*z*[1+z]^-1*sqrtxz1^-1
       - 2*ln(1 + sqrtxz1 - z)^2*z*sqrtxz1^-1
       + 8*ln(1 + sqrtxz1 - z)^2*x*[1+z]^-1*sqrtxz1^-1
       - 8*ln(1 + sqrtxz1 - z)^2*x*sqrtxz1^-1
       - 4*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)
       + 4*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^-1*[1+z]^-1
       + 4*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^-1*[1+z]^-1*sqrtxz1
       - 4*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^-1
       - 4*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^-1*sqrtxz1
       - 4*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*[1+z]^-1*sqrtxz1^-1
       + 8*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*sqrtxz1^-1
       + 3*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*[1-x]^-1
       - 4*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z*[1+z]^-1*sqrtxz1^-1
       + 4*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z*sqrtxz1^-1
       - 3*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z
       - 4*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*z^-1*[1+z]^-1
       - 4*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*z^-1*[1+z]^-1*sqrtxz1
       + 4*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*z^-1
       + 4*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*z^-1*sqrtxz1
       + 4*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*[1+z]^-1*sqrtxz1^-1
       - 8*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*sqrtxz1^-1
       + 4*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*z*[1+z]^-1*sqrtxz1^-1
       - 4*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*z*sqrtxz1^-1
       + 4*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*z
       + 16*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*x*[1+z]^-1*sqrtxz1^-1
       - 16*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*x*sqrtxz1^-1
       - 4*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^-1*[1+z]^-1*rln2
       - 4*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^-1*[1+z]^-1*sqrtxz1*rln2
       + 4*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^-1*rln2
       + 4*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^-1*sqrtxz1*rln2
       + 4*ln(1 + sqrtxz1 + z)*[1-x]^-1*[1+z]^-1*sqrtxz1^-1*rln2
       - 8*ln(1 + sqrtxz1 + z)*[1-x]^-1*sqrtxz1^-1*rln2
       - 3*ln(1 + sqrtxz1 + z)*[1-x]^-1*rln2
       + 4*ln(1 + sqrtxz1 + z)*[1-x]^-1*z*[1+z]^-1*sqrtxz1^-1*rln2
       - 4*ln(1 + sqrtxz1 + z)*[1-x]^-1*z*sqrtxz1^-1*rln2
       + 3*ln(1 + sqrtxz1 + z)*[1-x]^-1*z*rln2
       + 4*ln(1 + sqrtxz1 + z)*z^-1*[1+z]^-1*rln2
       + 4*ln(1 + sqrtxz1 + z)*z^-1*[1+z]^-1*sqrtxz1*rln2
       - 4*ln(1 + sqrtxz1 + z)*z^-1*rln2
       - 4*ln(1 + sqrtxz1 + z)*z^-1*sqrtxz1*rln2
       - 4*ln(1 + sqrtxz1 + z)*[1+z]^-1*sqrtxz1^-1*rln2
       + 8*ln(1 + sqrtxz1 + z)*sqrtxz1^-1*rln2
       + 4*ln(1 + sqrtxz1 + z)*rln2
       - 4*ln(1 + sqrtxz1 + z)*z*[1+z]^-1*sqrtxz1^-1*rln2
       + 4*ln(1 + sqrtxz1 + z)*z*sqrtxz1^-1*rln2
       - 4*ln(1 + sqrtxz1 + z)*z*rln2
       - 16*ln(1 + sqrtxz1 + z)*x*[1+z]^-1*sqrtxz1^-1*rln2
       + 16*ln(1 + sqrtxz1 + z)*x*sqrtxz1^-1*rln2
       + 1/2*ln(1 - 2*z + z^2 + 4*x*z)^2*[1-x]^-1*z^-1*[1+z]^-1*sqrtxz1
       - 1/2*ln(1 - 2*z + z^2 + 4*x*z)^2*[1-x]^-1*z^-1*sqrtxz1
       - 1/2*ln(1 - 2*z + z^2 + 4*x*z)^2*[1-x]^-1*[1+z]^-1*sqrtxz1^-1
       + ln(1 - 2*z + z^2 + 4*x*z)^2*[1-x]^-1*sqrtxz1^-1
       - 1/2*ln(1 - 2*z + z^2 + 4*x*z)^2*[1-x]^-1*z*[1+z]^-1*sqrtxz1^-1
       + 1/2*ln(1 - 2*z + z^2 + 4*x*z)^2*[1-x]^-1*z*sqrtxz1^-1
       - 1/2*ln(1 - 2*z + z^2 + 4*x*z)^2*z^-1*[1+z]^-1*sqrtxz1
       + 1/2*ln(1 - 2*z + z^2 + 4*x*z)^2*z^-1*sqrtxz1
       + 1/2*ln(1 - 2*z + z^2 + 4*x*z)^2*[1+z]^-1*sqrtxz1^-1
       - ln(1 - 2*z + z^2 + 4*x*z)^2*sqrtxz1^-1
       + 1/2*ln(1 - 2*z + z^2 + 4*x*z)^2*z*[1+z]^-1*sqrtxz1^-1
       - 1/2*ln(1 - 2*z + z^2 + 4*x*z)^2*z*sqrtxz1^-1
       + 2*ln(1 - 2*z + z^2 + 4*x*z)^2*x*[1+z]^-1*sqrtxz1^-1
       - 2*ln(1 - 2*z + z^2 + 4*x*z)^2*x*sqrtxz1^-1
       - 2*ln(x)
       - 4*ln(x)*x^-1*[1+x]^-1*z^-1*[1-z]^-1
       + 4*ln(x)*x^-1*[1+x]^-1*z^-1
       + 4*ln(x)*x^-1*[1+x]^-1
       + 4*ln(x)*x^-1*[1+x]^-1*z
       + 4*ln(x)*x^-1*z^-1*[1-z]^-1
       - 4*ln(x)*x^-1*z^-1
       - 4*ln(x)*x^-1
       - 4*ln(x)*x^-1*z
       + 2*ln(x)*[1-x]^-1*z^-1*[1+z]^-1*rln2
       + 4*ln(x)*[1-x]^-1*z^-1*[1+z]^-1*sqrtxz1
       + 3/4*ln(x)*[1-x]^-1*z^-1
       - 2*ln(x)*[1-x]^-1*z^-1*rln2
       - 4*ln(x)*[1-x]^-1*z^-1*sqrtxz1
       - 4*ln(x)*[1-x]^-1*[1+z]^-1*sqrtxz1^-1
       + 8*ln(x)*[1-x]^-1*sqrtxz1^-1
       - ln(x)*[1-x]^-1
       + 3/2*ln(x)*[1-x]^-1*rln2
       - 1/2*ln(x)*[1-x]^-1*sqrtxz1
       - 4*ln(x)*[1-x]^-1*z*[1+z]^-1*sqrtxz1^-1
       + 4*ln(x)*[1-x]^-1*z*sqrtxz1^-1
       + 1/4*ln(x)*[1-x]^-1*z
       - 3/2*ln(x)*[1-x]^-1*z*rln2
       - 4*ln(x)*[1+x]^-1*[1-z]^-1
       + 4*ln(x)*[1+x]^-1
       + 4*ln(x)*[1+x]^-1*z
       - 4*ln(x)*z^-1*[1-z]^-1
       - 2*ln(x)*z^-1*[1+z]^-1*rln2
       - 4*ln(x)*z^-1*[1+z]^-1*sqrtxz1
       + 9/2*ln(x)*z^-1
       + 2*ln(x)*z^-1*rln2
       + 4*ln(x)*z^-1*sqrtxz1
       + 4*ln(x)*[1-z]^-1
       + 4*ln(x)*[1+z]^-1*sqrtxz1^-1
       - 8*ln(x)*sqrtxz1^-1
       - 2*ln(x)*rln2
       + 4*ln(x)*z*[1+z]^-1*sqrtxz1^-1
       - 4*ln(x)*z*sqrtxz1^-1
       + 2*ln(x)*z*rln2
       + 1/2*ln(x)*x*z^-1
       + 16*ln(x)*x*[1+z]^-1*sqrtxz1^-1
       - 16*ln(x)*x*sqrtxz1^-1
       - 2*ln(x)*x
       + 2*ln(x)*ln(1 + sqrtxz1 - z)*[1-x]^-1*z^-1*[1+z]^-1*sqrtxz1
       - 2*ln(x)*ln(1 + sqrtxz1 - z)*[1-x]^-1*z^-1*sqrtxz1
       - 2*ln(x)*ln(1 + sqrtxz1 - z)*[1-x]^-1*[1+z]^-1*sqrtxz1^-1
       + 4*ln(x)*ln(1 + sqrtxz1 - z)*[1-x]^-1*sqrtxz1^-1
       - 2*ln(x)*ln(1 + sqrtxz1 - z)*[1-x]^-1*z*[1+z]^-1*sqrtxz1^-1
       + 2*ln(x)*ln(1 + sqrtxz1 - z)*[1-x]^-1*z*sqrtxz1^-1
       - 2*ln(x)*ln(1 + sqrtxz1 - z)*z^-1*[1+z]^-1*sqrtxz1
       + 2*ln(x)*ln(1 + sqrtxz1 - z)*z^-1*sqrtxz1
       + 2*ln(x)*ln(1 + sqrtxz1 - z)*[1+z]^-1*sqrtxz1^-1
       - 4*ln(x)*ln(1 + sqrtxz1 - z)*sqrtxz1^-1
       + 2*ln(x)*ln(1 + sqrtxz1 - z)*z*[1+z]^-1*sqrtxz1^-1
       - 2*ln(x)*ln(1 + sqrtxz1 - z)*z*sqrtxz1^-1
       + 8*ln(x)*ln(1 + sqrtxz1 - z)*x*[1+z]^-1*sqrtxz1^-1
       - 8*ln(x)*ln(1 + sqrtxz1 - z)*x*sqrtxz1^-1
       + 2*ln(x)*ln(1 + sqrtxz1 + z)
       - 2*ln(x)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^-1*[1+z]^-1
       - 2*ln(x)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^-1*[1+z]^-1*sqrtxz1
       + 2*ln(x)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^-1
       + 2*ln(x)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^-1*sqrtxz1
       + 2*ln(x)*ln(1 + sqrtxz1 + z)*[1-x]^-1*[1+z]^-1*sqrtxz1^-1
       - 4*ln(x)*ln(1 + sqrtxz1 + z)*[1-x]^-1*sqrtxz1^-1
       - 3/2*ln(x)*ln(1 + sqrtxz1 + z)*[1-x]^-1
       + 2*ln(x)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z*[1+z]^-1*sqrtxz1^-1
       - 2*ln(x)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z*sqrtxz1^-1
       + 3/2*ln(x)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z
       + 2*ln(x)*ln(1 + sqrtxz1 + z)*z^-1*[1+z]^-1
       + 2*ln(x)*ln(1 + sqrtxz1 + z)*z^-1*[1+z]^-1*sqrtxz1
       - 2*ln(x)*ln(1 + sqrtxz1 + z)*z^-1
       - 2*ln(x)*ln(1 + sqrtxz1 + z)*z^-1*sqrtxz1
       - 2*ln(x)*ln(1 + sqrtxz1 + z)*[1+z]^-1*sqrtxz1^-1
       + 4*ln(x)*ln(1 + sqrtxz1 + z)*sqrtxz1^-1
       - 2*ln(x)*ln(1 + sqrtxz1 + z)*z*[1+z]^-1*sqrtxz1^-1
       + 2*ln(x)*ln(1 + sqrtxz1 + z)*z*sqrtxz1^-1
       - 2*ln(x)*ln(1 + sqrtxz1 + z)*z
       - 8*ln(x)*ln(1 + sqrtxz1 + z)*x*[1+z]^-1*sqrtxz1^-1
       + 8*ln(x)*ln(1 + sqrtxz1 + z)*x*sqrtxz1^-1
       - ln(x)*ln(1 + x*z^-1)
       + ln(x)*ln(1 + x*z^-1)*x^-1*[1+x]^-1*z^-1*[1-z]^-1
       - ln(x)*ln(1 + x*z^-1)*x^-1*[1+x]^-1*z^-1
       - ln(x)*ln(1 + x*z^-1)*x^-1*[1+x]^-1
       - ln(x)*ln(1 + x*z^-1)*x^-1*[1+x]^-1*z
       - ln(x)*ln(1 + x*z^-1)*x^-1*z^-1*[1-z]^-1
       + ln(x)*ln(1 + x*z^-1)*x^-1*z^-1
       + ln(x)*ln(1 + x*z^-1)*x^-1
       + ln(x)*ln(1 + x*z^-1)*x^-1*z
       + ln(x)*ln(1 + x*z^-1)*z^-1*[1-z]^-1
       - ln(x)*ln(1 + x*z^-1)*z^-1
       - ln(x)*ln(1 + x*z^-1)*z
       + 2*ln(x)*ln(1 + x)*[1+x]^-1
       + 2*ln(x)*ln(1 + x)*[1+x]^-1*z
       - 2*ln(x)*ln(1 + x)*z
       + 2*ln(x)*ln(1 + x)*x
       - 2*ln(x)*ln(1 + x*z)
       + ln(x)*ln(1 + x*z)*x^-1*[1+x]^-1*z^-1*[1-z]^-1
       - ln(x)*ln(1 + x*z)*x^-1*[1+x]^-1*z^-1
       - ln(x)*ln(1 + x*z)*x^-1*[1+x]^-1
       - ln(x)*ln(1 + x*z)*x^-1*[1+x]^-1*z
       - ln(x)*ln(1 + x*z)*x^-1*z^-1*[1-z]^-1
       + ln(x)*ln(1 + x*z)*x^-1*z^-1
       + ln(x)*ln(1 + x*z)*x^-1
       + ln(x)*ln(1 + x*z)*x^-1*z
       + ln(x)*ln(1 + x*z)*[1-x]^-1*z^-1*[1+z]^-1
       - ln(x)*ln(1 + x*z)*[1-x]^-1*z^-1
       + ln(x)*ln(1 + x*z)*[1-x]^-1
       - ln(x)*ln(1 + x*z)*[1-x]^-1*z
       + ln(x)*ln(1 + x*z)*z^-1*[1-z]^-1
       - ln(x)*ln(1 + x*z)*z^-1*[1+z]^-1
       + ln(x)*ln(z + x)
       - ln(x)*ln(z + x)*[1-x]^-1*z^-1*[1+z]^-1
       + ln(x)*ln(z + x)*[1-x]^-1*z^-1
       - ln(x)*ln(z + x)*[1-x]^-1
       + ln(x)*ln(z + x)*[1-x]^-1*z
       + ln(x)*ln(z + x)*z^-1*[1+z]^-1
       - ln(x)*ln(z + x)*z^-1
       - ln(x)*ln(z + x)*z
       + 1/2*ln(x)^2
       + 3/2*ln(x)^2*x^-1*[1+x]^-1*z^-1*[1-z]^-1
       - 3/2*ln(x)^2*x^-1*[1+x]^-1*z^-1
       - 3/2*ln(x)^2*x^-1*[1+x]^-1
       - 3/2*ln(x)^2*x^-1*[1+x]^-1*z
       - 3/2*ln(x)^2*x^-1*z^-1*[1-z]^-1
       + 3/2*ln(x)^2*x^-1*z^-1
       + 3/2*ln(x)^2*x^-1
       + 3/2*ln(x)^2*x^-1*z
       - 3/2*ln(x)^2*[1-x]^-1*z^-1*[1+z]^-1*sqrtxz1
       + 1/2*ln(x)^2*[1-x]^-1*z^-1
       + 3/2*ln(x)^2*[1-x]^-1*z^-1*sqrtxz1
       + 3/2*ln(x)^2*[1-x]^-1*[1+z]^-1*sqrtxz1^-1
       - 3*ln(x)^2*[1-x]^-1*sqrtxz1^-1
       - ln(x)^2*[1-x]^-1
       + 3/2*ln(x)^2*[1-x]^-1*z*[1+z]^-1*sqrtxz1^-1
       - 3/2*ln(x)^2*[1-x]^-1*z*sqrtxz1^-1
       + 1/2*ln(x)^2*[1-x]^-1*z
       + 3/2*ln(x)^2*[1+x]^-1*[1-z]^-1
       - 2*ln(x)^2*[1+x]^-1
       - 2*ln(x)^2*[1+x]^-1*z
       + 3/2*ln(x)^2*z^-1*[1-z]^-1
       + 3/2*ln(x)^2*z^-1*[1+z]^-1*sqrtxz1
       - 7/4*ln(x)^2*z^-1
       - 3/2*ln(x)^2*z^-1*sqrtxz1
       - 3/2*ln(x)^2*[1-z]^-1
       - 3/2*ln(x)^2*[1+z]^-1*sqrtxz1^-1
       + 3*ln(x)^2*sqrtxz1^-1
       - 3/2*ln(x)^2*z*[1+z]^-1*sqrtxz1^-1
       + 3/2*ln(x)^2*z*sqrtxz1^-1
       - 1/4*ln(x)^2*x*z^-1
       - 6*ln(x)^2*x*[1+z]^-1*sqrtxz1^-1
       + 6*ln(x)^2*x*sqrtxz1^-1
       - ln(x)*ln([1-x])
       + 2*ln(x)*ln([1-x])*x^-1*[1+x]^-1*z^-1*[1-z]^-1
       - 2*ln(x)*ln([1-x])*x^-1*[1+x]^-1*z^-1
       - 2*ln(x)*ln([1-x])*x^-1*[1+x]^-1
       - 2*ln(x)*ln([1-x])*x^-1*[1+x]^-1*z
       - 2*ln(x)*ln([1-x])*x^-1*z^-1*[1-z]^-1
       + 2*ln(x)*ln([1-x])*x^-1*z^-1
       + 2*ln(x)*ln([1-x])*x^-1
       + 2*ln(x)*ln([1-x])*x^-1*z
       - ln(x)*ln([1-x])*[1-x]^-1*z^-1
       + 2*ln(x)*ln([1-x])*[1-x]^-1
       - ln(x)*ln([1-x])*[1-x]^-1*z
       + 2*ln(x)*ln([1-x])*[1+x]^-1*[1-z]^-1
       - 2*ln(x)*ln([1-x])*[1+x]^-1
       - 2*ln(x)*ln([1-x])*[1+x]^-1*z
       + 2*ln(x)*ln([1-x])*z^-1*[1-z]^-1
       - 3/2*ln(x)*ln([1-x])*z^-1
       - 2*ln(x)*ln([1-x])*[1-z]^-1
       + ln(x)*ln([1-x])*z
       + 1/2*ln(x)*ln([1-x])*x*z^-1
       - ln(x)*ln([1-x])*x
       + 2*ln(x)*ln([1+x])
       - 2*ln(x)*ln([1+x])*x^-1*[1+x]^-1*z^-1*[1-z]^-1
       + 2*ln(x)*ln([1+x])*x^-1*[1+x]^-1*z^-1
       + 2*ln(x)*ln([1+x])*x^-1*[1+x]^-1
       + 2*ln(x)*ln([1+x])*x^-1*[1+x]^-1*z
       + 2*ln(x)*ln([1+x])*x^-1*z^-1*[1-z]^-1
       - 2*ln(x)*ln([1+x])*x^-1*z^-1
       - 2*ln(x)*ln([1+x])*x^-1
       - 2*ln(x)*ln([1+x])*x^-1*z
       - 2*ln(x)*ln([1+x])*z^-1*[1-z]^-1
       + 2*ln(x)*ln([1+x])*z^-1
       + 2*ln(x)*ln([1+x])*z
       - ln(x)*ln(z)
       - ln(x)*ln(z)*x^-1*[1+x]^-1*z^-1*[1-z]^-1
       + ln(x)*ln(z)*x^-1*[1+x]^-1*z^-1
       + ln(x)*ln(z)*x^-1*[1+x]^-1
       + ln(x)*ln(z)*x^-1*[1+x]^-1*z
       + ln(x)*ln(z)*x^-1*z^-1*[1-z]^-1
       - ln(x)*ln(z)*x^-1*z^-1
       - ln(x)*ln(z)*x^-1
       - ln(x)*ln(z)*x^-1*z
       + ln(x)*ln(z)*[1-x]^-1*z^-1*[1+z]^-1
       - ln(x)*ln(z)*[1-x]^-1*z^-1*[1+z]^-1*sqrtxz1
       - ln(x)*ln(z)*[1-x]^-1*z^-1
       + ln(x)*ln(z)*[1-x]^-1*z^-1*sqrtxz1
       + ln(x)*ln(z)*[1-x]^-1*[1+z]^-1*sqrtxz1^-1
       - 2*ln(x)*ln(z)*[1-x]^-1*sqrtxz1^-1
       + ln(x)*ln(z)*[1-x]^-1
       + ln(x)*ln(z)*[1-x]^-1*z*[1+z]^-1*sqrtxz1^-1
       - ln(x)*ln(z)*[1-x]^-1*z*sqrtxz1^-1
       - ln(x)*ln(z)*[1-x]^-1*z
       - ln(x)*ln(z)*z^-1*[1-z]^-1
       - ln(x)*ln(z)*z^-1*[1+z]^-1
       + ln(x)*ln(z)*z^-1*[1+z]^-1*sqrtxz1
       + 2*ln(x)*ln(z)*z^-1
       - ln(x)*ln(z)*z^-1*sqrtxz1
       + 1/2*ln(x)*ln(z)*[1-z]^-1
       - ln(x)*ln(z)*[1+z]^-1*sqrtxz1^-1
       + 2*ln(x)*ln(z)*sqrtxz1^-1
       - ln(x)*ln(z)*z*[1+z]^-1*sqrtxz1^-1
       + ln(x)*ln(z)*z*sqrtxz1^-1
       + 2*ln(x)*ln(z)*z
       + 1/2*ln(x)*ln(z)*x*[1-z]^-1
       - 4*ln(x)*ln(z)*x*[1+z]^-1*sqrtxz1^-1
       + 4*ln(x)*ln(z)*x*sqrtxz1^-1
       - ln(x)*ln(z)*x
       + 2*ln(x)*ln([1-z])*[1-x]^-1*z^-1*[1+z]^-1*sqrtxz1
       - 2*ln(x)*ln([1-z])*[1-x]^-1*z^-1*sqrtxz1
       - 2*ln(x)*ln([1-z])*[1-x]^-1*[1+z]^-1*sqrtxz1^-1
       + 4*ln(x)*ln([1-z])*[1-x]^-1*sqrtxz1^-1
       - 2*ln(x)*ln([1-z])*[1-x]^-1*z*[1+z]^-1*sqrtxz1^-1
       + 2*ln(x)*ln([1-z])*[1-x]^-1*z*sqrtxz1^-1
       - 2*ln(x)*ln([1-z])*z^-1*[1+z]^-1*sqrtxz1
       + 2*ln(x)*ln([1-z])*z^-1*sqrtxz1
       + 2*ln(x)*ln([1-z])*[1+z]^-1*sqrtxz1^-1
       - 4*ln(x)*ln([1-z])*sqrtxz1^-1
       + 2*ln(x)*ln([1-z])*z*[1+z]^-1*sqrtxz1^-1
       - 2*ln(x)*ln([1-z])*z*sqrtxz1^-1
       + 8*ln(x)*ln([1-z])*x*[1+z]^-1*sqrtxz1^-1
       - 8*ln(x)*ln([1-z])*x*sqrtxz1^-1
       - 3/2*ln(z)
       + 6*ln(z)*[1-x]^-1*z^-1*[1+z]^-1*rln2
       + 4*ln(z)*[1-x]^-1*z^-1*[1+z]^-1*sqrtxz1
       + 8*ln(z)*[1-x]^-1*z^-1*[1+z]^-1*sqrtxz1*rln2
       - 6*ln(z)*[1-x]^-1*z^-1*rln2
       - 4*ln(z)*[1-x]^-1*z^-1*sqrtxz1
       - 8*ln(z)*[1-x]^-1*z^-1*sqrtxz1*rln2
       - 4*ln(z)*[1-x]^-1*[1+z]^-1*sqrtxz1^-1
       - 8*ln(z)*[1-x]^-1*[1+z]^-1*sqrtxz1^-1*rln2
       + 8*ln(z)*[1-x]^-1*sqrtxz1^-1
       + 16*ln(z)*[1-x]^-1*sqrtxz1^-1*rln2
       + 1/2*ln(z)*[1-x]^-1
       + 9/2*ln(z)*[1-x]^-1*rln2
       - 1/2*ln(z)*[1-x]^-1*sqrtxz1
       - 4*ln(z)*[1-x]^-1*z*[1+z]^-1*sqrtxz1^-1
       - 8*ln(z)*[1-x]^-1*z*[1+z]^-1*sqrtxz1^-1*rln2
       + 4*ln(z)*[1-x]^-1*z*sqrtxz1^-1
       + 8*ln(z)*[1-x]^-1*z*sqrtxz1^-1*rln2
       + 1/2*ln(z)*[1-x]^-1*z
       - 9/2*ln(z)*[1-x]^-1*z*rln2
       - 6*ln(z)*z^-1*[1+z]^-1*rln2
       - 4*ln(z)*z^-1*[1+z]^-1*sqrtxz1
       - 8*ln(z)*z^-1*[1+z]^-1*sqrtxz1*rln2
       + 6*ln(z)*z^-1*rln2
       + 4*ln(z)*z^-1*sqrtxz1
       + 8*ln(z)*z^-1*sqrtxz1*rln2
       + 1/2*ln(z)*[1-z]^-1
       + 4*ln(z)*[1+z]^-1*sqrtxz1^-1
       + 8*ln(z)*[1+z]^-1*sqrtxz1^-1*rln2
       - 8*ln(z)*sqrtxz1^-1
       - 16*ln(z)*sqrtxz1^-1*rln2
       - 6*ln(z)*rln2
       + 4*ln(z)*z*[1+z]^-1*sqrtxz1^-1
       + 8*ln(z)*z*[1+z]^-1*sqrtxz1^-1*rln2
       - 4*ln(z)*z*sqrtxz1^-1
       - 8*ln(z)*z*sqrtxz1^-1*rln2
       + 6*ln(z)*z*rln2
       - 1/2*ln(z)*x*[1-z]^-1
       + 16*ln(z)*x*[1+z]^-1*sqrtxz1^-1
       + 32*ln(z)*x*[1+z]^-1*sqrtxz1^-1*rln2
       - 16*ln(z)*x*sqrtxz1^-1
       - 32*ln(z)*x*sqrtxz1^-1*rln2
       + 3/2*ln(z)*x
       + 2*ln(z)*ln(1 + sqrtxz1 - z)
       - 2*ln(z)*ln(1 + sqrtxz1 - z)*[1-x]^-1*z^-1*[1+z]^-1
       - 4*ln(z)*ln(1 + sqrtxz1 - z)*[1-x]^-1*z^-1*[1+z]^-1*sqrtxz1
       + 2*ln(z)*ln(1 + sqrtxz1 - z)*[1-x]^-1*z^-1
       + 4*ln(z)*ln(1 + sqrtxz1 - z)*[1-x]^-1*z^-1*sqrtxz1
       + 4*ln(z)*ln(1 + sqrtxz1 - z)*[1-x]^-1*[1+z]^-1*sqrtxz1^-1
       - 8*ln(z)*ln(1 + sqrtxz1 - z)*[1-x]^-1*sqrtxz1^-1
       - 3/2*ln(z)*ln(1 + sqrtxz1 - z)*[1-x]^-1
       + 4*ln(z)*ln(1 + sqrtxz1 - z)*[1-x]^-1*z*[1+z]^-1*sqrtxz1^-1
       - 4*ln(z)*ln(1 + sqrtxz1 - z)*[1-x]^-1*z*sqrtxz1^-1
       + 3/2*ln(z)*ln(1 + sqrtxz1 - z)*[1-x]^-1*z
       + 2*ln(z)*ln(1 + sqrtxz1 - z)*z^-1*[1+z]^-1
       + 4*ln(z)*ln(1 + sqrtxz1 - z)*z^-1*[1+z]^-1*sqrtxz1
       - 2*ln(z)*ln(1 + sqrtxz1 - z)*z^-1
       - 4*ln(z)*ln(1 + sqrtxz1 - z)*z^-1*sqrtxz1
       - 4*ln(z)*ln(1 + sqrtxz1 - z)*[1+z]^-1*sqrtxz1^-1
       + 8*ln(z)*ln(1 + sqrtxz1 - z)*sqrtxz1^-1
       - 4*ln(z)*ln(1 + sqrtxz1 - z)*z*[1+z]^-1*sqrtxz1^-1
       + 4*ln(z)*ln(1 + sqrtxz1 - z)*z*sqrtxz1^-1
       - 2*ln(z)*ln(1 + sqrtxz1 - z)*z
       - 16*ln(z)*ln(1 + sqrtxz1 - z)*x*[1+z]^-1*sqrtxz1^-1
       + 16*ln(z)*ln(1 + sqrtxz1 - z)*x*sqrtxz1^-1
       + 4*ln(z)*ln(1 + sqrtxz1 + z)
       - 4*ln(z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^-1*[1+z]^-1
       - 4*ln(z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^-1*[1+z]^-1*sqrtxz1
       + 4*ln(z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^-1
       + 4*ln(z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^-1*sqrtxz1
       + 4*ln(z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*[1+z]^-1*sqrtxz1^-1
       - 8*ln(z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*sqrtxz1^-1
       - 3*ln(z)*ln(1 + sqrtxz1 + z)*[1-x]^-1
       + 4*ln(z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z*[1+z]^-1*sqrtxz1^-1
       - 4*ln(z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z*sqrtxz1^-1
       + 3*ln(z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z
       + 4*ln(z)*ln(1 + sqrtxz1 + z)*z^-1*[1+z]^-1
       + 4*ln(z)*ln(1 + sqrtxz1 + z)*z^-1*[1+z]^-1*sqrtxz1
       - 4*ln(z)*ln(1 + sqrtxz1 + z)*z^-1
       - 4*ln(z)*ln(1 + sqrtxz1 + z)*z^-1*sqrtxz1
       - 4*ln(z)*ln(1 + sqrtxz1 + z)*[1+z]^-1*sqrtxz1^-1
       + 8*ln(z)*ln(1 + sqrtxz1 + z)*sqrtxz1^-1
       - 4*ln(z)*ln(1 + sqrtxz1 + z)*z*[1+z]^-1*sqrtxz1^-1
       + 4*ln(z)*ln(1 + sqrtxz1 + z)*z*sqrtxz1^-1
       - 4*ln(z)*ln(1 + sqrtxz1 + z)*z
       - 16*ln(z)*ln(1 + sqrtxz1 + z)*x*[1+z]^-1*sqrtxz1^-1
       + 16*ln(z)*ln(1 + sqrtxz1 + z)*x*sqrtxz1^-1
       - 2*ln(z)*ln(1 + z)*[1-x]^-1*[1+z]^-1
       + ln(z)*ln(1 + z)*[1-x]^-1
       - ln(z)*ln(1 + z)*[1-x]^-1*z
       + ln(z)*ln(1 + x*z^-1)
       - ln(z)*ln(1 + x*z^-1)*x^-1*[1+x]^-1*z^-1*[1-z]^-1
       + ln(z)*ln(1 + x*z^-1)*x^-1*[1+x]^-1*z^-1
       + ln(z)*ln(1 + x*z^-1)*x^-1*[1+x]^-1
       + ln(z)*ln(1 + x*z^-1)*x^-1*[1+x]^-1*z
       + ln(z)*ln(1 + x*z^-1)*x^-1*z^-1*[1-z]^-1
       - ln(z)*ln(1 + x*z^-1)*x^-1*z^-1
       - ln(z)*ln(1 + x*z^-1)*x^-1
       - ln(z)*ln(1 + x*z^-1)*x^-1*z
       - ln(z)*ln(1 + x*z^-1)*z^-1*[1-z]^-1
       + ln(z)*ln(1 + x*z^-1)*z^-1
       + ln(z)*ln(1 + x*z^-1)*z
       - 2*ln(z)*ln(1 + x*z)
       + ln(z)*ln(1 + x*z)*x^-1*[1+x]^-1*z^-1*[1-z]^-1
       - ln(z)*ln(1 + x*z)*x^-1*[1+x]^-1*z^-1
       - ln(z)*ln(1 + x*z)*x^-1*[1+x]^-1
       - ln(z)*ln(1 + x*z)*x^-1*[1+x]^-1*z
       - ln(z)*ln(1 + x*z)*x^-1*z^-1*[1-z]^-1
       + ln(z)*ln(1 + x*z)*x^-1*z^-1
       + ln(z)*ln(1 + x*z)*x^-1
       + ln(z)*ln(1 + x*z)*x^-1*z
       + ln(z)*ln(1 + x*z)*[1-x]^-1*z^-1*[1+z]^-1
       - ln(z)*ln(1 + x*z)*[1-x]^-1*z^-1
       + ln(z)*ln(1 + x*z)*[1-x]^-1
       - ln(z)*ln(1 + x*z)*[1-x]^-1*z
       + ln(z)*ln(1 + x*z)*z^-1*[1-z]^-1
       - ln(z)*ln(1 + x*z)*z^-1*[1+z]^-1
       - ln(z)*ln(z + x)
       + ln(z)*ln(z + x)*[1-x]^-1*z^-1*[1+z]^-1
       - ln(z)*ln(z + x)*[1-x]^-1*z^-1
       + ln(z)*ln(z + x)*[1-x]^-1
       - ln(z)*ln(z + x)*[1-x]^-1*z
       - ln(z)*ln(z + x)*z^-1*[1+z]^-1
       + ln(z)*ln(z + x)*z^-1
       + ln(z)*ln(z + x)*z
       - 1/2*ln(z)^2
       - 1/2*ln(z)^2*x^-1*[1+x]^-1*z^-1*[1-z]^-1
       + 1/2*ln(z)^2*x^-1*[1+x]^-1*z^-1
       + 1/2*ln(z)^2*x^-1*[1+x]^-1
       + 1/2*ln(z)^2*x^-1*[1+x]^-1*z
       + 1/2*ln(z)^2*x^-1*z^-1*[1-z]^-1
       - 1/2*ln(z)^2*x^-1*z^-1
       - 1/2*ln(z)^2*x^-1
       - 1/2*ln(z)^2*x^-1*z
       + ln(z)^2*[1-x]^-1*z^-1*[1+z]^-1
       + 5/2*ln(z)^2*[1-x]^-1*z^-1*[1+z]^-1*sqrtxz1
       - ln(z)^2*[1-x]^-1*z^-1
       - 5/2*ln(z)^2*[1-x]^-1*z^-1*sqrtxz1
       - 5/2*ln(z)^2*[1-x]^-1*[1+z]^-1*sqrtxz1^-1
       + 1/2*ln(z)^2*[1-x]^-1*[1+z]^-1
       + 5*ln(z)^2*[1-x]^-1*sqrtxz1^-1
       + 1/4*ln(z)^2*[1-x]^-1
       - 5/2*ln(z)^2*[1-x]^-1*z*[1+z]^-1*sqrtxz1^-1
       + 5/2*ln(z)^2*[1-x]^-1*z*sqrtxz1^-1
       - 1/4*ln(z)^2*[1-x]^-1*z
       - 1/2*ln(z)^2*z^-1*[1-z]^-1
       - ln(z)^2*z^-1*[1+z]^-1
       - 5/2*ln(z)^2*z^-1*[1+z]^-1*sqrtxz1
       + 3/2*ln(z)^2*z^-1
       + 5/2*ln(z)^2*z^-1*sqrtxz1
       + 5/2*ln(z)^2*[1+z]^-1*sqrtxz1^-1
       - 5*ln(z)^2*sqrtxz1^-1
       + 5/2*ln(z)^2*z*[1+z]^-1*sqrtxz1^-1
       - 5/2*ln(z)^2*z*sqrtxz1^-1
       + 3/2*ln(z)^2*z
       + 10*ln(z)^2*x*[1+z]^-1*sqrtxz1^-1
       - 10*ln(z)^2*x*sqrtxz1^-1
       + 2*ln(z)*ln([1-z])*[1-x]^-1*z^-1*[1+z]^-1*sqrtxz1
       - 2*ln(z)*ln([1-z])*[1-x]^-1*z^-1*sqrtxz1
       - 2*ln(z)*ln([1-z])*[1-x]^-1*[1+z]^-1*sqrtxz1^-1
       + 4*ln(z)*ln([1-z])*[1-x]^-1*sqrtxz1^-1
       - 2*ln(z)*ln([1-z])*[1-x]^-1*z*[1+z]^-1*sqrtxz1^-1
       + 2*ln(z)*ln([1-z])*[1-x]^-1*z*sqrtxz1^-1
       - 2*ln(z)*ln([1-z])*z^-1*[1+z]^-1*sqrtxz1
       + 2*ln(z)*ln([1-z])*z^-1*sqrtxz1
       + 2*ln(z)*ln([1-z])*[1+z]^-1*sqrtxz1^-1
       - 4*ln(z)*ln([1-z])*sqrtxz1^-1
       + 2*ln(z)*ln([1-z])*z*[1+z]^-1*sqrtxz1^-1
       - 2*ln(z)*ln([1-z])*z*sqrtxz1^-1
       + 8*ln(z)*ln([1-z])*x*[1+z]^-1*sqrtxz1^-1
       - 8*ln(z)*ln([1-z])*x*sqrtxz1^-1
       + 4*ln([1-z])*[1-x]^-1*z^-1*[1+z]^-1*sqrtxz1*rln2
       - 4*ln([1-z])*[1-x]^-1*z^-1*sqrtxz1*rln2
       - 4*ln([1-z])*[1-x]^-1*[1+z]^-1*sqrtxz1^-1*rln2
       + 8*ln([1-z])*[1-x]^-1*sqrtxz1^-1*rln2
       - 4*ln([1-z])*[1-x]^-1*z*[1+z]^-1*sqrtxz1^-1*rln2
       + 4*ln([1-z])*[1-x]^-1*z*sqrtxz1^-1*rln2
       - 4*ln([1-z])*z^-1*[1+z]^-1*sqrtxz1*rln2
       + 4*ln([1-z])*z^-1*sqrtxz1*rln2
       + 4*ln([1-z])*[1+z]^-1*sqrtxz1^-1*rln2
       - 8*ln([1-z])*sqrtxz1^-1*rln2
       + 4*ln([1-z])*z*[1+z]^-1*sqrtxz1^-1*rln2
       - 4*ln([1-z])*z*sqrtxz1^-1*rln2
       + 16*ln([1-z])*x*[1+z]^-1*sqrtxz1^-1*rln2
       - 16*ln([1-z])*x*sqrtxz1^-1*rln2
       - 4*ln([1-z])*ln(1 + sqrtxz1 - z)*[1-x]^-1*z^-1*[1+z]^-1*sqrtxz1
       + 4*ln([1-z])*ln(1 + sqrtxz1 - z)*[1-x]^-1*z^-1*sqrtxz1
       + 4*ln([1-z])*ln(1 + sqrtxz1 - z)*[1-x]^-1*[1+z]^-1*sqrtxz1^-1
       - 8*ln([1-z])*ln(1 + sqrtxz1 - z)*[1-x]^-1*sqrtxz1^-1
       + 4*ln([1-z])*ln(1 + sqrtxz1 - z)*[1-x]^-1*z*[1+z]^-1*sqrtxz1^-1
       - 4*ln([1-z])*ln(1 + sqrtxz1 - z)*[1-x]^-1*z*sqrtxz1^-1
       + 4*ln([1-z])*ln(1 + sqrtxz1 - z)*z^-1*[1+z]^-1*sqrtxz1
       - 4*ln([1-z])*ln(1 + sqrtxz1 - z)*z^-1*sqrtxz1
       - 4*ln([1-z])*ln(1 + sqrtxz1 - z)*[1+z]^-1*sqrtxz1^-1
       + 8*ln([1-z])*ln(1 + sqrtxz1 - z)*sqrtxz1^-1
       - 4*ln([1-z])*ln(1 + sqrtxz1 - z)*z*[1+z]^-1*sqrtxz1^-1
       + 4*ln([1-z])*ln(1 + sqrtxz1 - z)*z*sqrtxz1^-1
       - 16*ln([1-z])*ln(1 + sqrtxz1 - z)*x*[1+z]^-1*sqrtxz1^-1
       + 16*ln([1-z])*ln(1 + sqrtxz1 - z)*x*sqrtxz1^-1
       - 2*ln([1-z])*ln(1 - 2*z + z^2 + 4*x*z)*[1-x]^-1*z^-1*[1+z]^-1*sqrtxz1
       + 2*ln([1-z])*ln(1 - 2*z + z^2 + 4*x*z)*[1-x]^-1*z^-1*sqrtxz1
       + 2*ln([1-z])*ln(1 - 2*z + z^2 + 4*x*z)*[1-x]^-1*[1+z]^-1*sqrtxz1^-1
       - 4*ln([1-z])*ln(1 - 2*z + z^2 + 4*x*z)*[1-x]^-1*sqrtxz1^-1
       + 2*ln([1-z])*ln(1 - 2*z + z^2 + 4*x*z)*[1-x]^-1*z*[1+z]^-1*sqrtxz1^-1
       - 2*ln([1-z])*ln(1 - 2*z + z^2 + 4*x*z)*[1-x]^-1*z*sqrtxz1^-1
       + 2*ln([1-z])*ln(1 - 2*z + z^2 + 4*x*z)*z^-1*[1+z]^-1*sqrtxz1
       - 2*ln([1-z])*ln(1 - 2*z + z^2 + 4*x*z)*z^-1*sqrtxz1
       - 2*ln([1-z])*ln(1 - 2*z + z^2 + 4*x*z)*[1+z]^-1*sqrtxz1^-1
       + 4*ln([1-z])*ln(1 - 2*z + z^2 + 4*x*z)*sqrtxz1^-1
       - 2*ln([1-z])*ln(1 - 2*z + z^2 + 4*x*z)*z*[1+z]^-1*sqrtxz1^-1
       + 2*ln([1-z])*ln(1 - 2*z + z^2 + 4*x*z)*z*sqrtxz1^-1
       - 8*ln([1-z])*ln(1 - 2*z + z^2 + 4*x*z)*x*[1+z]^-1*sqrtxz1^-1
       + 8*ln([1-z])*ln(1 - 2*z + z^2 + 4*x*z)*x*sqrtxz1^-1
       + 2*ln([1-z])^2*[1-x]^-1*z^-1*[1+z]^-1*sqrtxz1
       - 2*ln([1-z])^2*[1-x]^-1*z^-1*sqrtxz1
       - 2*ln([1-z])^2*[1-x]^-1*[1+z]^-1*sqrtxz1^-1
       + 4*ln([1-z])^2*[1-x]^-1*sqrtxz1^-1
       - 2*ln([1-z])^2*[1-x]^-1*z*[1+z]^-1*sqrtxz1^-1
       + 2*ln([1-z])^2*[1-x]^-1*z*sqrtxz1^-1
       - 2*ln([1-z])^2*z^-1*[1+z]^-1*sqrtxz1
       + 2*ln([1-z])^2*z^-1*sqrtxz1
       + 2*ln([1-z])^2*[1+z]^-1*sqrtxz1^-1
       - 4*ln([1-z])^2*sqrtxz1^-1
       + 2*ln([1-z])^2*z*[1+z]^-1*sqrtxz1^-1
       - 2*ln([1-z])^2*z*sqrtxz1^-1
       + 8*ln([1-z])^2*x*[1+z]^-1*sqrtxz1^-1
       - 8*ln([1-z])^2*x*sqrtxz1^-1
       - 2*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)
       + 2*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*z^-1*[1+z]^-1
       + 2*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*z^-1*[1+z]^-1*sqrtxz1
       - 2*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*z^-1
       - 2*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*z^-1*sqrtxz1
       - 2*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*[1+z]^-1*sqrtxz1^-1
       + 4*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*sqrtxz1^-1
       + 3/2*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1
       - 2*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*z*[1+z]^-1*sqrtxz1^-1
       + 2*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*z*sqrtxz1^-1
       - 3/2*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*z
       - 2*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z^-1*[1+z]^-1
       - 2*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z^-1*[1+z]^-1*sqrtxz1
       + 2*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z^-1
       + 2*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z^-1*sqrtxz1
       + 2*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1+z]^-1*sqrtxz1^-1
       - 4*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*sqrtxz1^-1
       + 2*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z*[1+z]^-1*sqrtxz1^-1
       - 2*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z*sqrtxz1^-1
       + 2*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z
       + 8*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x*[1+z]^-1*sqrtxz1^-1
       - 8*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x*sqrtxz1^-1
       + 2*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)
       - 2*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*z^-1*[1+z]^-1
       + 2*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*z^-1*[1+z]^-1*sqrtxz1
       + 2*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*z^-1
       - 2*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*z^-1*sqrtxz1
       - 2*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*[1+z]^-1*sqrtxz1^-1
       + 4*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*sqrtxz1^-1
       - 3/2*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1
       - 2*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*z*[1+z]^-1*sqrtxz1^-1
       + 2*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*z*sqrtxz1^-1
       + 3/2*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*z
       + 2*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z^-1*[1+z]^-1
       - 2*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z^-1*[1+z]^-1*sqrtxz1
       - 2*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z^-1
       + 2*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z^-1*sqrtxz1
       + 2*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1+z]^-1*sqrtxz1^-1
       - 4*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*sqrtxz1^-1
       + 2*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z*[1+z]^-1*sqrtxz1^-1
       - 2*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z*sqrtxz1^-1
       - 2*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z
       + 8*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x*[1+z]^-1*sqrtxz1^-1
       - 8*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x*sqrtxz1^-1
       + 2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)
       - 2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*[1-x]^-1*z^-1*[1+z]^-1
       - 2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*[1-x]^-1*z^-1*[1+z]^-1*sqrtxz1
       + 2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*[1-x]^-1*z^-1
       + 2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*[1-x]^-1*z^-1*sqrtxz1
       + 2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*[1-x]^-1*[1+z]^-1*sqrtxz1^-1
       - 4*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*[1-x]^-1*sqrtxz1^-1
       - 3/2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*[1-x]^-1
       + 2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*[1-x]^-1*z*[1+z]^-1*sqrtxz1^-1
       - 2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*[1-x]^-1*z*sqrtxz1^-1
       + 3/2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*[1-x]^-1*z
       + 2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*z^-1*[1+z]^-1
       + 2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*z^-1*[1+z]^-1*sqrtxz1
       - 2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*z^-1
       - 2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*z^-1*sqrtxz1
       - 2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*[1+z]^-1*sqrtxz1^-1
       + 4*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*sqrtxz1^-1
       - 2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*z*[1+z]^-1*sqrtxz1^-1
       + 2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*z*sqrtxz1^-1
       - 2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*z
       - 8*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*x*[1+z]^-1*sqrtxz1^-1
       + 8*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*x*sqrtxz1^-1
       - 2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)
       + 2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*[1-x]^-1*z^-1*[1+z]^-1
       - 2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*[1-x]^-1*z^-1*[1+z]^-1*sqrtxz1
       - 2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*[1-x]^-1*z^-1
       + 2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*[1-x]^-1*z^-1*sqrtxz1
       + 2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*[1-x]^-1*[1+z]^-1*sqrtxz1^-1
       - 4*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*[1-x]^-1*sqrtxz1^-1
       + 3/2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*[1-x]^-1
       + 2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*[1-x]^-1*z*[1+z]^-1*sqrtxz1^-1
       - 2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*[1-x]^-1*z*sqrtxz1^-1
       - 3/2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*[1-x]^-1*z
       - 2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*z^-1*[1+z]^-1
       + 2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*z^-1*[1+z]^-1*sqrtxz1
       + 2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*z^-1
       - 2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*z^-1*sqrtxz1
       - 2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*[1+z]^-1*sqrtxz1^-1
       + 4*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*sqrtxz1^-1
       - 2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*z*[1+z]^-1*sqrtxz1^-1
       + 2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*z*sqrtxz1^-1
       + 2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*z
       - 8*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*x*[1+z]^-1*sqrtxz1^-1
       + 8*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*x*sqrtxz1^-1
       + 4*Li2(sqrtxz1*[z-1]^-1)*[1-x]^-1*z^-1*[1+z]^-1*sqrtxz1
       - 4*Li2(sqrtxz1*[z-1]^-1)*[1-x]^-1*z^-1*sqrtxz1
       - 4*Li2(sqrtxz1*[z-1]^-1)*[1-x]^-1*[1+z]^-1*sqrtxz1^-1
       + 8*Li2(sqrtxz1*[z-1]^-1)*[1-x]^-1*sqrtxz1^-1
       - 4*Li2(sqrtxz1*[z-1]^-1)*[1-x]^-1*z*[1+z]^-1*sqrtxz1^-1
       + 4*Li2(sqrtxz1*[z-1]^-1)*[1-x]^-1*z*sqrtxz1^-1
       - 4*Li2(sqrtxz1*[z-1]^-1)*z^-1*[1+z]^-1*sqrtxz1
       + 4*Li2(sqrtxz1*[z-1]^-1)*z^-1*sqrtxz1
       + 4*Li2(sqrtxz1*[z-1]^-1)*[1+z]^-1*sqrtxz1^-1
       - 8*Li2(sqrtxz1*[z-1]^-1)*sqrtxz1^-1
       + 4*Li2(sqrtxz1*[z-1]^-1)*z*[1+z]^-1*sqrtxz1^-1
       - 4*Li2(sqrtxz1*[z-1]^-1)*z*sqrtxz1^-1
       + 16*Li2(sqrtxz1*[z-1]^-1)*x*[1+z]^-1*sqrtxz1^-1
       - 16*Li2(sqrtxz1*[z-1]^-1)*x*sqrtxz1^-1
       + 4*Li2([1-z]*sqrtxz1^-1)*[1-x]^-1*z^-1*[1+z]^-1*sqrtxz1
       - 4*Li2([1-z]*sqrtxz1^-1)*[1-x]^-1*z^-1*sqrtxz1
       - 4*Li2([1-z]*sqrtxz1^-1)*[1-x]^-1*[1+z]^-1*sqrtxz1^-1
       + 8*Li2([1-z]*sqrtxz1^-1)*[1-x]^-1*sqrtxz1^-1
       - 4*Li2([1-z]*sqrtxz1^-1)*[1-x]^-1*z*[1+z]^-1*sqrtxz1^-1
       + 4*Li2([1-z]*sqrtxz1^-1)*[1-x]^-1*z*sqrtxz1^-1
       - 4*Li2([1-z]*sqrtxz1^-1)*z^-1*[1+z]^-1*sqrtxz1
       + 4*Li2([1-z]*sqrtxz1^-1)*z^-1*sqrtxz1
       + 4*Li2([1-z]*sqrtxz1^-1)*[1+z]^-1*sqrtxz1^-1
       - 8*Li2([1-z]*sqrtxz1^-1)*sqrtxz1^-1
       + 4*Li2([1-z]*sqrtxz1^-1)*z*[1+z]^-1*sqrtxz1^-1
       - 4*Li2([1-z]*sqrtxz1^-1)*z*sqrtxz1^-1
       + 16*Li2([1-z]*sqrtxz1^-1)*x*[1+z]^-1*sqrtxz1^-1
       - 16*Li2([1-z]*sqrtxz1^-1)*x*sqrtxz1^-1
       - 2*Li2( - z)*[1-x]^-1*[1+z]^-1
       + Li2( - z)*[1-x]^-1
       - Li2( - z)*[1-x]^-1*z
       + Li2( - x*z^-1)*x^-1*[1+x]^-1*z^-1*[1-z]^-1
       - Li2( - x*z^-1)*x^-1*[1+x]^-1*z^-1
       - Li2( - x*z^-1)*x^-1*[1+x]^-1
       - Li2( - x*z^-1)*x^-1*[1+x]^-1*z
       - Li2( - x*z^-1)*x^-1*z^-1*[1-z]^-1
       + Li2( - x*z^-1)*x^-1*z^-1
       + Li2( - x*z^-1)*x^-1
       + Li2( - x*z^-1)*x^-1*z
       - Li2( - x*z^-1)*[1-x]^-1*z^-1*[1+z]^-1
       + Li2( - x*z^-1)*[1-x]^-1*z^-1
       - Li2( - x*z^-1)*[1-x]^-1
       + Li2( - x*z^-1)*[1-x]^-1*z
       + Li2( - x*z^-1)*z^-1*[1-z]^-1
       + Li2( - x*z^-1)*z^-1*[1+z]^-1
       - 2*Li2( - x*z^-1)*z^-1
       - 2*Li2( - x*z^-1)*z
       + 2*Li2( - x)
       - 2*Li2( - x)*x^-1*[1+x]^-1*z^-1*[1-z]^-1
       + 2*Li2( - x)*x^-1*[1+x]^-1*z^-1
       + 2*Li2( - x)*x^-1*[1+x]^-1
       + 2*Li2( - x)*x^-1*[1+x]^-1*z
       + 2*Li2( - x)*x^-1*z^-1*[1-z]^-1
       - 2*Li2( - x)*x^-1*z^-1
       - 2*Li2( - x)*x^-1
       - 2*Li2( - x)*x^-1*z
       + 2*Li2( - x)*[1+x]^-1
       + 2*Li2( - x)*[1+x]^-1*z
       - 2*Li2( - x)*z^-1*[1-z]^-1
       + 2*Li2( - x)*z^-1
       + 2*Li2( - x)*x
       - 2*Li2( - x*z)
       + Li2( - x*z)*x^-1*[1+x]^-1*z^-1*[1-z]^-1
       - Li2( - x*z)*x^-1*[1+x]^-1*z^-1
       - Li2( - x*z)*x^-1*[1+x]^-1
       - Li2( - x*z)*x^-1*[1+x]^-1*z
       - Li2( - x*z)*x^-1*z^-1*[1-z]^-1
       + Li2( - x*z)*x^-1*z^-1
       + Li2( - x*z)*x^-1
       + Li2( - x*z)*x^-1*z
       + Li2( - x*z)*[1-x]^-1*z^-1*[1+z]^-1
       - Li2( - x*z)*[1-x]^-1*z^-1
       + Li2( - x*z)*[1-x]^-1
       - Li2( - x*z)*[1-x]^-1*z
       + Li2( - x*z)*z^-1*[1-z]^-1
       - Li2( - x*z)*z^-1*[1+z]^-1
       + 2*Li2( - 1/(1 + sqrtxz1 - z) + 1/(1 + sqrtxz1 - z)*sqrtxz1 + 1/(1 + sqrtxz1 - z)*z)*[1-x]^-1*z^-1*[1+z]^-1*sqrtxz1
       - 2*Li2( - 1/(1 + sqrtxz1 - z) + 1/(1 + sqrtxz1 - z)*sqrtxz1 + 1/(1 + sqrtxz1 - z)*z)*[1-x]^-1*z^-1*sqrtxz1
       - 2*Li2( - 1/(1 + sqrtxz1 - z) + 1/(1 + sqrtxz1 - z)*sqrtxz1 + 1/(1 + sqrtxz1 - z)*z)*[1-x]^-1*[1+z]^-1*sqrtxz1^-1
       + 4*Li2( - 1/(1 + sqrtxz1 - z) + 1/(1 + sqrtxz1 - z)*sqrtxz1 + 1/(1 + sqrtxz1 - z)*z)*[1-x]^-1*sqrtxz1^-1
       - 2*Li2( - 1/(1 + sqrtxz1 - z) + 1/(1 + sqrtxz1 - z)*sqrtxz1 + 1/(1 + sqrtxz1 - z)*z)*[1-x]^-1*z*[1+z]^-1*sqrtxz1^-1
       + 2*Li2( - 1/(1 + sqrtxz1 - z) + 1/(1 + sqrtxz1 - z)*sqrtxz1 + 1/(1 + sqrtxz1 - z)*z)*[1-x]^-1*z*sqrtxz1^-1
       - 2*Li2( - 1/(1 + sqrtxz1 - z) + 1/(1 + sqrtxz1 - z)*sqrtxz1 + 1/(1 + sqrtxz1 - z)*z)*z^-1*[1+z]^-1*sqrtxz1
       + 2*Li2( - 1/(1 + sqrtxz1 - z) + 1/(1 + sqrtxz1 - z)*sqrtxz1 + 1/(1 + sqrtxz1 - z)*z)*z^-1*sqrtxz1
       + 2*Li2( - 1/(1 + sqrtxz1 - z) + 1/(1 + sqrtxz1 - z)*sqrtxz1 + 1/(1 + sqrtxz1 - z)*z)*[1+z]^-1*sqrtxz1^-1
       - 4*Li2( - 1/(1 + sqrtxz1 - z) + 1/(1 + sqrtxz1 - z)*sqrtxz1 + 1/(1 + sqrtxz1 - z)*z)*sqrtxz1^-1
       + 2*Li2( - 1/(1 + sqrtxz1 - z) + 1/(1 + sqrtxz1 - z)*sqrtxz1 + 1/(1 + sqrtxz1 - z)*z)*z*[1+z]^-1*sqrtxz1^-1
       - 2*Li2( - 1/(1 + sqrtxz1 - z) + 1/(1 + sqrtxz1 - z)*sqrtxz1 + 1/(1 + sqrtxz1 - z)*z)*z*sqrtxz1^-1
       + 8*Li2( - 1/(1 + sqrtxz1 - z) + 1/(1 + sqrtxz1 - z)*sqrtxz1 + 1/(1 + sqrtxz1 - z)*z)*x*[1+z]^-1*sqrtxz1^-1
       - 8*Li2( - 1/(1 + sqrtxz1 - z) + 1/(1 + sqrtxz1 - z)*sqrtxz1 + 1/(1 + sqrtxz1 - z)*z)*x*sqrtxz1^-1
       - Li2(x)
       + 2*Li2(x)*x^-1*[1+x]^-1*z^-1*[1-z]^-1
       - 2*Li2(x)*x^-1*[1+x]^-1*z^-1
       - 2*Li2(x)*x^-1*[1+x]^-1
       - 2*Li2(x)*x^-1*[1+x]^-1*z
       - 2*Li2(x)*x^-1*z^-1*[1-z]^-1
       + 2*Li2(x)*x^-1*z^-1
       + 2*Li2(x)*x^-1
       + 2*Li2(x)*x^-1*z
       - Li2(x)*[1-x]^-1*z^-1
       + 2*Li2(x)*[1-x]^-1
       - Li2(x)*[1-x]^-1*z
       + 2*Li2(x)*[1+x]^-1*[1-z]^-1
       - 2*Li2(x)*[1+x]^-1
       - 2*Li2(x)*[1+x]^-1*z
       + 2*Li2(x)*z^-1*[1-z]^-1
       - 3/2*Li2(x)*z^-1
       - 2*Li2(x)*[1-z]^-1
       + Li2(x)*z
       + 1/2*Li2(x)*x*z^-1
       - Li2(x)*x );

id DC2Q2QP1    = (  + NC^-1 * (  - 1/3
       + 5/72*z^-1
       - 37/24*z
       + 65/36*z^2
       - 19/72*x*z^-1
       - 5/6*x
       + 47/24*x*z
       - 31/36*x*z^2
       + 1/12*pi^2
       + 1/6*pi^2*z
       + 1/12*pi^2*x
       - 3/2*ln(x)
       + 2/3*ln(x)*[1-x]^-1*z^-1
       + 5/4*ln(x)*[1-x]^-1
       - 5/4*ln(x)*[1-x]^-1*z
       - 2/3*ln(x)*[1-x]^-1*z^2
       - 1/3*ln(x)*z^-1
       + 1/2*ln(x)*z
       + 4/3*ln(x)*z^2
       - 1/3*ln(x)*x*z^-1
       - 1/2*ln(x)*x
       + 3/2*ln(x)*x*z
       - 2/3*ln(x)*x*z^2
       - ln(x)*ln(z)
       + 3/2*ln(x)*ln(z)*[1-x]^-1
       + 3/2*ln(x)*ln(z)*[1-x]^-1*z
       - 2*ln(x)*ln(z)*z
       - ln(x)*ln(z)*x
       + 3/4*ln([1-x])
       + 1/6*ln([1-x])*z^-1
       - 1/4*ln([1-x])*z
       - 2/3*ln([1-x])*z^2
       + 1/6*ln([1-x])*x*z^-1
       + 1/4*ln([1-x])*x
       - 3/4*ln([1-x])*x*z
       + 1/3*ln([1-x])*x*z^2
       + 1/2*ln([1-x])*ln(z)
       + ln([1-x])*ln(z)*z
       + 1/2*ln([1-x])*ln(z)*x
       + 5/4*ln(z)
       + 1/6*ln(z)*z^-1
       - 3/4*ln(z)*z
       - 2/3*ln(z)*z^2
       + 1/6*ln(z)*x*z^-1
       - 1/4*ln(z)*x
       - 3/4*ln(z)*x*z
       + 1/3*ln(z)*x*z^2
       + 1/2*ln(z)^2
       + ln(z)^2*z
       + 1/2*ln(z)^2*x
       + 3/4*ln([1-z])
       + 1/6*ln([1-z])*z^-1
       - 1/4*ln([1-z])*z
       - 2/3*ln([1-z])*z^2
       + 1/6*ln([1-z])*x*z^-1
       + 1/4*ln([1-z])*x
       - 3/4*ln([1-z])*x*z
       + 1/3*ln([1-z])*x*z^2
       - 1/2*Li2(z)
       - Li2(z)*z
       - 1/2*Li2(z)*x
       )

       + NC * ( 1/3
       - 5/72*z^-1
       + 37/24*z
       - 65/36*z^2
       + 19/72*x*z^-1
       + 5/6*x
       - 47/24*x*z
       + 31/36*x*z^2
       - 1/12*pi^2
       - 1/6*pi^2*z
       - 1/12*pi^2*x
       + 3/2*ln(x)
       - 2/3*ln(x)*[1-x]^-1*z^-1
       - 5/4*ln(x)*[1-x]^-1
       + 5/4*ln(x)*[1-x]^-1*z
       + 2/3*ln(x)*[1-x]^-1*z^2
       + 1/3*ln(x)*z^-1
       - 1/2*ln(x)*z
       - 4/3*ln(x)*z^2
       + 1/3*ln(x)*x*z^-1
       + 1/2*ln(x)*x
       - 3/2*ln(x)*x*z
       + 2/3*ln(x)*x*z^2
       + ln(x)*ln(z)
       - 3/2*ln(x)*ln(z)*[1-x]^-1
       - 3/2*ln(x)*ln(z)*[1-x]^-1*z
       + 2*ln(x)*ln(z)*z
       + ln(x)*ln(z)*x
       - 3/4*ln([1-x])
       - 1/6*ln([1-x])*z^-1
       + 1/4*ln([1-x])*z
       + 2/3*ln([1-x])*z^2
       - 1/6*ln([1-x])*x*z^-1
       - 1/4*ln([1-x])*x
       + 3/4*ln([1-x])*x*z
       - 1/3*ln([1-x])*x*z^2
       - 1/2*ln([1-x])*ln(z)
       - ln([1-x])*ln(z)*z
       - 1/2*ln([1-x])*ln(z)*x
       - 5/4*ln(z)
       - 1/6*ln(z)*z^-1
       + 3/4*ln(z)*z
       + 2/3*ln(z)*z^2
       - 1/6*ln(z)*x*z^-1
       + 1/4*ln(z)*x
       + 3/4*ln(z)*x*z
       - 1/3*ln(z)*x*z^2
       - 1/2*ln(z)^2
       - ln(z)^2*z
       - 1/2*ln(z)^2*x
       - 3/4*ln([1-z])
       - 1/6*ln([1-z])*z^-1
       + 1/4*ln([1-z])*z
       + 2/3*ln([1-z])*z^2
       - 1/6*ln([1-z])*x*z^-1
       - 1/4*ln([1-z])*x
       + 3/4*ln([1-z])*x*z
       - 1/3*ln([1-z])*x*z^2
       + 1/2*Li2(z)
       + Li2(z)*z
       + 1/2*Li2(z)*x
       )

       + LMUA*NC^-1 * (  - 3/4
       - 1/6*z^-1
       + 1/4*z
       + 2/3*z^2
       - 1/6*x*z^-1
       - 1/4*x
       + 3/4*x*z
       - 1/3*x*z^2
       - 1/2*ln(z)
       - ln(z)*z
       - 1/2*ln(z)*x
       )

       + LMUA*NC * ( 3/4
       + 1/6*z^-1
       - 1/4*z
       - 2/3*z^2
       + 1/6*x*z^-1
       + 1/4*x
       - 3/4*x*z
       + 1/3*x*z^2
       + 1/2*ln(z)
       + ln(z)*z
       + 1/2*ln(z)*x
       )

       + Dd([1-x])*NC^-1 * ( 11/18
       + 107/216*z^-1
       + 2/9*z
       - 287/216*z^2
       - zeta3
       - zeta3*z
       + 1/12*pi^2
       + 1/8*pi^2*z
       + 5/2*ln(z)
       + 7/36*ln(z)*z^-1
       + 9/4*ln(z)*z
       + 31/36*ln(z)*z^2
       - 1/12*ln(z)*pi^2
       - 1/12*ln(z)*pi^2*z
       + 1/16*ln(z)^2
       - 1/6*ln(z)^2*z^-1
       + 5/16*ln(z)^2*z
       - 5/24*ln(z)^3
       - 5/24*ln(z)^3*z
       + 3/4*ln(z)*ln([1-z])^2
       + 3/4*ln(z)*ln([1-z])^2*z
       + 5/3*ln([1-z])
       + 7/36*ln([1-z])*z^-1
       - 7/6*ln([1-z])*z
       - 25/36*ln([1-z])*z^2
       - 1/6*ln([1-z])*pi^2
       - 1/6*ln([1-z])*pi^2*z
       - 1/8*ln([1-z])^2
       - 1/6*ln([1-z])^2*z^-1
       + 1/8*ln([1-z])^2*z
       + 1/6*ln([1-z])^2*z^2
       + 1/2*ln([1-z])*Li2(1 - z)
       + 1/2*ln([1-z])*Li2(1 - z)*z
       + ln([1-z])*Li2(z)
       + ln([1-z])*Li2(z)*z
       + 1/2*Li3(1 - z)
       + 1/2*Li3(1 - z)*z
       + Li3(z)
       + Li3(z)*z
       - 1/4*Li2(z)
       + 1/3*Li2(z)*z^-1
       - Li2(z)*z
       - 1/3*Li2(z)*z^2
       )

       + Dd([1-x])*NC * (  - 11/18
       - 107/216*z^-1
       - 2/9*z
       + 287/216*z^2
       + zeta3
       + zeta3*z
       - 1/12*pi^2
       - 1/8*pi^2*z
       - 5/2*ln(z)
       - 7/36*ln(z)*z^-1
       - 9/4*ln(z)*z
       - 31/36*ln(z)*z^2
       + 1/12*ln(z)*pi^2
       + 1/12*ln(z)*pi^2*z
       - 1/16*ln(z)^2
       + 1/6*ln(z)^2*z^-1
       - 5/16*ln(z)^2*z
       + 5/24*ln(z)^3
       + 5/24*ln(z)^3*z
       - 3/4*ln(z)*ln([1-z])^2
       - 3/4*ln(z)*ln([1-z])^2*z
       - 5/3*ln([1-z])
       - 7/36*ln([1-z])*z^-1
       + 7/6*ln([1-z])*z
       + 25/36*ln([1-z])*z^2
       + 1/6*ln([1-z])*pi^2
       + 1/6*ln([1-z])*pi^2*z
       + 1/8*ln([1-z])^2
       + 1/6*ln([1-z])^2*z^-1
       - 1/8*ln([1-z])^2*z
       - 1/6*ln([1-z])^2*z^2
       - 1/2*ln([1-z])*Li2(1 - z)
       - 1/2*ln([1-z])*Li2(1 - z)*z
       - ln([1-z])*Li2(z)
       - ln([1-z])*Li2(z)*z
       - 1/2*Li3(1 - z)
       - 1/2*Li3(1 - z)*z
       - Li3(z)
       - Li3(z)*z
       + 1/4*Li2(z)
       - 1/3*Li2(z)*z^-1
       + Li2(z)*z
       + 1/3*Li2(z)*z^2
       )

       + Dd([1-x])*LMUA*NC^-1 * (  - 5/3
       - 7/36*z^-1
       + 7/6*z
       + 25/36*z^2
       + 1/12*pi^2
       + 1/12*pi^2*z
       - 1/4*ln(z)
       + 1/3*ln(z)*z^-1
       - ln(z)*z
       - 1/3*ln(z)*z^2
       + 1/2*ln(z)^2
       + 1/2*ln(z)^2*z
       + 1/4*ln([1-z])
       + 1/3*ln([1-z])*z^-1
       - 1/4*ln([1-z])*z
       - 1/3*ln([1-z])*z^2
       - 1/2*Li2(z)
       - 1/2*Li2(z)*z
       )

       + Dd([1-x])*LMUA*NC * ( 5/3
       + 7/36*z^-1
       - 7/6*z
       - 25/36*z^2
       - 1/12*pi^2
       - 1/12*pi^2*z
       + 1/4*ln(z)
       - 1/3*ln(z)*z^-1
       + ln(z)*z
       + 1/3*ln(z)*z^2
       - 1/2*ln(z)^2
       - 1/2*ln(z)^2*z
       - 1/4*ln([1-z])
       - 1/3*ln([1-z])*z^-1
       + 1/4*ln([1-z])*z
       + 1/3*ln([1-z])*z^2
       + 1/2*Li2(z)
       + 1/2*Li2(z)*z
       )

       + Dd([1-x])*LMUA^2*NC^-1 * (  - 1/8
       - 1/6*z^-1
       + 1/8*z
       + 1/6*z^2
       - 1/4*ln(z)
       - 1/4*ln(z)*z
       )

       + Dd([1-x])*LMUA^2*NC * ( 1/8
       + 1/6*z^-1
       - 1/8*z
       - 1/6*z^2
       + 1/4*ln(z)
       + 1/4*ln(z)*z
       )

       + Dn(0,[1-x])*NC^-1 * ( 5/3
       + 7/36*z^-1
       - 7/6*z
       - 25/36*z^2
       - 1/12*pi^2
       - 1/12*pi^2*z
       + 1/4*ln(z)
       - 1/3*ln(z)*z^-1
       + ln(z)*z
       + 1/3*ln(z)*z^2
       - 1/2*ln(z)^2
       - 1/2*ln(z)^2*z
       - 1/4*ln([1-z])
       - 1/3*ln([1-z])*z^-1
       + 1/4*ln([1-z])*z
       + 1/3*ln([1-z])*z^2
       + 1/2*Li2(z)
       + 1/2*Li2(z)*z
       )

       + Dn(0,[1-x])*NC * (  - 5/3
       - 7/36*z^-1
       + 7/6*z
       + 25/36*z^2
       + 1/12*pi^2
       + 1/12*pi^2*z
       - 1/4*ln(z)
       + 1/3*ln(z)*z^-1
       - ln(z)*z
       - 1/3*ln(z)*z^2
       + 1/2*ln(z)^2
       + 1/2*ln(z)^2*z
       + 1/4*ln([1-z])
       + 1/3*ln([1-z])*z^-1
       - 1/4*ln([1-z])*z
       - 1/3*ln([1-z])*z^2
       - 1/2*Li2(z)
       - 1/2*Li2(z)*z
       )

       + Dn(0,[1-x])*LMUA*NC^-1 * ( 1/4
       + 1/3*z^-1
       - 1/4*z
       - 1/3*z^2
       + 1/2*ln(z)
       + 1/2*ln(z)*z
       )

       + Dn(0,[1-x])*LMUA*NC * (  - 1/4
       - 1/3*z^-1
       + 1/4*z
       + 1/3*z^2
       - 1/2*ln(z)
       - 1/2*ln(z)*z
       )

       + Dn(1,[1-x])*NC^-1 * (  - 1/4
       - 1/3*z^-1
       + 1/4*z
       + 1/3*z^2
       - 1/2*ln(z)
       - 1/2*ln(z)*z
       )

       + Dn(1,[1-x])*NC * ( 1/4
       + 1/3*z^-1
       - 1/4*z
       - 1/3*z^2
       + 1/2*ln(z)
       + 1/2*ln(z)*z ) );

id DC2Q2QP2    = (  + NC^-1 * (  - 5
       + 3/2*z^-1
       - 3/2*x*z^-1
       + 5*x
       - 1/12*pi^2*z^-1
       + 1/6*pi^2
       - 1/12*pi^2*x*z^-1
       + 1/6*pi^2*x
       - 4*ln(x)
       + 3/2*ln(x)*z^-1
       - 1/4*ln(x)*poly2^-1
       - 3/2*ln(x)*x*z^-1
       + 1/4*ln(x)*x*poly2^-1
       + 2*ln(x)*x
       + 1/4*ln(x)*x^2*poly2^-1
       - 1/4*ln(x)*x^3*poly2^-1
       - 1/8*ln(x)*ln(1 - sqrtxz2 + x)*sqrtxz2^-1*poly2^-1
       - 7/8*ln(x)*ln(1 - sqrtxz2 + x)*sqrtxz2^-1
       - 1/4*ln(x)*ln(1 - sqrtxz2 + x)*x*sqrtxz2^-1
       + 1/2*ln(x)*ln(1 - sqrtxz2 + x)*x*z*sqrtxz2^-1
       + 1/4*ln(x)*ln(1 - sqrtxz2 + x)*x^2*sqrtxz2^-1*poly2^-1
       - 7/8*ln(x)*ln(1 - sqrtxz2 + x)*x^2*sqrtxz2^-1
       - 1/8*ln(x)*ln(1 - sqrtxz2 + x)*x^4*sqrtxz2^-1*poly2^-1
       + 1/8*ln(x)*ln(1 + sqrtxz2 + x)*sqrtxz2^-1*poly2^-1
       + 7/8*ln(x)*ln(1 + sqrtxz2 + x)*sqrtxz2^-1
       + 1/4*ln(x)*ln(1 + sqrtxz2 + x)*x*sqrtxz2^-1
       - 1/2*ln(x)*ln(1 + sqrtxz2 + x)*x*z*sqrtxz2^-1
       - 1/4*ln(x)*ln(1 + sqrtxz2 + x)*x^2*sqrtxz2^-1*poly2^-1
       + 7/8*ln(x)*ln(1 + sqrtxz2 + x)*x^2*sqrtxz2^-1
       + 1/8*ln(x)*ln(1 + sqrtxz2 + x)*x^4*sqrtxz2^-1*poly2^-1
       - ln(x)^2
       + 1/2*ln(x)^2*z^-1
       + 1/2*ln(x)^2*x*z^-1
       - ln(x)^2*x
       + 1/2*ln(x)*ln(z)
       - 1/2*ln(x)*ln(z)*z^-1
       - 1/2*ln(x)*ln(z)*x*z^-1
       + 1/2*ln(x)*ln(z)*x
       + ln(x)*ln([1-z])
       - 1/2*ln(x)*ln([1-z])*z^-1
       - 1/2*ln(x)*ln([1-z])*x*z^-1
       + ln(x)*ln([1-z])*x
       + 5/2*ln([1-x])
       - 5/4*ln([1-x])*z^-1
       + 5/4*ln([1-x])*x*z^-1
       - 5/2*ln([1-x])*x
       + ln(z)
       - 5/4*ln(z)*z^-1
       - 1/2*ln(z)*[1-z]^-1
       - 1/4*ln(z)*poly2^-1
       + 5/4*ln(z)*x*z^-1
       + 1/2*ln(z)*x*[1-z]^-1
       - 1/4*ln(z)*x*poly2^-1
       - ln(z)*x
       + 1/4*ln(z)*x^2*poly2^-1
       + 1/4*ln(z)*x^3*poly2^-1
       + 5/2*ln([1-z])
       - 5/4*ln([1-z])*z^-1
       + 5/4*ln([1-z])*x*z^-1
       - 5/2*ln([1-z])*x
       - 1/8*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*sqrtxz2^-1*poly2^-1
       - 7/8*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*sqrtxz2^-1
       - 1/4*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x*sqrtxz2^-1
       + 1/2*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x*z*sqrtxz2^-1
       + 1/4*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x^2*sqrtxz2^-1*poly2^-1
       - 7/8*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x^2*sqrtxz2^-1
       - 1/8*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x^4*sqrtxz2^-1*poly2^-1
       + 1/8*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*sqrtxz2^-1*poly2^-1
       + 7/8*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*sqrtxz2^-1
       + 1/4*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x*sqrtxz2^-1
       - 1/2*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x*z*sqrtxz2^-1
       - 1/4*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x^2*sqrtxz2^-1*poly2^-1
       + 7/8*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x^2*sqrtxz2^-1
       + 1/8*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x^4*sqrtxz2^-1*poly2^-1
       + 1/8*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*sqrtxz2^-1*poly2^-1
       + 7/8*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*sqrtxz2^-1
       + 1/4*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x*sqrtxz2^-1
       - 1/2*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x*z*sqrtxz2^-1
       - 1/4*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x^2*sqrtxz2^-1*poly2^-1
       + 7/8*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x^2*sqrtxz2^-1
       + 1/8*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x^4*sqrtxz2^-1*poly2^-1
       - 1/8*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*sqrtxz2^-1*poly2^-1
       - 7/8*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*sqrtxz2^-1
       - 1/4*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x*sqrtxz2^-1
       + 1/2*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x*z*sqrtxz2^-1
       + 1/4*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x^2*sqrtxz2^-1*poly2^-1
       - 7/8*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x^2*sqrtxz2^-1
       - 1/8*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x^4*sqrtxz2^-1*poly2^-1
       - Li2(x)
       + 1/2*Li2(x)*z^-1
       + 1/2*Li2(x)*x*z^-1
       - Li2(x)*x
       )

       + NC * ( 5
       - 3/2*z^-1
       + 3/2*x*z^-1
       - 5*x
       + 1/12*pi^2*z^-1
       - 1/6*pi^2
       + 1/12*pi^2*x*z^-1
       - 1/6*pi^2*x
       + 4*ln(x)
       - 3/2*ln(x)*z^-1
       + 1/4*ln(x)*poly2^-1
       + 3/2*ln(x)*x*z^-1
       - 1/4*ln(x)*x*poly2^-1
       - 2*ln(x)*x
       - 1/4*ln(x)*x^2*poly2^-1
       + 1/4*ln(x)*x^3*poly2^-1
       + 1/8*ln(x)*ln(1 - sqrtxz2 + x)*sqrtxz2^-1*poly2^-1
       + 7/8*ln(x)*ln(1 - sqrtxz2 + x)*sqrtxz2^-1
       + 1/4*ln(x)*ln(1 - sqrtxz2 + x)*x*sqrtxz2^-1
       - 1/2*ln(x)*ln(1 - sqrtxz2 + x)*x*z*sqrtxz2^-1
       - 1/4*ln(x)*ln(1 - sqrtxz2 + x)*x^2*sqrtxz2^-1*poly2^-1
       + 7/8*ln(x)*ln(1 - sqrtxz2 + x)*x^2*sqrtxz2^-1
       + 1/8*ln(x)*ln(1 - sqrtxz2 + x)*x^4*sqrtxz2^-1*poly2^-1
       - 1/8*ln(x)*ln(1 + sqrtxz2 + x)*sqrtxz2^-1*poly2^-1
       - 7/8*ln(x)*ln(1 + sqrtxz2 + x)*sqrtxz2^-1
       - 1/4*ln(x)*ln(1 + sqrtxz2 + x)*x*sqrtxz2^-1
       + 1/2*ln(x)*ln(1 + sqrtxz2 + x)*x*z*sqrtxz2^-1
       + 1/4*ln(x)*ln(1 + sqrtxz2 + x)*x^2*sqrtxz2^-1*poly2^-1
       - 7/8*ln(x)*ln(1 + sqrtxz2 + x)*x^2*sqrtxz2^-1
       - 1/8*ln(x)*ln(1 + sqrtxz2 + x)*x^4*sqrtxz2^-1*poly2^-1
       + ln(x)^2
       - 1/2*ln(x)^2*z^-1
       - 1/2*ln(x)^2*x*z^-1
       + ln(x)^2*x
       - 1/2*ln(x)*ln(z)
       + 1/2*ln(x)*ln(z)*z^-1
       + 1/2*ln(x)*ln(z)*x*z^-1
       - 1/2*ln(x)*ln(z)*x
       - ln(x)*ln([1-z])
       + 1/2*ln(x)*ln([1-z])*z^-1
       + 1/2*ln(x)*ln([1-z])*x*z^-1
       - ln(x)*ln([1-z])*x
       - 5/2*ln([1-x])
       + 5/4*ln([1-x])*z^-1
       - 5/4*ln([1-x])*x*z^-1
       + 5/2*ln([1-x])*x
       - ln(z)
       + 5/4*ln(z)*z^-1
       + 1/2*ln(z)*[1-z]^-1
       + 1/4*ln(z)*poly2^-1
       - 5/4*ln(z)*x*z^-1
       - 1/2*ln(z)*x*[1-z]^-1
       + 1/4*ln(z)*x*poly2^-1
       + ln(z)*x
       - 1/4*ln(z)*x^2*poly2^-1
       - 1/4*ln(z)*x^3*poly2^-1
       - 5/2*ln([1-z])
       + 5/4*ln([1-z])*z^-1
       - 5/4*ln([1-z])*x*z^-1
       + 5/2*ln([1-z])*x
       + 1/8*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*sqrtxz2^-1*poly2^-1
       + 7/8*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*sqrtxz2^-1
       + 1/4*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x*sqrtxz2^-1
       - 1/2*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x*z*sqrtxz2^-1
       - 1/4*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x^2*sqrtxz2^-1*poly2^-1
       + 7/8*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x^2*sqrtxz2^-1
       + 1/8*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x^4*sqrtxz2^-1*poly2^-1
       - 1/8*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*sqrtxz2^-1*poly2^-1
       - 7/8*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*sqrtxz2^-1
       - 1/4*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x*sqrtxz2^-1
       + 1/2*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x*z*sqrtxz2^-1
       + 1/4*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x^2*sqrtxz2^-1*poly2^-1
       - 7/8*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x^2*sqrtxz2^-1
       - 1/8*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x^4*sqrtxz2^-1*poly2^-1
       - 1/8*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*sqrtxz2^-1*poly2^-1
       - 7/8*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*sqrtxz2^-1
       - 1/4*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x*sqrtxz2^-1
       + 1/2*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x*z*sqrtxz2^-1
       + 1/4*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x^2*sqrtxz2^-1*poly2^-1
       - 7/8*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x^2*sqrtxz2^-1
       - 1/8*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x^4*sqrtxz2^-1*poly2^-1
       + 1/8*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*sqrtxz2^-1*poly2^-1
       + 7/8*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*sqrtxz2^-1
       + 1/4*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x*sqrtxz2^-1
       - 1/2*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x*z*sqrtxz2^-1
       - 1/4*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x^2*sqrtxz2^-1*poly2^-1
       + 7/8*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x^2*sqrtxz2^-1
       + 1/8*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x^4*sqrtxz2^-1*poly2^-1
       + Li2(x)
       - 1/2*Li2(x)*z^-1
       - 1/2*Li2(x)*x*z^-1
       + Li2(x)*x
       )

       + LMUF*NC^-1 * (  - 5/2
       + 5/4*z^-1
       - 5/4*x*z^-1
       + 5/2*x
       - ln(x)
       + 1/2*ln(x)*z^-1
       + 1/2*ln(x)*x*z^-1
       - ln(x)*x
       )

       + LMUF*NC * ( 5/2
       - 5/4*z^-1
       + 5/4*x*z^-1
       - 5/2*x
       + ln(x)
       - 1/2*ln(x)*z^-1
       - 1/2*ln(x)*x*z^-1
       + ln(x)*x
       )

       + Dd([1-z])*NC^-1 * (  - 9/2
       + 9/2*x
       + 1/4*pi^2
       - 1/4*pi^2*x
       - 25/8*ln(x)
       + 3/8*ln(x)*x
       + 1/6*ln(x)*pi^2
       + 1/6*ln(x)*pi^2*x
       - 17/16*ln(x)^2
       + 15/16*ln(x)^2*x
       - 5/24*ln(x)^3
       - 5/24*ln(x)^3*x
       + 5/4*ln(x)*ln([1-x])
       - 5/4*ln(x)*ln([1-x])*x
       + 3/4*ln(x)*ln([1-x])^2
       + 3/4*ln(x)*ln([1-x])^2*x
       - 1/2*ln(x)*Li2(x)
       - 1/2*ln(x)*Li2(x)*x
       + 3/2*ln([1-x])
       - 3/2*ln([1-x])*x
       - 1/6*ln([1-x])*pi^2
       - 1/6*ln([1-x])*pi^2*x
       - 5/8*ln([1-x])^2
       + 5/8*ln([1-x])^2*x
       + 1/2*ln([1-x])*Li2(1 - x)
       + 1/2*ln([1-x])*Li2(1 - x)*x
       + ln([1-x])*Li2(x)
       + ln([1-x])*Li2(x)*x
       + 1/2*Li3(1 - x)
       + 1/2*Li3(1 - x)*x
       - 1/4*Li2(x)
       + 1/4*Li2(x)*x
       )

       + Dd([1-z])*NC * ( 9/2
       - 9/2*x
       - 1/4*pi^2
       + 1/4*pi^2*x
       + 25/8*ln(x)
       - 3/8*ln(x)*x
       - 1/6*ln(x)*pi^2
       - 1/6*ln(x)*pi^2*x
       + 17/16*ln(x)^2
       - 15/16*ln(x)^2*x
       + 5/24*ln(x)^3
       + 5/24*ln(x)^3*x
       - 5/4*ln(x)*ln([1-x])
       + 5/4*ln(x)*ln([1-x])*x
       - 3/4*ln(x)*ln([1-x])^2
       - 3/4*ln(x)*ln([1-x])^2*x
       + 1/2*ln(x)*Li2(x)
       + 1/2*ln(x)*Li2(x)*x
       - 3/2*ln([1-x])
       + 3/2*ln([1-x])*x
       + 1/6*ln([1-x])*pi^2
       + 1/6*ln([1-x])*pi^2*x
       + 5/8*ln([1-x])^2
       - 5/8*ln([1-x])^2*x
       - 1/2*ln([1-x])*Li2(1 - x)
       - 1/2*ln([1-x])*Li2(1 - x)*x
       - ln([1-x])*Li2(x)
       - ln([1-x])*Li2(x)*x
       - 1/2*Li3(1 - x)
       - 1/2*Li3(1 - x)*x
       + 1/4*Li2(x)
       - 1/4*Li2(x)*x
       )

       + Dd([1-z])*LMUF*NC^-1 * (  - 3/2
       + 3/2*x
       + 1/12*pi^2
       + 1/12*pi^2*x
       - 3/2*ln(x)
       + 3/2*ln(x)*x
       - 1/2*ln(x)^2
       - 1/2*ln(x)^2*x
       + 5/4*ln([1-x])
       - 5/4*ln([1-x])*x
       - 1/2*Li2(x)
       - 1/2*Li2(x)*x
       )

       + Dd([1-z])*LMUF*NC * ( 3/2
       - 3/2*x
       - 1/12*pi^2
       - 1/12*pi^2*x
       + 3/2*ln(x)
       - 3/2*ln(x)*x
       + 1/2*ln(x)^2
       + 1/2*ln(x)^2*x
       - 5/4*ln([1-x])
       + 5/4*ln([1-x])*x
       + 1/2*Li2(x)
       + 1/2*Li2(x)*x
       )

       + Dd([1-z])*LMUF^2*NC^-1 * (  - 5/8
       + 5/8*x
       - 1/4*ln(x)
       - 1/4*ln(x)*x
       )

       + Dd([1-z])*LMUF^2*NC * ( 5/8
       - 5/8*x
       + 1/4*ln(x)
       + 1/4*ln(x)*x
       )

       + Dn(0,[1-z])*NC^-1 * ( 3/2
       - 3/2*x
       - 1/12*pi^2
       - 1/12*pi^2*x
       + 3/2*ln(x)
       - 3/2*ln(x)*x
       + 1/2*ln(x)^2
       + 1/2*ln(x)^2*x
       - 5/4*ln([1-x])
       + 5/4*ln([1-x])*x
       + 1/2*Li2(x)
       + 1/2*Li2(x)*x
       )

       + Dn(0,[1-z])*NC * (  - 3/2
       + 3/2*x
       + 1/12*pi^2
       + 1/12*pi^2*x
       - 3/2*ln(x)
       + 3/2*ln(x)*x
       - 1/2*ln(x)^2
       - 1/2*ln(x)^2*x
       + 5/4*ln([1-x])
       - 5/4*ln([1-x])*x
       - 1/2*Li2(x)
       - 1/2*Li2(x)*x
       )

       + Dn(0,[1-z])*LMUF*NC^-1 * ( 5/4
       - 5/4*x
       + 1/2*ln(x)
       + 1/2*ln(x)*x
       )

       + Dn(0,[1-z])*LMUF*NC * (  - 5/4
       + 5/4*x
       - 1/2*ln(x)
       - 1/2*ln(x)*x
       )

       + Dn(1,[1-z])*NC^-1 * (  - 5/4
       + 5/4*x
       - 1/2*ln(x)
       - 1/2*ln(x)*x
       )

       + Dn(1,[1-z])*NC * ( 5/4
       - 5/4*x
       + 1/2*ln(x)
       + 1/2*ln(x)*x ) );

id DC2Q2QP3    = (  + NC^-1 * (  - 5/2
       - 3*[1-x]^-1*rln2^2
       - 3*[1-x]^-1*sqrtxz1*rln2
       - 2*[1-x]^-1*z*rln2^2
       - 8*[1-x]^-1*z^2*rln2^2
       + 2*rln2^2
       + 2*sqrtxz1*rln2
       + 2*z
       + 8*z^2*rln2^2
       + 4*x*z*rln2^2
       + 1/6*pi^2*x^-2*[1+x]^-1
       - 1/3*pi^2*x^-2*[1+x]^-1*z
       - 1/6*pi^2*x^-2
       + 1/3*pi^2*x^-2*z
       - 2/3*pi^2*x^-1*[1+x]^-1*z^2
       + 1/6*pi^2*x^-1
       - 1/3*pi^2*x^-1*z
       + 2/3*pi^2*x^-1*z^2
       + 5/12*pi^2*[1-x]^-1
       - 5/6*pi^2*[1-x]^-1*z
       + pi^2*[1-x]^-1*z^2
       - 1/3*pi^2*[1+x]^-1
       + 2/3*pi^2*[1+x]^-1*z
       - pi^2*[1+x]^-1*z^2
       - 1/6*pi^2
       + 1/3*pi^2*z
       - 2/3*pi^2*z^2
       + 1/6*pi^2*x
       - 1/3*pi^2*x*z
       + 4*ln(1 + sqrtxz1 - z)*[1-x]^-1*rln2
       + 3*ln(1 + sqrtxz1 - z)*[1-x]^-1*sqrtxz1
       + 12*ln(1 + sqrtxz1 - z)*[1-x]^-1*z^2*rln2
       - 3*ln(1 + sqrtxz1 - z)*rln2
       - 2*ln(1 + sqrtxz1 - z)*sqrtxz1
       + 2*ln(1 + sqrtxz1 - z)*z*rln2
       - 12*ln(1 + sqrtxz1 - z)*z^2*rln2
       + ln(1 + sqrtxz1 - z)*x*rln2
       - 6*ln(1 + sqrtxz1 - z)*x*z*rln2
       + ln(1 + sqrtxz1 - z)^2
       - ln(1 + sqrtxz1 - z)^2*[1-x]^-1
       + 2*ln(1 + sqrtxz1 - z)^2*[1-x]^-1*z
       - 4*ln(1 + sqrtxz1 - z)^2*[1-x]^-1*z^2
       - 2*ln(1 + sqrtxz1 - z)^2*z
       + 4*ln(1 + sqrtxz1 - z)^2*z^2
       - ln(1 + sqrtxz1 - z)^2*x
       + 2*ln(1 + sqrtxz1 - z)^2*x*z
       + ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)
       - 2*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*[1-x]^-1
       - 4*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z
       - 4*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^2
       + 2*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*z
       + 4*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*z^2
       + ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*x
       + 2*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*x*z
       + 2*ln(1 + sqrtxz1 + z)*[1-x]^-1*rln2
       + 4*ln(1 + sqrtxz1 + z)*[1-x]^-1*z*rln2
       + 4*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^2*rln2
       - ln(1 + sqrtxz1 + z)*rln2
       - 2*ln(1 + sqrtxz1 + z)*z*rln2
       - 4*ln(1 + sqrtxz1 + z)*z^2*rln2
       - ln(1 + sqrtxz1 + z)*x*rln2
       - 2*ln(1 + sqrtxz1 + z)*x*z*rln2
       + 4*ln(z*sqrtxz3)*ArcTan(z*sqrtxz3)*z*sqrtxz3
       - 3/4*ln(x)
       + 2*ln(x)*x^-2*[1+x]^-1
       - 4*ln(x)*x^-2*[1+x]^-1*z
       - 2*ln(x)*x^-2
       + 4*ln(x)*x^-2*z
       - 8*ln(x)*x^-1*[1+x]^-1*z^2
       + 2*ln(x)*x^-1
       - 4*ln(x)*x^-1*z
       + 8*ln(x)*x^-1*z^2
       + ln(x)*[1-x]^-1
       - 5/2*ln(x)*[1-x]^-1*rln2
       - 3/2*ln(x)*[1-x]^-1*sqrtxz1
       - 3/2*ln(x)*[1-x]^-1*z
       + ln(x)*[1-x]^-1*z*rln2
       - 8*ln(x)*[1-x]^-1*z^2*rln2
       - 2*ln(x)*[1+x]^-1
       + 4*ln(x)*[1+x]^-1*z
       - 8*ln(x)*[1+x]^-1*z^2
       + 1/4*ln(x)*poly2^-1
       + 2*ln(x)*rln2
       + ln(x)*sqrtxz1
       + ln(x)*z
       - 2*ln(x)*z*rln2
       + 8*ln(x)*z^2*rln2
       + 1/2*ln(x)*x*[x-z]^-1
       - 1/4*ln(x)*x*poly2^-1
       + 3/4*ln(x)*x
       - ln(x)*x*rln2
       + 4*ln(x)*x*z*rln2
       - ln(x)*x^2*[x-z]^-1
       - 1/4*ln(x)*x^2*poly2^-1
       + 1/4*ln(x)*x^3*poly2^-1
       + 1/8*ln(x)*ln(1 - sqrtxz2 + x)*sqrtxz2^-1*poly2^-1
       - 1/8*ln(x)*ln(1 - sqrtxz2 + x)*sqrtxz2^-1
       + 5/4*ln(x)*ln(1 - sqrtxz2 + x)*x*sqrtxz2^-1
       - 5/2*ln(x)*ln(1 - sqrtxz2 + x)*x*z*sqrtxz2^-1
       - 1/4*ln(x)*ln(1 - sqrtxz2 + x)*x^2*sqrtxz2^-1*poly2^-1
       - 1/8*ln(x)*ln(1 - sqrtxz2 + x)*x^2*sqrtxz2^-1
       + 1/8*ln(x)*ln(1 - sqrtxz2 + x)*x^4*sqrtxz2^-1*poly2^-1
       - 1/8*ln(x)*ln(1 + sqrtxz2 + x)*sqrtxz2^-1*poly2^-1
       + 1/8*ln(x)*ln(1 + sqrtxz2 + x)*sqrtxz2^-1
       - 5/4*ln(x)*ln(1 + sqrtxz2 + x)*x*sqrtxz2^-1
       + 5/2*ln(x)*ln(1 + sqrtxz2 + x)*x*z*sqrtxz2^-1
       + 1/4*ln(x)*ln(1 + sqrtxz2 + x)*x^2*sqrtxz2^-1*poly2^-1
       + 1/8*ln(x)*ln(1 + sqrtxz2 + x)*x^2*sqrtxz2^-1
       - 1/8*ln(x)*ln(1 + sqrtxz2 + x)*x^4*sqrtxz2^-1*poly2^-1
       - ln(x)*ln(1 + sqrtxz1 - z)
       + ln(x)*ln(1 + sqrtxz1 - z)*[1-x]^-1
       - 2*ln(x)*ln(1 + sqrtxz1 - z)*[1-x]^-1*z
       + 4*ln(x)*ln(1 + sqrtxz1 - z)*[1-x]^-1*z^2
       + 2*ln(x)*ln(1 + sqrtxz1 - z)*z
       - 4*ln(x)*ln(1 + sqrtxz1 - z)*z^2
       + ln(x)*ln(1 + sqrtxz1 - z)*x
       - 2*ln(x)*ln(1 + sqrtxz1 - z)*x*z
       - ln(x)*ln(1 + sqrtxz1 + z)
       + 3/2*ln(x)*ln(1 + sqrtxz1 + z)*[1-x]^-1
       + ln(x)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z
       + 4*ln(x)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^2
       - 4*ln(x)*ln(1 + sqrtxz1 + z)*z^2
       - 2*ln(x)*ln(1 + sqrtxz1 + z)*x*z
       + 1/2*ln(x)*ln(1 + x*z^-1)
       - 1/2*ln(x)*ln(1 + x*z^-1)*x^-2*[1+x]^-1
       + ln(x)*ln(1 + x*z^-1)*x^-2*[1+x]^-1*z
       + 1/2*ln(x)*ln(1 + x*z^-1)*x^-2
       - ln(x)*ln(1 + x*z^-1)*x^-2*z
       + 2*ln(x)*ln(1 + x*z^-1)*x^-1*[1+x]^-1*z^2
       - 1/2*ln(x)*ln(1 + x*z^-1)*x^-1
       + ln(x)*ln(1 + x*z^-1)*x^-1*z
       - 2*ln(x)*ln(1 + x*z^-1)*x^-1*z^2
       - 1/2*ln(x)*ln(1 + x*z^-1)*[1+x]^-1
       + ln(x)*ln(1 + x*z^-1)*[1+x]^-1*z
       - ln(x)*ln(1 + x*z^-1)*z
       + 2*ln(x)*ln(1 + x*z^-1)*z^2
       - 1/2*ln(x)*ln(1 + x*z^-1)*x
       + ln(x)*ln(1 + x*z^-1)*x*z
       + 3*ln(x)*ln(1 + x)
       - 2*ln(x)*ln(1 + x)*[1+x]^-1
       + 4*ln(x)*ln(1 + x)*[1+x]^-1*z
       - 4*ln(x)*ln(1 + x)*[1+x]^-1*z^2
       - 6*ln(x)*ln(1 + x)*z
       + 4*ln(x)*ln(1 + x)*z^2
       + ln(x)*ln(1 + x)*x
       - 2*ln(x)*ln(1 + x)*x*z
       + ln(x)*ln(1 + x*z)
       - 1/2*ln(x)*ln(1 + x*z)*x^-2*[1+x]^-1
       + ln(x)*ln(1 + x*z)*x^-2*[1+x]^-1*z
       + 1/2*ln(x)*ln(1 + x*z)*x^-2
       - ln(x)*ln(1 + x*z)*x^-2*z
       + 2*ln(x)*ln(1 + x*z)*x^-1*[1+x]^-1*z^2
       - 1/2*ln(x)*ln(1 + x*z)*x^-1
       + ln(x)*ln(1 + x*z)*x^-1*z
       - 2*ln(x)*ln(1 + x*z)*x^-1*z^2
       - ln(x)*ln(1 + x*z)*[1-x]^-1
       - 2*ln(x)*ln(1 + x*z)*[1-x]^-1*z
       - 2*ln(x)*ln(1 + x*z)*[1-x]^-1*z^2
       - 1/2*ln(x)*ln(1 + x*z)*[1+x]^-1
       + ln(x)*ln(1 + x*z)*[1+x]^-1*z
       + 4*ln(x)*ln(1 + x*z)*z^2
       + 2*ln(x)*ln(1 + x*z)*x*z
       - 1/2*ln(x)*ln(z + x)
       + ln(x)*ln(z + x)*[1-x]^-1
       + 2*ln(x)*ln(z + x)*[1-x]^-1*z
       + 2*ln(x)*ln(z + x)*[1-x]^-1*z^2
       - ln(x)*ln(z + x)*z
       - 2*ln(x)*ln(z + x)*z^2
       - 1/2*ln(x)*ln(z + x)*x
       - ln(x)*ln(z + x)*x*z
       - 3/2*ln(x)^2
       - 3/4*ln(x)^2*x^-2*[1+x]^-1
       + 3/2*ln(x)^2*x^-2*[1+x]^-1*z
       + 3/4*ln(x)^2*x^-2
       - 3/2*ln(x)^2*x^-2*z
       + 3*ln(x)^2*x^-1*[1+x]^-1*z^2
       - 3/4*ln(x)^2*x^-1
       + 3/2*ln(x)^2*x^-1*z
       - 3*ln(x)^2*x^-1*z^2
       + ln(x)^2*[1-x]^-1
       - 2*ln(x)^2*[1-x]^-1*z
       + ln(x)^2*[1-x]^-1*z^2
       + 5/4*ln(x)^2*[1+x]^-1
       - 5/2*ln(x)^2*[1+x]^-1*z
       + 4*ln(x)^2*[1+x]^-1*z^2
       + 3*ln(x)^2*z
       - 2*ln(x)^2*z^2
       - ln(x)^2*x
       + 2*ln(x)^2*x*z
       + ln(x)*ln([1-x])
       - ln(x)*ln([1-x])*x^-2*[1+x]^-1
       + 2*ln(x)*ln([1-x])*x^-2*[1+x]^-1*z
       + ln(x)*ln([1-x])*x^-2
       - 2*ln(x)*ln([1-x])*x^-2*z
       + 4*ln(x)*ln([1-x])*x^-1*[1+x]^-1*z^2
       - ln(x)*ln([1-x])*x^-1
       + 2*ln(x)*ln([1-x])*x^-1*z
       - 4*ln(x)*ln([1-x])*x^-1*z^2
       - 3/2*ln(x)*ln([1-x])*[1-x]^-1
       + 3*ln(x)*ln([1-x])*[1-x]^-1*z
       - 4*ln(x)*ln([1-x])*[1-x]^-1*z^2
       + ln(x)*ln([1-x])*[1+x]^-1
       - 2*ln(x)*ln([1-x])*[1+x]^-1*z
       + 4*ln(x)*ln([1-x])*[1+x]^-1*z^2
       - 2*ln(x)*ln([1-x])*z
       + 4*ln(x)*ln([1-x])*z^2
       - ln(x)*ln([1+x])
       + ln(x)*ln([1+x])*x^-2*[1+x]^-1
       - 2*ln(x)*ln([1+x])*x^-2*[1+x]^-1*z
       - ln(x)*ln([1+x])*x^-2
       + 2*ln(x)*ln([1+x])*x^-2*z
       - 4*ln(x)*ln([1+x])*x^-1*[1+x]^-1*z^2
       + ln(x)*ln([1+x])*x^-1
       - 2*ln(x)*ln([1+x])*x^-1*z
       + 4*ln(x)*ln([1+x])*x^-1*z^2
       + ln(x)*ln([1+x])*[1+x]^-1
       - 2*ln(x)*ln([1+x])*[1+x]^-1*z
       + 2*ln(x)*ln([1+x])*z
       - 4*ln(x)*ln([1+x])*z^2
       + ln(x)*ln([1+x])*x
       - 2*ln(x)*ln([1+x])*x*z
       + 1/2*ln(x)*ln(z)
       + 1/2*ln(x)*ln(z)*x^-2*[1+x]^-1
       - ln(x)*ln(z)*x^-2*[1+x]^-1*z
       - 1/2*ln(x)*ln(z)*x^-2
       + ln(x)*ln(z)*x^-2*z
       - 2*ln(x)*ln(z)*x^-1*[1+x]^-1*z^2
       + 1/2*ln(x)*ln(z)*x^-1
       - ln(x)*ln(z)*x^-1*z
       + 2*ln(x)*ln(z)*x^-1*z^2
       - 2*ln(x)*ln(z)*[1-x]^-1
       + 2*ln(x)*ln(z)*[1-x]^-1*z
       - 4*ln(x)*ln(z)*[1-x]^-1*z^2
       + 1/2*ln(x)*ln(z)*[1+x]^-1
       - ln(x)*ln(z)*[1+x]^-1*z
       - ln(x)*ln(z)*z
       + 2*ln(x)*ln(z)*z^2
       - 1/2*ln(x)*ln(z)*x
       + ln(x)*ln(z)*x*z
       - 9/4*ln(z)
       + 3/2*ln(z)*[1-x]^-1
       - 7/2*ln(z)*[1-x]^-1*rln2
       - 3/2*ln(z)*[1-x]^-1*sqrtxz1
       + 3/2*ln(z)*[1-x]^-1*z
       - 5*ln(z)*[1-x]^-1*z*rln2
       - 8*ln(z)*[1-x]^-1*z^2*rln2
       + 1/4*ln(z)*poly2^-1
       + 2*ln(z)*rln2
       + ln(z)*sqrtxz1
       - ln(z)*z
       + 2*ln(z)*z*rln2
       + 8*ln(z)*z^2*rln2
       - 1/2*ln(z)*x*[x-z]^-1
       + 1/4*ln(z)*x*poly2^-1
       - 3/4*ln(z)*x
       + ln(z)*x*rln2
       + 4*ln(z)*x*z*rln2
       + ln(z)*x^2*[x-z]^-1
       - 1/4*ln(z)*x^2*poly2^-1
       - 1/4*ln(z)*x^3*poly2^-1
       - ln(z)*ln(1 + sqrtxz1 - z)
       + 3/2*ln(z)*ln(1 + sqrtxz1 - z)*[1-x]^-1
       + ln(z)*ln(1 + sqrtxz1 - z)*[1-x]^-1*z
       + 4*ln(z)*ln(1 + sqrtxz1 - z)*[1-x]^-1*z^2
       - 4*ln(z)*ln(1 + sqrtxz1 - z)*z^2
       - 2*ln(z)*ln(1 + sqrtxz1 - z)*x*z
       - ln(z)*ln(1 + sqrtxz1 + z)
       + 2*ln(z)*ln(1 + sqrtxz1 + z)*[1-x]^-1
       + 4*ln(z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z
       + 4*ln(z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^2
       - 2*ln(z)*ln(1 + sqrtxz1 + z)*z
       - 4*ln(z)*ln(1 + sqrtxz1 + z)*z^2
       - ln(z)*ln(1 + sqrtxz1 + z)*x
       - 2*ln(z)*ln(1 + sqrtxz1 + z)*x*z
       - 1/2*ln(z)*ln(1 + x*z^-1)
       + 1/2*ln(z)*ln(1 + x*z^-1)*x^-2*[1+x]^-1
       - ln(z)*ln(1 + x*z^-1)*x^-2*[1+x]^-1*z
       - 1/2*ln(z)*ln(1 + x*z^-1)*x^-2
       + ln(z)*ln(1 + x*z^-1)*x^-2*z
       - 2*ln(z)*ln(1 + x*z^-1)*x^-1*[1+x]^-1*z^2
       + 1/2*ln(z)*ln(1 + x*z^-1)*x^-1
       - ln(z)*ln(1 + x*z^-1)*x^-1*z
       + 2*ln(z)*ln(1 + x*z^-1)*x^-1*z^2
       + 1/2*ln(z)*ln(1 + x*z^-1)*[1+x]^-1
       - ln(z)*ln(1 + x*z^-1)*[1+x]^-1*z
       + ln(z)*ln(1 + x*z^-1)*z
       - 2*ln(z)*ln(1 + x*z^-1)*z^2
       + 1/2*ln(z)*ln(1 + x*z^-1)*x
       - ln(z)*ln(1 + x*z^-1)*x*z
       + ln(z)*ln(1 + x*z)
       - 1/2*ln(z)*ln(1 + x*z)*x^-2*[1+x]^-1
       + ln(z)*ln(1 + x*z)*x^-2*[1+x]^-1*z
       + 1/2*ln(z)*ln(1 + x*z)*x^-2
       - ln(z)*ln(1 + x*z)*x^-2*z
       + 2*ln(z)*ln(1 + x*z)*x^-1*[1+x]^-1*z^2
       - 1/2*ln(z)*ln(1 + x*z)*x^-1
       + ln(z)*ln(1 + x*z)*x^-1*z
       - 2*ln(z)*ln(1 + x*z)*x^-1*z^2
       - ln(z)*ln(1 + x*z)*[1-x]^-1
       - 2*ln(z)*ln(1 + x*z)*[1-x]^-1*z
       - 2*ln(z)*ln(1 + x*z)*[1-x]^-1*z^2
       - 1/2*ln(z)*ln(1 + x*z)*[1+x]^-1
       + ln(z)*ln(1 + x*z)*[1+x]^-1*z
       + 4*ln(z)*ln(1 + x*z)*z^2
       + 2*ln(z)*ln(1 + x*z)*x*z
       + 1/2*ln(z)*ln(z + x)
       - ln(z)*ln(z + x)*[1-x]^-1
       - 2*ln(z)*ln(z + x)*[1-x]^-1*z
       - 2*ln(z)*ln(z + x)*[1-x]^-1*z^2
       + ln(z)*ln(z + x)*z
       + 2*ln(z)*ln(z + x)*z^2
       + 1/2*ln(z)*ln(z + x)*x
       + ln(z)*ln(z + x)*x*z
       - ln(z)^2
       + 1/4*ln(z)^2*x^-2*[1+x]^-1
       - 1/2*ln(z)^2*x^-2*[1+x]^-1*z
       - 1/4*ln(z)^2*x^-2
       + 1/2*ln(z)^2*x^-2*z
       - ln(z)^2*x^-1*[1+x]^-1*z^2
       + 1/4*ln(z)^2*x^-1
       - 1/2*ln(z)^2*x^-1*z
       + ln(z)^2*x^-1*z^2
       + 1/2*ln(z)^2*[1-x]^-1
       - ln(z)^2*[1-x]^-1*z
       + ln(z)^2*[1-x]^-1*z^2
       + 1/4*ln(z)^2*[1+x]^-1
       - 1/2*ln(z)^2*[1+x]^-1*z
       + 2*ln(z)^2*z
       - 2*ln(z)^2*z^2
       + 1/2*ln(z)^2*x
       - ln(z)^2*x*z
       + ln(z)*ln([1-z])
       - 1/2*ln(z)*ln([1-z])*[1-x]^-1
       + ln(z)*ln([1-z])*[1-x]^-1*z
       - 2*ln(z)*ln([1-z])*z
       - 4*ln(sqrtxz3)*ArcTan(sqrtxz3)*z*sqrtxz3
       + 1/8*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*sqrtxz2^-1*poly2^-1
       - 1/8*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*sqrtxz2^-1
       + 5/4*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x*sqrtxz2^-1
       - 5/2*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x*z*sqrtxz2^-1
       - 1/4*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x^2*sqrtxz2^-1*poly2^-1
       - 1/8*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x^2*sqrtxz2^-1
       + 1/8*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x^4*sqrtxz2^-1*poly2^-1
       - 1/8*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*sqrtxz2^-1*poly2^-1
       + 1/8*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*sqrtxz2^-1
       - 5/4*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x*sqrtxz2^-1
       + 5/2*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x*z*sqrtxz2^-1
       + 1/4*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x^2*sqrtxz2^-1*poly2^-1
       + 1/8*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x^2*sqrtxz2^-1
       - 1/8*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x^4*sqrtxz2^-1*poly2^-1
       - 1/2*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1
       - 3*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*z
       + 2*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z
       + Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x
       + 1/2*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1
       + 3*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*z
       - 2*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z
       - Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x
       - 1/8*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*sqrtxz2^-1*poly2^-1
       + 1/8*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*sqrtxz2^-1
       - 5/4*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x*sqrtxz2^-1
       + 5/2*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x*z*sqrtxz2^-1
       + 1/4*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x^2*sqrtxz2^-1*poly2^-1
       + 1/8*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x^2*sqrtxz2^-1
       - 1/8*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x^4*sqrtxz2^-1*poly2^-1
       + 1/8*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*sqrtxz2^-1*poly2^-1
       - 1/8*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*sqrtxz2^-1
       + 5/4*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x*sqrtxz2^-1
       - 5/2*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x*z*sqrtxz2^-1
       - 1/4*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x^2*sqrtxz2^-1*poly2^-1
       - 1/8*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x^2*sqrtxz2^-1
       + 1/8*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x^4*sqrtxz2^-1*poly2^-1
       - Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)
       + 3/2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*[1-x]^-1
       + Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*[1-x]^-1*z
       + 4*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*[1-x]^-1*z^2
       - 4*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*z^2
       - 2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*x*z
       + Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)
       - 3/2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*[1-x]^-1
       - Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*[1-x]^-1*z
       - 4*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*[1-x]^-1*z^2
       + 4*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*z^2
       + 2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*x*z
       - Li2(1 - x*z^-1)
       + 1/2*Li2(1 - x*z^-1)*[1-x]^-1
       - Li2(1 - x*z^-1)*[1-x]^-1*z
       + 2*Li2(1 - x*z^-1)*z
       - 1/2*Li2( - x*z^-1)*x^-2*[1+x]^-1
       + Li2( - x*z^-1)*x^-2*[1+x]^-1*z
       + 1/2*Li2( - x*z^-1)*x^-2
       - Li2( - x*z^-1)*x^-2*z
       + 2*Li2( - x*z^-1)*x^-1*[1+x]^-1*z^2
       - 1/2*Li2( - x*z^-1)*x^-1
       + Li2( - x*z^-1)*x^-1*z
       - 2*Li2( - x*z^-1)*x^-1*z^2
       + Li2( - x*z^-1)*[1-x]^-1
       + 2*Li2( - x*z^-1)*[1-x]^-1*z
       + 2*Li2( - x*z^-1)*[1-x]^-1*z^2
       - 1/2*Li2( - x*z^-1)*[1+x]^-1
       + Li2( - x*z^-1)*[1+x]^-1*z
       - 2*Li2( - x*z^-1)*z
       - Li2( - x*z^-1)*x
       + 2*Li2( - x)
       + Li2( - x)*x^-2*[1+x]^-1
       - 2*Li2( - x)*x^-2*[1+x]^-1*z
       - Li2( - x)*x^-2
       + 2*Li2( - x)*x^-2*z
       - 4*Li2( - x)*x^-1*[1+x]^-1*z^2
       + Li2( - x)*x^-1
       - 2*Li2( - x)*x^-1*z
       + 4*Li2( - x)*x^-1*z^2
       - Li2( - x)*[1+x]^-1
       + 2*Li2( - x)*[1+x]^-1*z
       - 4*Li2( - x)*[1+x]^-1*z^2
       - 4*Li2( - x)*z
       + 2*Li2( - x)*x
       - 4*Li2( - x)*x*z
       + Li2( - x*z)
       - 1/2*Li2( - x*z)*x^-2*[1+x]^-1
       + Li2( - x*z)*x^-2*[1+x]^-1*z
       + 1/2*Li2( - x*z)*x^-2
       - Li2( - x*z)*x^-2*z
       + 2*Li2( - x*z)*x^-1*[1+x]^-1*z^2
       - 1/2*Li2( - x*z)*x^-1
       + Li2( - x*z)*x^-1*z
       - 2*Li2( - x*z)*x^-1*z^2
       - Li2( - x*z)*[1-x]^-1
       - 2*Li2( - x*z)*[1-x]^-1*z
       - 2*Li2( - x*z)*[1-x]^-1*z^2
       - 1/2*Li2( - x*z)*[1+x]^-1
       + Li2( - x*z)*[1+x]^-1*z
       + 4*Li2( - x*z)*z^2
       + 2*Li2( - x*z)*x*z
       + Li2(x)
       - Li2(x)*x^-2*[1+x]^-1
       + 2*Li2(x)*x^-2*[1+x]^-1*z
       + Li2(x)*x^-2
       - 2*Li2(x)*x^-2*z
       + 4*Li2(x)*x^-1*[1+x]^-1*z^2
       - Li2(x)*x^-1
       + 2*Li2(x)*x^-1*z
       - 4*Li2(x)*x^-1*z^2
       - 3/2*Li2(x)*[1-x]^-1
       + 3*Li2(x)*[1-x]^-1*z
       - 4*Li2(x)*[1-x]^-1*z^2
       + Li2(x)*[1+x]^-1
       - 2*Li2(x)*[1+x]^-1*z
       + 4*Li2(x)*[1+x]^-1*z^2
       - 2*Li2(x)*z
       + 4*Li2(x)*z^2
       + Li2(z)
       - 1/2*Li2(z)*[1-x]^-1
       + Li2(z)*[1-x]^-1*z
       - 2*Li2(z)*z
       - 2*InvTanInt( - sqrtxz3)*z*sqrtxz3
       - 4*InvTanInt(z*sqrtxz3)*z*sqrtxz3
       + 2*InvTanInt(sqrtxz3)*z*sqrtxz3
       )

       + NC * ( 5/2
       + 3*[1-x]^-1*rln2^2
       + 3*[1-x]^-1*sqrtxz1*rln2
       + 2*[1-x]^-1*z*rln2^2
       + 8*[1-x]^-1*z^2*rln2^2
       - 2*rln2^2
       - 2*sqrtxz1*rln2
       - 2*z
       - 8*z^2*rln2^2
       - 4*x*z*rln2^2
       - 1/6*pi^2*x^-2*[1+x]^-1
       + 1/3*pi^2*x^-2*[1+x]^-1*z
       + 1/6*pi^2*x^-2
       - 1/3*pi^2*x^-2*z
       + 2/3*pi^2*x^-1*[1+x]^-1*z^2
       - 1/6*pi^2*x^-1
       + 1/3*pi^2*x^-1*z
       - 2/3*pi^2*x^-1*z^2
       - 5/12*pi^2*[1-x]^-1
       + 5/6*pi^2*[1-x]^-1*z
       - pi^2*[1-x]^-1*z^2
       + 1/3*pi^2*[1+x]^-1
       - 2/3*pi^2*[1+x]^-1*z
       + pi^2*[1+x]^-1*z^2
       + 1/6*pi^2
       - 1/3*pi^2*z
       + 2/3*pi^2*z^2
       - 1/6*pi^2*x
       + 1/3*pi^2*x*z
       - 4*ln(1 + sqrtxz1 - z)*[1-x]^-1*rln2
       - 3*ln(1 + sqrtxz1 - z)*[1-x]^-1*sqrtxz1
       - 12*ln(1 + sqrtxz1 - z)*[1-x]^-1*z^2*rln2
       + 3*ln(1 + sqrtxz1 - z)*rln2
       + 2*ln(1 + sqrtxz1 - z)*sqrtxz1
       - 2*ln(1 + sqrtxz1 - z)*z*rln2
       + 12*ln(1 + sqrtxz1 - z)*z^2*rln2
       - ln(1 + sqrtxz1 - z)*x*rln2
       + 6*ln(1 + sqrtxz1 - z)*x*z*rln2
       - ln(1 + sqrtxz1 - z)^2
       + ln(1 + sqrtxz1 - z)^2*[1-x]^-1
       - 2*ln(1 + sqrtxz1 - z)^2*[1-x]^-1*z
       + 4*ln(1 + sqrtxz1 - z)^2*[1-x]^-1*z^2
       + 2*ln(1 + sqrtxz1 - z)^2*z
       - 4*ln(1 + sqrtxz1 - z)^2*z^2
       + ln(1 + sqrtxz1 - z)^2*x
       - 2*ln(1 + sqrtxz1 - z)^2*x*z
       - ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)
       + 2*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*[1-x]^-1
       + 4*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z
       + 4*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^2
       - 2*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*z
       - 4*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*z^2
       - ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*x
       - 2*ln(1 + sqrtxz1 - z)*ln(1 + sqrtxz1 + z)*x*z
       - 2*ln(1 + sqrtxz1 + z)*[1-x]^-1*rln2
       - 4*ln(1 + sqrtxz1 + z)*[1-x]^-1*z*rln2
       - 4*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^2*rln2
       + ln(1 + sqrtxz1 + z)*rln2
       + 2*ln(1 + sqrtxz1 + z)*z*rln2
       + 4*ln(1 + sqrtxz1 + z)*z^2*rln2
       + ln(1 + sqrtxz1 + z)*x*rln2
       + 2*ln(1 + sqrtxz1 + z)*x*z*rln2
       - 4*ln(z*sqrtxz3)*ArcTan(z*sqrtxz3)*z*sqrtxz3
       + 3/4*ln(x)
       - 2*ln(x)*x^-2*[1+x]^-1
       + 4*ln(x)*x^-2*[1+x]^-1*z
       + 2*ln(x)*x^-2
       - 4*ln(x)*x^-2*z
       + 8*ln(x)*x^-1*[1+x]^-1*z^2
       - 2*ln(x)*x^-1
       + 4*ln(x)*x^-1*z
       - 8*ln(x)*x^-1*z^2
       - ln(x)*[1-x]^-1
       + 5/2*ln(x)*[1-x]^-1*rln2
       + 3/2*ln(x)*[1-x]^-1*sqrtxz1
       + 3/2*ln(x)*[1-x]^-1*z
       - ln(x)*[1-x]^-1*z*rln2
       + 8*ln(x)*[1-x]^-1*z^2*rln2
       + 2*ln(x)*[1+x]^-1
       - 4*ln(x)*[1+x]^-1*z
       + 8*ln(x)*[1+x]^-1*z^2
       - 1/4*ln(x)*poly2^-1
       - 2*ln(x)*rln2
       - ln(x)*sqrtxz1
       - ln(x)*z
       + 2*ln(x)*z*rln2
       - 8*ln(x)*z^2*rln2
       - 1/2*ln(x)*x*[x-z]^-1
       + 1/4*ln(x)*x*poly2^-1
       - 3/4*ln(x)*x
       + ln(x)*x*rln2
       - 4*ln(x)*x*z*rln2
       + ln(x)*x^2*[x-z]^-1
       + 1/4*ln(x)*x^2*poly2^-1
       - 1/4*ln(x)*x^3*poly2^-1
       - 1/8*ln(x)*ln(1 - sqrtxz2 + x)*sqrtxz2^-1*poly2^-1
       + 1/8*ln(x)*ln(1 - sqrtxz2 + x)*sqrtxz2^-1
       - 5/4*ln(x)*ln(1 - sqrtxz2 + x)*x*sqrtxz2^-1
       + 5/2*ln(x)*ln(1 - sqrtxz2 + x)*x*z*sqrtxz2^-1
       + 1/4*ln(x)*ln(1 - sqrtxz2 + x)*x^2*sqrtxz2^-1*poly2^-1
       + 1/8*ln(x)*ln(1 - sqrtxz2 + x)*x^2*sqrtxz2^-1
       - 1/8*ln(x)*ln(1 - sqrtxz2 + x)*x^4*sqrtxz2^-1*poly2^-1
       + 1/8*ln(x)*ln(1 + sqrtxz2 + x)*sqrtxz2^-1*poly2^-1
       - 1/8*ln(x)*ln(1 + sqrtxz2 + x)*sqrtxz2^-1
       + 5/4*ln(x)*ln(1 + sqrtxz2 + x)*x*sqrtxz2^-1
       - 5/2*ln(x)*ln(1 + sqrtxz2 + x)*x*z*sqrtxz2^-1
       - 1/4*ln(x)*ln(1 + sqrtxz2 + x)*x^2*sqrtxz2^-1*poly2^-1
       - 1/8*ln(x)*ln(1 + sqrtxz2 + x)*x^2*sqrtxz2^-1
       + 1/8*ln(x)*ln(1 + sqrtxz2 + x)*x^4*sqrtxz2^-1*poly2^-1
       + ln(x)*ln(1 + sqrtxz1 - z)
       - ln(x)*ln(1 + sqrtxz1 - z)*[1-x]^-1
       + 2*ln(x)*ln(1 + sqrtxz1 - z)*[1-x]^-1*z
       - 4*ln(x)*ln(1 + sqrtxz1 - z)*[1-x]^-1*z^2
       - 2*ln(x)*ln(1 + sqrtxz1 - z)*z
       + 4*ln(x)*ln(1 + sqrtxz1 - z)*z^2
       - ln(x)*ln(1 + sqrtxz1 - z)*x
       + 2*ln(x)*ln(1 + sqrtxz1 - z)*x*z
       + ln(x)*ln(1 + sqrtxz1 + z)
       - 3/2*ln(x)*ln(1 + sqrtxz1 + z)*[1-x]^-1
       - ln(x)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z
       - 4*ln(x)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^2
       + 4*ln(x)*ln(1 + sqrtxz1 + z)*z^2
       + 2*ln(x)*ln(1 + sqrtxz1 + z)*x*z
       - 1/2*ln(x)*ln(1 + x*z^-1)
       + 1/2*ln(x)*ln(1 + x*z^-1)*x^-2*[1+x]^-1
       - ln(x)*ln(1 + x*z^-1)*x^-2*[1+x]^-1*z
       - 1/2*ln(x)*ln(1 + x*z^-1)*x^-2
       + ln(x)*ln(1 + x*z^-1)*x^-2*z
       - 2*ln(x)*ln(1 + x*z^-1)*x^-1*[1+x]^-1*z^2
       + 1/2*ln(x)*ln(1 + x*z^-1)*x^-1
       - ln(x)*ln(1 + x*z^-1)*x^-1*z
       + 2*ln(x)*ln(1 + x*z^-1)*x^-1*z^2
       + 1/2*ln(x)*ln(1 + x*z^-1)*[1+x]^-1
       - ln(x)*ln(1 + x*z^-1)*[1+x]^-1*z
       + ln(x)*ln(1 + x*z^-1)*z
       - 2*ln(x)*ln(1 + x*z^-1)*z^2
       + 1/2*ln(x)*ln(1 + x*z^-1)*x
       - ln(x)*ln(1 + x*z^-1)*x*z
       - 3*ln(x)*ln(1 + x)
       + 2*ln(x)*ln(1 + x)*[1+x]^-1
       - 4*ln(x)*ln(1 + x)*[1+x]^-1*z
       + 4*ln(x)*ln(1 + x)*[1+x]^-1*z^2
       + 6*ln(x)*ln(1 + x)*z
       - 4*ln(x)*ln(1 + x)*z^2
       - ln(x)*ln(1 + x)*x
       + 2*ln(x)*ln(1 + x)*x*z
       - ln(x)*ln(1 + x*z)
       + 1/2*ln(x)*ln(1 + x*z)*x^-2*[1+x]^-1
       - ln(x)*ln(1 + x*z)*x^-2*[1+x]^-1*z
       - 1/2*ln(x)*ln(1 + x*z)*x^-2
       + ln(x)*ln(1 + x*z)*x^-2*z
       - 2*ln(x)*ln(1 + x*z)*x^-1*[1+x]^-1*z^2
       + 1/2*ln(x)*ln(1 + x*z)*x^-1
       - ln(x)*ln(1 + x*z)*x^-1*z
       + 2*ln(x)*ln(1 + x*z)*x^-1*z^2
       + ln(x)*ln(1 + x*z)*[1-x]^-1
       + 2*ln(x)*ln(1 + x*z)*[1-x]^-1*z
       + 2*ln(x)*ln(1 + x*z)*[1-x]^-1*z^2
       + 1/2*ln(x)*ln(1 + x*z)*[1+x]^-1
       - ln(x)*ln(1 + x*z)*[1+x]^-1*z
       - 4*ln(x)*ln(1 + x*z)*z^2
       - 2*ln(x)*ln(1 + x*z)*x*z
       + 1/2*ln(x)*ln(z + x)
       - ln(x)*ln(z + x)*[1-x]^-1
       - 2*ln(x)*ln(z + x)*[1-x]^-1*z
       - 2*ln(x)*ln(z + x)*[1-x]^-1*z^2
       + ln(x)*ln(z + x)*z
       + 2*ln(x)*ln(z + x)*z^2
       + 1/2*ln(x)*ln(z + x)*x
       + ln(x)*ln(z + x)*x*z
       + 3/2*ln(x)^2
       + 3/4*ln(x)^2*x^-2*[1+x]^-1
       - 3/2*ln(x)^2*x^-2*[1+x]^-1*z
       - 3/4*ln(x)^2*x^-2
       + 3/2*ln(x)^2*x^-2*z
       - 3*ln(x)^2*x^-1*[1+x]^-1*z^2
       + 3/4*ln(x)^2*x^-1
       - 3/2*ln(x)^2*x^-1*z
       + 3*ln(x)^2*x^-1*z^2
       - ln(x)^2*[1-x]^-1
       + 2*ln(x)^2*[1-x]^-1*z
       - ln(x)^2*[1-x]^-1*z^2
       - 5/4*ln(x)^2*[1+x]^-1
       + 5/2*ln(x)^2*[1+x]^-1*z
       - 4*ln(x)^2*[1+x]^-1*z^2
       - 3*ln(x)^2*z
       + 2*ln(x)^2*z^2
       + ln(x)^2*x
       - 2*ln(x)^2*x*z
       - ln(x)*ln([1-x])
       + ln(x)*ln([1-x])*x^-2*[1+x]^-1
       - 2*ln(x)*ln([1-x])*x^-2*[1+x]^-1*z
       - ln(x)*ln([1-x])*x^-2
       + 2*ln(x)*ln([1-x])*x^-2*z
       - 4*ln(x)*ln([1-x])*x^-1*[1+x]^-1*z^2
       + ln(x)*ln([1-x])*x^-1
       - 2*ln(x)*ln([1-x])*x^-1*z
       + 4*ln(x)*ln([1-x])*x^-1*z^2
       + 3/2*ln(x)*ln([1-x])*[1-x]^-1
       - 3*ln(x)*ln([1-x])*[1-x]^-1*z
       + 4*ln(x)*ln([1-x])*[1-x]^-1*z^2
       - ln(x)*ln([1-x])*[1+x]^-1
       + 2*ln(x)*ln([1-x])*[1+x]^-1*z
       - 4*ln(x)*ln([1-x])*[1+x]^-1*z^2
       + 2*ln(x)*ln([1-x])*z
       - 4*ln(x)*ln([1-x])*z^2
       + ln(x)*ln([1+x])
       - ln(x)*ln([1+x])*x^-2*[1+x]^-1
       + 2*ln(x)*ln([1+x])*x^-2*[1+x]^-1*z
       + ln(x)*ln([1+x])*x^-2
       - 2*ln(x)*ln([1+x])*x^-2*z
       + 4*ln(x)*ln([1+x])*x^-1*[1+x]^-1*z^2
       - ln(x)*ln([1+x])*x^-1
       + 2*ln(x)*ln([1+x])*x^-1*z
       - 4*ln(x)*ln([1+x])*x^-1*z^2
       - ln(x)*ln([1+x])*[1+x]^-1
       + 2*ln(x)*ln([1+x])*[1+x]^-1*z
       - 2*ln(x)*ln([1+x])*z
       + 4*ln(x)*ln([1+x])*z^2
       - ln(x)*ln([1+x])*x
       + 2*ln(x)*ln([1+x])*x*z
       - 1/2*ln(x)*ln(z)
       - 1/2*ln(x)*ln(z)*x^-2*[1+x]^-1
       + ln(x)*ln(z)*x^-2*[1+x]^-1*z
       + 1/2*ln(x)*ln(z)*x^-2
       - ln(x)*ln(z)*x^-2*z
       + 2*ln(x)*ln(z)*x^-1*[1+x]^-1*z^2
       - 1/2*ln(x)*ln(z)*x^-1
       + ln(x)*ln(z)*x^-1*z
       - 2*ln(x)*ln(z)*x^-1*z^2
       + 2*ln(x)*ln(z)*[1-x]^-1
       - 2*ln(x)*ln(z)*[1-x]^-1*z
       + 4*ln(x)*ln(z)*[1-x]^-1*z^2
       - 1/2*ln(x)*ln(z)*[1+x]^-1
       + ln(x)*ln(z)*[1+x]^-1*z
       + ln(x)*ln(z)*z
       - 2*ln(x)*ln(z)*z^2
       + 1/2*ln(x)*ln(z)*x
       - ln(x)*ln(z)*x*z
       + 9/4*ln(z)
       - 3/2*ln(z)*[1-x]^-1
       + 7/2*ln(z)*[1-x]^-1*rln2
       + 3/2*ln(z)*[1-x]^-1*sqrtxz1
       - 3/2*ln(z)*[1-x]^-1*z
       + 5*ln(z)*[1-x]^-1*z*rln2
       + 8*ln(z)*[1-x]^-1*z^2*rln2
       - 1/4*ln(z)*poly2^-1
       - 2*ln(z)*rln2
       - ln(z)*sqrtxz1
       + ln(z)*z
       - 2*ln(z)*z*rln2
       - 8*ln(z)*z^2*rln2
       + 1/2*ln(z)*x*[x-z]^-1
       - 1/4*ln(z)*x*poly2^-1
       + 3/4*ln(z)*x
       - ln(z)*x*rln2
       - 4*ln(z)*x*z*rln2
       - ln(z)*x^2*[x-z]^-1
       + 1/4*ln(z)*x^2*poly2^-1
       + 1/4*ln(z)*x^3*poly2^-1
       + ln(z)*ln(1 + sqrtxz1 - z)
       - 3/2*ln(z)*ln(1 + sqrtxz1 - z)*[1-x]^-1
       - ln(z)*ln(1 + sqrtxz1 - z)*[1-x]^-1*z
       - 4*ln(z)*ln(1 + sqrtxz1 - z)*[1-x]^-1*z^2
       + 4*ln(z)*ln(1 + sqrtxz1 - z)*z^2
       + 2*ln(z)*ln(1 + sqrtxz1 - z)*x*z
       + ln(z)*ln(1 + sqrtxz1 + z)
       - 2*ln(z)*ln(1 + sqrtxz1 + z)*[1-x]^-1
       - 4*ln(z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z
       - 4*ln(z)*ln(1 + sqrtxz1 + z)*[1-x]^-1*z^2
       + 2*ln(z)*ln(1 + sqrtxz1 + z)*z
       + 4*ln(z)*ln(1 + sqrtxz1 + z)*z^2
       + ln(z)*ln(1 + sqrtxz1 + z)*x
       + 2*ln(z)*ln(1 + sqrtxz1 + z)*x*z
       + 1/2*ln(z)*ln(1 + x*z^-1)
       - 1/2*ln(z)*ln(1 + x*z^-1)*x^-2*[1+x]^-1
       + ln(z)*ln(1 + x*z^-1)*x^-2*[1+x]^-1*z
       + 1/2*ln(z)*ln(1 + x*z^-1)*x^-2
       - ln(z)*ln(1 + x*z^-1)*x^-2*z
       + 2*ln(z)*ln(1 + x*z^-1)*x^-1*[1+x]^-1*z^2
       - 1/2*ln(z)*ln(1 + x*z^-1)*x^-1
       + ln(z)*ln(1 + x*z^-1)*x^-1*z
       - 2*ln(z)*ln(1 + x*z^-1)*x^-1*z^2
       - 1/2*ln(z)*ln(1 + x*z^-1)*[1+x]^-1
       + ln(z)*ln(1 + x*z^-1)*[1+x]^-1*z
       - ln(z)*ln(1 + x*z^-1)*z
       + 2*ln(z)*ln(1 + x*z^-1)*z^2
       - 1/2*ln(z)*ln(1 + x*z^-1)*x
       + ln(z)*ln(1 + x*z^-1)*x*z
       - ln(z)*ln(1 + x*z)
       + 1/2*ln(z)*ln(1 + x*z)*x^-2*[1+x]^-1
       - ln(z)*ln(1 + x*z)*x^-2*[1+x]^-1*z
       - 1/2*ln(z)*ln(1 + x*z)*x^-2
       + ln(z)*ln(1 + x*z)*x^-2*z
       - 2*ln(z)*ln(1 + x*z)*x^-1*[1+x]^-1*z^2
       + 1/2*ln(z)*ln(1 + x*z)*x^-1
       - ln(z)*ln(1 + x*z)*x^-1*z
       + 2*ln(z)*ln(1 + x*z)*x^-1*z^2
       + ln(z)*ln(1 + x*z)*[1-x]^-1
       + 2*ln(z)*ln(1 + x*z)*[1-x]^-1*z
       + 2*ln(z)*ln(1 + x*z)*[1-x]^-1*z^2
       + 1/2*ln(z)*ln(1 + x*z)*[1+x]^-1
       - ln(z)*ln(1 + x*z)*[1+x]^-1*z
       - 4*ln(z)*ln(1 + x*z)*z^2
       - 2*ln(z)*ln(1 + x*z)*x*z
       - 1/2*ln(z)*ln(z + x)
       + ln(z)*ln(z + x)*[1-x]^-1
       + 2*ln(z)*ln(z + x)*[1-x]^-1*z
       + 2*ln(z)*ln(z + x)*[1-x]^-1*z^2
       - ln(z)*ln(z + x)*z
       - 2*ln(z)*ln(z + x)*z^2
       - 1/2*ln(z)*ln(z + x)*x
       - ln(z)*ln(z + x)*x*z
       + ln(z)^2
       - 1/4*ln(z)^2*x^-2*[1+x]^-1
       + 1/2*ln(z)^2*x^-2*[1+x]^-1*z
       + 1/4*ln(z)^2*x^-2
       - 1/2*ln(z)^2*x^-2*z
       + ln(z)^2*x^-1*[1+x]^-1*z^2
       - 1/4*ln(z)^2*x^-1
       + 1/2*ln(z)^2*x^-1*z
       - ln(z)^2*x^-1*z^2
       - 1/2*ln(z)^2*[1-x]^-1
       + ln(z)^2*[1-x]^-1*z
       - ln(z)^2*[1-x]^-1*z^2
       - 1/4*ln(z)^2*[1+x]^-1
       + 1/2*ln(z)^2*[1+x]^-1*z
       - 2*ln(z)^2*z
       + 2*ln(z)^2*z^2
       - 1/2*ln(z)^2*x
       + ln(z)^2*x*z
       - ln(z)*ln([1-z])
       + 1/2*ln(z)*ln([1-z])*[1-x]^-1
       - ln(z)*ln([1-z])*[1-x]^-1*z
       + 2*ln(z)*ln([1-z])*z
       + 4*ln(sqrtxz3)*ArcTan(sqrtxz3)*z*sqrtxz3
       - 1/8*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*sqrtxz2^-1*poly2^-1
       + 1/8*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*sqrtxz2^-1
       - 5/4*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x*sqrtxz2^-1
       + 5/2*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x*z*sqrtxz2^-1
       + 1/4*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x^2*sqrtxz2^-1*poly2^-1
       + 1/8*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x^2*sqrtxz2^-1
       - 1/8*Li2(1/2 - 1/2*x^-1 - 1/2*x^-1*sqrtxz2)*x^4*sqrtxz2^-1*poly2^-1
       + 1/8*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*sqrtxz2^-1*poly2^-1
       - 1/8*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*sqrtxz2^-1
       + 5/4*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x*sqrtxz2^-1
       - 5/2*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x*z*sqrtxz2^-1
       - 1/4*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x^2*sqrtxz2^-1*poly2^-1
       - 1/8*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x^2*sqrtxz2^-1
       + 1/8*Li2(1/2 - 1/2*x^-1 + 1/2*x^-1*sqrtxz2)*x^4*sqrtxz2^-1*poly2^-1
       + 1/2*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1
       + 3*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*z
       - 2*Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z
       - Li2(1/2 - 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x
       - 1/2*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1
       - 3*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*[1-x]^-1*z
       + 2*Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*z
       + Li2(1/2 + 1/2*z^-1 - 1/2*z^-1*sqrtxz1)*x
       + 1/8*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*sqrtxz2^-1*poly2^-1
       - 1/8*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*sqrtxz2^-1
       + 5/4*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x*sqrtxz2^-1
       - 5/2*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x*z*sqrtxz2^-1
       - 1/4*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x^2*sqrtxz2^-1*poly2^-1
       - 1/8*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x^2*sqrtxz2^-1
       + 1/8*Li2(1/2 - 1/2*sqrtxz2 - 1/2*x)*x^4*sqrtxz2^-1*poly2^-1
       - 1/8*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*sqrtxz2^-1*poly2^-1
       + 1/8*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*sqrtxz2^-1
       - 5/4*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x*sqrtxz2^-1
       + 5/2*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x*z*sqrtxz2^-1
       + 1/4*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x^2*sqrtxz2^-1*poly2^-1
       + 1/8*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x^2*sqrtxz2^-1
       - 1/8*Li2(1/2 + 1/2*sqrtxz2 - 1/2*x)*x^4*sqrtxz2^-1*poly2^-1
       + Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)
       - 3/2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*[1-x]^-1
       - Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*[1-x]^-1*z
       - 4*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*[1-x]^-1*z^2
       + 4*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*z^2
       + 2*Li2(1/2 - 1/2*sqrtxz1 - 1/2*z)*x*z
       - Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)
       + 3/2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*[1-x]^-1
       + Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*[1-x]^-1*z
       + 4*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*[1-x]^-1*z^2
       - 4*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*z^2
       - 2*Li2(1/2 - 1/2*sqrtxz1 + 1/2*z)*x*z
       + Li2(1 - x*z^-1)
       - 1/2*Li2(1 - x*z^-1)*[1-x]^-1
       + Li2(1 - x*z^-1)*[1-x]^-1*z
       - 2*Li2(1 - x*z^-1)*z
       + 1/2*Li2( - x*z^-1)*x^-2*[1+x]^-1
       - Li2( - x*z^-1)*x^-2*[1+x]^-1*z
       - 1/2*Li2( - x*z^-1)*x^-2
       + Li2( - x*z^-1)*x^-2*z
       - 2*Li2( - x*z^-1)*x^-1*[1+x]^-1*z^2
       + 1/2*Li2( - x*z^-1)*x^-1
       - Li2( - x*z^-1)*x^-1*z
       + 2*Li2( - x*z^-1)*x^-1*z^2
       - Li2( - x*z^-1)*[1-x]^-1
       - 2*Li2( - x*z^-1)*[1-x]^-1*z
       - 2*Li2( - x*z^-1)*[1-x]^-1*z^2
       + 1/2*Li2( - x*z^-1)*[1+x]^-1
       - Li2( - x*z^-1)*[1+x]^-1*z
       + 2*Li2( - x*z^-1)*z
       + Li2( - x*z^-1)*x
       - 2*Li2( - x)
       - Li2( - x)*x^-2*[1+x]^-1
       + 2*Li2( - x)*x^-2*[1+x]^-1*z
       + Li2( - x)*x^-2
       - 2*Li2( - x)*x^-2*z
       + 4*Li2( - x)*x^-1*[1+x]^-1*z^2
       - Li2( - x)*x^-1
       + 2*Li2( - x)*x^-1*z
       - 4*Li2( - x)*x^-1*z^2
       + Li2( - x)*[1+x]^-1
       - 2*Li2( - x)*[1+x]^-1*z
       + 4*Li2( - x)*[1+x]^-1*z^2
       + 4*Li2( - x)*z
       - 2*Li2( - x)*x
       + 4*Li2( - x)*x*z
       - Li2( - x*z)
       + 1/2*Li2( - x*z)*x^-2*[1+x]^-1
       - Li2( - x*z)*x^-2*[1+x]^-1*z
       - 1/2*Li2( - x*z)*x^-2
       + Li2( - x*z)*x^-2*z
       - 2*Li2( - x*z)*x^-1*[1+x]^-1*z^2
       + 1/2*Li2( - x*z)*x^-1
       - Li2( - x*z)*x^-1*z
       + 2*Li2( - x*z)*x^-1*z^2
       + Li2( - x*z)*[1-x]^-1
       + 2*Li2( - x*z)*[1-x]^-1*z
       + 2*Li2( - x*z)*[1-x]^-1*z^2
       + 1/2*Li2( - x*z)*[1+x]^-1
       - Li2( - x*z)*[1+x]^-1*z
       - 4*Li2( - x*z)*z^2
       - 2*Li2( - x*z)*x*z
       - Li2(x)
       + Li2(x)*x^-2*[1+x]^-1
       - 2*Li2(x)*x^-2*[1+x]^-1*z
       - Li2(x)*x^-2
       + 2*Li2(x)*x^-2*z
       - 4*Li2(x)*x^-1*[1+x]^-1*z^2
       + Li2(x)*x^-1
       - 2*Li2(x)*x^-1*z
       + 4*Li2(x)*x^-1*z^2
       + 3/2*Li2(x)*[1-x]^-1
       - 3*Li2(x)*[1-x]^-1*z
       + 4*Li2(x)*[1-x]^-1*z^2
       - Li2(x)*[1+x]^-1
       + 2*Li2(x)*[1+x]^-1*z
       - 4*Li2(x)*[1+x]^-1*z^2
       + 2*Li2(x)*z
       - 4*Li2(x)*z^2
       - Li2(z)
       + 1/2*Li2(z)*[1-x]^-1
       - Li2(z)*[1-x]^-1*z
       + 2*Li2(z)*z
       + 2*InvTanInt( - sqrtxz3)*z*sqrtxz3
       + 4*InvTanInt(z*sqrtxz3)*z*sqrtxz3
       - 2*InvTanInt(sqrtxz3)*z*sqrtxz3 ) );

