******************************************************************************
**  hplog: a subroutine for the evaluation of harmonic polylogarithms
**         Version 1.0         12/07/2001
**  described in: 
**  T.Gehrmann and E.Remiddi: Numerical Evaluation of the Harmonic 
**                            Polylogarithms up to Weight 4
**                            (hep-ph/0107173; CERN-TH/2001/188)
**  the harmonic polylogarithms are defined in: 
**  E.Remiddi and J.Vermaseren: Harmonic Polylogarithms
**                            (hep-ph/9905237; Int.J.Mod.Phys. A15 (2000) 725)
**  email:
**  Thomas.Gehrmann@cern.ch and Ettore.Remiddi@bo.infn.it
**
******************************************************************************
c      subroutine hplog(x,nw,Hc1,Hc2,Hc3,Hc4, 
c     $                      Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
      subroutine hplog(x,nw,Hr1,Hr2,Hr3,Hr4,n1,n2) 
****** 
** x is the argument of the 1dHPL's (1 dimensional Harmonic PolyLogarithms) 
**   to be evaluated; 
** nw is the maximum weight of the required 1dHPL's; 
**    the maximum allowed value of nw of this implementation is 4; 
** Hc1,Hc2,Hc3,Hc4 are the complex*16 values of the 1dHPL; 
**    they must all be supplied in the arguments even if some of them 
**    are not to be evaluated; 
** Hr1,Hr2,Hr3,Hr4 are the double precision real parts of 
**    Hc1,Hc2,Hc3,Hc4; 
** Hi1,Hi2,Hi3,Hi4 are the double precision immaginary parts of 
**    Hc1,Hc2,Hc3,Hc4 divided by pi=3.114159.... 
** n1,n2 is the required range of indices, the allowed ranges are 
**    (0,1), (-1,0), (-1,1) ; 
****** 
      implicit double precision (a-h,o-z) 
      complex*16 Hc1,Hc2,Hc3,Hc4 
      dimension Hc1(n1:n2),Hc2(n1:n2,n1:n2),Hc3(n1:n2,n1:n2,n1:n2), 
     $          Hc4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension Hr1(n1:n2),Hr2(n1:n2,n1:n2),Hr3(n1:n2,n1:n2,n1:n2), 
     $          Hr4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2) 
      parameter (r2   = 1.4142135623730950488d0) 
** check on the weight nw 
      if ( (nw.lt.1).or.(nw.gt.4) ) then 
        print*, ' illegal call of eval1dhpl with second argument', 
     $          ' (the weight) = ',nw 
        print*, ' the allowed values of the weight are 1,2,3,4 ' 
        stop
      endif
** setting the immaginary parts equal to zero 
      call setzero(nw,Hi1,Hi2,Hi3,Hi4,n1,n2) 
** looking at the range of the argument 
*      r2 = sqrt(2.d0) 
      r2m1 = r2 - 1 
      r2p1 = r2 + 1 
      if ( ( x.gt.-r2m1 ).and.( x.le.r2m1) ) then 
*        print*, ' eval1dhpl:      x = ',x,', call eval1dhplat0 ' 
        call eval1dhplat0(x,nw,Hc1,Hc2,Hc3,Hc4, 
     $                          Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
        return 
      elseif ( x.eq.1d0 ) then
*        print*, ' eval1dhpl:      x = ',x,', call eval1dhplin1 ' 
        call eval1dhplin1(x,nw,Hc1,Hc2,Hc3,Hc4, 
     $                          Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
        return 
      elseif ( ( x.gt.r2m1 ).and.( x.le.r2p1) ) then 
*        print*, ' eval1dhpl:      x = ',x,', call eval1dhplat1 ' 
        call eval1dhplat1(x,nw,Hc1,Hc2,Hc3,Hc4, 
     $                          Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
        return 
      elseif ( ( x.gt.r2p1 ) ) then 
*        print*, ' eval1dhpl:      x = ',x,', call eval1dhplatinf ' 
        call eval1dhplatinf(x,nw,Hc1,Hc2,Hc3,Hc4, 
     $                          Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
        return 
      elseif ( ( x.le.-r2p1) ) then 
*        print*, ' eval1dhpl:      x = ',x,', call eval1dhplatminf ' 
        call eval1dhplatminf(x,nw,Hc1,Hc2,Hc3,Hc4, 
     $                          Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
        return 
      elseif ( x.eq.-1d0 ) then
*        print*, ' eval1dhpl:      x = ',x,', call eval1dhplinm1 ' 
        call eval1dhplinm1(x,nw,Hc1,Hc2,Hc3,Hc4, 
     $                          Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
        return 
      elseif ( ( x.gt.-r2p1 ).and.( x.le.-r2m1) ) then 
*        print*, ' eval1dhpl:      x = ',x,', call eval1dhplatm1 ' 
        call eval1dhplatm1(x,nw,Hc1,Hc2,Hc3,Hc4, 
     $                          Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
        return 
      endif 
** 
      end 
************************************************************************ 
      subroutine eval1dhplat0(y,nw,H1,H2,H3,H4, 
     $                          HY1,HY2,HY3,HY4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
** evaluates 1dhpl's in the 0-range  -(r2-1) < y <= (r2-1) 
** by direct series expansion (Bernoulli-accelerated) 
      implicit double precision (a-h,o-z) 
      complex*16 H1,H2,H3,H4 
      dimension H1(n1:n2),H2(n1:n2,n1:n2),H3(n1:n2,n1:n2,n1:n2), 
     $          H4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension HY1(n1:n2),HY2(n1:n2,n1:n2),HY3(n1:n2,n1:n2,n1:n2), 
     $          HY4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2) 
** evaluate the irreducible 1dHPL's first 
      call fillh1(y,H1,HY1,Hi1,n1,n2) 
      if ( nw.eq.1 ) return 
      call fillirr1dhplat0(y,nw,HY1,HY2,HY3,HY4,n1,n2) 
** then the reducible 1dHPL's 
      call fillred1dhpl(nw,H1,H2,H3,H4, 
     $                     HY1,HY2,HY3,HY4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
      return 
      end 
************************************************************************ 
      subroutine eval1dhplin1(y,nw,H1,H2,H3,H4, 
     $                          HY1,HY2,HY3,HY4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
** evaluates 1dhpl's for y=1 (explicit values are tabulated)
      implicit double precision (a-h,o-z) 
      complex*16 H1,H2,H3,H4 
      dimension H1(n1:n2),H2(n1:n2,n1:n2),H3(n1:n2,n1:n2,n1:n2), 
     $          H4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension HY1(n1:n2),HY2(n1:n2,n1:n2),HY3(n1:n2,n1:n2,n1:n2), 
     $          HY4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2) 
      parameter (pi   = 3.14159265358979324d0) 
** evaluate the irreducible 1dHPL's first 
      call fillh1(y,H1,HY1,Hi1,n1,n2) 
      if ( nw.eq.1 ) return 
      call fillirr1dhplin1(y,nw,HY1,HY2,HY3,HY4,n1,n2) 
** then the reducible 1dHPL's 
      call fillred1dhpl(nw,H1,H2,H3,H4, 
     $                     HY1,HY2,HY3,HY4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
      if (n2.eq.0) return
** correct the ill-defined entries
      HY2(1,0) = - HY2(0,1) 
      Hi2(1,0) = 0d0 
      H2(1,0) = dcmplx(HY2(1,0),Hi2(1,0)*pi)
      if ( nw.eq.2 ) return
      HY3(1,0,0) = HY3(0,0,1) 
      Hi3(1,0,0) = 0d0 
      H3(1,0,0) = dcmplx(HY3(1,0,0),Hi3(1,0,0)*pi)
      if ( nw.eq.3 ) return
      HY4(1,0,0,0) = -HY4(0,0,0,1) 
      Hi4(1,0,0,0) = 0d0 
      H4(1,0,0,0) = dcmplx(HY4(1,0,0,0),Hi4(1,0,0,0)*pi)
      return 
      end 
************************************************************************ 
      subroutine eval1dhplat1(y,nw,H1,H2,H3,H4, 
     $                          HY1,HY2,HY3,HY4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
** evaluates 1dhpl's in the 1-range  (r2-1) < y <= (r2+1) 
** evaluating first the H(..,r=(1-y)/(1+y)) by calling eval1dhplat0(r)  
** and then expressing H(..,y=(1-r)/(1+r)) in terms of H(..,r) 
      implicit double precision (a-h,o-z) 
      complex*16 H1,H2,H3,H4 
      dimension H1(n1:n2),H2(n1:n2,n1:n2),H3(n1:n2,n1:n2,n1:n2), 
     $          H4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension HY1(n1:n2),HY2(n1:n2,n1:n2),HY3(n1:n2,n1:n2,n1:n2), 
     $          HY4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2) 
** additional arrays required within this routine 
      dimension HR1(-1:1),HR2(-1:1,-1:1),HR3(-1:1,-1:1,-1:1), 
     $          HR4(-1:1,-1:1,-1:1,-1:1) 
** the nw = 1 case 
      call fillh1(y,H1,HY1,Hi1,n1,n2) 
      if ( nw.eq.1 ) return 
** the nw > 1 case 
      r = (1.d0-y)/(1.d0+y) 
*      print*,' eval1dhplat1: y = ',y,', r = ',r 
** the whole (-1,1) range is in general needed for any pair (n1,n2)
      call fillirr1dhplat0(r,nw,HR1,HR2,HR3,HR4,-1,1) 
** fillirr1dhplat1 takes care automatically of all the immaginary 
** parts as well as of the jump across y=1 
      call fillirr1dhplat1(r,nw,HR1,HR2,HR3,HR4, 
     $                          HY1,HY2,HY3,HY4, 
     $                          Hi1,Hi2,Hi3,Hi4,n1,n2) 
** then the reducible 1dHPL's 
      call fillred1dhpl(nw,H1,H2,H3,H4, 
     $                     HY1,HY2,HY3,HY4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
      return 
      end 
************************************************************************ 
      subroutine eval1dhplatinf(y,nw,H1,H2,H3,H4, 
     $                          HY1,HY2,HY3,HY4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
** evaluates 1dhpl's in the inf-range  (r2+1) < abs(y) 
** evaluating first the H(..,x=1/y) by calling eval1dhplat0(x)  
** and then expressing H(..,y=1/x) in terms of H(..,x) 
      implicit double precision (a-h,o-z) 
      complex*16 H1,H2,H3,H4 
      dimension H1(n1:n2),H2(n1:n2,n1:n2),H3(n1:n2,n1:n2,n1:n2), 
     $          H4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension HY1(n1:n2),HY2(n1:n2,n1:n2),HY3(n1:n2,n1:n2,n1:n2), 
     $          HY4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2) 
** additional arrays required within this routine 
      dimension HX1(n1:n2),HX2(n1:n2,n1:n2),HX3(n1:n2,n1:n2,n1:n2), 
     $          HX4(n1:n2,n1:n2,n1:n2,n1:n2) 
      parameter (pi   = 3.14159265358979324d0) 
** the nw = 1 case 
      call fillh1(y,H1,HY1,Hi1,n1,n2) 
      if ( nw.eq.1 ) return 
** the nw > 1 case 
      x = 1.d0/y 
*      print*,' eval1dhplatinf: y = ',y,', x = ',x 
      call fillirr1dhplat0(x,nw,HX1,HX2,HX3,HX4,n1,n2) 
** fillirr1dhplatinf takes care automatically of all the immaginary 
** parts as well as of the jump across y=1 
      call fillirr1dhplatinf(x,nw,HX1,HX2,HX3,HX4, 
     $                            HY1,HY2,HY3,HY4, 
     $                            Hi1,Hi2,Hi3,Hi4,n1,n2) 
** then the reducible 1dHPL's 
      call fillred1dhpl(nw,H1,H2,H3,H4, 
     $                     HY1,HY2,HY3,HY4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
      return 
      end 
************************************************************************ 
      subroutine eval1dhplinm1(y,nw,H1,H2,H3,H4, 
     $                           HY1,HY2,HY3,HY4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
** evaluates 1dhpl's for y=-1 (explicit values are tabulated)
      implicit double precision (a-h,o-z) 
      complex*16 H1,H2,H3,H4 
      complex*16 G1,G2,G3,G4  
      dimension H1(n1:n2),H2(n1:n2,n1:n2),H3(n1:n2,n1:n2,n1:n2), 
     $          H4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension HY1(n1:n2),HY2(n1:n2,n1:n2),HY3(n1:n2,n1:n2,n1:n2), 
     $          HY4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2) 
** additional arrays required within this routine 
      dimension G1(-n2:-n1),G2(-n2:-n1,-n2:-n1),
     $          G3(-n2:-n1,-n2:-n1,-n2:-n1), 
     $          G4(-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1) 
      dimension GY1(-n2:-n1),GY2(-n2:-n1,-n2:-n1),
     $          GY3(-n2:-n1,-n2:-n1,-n2:-n1), 
     $          GY4(-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1) 
      dimension Gi1(-n2:-n1),Gi2(-n2:-n1,-n2:-n1),
     $          Gi3(-n2:-n1,-n2:-n1,-n2:-n1), 
     $          Gi4(-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1) 
      dimension istorfill(3)
      dimension nphase(-1:1) 
      data nphase/-1,1,-1/ 
      parameter (pi   = 3.14159265358979324d0) 
*      print*,' eval1dhplatm1: y = ',y 
** evaluate H(...,-y) 
      call setzero(nw,Gi1,Gi2,Gi3,Gi4,-n2,-n1) 
      Gi1(0) = -1
      call eval1dhplin1(-y,nw,G1,G2,G3,G4, 
     $                        GY1,GY2,GY3,GY4,Gi1,Gi2,Gi3,Gi4,-n2,-n1) 
** fill the arrays H's 
      do k1=n1,n2 
        nph1 = nphase(k1) 
        HY1(k1) =   nph1*GY1(-k1) 
        Hi1(k1) = - nph1*Gi1(-k1) 
        H1(k1)  = dcmplx(HY1(k1),Hi1(k1)*pi) 
        if ( nw.gt.1 ) then 
          do k2=n1,n2 
            nph2 = nph1*nphase(k2) 
            HY2(k1,k2) =   nph2*GY2(-k1,-k2) 
            Hi2(k1,k2) = - nph2*Gi2(-k1,-k2) 
            H2(k1,k2)  = dcmplx(HY2(k1,k2),Hi2(k1,k2)*pi) 
            if ( nw.gt.2 ) then 
              do k3=n1,n2 
                nph3 = nph2*nphase(k3) 
                HY3(k1,k2,k3) =   nph3*GY3(-k1,-k2,-k3) 
                Hi3(k1,k2,k3) = - nph3*Gi3(-k1,-k2,-k3) 
                H3(k1,k2,k3)  = dcmplx(HY3(k1,k2,k3),Hi3(k1,k2,k3)*pi) 
                if ( nw.gt.3 ) then 
                  do k4=n1,n2 
                    nph4 = nph3*nphase(k4) 
                    HY4(k1,k2,k3,k4) =   nph4*GY4(-k1,-k2,-k3,-k4) 
                    Hi4(k1,k2,k3,k4) = - nph4*Gi4(-k1,-k2,-k3,-k4) 
                    H4(k1,k2,k3,k4)  = 
     $                    dcmplx(HY4(k1,k2,k3,k4),Hi4(k1,k2,k3,k4)*pi) 
                  enddo 
                endif 
              enddo 
            endif 
          enddo 
        endif 
      enddo 
      if (n1.eq.0) return
** correct the ill-defined entries
      HY2(-1,0) = - HY2(0,-1) 
      Hi2(-1,0) = Hi1(0)*HY1(-1)
      H2(-1,0) = dcmplx(HY2(-1,0),Hi2(-1,0)*pi)
      if ( nw.eq.2 ) return
      HY3(-1,0,0) = HY1(-1)*HY2(0,0)+HY3(0,0,-1) 
      Hi3(-1,0,0) = HY1(-1)*Hi2(0,0)-HY2(0,-1)*Hi1(0)
      H3(-1,0,0) = dcmplx(HY3(-1,0,0),Hi3(-1,0,0)*pi)
      if ( nw.eq.3 ) return
      HY4(-1,0,0,0) = -HY2(0,-1)*HY2(0,0)-HY4(0,0,0,-1) 
      Hi4(-1,0,0,0) = HY1(-1)*Hi3(0,0,0)+HY3(0,0,-1)*Hi1(0)
      H4(-1,0,0,0) = dcmplx(HY4(-1,0,0,0),Hi4(-1,0,0,0)*pi)
      return 
      end 
************************************************************************ 
      subroutine eval1dhplatm1(y,nw,H1,H2,H3,H4, 
     $                          HY1,HY2,HY3,HY4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
** evaluates 1dhpl's in the (-1)-range  -(r2+1) < y <= -(r2-1) 
** evaluating first the H(..,-y) by calling eval1dhplat1(-y), 
** and then expressing H(..,y) in terms of H(..,-y) 
      implicit double precision (a-h,o-z) 
      complex*16 H1,H2,H3,H4  
      complex*16 G1,G2,G3,G4  
      dimension H1(n1:n2),H2(n1:n2,n1:n2),H3(n1:n2,n1:n2,n1:n2), 
     $          H4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension HY1(n1:n2),HY2(n1:n2,n1:n2),HY3(n1:n2,n1:n2,n1:n2), 
     $          HY4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2) 
** additional arrays required within this routine 
      dimension G1(-n2:-n1),G2(-n2:-n1,-n2:-n1),
     $          G3(-n2:-n1,-n2:-n1,-n2:-n1), 
     $          G4(-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1) 
      dimension GY1(-n2:-n1),GY2(-n2:-n1,-n2:-n1),
     $          GY3(-n2:-n1,-n2:-n1,-n2:-n1), 
     $          GY4(-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1) 
      dimension Gi1(-n2:-n1),Gi2(-n2:-n1,-n2:-n1),
     $          Gi3(-n2:-n1,-n2:-n1,-n2:-n1), 
     $          Gi4(-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1) 
** 
      dimension istorfill(3)
      dimension nphase(-1:1) 
      data nphase/-1,1,-1/ 
      parameter (pi   = 3.14159265358979324d0) 
*      print*,' eval1dhplatm1: y = ',y 
** evaluate H(...,-y) 
      call setzero(nw,Gi1,Gi2,Gi3,Gi4,-n2,-n1) 
      Gi1(0) = -1
      call eval1dhplat1(-y,nw,G1,G2,G3,G4, 
     $                        GY1,GY2,GY3,GY4,Gi1,Gi2,Gi3,Gi4,-n2,-n1) 
** fill the arrays H's 
      do k1=n1,n2 
        nph1 = nphase(k1) 
        HY1(k1) =   nph1*GY1(-k1) 
        Hi1(k1) = - nph1*Gi1(-k1) 
        H1(k1)  = dcmplx(HY1(k1),Hi1(k1)*pi) 
        if ( nw.gt.1 ) then 
          do k2=n1,n2 
            nph2 = nph1*nphase(k2) 
            HY2(k1,k2) =   nph2*GY2(-k1,-k2) 
            Hi2(k1,k2) = - nph2*Gi2(-k1,-k2) 
            H2(k1,k2)  = dcmplx(HY2(k1,k2),Hi2(k1,k2)*pi) 
            if ( nw.gt.2 ) then 
              do k3=n1,n2 
                nph3 = nph2*nphase(k3) 
                HY3(k1,k2,k3) =   nph3*GY3(-k1,-k2,-k3) 
                Hi3(k1,k2,k3) = - nph3*Gi3(-k1,-k2,-k3) 
                H3(k1,k2,k3)  = dcmplx(HY3(k1,k2,k3),Hi3(k1,k2,k3)*pi) 
                if ( nw.gt.3 ) then 
                  do k4=n1,n2 
                    nph4 = nph3*nphase(k4) 
                    HY4(k1,k2,k3,k4) =   nph4*GY4(-k1,-k2,-k3,-k4) 
                    Hi4(k1,k2,k3,k4) = - nph4*Gi4(-k1,-k2,-k3,-k4) 
                    H4(k1,k2,k3,k4)  = 
     $                    dcmplx(HY4(k1,k2,k3,k4),Hi4(k1,k2,k3,k4)*pi) 
                  enddo 
                endif 
              enddo 
            endif 
          enddo 
        endif 
      enddo 
      return 
      end 


      subroutine eval1dhplatminf(y,nw,H1,H2,H3,H4, 
     $                          HY1,HY2,HY3,HY4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
** evaluates 1dhpl's in the (-1)-range y  <= -(r2+1) 
** evaluating first the H(..,-y) by calling eval1dhplatinf(-y), 
** and then expressing H(..,y) in terms of H(..,-y) 
      implicit double precision (a-h,o-z) 
      complex*16 H1,H2,H3,H4  
      complex*16 G1,G2,G3,G4  
      dimension H1(n1:n2),H2(n1:n2,n1:n2),H3(n1:n2,n1:n2,n1:n2), 
     $          H4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension HY1(n1:n2),HY2(n1:n2,n1:n2),HY3(n1:n2,n1:n2,n1:n2), 
     $          HY4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2) 
** additional arrays required within this routine 
      dimension G1(-n2:-n1),G2(-n2:-n1,-n2:-n1),
     $          G3(-n2:-n1,-n2:-n1,-n2:-n1), 
     $          G4(-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1) 
      dimension GY1(-n2:-n1),GY2(-n2:-n1,-n2:-n1),
     $          GY3(-n2:-n1,-n2:-n1,-n2:-n1), 
     $          GY4(-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1) 
      dimension Gi1(-n2:-n1),Gi2(-n2:-n1,-n2:-n1),
     $          Gi3(-n2:-n1,-n2:-n1,-n2:-n1), 
     $          Gi4(-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1) 
** 
      dimension istorfill(3)
      dimension nphase(-1:1) 
      data nphase/-1,1,-1/ 
      parameter (pi   = 3.14159265358979324d0) 
*      print*,' eval1dhplatm1: y = ',y 
** evaluate H(...,-y) 
      call setzero(nw,Gi1,Gi2,Gi3,Gi4,-n2,-n1) 
      Gi1(0) = -1
      call eval1dhplatinf(-y,nw,G1,G2,G3,G4, 
     $                        GY1,GY2,GY3,GY4,Gi1,Gi2,Gi3,Gi4,-n2,-n1) 
** fill the arrays H's 
      do k1=n1,n2 
        nph1 = nphase(k1) 
        HY1(k1) =   nph1*GY1(-k1) 
        Hi1(k1) = - nph1*Gi1(-k1) 
        H1(k1)  = dcmplx(HY1(k1),Hi1(k1)*pi) 
        if ( nw.gt.1 ) then 
          do k2=n1,n2 
            nph2 = nph1*nphase(k2) 
            HY2(k1,k2) =   nph2*GY2(-k1,-k2) 
            Hi2(k1,k2) = - nph2*Gi2(-k1,-k2) 
            H2(k1,k2)  = dcmplx(HY2(k1,k2),Hi2(k1,k2)*pi) 
            if ( nw.gt.2 ) then 
              do k3=n1,n2 
                nph3 = nph2*nphase(k3) 
                HY3(k1,k2,k3) =   nph3*GY3(-k1,-k2,-k3) 
                Hi3(k1,k2,k3) = - nph3*Gi3(-k1,-k2,-k3) 
                H3(k1,k2,k3)  = dcmplx(HY3(k1,k2,k3),Hi3(k1,k2,k3)*pi) 
                if ( nw.gt.3 ) then 
                  do k4=n1,n2 
                    nph4 = nph3*nphase(k4) 
                    HY4(k1,k2,k3,k4) =   nph4*GY4(-k1,-k2,-k3,-k4) 
                    Hi4(k1,k2,k3,k4) = - nph4*Gi4(-k1,-k2,-k3,-k4) 
                    H4(k1,k2,k3,k4)  = 
     $                    dcmplx(HY4(k1,k2,k3,k4),Hi4(k1,k2,k3,k4)*pi) 
                  enddo 
                endif 
              enddo 
            endif 
          enddo 
        endif 
      enddo 
      return 
      end 
************************************************************************ 
      subroutine setzero(nw,Hi1,Hi2,Hi3,Hi4,n1,n2) 
** initializes with 0 the elements of the arrays 
      implicit double precision (a-h,o-z) 
      dimension Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2) 
      do k1=n1,n2 
        Hi1(k1) = 0.d0 
        if ( nw.gt.1 ) then 
          do k2=n1,n2 
            Hi2(k1,k2) = 0.d0 
            if ( nw.gt.2 ) then 
              do k3=n1,n2 
                Hi3(k1,k2,k3) = 0.d0 
                if ( nw.gt.3 ) then 
                  do k4=n1,n2 
                    Hi4(k1,k2,k3,k4) = 0.d0 
                  enddo 
                endif 
              enddo 
            endif 
          enddo 
        endif 
      enddo 
      return 
      end 
************************************************************************ 
      subroutine fillred1dhpl(nw,H1,H2,H3,H4, 
     $                     HY1,HY2,HY3,HY4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
* fills the reducible 1dhpl from the irreducible set
      implicit double precision (a-h,o-z) 
      complex*16 H1,H2,H3,H4 
      dimension H1(n1:n2),H2(n1:n2,n1:n2),H3(n1:n2,n1:n2,n1:n2), 
     $          H4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension HY1(n1:n2),HY2(n1:n2,n1:n2),HY3(n1:n2,n1:n2,n1:n2), 
     $          HY4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2) 
      parameter (pinv = 0.318309886183790672d0) 
      parameter (pi   = 3.14159265358979324d0) 
      dimension infill(3)
** combining real and immaginary into the complex value 
      do k1=n1,n2 
      do k2=n1,n2 
        H2(k1,k2) = dcmplx(HY2(k1,k2),Hi2(k1,k2)*pi) 
        if ( nw.gt.2 ) then 
          do k3=n1,n2 
            H3(k1,k2,k3) = dcmplx(HY3(k1,k2,k3),Hi3(k1,k2,k3)*pi) 
            if ( nw.gt.3 ) then 
              do k4=n1,n2 
                H4(k1,k2,k3,k4) = 
     $               dcmplx(HY4(k1,k2,k3,k4),Hi4(k1,k2,k3,k4)*pi) 
              enddo 
            endif 
          enddo 
        endif 
      enddo 
      enddo 
** evaluating the reduced HPL's 
** iflag = 0 to suppress auxiliary printings of FILLREDHPLx 
      infill(1) =  0
      infill(2) = -1
      infill(3) =  1
      iflag = 0 
      do ia =  1,3 
      do ib = ia,3 
        call FILLREDHPL2(iflag,H1,H2,n1,n2,infill(ia),infill(ib)) 
        if ( nw.gt.2 ) then 
          do ic = ib,3 
            call FILLREDHPL3(iflag,H1,H2,H3,n1,n2, 
     $                          infill(ia),infill(ib),infill(ic)) 
            if ( nw.gt.3 ) then 
              do id = ic,3 
                call FILLREDHPL4(iflag,H1,H2,H3,H4,n1,n2, 
     $               infill(ia),infill(ib),infill(ic),infill(id)) 
              enddo 
            endif 
          enddo 
        endif 
      enddo 
      enddo 
** extractin real and immaginary parts from the complex value 
      do k1=n1,n2 
      do k2=n1,n2 
        HY2(k1,k2) =  dble(H2(k1,k2)) 
        Hi2(k1,k2) = aimag(H2(k1,k2))*pinv 
        if ( nw.gt.2 ) then 
          do k3=n1,n2 
            HY3(k1,k2,k3) =  dble(H3(k1,k2,k3)) 
            Hi3(k1,k2,k3) = aimag(H3(k1,k2,k3))*pinv 
            if ( nw.gt.3 ) then 
              do k4=n1,n2 
                HY4(k1,k2,k3,k4) =  dble(H4(k1,k2,k3,k4)) 
                Hi4(k1,k2,k3,k4) = aimag(H4(k1,k2,k3,k4))*pinv 
              enddo 
            endif 
          enddo 
        endif 
      enddo 
      enddo 
      return 
      end 
************************************************************************ 
      subroutine FILLREDHPL2(iflag,H1,H2,i1,i2,na,nb) 
      implicit double precision (a-h,o-z) 
      complex*16 H1,H2 
      dimension H1(i1:i2),H2(i1:i2,i1:i2) 
*23456789012345678901234567890123456789012345678901234567890123456789012 
* must be called with ordered indices na <= nb 
*      print*,' FILLREDHPL2, iflag =',iflag 
      if ( na.eq.nb ) then 
        H2(na,na) = 1.d0/2*( H1(na) )**2 
      else 
        H2(nb,na) = + H1(na)*H1(nb) - H2(na,nb) 
        if ( iflag.eq.1 ) then 
          call printer2(na,nb) 
        endif 
      endif 
      return 
      end 
************************************************************************ 
      subroutine FILLREDHPL3(iflag,H1,H2,H3,i1,i2,ia,ib,ic) 
      implicit double precision (a-h,o-z) 
      complex*16 H1,H2,H3 
      dimension H1(i1:i2),H2(i1:i2,i1:i2),H3(i1:i2,i1:i2,i1:i2) 
* must be called with "properly ordered" indices 
* note in particular the remapping of, say, (ia,ib,ic) into 
* (na,na,nb) of ReducerTest.out 
      na = ia 
      if ( (ia.eq.ib).and.(ib.eq.ic) ) then 
* case (na,na,na) 
        H3(na,na,na) = 1.d0/6*( H1(na) )**3 
* ic cannot be anymore equal to ia 
      else if ( ic.eq.ia ) then 
        print*,' FILLREDHPL3, error 1, called with arguments ' 
        print*,'               ',ia,ib,ic 
        stop 
      else if ( ia.eq.ib ) then 
* case (na,na,nb) 
        nb = ic 
        if ( iflag.eq.1 ) then 
          call printer3(na,na,nb) 
        endif 
        H3(na,nb,na) = + H1(na)*H2(na,nb) - 2*H3(na,na,nb) 
        H3(nb,na,na) = + 1.d0/2*H1(na)*H1(na)*H1(nb) 
     $                 - H1(na)*H2(na,nb) + H3(na,na,nb) 
* ib cannot be anymore equal to ia 
      else if ( ib.eq.ia ) then 
        print*,' FILLREDHPL3, error 2, called with arguments ' 
        print*,'               ',ia,ib,ic 
        stop 
      else if ( ib.eq.ic ) then 
* case (na,nb,nb) 
        nb = ib 
        if ( iflag.eq.1 ) then 
          call printer3(na,nb,nb) 
        endif 
        H3(nb,na,nb) = + H1(nb)*H2(na,nb) - 2*H3(na,nb,nb) 
        H3(nb,nb,na) = + 1.d0/2*H1(na)*H1(nb)*H1(nb) 
     $                 - H1(nb)*H2(na,nb) + H3(na,nb,nb) 
* no need to protect against ic.eq.ib 
* when arriving here all indices are different 
      else 
* case (na,nb,nc)    all indices are different 
        nb = ib 
        nc = ic 
        if ( iflag.eq.1 ) then 
          call printer3(na,nb,nc) 
          call printer3(na,nc,nb) 
        endif 
        H3(nb,na,nc) = + H1(nb)*H2(na,nc) 
     $                 - H3(na,nb,nc) - H3(na,nc,nb) 
        H3(nb,nc,na) = + H1(na)*H2(nb,nc) 
     $                 - H1(nb)*H2(na,nc) + H3(na,nc,nb) 
        H3(nc,na,nb) = + H1(nc)*H2(na,nb) 
     $                 - H3(na,nb,nc) - H3(na,nc,nb) 
        H3(nc,nb,na) = + H1(na)*H1(nb)*H1(nc) - H1(na)*H2(nb,nc) 
     $                 - H1(nc)*H2(na,nb) + H3(na,nb,nc) 
      endif 
*23456789012345678901234567890123456789012345678901234567890123456789012 
      return 
      end 
************************************************************************ 
      subroutine FILLREDHPL4(iflag,H1,H2,H3,H4,i1,i2,ia,ib,ic,id) 
      implicit double precision (a-h,o-z) 
      complex*16 H1,H2,H3,H4 
      dimension H1(i1:i2),H2(i1:i2,i1:i2),H3(i1:i2,i1:i2,i1:i2) 
      dimension H4(i1:i2,i1:i2,i1:i2,i1:i2) 
*23456789012345678901234567890123456789012345678901234567890123456789012 
* must be called with "properly ordered" indices 
* note in particular the remapping of, say, (ia,ib,ic) into 
* (na,na,nb) of ReducerTest.out 
      na = ia 
      if ( (ia.eq.ib).and.(ib.eq.ic).and.(ic.eq.id) ) then 
* case (na,na,na,na) 
        H4(na,na,na,na) = 1.d0/24*( H1(na) )**4 
* id cannot be anymore equal to ia 
      else if ( id.eq.ia ) then 
        print*,' FILLREDHPL4, error 1, called with arguments ' 
        print*,'               ',ia,ib,ic,id 
        stop 
      else if ( (ia.eq.ib).and.(ib.eq.ic) ) then 
* case (na,na,na,nb) 
        nb = id 
        H4(na,na,nb,na) = + H1(na)*H3(na,na,nb) - 3*H4(na,na,na,nb) 
        H4(na,nb,na,na) = + 1.d0/2*H1(na)*H1(na)*H2(na,nb) 
     $                    - 2*H1(na)*H3(na,na,nb) + 3*H4(na,na,na,nb) 
        H4(nb,na,na,na) = + 1.d0/6*H1(na)*H1(na)*H1(na)*H1(nb) 
     $                    - 1.d0/2*H1(na)*H1(na)*H2(na,nb) 
     $                    + H1(na)*H3(na,na,nb) - H4(na,na,na,nb) 
        if ( iflag.eq.1 ) then 
          call printer4(na,na,na,nb) 
        endif 
* ic cannot be anymore equal to ia 
      else if ( ic.eq.ia ) then 
        print*,' FILLREDHPL4, error 2, called with arguments ' 
        print*,'               ',ia,ib,ic,id 
        stop 
      else if ( (ia.eq.ib).and.(ic.eq.id) ) then 
* case (na,na,nb,nb) 
        nb = ic 
        H4(na,nb,na,nb) = + 1.d0/2*H2(na,nb)*H2(na,nb) 
     $                    - 2*H4(na,na,nb,nb) 
        H4(na,nb,nb,na) = + H1(na)*H3(na,nb,nb) 
     $                    - 1.d0/2*H2(na,nb)*H2(na,nb) 
        H4(nb,na,na,nb) = + H1(nb)*H3(na,na,nb) 
     $                    - 1.d0/2*H2(na,nb)*H2(na,nb) 
        H4(nb,na,nb,na) = + H1(na)*H1(nb)*H2(na,nb) 
     $                    - 2*H1(na)*H3(na,nb,nb) 
     $                    - 2*H1(nb)*H3(na,na,nb) 
     $                    + 1.d0/2*H2(na,nb)*H2(na,nb) 
     $                    + 2*H4(na,na,nb,nb) 
        H4(nb,nb,na,na) = + 1.d0/4*H1(na)*H1(na)*H1(nb)*H1(nb) 
     $                    - H1(na)*H1(nb)*H2(na,nb) 
     $                    + H1(na)*H3(na,nb,nb) 
     $                    + H1(nb)*H3(na,na,nb) - H4(na,na,nb,nb) 
        if ( iflag.eq.1 ) then 
          call printer4(na,na,nb,nb) 
        endif 
      else if ( ia.eq.ib ) then 
* case (na,na,nb,nc) 
        nb = ic 
        nc = id 
        H4(na,nb,nc,na) = + H1(na)*H3(na,nb,nc) - 2*H4(na,na,nb,nc) 
     $                    - H4(na,nb,na,nc) 
        H4(na,nc,na,nb) = + H2(na,nb)*H2(na,nc) - 2*H4(na,na,nb,nc) 
     $                    - 2*H4(na,na,nc,nb) - H4(na,nb,na,nc) 
        H4(na,nc,nb,na) = + H1(na)*H3(na,nc,nb) - H2(na,nb)*H2(na,nc) 
     $                    + 2*H4(na,na,nb,nc) + H4(na,nb,na,nc) 
        H4(nb,na,na,nc) = + H1(nb)*H3(na,na,nc) - H4(na,na,nb,nc) 
     $                    - H4(na,na,nc,nb) - H4(na,nb,na,nc) 
        H4(nb,na,nc,na) = + H1(na)*H1(nb)*H2(na,nc) 
     $                    - H1(na)*H3(na,nb,nc) - H1(na)*H3(na,nc,nb) 
     $                    - 2*H1(nb)*H3(na,na,nc) + 2*H4(na,na,nb,nc) 
     $                    + 2*H4(na,na,nc,nb) + H4(na,nb,na,nc) 
        H4(nb,nc,na,na) = + 1.d0/2*H1(na)*H1(na)*H2(nb,nc) 
     $                    - H1(na)*H1(nb)*H2(na,nc) 
     $                    + H1(na)*H3(na,nc,nb) + H1(nb)*H3(na,na,nc) 
     $                    - H4(na,na,nc,nb) 
        H4(nc,na,na,nb) = + H1(nc)*H3(na,na,nb) - H2(na,nb)*H2(na,nc) 
     $                    + H4(na,na,nb,nc) + H4(na,na,nc,nb) 
     $                    + H4(na,nb,na,nc) 
        H4(nc,na,nb,na) = + H1(na)*H1(nc)*H2(na,nb) 
     $                    - H1(na)*H3(na,nb,nc) - H1(na)*H3(na,nc,nb) 
     $                    - 2*H1(nc)*H3(na,na,nb) + H2(na,nb)*H2(na,nc) 
     $                    - H4(na,nb,na,nc) 
        H4(nc,nb,na,na) = + 1.d0/2*H1(na)*H1(na)*H1(nb)*H1(nc) 
     $                    - 1.d0/2*H1(na)*H1(na)*H2(nb,nc) 
     $                    - H1(na)*H1(nc)*H2(na,nb) 
     $                    + H1(na)*H3(na,nb,nc) + H1(nc)*H3(na,na,nb) 
     $                    - H4(na,na,nb,nc)  
        if ( iflag.eq.1 ) then 
          call printer4(na,na,nb,nc) 
          call printer4(na,na,nc,nb) 
          call printer4(na,nb,na,nc) 
        endif 
* ib cannot be anymore equal to ia 
      else if ( ib.eq.ia ) then 
        print*,' FILLREDHPL4, error 3, called with arguments ' 
        print*,'               ',ia,ib,ic,id 
        stop 
      else if ( (ib.eq.ic).and.(ic.eq.id) ) then 
* case (na,nb,nb,nb) 
        nb = ib 
        H4(nb,na,nb,nb) = + H1(nb)*H3(na,nb,nb) - 3*H4(na,nb,nb,nb) 
        H4(nb,nb,na,nb) = + 1.d0/2*H1(nb)*H1(nb)*H2(na,nb) 
     $                    - 2*H1(nb)*H3(na,nb,nb) + 3*H4(na,nb,nb,nb) 
        H4(nb,nb,nb,na) = + 1.d0/6*H1(na)*H1(nb)*H1(nb)*H1(nb) 
     $                    - 1.d0/2*H1(nb)*H1(nb)*H2(na,nb) 
     $                    + H1(nb)*H3(na,nb,nb) - H4(na,nb,nb,nb) 
        if ( iflag.eq.1 ) then 
          call printer4(na,nb,nb,nb) 
        endif 
* id cannot be anymore equal to ib 
      else if ( id.eq.ib ) then 
        print*,' FILLREDHPL4, error 4, called with arguments ' 
        print*,'               ',ia,ib,ic,id 
        stop 
      else if ( ib.eq.ic ) then 
* case (na,nb,nb,nc) 
        nb = ib 
        nc = id 
        H4(nb,na,nb,nc) = + H1(nb)*H3(na,nb,nc) 
     $                    - 2*H4(na,nb,nb,nc) - H4(na,nb,nc,nb) 
        H4(nb,na,nc,nb) = + H1(nb)*H3(na,nc,nb) - H4(na,nb,nc,nb) 
     $                    - 2*H4(na,nc,nb,nb) 
        H4(nb,nb,na,nc) = + 1.d0/2*H1(nb)*H1(nb)*H2(na,nc) 
     $                    - H1(nb)*H3(na,nb,nc) - H1(nb)*H3(na,nc,nb) 
     $                    + H4(na,nb,nb,nc) + H4(na,nb,nc,nb) 
     $                    + H4(na,nc,nb,nb) 
        H4(nb,nb,nc,na) = + H1(na)*H3(nb,nb,nc) 
     $                    - 1.d0/2*H1(nb)*H1(nb)*H2(na,nc) 
     $                    + H1(nb)*H3(na,nc,nb) - H4(na,nc,nb,nb) 
        H4(nb,nc,na,nb) = - H1(nb)*H3(na,nb,nc) - H1(nb)*H3(na,nc,nb) 
     $                    + H2(na,nb)*H2(nb,nc) + H4(na,nb,nc,nb) 
     $                    + 2*H4(na,nc,nb,nb) 
        H4(nb,nc,nb,na) = + H1(na)*H1(nb)*H2(nb,nc) 
     $                    - 2*H1(na)*H3(nb,nb,nc) 
     $                    + H1(nb)*H3(na,nb,nc) 
     $                    - H2(na,nb)*H2(nb,nc) - H4(na,nb,nc,nb) 
        H4(nc,na,nb,nb) = + H1(nc)*H3(na,nb,nb) - H4(na,nb,nb,nc) 
     $                    - H4(na,nb,nc,nb) - H4(na,nc,nb,nb) 
        H4(nc,nb,na,nb) = + H1(nb)*H1(nc)*H2(na,nb) 
     $                    - 2*H1(nc)*H3(na,nb,nb) 
     $                    - H2(na,nb)*H2(nb,nc) + 2*H4(na,nb,nb,nc) 
     $                    + H4(na,nb,nc,nb) 
        H4(nc,nb,nb,na) = + 1.d0/2*H1(na)*H1(nb)*H1(nb)*H1(nc) 
     $                    - H1(na)*H1(nb)*H2(nb,nc) 
     $                    + H1(na)*H3(nb,nb,nc) 
     $                    - H1(nb)*H1(nc)*H2(na,nb) 
     $                    + H1(nc)*H3(na,nb,nb) + H2(na,nb)*H2(nb,nc) 
     $                    - H4(na,nb,nb,nc) 
        if ( iflag.eq.1 ) then 
          call printer4(na,nb,nb,nc) 
          call printer4(na,nb,nc,nb) 
          call printer4(na,nc,nb,nb) 
        endif 
* ic cannot be anymore equal to ib 
      else if ( ic.eq.ib ) then 
        print*,' FILLREDHPL4, error 5, called with arguments ' 
        print*,'               ',ia,ib,ic,id 
        stop 
      else if ( ic.eq.id ) then 
* case (na,nb,nc,nc) 
        nb = ib 
        nc = ic 
        H4(nb,na,nc,nc) = + H1(nb)*H3(na,nc,nc) - H4(na,nb,nc,nc) 
     $                    - H4(na,nc,nb,nc) - H4(na,nc,nc,nb) 
        H4(nb,nc,na,nc) = - 2*H1(nb)*H3(na,nc,nc) + H2(na,nc)*H2(nb,nc) 
     $                    + H4(na,nc,nb,nc) + 2*H4(na,nc,nc,nb) 
        H4(nb,nc,nc,na) = + H1(na)*H3(nb,nc,nc) + H1(nb)*H3(na,nc,nc) 
     $                    - H2(na,nc)*H2(nb,nc) - H4(na,nc,nc,nb) 
        H4(nc,na,nb,nc) = + H1(nc)*H3(na,nb,nc) - 2*H4(na,nb,nc,nc) 
     $                    - H4(na,nc,nb,nc) 
        H4(nc,na,nc,nb) = + H1(nc)*H3(na,nc,nb) - H4(na,nc,nb,nc) 
     $                    - 2*H4(na,nc,nc,nb) 
        H4(nc,nb,na,nc) = + H1(nb)*H1(nc)*H2(na,nc) 
     $                    - H1(nc)*H3(na,nb,nc) - H1(nc)*H3(na,nc,nb) 
     $                    - H2(na,nc)*H2(nb,nc) + 2*H4(na,nb,nc,nc) 
     $                    + H4(na,nc,nb,nc) 
        H4(nc,nb,nc,na) = + H1(na)*H1(nc)*H2(nb,nc) 
     $                    - 2*H1(na)*H3(nb,nc,nc) 
     $                    - H1(nb)*H1(nc)*H2(na,nc) 
     $                    + H1(nc)*H3(na,nc,nb) + H2(na,nc)*H2(nb,nc) 
     $                    - H4(na,nc,nb,nc) 
        H4(nc,nc,na,nb) = + 1.d0/2*H1(nc)*H1(nc)*H2(na,nb) 
     $                    - H1(nc)*H3(na,nb,nc) - H1(nc)*H3(na,nc,nb) 
     $                    + H4(na,nb,nc,nc) + H4(na,nc,nb,nc) 
     $                    + H4(na,nc,nc,nb) 
        H4(nc,nc,nb,na) = + 1.d0/2*H1(na)*H1(nb)*H1(nc)*H1(nc) 
     $                    - H1(na)*H1(nc)*H2(nb,nc) 
     $                    + H1(na)*H3(nb,nc,nc) 
     $                    - 1.d0/2*H1(nc)*H1(nc)*H2(na,nb) 
     $                    + H1(nc)*H3(na,nb,nc) - H4(na,nb,nc,nc) 
        if ( iflag.eq.1 ) then 
          call printer4(na,nb,nc,nc) 
          call printer4(na,nc,nb,nc) 
          call printer4(na,nc,nc,nb) 
        endif 
* no need to protect against id.eq.ic 
* when arriving here all indices are different 
      else 
* case (na,nb,nc,nd) all indices are different 
        nb = ib 
        nc = ic 
        nd = id 
        H4(nb,na,nc,nd) = + H1(nb)*H3(na,nc,nd) - H4(na,nb,nc,nd) 
     $                    - H4(na,nc,nb,nd) - H4(na,nc,nd,nb) 
        H4(nb,na,nd,nc) = + H1(nb)*H3(na,nd,nc) - H4(na,nb,nd,nc) 
     $                    - H4(na,nd,nb,nc) - H4(na,nd,nc,nb) 
        H4(nb,nc,na,nd) = - H1(nb)*H3(na,nc,nd) - H1(nb)*H3(na,nd,nc) 
     $                    + H2(na,nd)*H2(nb,nc) + H4(na,nc,nb,nd) 
     $                    + H4(na,nc,nd,nb) + H4(na,nd,nc,nb) 
        H4(nb,nc,nd,na) = + H1(na)*H3(nb,nc,nd) + H1(nb)*H3(na,nd,nc) 
     $                    - H2(na,nd)*H2(nb,nc) - H4(na,nd,nc,nb) 
        H4(nb,nd,na,nc) = - H1(nb)*H3(na,nc,nd) - H1(nb)*H3(na,nd,nc) 
     $                    + H2(na,nc)*H2(nb,nd) + H4(na,nc,nd,nb) 
     $                    + H4(na,nd,nb,nc) + H4(na,nd,nc,nb) 
        H4(nb,nd,nc,na) = + H1(na)*H3(nb,nd,nc) + H1(nb)*H3(na,nc,nd) 
     $                    - H2(na,nc)*H2(nb,nd) - H4(na,nc,nd,nb) 
        H4(nc,na,nb,nd) = + H1(nc)*H3(na,nb,nd) - H4(na,nb,nc,nd) 
     $                    - H4(na,nb,nd,nc) - H4(na,nc,nb,nd) 
        H4(nc,na,nd,nb) = + H1(nc)*H3(na,nd,nb) - H4(na,nc,nd,nb) 
     $                    - H4(na,nd,nb,nc) - H4(na,nd,nc,nb) 
        H4(nc,nb,na,nd) = + H1(nb)*H1(nc)*H2(na,nd) 
     $                    - H1(nc)*H3(na,nb,nd) - H1(nc)*H3(na,nd,nb) 
     $                    - H2(na,nd)*H2(nb,nc) + H4(na,nb,nc,nd) 
     $                    + H4(na,nb,nd,nc) + H4(na,nd,nb,nc) 
        H4(nc,nb,nd,na) = + H1(na)*H1(nc)*H2(nb,nd) 
     $                    - H1(na)*H3(nb,nc,nd) - H1(na)*H3(nb,nd,nc) 
     $                    - H1(nb)*H1(nc)*H2(na,nd) 
     $                    + H1(nc)*H3(na,nd,nb) + H2(na,nd)*H2(nb,nc) 
     $                    - H4(na,nd,nb,nc) 
        H4(nc,nd,na,nb) = - H1(nc)*H3(na,nb,nd) - H1(nc)*H3(na,nd,nb) 
     $                    + H2(na,nb)*H2(nc,nd) + H4(na,nb,nd,nc) 
     $                    + H4(na,nd,nb,nc) + H4(na,nd,nc,nb) 
        H4(nc,nd,nb,na) = + H1(na)*H1(nb)*H2(nc,nd) 
     $                    - H1(na)*H1(nc)*H2(nb,nd) 
     $                    + H1(na)*H3(nb,nd,nc) + H1(nc)*H3(na,nb,nd) 
     $                    - H2(na,nb)*H2(nc,nd) - H4(na,nb,nd,nc) 
        H4(nd,na,nb,nc) = + H1(nd)*H3(na,nb,nc) - H4(na,nb,nc,nd) 
     $                    - H4(na,nb,nd,nc) - H4(na,nd,nb,nc) 
        H4(nd,na,nc,nb) = + H1(nd)*H3(na,nc,nb) - H4(na,nc,nb,nd) 
     $                    - H4(na,nc,nd,nb) - H4(na,nd,nc,nb) 
        H4(nd,nb,na,nc) = + H1(nb)*H1(nd)*H2(na,nc) 
     $                    - H1(nd)*H3(na,nb,nc) - H1(nd)*H3(na,nc,nb) 
     $                    - H2(na,nc)*H2(nb,nd) + H4(na,nb,nc,nd) 
     $                    + H4(na,nb,nd,nc) + H4(na,nc,nb,nd) 
        H4(nd,nb,nc,na) = + H1(na)*H1(nd)*H2(nb,nc) 
     $                    - H1(na)*H3(nb,nc,nd) - H1(na)*H3(nb,nd,nc) 
     $                    - H1(nb)*H1(nd)*H2(na,nc) 
     $                    + H1(nd)*H3(na,nc,nb) + H2(na,nc)*H2(nb,nd) 
     $                    - H4(na,nc,nb,nd) 
        H4(nd,nc,na,nb) = + H1(nc)*H1(nd)*H2(na,nb) 
     $                    - H1(nd)*H3(na,nb,nc) - H1(nd)*H3(na,nc,nb) 
     $                    - H2(na,nb)*H2(nc,nd) + H4(na,nb,nc,nd) 
     $                    + H4(na,nc,nb,nd) + H4(na,nc,nd,nb) 
        H4(nd,nc,nb,na) = + H1(na)*H1(nb)*H1(nc)*H1(nd) 
     $                    - H1(na)*H1(nb)*H2(nc,nd) 
     $                    - H1(na)*H1(nd)*H2(nb,nc) 
     $                    + H1(na)*H3(nb,nc,nd) 
     $                    - H1(nc)*H1(nd)*H2(na,nb) 
     $                    + H1(nd)*H3(na,nb,nc) 
     $                    + H2(na,nb)*H2(nc,nd) - H4(na,nb,nc,nd) 
        if ( iflag.eq.1 ) then 
          call printer4(na,nb,nc,nd) 
          call printer4(na,nb,nd,nc) 
          call printer4(na,nc,nb,nd) 
          call printer4(na,nc,nb,nd) 
          call printer4(na,nd,nb,nc) 
          call printer4(na,nd,nc,nb) 
        endif 
      endif 
*23456789012345678901234567890123456789012345678901234567890123456789012 
      return 
      end 
************************************************************************ 
      subroutine printer2(na,nb) 

      write(11,'(''g [H('',$)') 
      call subprint(11,na) 
      write(11,'('','',$)') 
      call subprint(11,nb) 
      write(11,'('',y)] = H('',$)') 
      call subprint(11,na) 
      write(11,'('','',$)') 
      call subprint(11,nb) 
      write(11,'('',y) ; '')') 

      write(12,'(''id H('',$)') 
      call subprint(12,na) 
      write(12,'('','',$)') 
      call subprint(12,nb) 
      write(12,'('',y) = H[('',$)') 
      call subprint(12,na) 
      write(12,'('','',$)') 
      call subprint(12,nb) 
      write(12,'('',y)] ; '')') 

      return 
      end 
*** 
      subroutine printer3(na,nb,nc) 

      write(11,'(''g [H('',$)') 
      call subprint(11,na) 
      write(11,'('','',$)') 
      call subprint(11,nb) 
      write(11,'('','',$)') 
      call subprint(11,nc) 
      write(11,'('',y)] = H('',$)') 
      call subprint(11,na) 
      write(11,'('','',$)') 
      call subprint(11,nb) 
      write(11,'('','',$)') 
      call subprint(11,nc) 
      write(11,'('',y) ; '')') 

      write(12,'(''id H('',$)') 
      call subprint(12,na) 
      write(12,'('','',$)') 
      call subprint(12,nb) 
      write(12,'('','',$)') 
      call subprint(12,nc) 
      write(12,'('',y) = H[('',$)') 
      call subprint(12,na) 
      write(12,'('','',$)') 
      call subprint(12,nb) 
      write(12,'('',y)] ; '')') 

      return 
      end 
*** 
      subroutine printer4(na,nb,nc,nd) 

      write(11,'(''g [H('',$)') 
      call subprint(11,na) 
      write(11,'('','',$)') 
      call subprint(11,nb) 
      write(11,'('','',$)') 
      call subprint(11,nc) 
      write(11,'('','',$)') 
      call subprint(11,nd) 
      write(11,'('',y)] = H('',$)') 
      call subprint(11,na) 
      write(11,'('','',$)') 
      call subprint(11,nb) 
      write(11,'('','',$)') 
      call subprint(11,nc) 
      write(11,'('','',$)') 
      call subprint(11,nd) 
      write(11,'('',y) ; '')') 

      write(12,'(''id H('',$)') 
      call subprint(12,na) 
      write(12,'('','',$)') 
      call subprint(12,nb) 
      write(12,'('','',$)') 
      call subprint(12,nc) 
      write(12,'('','',$)') 
      call subprint(12,nd) 
      write(12,'('',y) = H[('',$)') 
      call subprint(12,na) 
      write(12,'('','',$)') 
      call subprint(12,nb) 
      write(12,'('','',$)') 
      call subprint(12,nc) 
      write(12,'('','',$)') 
      call subprint(12,nd) 
      write(12,'('',y)] ; '')') 

      return 
      end 
*** 
      subroutine subprint(n,na) 
      if ( na.lt.0 ) then 
        write (n,102) na 
      else 
        write (n,101) na 
      endif 
      return 
  101 format(i1,$) 
  102 format(i2,$) 
      end 

************************************************************************
** the following routines contain th set of routines evaluating 
** irreducible 1dhpl's for various values of the arguments 
************************************************************************ 
      subroutine fillh1(y,H1,HY1,Hi1,n1,n2) 
** fillh1 evaluates the 1dhpl's of weight 1 
      implicit double precision (a-h,o-z) 
      complex*16 H1 
      dimension H1(n1:n2) 
      dimension HY1(n1:n2) 
      dimension Hi1(n1:n2) 
      parameter (pi   = 3.14159265358979324d0) 
      if ( n1.eq.-1) then 
        if ( y.ge.-1.d0 ) then 
          HY1(-1) = log(1.d0+y) 
          Hi1(-1) = 0.d0 
        elseif ( y.lt.-1.d0 ) then 
          HY1(-1) = log(-1.d0-y) 
          Hi1(-1) = 1.d0 
        endif 
        H1(-1) = dcmplx(HY1(-1),pi*Hi1(-1)) 
      endif 
      if ( y.ge.0.d0 ) then 
        HY1(0) = log(y) 
*        Hi1(0) = 0.d0 
      elseif ( y.lt.0.d0 ) then 
        HY1(0) = log(-y) 
        Hi1(0) = 1.d0 
      endif 
      H1(0) = dcmplx(HY1(0),pi*Hi1(0)) 
      if ( n2.eq.1 ) then 
        if ( y.ge.1.d0 ) then 
          HY1(1) = - log(-1.d0+y) 
          Hi1(1) = 1.d0 
        elseif ( y.lt.1.d0 ) then 
          HY1(1) = - log(1.d0-y) 
          Hi1(1) = 0.d0 
        endif 
        H1(1) = dcmplx(HY1(1),pi*Hi1(1)) 
      endif 
      return 
      end 
************************************************************************ 
      subroutine fillirr1dhplat0(y,nw,HY1,HY2,HY3,HY4,n1,n2) 
** evaluate the HPL from their power series expansions
** fillirr1dhplat0 is called by eval1dhplat0; 
** it is guaranteed that nw is in the range 1:4, and that (n1,n2) 
** take one of the pairs of values (0,1), (-1,0) or (-1,1) 
** 
** for y < 0 DOES NOT evaluates the immaginary part of H(0,y) = log(y) 
      implicit double precision (a-h,o-z) 
      dimension HY1(n1:n2),HY2(n1:n2,n1:n2),HY3(n1:n2,n1:n2,n1:n2), 
     $          HY4(n1:n2,n1:n2,n1:n2,n1:n2) 
** evaluating the required 1dHPL of weight 1 
      if ( n1.eq.-1) then 
** 1+y = (1+ep)/(1-ep), ep = y/(2+y) 
** log(1+y) = log((1+y)/(1-y)) = 2*ep*(1+ep^2/3+ep^4/5+.....) 
** at y= -(r2-1) = - 0.4142135624, ep = - 0.26120387496 
** ep2 = 0.068227464296, ep2^13 = 6.9 x 10^(-16) 
         ep = y/(2.d0+y) 
         e2 = ep*ep 
*         v = log(1.d0+y) 
         v = 2*ep*(1+e2*(1.d0/ 3+e2*(1.d0/ 5+e2*(1.d0/ 7+e2*(1.d0/ 9 
     $              +e2*(1.d0/11+e2*(1.d0/13+e2*(1.d0/15+e2*(1.d0/17    
     $              +e2*(1.d0/19+e2*(1.d0/21+e2*(1.d0/23+e2*(1.d0/25    
     $              )))))))))))))
         HY1(-1) = v 
      endif 
      if (y.ge.0d0) then 
         HY1(0) = log(y) 
      else 
         HY1(0) = log(-y) 
** the immaginary part is evaluated in the calling routine eval1dhplat0 
**       Hi1(0) = 1d0 
      endif 
      if ( n2.eq.1) then 
** 1-y = (1-ep)/(1+ep), ep = y/(2-y) 
         ep = y/(2.d0-y) 
         e2 = ep*ep 
*         u = - log(1.d0-y) 
         u = 2*ep*(1+e2*(1.d0/ 3+e2*(1.d0/ 5+e2*(1.d0/ 7+e2*(1.d0/ 9 
     $              +e2*(1.d0/11+e2*(1.d0/13+e2*(1.d0/15+e2*(1.d0/17    
     $              +e2*(1.d0/19+e2*(1.d0/21+e2*(1.d0/23+e2*(1.d0/25    
     $              )))))))))))))
         HY1(1) = u 
      endif 
      if ( nw.eq.1 ) return 
** from now on nw > 1 
** evaluating the Cebyshev polynomials for the expansions 
      ep = y 
      if ( n2.eq.1) then 
        tu01 = 20d0/11d0*u 
        tu02 = 2d0*tu01*tu01 - 1d0 
        tu03 = 2d0*tu01*tu02 - tu01 
        tu04 = 2d0*tu01*tu03 - tu02 
        tu05 = 2d0*tu01*tu04 - tu03 
        tu06 = 2d0*tu01*tu05 - tu04 
        tu07 = 2d0*tu01*tu06 - tu05 
        tu08 = 2d0*tu01*tu07 - tu06 
        tu09 = 2d0*tu01*tu08 - tu07 
        tu10 = 2d0*tu01*tu09 - tu08 
        tu11 = 2d0*tu01*tu10 - tu09 
        tu12 = 2d0*tu01*tu11 - tu10 
      endif 
      if ( n1.eq.-1 ) then 
        tv01 = 20d0/11d0*v 
        tv02 = 2d0*tv01*tv01 - 1d0 
        tv03 = 2d0*tv01*tv02 - tv01 
        tv04 = 2d0*tv01*tv03 - tv02 
        tv05 = 2d0*tv01*tv04 - tv03 
        tv06 = 2d0*tv01*tv05 - tv04 
        tv07 = 2d0*tv01*tv06 - tv05 
        tv08 = 2d0*tv01*tv07 - tv06 
        tv09 = 2d0*tv01*tv08 - tv07 
        tv10 = 2d0*tv01*tv09 - tv08 
        tv11 = 2d0*tv01*tv10 - tv09 
        tv12 = 2d0*tv01*tv11 - tv10 
      endif 
** evaluating the expansions 
** (n1,n2) = (0,1) or (-1,1) 
      if (    ( (n1.eq.0).and.(n2.eq.1) ) 
     $    .or.( (n1.eq.-1).and.(n2.eq.1) ) ) then 
      HY2(0,1) = 
     $  - 3.781250000000000d-02 
     $  + 5.534574473824441d-01*tu01 
     $  - 3.781250000000000d-02*tu02 
     $  + 1.151036617760703d-03*tu03 
     $  - 8.659502433858922d-07*tu05 
     $  + 1.109042494804544d-09*tu07 
     $  - 1.624415058184216d-12*tu09 
     $  + 2.528376460336939d-15*tu11 
**    it would be wrong to write 
**    if ( nw.eq.2 ) return 
**    because the (n1.eq.-1).and.(n2.eq.1) case is not yet complete 
      if ( nw.gt.2 ) then 
      HY3(0,0,1) = 
     $  - 5.701592410758114d-02 
     $  + 5.598247957892565d-01*tu01 
     $  - 5.711486614505007d-02*tu02 
     $  + 3.275603992203700d-03*tu03 
     $  - 9.887255877938583d-05*tu04 
     $  + 4.021153684652295d-07*tu05 
     $  + 6.939288687864526d-08*tu06 
     $  - 7.995347631322020d-10*tu07 
     $  - 8.567978673919505d-11*tu08 
     $  + 1.526387027481200d-12*tu09 
     $  + 1.226899454816980d-13*tu10 
     $  - 2.848614761014972d-15*tu11 
     $  - 1.880542777479446d-16*tu12 
      HY3(0,1,1) = 
     $  + 3.816894981500984d-02 
     $  - 1.039843750000000d-02*tu01 
     $  + 3.828760080995617d-02*tu02 
     $  - 3.466145833333333d-03*tu03 
     $  + 1.185518160084905d-04*tu04 
     $  - 9.904555648775859d-08*tu06 
     $  + 1.331803984518588d-10*tu08 
     $  - 2.006389465106708d-13*tu10 
     $  + 3.180731062055677d-16*tu12 
      endif 
      if ( nw.gt.3 ) then 
      HY4(0,0,0,1) = 
     $  - 6.685228257646101d-02 
     $  + 5.645990701998083d-01*tu01 
     $  - 6.707912936340146d-02*tu02 
     $  + 4.876429488624746d-03*tu03 
     $  - 2.268732672568699d-04*tu04 
     $  + 6.038494106229146d-06*tu05 
     $  - 2.642577015932576d-08*tu06 
     $  - 3.679843316593900d-09*tu07 
     $  + 5.444046563879984d-11*tu08 
     $  + 4.063821221202881d-12*tu09 
     $  - 1.055985864474070d-13*tu10 
     $  - 5.190408125225683d-15*tu11 
     $  + 1.985464489219049d-16*tu12 
      HY4(0,0,1,1) = 
     $  + 1.953236111099851d-02 
     $  - 8.741612828671381d-03*tu01 
     $  + 1.974116110893196d-02*tu02 
     $  - 2.926558492394004d-03*tu03 
     $  + 2.088576190269387d-04*tu04 
     $  - 7.604351107741397d-06*tu05 
     $  + 5.751031394942524d-08*tu06 
     $  + 5.832253077603139d-09*tu07 
     $  - 1.105713721511985d-10*tu08 
     $  - 7.453416210082473d-12*tu09 
     $  + 2.077758906032370d-13*tu10 
     $  + 1.085601519719514d-14*tu11 
     $  - 3.848312092795918d-16*tu12 
      HY4(0,1,1,1) = 
     $  - 7.148925781250000d-04 
     $  + 7.019393481825299d-03*tu01 
     $  - 9.531901041666666d-04*tu02 
     $  + 2.354287493676137d-03*tu03 
     $  - 2.382975260416666d-04*tu04 
     $  + 8.682904829408987d-06*tu05 
     $  - 7.768198634676578d-09*tu07 
     $  + 1.083130072188330d-11*tu09 
     $  - 1.668810490326842d-14*tu11 
      endif 
** nw > 3 endif 
      endif 
** (n1,n2) = (0,1) or (-1,1) endif 
************ 
** (n1,n2) = (-1,0) or (-1,1) 
      if (    ( (n1.eq.-1).and.(n2.eq.0) ) 
     $    .or.( (n1.eq.-1).and.(n2.eq.1) ) ) then 
      HY2(0,-1) = 
     $  + 3.781250000000000d-02 
     $  + 5.534574473824441d-01*tv01 
     $  + 3.781250000000000d-02*tv02 
     $  + 1.151036617760703d-03*tv03 
     $  - 8.659502433858922d-07*tv05 
     $  + 1.109042494804544d-09*tv07 
     $  - 1.624415058184216d-12*tv09 
     $  + 2.528376460336939d-15*tv11 
      if ( nw.gt.2 ) then 
      HY3(0,0,-1) = 
     $  + 5.701592410758114d-02 
     $  + 5.598247957892565d-01*tv01 
     $  + 5.711486614505007d-02*tv02 
     $  + 3.275603992203700d-03*tv03 
     $  + 9.887255877938583d-05*tv04 
     $  + 4.021153684652295d-07*tv05 
     $  - 6.939288687864526d-08*tv06 
     $  - 7.995347631322020d-10*tv07 
     $  + 8.567978673919505d-11*tv08 
     $  + 1.526387027481200d-12*tv09 
     $  - 1.226899454816980d-13*tv10 
     $  - 2.848614761014972d-15*tv11 
     $  + 1.880542777479446d-16*tv12 
      HY3(0,-1,-1) = 
     $  + 3.816894981500984d-02 
     $  + 1.039843750000000d-02*tv01 
     $  + 3.828760080995617d-02*tv02 
     $  + 3.466145833333333d-03*tv03 
     $  + 1.185518160084905d-04*tv04 
     $  - 9.904555648775859d-08*tv06 
     $  + 1.331803984518588d-10*tv08 
     $  - 2.006389465106708d-13*tv10 
     $  + 3.180731062055677d-16*tv12 
      endif 
      if ( nw.gt.3 ) then 
      HY4(0,0,0,-1) = 
     $  + 6.685228257646101d-02 
     $  + 5.645990701998083d-01*tv01 
     $  + 6.707912936340146d-02*tv02 
     $  + 4.876429488624746d-03*tv03 
     $  + 2.268732672568699d-04*tv04 
     $  + 6.038494106229146d-06*tv05 
     $  + 2.642577015932576d-08*tv06 
     $  - 3.679843316593900d-09*tv07 
     $  - 5.444046563879984d-11*tv08 
     $  + 4.063821221202881d-12*tv09 
     $  + 1.055985864474070d-13*tv10 
     $  - 5.190408125225683d-15*tv11 
     $  - 1.985464489219049d-16*tv12 
      HY4(0,0,-1,-1) = 
     $  + 1.953236111099851d-02 
     $  + 8.741612828671381d-03*tv01 
     $  + 1.974116110893196d-02*tv02 
     $  + 2.926558492394004d-03*tv03 
     $  + 2.088576190269387d-04*tv04 
     $  + 7.604351107741397d-06*tv05 
     $  + 5.751031394942524d-08*tv06 
     $  - 5.832253077603139d-09*tv07 
     $  - 1.105713721511985d-10*tv08 
     $  + 7.453416210082473d-12*tv09 
     $  + 2.077758906032370d-13*tv10 
     $  - 1.085601519719514d-14*tv11 
     $  - 3.848312092795918d-16*tv12 
      HY4(0,-1,-1,-1) = 
     $  + 7.148925781250000d-04 
     $  + 7.019393481825299d-03*tv01 
     $  + 9.531901041666666d-04*tv02 
     $  + 2.354287493676137d-03*tv03 
     $  + 2.382975260416666d-04*tv04 
     $  + 8.682904829408987d-06*tv05 
     $  - 7.768198634676578d-09*tv07 
     $  + 1.083130072188330d-11*tv09 
     $  - 1.668810490326842d-14*tv11 
      endif 
** nw > 3 endif 
      endif 
** (n1,n2) = (-1,0) or (-1,1) endif 
** (n1,n2) = (-1,1) -- completion 
      if ( (n1.eq.-1).and.(n2.eq.1) ) then 
      HY2(-1,1) = 
     $  - 2.924454241163343d-02 
     $  + 3.845279287117326d-01*tu01 
     $  - 2.925485694830038d-02*tu02 
     $  + 1.097780471057338d-03*tu03 
     $  - 1.029703135442673d-05*tu04 
     $  - 7.265175511511970d-07*tu05 
     $  + 1.747461299829753d-08*tu06 
     $  + 7.707353556013722d-10*tu07 
     $  - 3.064611747990741d-11*tu08 
     $  - 8.531228176305706d-13*tu09 
     $  + 5.331187822989144d-14*tu10 
     $  + 8.500141365188675d-16*tu11 
     $  - 6.931471805599453d-01*HY1(-1) 
      if ( nw.gt.2 ) then 
      HY3(0,-1,1) = 
     $  - 4.107537580582269d-02 
     $  + 3.887609555197323d-01*tu01 
     $  - 4.116162793629221d-02*tu02 
     $  + 2.511526558054413d-03*tu03 
     $  - 8.620496933228561d-05*tu04 
     $  + 9.128023201466990d-07*tu05 
     $  + 4.711634663963971d-08*tu06 
     $  - 1.347359673414334d-09*tu07 
     $  - 4.474345520888852d-11*tu08 
     $  + 2.138249646727980d-12*tu09 
     $  + 4.709915818801180d-14*tu10 
     $  - 3.454431385666621d-15*tu11 
     $  - 6.931471805599453d-01*HY2(0,-1) 
      HY3(0,1,-1) = 
     $  - 4.107537580582269d-02 
     $  - 3.887609555197323d-01*tv01 
     $  - 4.116162793629221d-02*tv02 
     $  - 2.511526558054413d-03*tv03 
     $  - 8.620496933228561d-05*tv04 
     $  - 9.128023201466990d-07*tv05 
     $  + 4.711634663963971d-08*tv06 
     $  + 1.347359673414334d-09*tv07 
     $  - 4.474345520888852d-11*tv08 
     $  - 2.138249646727980d-12*tv09 
     $  + 4.709915818801180d-14*tv10 
     $  + 3.454431385666621d-15*tv11 
     $  + 6.931471805599453d-01*HY2(0,1) 
      HY3(-1,-1,1) = 
     $  - 3.590863871372201d-02 
     $  + 3.272029419300922d-01*tu01 
     $  - 3.599657175069328d-02*tu02 
     $  + 2.325685169395631d-03*tu03 
     $  - 8.788997314012583d-05*tu04 
     $  + 1.277831858501559d-06*tu05 
     $  + 4.303730428865162d-08*tu06 
     $  - 1.992295216809703d-09*tu07 
     $  - 2.652932076676834d-11*tu08 
     $  + 3.159865930142703d-12*tu09 
     $  - 2.395589527593406d-15*tu10 
     $  - 4.870947810519399d-15*tu11 
     $  - 5.822405264650125d-01*HY1(-1) 
     $  - 3.465735902799726d-01*HY1(-1)*HY1(-1) 
      HY3(-1,1,1) = 
     $  + 3.668493142404161d-02 
     $  - 1.413123104773291d-01*tu01 
     $  + 3.680167312678666d-02*tu02 
     $  - 3.064044728536094d-03*tu03 
     $  + 1.166524199994130d-04*tu04 
     $  - 8.779983417383380d-07*tu05 
     $  - 8.917940330502000d-08*tu06 
     $  + 1.787575622706040d-09*tu07 
     $  + 1.032182649980912d-10*tu08 
     $  - 3.441821872732193d-12*tu09 
     $  - 1.239218730863368d-13*tu10 
     $  + 6.355731482672869d-15*tu11 
     $  + 1.386175839607904d-16*tu12 
     $  + 2.402265069591007d-01*HY1(-1)       
      endif 
      if ( nw.gt.3 ) then 
      HY4(0,0,-1,1) = 
     $  - 4.713463351559199d-02 
     $  + 3.918037828258655d-01*tu01 
     $  - 4.730698763577787d-02*tu02 
     $  + 3.532784273601097d-03*tu03 
     $  - 1.724036773635937d-04*tu04 
     $  + 5.100573466380115d-06*tu05 
     $  - 4.948996960052575d-08*tu06 
     $  - 2.345390965359666d-09*tu07 
     $  + 6.710522628543514d-11*tu08 
     $  + 1.979867116023822d-12*tu09 
     $  - 1.027163441987459d-13*tu10 
     $  - 1.836436639605094d-15*tu11 
     $  + 1.633620651699784d-16*tu12 
     $  - 6.931471805599453d-01*HY3(0,0,-1) 
      HY4(0,0,1,-1) = 
     $  - 4.713463351559199d-02 
     $  - 3.918037828258655d-01*tv01 
     $  - 4.730698763577787d-02*tv02 
     $  - 3.532784273601097d-03*tv03 
     $  - 1.724036773635937d-04*tv04 
     $  - 5.100573466380115d-06*tv05 
     $  - 4.948996960052575d-08*tv06 
     $  + 2.345390965359666d-09*tv07 
     $  + 6.710522628543514d-11*tv08 
     $  - 1.979867116023822d-12*tv09 
     $  - 1.027163441987459d-13*tv10  
     $  + 1.836436639605094d-15*tv11 
     $  + 1.633620651699784d-16*tv12 
     $  + 6.931471805599453d-01*HY3(0,0,1) 
      HY4(0,-1,0,1) = 
     $  - 5.610575179941452d-02 
     $  + 4.649892609082033d-01*tu01 
     $  - 5.631239161843284d-02*tu02 
     $  + 4.220972769653239d-03*tu03 
     $  - 2.066940413626322d-04*tu04 
     $  + 6.100628682175971d-06*tu05 
     $  - 5.412969106099992d-08*tu06 
     $  - 3.230915912784154d-09*tu07 
     $  + 9.249866333323043d-11*tu08 
     $  + 2.685990764581699d-12*tu09 
     $  - 1.543312114608473d-13*tu10 
     $  - 2.036971731594398d-15*tu11 
     $  + 2.517450307574790d-16*tu12 
     $  - 8.224670334241132d-01*HY2(0,-1) 
      HY4(0,-1,-1,1) = 
     $  - 4.031271939759038d-02 
     $  + 3.295217254379970d-01*tu01 
     $  - 4.047097737450547d-02*tu02 
     $  + 3.104955391145708d-03*tu03 
     $  - 1.583251510732719d-04*tu04 
     $  + 5.083334568184305d-06*tu05 
     $  - 6.708598619683341d-08*tu06 
     $  - 1.944278941559733d-09*tu07 
     $  + 8.804863765356287d-11*tu08 
     $  + 9.341312729419985d-13*tu09 
     $  - 1.231746977889946d-13*tu10 
     $  + 3.370647349658755d-16*tu11 
     $  + 1.718647072955689d-16*tu12 
     $  - 5.822405264650125d-01*HY2(0,-1) 
     $  - 6.931471805599453d-01*HY3(0,-1,-1) 
      HY4(0,-1,1,-1) = 
     $  - 4.495764739674318d-02 
     $  - 2.758514579198452d-01*tv01 
     $  - 4.515130668959398d-02*tv02 
     $  - 3.875995092451054d-03*tv03 
     $  - 1.936768370518385d-04*tv04 
     $  - 5.133195476137788d-06*tv05 
     $  - 1.752786900562004d-08*tv06 
     $  + 2.715518363893619d-09*tv07 
     $  + 1.631155670579918d-11*tv08 
     $  - 2.940721244025822d-12*tv09 
     $  - 2.045219059123054d-14*tv10 
     $  + 3.895696592051861d-15*tv11 
     $  + 4.804530139182014d-01*HY2(0,-1) 
     $  + 6.931471805599453d-01*HY3(0,-1,1) 
      HY4(0,1,-1,-1) = 
     $  - 2.782664607935622d-02 
     $  - 1.410831481728889d-01*tv01 
     $  - 2.801876266982354d-02*tv02 
     $  - 2.997894208020603d-03*tv03 
     $  - 1.921960113936824d-04*tv04 
     $  - 7.016503666427137d-06*tv05 
     $  - 7.928257765061337d-08*tv06 
     $  + 4.388745575295455d-09*tv07 
     $  + 1.381107719492586d-10*tv08 
     $  - 4.341921500497716d-12*tv09 
     $  - 2.375364913875066d-13*tv10 
     $  + 4.522044546598701d-15*tv11 
     $  + 4.033357472727688d-16*tv12 
     $  + 2.402265069591007d-01*HY2(0,1) 
      HY4(0,-1,1,1) = 
     $  + 2.782664607935622d-02 
     $  - 1.410831481728889d-01*tu01 
     $  + 2.801876266982354d-02*tu02 
     $  - 2.997894208020603d-03*tu03 
     $  + 1.921960113936824d-04*tu04 
     $  - 7.016503666427137d-06*tu05 
     $  + 7.928257765061337d-08*tu06 
     $  + 4.388745575295455d-09*tu07 
     $  - 1.381107719492586d-10*tu08 
     $  - 4.341921500497716d-12*tu09 
     $  + 2.375364913875066d-13*tu10 
     $  + 4.522044546598701d-15*tu11 
     $  - 4.033357472727688d-16*tu12 
     $  + 2.402265069591007d-01*HY2(0,-1) 
      HY4(0,1,-1,1) = 
     $  + 4.495764739674318d-02 
     $  - 2.758514579198452d-01*tu01 
     $  + 4.515130668959398d-02*tu02 
     $  - 3.875995092451054d-03*tu03 
     $  + 1.936768370518385d-04*tu04 
     $  - 5.133195476137788d-06*tu05 
     $  + 1.752786900562004d-08*tu06 
     $  + 2.715518363893619d-09*tu07 
     $  - 1.631155670579918d-11*tu08 
     $  - 2.940721244025822d-12*tu09 
     $  + 2.045219059123054d-14*tu10 
     $  + 3.895696592051861d-15*tu11 
     $  + 4.804530139182014d-01*HY2(0,1) 
     $  - 6.931471805599453d-01*HY3(0,1,-1) 
      HY4(0,1,1,-1) = 
     $  + 4.031271939759038d-02 
     $  + 3.295217254379970d-01*tv01 
     $  + 4.047097737450547d-02*tv02 
     $  + 3.104955391145708d-03*tv03 
     $  + 1.583251510732719d-04*tv04 
     $  + 5.083334568184305d-06*tv05 
     $  + 6.708598619683341d-08*tv06 
     $  - 1.944278941559733d-09*tv07 
     $  - 8.804863765356287d-11*tv08 
     $  + 9.341312729419985d-13*tv09 
     $  + 1.231746977889946d-13*tv10 
     $  + 3.370647349658755d-16*tv11 
     $  - 1.718647072955689d-16*tv12 
     $  - 5.822405264650125d-01*HY2(0,1) 
     $  + 6.931471805599453d-01*HY3(0,1,1) 
      HY4(-1,-1,-1,1) = 
     $  - 3.768651335815766d-02 
     $  + 3.043162147119780d-01*tu01 
     $  - 3.784162844891144d-02*tu02 
     $  + 2.958351024362477d-03*tu03 
     $  - 1.551924666783514d-04*tu04 
     $  + 5.216293832777793d-06*tu05 
     $  - 7.726843592398867d-08*tu06 
     $  - 1.910379383726989d-09*tu07 
     $  + 1.073377838077624d-10*tu08 
     $  + 4.147979000313175d-13*tu09 
     $  - 1.506593045440627d-13*tu10 
     $  + 1.921276747438603d-15*tu11 
     $  + 1.977332880766160d-16*tu12 
     $  - 5.372131936080402d-01*HY1(-1) 
     $  - 2.911202632325062d-01*HY1(-1)*HY1(-1) 
     $  - 1.155245300933242d-01*HY1(-1)*HY1(-1)*HY1(-1) 
      HY4(-1,-1,1,1) = 
     $  + 2.908893189635991d-02 
     $  - 1.784837106345115d-01*tu01 
     $  + 2.927117884632272d-02*tu02 
     $  - 2.888221776586007d-03*tu03 
     $  + 1.823501630828519d-04*tu04 
     $  - 6.976883920991888d-06*tu05 
     $  + 1.030302948541690d-07*tu06 
     $  + 3.794029548474434d-09*tu07 
     $  - 1.825184393299693d-10*tu08 
     $  - 2.300206200729610d-12*tu09 
     $  + 3.062629564489397d-13*tu10 
     $  - 7.629393984387632d-16*tu11 
     $  - 4.860728618463296d-16*tu12 
     $  + 3.088253750968339d-01*HY1(-1) 
     $  + 1.201132534795503d-01*HY1(-1)*HY1(-1) 
      HY4(-1,1,1,1) = 
     $  - 9.029205146496301d-03 
     $  + 3.753824045412342d-02*tu01 
     $  - 9.240717745810759d-03*tu02 
     $  + 2.351153976182453d-03*tu03 
     $  - 2.115782190216214d-04*tu04 
     $  + 8.486524807740892d-06*tu05 
     $  - 6.547885807612483d-08*tu06 
     $  - 6.934422754020238d-09*tu07 
     $  + 1.405695202725693d-10*tu08 
     $  + 8.329441237576153d-12*tu09 
     $  - 2.790404594803712d-13*tu10 
     $  - 1.024489568815216d-14*tu11 
     $  + 5.256388245544115d-16*tu12 
     $  - 5.550410866482157d-02*HY1(-1) 
      endif 
** nw > 3 endif 
      endif 
** (n1,n2) = (-1,1) -- completion endif 
      return 
      end 
************************************************************************ 
      subroutine fillirr1dhplat1(r,nw,HR1,HR2,HR3,HR4, 
     $                                HY1,HY2,HY3,HY4, 
     $                                Hi1,Hi2,Hi3,Hi4,n1,n2) 
** evaluates the HPL for r2m1 < y < r2p1
** fillirr1dhplat1 is called by eval1dhplat1 after calling 
** fillirr1dhplat0 with argument r=(1-y)/(1+y) 
** it is guaranteed that nw is in the range 2:4, and that (n1,n2) 
** take one of the pairs of values (0,1), (-1,0) or (-1,1) 
      implicit double precision (a-h,o-z) 
      dimension HR1(-1:1),HR2(-1:1,-1:1),HR3(-1:1,-1:1,-1:1), 
     $          HR4(-1:1,-1:1,-1:1,-1:1) 
      dimension HY1(n1:n2),HY2(n1:n2,n1:n2),HY3(n1:n2,n1:n2,n1:n2), 
     $          HY4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2) 
** (n1,n2) = (0,1) or (-1,1) 
      if (    ( (n1.eq.0).and.(n2.eq.1) ) 
     $    .or.( (n1.eq.-1).and.(n2.eq.1) ) ) then 
      HY2(0,1) = 
     $  + 1.6449340668482264d+00
     $  + 6.9314718055994530d-01*HR1(-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)
     $  + HR1( -1)*HR1(0)
     $  - HR1( -1)*HR1(1)
     $  + HR1(0) *HR1(1)
     $  + 6.9314718055994530d-01*HR1(1)
     $  + HR2( -1,1)
     $  - HR2(0, -1)
     $  - HR2(0,1) 
      if (r.lt.0d0) then 
      Hi2(0,1) = 
     $  - HR1( -1) 
     $  - HR1(1) 
      endif 
      if ( nw.gt.2 ) then 
      HY3(0,0,1) = 
     $  + 1.2020569031595942d+00
     $  - 1.6449340668482264d+00*HR1(-1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - HR1( -1)*HR1(0)*HR1(1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(1)*HR1(1)
     $  + HR1( -1)*HR2(0,-1)
     $  + HR1( -1)*HR2(0,1)
     $  - 5.0000000000000000d-01*HR1(0)*HR1(1)*HR1(1)
     $  - 1.6449340668482264d+00*HR1(1)
     $  - 3.4657359027997265d-01*HR1(1)*HR1(1)
     $  - HR1(1) *HR2(-1,1)
     $  + HR1(1) *HR2(0,-1)
     $  + HR1(1) *HR2(0,1)
     $  - HR3( -1,-1,1)
     $  + HR3( -1,1,1)
     $  - HR3(0, -1,-1)
     $  - HR3(0, -1,1)
     $  - HR3(0,1, -1)
     $  - HR3(0,1,1) 
      HY3(0,1,1) = 
     $  + 1.2020569031595942d+00
     $  - 2.4022650695910071d-01*HR1(-1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(0)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(0)*HR1(0)
     $  + HR1( -1)*HR1(0)*HR1(1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(1)
     $  + HR1( -1)*HR2(-1,1)
     $  - 5.0000000000000000d-01*HR1(0)*HR1(0)*HR1(1)
     $  - 6.9314718055994530d-01*HR1(0)*HR1(1)
     $  - HR1(0) *HR2(-1,1)
     $  + HR1(0) *HR2(0,-1)
     $  + HR1(0) *HR2(0,1)
     $  - 2.4022650695910071d-01*HR1(1)
     $  - 6.9314718055994530d-01*HR2(-1,1)
     $  + 6.9314718055994530d-01*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR2(0,1)
     $  - HR3( -1,-1,1)
     $  - HR3(0, -1,-1)
     $  - HR3(0,0, -1)
     $  - HR3(0,0,1) 
     $  - HR3(0,1, -1)
      if (r.lt.0d0) then 
      HY3(0,1,1) = HY3(0,1,1) 
     $  + 4.9348022005446793d+00*HR1(-1)
     $  + 4.9348022005446793d+00*HR1(1)
      Hi3(0,0,1) = 
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1) 
     $  + HR1( -1)*HR1(1) 
     $  + 5.0000000000000000d-01*HR1(1)*HR1(1) 
      Hi3(0,1,1) = 
     $  + 6.9314718055994530d-01*HR1(-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)
     $  + HR1( -1)*HR1(0)
     $  - HR1( -1)*HR1(1)
     $  + HR1(0) *HR1(1)
     $  + 6.9314718055994530d-01*HR1(1)
     $  + HR2( -1,1)
     $  - HR2(0, -1)
     $  - HR2(0,1) 
      endif
      endif
      if ( nw.gt.3 ) then 
      HY4(0,0,0,1) = 
     $  + 1.0823232337111381d+00
     $  - 1.2020569031595942d+00*HR1(-1)
     $  + 8.2246703342411321d-01*HR1(-1)*HR1(-1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - 2.5000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(0,1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(0)*HR1(1)*HR1(1)
     $  + 1.6449340668482264d+00*HR1(-1)*HR1(1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(1)*HR1(1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(1)*HR1(1)*HR1(1)
     $  - HR1( -1)*HR1(1)*HR2(0,-1)
     $  - HR1( -1)*HR1(1)*HR2(0,1)
     $  + HR1( -1)*HR3(0,-1,-1)
     $  + HR1( -1)*HR3(0,-1,1)
     $  + HR1( -1)*HR3(0,1,-1)
     $  + HR1( -1)*HR3(0,1,1)
     $  + 1.6666666666666666d-01*HR1(0)*HR1(1)*HR1(1)*HR1(1)
     $  - 1.2020569031595942d+00*HR1(1)
     $  + 8.2246703342411321d-01*HR1(1)*HR1(1)
     $  + 1.1552453009332421d-01*HR1(1)*HR1(1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(1)*HR1(1)*HR2(-1,1)
     $  - 5.0000000000000000d-01*HR1(1)*HR1(1)*HR2(0,-1)
     $  - 5.0000000000000000d-01*HR1(1)*HR1(1)*HR2(0,1)
     $  + HR1(1) *HR3(-1,-1,1)
     $  - HR1(1) *HR3(-1,1,1)
     $  + HR1(1) *HR3(0,-1,-1)
     $  + HR1(1) *HR3(0,-1,1)
     $  + HR1(1) *HR3(0,1,-1)
     $  + HR1(1) *HR3(0,1,1)
     $  + HR4( -1,-1,-1,1)
     $  - HR4( -1,-1,1,1)
     $  + HR4( -1,1,1,1)
     $  - HR4(0, -1,-1,-1)
     $  - HR4(0, -1,-1,1)
     $  - HR4(0, -1,1,-1)
     $  - HR4(0, -1,1,1)
     $  - HR4(0,1, -1,-1)
     $  - HR4(0,1, -1,1)
     $  - HR4(0,1,1, -1)
     $  - HR4(0,1,1,1) 
      HY4(0,0,1,1) = 
     $  + 2.7058080842778454d-01
     $  - 1.2020569031595942d+00*HR1(-1)
     $  + 1.2011325347955035d-01*HR1(-1)*HR1(-1)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  + 2.5000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(0)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  + 2.5000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(0)*HR1(0)*HR1(1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(0)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(0)*HR1(1)*HR1(1)
     $  - HR1( -1)*HR1(0)*HR2(0,-1)
     $  - HR1( -1)*HR1(0)*HR2(0,1)
     $  + 2.4022650695910071d-01*HR1(-1)*HR1(1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(1)*HR1(1)
     $  - HR1( -1)*HR1(1)*HR2(-1,1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR2(0,1)
     $  - HR1( -1)*HR3(-1,-1,1)
     $  + HR1( -1)*HR3(-1,1,1)
     $  + HR1( -1)*HR3(0,-1,-1)
     $  + HR1( -1)*HR3(0,0,-1)
     $  + HR1( -1)*HR3(0,0,1)
     $  + HR1( -1)*HR3(0,1,-1)
     $  + 2.5000000000000000d-01*HR1(0)*HR1(0)*HR1(1)*HR1(1)
     $  + 3.4657359027997265d-01*HR1(0)*HR1(1)*HR1(1)
     $  + HR1(0) *HR1(1)*HR2(-1,1)
     $  - HR1(0) *HR1(1)*HR2(0,-1)
     $  - HR1(0) *HR1(1)*HR2(0,1)
     $  + HR1(0) *HR3(-1,-1,1)
     $  - HR1(0) *HR3(-1,1,1)
     $  + HR1(0) *HR3(0,-1,-1)
     $  + HR1(0) *HR3(0,-1,1)
     $  + HR1(0) *HR3(0,1,-1)
     $  + HR1(0) *HR3(0,1,1)
     $  - 1.2020569031595942d+00*HR1(1)
     $  + 1.2011325347955035d-01*HR1(1)*HR1(1)
     $  + 6.9314718055994530d-01*HR1(1)*HR2(-1,1)
     $  - 6.9314718055994530d-01*HR1(1)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(1)*HR2(0,1)
     $  + HR1(1) *HR3(-1,-1,1)
     $  + HR1(1) *HR3(0,-1,-1)
     $  + HR1(1) *HR3(0,0,-1)
     $  + HR1(1) *HR3(0,0,1)
     $  + HR1(1) *HR3(0,1,-1)
     $  + 6.9314718055994530d-01*HR3(-1,-1,1)
     $  - 6.9314718055994530d-01*HR3(-1,1,1)
     $  + 6.9314718055994530d-01*HR3(0,-1,-1)
     $  + 6.9314718055994530d-01*HR3(0,-1,1)
     $  + 6.9314718055994530d-01*HR3(0,1,-1)
     $  + 6.9314718055994530d-01*HR3(0,1,1)
     $  + 2.0000000000000000d+00*HR4(-1,-1,-1,1)
     $  - HR4( -1,-1,1,1)
     $  - 2.0000000000000000d+00*HR4(0,-1,-1,-1)
     $  - HR4(0, -1,-1,1)
     $  - HR4(0, -1,1,-1)
     $  - HR4(0,0, -1,-1)
     $  - HR4(0,0, -1,1)
     $  - HR4(0,0,1, -1)
     $  - HR4(0,0,1,1) 
     $  - 2.0000000000000000d+00*HR4(0,1,-1,-1)
     $  - HR4(0,1, -1,1)
     $  - HR4(0,1,1, -1)
      HY4(0,1,1,1) = 
     $  + 1.0823232337111381d+00
     $  + 5.5504108664821579d-02*HR1(-1)
     $  - 1.2011325347955035d-01*HR1(-1)*HR1(-1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  - 2.5000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(0)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  + 2.4022650695910071d-01*HR1(-1)*HR1(0)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(0)*HR1(0)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(0)*HR1(0)*HR1(0)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(0)*HR1(0)*HR1(1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(0)*HR1(1)
     $  - HR1( -1)*HR1(0)*HR2(-1,1)
     $  - 2.4022650695910071d-01*HR1(-1)*HR1(1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR2(-1,1)
     $  - HR1( -1)*HR3(-1,-1,1)
     $  + 1.6666666666666666d-01*HR1(0)*HR1(0)*HR1(0)*HR1(1)
     $  + 3.4657359027997265d-01*HR1(0)*HR1(0)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(0)*HR1(0)*HR2(-1,1)
     $  - 5.0000000000000000d-01*HR1(0)*HR1(0)*HR2(0,-1)
     $  - 5.0000000000000000d-01*HR1(0)*HR1(0)*HR2(0,1)
     $  + 2.4022650695910071d-01*HR1(0)*HR1(1)
     $  + 6.9314718055994530d-01*HR1(0)*HR2(-1,1)
     $  - 6.9314718055994530d-01*HR1(0)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(0)*HR2(0,1)
     $  + HR1(0) *HR3(-1,-1,1)
     $  + HR1(0) *HR3(0,-1,-1)
     $  + HR1(0) *HR3(0,0,-1)
     $  + HR1(0) *HR3(0,0,1)
     $  + HR1(0) *HR3(0,1,-1)
     $  + 5.5504108664821579d-02*HR1(1)
     $  + 2.4022650695910071d-01*HR2(-1,1)
     $  - 2.4022650695910071d-01*HR2(0,-1)
     $  - 2.4022650695910071d-01*HR2(0,1)
     $  + 6.9314718055994530d-01*HR3(-1,-1,1)
     $  + 6.9314718055994530d-01*HR3(0,-1,-1)
     $  + 6.9314718055994530d-01*HR3(0,0,-1)
     $  + 6.9314718055994530d-01*HR3(0,0,1)
     $  + 6.9314718055994530d-01*HR3(0,1,-1)
     $  + HR4( -1,-1,-1,1)
     $  - HR4(0, -1,-1,-1)
     $  - HR4(0,0, -1,-1)
     $  - HR4(0,0,0, -1)
     $  - HR4(0,0,0,1) 
     $  - HR4(0,0,1, -1)
     $  - HR4(0,1, -1,-1)
      if (r.lt.0d0) then 
      HY4(0,0,1,1) = HY4(0,0,1,1) 
     $  - 2.4674011002723396d+00*HR1(-1)*HR1(-1)
     $  - 4.9348022005446793d+00*HR1(-1)*HR1(1)
     $  - 2.4674011002723396d+00*HR1(1)*HR1(1)
      HY4(0,1,1,1) = HY4(0,1,1,1) 
     $  - 3.4205442319285582d+00*HR1(-1)
     $  + 2.4674011002723396d+00*HR1(-1)*HR1(-1)
     $  - 4.9348022005446793d+00*HR1(-1)*HR1(0)
     $  + 4.9348022005446793d+00*HR1(-1)*HR1(1)
     $  - 4.9348022005446793d+00*HR1(0)*HR1(1)
     $  - 3.4205442319285582d+00*HR1(1)
     $  - 4.9348022005446793d+00*HR2(-1,1)
     $  + 4.9348022005446793d+00*HR2(0,-1)
     $  + 4.9348022005446793d+00*HR2(0,1)
      Hi4(0,0,0,1) = 
     $  - 1.666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1) 
     $  - 5.000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1) 
     $  - 5.000000000000000d-01*HR1(-1)*HR1(1)*HR1(1) 
     $  - 1.666666666666666d-01*HR1(1)*HR1(1)*HR1(1) 
      Hi4(0,0,1,1) = 
     $  - 3.465735902799726d-01*HR1(-1)*HR1(-1) 
     $  + 1.666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1) 
     $  - 5.000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0) 
     $  + 5.000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1) 
     $  - HR1( -1)*HR1(0)*HR1(1) 
     $  - 6.931471805599453d-01*HR1(-1)*HR1(1) 
     $  + 5.000000000000000d-01*HR1(-1)*HR1(1)*HR1(1) 
     $  + HR1( -1)*HR2(0,-1) 
     $  + HR1( -1)*HR2(0,1) 
     $  - 5.000000000000000d-01*HR1(0)*HR1(1)*HR1(1) 
     $  - 3.465735902799726d-01*HR1(1)*HR1(1) 
     $  - HR1(1) *HR2(-1,1) 
     $  + HR1(1) *HR2(0,-1) 
     $  + HR1(1) *HR2(0,1) 
     $  - HR3( -1,-1,1) 
     $  + HR3( -1,1,1) 
     $  - HR3(0, -1,-1) 
     $  - HR3(0, -1,1) 
     $  - HR3(0,1, -1) 
     $  - HR3(0,1,1) 
      Hi4(0,1,1,1) = 
     $  + 1.404707559889125d+00*HR1(-1) 
     $  + 3.465735902799726d-01*HR1(-1)*HR1(-1) 
     $  - 1.666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1) 
     $  + 5.000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0) 
     $  - 5.000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1) 
     $  - 6.931471805599453d-01*HR1(-1)*HR1(0) 
     $  - 5.000000000000000d-01*HR1(-1)*HR1(0)*HR1(0) 
     $  + HR1( -1)*HR1(0)*HR1(1) 
     $  + 6.931471805599453d-01*HR1(-1)*HR1(1) 
     $  + HR1( -1)*HR2(-1,1) 
     $  - 5.000000000000000d-01*HR1(0)*HR1(0)*HR1(1) 
     $  - 6.931471805599453d-01*HR1(0)*HR1(1) 
     $  - HR1(0) *HR2(-1,1) 
     $  + HR1(0) *HR2(0,-1) 
     $  + HR1(0) *HR2(0,1) 
     $  + 1.404707559889125d+00*HR1(1) 
     $  - 6.931471805599453d-01*HR2(-1,1) 
     $  + 6.931471805599453d-01*HR2(0,-1) 
     $  + 6.931471805599453d-01*HR2(0,1) 
     $  - HR3( -1,-1,1) 
     $  - HR3(0, -1,-1) 
     $  - HR3(0,0, -1) 
     $  - HR3(0,0,1)  
     $  - HR3(0,1, -1) 
      endif 
      endif 
** nw > 3 endif 
      endif 
** (n1,n2) = (0,1) or (-1,1) endif 
************ 
** (n1,n2) = (-1,0) or (-1,1) 
      if (    ( (n1.eq.-1).and.(n2.eq.0) ) 
     $    .or.( (n1.eq.-1).and.(n2.eq.1) ) ) then 
       HY2(0,-1) = 
     $  + 8.2246703342411321d-01
     $  - 6.9314718055994530d-01*HR1(-1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)
     $  + HR1( -1)*HR1(1)
     $  - 6.9314718055994530d-01*HR1(1)
     $  - HR2( -1,1)
      if ( nw.gt.2 ) then 
      HY3(0,0,-1) = 
     $  + 9.0154267736969571d-01
     $  - 8.2246703342411321d-01*HR1(-1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(1)*HR1(1)
     $  - 8.2246703342411321d-01*HR1(1)
     $  + 3.4657359027997265d-01*HR1(1)*HR1(1)
     $  + HR1(1) *HR2(-1,1)
     $  + HR3( -1,-1,1)
     $  - HR3( -1,1,1)
      HY3(0,-1,-1) = 
     $  + 1.5025711289494928d-01
     $  - 2.4022650695910071d-01*HR1(-1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(1)
     $  + HR1( -1)*HR2(-1,1)
     $  - 2.4022650695910071d-01*HR1(1)
     $  - 6.9314718055994530d-01*HR2(-1,1)
     $  - HR3( -1,-1,1)
      endif 
      if ( nw.gt.3 ) then 
      HY4(0,0,0,-1) = 
     $  + 9.4703282949724591d-01
     $  - 9.0154267736969571d-01*HR1(-1)
     $  + 4.1123351671205660d-01*HR1(-1)*HR1(-1)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  + 2.5000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  + 8.2246703342411321d-01*HR1(-1)*HR1(1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(1)*HR1(1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(1)*HR1(1)*HR1(1)
     $  - 9.0154267736969571d-01*HR1(1)
     $  + 4.1123351671205660d-01*HR1(1)*HR1(1)
     $  - 1.1552453009332421d-01*HR1(1)*HR1(1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(1)*HR1(1)*HR2(-1,1)
     $  - HR1(1) *HR3(-1,-1,1)
     $  + HR1(1) *HR3(-1,1,1)
     $  - HR4( -1,-1,-1,1)
     $  + HR4( -1,-1,1,1)
     $  - HR4( -1,1,1,1)
      HY4(0,0,-1,-1) = 
     $  + 8.7785671568655302d-02
     $  - 1.5025711289494928d-01*HR1(-1)
     $  + 1.2011325347955035d-01*HR1(-1)*HR1(-1)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  + 2.5000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  + 2.4022650695910071d-01*HR1(-1)*HR1(1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(1)*HR1(1)
     $  - HR1( -1)*HR1(1)*HR2(-1,1)
     $  - HR1( -1)*HR3(-1,-1,1)
     $  + HR1( -1)*HR3(-1,1,1)
     $  - 1.5025711289494928d-01*HR1(1)
     $  + 1.2011325347955035d-01*HR1(1)*HR1(1)
     $  + 6.9314718055994530d-01*HR1(1)*HR2(-1,1)
     $  + HR1(1) *HR3(-1,-1,1)
     $  + 6.9314718055994530d-01*HR3(-1,-1,1)
     $  - 6.9314718055994530d-01*HR3(-1,1,1)
     $  + 2.0000000000000000d+00*HR4(-1,-1,-1,1)
     $  - HR4( -1,-1,1,1)
      HY4(0,-1,-1,-1) = 
     $  + 2.3752366322618485d-02
     $  - 5.5504108664821579d-02*HR1(-1)
     $  + 1.2011325347955035d-01*HR1(-1)*HR1(-1)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  + 2.4022650695910071d-01*HR1(-1)*HR1(1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR2(-1,1)
     $  + HR1( -1)*HR3(-1,-1,1)
     $  - 5.5504108664821579d-02*HR1(1)
     $  - 2.4022650695910071d-01*HR2(-1,1)
     $  - 6.9314718055994530d-01*HR3(-1,-1,1)
     $  - HR4( -1,-1,-1,1)
      endif 
** nw > 3 endif 
      endif 
** (n1,n2) = (-1,0) or (-1,1) endif 
** (n1,n2) = (-1,1) -- completion 
      if ( (n1.eq.-1).and.(n2.eq.1) ) then 
      HY2(-1,1) = 
     $  + 5.8224052646501250d-01
     $  + 6.9314718055994530d-01*HR1(-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)
     $  + HR1( -1)*HR1(0)
     $  - HR2(0, -1)
      if (r.lt.0d0) then 
      Hi2(-1,1) = 
     $  - HR1( -1) 
      endif 
      if ( nw.gt.2 ) then 
      HY3(0,-1,1) = 
     $  + 2.4307035167006157d-01
     $  - 5.8224052646501250d-01*HR1(-1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - HR1( -1)*HR1(0)*HR1(1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(1)
     $  - HR1( -1)*HR2(-1,1)
     $  + HR1( -1)*HR2(0,-1)
     $  + HR1(0) *HR2(-1,1)
     $  - 5.8224052646501250d-01*HR1(1)
     $  + HR1(1) *HR2(0,-1)
     $  + 6.9314718055994530d-01*HR2(-1,1)
     $  + HR3( -1,-1,1)
     $  - HR3(0, -1,-1)
     $  - HR3(0, -1,1)
      HY3(0,1,-1) = 
     $  + 5.0821521280468485d-01
     $  + 1.0626935403832139d+00*HR1(-1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(0)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(1)
     $  - HR1( -1)*HR2(-1,1)
     $  - HR1( -1)*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR1(0)*HR1(1)
     $  + 1.0626935403832139d+00*HR1(1)
     $  - HR1(1) *HR2(0,-1)
     $  + 6.9314718055994530d-01*HR2(-1,1)
     $  - 6.9314718055994530d-01*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR2(0,1)
     $  + HR3( -1,-1,1)
     $  + 2.0000000000000000d+00*HR3(0,-1,-1)
     $  + HR3(0, -1,1)
     $  + HR3(0,1, -1)
      HY3(-1,-1,1) = 
     $  + 9.4753004230127705d-02
     $  - 5.8224052646501250d-01*HR1(-1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  + HR1( -1)*HR2(0,-1)
     $  - HR3(0, -1,-1)
      HY3(-1,1,1) = 
     $  + 5.3721319360804020d-01
     $  - 2.4022650695910071d-01*HR1(-1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(0)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(0)*HR1(0)
     $  + HR1(0) *HR2(0,-1)
     $  + 6.9314718055994530d-01*HR2(0,-1)
     $  - HR3(0, -1,-1)
     $  - HR3(0,0, -1)
      if (r.lt.0d0) then 
      HY3(-1,1,1) = HY3(-1,1,1) 
     $  + 4.9348022005446793d+00*HR1(-1)
      Hi3(0,-1,1) = 
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1) 
     $  + HR1( -1)*HR1(1) 
     $  - HR2( -1,1) 
      Hi3(0,1,-1) = 
     $  - 6.9314718055994530d-01*HR1(-1) 
     $  - 6.9314718055994530d-01*HR1(1) 
      Hi3(-1,-1,1) = 
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1) 
      Hi3(-1,1,1) = 
     $  + 6.9314718055994530d-01*HR1(-1) 
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1) 
     $  + HR1( -1)*HR1(0) 
     $  - HR2(0, -1) 
      endif 
      endif 
      if ( nw.gt.3 ) then 
      HY4(0,0,-1,1) = 
     $  + 1.1787599965050932d-01
     $  - 2.4307035167006157d-01*HR1(-1)
     $  + 2.9112026323250625d-01*HR1(-1)*HR1(-1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - 2.5000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(0)*HR1(1)*HR1(1)
     $  + 5.8224052646501250d-01*HR1(-1)*HR1(1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(1)*HR1(1)
     $  + HR1( -1)*HR1(1)*HR2(-1,1)
     $  - HR1( -1)*HR1(1)*HR2(0,-1)
     $  + HR1( -1)*HR3(-1,-1,1)
     $  - HR1( -1)*HR3(-1,1,1)
     $  + HR1( -1)*HR3(0,-1,-1)
     $  + HR1( -1)*HR3(0,-1,1)
     $  - HR1(0) *HR1(1)*HR2(-1,1)
     $  - HR1(0) *HR3(-1,-1,1)
     $  + HR1(0) *HR3(-1,1,1)
     $  - 2.4307035167006157d-01*HR1(1)
     $  + 2.9112026323250625d-01*HR1(1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(1)*HR1(1)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(1)*HR2(-1,1)
     $  - HR1(1) *HR3(-1,-1,1)
     $  + HR1(1) *HR3(0,-1,-1)
     $  + HR1(1) *HR3(0,-1,1)
     $  - 6.9314718055994530d-01*HR3(-1,-1,1)
     $  + 6.9314718055994530d-01*HR3(-1,1,1)
     $  - 2.0000000000000000d+00*HR4(-1,-1,-1,1)
     $  + HR4( -1,-1,1,1)
     $  - HR4(0, -1,-1,-1)
     $  - HR4(0, -1,-1,1)
     $  - HR4(0, -1,1,-1)
     $  - HR4(0, -1,1,1)
      HY4(0,0,1,-1) = 
     $  + 1.7284527823898438d-01
     $  - 5.0821521280468485d-01*HR1(-1)
     $  - 5.3134677019160696d-01*HR1(-1)*HR1(-1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - 2.5000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(0)*HR1(1)
     $  - 1.0626935403832139d+00*HR1(-1)*HR1(1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(1)*HR1(1)
     $  + HR1( -1)*HR1(1)*HR2(-1,1)
     $  + HR1( -1)*HR1(1)*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR2(0,1)
     $  + HR1( -1)*HR3(-1,-1,1)
     $  - HR1( -1)*HR3(-1,1,1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR3(0,-1,-1)
     $  - HR1( -1)*HR3(0,-1,1)
     $  - HR1( -1)*HR3(0,1,-1)
     $  - 3.4657359027997265d-01*HR1(0)*HR1(1)*HR1(1)
     $  - 5.0821521280468485d-01*HR1(1)
     $  - 5.3134677019160696d-01*HR1(1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(1)*HR1(1)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(1)*HR2(-1,1)
     $  + 6.9314718055994530d-01*HR1(1)*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR1(1)*HR2(0,1)
     $  - HR1(1) *HR3(-1,-1,1)
     $  - 2.0000000000000000d+00*HR1(1)*HR3(0,-1,-1)
     $  - HR1(1) *HR3(0,-1,1)
     $  - HR1(1) *HR3(0,1,-1)
     $  - 6.9314718055994530d-01*HR3(-1,-1,1)
     $  + 6.9314718055994530d-01*HR3(-1,1,1)
     $  - 6.9314718055994530d-01*HR3(0,-1,-1)
     $  - 6.9314718055994530d-01*HR3(0,-1,1)
     $  - 6.9314718055994530d-01*HR3(0,1,-1)
     $  - 6.9314718055994530d-01*HR3(0,1,1)
     $  - 2.0000000000000000d+00*HR4(-1,-1,-1,1)
     $  + HR4( -1,-1,1,1)
     $  + 3.0000000000000000d+00*HR4(0,-1,-1,-1)
     $  + 2.0000000000000000d+00*HR4(0,-1,-1,1)
     $  + 2.0000000000000000d+00*HR4(0,-1,1,-1)
     $  + HR4(0, -1,1,1)
     $  + 2.0000000000000000d+00*HR4(0,1,-1,-1)
     $  + HR4(0,1, -1,1)
     $  + HR4(0,1,1, -1)
      HY4(0,-1,0,1) = 
     $  + 2.0293560632083841d-01
     $  - 3.8889584616810632d-01*HR1(-1)
     $  + 8.2246703342411321d-01*HR1(-1)*HR1(-1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(0,1)
     $  - HR1( -1)*HR1(0)*HR2(-1,1)
     $  + 1.6449340668482264d+00*HR1(-1)*HR1(1)
     $  - HR1( -1)*HR1(1)*HR2(-1,1)
     $  - HR1( -1)*HR1(1)*HR2(0,-1)
     $  - HR1( -1)*HR1(1)*HR2(0,1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR2(-1,1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR3(-1,-1,1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR3(-1,1,1)
     $  + HR1( -1)*HR3(0,-1,-1)
     $  + HR1( -1)*HR3(0,1,-1)
     $  + HR1(0) *HR1(1)*HR2(-1,1)
     $  + 2.0000000000000000d+00*HR1(0)*HR3(-1,-1,1)
     $  - 2.0000000000000000d+00*HR1(0)*HR3(-1,1,1)
     $  - 3.8889584616810632d-01*HR1(1)
     $  + 6.9314718055994530d-01*HR1(1)*HR2(-1,1)
     $  + 2.0000000000000000d+00*HR1(1)*HR3(-1,-1,1)
     $  + HR1(1) *HR3(0,-1,-1)
     $  + HR1(1) *HR3(0,1,-1)
     $  - 1.6449340668482264d+00*HR2(-1,1)
     $  - 5.0000000000000000d-01*HR2(-1,1)*HR2(-1,1)
     $  + HR2( -1,1)*HR2(0,-1)
     $  + HR2( -1,1)*HR2(0,1)
     $  + 1.3862943611198906d+00*HR3(-1,-1,1)
     $  - 1.3862943611198906d+00*HR3(-1,1,1)
     $  + 4.0000000000000000d+00*HR4(-1,-1,-1,1)
     $  - 2.0000000000000000d+00*HR4(-1,-1,1,1)
     $  - HR4(0, -1,-1,-1)
     $  - HR4(0, -1,-1,1)
     $  - HR4(0,1, -1,-1)
     $  - HR4(0,1, -1,1)
      HY4(0,-1,-1,1) = 
     $  + 3.4159126166513913d-02
     $  - 9.4753004230127705d-02*HR1(-1)
     $  + 2.9112026323250625d-01*HR1(-1)*HR1(-1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  - HR1( -1)*HR1(0)*HR2(-1,1)
     $  + 5.8224052646501250d-01*HR1(-1)*HR1(1)
     $  - HR1( -1)*HR1(1)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR2(-1,1)
     $  - HR1( -1)*HR3(-1,-1,1)
     $  + HR1( -1)*HR3(0,-1,-1)
     $  + HR1(0) *HR3(-1,-1,1)
     $  - 9.4753004230127705d-02*HR1(1)
     $  + HR1(1) *HR3(0,-1,-1)
     $  - 5.8224052646501250d-01*HR2(-1,1)
     $  + HR2( -1,1)*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR3(-1,-1,1)
     $  + HR4( -1,-1,-1,1)
     $  - HR4(0, -1,-1,-1)
     $  - HR4(0, -1,-1,1)
      HY4(0,-1,1,-1) = 
     $  + 5.4653052738263652d-02
     $  - 2.1407237086670622d-01*HR1(-1)
     $  - 5.3134677019160696d-01*HR1(-1)*HR1(-1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(0)*HR1(1)
     $  - 1.0626935403832139d+00*HR1(-1)*HR1(1)
     $  + HR1( -1)*HR1(1)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR2(-1,1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR2(0,-1)
     $  - HR1( -1)*HR3(-1,-1,1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR3(0,-1,-1)
     $  + 6.9314718055994530d-01*HR1(0)*HR2(-1,1)
     $  - 2.1407237086670622d-01*HR1(1)
     $  + 6.9314718055994530d-01*HR1(1)*HR2(0,-1)
     $  - 2.0000000000000000d+00*HR1(1)*HR3(0,-1,-1)
     $  + 1.0626935403832139d+00*HR2(-1,1)
     $  - HR2( -1,1)*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR3(-1,-1,1)
     $  - 6.9314718055994530d-01*HR3(0,-1,-1)
     $  - 6.9314718055994530d-01*HR3(0,-1,1)
     $  + HR4( -1,-1,-1,1)
     $  + 3.0000000000000000d+00*HR4(0,-1,-1,-1)
     $  + 2.0000000000000000d+00*HR4(0,-1,-1,1)
     $  + HR4(0, -1,1,-1)
      HY4(0,1,-1,-1) = 
     $  + 1.1412342741606084d-01
     $  + 4.7533770109129867d-01*HR1(-1)
     $  - 1.2011325347955035d-01*HR1(-1)*HR1(-1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  + 2.4022650695910071d-01*HR1(-1)*HR1(0)
     $  - 2.4022650695910071d-01*HR1(-1)*HR1(1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR2(-1,1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR2(0,-1)
     $  - HR1( -1)*HR3(-1,-1,1)
     $  + HR1( -1)*HR3(0,-1,-1)
     $  + 2.4022650695910071d-01*HR1(0)*HR1(1)
     $  + 4.7533770109129867d-01*HR1(1)
     $  - 6.9314718055994530d-01*HR1(1)*HR2(0,-1)
     $  + HR1(1) *HR3(0,-1,-1)
     $  + 2.4022650695910071d-01*HR2(-1,1)
     $  - 2.4022650695910071d-01*HR2(0,-1)
     $  - 2.4022650695910071d-01*HR2(0,1)
     $  + 6.9314718055994530d-01*HR3(-1,-1,1)
     $  + 1.3862943611198906d+00*HR3(0,-1,-1)
     $  + 6.9314718055994530d-01*HR3(0,-1,1)
     $  + 6.9314718055994530d-01*HR3(0,1,-1)
     $  + HR4( -1,-1,-1,1)
     $  - 3.0000000000000000d+00*HR4(0,-1,-1,-1)
     $  - HR4(0, -1,-1,1)
     $  - HR4(0, -1,1,-1)
     $  - HR4(0,1, -1,-1)
      HY4(0,-1,1,1) = 
     $  + 9.3097125991768577d-02
     $  - 5.3721319360804020d-01*HR1(-1)
     $  + 1.2011325347955035d-01*HR1(-1)*HR1(-1)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  + 2.5000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(0)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(0)*HR1(0)*HR1(1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(0)*HR1(1)
     $  + HR1( -1)*HR1(0)*HR2(-1,1)
     $  - HR1( -1)*HR1(0)*HR2(0,-1)
     $  + 2.4022650695910071d-01*HR1(-1)*HR1(1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR2(-1,1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR2(0,-1)
     $  + HR1( -1)*HR3(-1,-1,1)
     $  + HR1( -1)*HR3(0,-1,-1)
     $  + HR1( -1)*HR3(0,0,-1)
     $  - 5.0000000000000000d-01*HR1(0)*HR1(0)*HR2(-1,1)
     $  - HR1(0) *HR1(1)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(0)*HR2(-1,1)
     $  - HR1(0) *HR3(-1,-1,1)
     $  + HR1(0) *HR3(0,-1,-1)
     $  + HR1(0) *HR3(0,-1,1)
     $  - 5.3721319360804020d-01*HR1(1)
     $  - 6.9314718055994530d-01*HR1(1)*HR2(0,-1)
     $  + HR1(1) *HR3(0,-1,-1)
     $  + HR1(1) *HR3(0,0,-1)
     $  - 2.4022650695910071d-01*HR2(-1,1)
     $  - 6.9314718055994530d-01*HR3(-1,-1,1)
     $  + 6.9314718055994530d-01*HR3(0,-1,-1)
     $  + 6.9314718055994530d-01*HR3(0,-1,1)
     $  - HR4( -1,-1,-1,1)
     $  - 2.0000000000000000d+00*HR4(0,-1,-1,-1)
     $  - HR4(0, -1,-1,1)
     $  - HR4(0, -1,1,-1)
     $  - HR4(0,0, -1,-1)
     $  - HR4(0,0, -1,1)
      HY4(0,1,-1,1) = 
     $  + 1.9355535381306524d-01
     $  + 1.4780047665430420d+00*HR1(-1)
     $  - 2.9112026323250625d-01*HR1(-1)*HR1(-1)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  + 5.8224052646501250d-01*HR1(-1)*HR1(0)
     $  + HR1( -1)*HR1(0)*HR2(-1,1)
     $  + HR1( -1)*HR1(0)*HR2(0,-1)
     $  - 5.8224052646501250d-01*HR1(-1)*HR1(1)
     $  + HR1( -1)*HR1(1)*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR2(-1,1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR2(0,-1)
     $  + HR1( -1)*HR3(-1,-1,1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR3(0,-1,-1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR3(0,0,-1)
     $  + 5.8224052646501250d-01*HR1(0)*HR1(1)
     $  + HR1(0) *HR1(1)*HR2(0,-1)
     $  - HR1(0) *HR3(-1,-1,1)
     $  - 2.0000000000000000d+00*HR1(0)*HR3(0,-1,-1)
     $  - HR1(0) *HR3(0,-1,1)
     $  - HR1(0) *HR3(0,1,-1)
     $  + 1.4780047665430420d+00*HR1(1)
     $  + 6.9314718055994530d-01*HR1(1)*HR2(0,-1)
     $  - 2.0000000000000000d+00*HR1(1)*HR3(0,-1,-1)
     $  - 2.0000000000000000d+00*HR1(1)*HR3(0,0,-1)
     $  + 5.8224052646501250d-01*HR2(-1,1)
     $  - HR2( -1,1)*HR2(0,-1)
     $  - 5.8224052646501250d-01*HR2(0,-1)
     $  + 5.0000000000000000d-01*HR2(0,-1)*HR2(0,-1)
     $  + HR2(0, -1)*HR2(0,1)
     $  - 5.8224052646501250d-01*HR2(0,1)
     $  - 6.9314718055994530d-01*HR3(-1,-1,1)
     $  - 1.3862943611198906d+00*HR3(0,-1,-1)
     $  - 6.9314718055994530d-01*HR3(0,-1,1)
     $  - 6.9314718055994530d-01*HR3(0,1,-1)
     $  - HR4( -1,-1,-1,1)
     $  + 4.0000000000000000d+00*HR4(0,-1,-1,-1)
     $  + 2.0000000000000000d+00*HR4(0,-1,-1,1)
     $  - HR4(0, -1,0,1)
     $  + HR4(0, -1,1,-1)
     $  + 2.0000000000000000d+00*HR4(0,0,-1,-1)
     $  + HR4(0,1, -1,-1)
      HY4(0,1,1,-1) = 
     $  + 4.3369237704895519d-01
     $  - 1.1073038989294665d+00*HR1(-1)
     $  + 5.3134677019160696d-01*HR1(-1)*HR1(-1)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  - 1.0626935403832139d+00*HR1(-1)*HR1(0)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(0)*HR1(0)
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(0)*HR1(1)
     $  + 1.0626935403832139d+00*HR1(-1)*HR1(1)
     $  - HR1( -1)*HR1(1)*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR2(-1,1)
     $  + HR1( -1)*HR3(-1,-1,1)
     $  + HR1( -1)*HR3(0,-1,-1)
     $  + HR1( -1)*HR3(0,0,-1)
     $  - 3.4657359027997265d-01*HR1(0)*HR1(0)*HR1(1)
     $  - 1.0626935403832139d+00*HR1(0)*HR1(1)
     $  - 6.9314718055994530d-01*HR1(0)*HR2(-1,1)
     $  + 6.9314718055994530d-01*HR1(0)*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR1(0)*HR2(0,1)
     $  - 1.1073038989294665d+00*HR1(1)
     $  + HR1(1) *HR3(0,-1,-1)
     $  + HR1(1) *HR3(0,0,-1)
     $  - 1.0626935403832139d+00*HR2(-1,1)
     $  + HR2( -1,1)*HR2(0,-1)
     $  + 1.0626935403832139d+00*HR2(0,-1)
     $  - 5.0000000000000000d-01*HR2(0,-1)*HR2(0,-1)
     $  - HR2(0, -1)*HR2(0,1)
     $  + 1.0626935403832139d+00*HR2(0,1)
     $  - 6.9314718055994530d-01*HR3(-1,-1,1)
     $  - 6.9314718055994530d-01*HR3(0,-1,-1)
     $  - 6.9314718055994530d-01*HR3(0,0,-1)
     $  - 6.9314718055994530d-01*HR3(0,0,1)
     $  - 6.9314718055994530d-01*HR3(0,1,-1)
     $  - HR4( -1,-1,-1,1)
     $  - HR4(0, -1,-1,1)
     $  + HR4(0, -1,0,1)
     $  + HR4(0,0, -1,1)
     $  + HR4(0,0,1, -1)
     $  + HR4(0,1, -1,-1)
      HY4(-1,-1,-1,1) = 
     $  + 1.4134237214990008d-02
     $  - 9.4753004230127705d-02*HR1(-1)
     $  + 2.9112026323250625d-01*HR1(-1)*HR1(-1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  + HR1( -1)*HR3(0,-1,-1)
     $  - HR4(0, -1,-1,-1)
      HY4(-1,-1,1,1) = 
     $  + 4.0758239159309251d-02
     $  - 5.3721319360804020d-01*HR1(-1)
     $  + 1.2011325347955035d-01*HR1(-1)*HR1(-1)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  + 2.5000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(0)
     $  - HR1( -1)*HR1(0)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR2(0,-1)
     $  + HR1( -1)*HR3(0,-1,-1)
     $  + HR1( -1)*HR3(0,0,-1)
     $  + HR1(0) *HR3(0,-1,-1)
     $  + 6.9314718055994530d-01*HR3(0,-1,-1)
     $  - 2.0000000000000000d+00*HR4(0,-1,-1,-1)
     $  - HR4(0,0, -1,-1)
      HY4(-1,1,1,1) = 
     $  + 5.1747906167389938d-01
     $  + 5.5504108664821579d-02*HR1(-1)
     $  - 1.2011325347955035d-01*HR1(-1)*HR1(-1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  - 2.5000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(0)
     $  + 2.4022650695910071d-01*HR1(-1)*HR1(0)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(0)*HR1(0)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(0)*HR1(0)*HR1(0)
     $  - 5.0000000000000000d-01*HR1(0)*HR1(0)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(0)*HR2(0,-1)
     $  + HR1(0) *HR3(0,-1,-1)
     $  + HR1(0) *HR3(0,0,-1)
     $  - 2.4022650695910071d-01*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR3(0,-1,-1)
     $  + 6.9314718055994530d-01*HR3(0,0,-1)
     $  - HR4(0, -1,-1,-1)
     $  - HR4(0,0, -1,-1)
     $  - HR4(0,0,0, -1)
      if (r.lt.0d0) then 
      HY4(0,-1,1,1) = HY4(0,-1,1,1) 
     $  - 2.4674011002723396d+00*HR1(-1)*HR1(-1)
     $  - 4.9348022005446793d+00*HR1(-1)*HR1(1)
     $  + 4.9348022005446793d+00*HR2(-1,1)
      HY4(0,1,1,-1) = HY4(0,1,1,-1) 
     $  + 3.4205442319285582d+00*HR1(-1)
     $  + 3.4205442319285582d+00*HR1(1)
      HY4(-1,-1,1,1) = HY4(-1,-1,1,1) 
     $  - 2.4674011002723396d+00*HR1(-1)*HR1(-1) 
      HY4(-1,1,1,1) = HY4(-1,1,1,1)
     $  - 3.4205442319285582d+00*HR1(-1)
     $  + 2.4674011002723396d+00*HR1(-1)*HR1(-1)
     $  - 4.9348022005446793d+00*HR1(-1)*HR1(0)
     $  + 4.9348022005446793d+00*HR2(0,-1)
      Hi4(0,0,-1,1) = 
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1) 
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1) 
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(1)*HR1(1) 
     $  + HR1(1) *HR2(-1,1) 
     $  + HR3( -1,-1,1) 
     $  - HR3( -1,1,1) 
      Hi4(0,0,1,-1) = 
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1) 
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(1) 
     $  + 3.4657359027997265d-01*HR1(1)*HR1(1) 
      Hi4(0,-1,0,1) = 
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1) 
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1) 
     $  + HR1( -1)*HR2(-1,1) 
     $  - HR1(1) *HR2(-1,1) 
     $  - 2.0000000000000000d+00*HR3(-1,-1,1) 
     $  + 2.0000000000000000d+00*HR3(-1,1,1) 
      Hi4(0,-1,-1,1) = 
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1) 
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1) 
     $  + HR1( -1)*HR2(-1,1) 
     $  - HR3( -1,-1,1) 
      Hi4(0,-1,1,-1) = 
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1) 
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(1) 
     $  - 6.9314718055994530d-01*HR2(-1,1) 
      Hi4(0,1,-1,-1) =  
     $  - 2.4022650695910071d-01*HR1(-1) 
     $  - 2.4022650695910071d-01*HR1(1) 
      Hi4(0,-1,1,1) = 
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - HR1( -1)*HR1(0)*HR1(1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(1)
     $  - HR1( -1)*HR2(-1,1)
     $  + HR1( -1)*HR2(0,-1)
     $  + HR1(0) *HR2(-1,1)
     $  + HR1(1) *HR2(0,-1)
     $  + 6.9314718055994530d-01*HR2(-1,1)
     $  + HR3( -1,-1,1)
     $  - HR3(0, -1,-1)
     $  - HR3(0, -1,1)
      Hi4(0,1,-1,1) = 
     $  - 5.8224052646501250d-01*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - HR1( -1)*HR2(-1,1)
     $  - HR1( -1)*HR2(0,-1)
     $  - 5.8224052646501250d-01*HR1(1)
     $  - HR1(1) *HR2(0,-1)
     $  + HR3( -1,-1,1)
     $  + 2.0000000000000000d+00*HR3(0,-1,-1)
     $  + HR3(0, -1,1)
     $  + HR3(0,1, -1)
      Hi4(0,1,1,-1) = 
     $  + 1.0626935403832139d+00*HR1(-1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(0)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(1)
     $  + 6.9314718055994530d-01*HR1(0)*HR1(1)
     $  + 1.0626935403832139d+00*HR1(1)
     $  + 6.9314718055994530d-01*HR2(-1,1)
     $  - 6.9314718055994530d-01*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR2(0,1)
      Hi4(-1,-1,-1,1) = 
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1) 
      Hi4(-1,-1,1,1) = 
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  + HR1( -1)*HR2(0,-1)
     $  - HR3(0, -1,-1)
      Hi4(-1,1,1,1) = 
     $  + 1.4047075598891257d+00*HR1(-1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(0)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(0)*HR1(0)
     $  + HR1(0) *HR2(0,-1)
     $  + 6.9314718055994530d-01*HR2(0,-1)
     $  - HR3(0, -1,-1)
     $  - HR3(0,0, -1)
      endif  
      endif 
** nw > 3 endif 
      endif 
** (n1,n2) = (-1,1) -- completion endif 
      return 
      end 
************************************************************************ 
      subroutine fillirr1dhplatinf(x,nw,HX1,HX2,HX3,HX4, 
     $                                HY1,HY2,HY3,HY4, 
     $                                Hi1,Hi2,Hi3,Hi4,n1,n2) 
** evaluates the HPL for y > r2p1
** fillirr1dhplatinf is called by eval1dhplatinf after calling 
** fillirr1dhplat0 with argument r=1/y 
** it is guaranteed that nw is in the range 2:4, and that (n1,n2) 
** take one of the pairs of values (0,1), (-1,0) or (-1,1) 
      implicit double precision (a-h,o-z) 
      dimension HX1(n1:n2),HX2(n1:n2,n1:n2),HX3(n1:n2,n1:n2,n1:n2), 
     $          HX4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension HY1(n1:n2),HY2(n1:n2,n1:n2),HY3(n1:n2,n1:n2,n1:n2), 
     $          HY4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2) 
** (n1,n2) = (0,1) or (-1,1) 
      if (    ( (n1.eq.0).and.(n2.eq.1) ) 
     $    .or.( (n1.eq.-1).and.(n2.eq.1) ) ) then 
      HY2(0,1) = 
     $  + 3.2898681336964528d+00 
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0) 
     $  - HX2(0,1) 
      Hi2(0,1) = 
     $  - HX1(0) 
      if ( nw.gt.2 ) then 
      HY3(0,0,1) = 
     $  - 3.2898681336964528d+00*HX1(0) 
     $  + 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0) 
     $  + HX3(0,0,1) 
      HY3(0,1,1) = 
     $  + 1.2020569031595942d+00 
     $  + 4.9348022005446793d+00*HX1(0) 
     $  - 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0) 
     $  - HX1(0) *HX2(0,1) 
     $  + HX3(0,0,1) 
     $  - HX3(0,1,1) 
      Hi3(0,0,1) = 
     $  + 5.000000000000000d-01*HX1(0)*HX1(0) 
      Hi3(0,1,1) = 
     $  + 1.6449340668482264d+00 
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0) 
     $  - HX2(0,1) 
      endif 
      if ( nw.gt.3 ) then 
      HY4(0,0,0,1) = 
     $  + 2.1646464674222763d+00 
     $  + 1.6449340668482264d+00*HX1(0)*HX1(0) 
     $  - 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0) 
     $  - HX4(0,0,0,1) 
      HY4(0,0,1,1) = 
     $  + 2.1646464674222763d+00 
     $  - 1.2020569031595942d+00*HX1(0) 
     $  - 2.4674011002723396d+00*HX1(0)*HX1(0) 
     $  + 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0) 
     $  + HX1(0) *HX3(0,0,1) 
     $  - 2.0000000000000000d+00*HX4(0,0,0,1) 
     $  + HX4(0,0,1,1) 
      HY4(0,1,1,1) = 
     $  - 5.1410353601279064d+00 
     $  + 2.4674011002723396d+00*HX1(0)*HX1(0) 
     $  - 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0) 
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0)*HX2(0,1) 
     $  + HX1(0) *HX3(0,0,1) 
     $  - HX1(0) *HX3(0,1,1) 
     $  + 4.9348022005446793d+00*HX2(0,1) 
     $  - HX4(0,0,0,1) 
     $  + HX4(0,0,1,1) 
     $  - HX4(0,1,1,1) 
      Hi4(0,0,0,1) = 
     $  - 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0) 
      Hi4(0,0,1,1) = 
     $  - 1.2020569031595942d+00 
     $  - 1.6449340668482264d+00*HX1(0) 
     $  + 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0) 
     $  + HX3(0,0,1) 
      Hi4(0,1,1,1) = 
     $  + 1.6449340668482264d+00*HX1(0) 
     $  - 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0) 
     $  - HX1(0) *HX2(0,1) 
     $  + HX3(0,0,1) 
     $  - HX3(0,1,1) 
      endif 
** nw > 3 endif 
      endif 
** (n1,n2) = (0,1) or (-1,1) endif 
************ 
** (n1,n2) = (-1,0) or (-1,1) 
      if (    ( (n1.eq.-1).and.(n2.eq.0) ) 
     $    .or.( (n1.eq.-1).and.(n2.eq.1) ) ) then 
      HY2(0,-1) = 
     $  + 1.6449340668482264d+00 
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0) 
     $  - HX2(0, -1) 
      if ( nw.gt.2 ) then 
      HY3(0,0,-1) = 
     $  - 1.6449340668482264d+00*HX1(0) 
     $  - 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0) 
     $  + HX3(0,0, -1) 
      HY3(0,-1,-1) = 
     $  + 1.2020569031595942d+00 
     $  - 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0) 
     $  + HX1(0) *HX2(0,-1) 
     $  - HX3(0, -1,-1) 
     $  - HX3(0,0, -1) 
      endif 
      if ( nw.gt.3 ) then 
      HY4(0,0,0,-1) = 
     $  + 1.8940656589944918d+00 
     $  + 8.2246703342411321d-01*HX1(0)*HX1(0) 
     $  + 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0) 
     $  - HX4(0,0,0, -1) 
      HY4(0,0,-1,-1) = 
     $  - 1.8940656589944918d+00 
     $  - 1.2020569031595942d+00*HX1(0) 
     $  + 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0) 
     $  - HX1(0) *HX3(0,0,-1) 
     $  + HX4(0,0, -1,-1) 
     $  + 2.0000000000000000d+00*HX4(0,0,0,-1) 
      HY4(0,-1,-1,-1) = 
     $  + 1.0823232337111381d+00 
     $  + 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0) 
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0)*HX2(0,-1) 
     $  + HX1(0) *HX3(0,-1,-1) 
     $  + HX1(0) *HX3(0,0,-1) 
     $  - HX4(0, -1,-1,-1) 
     $  - HX4(0,0, -1,-1) 
     $  - HX4(0,0,0, -1) 
      endif 
** nw > 3 endif 
      endif 
** (n1,n2) = (-1,0) or (-1,1) endif 
** (n1,n2) = (-1,1) -- completion 
      if ( (n1.eq.-1).and.(n2.eq.1) ) then 
      HY2(-1,1) = 
     $  + 2.4674011002723396d+00 
     $  + HX1( -1)*HX1(0) 
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0) 
     $  + HX2( -1,1) 
     $  - HX2(0, -1) 
     $  - HX2(0,1) 
      Hi2(-1,1) = 
     $  - 6.9314718055994530d-01 
     $  + HX1( -1) 
     $  - HX1(0) 
      if ( nw.gt.2 ) then 
      HY3(0,-1,1) = 
     $  - 2.5190015545588625d+00 
     $  - 2.4674011002723396d+00*HX1(0) 
     $  + 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0) 
     $  - HX1(0) *HX2(0,-1) 
     $  - HX3(0, -1,1) 
     $  + 2.0000000000000000d+00*HX3(0,0,-1) 
     $  + HX3(0,0,1) 
      HY3(0,1,-1) = 
     $  + 4.3220869092982539d+00 
     $  + 2.4674011002723396d+00*HX1(0) 
     $  + 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0) 
     $  + HX1(0) *HX2(0,1) 
     $  - HX3(0,0, -1) 
     $  - 2.0000000000000000d+00*HX3(0,0,1) 
     $  - HX3(0,1, -1) 
      HY3(-1,-1,1) = 
     $  - 2.7620719062289241d+00 
     $  + 2.4674011002723396d+00*HX1(-1) 
     $  + 5.0000000000000000d-01*HX1(-1)*HX1(-1)*HX1(0) 
     $  - 5.0000000000000000d-01*HX1(-1)*HX1(0)*HX1(0) 
     $  - HX1( -1)*HX2(0,-1) 
     $  - HX1( -1)*HX2(0,1) 
     $  - 2.4674011002723396d+00*HX1(0) 
     $  + 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0) 
     $  + HX3( -1,-1,1) 
     $  + HX3(0, -1,-1) 
     $  + HX3(0,0, -1) 
     $  + HX3(0,0,1) 
     $  + HX3(0,1, -1) 
      HY3(-1,1,1) = 
     $  + 2.7620719062289241d+00 
     $  - 4.9348022005446793d+00*HX1(-1) 
     $  + 5.0000000000000000d-01*HX1(-1)*HX1(0)*HX1(0) 
     $  + 4.9348022005446793d+00*HX1(0) 
     $  - 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0) 
     $  + HX1(0) *HX2(-1,1) 
     $  - HX1(0) *HX2(0,-1) 
     $  - HX1(0) *HX2(0,1) 
     $  + HX3( -1,1,1) 
     $  - HX3(0, -1,1) 
     $  + HX3(0,0, -1) 
     $  + HX3(0,0,1) 
     $  - HX3(0,1,1) 
      Hi3(0,-1,1) = 
     $  + 8.2246703342411321d-01 
     $  + 6.9314718055994530d-01*HX1(0) 
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0) 
     $  - HX2(0, -1) 
      Hi3(0,1,-1) = 
     $  - 6.9314718055994530d-01*HX1(0) 
      Hi3(-1,-1,1) = 
     $  + 2.4022650695910071d-01 
     $  - 6.9314718055994530d-01*HX1(-1) 
     $  + 5.0000000000000000d-01*HX1(-1)*HX1(-1) 
     $  - HX1( -1)*HX1(0) 
     $  + 6.9314718055994530d-01*HX1(0) 
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0) 
      Hi3(-1,1,1) = 
     $  + 1.8851605738073271d+00 
     $  + HX1( -1)*HX1(0) 
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0) 
     $  + HX2( -1,1) 
     $  - HX2(0, -1) 
     $  - HX2(0,1) 
      endif 
      if ( nw.gt.3 ) then 
      HY4(0,0,-1,1) = 
     $  + 3.9234217222028759d+00
     $  + 2.5190015545588625d+00*HX1(0)
     $  + 1.2337005501361698d+00*HX1(0)*HX1(0)
     $  - 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + HX1(0) *HX3(0,0,-1)
     $  + HX4(0,0, -1,1)
     $  - 3.0000000000000000d+00*HX4(0,0,0,-1)
     $  - HX4(0,0,0,1) 
      HY4(0,0,1,-1) = 
     $  - 4.1940025306306604d+00
     $  - 4.3220869092982539d+00*HX1(0)
     $  - 1.2337005501361698d+00*HX1(0)*HX1(0)
     $  - 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  - HX1(0) *HX3(0,0,1)
     $  + HX4(0,0,0, -1)
     $  + 3.0000000000000000d+00*HX4(0,0,0,1)
     $  + HX4(0,0,1, -1)
      HY4(0,-1,0,1) = 
     $  + 9.4703282949724591d-01
     $  + 1.8030853547393914d+00*HX1(0)
     $  + 1.6449340668482264d+00*HX1(0)*HX1(0)
     $  - 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0)*HX2(0,-1)
     $  - 2.0000000000000000d+00*HX1(0)*HX3(0,0,-1)
     $  - 3.2898681336964528d+00*HX2(0,-1)
     $  + HX4(0, -1,0,1)
     $  + 3.0000000000000000d+00*HX4(0,0,0,-1)
     $  - HX4(0,0,0,1) 
      HY4(0,-1,-1,1) = 
     $  + 2.5209599327464717d+00
     $  + 2.7620719062289241d+00*HX1(0)
     $  + 1.2337005501361698d+00*HX1(0)*HX1(0)
     $  - 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0)*HX2(0,-1)
     $  - HX1(0) *HX3(0,-1,-1)
     $  - HX1(0) *HX3(0,0,-1)
     $  - 2.4674011002723396d+00*HX2(0,-1)
     $  + 5.0000000000000000d-01*HX2(0,-1)*HX2(0,-1)
     $  - HX4(0, -1,-1,1)
     $  + HX4(0, -1,0,1)
     $  + HX4(0,0, -1,1)
     $  - HX4(0,0,0,1) 
      HY4(0,-1,1,-1) = 
     $  - 8.5266539820739622d+00
     $  - 5.5241438124578482d+00*HX1(0)
     $  - 1.2337005501361698d+00*HX1(0)*HX1(0)
     $  - 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0)*HX2(0,-1)
     $  + HX1(0) *HX3(0,-1,1)
     $  - 2.0000000000000000d+00*HX1(0)*HX3(0,0,-1)
     $  - HX1(0) *HX3(0,0,1)
     $  + 2.4674011002723396d+00*HX2(0,-1)
     $  - 5.0000000000000000d-01*HX2(0,-1)*HX2(0,-1)
     $  - HX4(0, -1,0,1)
     $  - HX4(0, -1,1,-1)
     $  + 2.0000000000000000d+00*HX4(0,0,-1,-1)
     $  - 2.0000000000000000d+00*HX4(0,0,-1,1)
     $  + 4.0000000000000000d+00*HX4(0,0,0,-1)
     $  + 3.0000000000000000d+00*HX4(0,0,0,1)
     $  + HX4(0,0,1, -1)
      HY4(0,1,-1,-1) = 
     $  + 5.8027584430066521d+00
     $  + 2.7620719062289241d+00*HX1(0)
     $  - 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0)*HX2(0,1)
     $  + HX1(0) *HX3(0,0,-1)
     $  + 2.0000000000000000d+00*HX1(0)*HX3(0,0,1)
     $  + HX1(0) *HX3(0,1,-1)
     $  - HX4(0,0, -1,-1)
     $  - 2.0000000000000000d+00*HX4(0,0,0,-1)
     $  - 3.0000000000000000d+00*HX4(0,0,0,1)
     $  - 2.0000000000000000d+00*HX4(0,0,1,-1)
     $  - HX4(0,1, -1,-1)
      HY4(0,-1,1,1) = 
     $  + 6.2689427375197987d-01
     $  - 2.7620719062289241d+00*HX1(0)
     $  - 2.4674011002723396d+00*HX1(0)*HX1(0)
     $  + 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0)*HX2(0,-1)
     $  - HX1(0) *HX3(0,-1,1)
     $  + 2.0000000000000000d+00*HX1(0)*HX3(0,0,-1)
     $  + HX1(0) *HX3(0,0,1)
     $  + 4.9348022005446793d+00*HX2(0,-1)
     $  - HX4(0, -1,1,1)
     $  + 2.0000000000000000d+00*HX4(0,0,-1,1)
     $  - 3.0000000000000000d+00*HX4(0,0,0,-1)
     $  - 2.0000000000000000d+00*HX4(0,0,0,1)
     $  + HX4(0,0,1,1) 
      HY4(0,1,-1,1) = 
     $  - 4.3326514514433017d+00
     $  - 1.3169446513992682d+00*HX1(0)
     $  - 1.2337005501361698d+00*HX1(0)*HX1(0)
     $  + 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0)*HX2(0,1)
     $  - HX1(0) *HX3(0,0,-1)
     $  - 2.0000000000000000d+00*HX1(0)*HX3(0,0,1)
     $  - HX1(0) *HX3(0,1,-1)
     $  + HX2(0, -1)*HX2(0,1)
     $  - 2.4674011002723396d+00*HX2(0,1)
     $  + 5.0000000000000000d-01*HX2(0,1)*HX2(0,1)
     $  - HX4(0, -1,0,1)
     $  - 3.0000000000000000d+00*HX4(0,0,-1,1)
     $  + 3.0000000000000000d+00*HX4(0,0,0,-1)
     $  + 4.0000000000000000d+00*HX4(0,0,0,1)
     $  - 2.0000000000000000d+00*HX4(0,0,1,1)
     $  - HX4(0,1, -1,1)
      HY4(0,1,1,-1) = 
     $  - 1.5001934240460787d-01
     $  + 4.0790165576281924d+00*HX1(0)
     $  + 1.2337005501361698d+00*HX1(0)*HX1(0)
     $  + 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0)*HX2(0,1)
     $  - HX1(0) *HX3(0,0,1)
     $  + HX1(0) *HX3(0,1,1)
     $  - HX2(0, -1)*HX2(0,1)
     $  + 2.4674011002723396d+00*HX2(0,1)
     $  - 5.0000000000000000d-01*HX2(0,1)*HX2(0,1)
     $  + HX4(0, -1,0,1)
     $  + 2.0000000000000000d+00*HX4(0,0,-1,1)
     $  - HX4(0,0,0, -1)
     $  + HX4(0,0,1, -1)
     $  - HX4(0,1,1, -1)
      HY4(-1,-1,-1,1) = 
     $  + 2.4278628067547031d+00
     $  - 2.7620719062289241d+00*HX1(-1)
     $  + 1.2337005501361698d+00*HX1(-1)*HX1(-1)
     $  + 1.6666666666666666d-01*HX1(-1)*HX1(-1)*HX1(-1)*HX1(0)
     $  - 2.5000000000000000d-01*HX1(-1)*HX1(-1)*HX1(0)*HX1(0)
     $  - 5.0000000000000000d-01*HX1(-1)*HX1(-1)*HX2(0,-1)
     $  - 5.0000000000000000d-01*HX1(-1)*HX1(-1)*HX2(0,1)
     $  - 2.4674011002723396d+00*HX1(-1)*HX1(0)
     $  + 1.6666666666666666d-01*HX1(-1)*HX1(0)*HX1(0)*HX1(0)
     $  + HX1( -1)*HX3(0,-1,-1)
     $  + HX1( -1)*HX3(0,0,-1)
     $  + HX1( -1)*HX3(0,0,1)
     $  + HX1( -1)*HX3(0,1,-1)
     $  + 2.7620719062289241d+00*HX1(0)
     $  + 1.2337005501361698d+00*HX1(0)*HX1(0)
     $  - 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + HX4( -1,-1,-1,1)
     $  - HX4(0, -1,-1,-1)
     $  - HX4(0,0, -1,-1)
     $  - HX4(0,0,0, -1)
     $  - HX4(0,0,0,1) 
     $  - HX4(0,0,1, -1)
     $  - HX4(0,1, -1,-1)
      HY4(-1,-1,1,1) = 
     $  + 2.0293560632083841d+00
     $  + 2.7620719062289241d+00*HX1(-1)
     $  - 2.4674011002723396d+00*HX1(-1)*HX1(-1)
     $  + 2.5000000000000000d-01*HX1(-1)*HX1(-1)*HX1(0)*HX1(0)
     $  + 4.9348022005446793d+00*HX1(-1)*HX1(0)
     $  - 1.6666666666666666d-01*HX1(-1)*HX1(0)*HX1(0)*HX1(0)
     $  - HX1( -1)*HX1(0)*HX2(0,-1)
     $  - HX1( -1)*HX1(0)*HX2(0,1)
     $  - HX1( -1)*HX3(0,-1,1)
     $  + HX1( -1)*HX3(0,0,-1)
     $  + HX1( -1)*HX3(0,0,1)
     $  - HX1( -1)*HX3(0,1,1)
     $  - 2.7620719062289241d+00*HX1(0)
     $  - 2.4674011002723396d+00*HX1(0)*HX1(0)
     $  + 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + HX1(0) *HX3(-1,-1,1)
     $  + HX1(0) *HX3(0,-1,-1)
     $  + HX1(0) *HX3(0,0,-1)
     $  + HX1(0) *HX3(0,0,1)
     $  + HX1(0) *HX3(0,1,-1)
     $  + HX4( -1,-1,1,1)
     $  + HX4(0, -1,-1,1)
     $  + HX4(0, -1,1,-1)
     $  - HX4(0,0, -1,-1)
     $  + HX4(0,0, -1,1)
     $  - 2.0000000000000000d+00*HX4(0,0,0,-1)
     $  - 2.0000000000000000d+00*HX4(0,0,0,1)
     $  - HX4(0,0,1, -1)
     $  + HX4(0,0,1,1) 
     $  + HX4(0,1, -1,1)
     $  + HX4(0,1,1, -1)
      HY4(-1,1,1,1) = 
     $  - 6.4865749331714713d+00
     $  - 4.9348022005446793d+00*HX1(-1)*HX1(0)
     $  + 1.6666666666666666d-01*HX1(-1)*HX1(0)*HX1(0)*HX1(0)
     $  + 2.4674011002723396d+00*HX1(0)*HX1(0)
     $  - 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0)*HX2(-1,1)
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0)*HX2(0,-1)
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0)*HX2(0,1)
     $  + HX1(0) *HX3(-1,1,1)
     $  - HX1(0) *HX3(0,-1,1)
     $  + HX1(0) *HX3(0,0,-1)
     $  + HX1(0) *HX3(0,0,1)
     $  - HX1(0) *HX3(0,1,1)
     $  - 4.9348022005446793d+00*HX2(-1,1)
     $  + 4.9348022005446793d+00*HX2(0,-1)
     $  + 4.9348022005446793d+00*HX2(0,1)
     $  + HX4( -1,1,1,1)
     $  - HX4(0, -1,1,1)
     $  + HX4(0,0, -1,1)
     $  - HX4(0,0,0, -1)
     $  - HX4(0,0,0,1) 
     $  + HX4(0,0,1,1) 
     $  - HX4(0,1,1,1) 
      Hi4(0,0,-1,1) = 
     $  - 9.0154267736969571d-01 
     $  - 8.2246703342411321d-01*HX1(0) 
     $  - 3.4657359027997265d-01*HX1(0)*HX1(0) 
     $  - 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0) 
     $  + HX3(0,0, -1) 
      Hi4(0,0,1,-1) = 
     $  + 3.4657359027997265d-01*HX1(0)*HX1(0) 
      Hi4(0,-1,0,1) = 
     $  + 1.8030853547393914d+00 
     $  + 8.2246703342411321d-01*HX1(0) 
     $  - 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0) 
     $  + HX1(0) *HX2(0,-1) 
     $  - 2.0000000000000000d+00*HX3(0,0,-1) 
      Hi4(0,-1,-1,1) = 
     $  + 4.8170908494321862d-01 
     $  - 2.4022650695910071d-01*HX1(0) 
     $  - 3.4657359027997265d-01*HX1(0)*HX1(0) 
     $  - 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0) 
     $  + HX1(0) *HX2(0,-1) 
     $  + 6.9314718055994530d-01*HX2(0,-1) 
     $  - HX3(0, -1,-1) 
     $  - HX3(0,0, -1) 
      Hi4(0,-1,1,-1) = 
     $  + 5.7009070532142637d-01 
     $  + 4.8045301391820142d-01*HX1(0) 
     $  + 3.4657359027997265d-01*HX1(0)*HX1(0) 
     $  - 6.9314718055994530d-01*HX2(0,-1) 
      Hi4(0,1,-1,-1) = 
     $  - 2.4022650695910071d-01*HX1(0) 
      Hi4(0,-1,1,1) = 
     $  - 2.7620719062289241d+00 
     $  - 1.8851605738073271d+00*HX1(0) 
     $  + 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0) 
     $  - HX1(0) *HX2(0,-1) 
     $  - HX3(0, -1,1) 
     $  + 2.0000000000000000d+00*HX3(0,0,-1) 
     $  + HX3(0,0,1) 
      Hi4(0,1,-1,1) = 
     $  + 2.6736902858507163d+00 
     $  + 1.3029200473423146d+00*HX1(0) 
     $  + 3.4657359027997265d-01*HX1(0)*HX1(0) 
     $  + 1.6666666666666665d-01*HX1(0)*HX1(0)*HX1(0) 
     $  + HX1(0) *HX2(0,1) 
     $  + 6.9314718055994530d-01*HX2(0,1) 
     $  - HX3(0,0, -1) 
     $  - 2.0000000000000000d+00*HX3(0,0,1) 
     $  - HX3(0,1, -1) 
      Hi4(0,1,1,-1) =  
     $  + 1.1401814106428527d+00 
     $  + 5.8224052646501250d-01*HX1(0) 
     $  - 3.4657359027997265d-01*HX1(0)*HX1(0) 
     $  - 6.9314718055994530d-01*HX2(0,1) 
      Hi4(-1,-1,-1,1) = 
     $  - 5.5504108664821579d-02
     $  + 2.4022650695910071d-01*HX1(-1)
     $  - 3.4657359027997265d-01*HX1(-1)*HX1(-1)
     $  + 1.6666666666666666d-01*HX1(-1)*HX1(-1)*HX1(-1)
     $  - 5.0000000000000000d-01*HX1(-1)*HX1(-1)*HX1(0)
     $  + 6.9314718055994530d-01*HX1(-1)*HX1(0)
     $  + 5.0000000000000000d-01*HX1(-1)*HX1(0)*HX1(0)
     $  - 2.4022650695910071d-01*HX1(0)
     $  - 3.4657359027997265d-01*HX1(0)*HX1(0)
     $  - 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0)
      Hi4(-1,-1,1,1) = 
     $  - 2.4532465311320902d+00
     $  + 1.8851605738073271d+00*HX1(-1)
     $  + 5.0000000000000000d-01*HX1(-1)*HX1(-1)*HX1(0)
     $  - 5.0000000000000000d-01*HX1(-1)*HX1(0)*HX1(0)
     $  - HX1( -1)*HX2(0,-1)
     $  - HX1( -1)*HX2(0,1)
     $  - 1.8851605738073271d+00*HX1(0)
     $  + 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0)
     $  + HX3( -1,-1,1)
     $  + HX3(0, -1,-1)
     $  + HX3(0,0, -1)
     $  + HX3(0,0,1) 
     $  + HX3(0,1, -1)
      Hi4(-1,1,1,1) = 
     $  - 5.5504108664821579d-02
     $  - 1.6449340668482264d+00*HX1(-1)
     $  + 5.0000000000000000d-01*HX1(-1)*HX1(0)*HX1(0)
     $  + 1.6449340668482264d+00*HX1(0)
     $  - 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0)
     $  + HX1(0) *HX2(-1,1)
     $  - HX1(0) *HX2(0,-1)
     $  - HX1(0) *HX2(0,1)
     $  + HX3( -1,1,1)
     $  - HX3(0, -1,1)
     $  + HX3(0,0, -1)
     $  + HX3(0,0,1) 
     $  - HX3(0,1,1) 
      endif 
** nw > 3 endif 
      endif 
** (n1,n2) = (-1,1) -- completion endif 
      return 
      end 
************************************************************************ 
      subroutine fillirr1dhplin1(y,nw,HY1,HY2,HY3,HY4,n1,n2) 
** evaluates the irreducible HPL for y =1
** it is guaranteed that nw is in the range 2:4, and that (n1,n2) 
** take one of the pairs of values (0,1), (-1,0) or (-1,1) 
      implicit double precision (a-h,o-z) 
      dimension HY1(n1:n2),HY2(n1:n2,n1:n2),HY3(n1:n2,n1:n2,n1:n2), 
     $          HY4(n1:n2,n1:n2,n1:n2,n1:n2) 
** (n1,n2) = (0,1) or (-1,1) 
      if (    ( (n1.eq.0).and.(n2.eq.1) ) 
     $    .or.( (n1.eq.-1).and.(n2.eq.1) ) ) then 
      HY2(0,1) = 
     $  + 1.6449340668482264d+00
      if (nw.gt.2) then 
      HY3(0,0,1) = 
     $  + 1.2020569031595942d+00
      HY3(0,1,1) = 
     $  + 1.2020569031595942d+00
      endif
      if (nw.gt.3) then 
      HY4(0,0,0,1) = 
     $  + 1.0823232337111381d+00
      HY4(0,0,1,1) = 
     $  + 2.7058080842778454d-01
      HY4(0,1,1,1) = 
     $  + 1.0823232337111381d+00
      endif
      endif
** (n1,n2) = (0,1) or (-1,1) endif 
************ 
** (n1,n2) = (-1,0) or (-1,1) 
      if (    ( (n1.eq.-1).and.(n2.eq.0) ) 
     $    .or.( (n1.eq.-1).and.(n2.eq.1) ) ) then 
      HY2(0,-1) = 
     $  + 8.2246703342411321d-01
      if (nw.gt.2) then
      HY3(0,-1,-1) = 
     $  + 1.5025711289494928d-01
      HY3(0,0,-1) = 
     $  + 9.0154267736969571d-01
      endif
      if (nw.gt.3) then
      HY4(0,-1,-1,-1) = 
     $  + 2.3752366322618485d-02
      HY4(0,0,-1,-1) = 
     $  + 8.7785671568655302d-02
      HY4(0,0,0,-1) = 
     $  + 9.4703282949724591d-01
      endif
      endif 
** (n1,n2) = (-1,0) or (-1,1) endif 
** (n1,n2) = (-1,1) -- completion 
      if ( (n1.eq.-1).and.(n2.eq.1) ) then 
      HY2(-1,1) = 
     $  + 5.8224052646501250d-01
      if (nw.gt.2) then
      HY3(0,-1,1) = 
     $  + 2.4307035167006157d-01
      HY3(0,1,-1) = 
     $  + 5.0821521280468485d-01
      HY3(-1,-1,1) = 
     $  + 9.4753004230127705d-02
      HY3(-1,1,1) = 
     $  + 5.3721319360804020d-01
      endif
      if (nw.gt.3) then 
      HY4(0,0,-1,1) = 
     $  + 1.1787599965050932d-01
      HY4(0,0,1,-1) = 
     $  + 1.7284527823898438d-01
      HY4(0,-1,0,1) = 
     $  + 2.0293560632083841d-01
      HY4(0,-1,-1,1) = 
     $  + 3.4159126166513913d-02
      HY4(0,-1,1,-1) = 
     $  + 5.4653052738263652d-02
      HY4(0,1,-1,-1) = 
     $  + 1.1412342741606084d-01
      HY4(0,-1,1,1) = 
     $  + 9.3097125991768577d-02
      HY4(0,1,-1,1) = 
     $  + 1.9355535381306524d-01
      HY4(0,1,1,-1) = 
     $  + 4.3369237704895519d-01
      HY4(-1,-1,-1,1) = 
     $  + 1.4134237214990008d-02
      HY4(-1,-1,1,1) = 
     $  + 4.0758239159309251d-02
      HY4(-1,1,1,1) = 
     $  + 5.1747906167389938d-01
      endif
      endif
** (n1,n2) = (-1,1) -- completion endif 
      return
      end
