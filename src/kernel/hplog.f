******************************************************************************
**  hplog5: a subroutine for the evaluation of harmonic polylogarithms
**          Version 1.0         15/11/2001
**  upgraded to w=5 from hplog as described in: 
**  T.Gehrmann and E.Remiddi: Numerical Evaluation of the Harmonic 
**                            Polylogarithms up to Weight 4
**                            (hep-ph/0107173; Comp.Phys.Comm. 141 (2001) 296)
**  the harmonic polylogarithms are defined in: 
**  E.Remiddi and J.Vermaseren: Harmonic Polylogarithms
**                            (hep-ph/9905237; Int.J.Mod.Phys. A15 (2000) 725)
**  email:
**  Thomas.Gehrmann@cern.ch and Ettore.Remiddi@bo.infn.it
**
******************************************************************************
c      subroutine hplog5(x,nw,Hc1,Hc2,Hc3,Hc4,Hc5, 
c     $                       Hr1,Hr2,Hr3,Hr4,Hr5,
c     $                       Hi1,Hi2,Hi3,Hi4,Hi5,n1,n2) 
      subroutine hplog(x,nw,Hr1,Hr2,Hr3,Hr4,Hr5,n1,n2)
****** 
** x is the argument of the 1dHPL's (1 dimensional Harmonic PolyLogarithms) 
**   to be evaluated; 
** nw is the maximum weight of the required 1dHPL's; 
**    the maximum allowed value of nw of this implementation is 5; 
** Hc1,Hc2,Hc3,Hc4,Hc5 are the complex*16 values of the 1dHPL; 
**    they must all be supplied in the arguments even if some of them 
**    are not to be evaluated; 
** Hr1,Hr2,Hr3,Hr4,Hr5 are the double precision real parts of 
**    Hc1,Hc2,Hc3,Hc4,Hc5; 
** Hi1,Hi2,Hi3,Hi4,Hi5 are the double precision immaginary parts of 
**    Hc1,Hc2,Hc3,Hc4,Hc5 divided by pi=3.114159.... 
** n1,n2 is the required range of indices, the allowed ranges are 
**    (0,1), (-1,0), (-1,1) ; 
****** 
      implicit double precision (a-h,o-z) 
      complex*16 Hc1,Hc2,Hc3,Hc4,Hc5 
      dimension Hc1(n1:n2),Hc2(n1:n2,n1:n2),Hc3(n1:n2,n1:n2,n1:n2), 
     $          Hc4(n1:n2,n1:n2,n1:n2,n1:n2), 
     $          Hc5(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension Hr1(n1:n2),Hr2(n1:n2,n1:n2),Hr3(n1:n2,n1:n2,n1:n2), 
     $          Hr4(n1:n2,n1:n2,n1:n2,n1:n2), 
     $          Hr5(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2), 
     $          Hi5(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2) 
      common /fillred/infilldim,infill(3) 
      parameter (r2   = 1.4142135623730950488d0) 
** check on the weight nw 
      if ( (nw.lt.1).or.(nw.gt.5) ) then 
        print*, ' illegal call of eval1dhpl with second argument', 
     $          ' (the weight) = ',nw 
        print*, ' the allowed values of the weight are 1,...,5 ' 
        stop
      endif
** check on the range n1:n2 
      if ( (n1.eq.-1).and.(n2.eq.0) ) then 
        infilldim =  2 
        infill(1) =  0 
        infill(2) = -1  
      elseif ( (n1.eq.0).and.(n2.eq.1) ) then 
        infilldim =  2 
        infill(1) =  0 
        infill(2) =  1  
      elseif ( (n1.eq.-1).and.(n2.eq.1) ) then 
        infilldim =  3 
        infill(1) =  0 
        infill(2) = -1  
        infill(3) =  1  
      else 
        print*, ' illegal call of eval1dhpl with the two last ', 
     $          'arguments = (',n1,',',n2,')' 
        print*, ' the allowed values are (-1,0), (0,1), (-1,1) ' 
        stop 
      endif 
** setting the immaginary parts equal to zero 
      call setzero(nw,Hi1,Hi2,Hi3,Hi4,Hi5,n1,n2) 
** looking at the range of the argument 
*      r2 = sqrt(2.d0) 
      r2m1 = r2 - 1 
      r2p1 = r2 + 1 
      if ( ( x.gt.-r2m1 ).and.( x.le.r2m1) ) then 
*       print*, ' eval1dhpl:      x = ',x,', call eval1dhplat0 ' 
        call eval1dhplat0(x,nw,Hc1,Hc2,Hc3,Hc4,Hc5, 
     $                 Hr1,Hr2,Hr3,Hr4,Hr5,Hi1,Hi2,Hi3,Hi4,Hi5,n1,n2) 
        return 
      elseif ( x.eq.1d0 ) then
*       print*, ' eval1dhpl:      x = ',x,', call eval1dhplin1 ' 
        call eval1dhplin1(x,nw,Hc1,Hc2,Hc3,Hc4,Hc5, 
     $                 Hr1,Hr2,Hr3,Hr4,Hr5,Hi1,Hi2,Hi3,Hi4,Hi5,n1,n2) 
        return 
      elseif ( ( x.gt.r2m1 ).and.( x.le.r2p1) ) then 
*       print*, ' eval1dhpl:      x = ',x,', call eval1dhplat1 ' 
        call eval1dhplat1(x,nw,Hc1,Hc2,Hc3,Hc4,Hc5, 
     $                 Hr1,Hr2,Hr3,Hr4,Hr5,Hi1,Hi2,Hi3,Hi4,Hi5,n1,n2) 
        return 
      elseif ( ( x.gt.r2p1 ) ) then 
*       print*, ' eval1dhpl:      x = ',x,', call eval1dhplatinf ' 
        call eval1dhplatinf(x,nw,Hc1,Hc2,Hc3,Hc4,Hc5,
     $                   Hr1,Hr2,Hr3,Hr4,Hr5,Hi1,Hi2,Hi3,Hi4,Hi5,n1,n2) 
        return 
      elseif ( ( x.le.-r2p1) ) then 
*       print*, ' eval1dhpl:      x = ',x,', call eval1dhplatminf ' 
        call eval1dhplatminf(x,nw,Hc1,Hc2,Hc3,Hc4,Hc5, 
     $                    Hr1,Hr2,Hr3,Hr4,Hr5,Hi1,Hi2,Hi3,Hi4,Hi5,n1,n2) 
        return 
      elseif ( x.eq.-1d0 ) then
*       print*, ' eval1dhpl:      x = ',x,', call eval1dhplinm1 ' 
        call eval1dhplinm1(x,nw,Hc1,Hc2,Hc3,Hc4,Hc5,
     $                  Hr1,Hr2,Hr3,Hr4,Hr5,Hi1,Hi2,Hi3,Hi4,Hi5,n1,n2) 
        return 
      elseif ( ( x.gt.-r2p1 ).and.( x.le.-r2m1) ) then 
*       print*, ' eval1dhpl:      x = ',x,', call eval1dhplatm1 ' 
        call eval1dhplatm1(x,nw,Hc1,Hc2,Hc3,Hc4,Hc5, 
     $                  Hr1,Hr2,Hr3,Hr4,Hr5,Hi1,Hi2,Hi3,Hi4,Hi5,n1,n2) 
        return 
      endif 
** 
      end 
************************************************************************ 
      subroutine eval1dhplat0(y,nw,H1,H2,H3,H4,H5, 
     $                          HY1,HY2,HY3,HY4,HY5,
     $                          Hi1,Hi2,Hi3,Hi4,Hi5,n1,n2) 
** evaluates 1dhpl's in the 0-range  -(r2-1) < y <= (r2-1) 
** by direct series expansion (Bernoulli-accelerated) 
      implicit double precision (a-h,o-z) 
      complex*16 H1,H2,H3,H4,H5 
      dimension H1(n1:n2),H2(n1:n2,n1:n2),H3(n1:n2,n1:n2,n1:n2), 
     $          H4(n1:n2,n1:n2,n1:n2,n1:n2), 
     $          H5(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension HY1(n1:n2),HY2(n1:n2,n1:n2),HY3(n1:n2,n1:n2,n1:n2), 
     $          HY4(n1:n2,n1:n2,n1:n2,n1:n2), 
     $          HY5(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2), 
     $          Hi5(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2) 
** evaluate the irreducible 1dHPL's first 
      call fillh1(y,H1,HY1,Hi1,n1,n2) 
      if ( nw.eq.1 ) return 
      call fillirr1dhplat0(y,nw,HY1,HY2,HY3,HY4,HY5,n1,n2) 
** then the reducible 1dHPL's 
      call fillred1dhpl(nw,H1,H2,H3,H4,H5, 
     $              HY1,HY2,HY3,HY4,HY5,Hi1,Hi2,Hi3,Hi4,Hi5,n1,n2) 
      return 
      end 
************************************************************************ 
      subroutine eval1dhplin1(y,nw,H1,H2,H3,H4,H5, 
     $                   HY1,HY2,HY3,HY4,HY5,Hi1,Hi2,Hi3,Hi4,Hi5,n1,n2) 
** evaluates 1dhpl's for y=1 (explicit values are tabulated)
      implicit double precision (a-h,o-z) 
      complex*16 H1,H2,H3,H4,H5 
      dimension H1(n1:n2),H2(n1:n2,n1:n2),H3(n1:n2,n1:n2,n1:n2), 
     $          H4(n1:n2,n1:n2,n1:n2,n1:n2), 
     $          H5(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension HY1(n1:n2),HY2(n1:n2,n1:n2),HY3(n1:n2,n1:n2,n1:n2), 
     $          HY4(n1:n2,n1:n2,n1:n2,n1:n2), 
     $          HY5(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2), 
     $          Hi5(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2) 
      parameter (pi   = 3.14159265358979324d0) 
** evaluate the irreducible 1dHPL's first 
      call fillh1(y,H1,HY1,Hi1,n1,n2) 
      if ( nw.eq.1 ) return 
      call fillirr1dhplin1(y,nw,HY1,HY2,HY3,HY4,HY5,n1,n2) 
** then the reducible 1dHPL's 
      call fillred1dhpl(nw,H1,H2,H3,H4,H5, 
     $                HY1,HY2,HY3,HY4,HY5,Hi1,Hi2,Hi3,Hi4,Hi5,n1,n2) 
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
      if ( nw.eq.4 ) return
      HY5(1,0,0,0,0) = HY5(0,0,0,0,1) 
      Hi5(1,0,0,0,0) = 0d0 
      H5(1,0,0,0,0) = dcmplx(HY5(1,0,0,0,0),Hi5(1,0,0,0,0)*pi)
      return 
      end 
************************************************************************ 
      subroutine eval1dhplat1(y,nw,H1,H2,H3,H4,H5, 
     $                 HY1,HY2,HY3,HY4,HY5,Hi1,Hi2,Hi3,Hi4,Hi5,n1,n2) 
** evaluates 1dhpl's in the 1-range  (r2-1) < y <= (r2+1) 
** evaluating first the H(..,r=(1-y)/(1+y)) by calling eval1dhplat0(r)  
** and then expressing H(..,y=(1-r)/(1+r)) in terms of H(..,r) 
      implicit double precision (a-h,o-z) 
      complex*16 H1,H2,H3,H4,H5 
      dimension H1(n1:n2),H2(n1:n2,n1:n2),H3(n1:n2,n1:n2,n1:n2), 
     $          H4(n1:n2,n1:n2,n1:n2,n1:n2), 
     $          H5(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension HY1(n1:n2),HY2(n1:n2,n1:n2),HY3(n1:n2,n1:n2,n1:n2), 
     $          HY4(n1:n2,n1:n2,n1:n2,n1:n2), 
     $          HY5(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2), 
     $          Hi5(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2) 
** additional arrays required within this routine 
      dimension HR1(-1:1),HR2(-1:1,-1:1),HR3(-1:1,-1:1,-1:1), 
     $          HR4(-1:1,-1:1,-1:1,-1:1), 
     $          HR5(-1:1,-1:1,-1:1,-1:1,-1:1) 
** the nw = 1 case 
      call fillh1(y,H1,HY1,Hi1,n1,n2) 
      if ( nw.eq.1 ) return 
** the nw > 1 case 
      r = (1.d0-y)/(1.d0+y) 
*      print*,' eval1dhplat1: y = ',y,', r = ',r 
** the whole (-1,1) range is in general needed for any pair (n1,n2)
      call fillirr1dhplat0(r,nw,HR1,HR2,HR3,HR4,HR5,-1,1) 
** fillirr1dhplat1 takes care automatically of all the immaginary 
** parts as well as of the jump across y=1 
      call fillirr1dhplat1(r,nw,HR1,HR2,HR3,HR4,HR5, 
     $                          HY1,HY2,HY3,HY4,HY5, 
     $                          Hi1,Hi2,Hi3,Hi4,Hi5,n1,n2) 
** then the reducible 1dHPL's 
      call fillred1dhpl(nw,H1,H2,H3,H4,H5, 
     $                HY1,HY2,HY3,HY4,HY5,Hi1,Hi2,Hi3,Hi4,Hi5,n1,n2) 
      return 
      end 
************************************************************************ 
      subroutine eval1dhplatinf(y,nw,H1,H2,H3,H4,H5, 
     $                  HY1,HY2,HY3,HY4,HY5,Hi1,Hi2,Hi3,Hi4,Hi5,n1,n2) 
** evaluates 1dhpl's in the inf-range  (r2+1) < abs(y) 
** evaluating first the H(..,x=1/y) by calling eval1dhplat0(x)  
** and then expressing H(..,y=1/x) in terms of H(..,x) 
      implicit double precision (a-h,o-z) 
      complex*16 H1,H2,H3,H4,H5 
      dimension H1(n1:n2),H2(n1:n2,n1:n2),H3(n1:n2,n1:n2,n1:n2), 
     $          H4(n1:n2,n1:n2,n1:n2,n1:n2), 
     $          H5(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension HY1(n1:n2),HY2(n1:n2,n1:n2),HY3(n1:n2,n1:n2,n1:n2), 
     $          HY4(n1:n2,n1:n2,n1:n2,n1:n2), 
     $          HY5(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2), 
     $          Hi5(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2) 
** additional arrays required within this routine 
      dimension HX1(n1:n2),HX2(n1:n2,n1:n2),HX3(n1:n2,n1:n2,n1:n2), 
     $          HX4(n1:n2,n1:n2,n1:n2,n1:n2), 
     $          HX5(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2) 
      parameter (pi   = 3.14159265358979324d0) 
** the nw = 1 case 
      call fillh1(y,H1,HY1,Hi1,n1,n2) 
      if ( nw.eq.1 ) return 
** the nw > 1 case 
      x = 1.d0/y 
*      print*,' eval1dhplatinf: y = ',y,', x = ',x 
      call fillirr1dhplat0(x,nw,HX1,HX2,HX3,HX4,HX5,n1,n2) 
** fillirr1dhplatinf takes care automatically of all the immaginary 
** parts as well as of the jump across y=1 
      call fillirr1dhplatinf(x,nw,HX1,HX2,HX3,HX4,HX5, 
     $                            HY1,HY2,HY3,HY4,HY5, 
     $                            Hi1,Hi2,Hi3,Hi4,Hi5,n1,n2) 
** then the reducible 1dHPL's 
      call fillred1dhpl(nw,H1,H2,H3,H4,H5, 
     $                 HY1,HY2,HY3,HY4,HY5,Hi1,Hi2,Hi3,Hi4,Hi5,n1,n2) 
      return 
      end 
************************************************************************ 
      subroutine eval1dhplinm1(y,nw,H1,H2,H3,H4,H5, 
     $                 HY1,HY2,HY3,HY4,HY5,Hi1,Hi2,Hi3,Hi4,Hi5,n1,n2) 
** evaluates 1dhpl's for y=-1 (explicit values are tabulated)
      implicit double precision (a-h,o-z) 
      complex*16 H1,H2,H3,H4,H5 
      complex*16 G1,G2,G3,G4,G5  
      dimension H1(n1:n2),H2(n1:n2,n1:n2),H3(n1:n2,n1:n2,n1:n2), 
     $          H4(n1:n2,n1:n2,n1:n2,n1:n2), 
     $          H5(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension HY1(n1:n2),HY2(n1:n2,n1:n2),HY3(n1:n2,n1:n2,n1:n2), 
     $          HY4(n1:n2,n1:n2,n1:n2,n1:n2), 
     $          HY5(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2), 
     $          Hi5(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2) 
** additional arrays required within this routine 
      dimension G1(-n2:-n1),G2(-n2:-n1,-n2:-n1),
     $          G3(-n2:-n1,-n2:-n1,-n2:-n1), 
     $          G4(-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1), 
     $          G5(-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1) 
      dimension GY1(-n2:-n1),GY2(-n2:-n1,-n2:-n1),
     $          GY3(-n2:-n1,-n2:-n1,-n2:-n1), 
     $          GY4(-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1), 
     $          GY5(-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1) 
      dimension Gi1(-n2:-n1),Gi2(-n2:-n1,-n2:-n1),
     $          Gi3(-n2:-n1,-n2:-n1,-n2:-n1), 
     $          Gi4(-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1), 
     $          Gi5(-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1) 
      common /fillred/infilldim,infill(3) 
      dimension istorfill(3)
      dimension nphase(-1:1) 
      data nphase/-1,1,-1/ 
      parameter (pi   = 3.14159265358979324d0) 
*      print*,' eval1dhplatm1: y = ',y 
      if (infilldim.eq.2) then
         do i=1,2
            istorfill(i) = infill(i)
            infill(i) = -istorfill(i)
         enddo
      endif
** evaluate H(...,-y) 
      call setzero(nw,Gi1,Gi2,Gi3,Gi4,Gi5,-n2,-n1) 
      Gi1(0) = -1
      call eval1dhplin1(-y,nw,G1,G2,G3,G4,G5, 
     $                  GY1,GY2,GY3,GY4,GY5,Gi1,Gi2,Gi3,Gi4,Gi5,-n2,-n1) 
      if (infilldim.eq.2) then
         do i=1,2
            infill(i) = istorfill(i)
         enddo
      endif
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
                if ( nw.gt.4 ) then 
                  do k5=n1,n2 
                    nph5 = nph4*nphase(k5) 
                    HY5(k1,k2,k3,k4,k5) = nph5*GY5(-k1,-k2,-k3,-k4,-k5) 
                    Hi5(k1,k2,k3,k4,k5) =-nph5*Gi5(-k1,-k2,-k3,-k4,-k5) 
                    H5(k1,k2,k3,k4,k5)  = 
     $                dcmplx(HY5(k1,k2,k3,k4,k5),Hi5(k1,k2,k3,k4,k5)*pi) 
                  enddo 
                endif 
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
      if ( nw.eq.4 ) return
      HY5(-1,0,0,0,0) = HY4(0,0,0,0)*HY1(-1)+HY2(0,0)*HY3(0,0,-1)
     $                  + HY5(0,0,0,0,-1)
      Hi5(-1,0,0,0,0) = -HY2(0,-1)*Hi3(0,0,0)-HY4(0,0,0,-1)*Hi1(0)
      H5(-1,0,0,0,0) = dcmplx(HY5(-1,0,0,0,0),Hi5(-1,0,0,0,0)*pi)
      return 
      end 
************************************************************************ 
      subroutine eval1dhplatm1(y,nw,H1,H2,H3,H4,H5, 
     $                    HY1,HY2,HY3,HY4,HY5,Hi1,Hi2,Hi3,Hi4,Hi5,n1,n2) 
** evaluates 1dhpl's in the (-1)-range  -(r2+1) < y <= -(r2-1) 
** evaluating first the H(..,-y) by calling eval1dhplat1(-y), 
** and then expressing H(..,y) in terms of H(..,-y) 
      implicit double precision (a-h,o-z) 
      complex*16 H1,H2,H3,H4,H5 
      complex*16 G1,G2,G3,G4,G5  
      dimension H1(n1:n2),H2(n1:n2,n1:n2),H3(n1:n2,n1:n2,n1:n2), 
     $          H4(n1:n2,n1:n2,n1:n2,n1:n2), 
     $          H5(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension HY1(n1:n2),HY2(n1:n2,n1:n2),HY3(n1:n2,n1:n2,n1:n2), 
     $          HY4(n1:n2,n1:n2,n1:n2,n1:n2), 
     $          HY5(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2), 
     $          Hi5(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2) 
** additional arrays required within this routine 
      dimension G1(-n2:-n1),G2(-n2:-n1,-n2:-n1),
     $          G3(-n2:-n1,-n2:-n1,-n2:-n1), 
     $          G4(-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1), 
     $          G5(-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1) 
      dimension GY1(-n2:-n1),GY2(-n2:-n1,-n2:-n1),
     $          GY3(-n2:-n1,-n2:-n1,-n2:-n1), 
     $          GY4(-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1), 
     $          GY5(-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1) 
      dimension Gi1(-n2:-n1),Gi2(-n2:-n1,-n2:-n1),
     $          Gi3(-n2:-n1,-n2:-n1,-n2:-n1), 
     $          Gi4(-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1), 
     $          Gi5(-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1) 
** 
      common /fillred/infilldim,infill(3) 
      dimension istorfill(3)
      dimension nphase(-1:1) 
      data nphase/-1,1,-1/ 
      parameter (pi   = 3.14159265358979324d0) 
*      print*,' eval1dhplatm1: y = ',y 
      if (infilldim.eq.2) then
         do i=1,2
            istorfill(i) = infill(i)
            infill(i) = -istorfill(i)
         enddo
      endif
** evaluate H(...,-y) 
      call setzero(nw,Gi1,Gi2,Gi3,Gi4,Gi5,-n2,-n1) 
      Gi1(0) = -1
      call eval1dhplat1(-y,nw,G1,G2,G3,G4,G5, 
     $                  GY1,GY2,GY3,GY4,GY5,Gi1,Gi2,Gi3,Gi4,Gi5,-n2,-n1) 
      if (infilldim.eq.2) then
         do i=1,2
            infill(i) = istorfill(i)
         enddo
      endif
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
                if ( nw.gt.4 ) then 
                  do k5=n1,n2 
                    nph5 = nph4*nphase(k5) 
                    HY5(k1,k2,k3,k4,k5) = nph5*GY5(-k1,-k2,-k3,-k4,-k5) 
                    Hi5(k1,k2,k3,k4,k5) =-nph5*Gi5(-k1,-k2,-k3,-k4,-k5) 
                    H5(k1,k2,k3,k4,k5)  = 
     $                dcmplx(HY5(k1,k2,k3,k4,k5),Hi5(k1,k2,k3,k4,k5)*pi) 
                  enddo 
                endif 
                  enddo 
                endif 
              enddo 
            endif 
          enddo 
        endif 
      enddo 
      return 
      end 


      subroutine eval1dhplatminf(y,nw,H1,H2,H3,H4,H5, 
     $                   HY1,HY2,HY3,HY4,HY5,Hi1,Hi2,Hi3,Hi4,Hi5,n1,n2) 
** evaluates 1dhpl's in the (-1)-range y  <= -(r2+1) 
** evaluating first the H(..,-y) by calling eval1dhplatinf(-y), 
** and then expressing H(..,y) in terms of H(..,-y) 
      implicit double precision (a-h,o-z) 
      complex*16 H1,H2,H3,H4,H5 
      complex*16 G1,G2,G3,G4,G5  
      dimension H1(n1:n2),H2(n1:n2,n1:n2),H3(n1:n2,n1:n2,n1:n2), 
     $          H4(n1:n2,n1:n2,n1:n2,n1:n2), 
     $          H5(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension HY1(n1:n2),HY2(n1:n2,n1:n2),HY3(n1:n2,n1:n2,n1:n2), 
     $          HY4(n1:n2,n1:n2,n1:n2,n1:n2), 
     $          HY5(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2), 
     $          Hi5(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2) 
** additional arrays required within this routine 
      dimension G1(-n2:-n1),G2(-n2:-n1,-n2:-n1),
     $          G3(-n2:-n1,-n2:-n1,-n2:-n1), 
     $          G4(-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1), 
     $          G5(-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1) 
      dimension GY1(-n2:-n1),GY2(-n2:-n1,-n2:-n1),
     $          GY3(-n2:-n1,-n2:-n1,-n2:-n1), 
     $          GY4(-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1), 
     $          GY5(-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1) 
      dimension Gi1(-n2:-n1),Gi2(-n2:-n1,-n2:-n1),
     $          Gi3(-n2:-n1,-n2:-n1,-n2:-n1), 
     $          Gi4(-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1), 
     $          Gi5(-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1) 
** 
      common /fillred/infilldim,infill(3) 
      dimension istorfill(3)
      dimension nphase(-1:1) 
      data nphase/-1,1,-1/ 
      parameter (pi   = 3.14159265358979324d0) 
*      print*,' eval1dhplatm1: y = ',y 
      if (infilldim.eq.2) then
         do i=1,2
            istorfill(i) = infill(i)
            infill(i) = -istorfill(i)
         enddo
      endif
** evaluate H(...,-y) 
      call setzero(nw,Gi1,Gi2,Gi3,Gi4,Gi5,-n2,-n1) 
      Gi1(0) = -1
      call eval1dhplatinf(-y,nw,G1,G2,G3,G4,G5, 
     $                  GY1,GY2,GY3,GY4,GY5,Gi1,Gi2,Gi3,Gi4,Gi5,-n2,-n1) 
      if (infilldim.eq.2) then
         do i=1,2
            infill(i) = istorfill(i)
         enddo
      endif
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
                if ( nw.gt.4 ) then 
                  do k5=n1,n2 
                    nph5 = nph4*nphase(k5) 
                    HY5(k1,k2,k3,k4,k5) = nph5*GY5(-k1,-k2,-k3,-k4,-k5) 
                    Hi5(k1,k2,k3,k4,k5) =-nph5*Gi5(-k1,-k2,-k3,-k4,-k5) 
                    H5(k1,k2,k3,k4,k5)  = 
     $                dcmplx(HY5(k1,k2,k3,k4,k5),Hi5(k1,k2,k3,k4,k5)*pi) 
                  enddo 
                endif 
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
      subroutine setzero(nw,Hi1,Hi2,Hi3,Hi4,Hi5,n1,n2) 
** initializes with 0 the elements of the arrays 
      implicit double precision (a-h,o-z) 
      dimension Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2), 
     $          Hi5(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2) 
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
                    if ( nw.gt.4 ) then 
                      do k5=n1,n2 
                        Hi5(k1,k2,k3,k4,k5) = 0.d0 
                      enddo 
                    endif 
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
      subroutine fillred1dhpl(nw,H1,H2,H3,H4,H5, 
     $                 HY1,HY2,HY3,HY4,HY5,Hi1,Hi2,Hi3,Hi4,Hi5,n1,n2) 
* fills the reducible 1dhpl from the irreducible set
      implicit double precision (a-h,o-z) 
      complex*16 H1,H2,H3,H4,H5 
      dimension H1(n1:n2),H2(n1:n2,n1:n2),H3(n1:n2,n1:n2,n1:n2), 
     $          H4(n1:n2,n1:n2,n1:n2,n1:n2), 
     $          H5(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension HY1(n1:n2),HY2(n1:n2,n1:n2),HY3(n1:n2,n1:n2,n1:n2), 
     $          HY4(n1:n2,n1:n2,n1:n2,n1:n2), 
     $          HY5(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2), 
     $          Hi5(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2) 
      common /fillred/infilldim,infill(3) 
      parameter (pinv = 0.318309886183790672d0) 
      parameter (pi   = 3.14159265358979324d0) 
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
            if ( nw.gt.4 ) then 
              do k5=n1,n2 
                H5(k1,k2,k3,k4,k5) = 
     $             dcmplx(HY5(k1,k2,k3,k4,k5),Hi5(k1,k2,k3,k4,k5)*pi) 
              enddo 
            endif 
              enddo 
            endif 
          enddo 
        endif 
      enddo 
      enddo 
** evaluating the reduced HPL's 
** iflag = 0 to suppress auxiliary printings of FILLREDHPLx 
      iflag = 0 
      do ia =  1,infilldim 
      do ib = ia,infilldim 
        call FILLREDHPL2(iflag,H1,H2,n1,n2,infill(ia),infill(ib)) 
        if ( nw.gt.2 ) then 
          do ic = ib,infilldim 
            call FILLREDHPL3(iflag,H1,H2,H3,n1,n2, 
     $                          infill(ia),infill(ib),infill(ic)) 
            if ( nw.gt.3 ) then 
              do id = ic,infilldim 
                call FILLREDHPL4(iflag,H1,H2,H3,H4,n1,n2, 
     $               infill(ia),infill(ib),infill(ic),infill(id)) 
              enddo 
            endif 
          enddo 
        endif 
      enddo 
      enddo 
      if (nw.gt.4) call FILLREDHPL5(iflag,H1,H2,H3,H4,H5,n1,n2)
** extractin real and immaginary parts from the complex value 
      do k1=n1,n2 
      do k2=n1,n2 
        HY2(k1,k2) =  dble(H2(k1,k2)) 
        Hi2(k1,k2) = dimag(H2(k1,k2))*pinv 
        if ( nw.gt.2 ) then 
          do k3=n1,n2 
            HY3(k1,k2,k3) =  dble(H3(k1,k2,k3)) 
            Hi3(k1,k2,k3) = dimag(H3(k1,k2,k3))*pinv 
            if ( nw.gt.3 ) then 
              do k4=n1,n2 
                HY4(k1,k2,k3,k4) =  dble(H4(k1,k2,k3,k4)) 
                Hi4(k1,k2,k3,k4) = dimag(H4(k1,k2,k3,k4))*pinv 
            if ( nw.gt.4 ) then 
              do k5=n1,n2 
                HY5(k1,k2,k3,k4,k5) =  dble(H5(k1,k2,k3,k4,k5)) 
                Hi5(k1,k2,k3,k4,k5) = dimag(H5(k1,k2,k3,k4,k5))*pinv 
              enddo 
            endif 
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


      subroutine FILLREDHPL5(iflag,HZ1,HZ2,HZ3,HZ4,HZ5,n1,n2)
      implicit double precision (a-h,o-z) 
      complex*16 HZ1,HZ2,HZ3,HZ4,HZ5 
      dimension HZ1(n1:n2),HZ2(n1:n2,n1:n2),HZ3(n1:n2,n1:n2,n1:n2) 
      dimension HZ4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension HZ5(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2) 

** evaluating the expansions 

      HZ5(0,0,0,0,0) = 
     $  + 8.3333333333333333d-03*HZ1(0)*HZ1(0)*HZ1(0)*HZ1(0)*HZ1(0)
	
** (n1,n2) = (0,1) or (-1,1) 
      if (    ( (n1.eq.0).and.(n2.eq.1) ) 
     $    .or.( (n1.eq.-1).and.(n2.eq.1) ) ) then 
      
      HZ5(0,0,0,1,0) = 
     $  + HZ1(0) *HZ4(0,0,0,1)
     $  - 4.0000000000000000d+00*HZ5(0,0,0,0,1)
      HZ5(0,0,1,0,0) = 
     $  + 5.0000000000000000d-01*HZ1(0)*HZ1(0)*HZ3(0,0,1)
     $  - 3.0000000000000000d+00*HZ1(0)*HZ4(0,0,0,1)
     $  + 6.0000000000000000d+00*HZ5(0,0,0,0,1)
      HZ5(0,0,1,1,0) = 
     $  + HZ1(0) *HZ4(0,0,1,1)
     $  - 3.0000000000000000d+00*HZ5(0,0,0,1,1)
     $  - HZ5(0,0,1,0,1) 
      HZ5(0,1,0,0,0) = 
     $  + 1.6666666666666666d-01*HZ1(0)*HZ1(0)*HZ1(0)*HZ2(0,1)
     $  - HZ1(0) *HZ1(0)*HZ3(0,0,1)
     $  + 3.0000000000000000d+00*HZ1(0)*HZ4(0,0,0,1)
     $  - 4.0000000000000000d+00*HZ5(0,0,0,0,1)
      HZ5(0,1,0,0,1) = 
     $  + HZ2(0,1) *HZ3(0,0,1)
     $  - 6.0000000000000000d+00*HZ5(0,0,0,1,1)
     $  - 3.0000000000000000d+00*HZ5(0,0,1,0,1)
      HZ5(0,1,0,1,0) = 
     $  + 5.0000000000000000d-01*HZ1(0)*HZ2(0,1)*HZ2(0,1)
     $  - 2.0000000000000000d+00*HZ1(0)*HZ4(0,0,1,1)
     $  - 2.0000000000000000d+00*HZ2(0,1)*HZ3(0,0,1)
     $  + 1.2000000000000000d+01*HZ5(0,0,0,1,1)
     $  + 4.0000000000000000d+00*HZ5(0,0,1,0,1)
      HZ5(0,1,1,0,0) = 
     $  + 5.0000000000000000d-01*HZ1(0)*HZ1(0)*HZ3(0,1,1)
     $  - 5.0000000000000000d-01*HZ1(0)*HZ2(0,1)*HZ2(0,1)
     $  + HZ2(0,1) *HZ3(0,0,1)
     $  - 3.0000000000000000d+00*HZ5(0,0,0,1,1)
     $  - HZ5(0,0,1,0,1) 
      HZ5(0,1,1,0,1) = 
     $  + HZ2(0,1) *HZ3(0,1,1)
     $  - 6.0000000000000000d+00*HZ5(0,0,1,1,1)
     $  - 3.0000000000000000d+00*HZ5(0,1,0,1,1)
      HZ5(0,1,1,1,0) = 
     $  + HZ1(0) *HZ4(0,1,1,1)
     $  - HZ2(0,1) *HZ3(0,1,1)
      HZ5(1,0,0,0,1) = 
     $  + HZ1(1) *HZ4(0,0,0,1)
     $  - HZ2(0,1) *HZ3(0,0,1)
     $  + 4.0000000000000000d+00*HZ5(0,0,0,1,1)
     $  + 2.0000000000000000d+00*HZ5(0,0,1,0,1)
      HZ5(1,0,0,1,0) = 
     $  + HZ1(0) *HZ1(1)*HZ3(0,0,1)
     $  - 5.0000000000000000d-01*HZ1(0)*HZ2(0,1)*HZ2(0,1)
     $  - 3.0000000000000000d+00*HZ1(1)*HZ4(0,0,0,1)
     $  + 2.0000000000000000d+00*HZ2(0,1)*HZ3(0,0,1)
     $  - 6.0000000000000000d+00*HZ5(0,0,0,1,1)
     $  - 3.0000000000000000d+00*HZ5(0,0,1,0,1)
      HZ5(1,0,0,1,1) = 
     $  + HZ1(1) *HZ4(0,0,1,1)
     $  - 3.0000000000000000d+00*HZ5(0,0,1,1,1)
     $  - HZ5(0,1,0,1,1) 
      HZ5(1,0,1,0,0) = 
     $  + 5.0000000000000000d-01*HZ1(0)*HZ1(0)*HZ1(1)*HZ2(0,1)
     $  - HZ1(0) *HZ1(0)*HZ3(0,1,1)
     $  - 2.0000000000000000d+00*HZ1(0)*HZ1(1)*HZ3(0,0,1)
     $  + 5.0000000000000000d-01*HZ1(0)*HZ2(0,1)*HZ2(0,1)
     $  + 2.0000000000000000d+00*HZ1(0)*HZ4(0,0,1,1)
     $  + 3.0000000000000000d+00*HZ1(1)*HZ4(0,0,0,1)
     $  - HZ2(0,1) *HZ3(0,0,1)
     $  + HZ5(0,0,1,0,1) 
      HZ5(1,0,1,0,1) = 
     $  + 5.0000000000000000d-01*HZ1(1)*HZ2(0,1)*HZ2(0,1)
     $  - 2.0000000000000000d+00*HZ1(1)*HZ4(0,0,1,1)
     $  - 2.0000000000000000d+00*HZ2(0,1)*HZ3(0,1,1)
     $  + 1.2000000000000000d+01*HZ5(0,0,1,1,1)
     $  + 4.0000000000000000d+00*HZ5(0,1,0,1,1)
      HZ5(1,0,1,1,0) = 
     $  + HZ1(0) *HZ1(1)*HZ3(0,1,1)
     $  - 3.0000000000000000d+00*HZ1(0)*HZ4(0,1,1,1)
     $  - 5.0000000000000000d-01*HZ1(1)*HZ2(0,1)*HZ2(0,1)
     $  + 2.0000000000000000d+00*HZ2(0,1)*HZ3(0,1,1)
     $  - 6.0000000000000000d+00*HZ5(0,0,1,1,1)
     $  - 3.0000000000000000d+00*HZ5(0,1,0,1,1)
      HZ5(1,0,1,1,1) = 
     $  + HZ1(1) *HZ4(0,1,1,1)
     $  - 4.0000000000000000d+00*HZ5(0,1,1,1,1)
      HZ5(1,1,0,0,0) = 
     $  + 8.3333333333333333d-02*HZ1(0)*HZ1(0)*HZ1(0)*HZ1(1)*HZ1(1)
     $  - 5.0000000000000000d-01*HZ1(0)*HZ1(0)*HZ1(1)*HZ2(0,1)
     $  + 5.0000000000000000d-01*HZ1(0)*HZ1(0)*HZ3(0,1,1)
     $  + HZ1(0) *HZ1(1)*HZ3(0,0,1)
     $  - HZ1(0) *HZ4(0,0,1,1)
     $  - HZ1(1) *HZ4(0,0,0,1)
     $  + HZ5(0,0,0,1,1) 
      HZ5(1,1,0,0,1) = 
     $  + 5.0000000000000000d-01*HZ1(1)*HZ1(1)*HZ3(0,0,1)
     $  - 5.0000000000000000d-01*HZ1(1)*HZ2(0,1)*HZ2(0,1)
     $  + HZ2(0,1) *HZ3(0,1,1)
     $  - 3.0000000000000000d+00*HZ5(0,0,1,1,1)
     $  - HZ5(0,1,0,1,1) 
      HZ5(1,1,0,1,0) = 
     $  + 5.0000000000000000d-01*HZ1(0)*HZ1(1)*HZ1(1)*HZ2(0,1)
     $  - 2.0000000000000000d+00*HZ1(0)*HZ1(1)*HZ3(0,1,1)
     $  + 3.0000000000000000d+00*HZ1(0)*HZ4(0,1,1,1)
     $  - HZ1(1) *HZ1(1)*HZ3(0,0,1)
     $  + 5.0000000000000000d-01*HZ1(1)*HZ2(0,1)*HZ2(0,1)
     $  + 2.0000000000000000d+00*HZ1(1)*HZ4(0,0,1,1)
     $  - HZ2(0,1) *HZ3(0,1,1)
     $  + HZ5(0,1,0,1,1) 
      HZ5(1,1,0,1,1) = 
     $  + 5.0000000000000000d-01*HZ1(1)*HZ1(1)*HZ3(0,1,1)
     $  - 3.0000000000000000d+00*HZ1(1)*HZ4(0,1,1,1)
     $  + 6.0000000000000000d+00*HZ5(0,1,1,1,1)
      HZ5(1,1,1,0,0) = 
     $  + 8.3333333333333333d-02*HZ1(0)*HZ1(0)*HZ1(1)*HZ1(1)*HZ1(1)
     $  - 5.0000000000000000d-01*HZ1(0)*HZ1(1)*HZ1(1)*HZ2(0,1)
     $  + HZ1(0) *HZ1(1)*HZ3(0,1,1)
     $  - HZ1(0) *HZ4(0,1,1,1)
     $  + 5.0000000000000000d-01*HZ1(1)*HZ1(1)*HZ3(0,0,1)
     $  - HZ1(1) *HZ4(0,0,1,1)
     $  + HZ5(0,0,1,1,1) 
      HZ5(1,1,1,0,1) = 
     $  + 1.6666666666666666d-01*HZ1(1)*HZ1(1)*HZ1(1)*HZ2(0,1)
     $  - HZ1(1) *HZ1(1)*HZ3(0,1,1)
     $  + 3.0000000000000000d+00*HZ1(1)*HZ4(0,1,1,1)
     $  - 4.0000000000000000d+00*HZ5(0,1,1,1,1)
      HZ5(1,1,1,1,0) = 
     $  + 4.1666666666666666d-02*HZ1(0)*HZ1(1)*HZ1(1)*HZ1(1)*HZ1(1)
     $  - 1.6666666666666666d-01*HZ1(1)*HZ1(1)*HZ1(1)*HZ2(0,1)
     $  + 5.0000000000000000d-01*HZ1(1)*HZ1(1)*HZ3(0,1,1)
     $  - HZ1(1) *HZ4(0,1,1,1)
     $  + HZ5(0,1,1,1,1) 
      HZ5(1,1,1,1,1) = 
     $  + 8.3333333333333333d-03*HZ1(1)*HZ1(1)*HZ1(1)*HZ1(1)*HZ1(1)
      endif

      if (    ( (n1.eq.-1).and.(n2.eq.0) ) 
     $    .or.( (n1.eq.-1).and.(n2.eq.1) ) ) then 
      HZ5(-1,-1,-1,-1,-1) = 
     $  + 8.3333333333333333d-03*HZ1(-1)*HZ1(-1)*HZ1(-1)*HZ1(-1)*HZ1(-1)
      HZ5(-1,-1,-1,-1,0) = 
     $  + 4.1666666666666666d-02*HZ1(-1)*HZ1(-1)*HZ1(-1)*HZ1(-1)*HZ1(0)
     $  - 1.6666666666666666d-01*HZ1(-1)*HZ1(-1)*HZ1(-1)*HZ2(0,-1)
     $  + 5.0000000000000000d-01*HZ1(-1)*HZ1(-1)*HZ3(0,-1,-1)
     $  - HZ1( -1)*HZ4(0,-1,-1,-1)
     $  + HZ5(0, -1,-1,-1,-1)
      HZ5(-1,-1,-1,0,-1) = 
     $  + 1.6666666666666666d-01*HZ1(-1)*HZ1(-1)*HZ1(-1)*HZ2(0,-1)
     $  - HZ1( -1)*HZ1(-1)*HZ3(0,-1,-1)
     $  + 3.0000000000000000d+00*HZ1(-1)*HZ4(0,-1,-1,-1)
     $  - 4.0000000000000000d+00*HZ5(0,-1,-1,-1,-1)
      HZ5(-1,-1,-1,0,0) = 
     $  + 8.3333333333333333d-02*HZ1(-1)*HZ1(-1)*HZ1(-1)*HZ1(0)*HZ1(0)
     $  - 5.0000000000000000d-01*HZ1(-1)*HZ1(-1)*HZ1(0)*HZ2(0,-1)
     $  + 5.0000000000000000d-01*HZ1(-1)*HZ1(-1)*HZ3(0,0,-1)
     $  + HZ1( -1)*HZ1(0)*HZ3(0,-1,-1)
     $  - HZ1( -1)*HZ4(0,0,-1,-1)
     $  - HZ1(0) *HZ4(0,-1,-1,-1)
     $  + HZ5(0,0, -1,-1,-1)
      HZ5(-1,-1,0,-1,-1) = 
     $  + 5.0000000000000000d-01*HZ1(-1)*HZ1(-1)*HZ3(0,-1,-1)
     $  - 3.0000000000000000d+00*HZ1(-1)*HZ4(0,-1,-1,-1)
     $  + 6.0000000000000000d+00*HZ5(0,-1,-1,-1,-1)
      HZ5(-1,-1,0,-1,0) = 
     $  + 5.0000000000000000d-01*HZ1(-1)*HZ1(-1)*HZ1(0)*HZ2(0,-1)
     $  - HZ1( -1)*HZ1(-1)*HZ3(0,0,-1)
     $  - 2.0000000000000000d+00*HZ1(-1)*HZ1(0)*HZ3(0,-1,-1)
     $  + 5.0000000000000000d-01*HZ1(-1)*HZ2(0,-1)*HZ2(0,-1)
     $  + 2.0000000000000000d+00*HZ1(-1)*HZ4(0,0,-1,-1)
     $  + 3.0000000000000000d+00*HZ1(0)*HZ4(0,-1,-1,-1)
     $  - HZ2(0, -1)*HZ3(0,-1,-1)
     $  + HZ5(0, -1,0,-1,-1)
      HZ5(-1,-1,0,0,-1) = 
     $  + 5.0000000000000000d-01*HZ1(-1)*HZ1(-1)*HZ3(0,0,-1)
     $  - 5.0000000000000000d-01*HZ1(-1)*HZ2(0,-1)*HZ2(0,-1)
     $  + HZ2(0, -1)*HZ3(0,-1,-1)
     $  - HZ5(0, -1,0,-1,-1)
     $  - 3.0000000000000000d+00*HZ5(0,0,-1,-1,-1)
      HZ5(-1,-1,0,0,0) = 
     $  + 8.3333333333333333d-02*HZ1(-1)*HZ1(-1)*HZ1(0)*HZ1(0)*HZ1(0)
     $  - 5.0000000000000000d-01*HZ1(-1)*HZ1(0)*HZ1(0)*HZ2(0,-1)
     $  + HZ1( -1)*HZ1(0)*HZ3(0,0,-1)
     $  - HZ1( -1)*HZ4(0,0,0,-1)
     $  + 5.0000000000000000d-01*HZ1(0)*HZ1(0)*HZ3(0,-1,-1)
     $  - HZ1(0) *HZ4(0,0,-1,-1)
     $  + HZ5(0,0,0, -1,-1)
      HZ5(-1,0,-1,-1,-1) = 
     $  + HZ1( -1)*HZ4(0,-1,-1,-1)
     $  - 4.0000000000000000d+00*HZ5(0,-1,-1,-1,-1)
      HZ5(-1,0,-1,-1,0) = 
     $  + HZ1( -1)*HZ1(0)*HZ3(0,-1,-1)
     $  - 5.0000000000000000d-01*HZ1(-1)*HZ2(0,-1)*HZ2(0,-1)
     $  - 3.0000000000000000d+00*HZ1(0)*HZ4(0,-1,-1,-1)
     $  + 2.0000000000000000d+00*HZ2(0,-1)*HZ3(0,-1,-1)
     $  - 3.0000000000000000d+00*HZ5(0,-1,0,-1,-1)
     $  - 6.0000000000000000d+00*HZ5(0,0,-1,-1,-1)
      HZ5(-1,0,-1,0,-1) = 
     $  + 5.0000000000000000d-01*HZ1(-1)*HZ2(0,-1)*HZ2(0,-1)
     $  - 2.0000000000000000d+00*HZ1(-1)*HZ4(0,0,-1,-1)
     $  - 2.0000000000000000d+00*HZ2(0,-1)*HZ3(0,-1,-1)
     $  + 4.0000000000000000d+00*HZ5(0,-1,0,-1,-1)
     $  + 1.2000000000000000d+01*HZ5(0,0,-1,-1,-1)
      HZ5(-1,0,-1,0,0) = 
     $  + 5.0000000000000000d-01*HZ1(-1)*HZ1(0)*HZ1(0)*HZ2(0,-1)
     $  - 2.0000000000000000d+00*HZ1(-1)*HZ1(0)*HZ3(0,0,-1)
     $  + 3.0000000000000000d+00*HZ1(-1)*HZ4(0,0,0,-1)
     $  - HZ1(0) *HZ1(0)*HZ3(0,-1,-1)
     $  + 5.0000000000000000d-01*HZ1(0)*HZ2(0,-1)*HZ2(0,-1)
     $  + 2.0000000000000000d+00*HZ1(0)*HZ4(0,0,-1,-1)
     $  - HZ2(0, -1)*HZ3(0,0,-1)
     $  + HZ5(0,0, -1,0,-1)
      HZ5(-1,0,0,-1,-1) = 
     $  + HZ1( -1)*HZ4(0,0,-1,-1)
     $  - HZ5(0, -1,0,-1,-1)
     $  - 3.0000000000000000d+00*HZ5(0,0,-1,-1,-1)
      HZ5(-1,0,0,-1,0) = 
     $  + HZ1( -1)*HZ1(0)*HZ3(0,0,-1)
     $  - 3.0000000000000000d+00*HZ1(-1)*HZ4(0,0,0,-1)
     $  - 5.0000000000000000d-01*HZ1(0)*HZ2(0,-1)*HZ2(0,-1)
     $  + 2.0000000000000000d+00*HZ2(0,-1)*HZ3(0,0,-1)
     $  - 3.0000000000000000d+00*HZ5(0,0,-1,0,-1)
     $  - 6.0000000000000000d+00*HZ5(0,0,0,-1,-1)
      HZ5(-1,0,0,0,-1) = 
     $  + HZ1( -1)*HZ4(0,0,0,-1)
     $  - HZ2(0, -1)*HZ3(0,0,-1)
     $  + 2.0000000000000000d+00*HZ5(0,0,-1,0,-1)
     $  + 4.0000000000000000d+00*HZ5(0,0,0,-1,-1)
      HZ5(-1,0,0,0,0) = 
     $  + 4.1666666666666666d-02*HZ1(-1)*HZ1(0)*HZ1(0)*HZ1(0)*HZ1(0)
     $  - 1.6666666666666666d-01*HZ1(0)*HZ1(0)*HZ1(0)*HZ2(0,-1)
     $  + 5.0000000000000000d-01*HZ1(0)*HZ1(0)*HZ3(0,0,-1)
     $  - HZ1(0) *HZ4(0,0,0,-1)
     $  + HZ5(0,0,0,0, -1)
      HZ5(0,-1,-1,-1,0) = 
     $  + HZ1(0) *HZ4(0,-1,-1,-1)
     $  - HZ2(0, -1)*HZ3(0,-1,-1)
     $  + 2.0000000000000000d+00*HZ5(0,-1,0,-1,-1)
     $  + 4.0000000000000000d+00*HZ5(0,0,-1,-1,-1)
      HZ5(0,-1,-1,0,-1) = 
     $  + HZ2(0, -1)*HZ3(0,-1,-1)
     $  - 3.0000000000000000d+00*HZ5(0,-1,0,-1,-1)
     $  - 6.0000000000000000d+00*HZ5(0,0,-1,-1,-1)
      HZ5(0,-1,-1,0,0) = 
     $  + 5.0000000000000000d-01*HZ1(0)*HZ1(0)*HZ3(0,-1,-1)
     $  - 5.0000000000000000d-01*HZ1(0)*HZ2(0,-1)*HZ2(0,-1)
     $  + HZ2(0, -1)*HZ3(0,0,-1)
     $  - HZ5(0,0, -1,0,-1)
     $  - 3.0000000000000000d+00*HZ5(0,0,0,-1,-1)
      HZ5(0,-1,0,-1,0) = 
     $  + 5.0000000000000000d-01*HZ1(0)*HZ2(0,-1)*HZ2(0,-1)
     $  - 2.0000000000000000d+00*HZ1(0)*HZ4(0,0,-1,-1)
     $  - 2.0000000000000000d+00*HZ2(0,-1)*HZ3(0,0,-1)
     $  + 4.0000000000000000d+00*HZ5(0,0,-1,0,-1)
     $  + 1.2000000000000000d+01*HZ5(0,0,0,-1,-1)
      HZ5(0,-1,0,0,-1) = 
     $  + HZ2(0, -1)*HZ3(0,0,-1)
     $  - 3.0000000000000000d+00*HZ5(0,0,-1,0,-1)
     $  - 6.0000000000000000d+00*HZ5(0,0,0,-1,-1)
      HZ5(0,-1,0,0,0) = 
     $  + 1.6666666666666666d-01*HZ1(0)*HZ1(0)*HZ1(0)*HZ2(0,-1)
     $  - HZ1(0) *HZ1(0)*HZ3(0,0,-1)
     $  + 3.0000000000000000d+00*HZ1(0)*HZ4(0,0,0,-1)
     $  - 4.0000000000000000d+00*HZ5(0,0,0,0,-1)
      HZ5(0,0,-1,-1,0) = 
     $  + HZ1(0) *HZ4(0,0,-1,-1)
     $  - HZ5(0,0, -1,0,-1)
     $  - 3.0000000000000000d+00*HZ5(0,0,0,-1,-1)
      HZ5(0,0,-1,0,0) = 
     $  + 5.0000000000000000d-01*HZ1(0)*HZ1(0)*HZ3(0,0,-1)
     $  - 3.0000000000000000d+00*HZ1(0)*HZ4(0,0,0,-1)
     $  + 6.0000000000000000d+00*HZ5(0,0,0,0,-1)
      HZ5(0,0,0,-1,0) = 
     $  + HZ1(0) *HZ4(0,0,0,-1)
     $  - 4.0000000000000000d+00*HZ5(0,0,0,0,-1)
      endif

      if ( (n1.eq.-1).and.(n2.eq.1) ) then 
      HZ5(-1,-1,-1,0,1) = 
     $  + 1.6666666666666666d-01*HZ1(-1)*HZ1(-1)*HZ1(-1)*HZ2(0,1)
     $  - 5.0000000000000000d-01*HZ1(-1)*HZ1(-1)*HZ3(-1,0,1)
     $  + HZ1( -1)*HZ4(-1,-1,0,1)
     $  - HZ5(0, -1,-1,-1,1)
     $  - HZ5(0, -1,-1,1,-1)
     $  - HZ5(0, -1,1,-1,-1)
     $  - HZ5(0,1, -1,-1,-1)
      HZ5(-1,-1,-1,1,-1) = 
     $  + HZ1( -1)*HZ4(-1,-1,-1,1)
     $  - 4.0000000000000000d+00*HZ5(-1,-1,-1,-1,1)
      HZ5(-1,-1,-1,1,0) = 
     $  - 1.6666666666666666d-01*HZ1(-1)*HZ1(-1)*HZ1(-1)*HZ2(0,1)
     $  + 5.0000000000000000d-01*HZ1(-1)*HZ1(-1)*HZ3(-1,0,1)
     $  + 5.0000000000000000d-01*HZ1(-1)*HZ1(-1)*HZ3(0,-1,1)
     $  - HZ1( -1)*HZ4(-1,-1,0,1)
     $  - HZ1( -1)*HZ4(-1,0,-1,1)
     $  - HZ1( -1)*HZ4(0,-1,-1,1)
     $  + HZ1(0) *HZ4(-1,-1,-1,1)
     $  + HZ5(0,1, -1,-1,-1)
      HZ5(-1,-1,0,-1,1) = 
     $  - 5.0000000000000000d-01*HZ1(-1)*HZ1(-1)*HZ3(0,-1,1)
     $  + HZ1( -1)*HZ4(-1,0,-1,1)
     $  + 3.0000000000000000d+00*HZ5(0,-1,-1,-1,1)
     $  + 2.0000000000000000d+00*HZ5(0,-1,-1,1,-1)
     $  + HZ5(0, -1,1,-1,-1)
      HZ5(-1,-1,0,0,1) = 
     $  - 5.0000000000000000d-01*HZ1(-1)*HZ1(-1)*HZ3(0,0,1)
     $  + HZ1( -1)*HZ4(-1,0,0,1)
     $  + HZ5(0, -1,-1,0,1)
     $  + HZ5(0, -1,0,-1,1)
     $  + HZ5(0, -1,0,1,-1)
     $  + HZ5(0,0, -1,-1,1)
     $  + HZ5(0,0, -1,1,-1)
     $  + HZ5(0,0,1, -1,-1)
      HZ5(-1,-1,0,1,-1) = 
     $  - 5.0000000000000000d-01*HZ1(-1)*HZ1(-1)*HZ1(-1)*HZ2(0,1)
     $  + 1.5000000000000000d+00*HZ1(-1)*HZ1(-1)*HZ3(-1,0,1)
     $  + 5.0000000000000000d-01*HZ1(-1)*HZ1(-1)*HZ3(0,-1,1)
     $  - 2.0000000000000000d+00*HZ1(-1)*HZ4(-1,-1,0,1)
     $  - HZ1( -1)*HZ4(-1,0,-1,1)
     $  + HZ5(0, -1,-1,1,-1)
     $  + 2.0000000000000000d+00*HZ5(0,-1,1,-1,-1)
     $  + 3.0000000000000000d+00*HZ5(0,1,-1,-1,-1)
      HZ5(-1,-1,0,1,0) = 
     $  + HZ1( -1)*HZ1(-1)*HZ3(0,0,1)
     $  - 2.0000000000000000d+00*HZ1(-1)*HZ4(-1,0,0,1)
     $  - HZ1( -1)*HZ4(0,-1,0,1)
     $  + HZ1(0) *HZ4(-1,-1,0,1)
     $  - HZ5(0, -1,-1,0,1)
     $  - HZ5(0, -1,0,-1,1)
     $  - HZ5(0, -1,0,1,-1)
     $  - 2.0000000000000000d+00*HZ5(0,0,-1,-1,1)
     $  - 2.0000000000000000d+00*HZ5(0,0,-1,1,-1)
     $  - 2.0000000000000000d+00*HZ5(0,0,1,-1,-1)
      HZ5(-1,-1,0,1,1) = 
     $  - 5.0000000000000000d-01*HZ1(-1)*HZ1(-1)*HZ3(0,1,1)
     $  + HZ1( -1)*HZ4(-1,0,1,1)
     $  + HZ5(0, -1,-1,1,1)
     $  + HZ5(0, -1,1,-1,1)
     $  + HZ5(0, -1,1,1,-1)
     $  + HZ5(0,1, -1,-1,1)
     $  + HZ5(0,1, -1,1,-1)
     $  + HZ5(0,1,1, -1,-1)
      HZ5(-1,-1,1,-1,-1) = 
     $  + 5.0000000000000000d-01*HZ1(-1)*HZ1(-1)*HZ3(-1,-1,1)
     $  - 3.0000000000000000d+00*HZ1(-1)*HZ4(-1,-1,-1,1)
     $  + 6.0000000000000000d+00*HZ5(-1,-1,-1,-1,1)
      HZ5(-1,-1,1,-1,0) = 
     $  - 5.0000000000000000d-01*HZ1(-1)*HZ1(-1)*HZ3(0,-1,1)
     $  + HZ1( -1)*HZ1(0)*HZ3(-1,-1,1)
     $  + HZ1( -1)*HZ4(-1,0,-1,1)
     $  + 2.0000000000000000d+00*HZ1(-1)*HZ4(0,-1,-1,1)
     $  - 3.0000000000000000d+00*HZ1(0)*HZ4(-1,-1,-1,1)
     $  - HZ2(0, -1)*HZ3(-1,-1,1)
     $  + HZ5(0, -1,1,-1,-1)
      HZ5(-1,-1,1,0,-1) = 
     $  + 5.0000000000000000d-01*HZ1(-1)*HZ1(-1)*HZ1(-1)*HZ2(0,1)
     $  - 1.5000000000000000d+00*HZ1(-1)*HZ1(-1)*HZ3(-1,0,1)
     $  - HZ1( -1)*HZ1(-1)*HZ3(0,-1,1)
     $  + 2.0000000000000000d+00*HZ1(-1)*HZ4(-1,-1,0,1)
     $  + HZ1( -1)*HZ4(-1,0,-1,1)
     $  + HZ2(0, -1)*HZ3(-1,-1,1)
     $  - HZ5(0, -1,1,-1,-1)
     $  - 3.0000000000000000d+00*HZ5(0,1,-1,-1,-1)
      HZ5(-1,-1,1,0,0) = 
     $  - 5.0000000000000000d-01*HZ1(-1)*HZ1(-1)*HZ3(0,0,1)
     $  + HZ1( -1)*HZ4(-1,0,0,1)
     $  + HZ1( -1)*HZ4(0,-1,0,1)
     $  + HZ1( -1)*HZ4(0,0,-1,1)
     $  + 5.0000000000000000d-01*HZ1(0)*HZ1(0)*HZ3(-1,-1,1)
     $  - HZ1(0) *HZ4(-1,-1,0,1)
     $  - HZ1(0) *HZ4(-1,0,-1,1)
     $  - HZ1(0) *HZ4(0,-1,-1,1)
     $  + HZ5(0,0,1, -1,-1)
      HZ5(-1,-1,1,0,1) = 
     $  + HZ1( -1)*HZ1(-1)*HZ3(0,1,1)
     $  - 2.0000000000000000d+00*HZ1(-1)*HZ4(-1,0,1,1)
     $  - 2.0000000000000000d+00*HZ1(-1)*HZ4(0,-1,1,1)
     $  - HZ1( -1)*HZ4(0,1,-1,1)
     $  + HZ2(0,1) *HZ3(-1,-1,1)
     $  - HZ5(0,1, -1,-1,1)
     $  - HZ5(0,1, -1,1,-1)
     $  - 2.0000000000000000d+00*HZ5(0,1,1,-1,-1)
      HZ5(-1,-1,1,1,-1) = 
     $  + HZ1( -1)*HZ4(-1,-1,1,1)
     $  - 3.0000000000000000d+00*HZ5(-1,-1,-1,1,1)
     $  - HZ5( -1,-1,1,-1,1)
      HZ5(-1,-1,1,1,0) = 
     $  - 5.0000000000000000d-01*HZ1(-1)*HZ1(-1)*HZ3(0,1,1)
     $  + HZ1( -1)*HZ4(-1,0,1,1)
     $  + HZ1( -1)*HZ4(0,-1,1,1)
     $  + HZ1( -1)*HZ4(0,1,-1,1)
     $  + HZ1(0) *HZ4(-1,-1,1,1)
     $  - HZ2(0,1) *HZ3(-1,-1,1)
     $  + HZ5(0,1,1, -1,-1)
      HZ5(-1,0,-1,-1,1) = 
     $  + HZ1( -1)*HZ4(0,-1,-1,1)
     $  - 3.0000000000000000d+00*HZ5(0,-1,-1,-1,1)
     $  - HZ5(0, -1,-1,1,-1)
      HZ5(-1,0,-1,0,1) = 
     $  + HZ1( -1)*HZ4(0,-1,0,1)
     $  - 2.0000000000000000d+00*HZ5(0,-1,-1,0,1)
     $  - HZ5(0, -1,0,-1,1)
     $  - HZ5(0, -1,0,1,-1)
      HZ5(-1,0,-1,1,-1) = 
     $  + HZ1( -1)*HZ1(-1)*HZ3(0,-1,1)
     $  - HZ1( -1)*HZ4(-1,0,-1,1)
     $  - 2.0000000000000000d+00*HZ1(-1)*HZ4(0,-1,-1,1)
     $  - 2.0000000000000000d+00*HZ5(0,-1,-1,1,-1)
     $  - 2.0000000000000000d+00*HZ5(0,-1,1,-1,-1)
      HZ5(-1,0,-1,1,0) = 
     $  - HZ1( -1)*HZ4(0,-1,0,1)
     $  - 2.0000000000000000d+00*HZ1(-1)*HZ4(0,0,-1,1)
     $  + HZ1(0) *HZ4(-1,0,-1,1)
     $  + 2.0000000000000000d+00*HZ5(0,-1,-1,0,1)
     $  + 2.0000000000000000d+00*HZ5(0,-1,0,-1,1)
     $  + HZ5(0, -1,0,1,-1)
     $  + 4.0000000000000000d+00*HZ5(0,0,-1,-1,1)
     $  + 2.0000000000000000d+00*HZ5(0,0,-1,1,-1)
      HZ5(-1,0,-1,1,1) = 
     $  + HZ1( -1)*HZ4(0,-1,1,1)
     $  - 2.0000000000000000d+00*HZ5(0,-1,-1,1,1)
     $  - HZ5(0, -1,1,-1,1)
     $  - HZ5(0, -1,1,1,-1)
      HZ5(-1,0,0,-1,1) = 
     $  + HZ1( -1)*HZ4(0,0,-1,1)
     $  - HZ5(0, -1,0,-1,1)
     $  - 2.0000000000000000d+00*HZ5(0,0,-1,-1,1)
     $  - HZ5(0,0, -1,1,-1)
      HZ5(-1,0,0,0,1) = 
     $  + HZ1( -1)*HZ4(0,0,0,1)
     $  - HZ2(0, -1)*HZ3(0,0,1)
     $  + HZ5(0,0, -1,0,1)
     $  + 2.0000000000000000d+00*HZ5(0,0,0,-1,1)
     $  + 2.0000000000000000d+00*HZ5(0,0,0,1,-1)
     $  + HZ5(0,0,1,0, -1)
      HZ5(-1,0,0,1,-1) = 
     $  + HZ1( -1)*HZ1(-1)*HZ3(0,0,1)
     $  - HZ1( -1)*HZ4(-1,0,0,1)
     $  - HZ1( -1)*HZ4(0,-1,0,1)
     $  - HZ1( -1)*HZ4(0,0,-1,1)
     $  - HZ5(0, -1,0,1,-1)
     $  - HZ5(0,0, -1,1,-1)
     $  - 2.0000000000000000d+00*HZ5(0,0,1,-1,-1)
      HZ5(-1,0,0,1,0) = 
     $  - 3.0000000000000000d+00*HZ1(-1)*HZ4(0,0,0,1)
     $  + HZ1(0) *HZ4(-1,0,0,1)
     $  + 2.0000000000000000d+00*HZ2(0,-1)*HZ3(0,0,1)
     $  - HZ5(0,0, -1,0,1)
     $  - 3.0000000000000000d+00*HZ5(0,0,0,-1,1)
     $  - 3.0000000000000000d+00*HZ5(0,0,0,1,-1)
     $  - 2.0000000000000000d+00*HZ5(0,0,1,0,-1)
      HZ5(-1,0,0,1,1) = 
     $  + HZ1( -1)*HZ4(0,0,1,1)
     $  - HZ5(0, -1,0,1,1)
     $  - HZ5(0,0, -1,1,1)
     $  - HZ5(0,0,1, -1,1)
     $  - HZ5(0,0,1,1, -1)
      HZ5(-1,0,1,-1,-1) = 
     $  + 5.0000000000000000d-01*HZ1(-1)*HZ1(-1)*HZ1(-1)*HZ2(0,1)
     $  - HZ1( -1)*HZ1(-1)*HZ3(-1,0,1)
     $  - HZ1( -1)*HZ1(-1)*HZ3(0,-1,1)
     $  + HZ1( -1)*HZ4(-1,-1,0,1)
     $  + HZ1( -1)*HZ4(-1,0,-1,1)
     $  + HZ1( -1)*HZ4(0,-1,-1,1)
     $  - HZ5(0, -1,1,-1,-1)
     $  - 3.0000000000000000d+00*HZ5(0,1,-1,-1,-1)
      HZ5(-1,0,1,-1,0) = 
     $  + HZ1( -1)*HZ1(0)*HZ3(-1,0,1)
     $  + HZ1( -1)*HZ4(0,-1,0,1)
     $  + 2.0000000000000000d+00*HZ1(-1)*HZ4(0,0,-1,1)
     $  - 2.0000000000000000d+00*HZ1(0)*HZ4(-1,-1,0,1)
     $  - HZ1(0) *HZ4(-1,0,-1,1)
     $  - HZ2(0, -1)*HZ3(-1,0,1)
     $  - 2.0000000000000000d+00*HZ5(0,-1,0,-1,1)
     $  - HZ5(0, -1,0,1,-1)
     $  - 4.0000000000000000d+00*HZ5(0,0,-1,-1,1)
     $  - 2.0000000000000000d+00*HZ5(0,0,-1,1,-1)
      HZ5(-1,0,1,-1,1) = 
     $  + HZ1( -1)*HZ4(0,1,-1,1)
     $  - HZ5(0, -1,1,-1,1)
     $  - 2.0000000000000000d+00*HZ5(0,1,-1,-1,1)
     $  - HZ5(0,1, -1,1,-1)
      HZ5(-1,0,1,0,-1) = 
     $  - 2.0000000000000000d+00*HZ1(-1)*HZ1(-1)*HZ3(0,0,1)
     $  + 2.0000000000000000d+00*HZ1(-1)*HZ4(-1,0,0,1)
     $  + HZ1( -1)*HZ4(0,-1,0,1)
     $  + HZ2(0, -1)*HZ3(-1,0,1)
     $  + 2.0000000000000000d+00*HZ5(0,-1,0,-1,1)
     $  + 2.0000000000000000d+00*HZ5(0,-1,0,1,-1)
     $  + 4.0000000000000000d+00*HZ5(0,0,-1,-1,1)
     $  + 4.0000000000000000d+00*HZ5(0,0,-1,1,-1)
     $  + 4.0000000000000000d+00*HZ5(0,0,1,-1,-1)
      HZ5(-1,0,1,0,0) = 
     $  + 3.0000000000000000d+00*HZ1(-1)*HZ4(0,0,0,1)
     $  + 5.0000000000000000d-01*HZ1(0)*HZ1(0)*HZ3(-1,0,1)
     $  - 2.0000000000000000d+00*HZ1(0)*HZ4(-1,0,0,1)
     $  - HZ1(0) *HZ4(0,-1,0,1)
     $  - HZ2(0, -1)*HZ3(0,0,1)
     $  + HZ5(0,0,1,0, -1)
      HZ5(-1,0,1,0,1) = 
     $  - 5.0000000000000000d-01*HZ1(-1)*HZ2(0,1)*HZ2(0,1)
     $  - 2.0000000000000000d+00*HZ1(-1)*HZ4(0,0,1,1)
     $  + HZ2(0,1) *HZ3(-1,0,1)
     $  + HZ2(0,1) *HZ3(0,-1,1)
     $  - HZ5(0, -1,1,0,1)
     $  + 2.0000000000000000d+00*HZ5(0,0,1,-1,1)
     $  + 4.0000000000000000d+00*HZ5(0,0,1,1,-1)
     $  + HZ5(0,1,0,1, -1)
      HZ5(-1,0,1,1,-1) = 
     $  + HZ1( -1)*HZ1(-1)*HZ3(0,1,1)
     $  - HZ1( -1)*HZ4(-1,0,1,1)
     $  - HZ1( -1)*HZ4(0,-1,1,1)
     $  - HZ1( -1)*HZ4(0,1,-1,1)
     $  - HZ5(0, -1,1,1,-1)
     $  - HZ5(0,1, -1,1,-1)
     $  - 2.0000000000000000d+00*HZ5(0,1,1,-1,-1)
      HZ5(-1,0,1,1,0) = 
     $  + 5.0000000000000000d-01*HZ1(-1)*HZ2(0,1)*HZ2(0,1)
     $  + HZ1(0) *HZ4(-1,0,1,1)
     $  - HZ2(0,1) *HZ3(-1,0,1)
     $  - HZ2(0,1) *HZ3(0,-1,1)
     $  + HZ5(0, -1,0,1,1)
     $  + HZ5(0, -1,1,0,1)
     $  + 2.0000000000000000d+00*HZ5(0,0,-1,1,1)
     $  - 2.0000000000000000d+00*HZ5(0,0,1,1,-1)
     $  - HZ5(0,1,0,1, -1)
      HZ5(-1,0,1,1,1) = 
     $  + HZ1( -1)*HZ4(0,1,1,1)
     $  - HZ5(0, -1,1,1,1)
     $  - HZ5(0,1, -1,1,1)
     $  - HZ5(0,1,1, -1,1)
     $  - HZ5(0,1,1,1, -1)
      HZ5(-1,1,-1,-1,-1) = 
     $  + 1.6666666666666666d-01*HZ1(-1)*HZ1(-1)*HZ1(-1)*HZ2(-1,1)
     $  - HZ1( -1)*HZ1(-1)*HZ3(-1,-1,1)
     $  + 3.0000000000000000d+00*HZ1(-1)*HZ4(-1,-1,-1,1)
     $  - 4.0000000000000000d+00*HZ5(-1,-1,-1,-1,1)
      HZ5(-1,1,-1,-1,0) = 
     $  + 5.0000000000000000d-01*HZ1(-1)*HZ1(-1)*HZ1(0)*HZ2(-1,1)
     $  - 2.0000000000000000d+00*HZ1(-1)*HZ1(0)*HZ3(-1,-1,1)
     $  - HZ1( -1)*HZ2(-1,1)*HZ2(0,-1)
     $  - HZ1( -1)*HZ4(0,-1,-1,1)
     $  + 3.0000000000000000d+00*HZ1(0)*HZ4(-1,-1,-1,1)
     $  + HZ2( -1,1)*HZ3(0,-1,-1)
     $  + 2.0000000000000000d+00*HZ2(0,-1)*HZ3(-1,-1,1)
     $  + HZ5(0, -1,-1,1,-1)
      HZ5(-1,1,-1,-1,1) = 
     $  + HZ2( -1,1)*HZ3(-1,-1,1)
     $  - 6.0000000000000000d+00*HZ5(-1,-1,-1,1,1)
     $  - 3.0000000000000000d+00*HZ5(-1,-1,1,-1,1)
      HZ5(-1,1,-1,0,-1) = 
     $  + HZ1( -1)*HZ1(-1)*HZ3(0,-1,1)
     $  + HZ1( -1)*HZ2(-1,1)*HZ2(0,-1)
     $  - HZ1( -1)*HZ4(-1,0,-1,1)
     $  - 2.0000000000000000d+00*HZ2(-1,1)*HZ3(0,-1,-1)
     $  - 2.0000000000000000d+00*HZ2(0,-1)*HZ3(-1,-1,1)
     $  - 2.0000000000000000d+00*HZ5(0,-1,-1,1,-1)
     $  - 2.0000000000000000d+00*HZ5(0,-1,1,-1,-1)
      HZ5(-1,1,-1,0,0) = 
     $  + 5.0000000000000000d-01*HZ1(-1)*HZ1(0)*HZ1(0)*HZ2(-1,1)
     $  - HZ1( -1)*HZ4(0,0,-1,1)
     $  - HZ1(0) *HZ1(0)*HZ3(-1,-1,1)
     $  - HZ1(0) *HZ2(-1,1)*HZ2(0,-1)
     $  + HZ1(0) *HZ4(-1,0,-1,1)
     $  + 2.0000000000000000d+00*HZ1(0)*HZ4(0,-1,-1,1)
     $  + HZ2( -1,1)*HZ3(0,0,-1)
     $  + HZ5(0,0, -1,1,-1)
      HZ5(-1,1,-1,0,1) = 
     $  + 2.0000000000000000d+00*HZ1(-1)*HZ4(0,-1,1,1)
     $  + HZ1( -1)*HZ4(0,1,-1,1)
     $  + HZ2( -1,1)*HZ3(-1,0,1)
     $  - 2.0000000000000000d+00*HZ2(0,1)*HZ3(-1,-1,1)
     $  - HZ5(0, -1,1,-1,1)
     $  - 2.0000000000000000d+00*HZ5(0,-1,1,1,-1)
     $  - HZ5(0,1, -1,1,-1)
      HZ5(-1,1,-1,1,-1) = 
     $  + 5.0000000000000000d-01*HZ1(-1)*HZ2(-1,1)*HZ2(-1,1)
     $  - 2.0000000000000000d+00*HZ1(-1)*HZ4(-1,-1,1,1)
     $  - 2.0000000000000000d+00*HZ2(-1,1)*HZ3(-1,-1,1)
     $  + 1.2000000000000000d+01*HZ5(-1,-1,-1,1,1)
     $  + 4.0000000000000000d+00*HZ5(-1,-1,1,-1,1)
      HZ5(-1,1,-1,1,0) = 
     $  - HZ1( -1)*HZ4(0,1,-1,1)
     $  + 5.0000000000000000d-01*HZ1(0)*HZ2(-1,1)*HZ2(-1,1)
     $  - 2.0000000000000000d+00*HZ1(0)*HZ4(-1,-1,1,1)
     $  - HZ2( -1,1)*HZ3(-1,0,1)
     $  - HZ2( -1,1)*HZ3(0,-1,1)
     $  + 2.0000000000000000d+00*HZ2(0,1)*HZ3(-1,-1,1)
     $  + HZ5(0,1, -1,1,-1)
      HZ5(-1,1,0,-1,-1) = 
     $  - 5.0000000000000000d-01*HZ1(-1)*HZ1(-1)*HZ1(-1)*HZ2(0,1)
     $  + HZ1( -1)*HZ1(-1)*HZ3(-1,0,1)
     $  - HZ1( -1)*HZ4(-1,-1,0,1)
     $  + HZ2( -1,1)*HZ3(0,-1,-1)
     $  + HZ5(0, -1,-1,1,-1)
     $  + 2.0000000000000000d+00*HZ5(0,-1,1,-1,-1)
     $  + 3.0000000000000000d+00*HZ5(0,1,-1,-1,-1)
      HZ5(-1,1,0,-1,0) = 
     $  - HZ1( -1)*HZ1(0)*HZ3(-1,0,1)
     $  - HZ1( -1)*HZ1(0)*HZ3(0,-1,1)
     $  - HZ1( -1)*HZ4(0,-1,0,1)
     $  + HZ1(0) *HZ2(-1,1)*HZ2(0,-1)
     $  + 2.0000000000000000d+00*HZ1(0)*HZ4(-1,-1,0,1)
     $  + HZ1(0) *HZ4(-1,0,-1,1)
     $  - 2.0000000000000000d+00*HZ2(-1,1)*HZ3(0,0,-1)
     $  + HZ2(0, -1)*HZ3(-1,0,1)
     $  + HZ2(0, -1)*HZ3(0,-1,1)
     $  + HZ5(0, -1,0,1,-1)
      HZ5(-1,1,0,-1,1) = 
     $  - 2.0000000000000000d+00*HZ1(-1)*HZ4(0,-1,1,1)
     $  - HZ1( -1)*HZ4(0,1,-1,1)
     $  + HZ2( -1,1)*HZ3(0,-1,1)
     $  + HZ5(0, -1,1,-1,1)
     $  + 2.0000000000000000d+00*HZ5(0,-1,1,1,-1)
     $  + 2.0000000000000000d+00*HZ5(0,1,-1,-1,1)
     $  + HZ5(0,1, -1,1,-1)
      HZ5(-1,1,0,0,-1) = 
     $  + HZ1( -1)*HZ1(-1)*HZ3(0,0,1)
     $  - HZ1( -1)*HZ4(-1,0,0,1)
     $  + HZ2( -1,1)*HZ3(0,0,-1)
     $  - HZ2(0, -1)*HZ3(-1,0,1)
     $  - HZ2(0, -1)*HZ3(0,-1,1)
     $  - HZ5(0, -1,0,1,-1)
     $  - HZ5(0,0, -1,1,-1)
     $  - 2.0000000000000000d+00*HZ5(0,0,1,-1,-1)
      HZ5(-1,1,0,0,0) = 
     $  - HZ1( -1)*HZ4(0,0,0,1)
     $  + 1.6666666666666666d-01*HZ1(0)*HZ1(0)*HZ1(0)*HZ2(-1,1)
     $  - 5.0000000000000000d-01*HZ1(0)*HZ1(0)*HZ3(-1,0,1)
     $  - 5.0000000000000000d-01*HZ1(0)*HZ1(0)*HZ3(0,-1,1)
     $  + HZ1(0) *HZ4(-1,0,0,1)
     $  + HZ1(0) *HZ4(0,-1,0,1)
     $  + HZ1(0) *HZ4(0,0,-1,1)
     $  + HZ5(0,0,0,1, -1)
      HZ5(-1,1,0,0,1) = 
     $  + 5.0000000000000000d-01*HZ1(-1)*HZ2(0,1)*HZ2(0,1)
     $  + HZ2( -1,1)*HZ3(0,0,1)
     $  - HZ2(0,1) *HZ3(-1,0,1)
     $  - HZ2(0,1) *HZ3(0,-1,1)
     $  - HZ5(0,0,1, -1,1)
     $  - 2.0000000000000000d+00*HZ5(0,0,1,1,-1)
     $  - HZ5(0,1,0,1, -1)
      HZ5(-1,1,0,1,-1) = 
     $  - 2.0000000000000000d+00*HZ1(-1)*HZ1(-1)*HZ3(0,1,1)
     $  + HZ1( -1)*HZ2(-1,1)*HZ2(0,1)
     $  + 2.0000000000000000d+00*HZ1(-1)*HZ4(-1,0,1,1)
     $  + 2.0000000000000000d+00*HZ1(-1)*HZ4(0,-1,1,1)
     $  + HZ1( -1)*HZ4(0,1,-1,1)
     $  - HZ2( -1,1)*HZ3(-1,0,1)
     $  - HZ2( -1,1)*HZ3(0,-1,1)
     $  + 2.0000000000000000d+00*HZ5(0,1,-1,1,-1)
     $  + 4.0000000000000000d+00*HZ5(0,1,1,-1,-1)
      HZ5(-1,1,0,1,0) = 
     $  - 5.0000000000000000d-01*HZ1(-1)*HZ2(0,1)*HZ2(0,1)
     $  + 2.0000000000000000d+00*HZ1(-1)*HZ4(0,0,1,1)
     $  + HZ1(0) *HZ2(-1,1)*HZ2(0,1)
     $  - 2.0000000000000000d+00*HZ1(0)*HZ4(-1,0,1,1)
     $  - 2.0000000000000000d+00*HZ1(0)*HZ4(0,-1,1,1)
     $  - HZ1(0) *HZ4(0,1,-1,1)
     $  - 2.0000000000000000d+00*HZ2(-1,1)*HZ3(0,0,1)
     $  + HZ2(0,1) *HZ3(-1,0,1)
     $  + HZ2(0,1) *HZ3(0,-1,1)
     $  + HZ5(0,1,0,1, -1)
      HZ5(-1,1,0,1,1) = 
     $  - 3.0000000000000000d+00*HZ1(-1)*HZ4(0,1,1,1)
     $  + HZ2( -1,1)*HZ3(0,1,1)
     $  + HZ5(0,1, -1,1,1)
     $  + 2.0000000000000000d+00*HZ5(0,1,1,-1,1)
     $  + 3.0000000000000000d+00*HZ5(0,1,1,1,-1)
      HZ5(-1,1,1,-1,-1) = 
     $  + 5.0000000000000000d-01*HZ1(-1)*HZ1(-1)*HZ3(-1,1,1)
     $  - 5.0000000000000000d-01*HZ1(-1)*HZ2(-1,1)*HZ2(-1,1)
     $  + HZ2( -1,1)*HZ3(-1,-1,1)
     $  - 3.0000000000000000d+00*HZ5(-1,-1,-1,1,1)
     $  - HZ5( -1,-1,1,-1,1)
      HZ5(-1,1,1,-1,0) = 
     $  + HZ1( -1)*HZ1(0)*HZ3(-1,1,1)
     $  - HZ1( -1)*HZ4(0,-1,1,1)
     $  - 5.0000000000000000d-01*HZ1(0)*HZ2(-1,1)*HZ2(-1,1)
     $  + HZ2( -1,1)*HZ3(0,-1,1)
     $  - HZ2(0, -1)*HZ3(-1,1,1)
     $  + HZ5(0, -1,1,1,-1)
      HZ5(-1,1,1,-1,1) = 
     $  + HZ2( -1,1)*HZ3(-1,1,1)
     $  - 6.0000000000000000d+00*HZ5(-1,-1,1,1,1)
     $  - 3.0000000000000000d+00*HZ5(-1,1,-1,1,1)
      HZ5(-1,1,1,0,-1) = 
     $  + HZ1( -1)*HZ1(-1)*HZ3(0,1,1)
     $  - HZ1( -1)*HZ2(-1,1)*HZ2(0,1)
     $  - HZ1( -1)*HZ4(-1,0,1,1)
     $  + HZ2( -1,1)*HZ3(-1,0,1)
     $  + HZ2(0, -1)*HZ3(-1,1,1)
     $  - HZ5(0, -1,1,1,-1)
     $  - HZ5(0,1, -1,1,-1)
     $  - 2.0000000000000000d+00*HZ5(0,1,1,-1,-1)
      HZ5(-1,1,1,0,0) = 
     $  - HZ1( -1)*HZ4(0,0,1,1)
     $  + 5.0000000000000000d-01*HZ1(0)*HZ1(0)*HZ3(-1,1,1)
     $  - HZ1(0) *HZ2(-1,1)*HZ2(0,1)
     $  + HZ1(0) *HZ4(-1,0,1,1)
     $  + HZ1(0) *HZ4(0,-1,1,1)
     $  + HZ1(0) *HZ4(0,1,-1,1)
     $  + HZ2( -1,1)*HZ3(0,0,1)
     $  + HZ5(0,0,1,1, -1)
      HZ5(-1,1,1,0,1) = 
     $  + 3.0000000000000000d+00*HZ1(-1)*HZ4(0,1,1,1)
     $  - 2.0000000000000000d+00*HZ2(-1,1)*HZ3(0,1,1)
     $  + HZ2(0,1) *HZ3(-1,1,1)
     $  - HZ5(0,1,1, -1,1)
     $  - 3.0000000000000000d+00*HZ5(0,1,1,1,-1)
      HZ5(-1,1,1,1,-1) = 
     $  + HZ1( -1)*HZ4(-1,1,1,1)
     $  - HZ2( -1,1)*HZ3(-1,1,1)
     $  + 4.0000000000000000d+00*HZ5(-1,-1,1,1,1)
     $  + 2.0000000000000000d+00*HZ5(-1,1,-1,1,1)
      HZ5(-1,1,1,1,0) = 
     $  - HZ1( -1)*HZ4(0,1,1,1)
     $  + HZ1(0) *HZ4(-1,1,1,1)
     $  + HZ2( -1,1)*HZ3(0,1,1)
     $  - HZ2(0,1) *HZ3(-1,1,1)
     $  + HZ5(0,1,1,1, -1)
      HZ5(0,-1,-1,1,0) = 
     $  + HZ1(0) *HZ4(0,-1,-1,1)
     $  - HZ5(0, -1,-1,0,1)
     $  - HZ5(0, -1,0,-1,1)
     $  - 2.0000000000000000d+00*HZ5(0,0,-1,-1,1)
      HZ5(0,-1,0,0,1) = 
     $  + HZ2(0, -1)*HZ3(0,0,1)
     $  - 2.0000000000000000d+00*HZ5(0,0,-1,0,1)
     $  - 3.0000000000000000d+00*HZ5(0,0,0,-1,1)
     $  - 3.0000000000000000d+00*HZ5(0,0,0,1,-1)
     $  - HZ5(0,0,1,0, -1)
      HZ5(0,-1,0,1,0) = 
     $  + HZ1(0) *HZ4(0,-1,0,1)
     $  - 2.0000000000000000d+00*HZ2(0,-1)*HZ3(0,0,1)
     $  + 2.0000000000000000d+00*HZ5(0,0,-1,0,1)
     $  + 6.0000000000000000d+00*HZ5(0,0,0,-1,1)
     $  + 6.0000000000000000d+00*HZ5(0,0,0,1,-1)
     $  + 2.0000000000000000d+00*HZ5(0,0,1,0,-1)
      HZ5(0,-1,1,-1,0) = 
     $  + HZ1( -1)*HZ1(0)*HZ3(0,-1,1)
     $  - HZ1(0) *HZ4(-1,0,-1,1)
     $  - 2.0000000000000000d+00*HZ1(0)*HZ4(0,-1,-1,1)
     $  - HZ2(0, -1)*HZ3(0,-1,1)
     $  + 2.0000000000000000d+00*HZ5(0,-1,0,-1,1)
     $  + 4.0000000000000000d+00*HZ5(0,0,-1,-1,1)
      HZ5(0,-1,1,0,-1) = 
     $  + HZ2(0, -1)*HZ3(0,-1,1)
     $  - 2.0000000000000000d+00*HZ5(0,-1,0,-1,1)
     $  - HZ5(0, -1,0,1,-1)
     $  - 4.0000000000000000d+00*HZ5(0,0,-1,-1,1)
     $  - 2.0000000000000000d+00*HZ5(0,0,-1,1,-1)
      HZ5(0,-1,1,0,0) = 
     $  + 5.0000000000000000d-01*HZ1(0)*HZ1(0)*HZ3(0,-1,1)
     $  - HZ1(0) *HZ4(0,-1,0,1)
     $  - 2.0000000000000000d+00*HZ1(0)*HZ4(0,0,-1,1)
     $  + HZ2(0, -1)*HZ3(0,0,1)
     $  - 3.0000000000000000d+00*HZ5(0,0,0,1,-1)
     $  - HZ5(0,0,1,0, -1)
      HZ5(0,-1,1,1,0) = 
     $  + HZ1(0) *HZ4(0,-1,1,1)
     $  - HZ5(0, -1,0,1,1)
     $  - HZ5(0, -1,1,0,1)
     $  - 2.0000000000000000d+00*HZ5(0,0,-1,1,1)
      HZ5(0,0,-1,1,0) = 
     $  + HZ1(0) *HZ4(0,0,-1,1)
     $  - HZ5(0,0, -1,0,1)
     $  - 3.0000000000000000d+00*HZ5(0,0,0,-1,1)
      HZ5(0,0,1,-1,0) = 
     $  + HZ1( -1)*HZ1(0)*HZ3(0,0,1)
     $  - HZ1(0) *HZ4(-1,0,0,1)
     $  - HZ1(0) *HZ4(0,-1,0,1)
     $  - HZ1(0) *HZ4(0,0,-1,1)
     $  - 3.0000000000000000d+00*HZ5(0,0,0,1,-1)
     $  - HZ5(0,0,1,0, -1)
      HZ5(0,1,-1,-1,0) = 
     $  + 5.0000000000000000d-01*HZ1(-1)*HZ1(-1)*HZ1(0)*HZ2(0,1)
     $  - HZ1( -1)*HZ1(0)*HZ3(-1,0,1)
     $  - HZ1( -1)*HZ1(0)*HZ3(0,-1,1)
     $  - HZ1( -1)*HZ2(0,-1)*HZ2(0,1)
     $  + HZ1(0) *HZ4(-1,-1,0,1)
     $  + HZ1(0) *HZ4(-1,0,-1,1)
     $  + HZ1(0) *HZ4(0,-1,-1,1)
     $  + HZ2(0, -1)*HZ3(-1,0,1)
     $  + HZ2(0, -1)*HZ3(0,-1,1)
     $  + HZ2(0,1) *HZ3(0,-1,-1)
     $  - HZ5(0, -1,-1,0,1)
     $  - HZ5(0, -1,0,-1,1)
     $  - 2.0000000000000000d+00*HZ5(0,0,-1,-1,1)
      HZ5(0,1,-1,0,-1) = 
     $  + HZ1( -1)*HZ2(0,-1)*HZ2(0,1)
     $  - HZ2(0, -1)*HZ3(-1,0,1)
     $  - HZ2(0, -1)*HZ3(0,-1,1)
     $  - 2.0000000000000000d+00*HZ2(0,1)*HZ3(0,-1,-1)
     $  + 2.0000000000000000d+00*HZ5(0,-1,-1,0,1)
     $  + 2.0000000000000000d+00*HZ5(0,-1,0,-1,1)
     $  + HZ5(0, -1,0,1,-1)
     $  + 4.0000000000000000d+00*HZ5(0,0,-1,-1,1)
     $  + 2.0000000000000000d+00*HZ5(0,0,-1,1,-1)
      HZ5(0,1,-1,0,0) = 
     $  + 5.0000000000000000d-01*HZ1(-1)*HZ1(0)*HZ1(0)*HZ2(0,1)
     $  - 5.0000000000000000d-01*HZ1(0)*HZ1(0)*HZ3(-1,0,1)
     $  - 5.0000000000000000d-01*HZ1(0)*HZ1(0)*HZ3(0,-1,1)
     $  - HZ1(0) *HZ2(0,-1)*HZ2(0,1)
     $  + HZ1(0) *HZ4(0,-1,0,1)
     $  + 2.0000000000000000d+00*HZ1(0)*HZ4(0,0,-1,1)
     $  + HZ2(0,1) *HZ3(0,0,-1)
     $  - HZ5(0,0, -1,0,1)
     $  - 3.0000000000000000d+00*HZ5(0,0,0,-1,1)
      HZ5(0,1,-1,0,1) = 
     $  + HZ1( -1)*HZ2(0,1)*HZ2(0,1)
     $  - HZ2(0,1) *HZ3(-1,0,1)
     $  - 2.0000000000000000d+00*HZ2(0,1)*HZ3(0,-1,1)
     $  + 2.0000000000000000d+00*HZ5(0,-1,0,1,1)
     $  + HZ5(0, -1,1,0,1)
     $  + 4.0000000000000000d+00*HZ5(0,0,-1,1,1)
     $  - 4.0000000000000000d+00*HZ5(0,0,1,1,-1)
     $  - 2.0000000000000000d+00*HZ5(0,1,0,1,-1)
      HZ5(0,1,-1,1,0) = 
     $  - HZ1( -1)*HZ2(0,1)*HZ2(0,1)
     $  + HZ1(0) *HZ4(0,1,-1,1)
     $  + HZ2(0,1) *HZ3(-1,0,1)
     $  + HZ2(0,1) *HZ3(0,-1,1)
     $  + 4.0000000000000000d+00*HZ5(0,0,1,1,-1)
     $  + 2.0000000000000000d+00*HZ5(0,1,0,1,-1)
      HZ5(0,1,0,-1,-1) = 
     $  + HZ2(0,1) *HZ3(0,-1,-1)
     $  - HZ5(0, -1,-1,0,1)
     $  - HZ5(0, -1,0,-1,1)
     $  - HZ5(0, -1,0,1,-1)
     $  - 2.0000000000000000d+00*HZ5(0,0,-1,-1,1)
     $  - 2.0000000000000000d+00*HZ5(0,0,-1,1,-1)
     $  - 2.0000000000000000d+00*HZ5(0,0,1,-1,-1)
      HZ5(0,1,0,-1,0) = 
     $  - 2.0000000000000000d+00*HZ1(-1)*HZ1(0)*HZ3(0,0,1)
     $  + HZ1(0) *HZ2(0,-1)*HZ2(0,1)
     $  + 2.0000000000000000d+00*HZ1(0)*HZ4(-1,0,0,1)
     $  + HZ1(0) *HZ4(0,-1,0,1)
     $  - 2.0000000000000000d+00*HZ2(0,1)*HZ3(0,0,-1)
     $  + 2.0000000000000000d+00*HZ5(0,0,-1,0,1)
     $  + 6.0000000000000000d+00*HZ5(0,0,0,-1,1)
     $  + 6.0000000000000000d+00*HZ5(0,0,0,1,-1)
     $  + 2.0000000000000000d+00*HZ5(0,0,1,0,-1)
      HZ5(0,1,0,-1,1) = 
     $  + HZ2(0,1) *HZ3(0,-1,1)
     $  - 2.0000000000000000d+00*HZ5(0,-1,0,1,1)
     $  - HZ5(0, -1,1,0,1)
     $  - 4.0000000000000000d+00*HZ5(0,0,-1,1,1)
     $  - 2.0000000000000000d+00*HZ5(0,0,1,-1,1)
      HZ5(0,1,0,0,-1) = 
     $  + HZ2(0,1) *HZ3(0,0,-1)
     $  - HZ5(0,0, -1,0,1)
     $  - 3.0000000000000000d+00*HZ5(0,0,0,-1,1)
     $  - 3.0000000000000000d+00*HZ5(0,0,0,1,-1)
     $  - 2.0000000000000000d+00*HZ5(0,0,1,0,-1)
      HZ5(0,1,1,-1,0) = 
     $  + HZ1( -1)*HZ1(0)*HZ3(0,1,1)
     $  - HZ1(0) *HZ4(-1,0,1,1)
     $  - HZ1(0) *HZ4(0,-1,1,1)
     $  - HZ1(0) *HZ4(0,1,-1,1)
     $  - HZ2(0, -1)*HZ3(0,1,1)
     $  + HZ2(0,1) *HZ3(0,-1,1)
     $  - HZ5(0, -1,0,1,1)
     $  - HZ5(0, -1,1,0,1)
     $  - 2.0000000000000000d+00*HZ5(0,0,-1,1,1)
      HZ5(0,1,1,0,-1) = 
     $  + HZ2(0, -1)*HZ3(0,1,1)
     $  - HZ2(0,1) *HZ3(0,-1,1)
     $  + HZ5(0, -1,0,1,1)
     $  + HZ5(0, -1,1,0,1)
     $  + 2.0000000000000000d+00*HZ5(0,0,-1,1,1)
     $  - 2.0000000000000000d+00*HZ5(0,0,1,1,-1)
     $  - HZ5(0,1,0,1, -1)
      HZ5(1,-1,-1,-1,-1) = 
     $  + 4.1666666666666666d-02*HZ1(-1)*HZ1(-1)*HZ1(-1)*HZ1(-1)*HZ1(1)
     $  - 1.6666666666666666d-01*HZ1(-1)*HZ1(-1)*HZ1(-1)*HZ2(-1,1)
     $  + 5.0000000000000000d-01*HZ1(-1)*HZ1(-1)*HZ3(-1,-1,1)
     $  - HZ1( -1)*HZ4(-1,-1,-1,1)
     $  + HZ5( -1,-1,-1,-1,1)
      HZ5(1,-1,-1,-1,0) = 
     $  + 1.6666666666666666d-01*HZ1(-1)*HZ1(-1)*HZ1(-1)*HZ1(0)*HZ1(1)
     $  - 5.0000000000000000d-01*HZ1(-1)*HZ1(-1)*HZ1(0)*HZ2(-1,1)
     $  - 5.0000000000000000d-01*HZ1(-1)*HZ1(-1)*HZ1(1)*HZ2(0,-1)
     $  + HZ1( -1)*HZ1(0)*HZ3(-1,-1,1)
     $  + HZ1( -1)*HZ1(1)*HZ3(0,-1,-1)
     $  + HZ1( -1)*HZ2(-1,1)*HZ2(0,-1)
     $  - HZ1(0) *HZ4(-1,-1,-1,1)
     $  - HZ1(1) *HZ4(0,-1,-1,-1)
     $  - HZ2( -1,1)*HZ3(0,-1,-1)
     $  - HZ2(0, -1)*HZ3(-1,-1,1)
     $  + HZ5(0, -1,-1,-1,1)
      HZ5(1,-1,-1,-1,1) = 
     $  + HZ1(1) *HZ4(-1,-1,-1,1)
     $  - HZ2( -1,1)*HZ3(-1,-1,1)
     $  + 4.0000000000000000d+00*HZ5(-1,-1,-1,1,1)
     $  + 2.0000000000000000d+00*HZ5(-1,-1,1,-1,1)
      HZ5(1,-1,-1,0,-1) = 
     $  + 5.0000000000000000d-01*HZ1(-1)*HZ1(-1)*HZ1(1)*HZ2(0,-1)
     $  - 2.0000000000000000d+00*HZ1(-1)*HZ1(1)*HZ3(0,-1,-1)
     $  - HZ1( -1)*HZ2(-1,1)*HZ2(0,-1)
     $  + 3.0000000000000000d+00*HZ1(1)*HZ4(0,-1,-1,-1)
     $  + 2.0000000000000000d+00*HZ2(-1,1)*HZ3(0,-1,-1)
     $  + HZ2(0, -1)*HZ3(-1,-1,1)
     $  - 3.0000000000000000d+00*HZ5(0,-1,-1,-1,1)
     $  - HZ5(0, -1,-1,1,-1)
      HZ5(1,-1,-1,0,0) = 
     $  + 2.5000000000000000d-01*HZ1(-1)*HZ1(-1)*HZ1(0)*HZ1(0)*HZ1(1)
     $  - 5.0000000000000000d-01*HZ1(-1)*HZ1(0)*HZ1(0)*HZ2(-1,1)
     $  - HZ1( -1)*HZ1(0)*HZ1(1)*HZ2(0,-1)
     $  + HZ1( -1)*HZ1(1)*HZ3(0,0,-1)
     $  + 5.0000000000000000d-01*HZ1(0)*HZ1(0)*HZ3(-1,-1,1)
     $  + HZ1(0) *HZ1(1)*HZ3(0,-1,-1)
     $  + HZ1(0) *HZ2(-1,1)*HZ2(0,-1)
     $  - HZ1(0) *HZ4(0,-1,-1,1)
     $  - HZ1(1) *HZ4(0,0,-1,-1)
     $  - HZ2( -1,1)*HZ3(0,0,-1)
     $  + HZ5(0,0, -1,-1,1)
      HZ5(1,-1,-1,0,1) = 
     $  + HZ1(1) *HZ4(-1,-1,0,1)
     $  - HZ2( -1,1)*HZ3(-1,0,1)
     $  + HZ2(0,1) *HZ3(-1,-1,1)
     $  - 2.0000000000000000d+00*HZ5(0,-1,-1,1,1)
     $  - HZ5(0, -1,1,-1,1)
     $  - HZ5(0,1, -1,-1,1)
      HZ5(1,-1,-1,1,-1) = 
     $  + HZ1( -1)*HZ1(1)*HZ3(-1,-1,1)
     $  - 5.0000000000000000d-01*HZ1(-1)*HZ2(-1,1)*HZ2(-1,1)
     $  - 3.0000000000000000d+00*HZ1(1)*HZ4(-1,-1,-1,1)
     $  + 2.0000000000000000d+00*HZ2(-1,1)*HZ3(-1,-1,1)
     $  - 6.0000000000000000d+00*HZ5(-1,-1,-1,1,1)
     $  - 3.0000000000000000d+00*HZ5(-1,-1,1,-1,1)
      HZ5(1,-1,-1,1,0) = 
     $  + HZ1(0) *HZ1(1)*HZ3(-1,-1,1)
     $  - 5.0000000000000000d-01*HZ1(0)*HZ2(-1,1)*HZ2(-1,1)
     $  - HZ1(1) *HZ4(-1,-1,0,1)
     $  - HZ1(1) *HZ4(-1,0,-1,1)
     $  - HZ1(1) *HZ4(0,-1,-1,1)
     $  + HZ2( -1,1)*HZ3(-1,0,1)
     $  + HZ2( -1,1)*HZ3(0,-1,1)
     $  - HZ2(0,1) *HZ3(-1,-1,1)
     $  + HZ5(0,1, -1,-1,1)
      HZ5(1,-1,-1,1,1) = 
     $  + HZ1(1) *HZ4(-1,-1,1,1)
     $  - 3.0000000000000000d+00*HZ5(-1,-1,1,1,1)
     $  - HZ5( -1,1,-1,1,1)
      HZ5(1,-1,0,-1,-1) = 
     $  + HZ1( -1)*HZ1(1)*HZ3(0,-1,-1)
     $  - 3.0000000000000000d+00*HZ1(1)*HZ4(0,-1,-1,-1)
     $  - HZ2( -1,1)*HZ3(0,-1,-1)
     $  + 3.0000000000000000d+00*HZ5(0,-1,-1,-1,1)
     $  + 2.0000000000000000d+00*HZ5(0,-1,-1,1,-1)
     $  + HZ5(0, -1,1,-1,-1)
      HZ5(1,-1,0,-1,0) = 
     $  + HZ1( -1)*HZ1(0)*HZ1(1)*HZ2(0,-1)
     $  + HZ1( -1)*HZ1(0)*HZ3(0,-1,1)
     $  - 2.0000000000000000d+00*HZ1(-1)*HZ1(1)*HZ3(0,0,-1)
     $  - 2.0000000000000000d+00*HZ1(0)*HZ1(1)*HZ3(0,-1,-1)
     $  - HZ1(0) *HZ2(-1,1)*HZ2(0,-1)
     $  - HZ1(0) *HZ4(-1,0,-1,1)
     $  + 5.0000000000000000d-01*HZ1(1)*HZ2(0,-1)*HZ2(0,-1)
     $  + 2.0000000000000000d+00*HZ1(1)*HZ4(0,0,-1,-1)
     $  + 2.0000000000000000d+00*HZ2(-1,1)*HZ3(0,0,-1)
     $  - HZ2(0, -1)*HZ3(0,-1,1)
     $  + HZ5(0, -1,0,-1,1)
      HZ5(1,-1,0,-1,1) = 
     $  + HZ1(1) *HZ4(-1,0,-1,1)
     $  - HZ2( -1,1)*HZ3(0,-1,1)
     $  + 4.0000000000000000d+00*HZ5(0,-1,-1,1,1)
     $  + 2.0000000000000000d+00*HZ5(0,-1,1,-1,1)
      HZ5(1,-1,0,0,-1) = 
     $  + HZ1( -1)*HZ1(1)*HZ3(0,0,-1)
     $  - 5.0000000000000000d-01*HZ1(1)*HZ2(0,-1)*HZ2(0,-1)
     $  - HZ2( -1,1)*HZ3(0,0,-1)
     $  + HZ2(0, -1)*HZ3(0,-1,1)
     $  - HZ5(0, -1,0,-1,1)
     $  - 2.0000000000000000d+00*HZ5(0,0,-1,-1,1)
     $  - HZ5(0,0, -1,1,-1)
      HZ5(1,-1,0,0,0) = 
     $  + 1.6666666666666666d-01*HZ1(-1)*HZ1(0)*HZ1(0)*HZ1(0)*HZ1(1)
     $  - 1.6666666666666666d-01*HZ1(0)*HZ1(0)*HZ1(0)*HZ2(-1,1)
     $  - 5.0000000000000000d-01*HZ1(0)*HZ1(0)*HZ1(1)*HZ2(0,-1)
     $  + 5.0000000000000000d-01*HZ1(0)*HZ1(0)*HZ3(0,-1,1)
     $  + HZ1(0) *HZ1(1)*HZ3(0,0,-1)
     $  - HZ1(0) *HZ4(0,0,-1,1)
     $  - HZ1(1) *HZ4(0,0,0,-1)
     $  + HZ5(0,0,0, -1,1)
      HZ5(1,-1,0,0,1) = 
     $  + HZ1(1) *HZ4(-1,0,0,1)
     $  - HZ2( -1,1)*HZ3(0,0,1)
     $  + 2.0000000000000000d+00*HZ5(0,-1,0,1,1)
     $  + HZ5(0, -1,1,0,1)
     $  + 2.0000000000000000d+00*HZ5(0,0,-1,1,1)
     $  + HZ5(0,0,1, -1,1)
      HZ5(1,-1,0,1,-1) = 
     $  + HZ1( -1)*HZ1(1)*HZ3(-1,0,1)
     $  - HZ1( -1)*HZ2(-1,1)*HZ2(0,1)
     $  - 2.0000000000000000d+00*HZ1(1)*HZ4(-1,-1,0,1)
     $  - HZ1(1) *HZ4(-1,0,-1,1)
     $  + HZ2( -1,1)*HZ3(-1,0,1)
     $  + HZ2( -1,1)*HZ3(0,-1,1)
     $  + HZ5(0, -1,1,-1,1)
     $  + 2.0000000000000000d+00*HZ5(0,-1,1,1,-1)
     $  + 2.0000000000000000d+00*HZ5(0,1,-1,-1,1)
     $  + HZ5(0,1, -1,1,-1)
      HZ5(1,-1,0,1,0) = 
     $  + HZ1(0) *HZ1(1)*HZ3(-1,0,1)
     $  - HZ1(0) *HZ2(-1,1)*HZ2(0,1)
     $  + 2.0000000000000000d+00*HZ1(0)*HZ4(0,-1,1,1)
     $  + HZ1(0) *HZ4(0,1,-1,1)
     $  - 2.0000000000000000d+00*HZ1(1)*HZ4(-1,0,0,1)
     $  - HZ1(1) *HZ4(0,-1,0,1)
     $  + 2.0000000000000000d+00*HZ2(-1,1)*HZ3(0,0,1)
     $  - 2.0000000000000000d+00*HZ5(0,-1,0,1,1)
     $  - HZ5(0, -1,1,0,1)
     $  - 4.0000000000000000d+00*HZ5(0,0,-1,1,1)
     $  - 2.0000000000000000d+00*HZ5(0,0,1,-1,1)
      HZ5(1,-1,0,1,1) = 
     $  + HZ1(1) *HZ4(-1,0,1,1)
     $  - HZ2( -1,1)*HZ3(0,1,1)
     $  + 3.0000000000000000d+00*HZ5(0,-1,1,1,1)
     $  + 2.0000000000000000d+00*HZ5(0,1,-1,1,1)
     $  + HZ5(0,1,1, -1,1)
      HZ5(1,-1,1,-1,-1) = 
     $  + 5.0000000000000000d-01*HZ1(-1)*HZ1(-1)*HZ1(1)*HZ2(-1,1)
     $  - HZ1( -1)*HZ1(-1)*HZ3(-1,1,1)
     $  - 2.0000000000000000d+00*HZ1(-1)*HZ1(1)*HZ3(-1,-1,1)
     $  + 5.0000000000000000d-01*HZ1(-1)*HZ2(-1,1)*HZ2(-1,1)
     $  + 2.0000000000000000d+00*HZ1(-1)*HZ4(-1,-1,1,1)
     $  + 3.0000000000000000d+00*HZ1(1)*HZ4(-1,-1,-1,1)
     $  - HZ2( -1,1)*HZ3(-1,-1,1)
     $  + HZ5( -1,-1,1,-1,1)
      HZ5(1,-1,1,-1,0) = 
     $  + HZ1( -1)*HZ1(0)*HZ1(1)*HZ2(-1,1)
     $  - 2.0000000000000000d+00*HZ1(-1)*HZ1(0)*HZ3(-1,1,1)
     $  - 2.0000000000000000d+00*HZ1(0)*HZ1(1)*HZ3(-1,-1,1)
     $  + 5.0000000000000000d-01*HZ1(0)*HZ2(-1,1)*HZ2(-1,1)
     $  + 2.0000000000000000d+00*HZ1(0)*HZ4(-1,-1,1,1)
     $  - HZ1(1) *HZ2(-1,1)*HZ2(0,-1)
     $  + HZ1(1) *HZ4(-1,0,-1,1)
     $  + 2.0000000000000000d+00*HZ1(1)*HZ4(0,-1,-1,1)
     $  - HZ2( -1,1)*HZ3(0,-1,1)
     $  + 2.0000000000000000d+00*HZ2(0,-1)*HZ3(-1,1,1)
     $  + HZ5(0, -1,1,-1,1)
      HZ5(1,-1,1,-1,1) = 
     $  + 5.0000000000000000d-01*HZ1(1)*HZ2(-1,1)*HZ2(-1,1)
     $  - 2.0000000000000000d+00*HZ1(1)*HZ4(-1,-1,1,1)
     $  - 2.0000000000000000d+00*HZ2(-1,1)*HZ3(-1,1,1)
     $  + 1.2000000000000000d+01*HZ5(-1,-1,1,1,1)
     $  + 4.0000000000000000d+00*HZ5(-1,1,-1,1,1)
      HZ5(1,-1,1,0,-1) = 
     $  - HZ1( -1)*HZ1(1)*HZ3(-1,0,1)
     $  - HZ1( -1)*HZ1(1)*HZ3(0,-1,1)
     $  + HZ1( -1)*HZ2(-1,1)*HZ2(0,1)
     $  + HZ1(1) *HZ2(-1,1)*HZ2(0,-1)
     $  + 2.0000000000000000d+00*HZ1(1)*HZ4(-1,-1,0,1)
     $  + HZ1(1) *HZ4(-1,0,-1,1)
     $  - HZ2( -1,1)*HZ3(-1,0,1)
     $  - 2.0000000000000000d+00*HZ2(0,-1)*HZ3(-1,1,1)
     $  - HZ5(0, -1,1,-1,1)
     $  - 2.0000000000000000d+00*HZ5(0,1,-1,-1,1)
     $  - HZ5(0,1, -1,1,-1)
      HZ5(1,-1,1,0,0) = 
     $  + 5.0000000000000000d-01*HZ1(0)*HZ1(0)*HZ1(1)*HZ2(-1,1)
     $  - HZ1(0) *HZ1(0)*HZ3(-1,1,1)
     $  - HZ1(0) *HZ1(1)*HZ3(-1,0,1)
     $  - HZ1(0) *HZ1(1)*HZ3(0,-1,1)
     $  + HZ1(0) *HZ2(-1,1)*HZ2(0,1)
     $  - HZ1(0) *HZ4(0,1,-1,1)
     $  + HZ1(1) *HZ4(-1,0,0,1)
     $  + HZ1(1) *HZ4(0,-1,0,1)
     $  + HZ1(1) *HZ4(0,0,-1,1)
     $  - HZ2( -1,1)*HZ3(0,0,1)
     $  + HZ5(0,0,1, -1,1)
      HZ5(1,-1,1,0,1) = 
     $  + HZ1(1) *HZ2(-1,1)*HZ2(0,1)
     $  - 2.0000000000000000d+00*HZ1(1)*HZ4(-1,0,1,1)
     $  - 2.0000000000000000d+00*HZ1(1)*HZ4(0,-1,1,1)
     $  - HZ1(1) *HZ4(0,1,-1,1)
     $  + 2.0000000000000000d+00*HZ2(-1,1)*HZ3(0,1,1)
     $  - 2.0000000000000000d+00*HZ2(0,1)*HZ3(-1,1,1)
     $  - 2.0000000000000000d+00*HZ5(0,1,-1,1,1)
     $  - 2.0000000000000000d+00*HZ5(0,1,1,-1,1)
      HZ5(1,-1,1,1,-1) = 
     $  + HZ1( -1)*HZ1(1)*HZ3(-1,1,1)
     $  - 3.0000000000000000d+00*HZ1(-1)*HZ4(-1,1,1,1)
     $  - 5.0000000000000000d-01*HZ1(1)*HZ2(-1,1)*HZ2(-1,1)
     $  + 2.0000000000000000d+00*HZ2(-1,1)*HZ3(-1,1,1)
     $  - 6.0000000000000000d+00*HZ5(-1,-1,1,1,1)
     $  - 3.0000000000000000d+00*HZ5(-1,1,-1,1,1)
      HZ5(1,-1,1,1,0) = 
     $  + HZ1(0) *HZ1(1)*HZ3(-1,1,1)
     $  - 3.0000000000000000d+00*HZ1(0)*HZ4(-1,1,1,1)
     $  - HZ1(1) *HZ2(-1,1)*HZ2(0,1)
     $  + HZ1(1) *HZ4(-1,0,1,1)
     $  + HZ1(1) *HZ4(0,-1,1,1)
     $  + HZ1(1) *HZ4(0,1,-1,1)
     $  - HZ2( -1,1)*HZ3(0,1,1)
     $  + 2.0000000000000000d+00*HZ2(0,1)*HZ3(-1,1,1)
     $  + HZ5(0,1,1, -1,1)
      HZ5(1,-1,1,1,1) = 
     $  + HZ1(1) *HZ4(-1,1,1,1)
     $  - 4.0000000000000000d+00*HZ5(-1,1,1,1,1)
      HZ5(1,0,-1,-1,-1) = 
     $  + HZ1(1) *HZ4(0,-1,-1,-1)
     $  - HZ5(0, -1,-1,-1,1)
     $  - HZ5(0, -1,-1,1,-1)
     $  - HZ5(0, -1,1,-1,-1)
     $  - HZ5(0,1, -1,-1,-1)
      HZ5(1,0,-1,-1,0) = 
     $  - 5.0000000000000000d-01*HZ1(-1)*HZ1(-1)*HZ1(0)*HZ2(0,1)
     $  + HZ1( -1)*HZ1(0)*HZ3(-1,0,1)
     $  + HZ1( -1)*HZ2(0,-1)*HZ2(0,1)
     $  + HZ1(0) *HZ1(1)*HZ3(0,-1,-1)
     $  - HZ1(0) *HZ4(-1,-1,0,1)
     $  - 5.0000000000000000d-01*HZ1(1)*HZ2(0,-1)*HZ2(0,-1)
     $  - HZ2(0, -1)*HZ3(-1,0,1)
     $  - HZ2(0,1) *HZ3(0,-1,-1)
     $  + HZ5(0, -1,-1,0,1)
      HZ5(1,0,-1,-1,1) = 
     $  + HZ1(1) *HZ4(0,-1,-1,1)
     $  - 2.0000000000000000d+00*HZ5(0,-1,-1,1,1)
     $  - HZ5(0, -1,1,-1,1)
     $  - HZ5(0,1, -1,-1,1)
      HZ5(1,0,-1,0,-1) = 
     $  - HZ1( -1)*HZ2(0,-1)*HZ2(0,1)
     $  + 5.0000000000000000d-01*HZ1(1)*HZ2(0,-1)*HZ2(0,-1)
     $  - 2.0000000000000000d+00*HZ1(1)*HZ4(0,0,-1,-1)
     $  + HZ2(0, -1)*HZ3(-1,0,1)
     $  + 2.0000000000000000d+00*HZ2(0,1)*HZ3(0,-1,-1)
     $  - 2.0000000000000000d+00*HZ5(0,-1,-1,0,1)
     $  - HZ5(0, -1,0,-1,1)
     $  - HZ5(0, -1,0,1,-1)
      HZ5(1,0,-1,0,0) = 
     $  - 5.0000000000000000d-01*HZ1(-1)*HZ1(0)*HZ1(0)*HZ2(0,1)
     $  + 5.0000000000000000d-01*HZ1(0)*HZ1(0)*HZ1(1)*HZ2(0,-1)
     $  + 5.0000000000000000d-01*HZ1(0)*HZ1(0)*HZ3(-1,0,1)
     $  - 2.0000000000000000d+00*HZ1(0)*HZ1(1)*HZ3(0,0,-1)
     $  + HZ1(0) *HZ2(0,-1)*HZ2(0,1)
     $  - HZ1(0) *HZ4(0,-1,0,1)
     $  + 3.0000000000000000d+00*HZ1(1)*HZ4(0,0,0,-1)
     $  - HZ2(0,1) *HZ3(0,0,-1)
     $  + HZ5(0,0, -1,0,1)
      HZ5(1,0,-1,0,1) = 
     $  - HZ1( -1)*HZ2(0,1)*HZ2(0,1)
     $  + HZ1(1) *HZ4(0,-1,0,1)
     $  + HZ2(0,1) *HZ3(-1,0,1)
     $  + 2.0000000000000000d+00*HZ2(0,1)*HZ3(0,-1,1)
     $  - 4.0000000000000000d+00*HZ5(0,-1,0,1,1)
     $  - 2.0000000000000000d+00*HZ5(0,-1,1,0,1)
     $  - 4.0000000000000000d+00*HZ5(0,0,-1,1,1)
     $  + 4.0000000000000000d+00*HZ5(0,0,1,1,-1)
     $  + 2.0000000000000000d+00*HZ5(0,1,0,1,-1)
      HZ5(1,0,-1,1,-1) = 
     $  + HZ1( -1)*HZ1(1)*HZ3(0,-1,1)
     $  - HZ1(1) *HZ4(-1,0,-1,1)
     $  - 2.0000000000000000d+00*HZ1(1)*HZ4(0,-1,-1,1)
     $  - HZ5(0, -1,1,-1,1)
     $  - 2.0000000000000000d+00*HZ5(0,-1,1,1,-1)
     $  - HZ5(0,1, -1,1,-1)
      HZ5(1,0,-1,1,0) = 
     $  + HZ1( -1)*HZ2(0,1)*HZ2(0,1)
     $  + HZ1(0) *HZ1(1)*HZ3(0,-1,1)
     $  - 2.0000000000000000d+00*HZ1(0)*HZ4(0,-1,1,1)
     $  - HZ1(0) *HZ4(0,1,-1,1)
     $  - HZ1(1) *HZ4(0,-1,0,1)
     $  - 2.0000000000000000d+00*HZ1(1)*HZ4(0,0,-1,1)
     $  - HZ2(0,1) *HZ3(-1,0,1)
     $  - HZ2(0,1) *HZ3(0,-1,1)
     $  + 2.0000000000000000d+00*HZ5(0,-1,0,1,1)
     $  + HZ5(0, -1,1,0,1)
     $  + 4.0000000000000000d+00*HZ5(0,0,-1,1,1)
     $  - 4.0000000000000000d+00*HZ5(0,0,1,1,-1)
     $  - 2.0000000000000000d+00*HZ5(0,1,0,1,-1)
      HZ5(1,0,-1,1,1) = 
     $  + HZ1(1) *HZ4(0,-1,1,1)
     $  - 3.0000000000000000d+00*HZ5(0,-1,1,1,1)
     $  - HZ5(0,1, -1,1,1)
      HZ5(1,0,0,-1,-1) = 
     $  + HZ1(1) *HZ4(0,0,-1,-1)
     $  - HZ2(0,1) *HZ3(0,-1,-1)
     $  + HZ5(0, -1,-1,0,1)
     $  + HZ5(0, -1,0,-1,1)
     $  + HZ5(0, -1,0,1,-1)
     $  + HZ5(0,0, -1,-1,1)
     $  + HZ5(0,0, -1,1,-1)
     $  + HZ5(0,0,1, -1,-1)
      HZ5(1,0,0,-1,0) = 
     $  + HZ1( -1)*HZ1(0)*HZ3(0,0,1)
     $  + HZ1(0) *HZ1(1)*HZ3(0,0,-1)
     $  - HZ1(0) *HZ2(0,-1)*HZ2(0,1)
     $  - HZ1(0) *HZ4(-1,0,0,1)
     $  - 3.0000000000000000d+00*HZ1(1)*HZ4(0,0,0,-1)
     $  + 2.0000000000000000d+00*HZ2(0,1)*HZ3(0,0,-1)
     $  - 2.0000000000000000d+00*HZ5(0,0,-1,0,1)
     $  - 3.0000000000000000d+00*HZ5(0,0,0,-1,1)
     $  - 3.0000000000000000d+00*HZ5(0,0,0,1,-1)
     $  - HZ5(0,0,1,0, -1)
      HZ5(1,0,0,-1,1) = 
     $  + HZ1(1) *HZ4(0,0,-1,1)
     $  - HZ2(0,1) *HZ3(0,-1,1)
     $  + 2.0000000000000000d+00*HZ5(0,-1,0,1,1)
     $  + HZ5(0, -1,1,0,1)
     $  + 2.0000000000000000d+00*HZ5(0,0,-1,1,1)
     $  + HZ5(0,0,1, -1,1)
      HZ5(1,0,0,0,-1) = 
     $  + HZ1(1) *HZ4(0,0,0,-1)
     $  - HZ2(0,1) *HZ3(0,0,-1)
     $  + HZ5(0,0, -1,0,1)
     $  + 2.0000000000000000d+00*HZ5(0,0,0,-1,1)
     $  + 2.0000000000000000d+00*HZ5(0,0,0,1,-1)
     $  + HZ5(0,0,1,0, -1)
      HZ5(1,0,0,1,-1) = 
     $  + HZ1( -1)*HZ1(1)*HZ3(0,0,1)
     $  - HZ1(1) *HZ4(-1,0,0,1)
     $  - HZ1(1) *HZ4(0,-1,0,1)
     $  - HZ1(1) *HZ4(0,0,-1,1)
     $  - HZ5(0,0,1, -1,1)
     $  - 2.0000000000000000d+00*HZ5(0,0,1,1,-1)
     $  - HZ5(0,1,0,1, -1)
      HZ5(1,0,1,-1,-1) = 
     $  + 5.0000000000000000d-01*HZ1(-1)*HZ1(-1)*HZ1(1)*HZ2(0,1)
     $  - HZ1( -1)*HZ1(1)*HZ3(-1,0,1)
     $  - HZ1( -1)*HZ1(1)*HZ3(0,-1,1)
     $  + HZ1(1) *HZ4(-1,-1,0,1)
     $  + HZ1(1) *HZ4(-1,0,-1,1)
     $  + HZ1(1) *HZ4(0,-1,-1,1)
     $  - HZ5(0,1, -1,-1,1)
     $  - HZ5(0,1, -1,1,-1)
     $  - 2.0000000000000000d+00*HZ5(0,1,1,-1,-1)
      HZ5(1,0,1,-1,0) = 
     $  + HZ1( -1)*HZ1(0)*HZ1(1)*HZ2(0,1)
     $  - 2.0000000000000000d+00*HZ1(-1)*HZ1(0)*HZ3(0,1,1)
     $  - HZ1(0) *HZ1(1)*HZ3(-1,0,1)
     $  - HZ1(0) *HZ1(1)*HZ3(0,-1,1)
     $  + 2.0000000000000000d+00*HZ1(0)*HZ4(-1,0,1,1)
     $  + 2.0000000000000000d+00*HZ1(0)*HZ4(0,-1,1,1)
     $  + HZ1(0) *HZ4(0,1,-1,1)
     $  - HZ1(1) *HZ2(0,-1)*HZ2(0,1)
     $  + HZ1(1) *HZ4(0,-1,0,1)
     $  + 2.0000000000000000d+00*HZ1(1)*HZ4(0,0,-1,1)
     $  + 2.0000000000000000d+00*HZ2(0,-1)*HZ3(0,1,1)
     $  - HZ2(0,1) *HZ3(0,-1,1)
     $  + HZ5(0, -1,1,0,1)
      HZ5(1,0,1,-1,1) = 
     $  + HZ1(1) *HZ4(0,1,-1,1)
     $  - 2.0000000000000000d+00*HZ5(0,1,-1,1,1)
     $  - 2.0000000000000000d+00*HZ5(0,1,1,-1,1)
      HZ5(1,0,1,0,-1) = 
     $  - 2.0000000000000000d+00*HZ1(-1)*HZ1(1)*HZ3(0,0,1)
     $  + HZ1(1) *HZ2(0,-1)*HZ2(0,1)
     $  + 2.0000000000000000d+00*HZ1(1)*HZ4(-1,0,0,1)
     $  + HZ1(1) *HZ4(0,-1,0,1)
     $  - 2.0000000000000000d+00*HZ2(0,-1)*HZ3(0,1,1)
     $  + HZ2(0,1) *HZ3(0,-1,1)
     $  - HZ5(0, -1,1,0,1)
     $  + 2.0000000000000000d+00*HZ5(0,0,1,-1,1)
     $  + 4.0000000000000000d+00*HZ5(0,0,1,1,-1)
     $  + HZ5(0,1,0,1, -1)
      HZ5(1,0,1,1,-1) = 
     $  + HZ1( -1)*HZ1(1)*HZ3(0,1,1)
     $  - HZ1(1) *HZ4(-1,0,1,1)
     $  - HZ1(1) *HZ4(0,-1,1,1)
     $  - HZ1(1) *HZ4(0,1,-1,1)
     $  - HZ5(0,1,1, -1,1)
     $  - 3.0000000000000000d+00*HZ5(0,1,1,1,-1)
      HZ5(1,1,-1,-1,-1) = 
     $  + 8.3333333333333333d-02*HZ1(-1)*HZ1(-1)*HZ1(-1)*HZ1(1)*HZ1(1)
     $  - 5.0000000000000000d-01*HZ1(-1)*HZ1(-1)*HZ1(1)*HZ2(-1,1)
     $  + 5.0000000000000000d-01*HZ1(-1)*HZ1(-1)*HZ3(-1,1,1)
     $  + HZ1( -1)*HZ1(1)*HZ3(-1,-1,1)
     $  - HZ1( -1)*HZ4(-1,-1,1,1)
     $  - HZ1(1) *HZ4(-1,-1,-1,1)
     $  + HZ5( -1,-1,-1,1,1)
      HZ5(1,1,-1,-1,0) = 
     $  + 2.5000000000000000d-01*HZ1(-1)*HZ1(-1)*HZ1(0)*HZ1(1)*HZ1(1)
     $  - HZ1( -1)*HZ1(0)*HZ1(1)*HZ2(-1,1)
     $  + HZ1( -1)*HZ1(0)*HZ3(-1,1,1)
     $  - 5.0000000000000000d-01*HZ1(-1)*HZ1(1)*HZ1(1)*HZ2(0,-1)
     $  + HZ1(0) *HZ1(1)*HZ3(-1,-1,1)
     $  - HZ1(0) *HZ4(-1,-1,1,1)
     $  + 5.0000000000000000d-01*HZ1(1)*HZ1(1)*HZ3(0,-1,-1)
     $  + HZ1(1) *HZ2(-1,1)*HZ2(0,-1)
     $  - HZ1(1) *HZ4(0,-1,-1,1)
     $  - HZ2(0, -1)*HZ3(-1,1,1)
     $  + HZ5(0, -1,-1,1,1)
      HZ5(1,1,-1,-1,1) = 
     $  + 5.0000000000000000d-01*HZ1(1)*HZ1(1)*HZ3(-1,-1,1)
     $  - 5.0000000000000000d-01*HZ1(1)*HZ2(-1,1)*HZ2(-1,1)
     $  + HZ2( -1,1)*HZ3(-1,1,1)
     $  - 3.0000000000000000d+00*HZ5(-1,-1,1,1,1)
     $  - HZ5( -1,1,-1,1,1)
      HZ5(1,1,-1,0,-1) = 
     $  + 5.0000000000000000d-01*HZ1(-1)*HZ1(1)*HZ1(1)*HZ2(0,-1)
     $  + HZ1( -1)*HZ1(1)*HZ3(0,-1,1)
     $  - HZ1(1) *HZ1(1)*HZ3(0,-1,-1)
     $  - HZ1(1) *HZ2(-1,1)*HZ2(0,-1)
     $  - HZ1(1) *HZ4(-1,0,-1,1)
     $  + HZ2(0, -1)*HZ3(-1,1,1)
     $  - 2.0000000000000000d+00*HZ5(0,-1,-1,1,1)
     $  - HZ5(0, -1,1,-1,1)
     $  - HZ5(0, -1,1,1,-1)
      HZ5(1,1,-1,0,0) = 
     $  + 2.5000000000000000d-01*HZ1(-1)*HZ1(0)*HZ1(0)*HZ1(1)*HZ1(1)
     $  - 5.0000000000000000d-01*HZ1(0)*HZ1(0)*HZ1(1)*HZ2(-1,1)
     $  + 5.0000000000000000d-01*HZ1(0)*HZ1(0)*HZ3(-1,1,1)
     $  - 5.0000000000000000d-01*HZ1(0)*HZ1(1)*HZ1(1)*HZ2(0,-1)
     $  + HZ1(0) *HZ1(1)*HZ3(0,-1,1)
     $  - HZ1(0) *HZ4(0,-1,1,1)
     $  + 5.0000000000000000d-01*HZ1(1)*HZ1(1)*HZ3(0,0,-1)
     $  - HZ1(1) *HZ4(0,0,-1,1)
     $  + HZ5(0,0, -1,1,1)
      HZ5(1,1,-1,0,1) = 
     $  + 5.0000000000000000d-01*HZ1(1)*HZ1(1)*HZ3(-1,0,1)
     $  - HZ1(1) *HZ2(-1,1)*HZ2(0,1)
     $  + 2.0000000000000000d+00*HZ1(1)*HZ4(0,-1,1,1)
     $  + HZ1(1) *HZ4(0,1,-1,1)
     $  + HZ2(0,1) *HZ3(-1,1,1)
     $  - 3.0000000000000000d+00*HZ5(0,-1,1,1,1)
     $  - HZ5(0,1, -1,1,1)
      HZ5(1,1,-1,1,-1) = 
     $  + 5.0000000000000000d-01*HZ1(-1)*HZ1(1)*HZ1(1)*HZ2(-1,1)
     $  - 2.0000000000000000d+00*HZ1(-1)*HZ1(1)*HZ3(-1,1,1)
     $  + 3.0000000000000000d+00*HZ1(-1)*HZ4(-1,1,1,1)
     $  - HZ1(1) *HZ1(1)*HZ3(-1,-1,1)
     $  + 5.0000000000000000d-01*HZ1(1)*HZ2(-1,1)*HZ2(-1,1)
     $  + 2.0000000000000000d+00*HZ1(1)*HZ4(-1,-1,1,1)
     $  - HZ2( -1,1)*HZ3(-1,1,1)
     $  + HZ5( -1,1,-1,1,1)
      HZ5(1,1,-1,1,0) = 
     $  + 5.0000000000000000d-01*HZ1(0)*HZ1(1)*HZ1(1)*HZ2(-1,1)
     $  - 2.0000000000000000d+00*HZ1(0)*HZ1(1)*HZ3(-1,1,1)
     $  + 3.0000000000000000d+00*HZ1(0)*HZ4(-1,1,1,1)
     $  - 5.0000000000000000d-01*HZ1(1)*HZ1(1)*HZ3(-1,0,1)
     $  - 5.0000000000000000d-01*HZ1(1)*HZ1(1)*HZ3(0,-1,1)
     $  + HZ1(1) *HZ2(-1,1)*HZ2(0,1)
     $  - HZ1(1) *HZ4(0,1,-1,1)
     $  - HZ2(0,1) *HZ3(-1,1,1)
     $  + HZ5(0,1, -1,1,1)
      HZ5(1,1,-1,1,1) = 
     $  + 5.0000000000000000d-01*HZ1(1)*HZ1(1)*HZ3(-1,1,1)
     $  - 3.0000000000000000d+00*HZ1(1)*HZ4(-1,1,1,1)
     $  + 6.0000000000000000d+00*HZ5(-1,1,1,1,1)
      HZ5(1,1,0,-1,-1) = 
     $  - 5.0000000000000000d-01*HZ1(-1)*HZ1(-1)*HZ1(1)*HZ2(0,1)
     $  + HZ1( -1)*HZ1(1)*HZ3(-1,0,1)
     $  + 5.0000000000000000d-01*HZ1(1)*HZ1(1)*HZ3(0,-1,-1)
     $  - HZ1(1) *HZ4(-1,-1,0,1)
     $  + HZ5(0, -1,-1,1,1)
     $  + HZ5(0, -1,1,-1,1)
     $  + HZ5(0, -1,1,1,-1)
     $  + HZ5(0,1, -1,-1,1)
     $  + HZ5(0,1, -1,1,-1)
     $  + HZ5(0,1,1, -1,-1)
      HZ5(1,1,0,-1,0) = 
     $  - HZ1( -1)*HZ1(0)*HZ1(1)*HZ2(0,1)
     $  + HZ1( -1)*HZ1(0)*HZ3(0,1,1)
     $  + 5.0000000000000000d-01*HZ1(0)*HZ1(1)*HZ1(1)*HZ2(0,-1)
     $  + HZ1(0) *HZ1(1)*HZ3(-1,0,1)
     $  - HZ1(0) *HZ4(-1,0,1,1)
     $  - HZ1(1) *HZ1(1)*HZ3(0,0,-1)
     $  + HZ1(1) *HZ2(0,-1)*HZ2(0,1)
     $  - HZ1(1) *HZ4(0,-1,0,1)
     $  - HZ2(0, -1)*HZ3(0,1,1)
     $  + HZ5(0, -1,0,1,1)
      HZ5(1,1,0,-1,1) = 
     $  + 5.0000000000000000d-01*HZ1(1)*HZ1(1)*HZ3(0,-1,1)
     $  - 2.0000000000000000d+00*HZ1(1)*HZ4(0,-1,1,1)
     $  - HZ1(1) *HZ4(0,1,-1,1)
     $  + 3.0000000000000000d+00*HZ5(0,-1,1,1,1)
     $  + 2.0000000000000000d+00*HZ5(0,1,-1,1,1)
     $  + HZ5(0,1,1, -1,1)
      HZ5(1,1,0,0,-1) = 
     $  + HZ1( -1)*HZ1(1)*HZ3(0,0,1)
     $  + 5.0000000000000000d-01*HZ1(1)*HZ1(1)*HZ3(0,0,-1)
     $  - HZ1(1) *HZ2(0,-1)*HZ2(0,1)
     $  - HZ1(1) *HZ4(-1,0,0,1)
     $  + HZ2(0, -1)*HZ3(0,1,1)
     $  - HZ5(0, -1,0,1,1)
     $  - HZ5(0,0, -1,1,1)
     $  - HZ5(0,0,1, -1,1)
     $  - HZ5(0,0,1,1, -1)
      HZ5(1,1,0,1,-1) = 
     $  + 5.0000000000000000d-01*HZ1(-1)*HZ1(1)*HZ1(1)*HZ2(0,1)
     $  - 2.0000000000000000d+00*HZ1(-1)*HZ1(1)*HZ3(0,1,1)
     $  - 5.0000000000000000d-01*HZ1(1)*HZ1(1)*HZ3(-1,0,1)
     $  - 5.0000000000000000d-01*HZ1(1)*HZ1(1)*HZ3(0,-1,1)
     $  + 2.0000000000000000d+00*HZ1(1)*HZ4(-1,0,1,1)
     $  + 2.0000000000000000d+00*HZ1(1)*HZ4(0,-1,1,1)
     $  + HZ1(1) *HZ4(0,1,-1,1)
     $  + HZ5(0,1, -1,1,1)
     $  + 2.0000000000000000d+00*HZ5(0,1,1,-1,1)
     $  + 3.0000000000000000d+00*HZ5(0,1,1,1,-1)
      HZ5(1,1,1,-1,-1) = 
     $  + 8.3333333333333333d-02*HZ1(-1)*HZ1(-1)*HZ1(1)*HZ1(1)*HZ1(1)
     $  - 5.0000000000000000d-01*HZ1(-1)*HZ1(1)*HZ1(1)*HZ2(-1,1)
     $  + HZ1( -1)*HZ1(1)*HZ3(-1,1,1)
     $  - HZ1( -1)*HZ4(-1,1,1,1)
     $  + 5.0000000000000000d-01*HZ1(1)*HZ1(1)*HZ3(-1,-1,1)
     $  - HZ1(1) *HZ4(-1,-1,1,1)
     $  + HZ5( -1,-1,1,1,1)
      HZ5(1,1,1,-1,0) = 
     $  + 1.6666666666666666d-01*HZ1(-1)*HZ1(0)*HZ1(1)*HZ1(1)*HZ1(1)
     $  - 5.0000000000000000d-01*HZ1(0)*HZ1(1)*HZ1(1)*HZ2(-1,1)
     $  + HZ1(0) *HZ1(1)*HZ3(-1,1,1)
     $  - HZ1(0) *HZ4(-1,1,1,1)
     $  - 1.6666666666666666d-01*HZ1(1)*HZ1(1)*HZ1(1)*HZ2(0,-1)
     $  + 5.0000000000000000d-01*HZ1(1)*HZ1(1)*HZ3(0,-1,1)
     $  - HZ1(1) *HZ4(0,-1,1,1)
     $  + HZ5(0, -1,1,1,1)
      HZ5(1,1,1,-1,1) = 
     $  + 1.6666666666666666d-01*HZ1(1)*HZ1(1)*HZ1(1)*HZ2(-1,1)
     $  - HZ1(1) *HZ1(1)*HZ3(-1,1,1)
     $  + 3.0000000000000000d+00*HZ1(1)*HZ4(-1,1,1,1)
     $  - 4.0000000000000000d+00*HZ5(-1,1,1,1,1)
      HZ5(1,1,1,0,-1) = 
     $  - 5.0000000000000000d-01*HZ1(-1)*HZ1(1)*HZ1(1)*HZ2(0,1)
     $  + HZ1( -1)*HZ1(1)*HZ3(0,1,1)
     $  + 1.6666666666666666d-01*HZ1(1)*HZ1(1)*HZ1(1)*HZ2(0,-1)
     $  + 5.0000000000000000d-01*HZ1(1)*HZ1(1)*HZ3(-1,0,1)
     $  - HZ1(1) *HZ4(-1,0,1,1)
     $  - HZ5(0, -1,1,1,1)
     $  - HZ5(0,1, -1,1,1)
     $  - HZ5(0,1,1, -1,1)
     $  - HZ5(0,1,1,1, -1)
      HZ5(1,1,1,1,-1) = 
     $  + 4.1666666666666666d-02*HZ1(-1)*HZ1(1)*HZ1(1)*HZ1(1)*HZ1(1)
     $  - 1.6666666666666666d-01*HZ1(1)*HZ1(1)*HZ1(1)*HZ2(-1,1)
     $  + 5.0000000000000000d-01*HZ1(1)*HZ1(1)*HZ3(-1,1,1)
     $  - HZ1(1) *HZ4(-1,1,1,1)
     $  + HZ5( -1,1,1,1,1)
      endif

      return
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
      subroutine fillirr1dhplat0(y,nw,HY1,HY2,HY3,HY4,HY5,n1,n2) 
** evaluate the HPL from their power series expansions
** fillirr1dhplat0 is called by eval1dhplat0; 
** it is guaranteed that nw is in the range 1:4, and that (n1,n2) 
** take one of the pairs of values (0,1), (-1,0) or (-1,1) 
** 
** for y < 0 DOES NOT evaluates the immaginary part of H(0,y) = log(y) 
      implicit double precision (a-h,o-z) 
      dimension HY1(n1:n2),HY2(n1:n2,n1:n2),HY3(n1:n2,n1:n2,n1:n2), 
     $          HY4(n1:n2,n1:n2,n1:n2,n1:n2), 
     $          HY5(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2) 
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
        tu13 = 2d0*tu01*tu12 - tu11 
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
        tv13 = 2d0*tv01*tv12 - tv11
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
      if (nw.gt.4) then
      HY5(0,0,0,0,1) = 
     $  - 7.1883935399674661d-02
     $  + 5.6753782109711975d-01*tu01
     $  - 7.2212923207857999d-02*tu02
     $  + 5.8670869471466117d-03*tu03
     $  - 3.2928069579812269d-04*tu04
     $  + 1.2689437803900269d-05*tu05
     $  - 2.9273255639591246d-07*tu06
     $  + 1.0606320330117018d-09*tu07
     $  + 1.5490385851360920d-10*tu08
     $  - 2.2474559413025843d-12*tu09
     $  - 1.5435486840448295d-13*tu10
     $  + 4.3623305832897835d-15*tu11
     $  + 1.7529593826254750d-16*tu12
      HY5(0,0,0,1,1) = 
     $  + 1.0000521958818753d-02
     $  - 5.5972378484058781d-03*tu01
     $  + 1.0183616888857122d-02*tu02
     $  - 1.8839550865568384d-03*tu03
     $  + 1.8347215924120138d-04*tu04
     $  - 1.0930879320715947d-05*tu05
     $  + 3.7698877255509719d-07*tu06
     $  - 3.8459944472218428d-09*tu07
     $  - 2.4017837972601217d-10*tu08
     $  + 7.0263208823587761d-12*tu09
     $  + 2.5160915023993803d-13*tu10
     $  - 1.2783871754086183d-14*tu11
     $  - 2.8797679566416130d-16*tu12
     $  + 2.3134706555649888d-17*tu13
      HY5(0,0,1,0,1) = 
     $  + 1.9758848802070053d-02
     $  - 9.9630794104335436d-03*tu01
     $  + 2.0043658059117399d-02*tu02
     $  - 3.3447469330099638d-03*tu03
     $  + 2.8517465661196517d-04*tu04
     $  - 1.4231602601165430d-05*tu05
     $  + 3.6511844421403422d-07*tu06
     $  + 4.8298298252332493d-10*tu07
     $  - 2.8074265749907893d-10*tu08
     $  + 5.8456103354109943d-13*tu09
     $  + 3.7717268411948233d-13*tu10
     $  - 2.6564037404227836d-15*tu11
     $  - 5.7368628737502723d-16*tu12
      HY5(0,0,1,1,1) = 
     $  - 4.2151303147672207d-04
     $  + 2.4230392899354007d-03*tu01
     $  - 5.6276481104680303d-04*tu02
     $  + 8.2640657870678230d-04*tu03
     $  - 1.4169972575483287d-04*tu04
     $  + 1.1241768997510332d-05*tu05
     $  - 4.4757876435236780d-07*tu06
     $  + 4.0464701571442372d-09*tu07
     $  + 3.6693605986548572d-10*tu08
     $  - 8.1486736073711877d-12*tu09
     $  - 4.8362219356601919d-13*tu10
     $  + 1.5732397614617162d-14*tu11
     $  + 7.1637443194486341d-16*tu12
     $  - 2.9682648694764162d-17*tu13
      HY5(0,1,0,1,1) = 
     $  - 5.4016355841362305d-04
     $  + 3.5885207401063810d-03*tu01
     $  - 7.2088270891443272d-04*tu02
     $  + 1.2165742608270584d-03*tu03
     $  - 1.8111750748505087d-04*tu04
     $  + 1.2242347596998815d-05*tu05
     $  - 3.9805632127123734d-07*tu06
     $  + 1.3818410996670467d-09*tu07
     $  + 3.0027181656102247d-10*tu08
     $  - 2.5190497445473117d-12*tu09
     $  - 3.9056998701282666d-13*tu10
     $  + 4.6313279849969705d-15*tu11
     $  + 5.8249180278623573d-16*tu12
      HY5(0,1,1,1,1) = 
     $  + 3.6243571484950436d-04
     $  - 1.3106363932291666d-04*tu01
     $  + 4.8407754872309633d-04*tu02
     $  - 6.5531819661458333d-05*tu03
     $  + 1.2213913602645317d-04*tu04
     $  - 1.3106363932291666d-05*tu05
     $  + 4.9683501491940088d-07*tu06
     $  - 4.6646796570171982d-10*tu08
     $  + 6.6892469436720110d-13*tu10
     $  - 1.0496975350928139d-15*tu12
      endif 
** nw > 4 endif 
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
      if ( nw.gt.4 ) then 
      HY5(0,0,0,0,-1) = 
     $  + 7.1883935399674661d-02
     $  + 5.6753782109711975d-01*tv01
     $  + 7.2212923207857999d-02*tv02
     $  + 5.8670869471466117d-03*tv03
     $  + 3.2928069579812269d-04*tv04
     $  + 1.2689437803900269d-05*tv05
     $  + 2.9273255639591246d-07*tv06
     $  + 1.0606320330117018d-09*tv07
     $  - 1.5490385851360920d-10*tv08
     $  - 2.2474559413025843d-12*tv09
     $  + 1.5435486840448295d-13*tv10
     $  + 4.3623305832897835d-15*tv11
     $  - 1.7529593826254750d-16*tv12
      HY5(0,0,0,-1,-1) = 
     $  + 1.0000521958818753d-02
     $  + 5.5972378484058781d-03*tv01
     $  + 1.0183616888857122d-02*tv02
     $  + 1.8839550865568384d-03*tv03
     $  + 1.8347215924120138d-04*tv04
     $  + 1.0930879320715947d-05*tv05
     $  + 3.7698877255509719d-07*tv06
     $  + 3.8459944472218428d-09*tv07
     $  - 2.4017837972601217d-10*tv08
     $  - 7.0263208823587761d-12*tv09
     $  + 2.5160915023993803d-13*tv10
     $  + 1.2783871754086183d-14*tv11
     $  - 2.8797679566416130d-16*tv12
     $  - 2.3134706555649888d-17*tv13
      HY5(0,0,-1,0,-1) = 
     $  + 1.9758848802070053d-02
     $  + 9.9630794104335436d-03*tv01
     $  + 2.0043658059117399d-02*tv02
     $  + 3.3447469330099638d-03*tv03
     $  + 2.8517465661196517d-04*tv04
     $  + 1.4231602601165430d-05*tv05
     $  + 3.6511844421403422d-07*tv06
     $  - 4.8298298252332493d-10*tv07
     $  - 2.8074265749907893d-10*tv08
     $  - 5.8456103354109943d-13*tv09
     $  + 3.7717268411948233d-13*tv10
     $  + 2.6564037404227836d-15*tv11
     $  - 5.7368628737502723d-16*tv12
      HY5(0,0,-1,-1,-1) = 
     $  + 4.2151303147672207d-04
     $  + 2.4230392899354007d-03*tv01
     $  + 5.6276481104680303d-04*tv02
     $  + 8.2640657870678230d-04*tv03
     $  + 1.4169972575483287d-04*tv04
     $  + 1.1241768997510332d-05*tv05
     $  + 4.4757876435236780d-07*tv06
     $  + 4.0464701571442372d-09*tv07
     $  - 3.6693605986548572d-10*tv08
     $  - 8.1486736073711877d-12*tv09
     $  + 4.8362219356601919d-13*tv10
     $  + 1.5732397614617162d-14*tv11
     $  - 7.1637443194486341d-16*tv12
     $  - 2.9682648694764162d-17*tv13
      HY5(0,-1,0,-1,-1) = 
     $  + 5.4016355841362305d-04
     $  + 3.5885207401063810d-03*tv01
     $  + 7.2088270891443272d-04*tv02
     $  + 1.2165742608270584d-03*tv03
     $  + 1.8111750748505087d-04*tv04
     $  + 1.2242347596998815d-05*tv05
     $  + 3.9805632127123734d-07*tv06
     $  + 1.3818410996670467d-09*tv07
     $  - 3.0027181656102247d-10*tv08
     $  - 2.5190497445473117d-12*tv09
     $  + 3.9056998701282666d-13*tv10
     $  + 4.6313279849969705d-15*tv11
     $  - 5.8249180278623573d-16*tv12
      HY5(0,-1,-1,-1,-1) = 
     $  + 3.6243571484950436d-04
     $  + 1.3106363932291666d-04*tv01
     $  + 4.8407754872309633d-04*tv02
     $  + 6.5531819661458333d-05*tv03
     $  + 1.2213913602645317d-04*tv04
     $  + 1.3106363932291666d-05*tv05
     $  + 4.9683501491940088d-07*tv06
     $  - 4.6646796570171982d-10*tv08
     $  + 6.6892469436720110d-13*tv10
     $  - 1.0496975350928139d-15*tv12
      endif 
** nw > 4 endif 
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
      if ( nw.gt.4 ) then 
      HY5(-1,-1,-1,-1,1) = 
     $  - 3.8228565672333040d-02
     $  + 2.9434499414685087d-01*tu01
     $  - 3.8424518843646597d-02*tu02
     $  + 3.2580449034204883d-03*tu03
     $  - 1.9620389859663782d-04*tu04
     $  + 8.5297225882767261d-06*tu05
     $  - 2.5065346065715455d-07*tu06
     $  + 3.4415061242343063d-09*tu07
     $  + 7.3821799422249131d-11*tu08
     $  - 4.2538624555519857d-12*tu09
     $  - 7.2315211312081768d-16*tu10
     $  + 5.3108417435543907d-15*tu11
     $  - 9.9323734178161058d-17*tu12
     $  - 5.1747906167389938d-01*HY1(-1)
     $  - 2.6860659680402010d-01*HY1(-1)*HY1(-1)
     $  - 9.7040087744168750d-02*HY1(-1)*HY1(-1)*HY1(-1)
     $  - 2.8881132523331054d-02*HY1(-1)*HY1(-1)*HY1(-1)*HY1(-1)
      HY5(-1,-1,-1,1,1) = 
     $  + 2.7357462427685970d-02
     $  - 1.9031941959865697d-01*tu01
     $  + 2.7543492667729091d-02*tu02
     $  - 2.6605966861706822d-03*tu03
     $  + 1.8636146557007899d-04*tu04
     $  - 9.5345987092093788d-06*tu05
     $  + 3.3110554578258299d-07*tu06
     $  - 5.4595158175118823d-09*tu07
     $  - 1.1999359443095388d-10*tu08
     $  + 8.2902112001115834d-12*tu09
     $  - 1.2181518159922712d-14*tu10
     $  - 1.1713353633997384d-14*tu11
     $  + 2.3651675075246772d-16*tu12
     $  + 3.3160957134970784d-01*HY1(-1)
     $  + 1.5441268754841696d-01*HY1(-1)*HY1(-1)
     $  + 4.0037751159850118d-02*HY1(-1)*HY1(-1)*HY1(-1)
      HY5(-1,-1,1,-1,1) = 
     $  + 4.1458814971381822d-02
     $  - 3.0303763301003844d-01*tu01
     $  + 4.1691982869997716d-02*tu02
     $  - 3.7240373018258494d-03*tu03
     $  + 2.3344990525159924d-04*tu04
     $  - 1.0247007230046718d-05*tu05
     $  + 2.8185771897077629d-07*tu06
     $  - 2.2751673303889503d-09*tu07
     $  - 1.4878969225475077d-10*tu08
     $  + 4.1341440052910325d-12*tu09
     $  + 1.2696091792274665d-13*tu10
     $  - 7.5007742483794021d-15*tu11
     $  - 8.1903610410494478d-17*tu12
     $  + 5.3075770941318154d-01*HY1(-1)
     $  + 2.6126533022459244d-01*HY1(-1)*HY1(-1)
     $  - 6.9314718055994530d-01*HY1(-1)*HY3(-1,-1,1)
     $  + 4.8045301391820142d-01*HY3(-1,-1,1)
     $  + 2.0794415416798359d+00*HY4(-1,-1,-1,1)
      HY5(-1,-1,1,1,1) = 
     $  - 9.5924489184207016d-03
     $  + 5.2772883655818822d-02*tu01
     $  - 9.7366715772116844d-03*tu02
     $  + 1.4142029020952310d-03*tu03
     $  - 1.4463004412432978d-04*tu04
     $  + 9.8218990989160592d-06*tu05
     $  - 4.0715654417888680d-07*tu06
     $  + 6.7544293368802512d-09*tu07
     $  + 2.2866755156107150d-10*tu08
     $  - 1.2497892086833029d-11*tu09
     $  - 1.2173408174458276d-13*tu10
     $  + 2.1387606509664189d-14*tu11
     $  - 1.1689492941664195d-16*tu12
     $  - 3.4168851038812753d-17*tu13
     $  - 8.8326067366010207d-02*HY1(-1)
     $  - 2.7752054332410789d-02*HY1(-1)*HY1(-1)
      HY5(-1,1,-1,1,1) = 
     $  - 1.8886299866599669d-02
     $  + 1.1916393884386234d-01*tu01
     $  - 1.9076154105183607d-02*tu02
     $  + 2.2711434786163337d-03*tu03
     $  - 1.9021198135197550d-04*tu04
     $  + 1.0646966480063299d-05*tu05
     $  - 3.5753730456224019d-07*tu06
     $  + 3.8528188610744423d-09*tu07
     $  + 2.0528377803607047d-10*tu08
     $  - 6.6747997821119625d-12*tu09
     $  - 1.7956391666197610d-13*tu10
     $  + 1.1491009542808092d-14*tu11
     $  + 1.3285418003675663d-16*tu12
     $  - 2.0437039310996688d-01*HY1(-1)
     $  + 2.4022650695910071d-01*HY1(-1)*HY2(-1,1)
     $  - 1.6651232599446473d-01*HY2(-1,1)
     $  - 4.8045301391820142d-01*HY3(-1,-1,1)
      HY5(-1,1,1,1,1) = 
     $  + 1.8171632147206520d-03
     $  - 6.6063186750176729d-03*tu01
     $  + 1.9387809591623771d-03*tu02
     $  - 4.5816039839889645d-04*tu03
     $  + 1.2210356950937828d-04*tu04
     $  - 1.1632122874703951d-05*tu05
     $  + 4.8540847374707648d-07*tu06
     $  - 3.8582255619059790d-09*tu07
     $  - 4.1607939221717784d-10*tu08
     $  + 8.5808819391347779d-12*tu09
     $  + 5.1386964924538773d-13*tu10
     $  - 1.7412399897832375d-14*tu11
     $  - 6.4350879750488940d-16*tu12
     $  + 3.3295304169373632d-17*tu13
     $  + 9.6181291076284771d-03*HY1(-1)
      HY5(0,-1,-1,-1,1) = 
     $  - 3.9520972886442092d-02
     $  + 3.0546490933012253d-01*tu01
     $  - 3.9720322777238976d-02*tu02
     $  + 3.3467933739052856d-03*tu03
     $  - 1.9959549584313573d-04*tu04
     $  + 8.5499586955813381d-06*tu05
     $  - 2.4553310777499893d-07*tu06
     $  + 3.2119489955423060d-09*tu07
     $  + 7.1927096749843198d-11*tu08
     $  - 3.7520798598557195d-12*tu09
     $  - 1.1440674588088539d-14*tu10
     $  + 4.6108325614199602d-15*tu11
     $  - 6.0832641844878483d-17*tu12
     $  - 5.3721319360804020d-01*HY2(0,-1)
     $  - 5.8224052646501250d-01*HY3(0,-1,-1)
     $  - 6.9314718055994530d-01*HY4(0,-1,-1,-1)
      HY5(0,-1,-1,0,1) = 
     $  - 5.4931882964030792d-02
     $  + 4.2696670756984732d-01*tu01
     $  - 5.5201538838609055d-02*tu02
     $  + 4.6052217228062647d-03*tu03
     $  - 2.6995693356946050d-04*tu04
     $  + 1.1207715314455869d-05*tu05
     $  - 3.0093286122242289d-07*tu06
     $  + 2.9252018980885367d-09*tu07
     $  + 1.2606407378756556d-10*tu08
     $  - 4.4893551397064275d-12*tu09
     $  - 6.5916974412553446d-14*tu10
     $  + 6.6066633977237873d-15*tu11
     $  - 7.5128556447474642d-01*HY2(0,-1)
     $  - 8.2246703342411321d-01*HY3(0,-1,-1)
      HY5(0,-1,-1,1,-1) = 
     $  - 3.4235225044469286d-02
     $  - 2.3328399380319120d-01*tv01
     $  - 3.4504895043303540d-02*tv02
     $  - 3.7912859711035303d-03*tv03
     $  - 2.6991423023257518d-04*tv04
     $  - 1.1593278610997365d-05*tv05
     $  - 2.4410984754189600d-07*tv06
     $  + 9.0404989403384568d-10*tv07
     $  + 1.2141001239198240d-10*tv08
     $  - 2.0560625998471416d-12*tv09
     $  - 1.4056423603087806d-13*tv10
     $  + 3.5816535397430763d-15*tv11
     $  + 2.0254744342336250d-16*tv12
     $  + 4.0357837932696163d-01*HY2(0,-1)
     $  + 4.8045301391820142d-01*HY3(0,-1,-1)
     $  + 6.9314718055994530d-01*HY4(0,-1,-1,1)
      HY5(0,-1,-1,1,1) = 
     $  + 2.6517569052959506d-02
     $  - 1.7779971519631702d-01*tu01
     $  + 2.6707704581331518d-02*tu02
     $  - 2.6648397927870009d-03*tu03
     $  + 1.9046558373551644d-04*tu04
     $  - 9.7589088614312887d-06*tu05
     $  + 3.2991992381966999d-07*tu06
     $  - 4.8561819501772588d-09*tu07
     $  - 1.3539320891044926d-10*tu08
     $  + 7.2902837287752652d-12*tu09
     $  + 4.6565266727231378d-14*tu10
     $  - 1.0737158384822919d-14*tu11
     $  + 8.8429222559013066d-17*tu12
     $  + 3.0882537509683393d-01*HY2(0,-1)
     $  + 2.4022650695910071d-01*HY3(0,-1,-1)
      HY5(0,-1,0,-1,1) = 
     $  - 4.6155088778537896d-02
     $  + 3.5912469692896398d-01*tu01
     $  - 4.6380859407450647d-02*tu02
     $  + 3.8633681101991486d-03*tu03
     $  - 2.2602267773980855d-04*tu04
     $  + 9.3668504249552797d-06*tu05
     $  - 2.5195369443104781d-07*tu06
     $  + 2.5727032623578792d-09*tu07
     $  + 9.5079191917240018d-11*tu08
     $  - 3.3913895763562877d-12*tu09
     $  - 5.3430708391974222d-14*tu10
     $  + 4.8668530674676388d-15*tu11
     $  - 6.3196619783816790d-01*HY2(0,-1)
     $  - 3.4657359027997265d-01*HY2(0,-1)*HY2(0,-1)
     $  + 1.3862943611198906d+00*HY4(0,0,-1,-1)
      HY5(0,-1,0,1,-1) = 
     $  - 4.8531326139900090d-02
     $  - 3.2657418615138227d-01*tv01
     $  - 4.8786752720619143d-02*tv02
     $  - 4.3579589418437103d-03*tv03
     $  - 2.5565369900267449d-04*tv04
     $  - 9.9174468231600028d-06*tv05
     $  - 2.2703683570409072d-07*tv06
     $  - 1.2319581804863699d-09*tv07
     $  + 8.1382125029078829d-11*tv08
     $  + 1.0513084104970317d-12*tv09
     $  - 6.5723688498545463d-14*tv10
     $  - 1.2162924878946765d-15*tv11
     $  + 6.9118041556917233d-17*tv12
     $  + 5.7009070532142637d-01*HY2(0,-1)
     $  + 6.9314718055994530d-01*HY4(0,-1,0,1)
      HY5(0,-1,0,1,1) = 
     $  + 1.7588255219497996d-02
     $  - 8.9173088210648816d-02*tu01
     $  + 1.7770012282463661d-02*tu02
     $  - 2.1943047496140125d-03*tu03
     $  + 1.8211824241696307d-04*tu04
     $  - 1.0254443539590487d-05*tu05
     $  + 3.6099159031901396d-07*tu06
     $  - 4.8590233971583916d-09*tu07
     $  - 1.8773473287955894d-10*tu08
     $  + 8.2188353537397648d-12*tu09
     $  + 1.2622338876212118d-13*tu10
     $  - 1.3549277713840031d-14*tu11
     $  - 2.2995013881263894d-17*tu12
     $  + 2.1790299145411898d-17*tu13
     $  + 1.5025711289494928d-01*HY2(0,-1)
      HY5(0,-1,1,-1,-1) = 
     $  - 1.5970421220241350d-02
     $  - 9.8098905332736657d-02*tv01
     $  - 1.6167874610049324d-02*tv02
     $  - 2.1909422567422571d-03*tv03
     $  - 1.9780330842002595d-04*tv04
     $  - 1.1143616753061323d-05*tv05
     $  - 3.4970788285258598d-07*tv06
     $  - 2.4737298345971262d-09*tv07
     $  + 2.1049905271699688d-10*tv08
     $  + 3.6829751452188057d-12*tv09
     $  - 2.2984898297836793d-13*tv10
     $  - 5.8041435355729247d-15*tv11
     $  + 2.9753327209204345d-16*tv12
     $  + 1.6651232599446473d-01*HY2(0,-1)
     $  + 2.4022650695910071d-01*HY3(0,-1,1)
      HY5(0,-1,1,-1,1) = 
     $  + 4.1824168833514614d-02
     $  - 2.9878381943835587d-01*tu01
     $  + 4.2063530444775636d-02*tu02
     $  - 3.8145573036522436d-03*tu03
     $  + 2.3963305636044003d-04*tu04
     $  - 1.0345567145401369d-05*tu05
     $  + 2.7130517497620301d-07*tu06
     $  - 1.7273255620176561d-09*tu07
     $  - 1.3978876407646153d-10*tu08
     $  + 2.7492226699257445d-12*tu09
     $  + 1.3553332679630927d-13*tu10
     $  - 4.7656087734784391d-15*tu11
     $  - 1.4423054189885051d-16*tu12
     $  + 5.2253066044918489d-01*HY2(0,-1)
     $  + 4.8045301391820142d-01*HY3(0,-1,1)
     $  - 6.9314718055994530d-01*HY4(0,-1,1,-1)
      HY5(0,-1,1,0,1) = 
     $  + 3.1399803193734681d-02
     $  - 1.5922070722448512d-01*tu01
     $  + 3.1673651121964784d-02*tu02
     $  - 3.6733024485531673d-03*tu03
     $  + 2.7421207386188742d-04*tu04
     $  - 1.3254697744560554d-05*tu05
     $  + 3.6389461745542725d-07*tu06
     $  - 1.3701468350861509d-09*tu07
     $  - 2.5073096543709732d-10*tu08
     $  + 3.7881257959182565d-12*tu09
     $  + 2.8301894299500673d-13*tu10
     $  - 7.9052198960553283d-15*tu11
     $  - 3.4475230007347254d-16*tu12
     $  + 2.6957647953152780d-01*HY2(0,-1)
      HY5(0,-1,1,1,-1) = 
     $  + 4.4239483860061491d-02
     $  + 3.2493754537326350d-01*tv01
     $  + 4.4462186515298750d-02*tv02
     $  + 3.8108642536456408d-03*tv03
     $  + 2.2293231798666157d-04*tv04
     $  + 8.9901447113222695d-06*tv05
     $  + 2.2959225139995976d-07*tv06
     $  + 2.1985163564934731d-09*tv07
     $  - 7.0460025669936385d-11*tv08
     $  - 2.1618428565949139d-12*tv09
     $  + 3.7961815108691857d-14*tv10
     $  + 2.5062265861081495d-15*tv11
     $  - 5.7009070532142637d-01*HY2(0,-1)
     $  - 5.8224052646501250d-01*HY3(0,-1,1)
     $  + 6.9314718055994530d-01*HY4(0,-1,1,1)
      HY5(0,-1,1,1,1) = 
     $  - 6.7365547952682237d-03
     $  + 3.4286529141621258d-02*tu01
     $  - 6.8840170170275732d-03*tu02
     $  + 1.2704332689803235d-03*tu03
     $  - 1.4787767126158972d-04*tu04
     $  + 1.0413585642621615d-05*tu05
     $  - 4.1517140759699248d-07*tu06
     $  + 5.3439659402510695d-09*tu07
     $  + 2.7781700180675136d-10*tu08
     $  - 9.9205586909926692d-12*tu09
     $  - 2.7736473454099143d-13*tu10
     $  + 1.7683378840938849d-14*tu11
     $  + 2.7645451246534933d-16*tu12
     $  - 3.0721858829145982d-17*tu13
     $  - 5.5504108664821579d-02*HY2(0,-1)
      HY5(0,0,-1,-1,1) = 
     $  - 4.2561833953883235d-02
     $  + 3.3089514249958932d-01*tu01
     $  - 4.2771142986440044d-02*tu02
     $  + 3.5688685742597137d-03*tu03
     $  - 2.0954973101068707d-04*tu04
     $  + 8.7544434158133471d-06*tu05
     $  - 2.4062186086108688d-07*tu06
     $  + 2.7728194984089187d-09*tu07
     $  + 7.6559566355403117d-11*tu08
     $  - 3.1501612342558870d-12*tu09
     $  - 3.3459877468585252d-14*tu10
     $  + 4.0333595766267801d-15*tu11
     $  - 5.8224052646501250d-01*HY3(0,0,-1)
     $  - 6.9314718055994530d-01*HY4(0,0,-1,-1)
      HY5(0,0,-1,0,1) = 
     $  - 5.9694681535720023d-02
     $  + 4.6714213068105458d-01*tu01
     $  - 5.9979199873202721d-02*tu02
     $  + 4.9475467229010258d-03*tu03
     $  - 2.8481020227431533d-04*tu04
     $  + 1.1478749310592714d-05*tu05
     $  - 2.9173546427111336d-07*tu06
     $  + 2.2634918139358350d-09*tu07
     $  + 1.2923277934839262d-10*tu08
     $  - 3.4600488918561596d-12*tu09
     $  - 9.4505709251699823d-14*tu10
     $  + 5.3596484866903464d-15*tu11
     $  + 6.0838246089431129d-17*tu12
     $  - 8.2246703342411321d-01*HY3(0,0,-1)
      HY5(0,0,-1,1,-1) = 
     $  - 4.1110148228889138d-02
     $  - 2.7597446804002515d-01*tv01
     $  - 4.1362136148261452d-02*tv02
     $  - 3.9264219532611858d-03*tv03
     $  - 2.5225247282294274d-04*tv04
     $  - 1.0792153931365190d-05*tv05
     $  - 2.6440687189609348d-07*tv06
     $  - 7.5507087052363252d-10*tv07
     $  + 1.4640592504212342d-10*tv08
     $  + 1.0410408318930554d-12*tv09
     $  - 1.7256282627937244d-13*tv10
     $  - 2.0115212902641882d-15*tv11
     $  + 2.4492838656439672d-16*tv12
     $  + 4.8045301391820142d-01*HY3(0,0,-1)
     $  + 6.9314718055994530d-01*HY4(0,0,-1,1)
      HY5(0,0,-1,1,1) = 
     $  + 2.3290431799544722d-02
     $  - 1.3976916667674197d-01*tu01
     $  + 2.3484837128475192d-02*tu02
     $  - 2.5652450331604480d-03*tu03
     $  + 1.9474432423788960d-04*tu04
     $  - 1.0235257714295329d-05*tu05
     $  + 3.3882452171323895d-07*tu06
     $  - 4.1398941829422872d-09*tu07
     $  - 1.7065390810468420d-10*tu08
     $  + 6.5500661689981377d-12*tu09
     $  + 1.3171293231469574d-13*tu10
     $  - 1.0601060206261924d-14*tu11
     $  - 8.5654606012582514d-17*tu12
     $  + 2.4022650695910071d-01*HY3(0,0,-1)
      HY5(0,0,0,-1,1) = 
     $  - 5.0232086741389027d-02
     $  + 3.9364648125763205d-01*tu01
     $  - 5.0470343357276722d-02*tu02
     $  + 4.1544817848152132d-03*tu03
     $  - 2.3850059599743822d-04*tu04
     $  + 9.5855250288941914d-06*tu05
     $  - 2.4388351988986297d-07*tu06
     $  + 2.0282340296972451d-09*tu07
     $  + 9.6515154204163078d-11*tu08
     $  - 2.5201744413354707d-12*tu09
     $  - 7.4636488424158657d-14*tu10
     $  + 3.7370382545755637d-15*tu11
     $  + 6.2133891701574733d-17*tu12
     $  - 6.9314718055994530d-01*HY4(0,0,0,-1)
      HY5(0,0,0,1,-1) = 
     $  - 5.0232086741389027d-02
     $  - 3.9364648125763205d-01*tv01
     $  - 5.0470343357276722d-02*tv02
     $  - 4.1544817848152132d-03*tv03
     $  - 2.3850059599743822d-04*tv04
     $  - 9.5855250288941914d-06*tv05
     $  - 2.4388351988986297d-07*tv06
     $  - 2.0282340296972451d-09*tv07
     $  + 9.6515154204163078d-11*tv08
     $  + 2.5201744413354707d-12*tv09
     $  - 7.4636488424158657d-14*tv10
     $  - 3.7370382545755637d-15*tv11
     $  + 6.2133891701574733d-17*tv12
     $  + 6.9314718055994530d-01*HY4(0,0,0,1)
      HY5(0,0,1,-1,-1) = 
     $  - 2.3290431799544722d-02
     $  - 1.3976916667674197d-01*tv01
     $  - 2.3484837128475192d-02*tv02
     $  - 2.5652450331604480d-03*tv03
     $  - 1.9474432423788960d-04*tv04
     $  - 1.0235257714295329d-05*tv05
     $  - 3.3882452171323895d-07*tv06
     $  - 4.1398941829422872d-09*tv07
     $  + 1.7065390810468420d-10*tv08
     $  + 6.5500661689981377d-12*tv09
     $  - 1.3171293231469574d-13*tv10
     $  - 1.0601060206261924d-14*tv11
     $  + 8.5654606012582514d-17*tv12
     $  + 2.4022650695910071d-01*HY3(0,0,1)
      HY5(0,0,1,-1,1) = 
     $  + 4.1110148228889138d-02
     $  - 2.7597446804002515d-01*tu01
     $  + 4.1362136148261452d-02*tu02
     $  - 3.9264219532611858d-03*tu03
     $  + 2.5225247282294274d-04*tu04
     $  - 1.0792153931365190d-05*tu05
     $  + 2.6440687189609348d-07*tu06
     $  - 7.5507087052363252d-10*tu07
     $  - 1.4640592504212342d-10*tu08
     $  + 1.0410408318930554d-12*tu09
     $  + 1.7256282627937244d-13*tu10
     $  - 2.0115212902641882d-15*tu11
     $  - 2.4492838656439672d-16*tu12
     $  + 4.8045301391820142d-01*HY3(0,0,1)
     $  - 6.9314718055994530d-01*HY4(0,0,1,-1)
      HY5(0,0,1,0,-1) = 
     $  - 5.9694681535720023d-02
     $  - 4.6714213068105458d-01*tv01
     $  - 5.9979199873202721d-02*tv02
     $  - 4.9475467229010258d-03*tv03
     $  - 2.8481020227431533d-04*tv04
     $  - 1.1478749310592714d-05*tv05
     $  - 2.9173546427111336d-07*tv06
     $  - 2.2634918139358350d-09*tv07
     $  + 1.2923277934839262d-10*tv08
     $  + 3.4600488918561596d-12*tv09
     $  - 9.4505709251699823d-14*tv10
     $  - 5.3596484866903464d-15*tv11
     $  + 6.0838246089431129d-17*tv12
     $  + 8.2246703342411321d-01*HY3(0,0,1)
      HY5(0,0,1,1,-1) = 
     $  + 4.2561833953883235d-02
     $  + 3.3089514249958932d-01*tv01
     $  + 4.2771142986440044d-02*tv02
     $  + 3.5688685742597137d-03*tv03
     $  + 2.0954973101068707d-04*tv04
     $  + 8.7544434158133471d-06*tv05
     $  + 2.4062186086108688d-07*tv06
     $  + 2.7728194984089187d-09*tv07
     $  - 7.6559566355403117d-11*tv08
     $  - 3.1501612342558870d-12*tv09
     $  + 3.3459877468585252d-14*tv10
     $  + 4.0333595766267801d-15*tv11
     $  - 5.8224052646501250d-01*HY3(0,0,1)
     $  + 6.9314718055994530d-01*HY4(0,0,1,1)
      HY5(0,1,-1,-1,-1) = 
     $  - 6.7365547952682237d-03
     $  - 3.4286529141621258d-02*tv01
     $  - 6.8840170170275732d-03*tv02
     $  - 1.2704332689803235d-03*tv03
     $  - 1.4787767126158972d-04*tv04
     $  - 1.0413585642621615d-05*tv05
     $  - 4.1517140759699248d-07*tv06
     $  - 5.3439659402510695d-09*tv07
     $  + 2.7781700180675136d-10*tv08
     $  + 9.9205586909926692d-12*tv09
     $  - 2.7736473454099143d-13*tv10
     $  - 1.7683378840938849d-14*tv11
     $  + 2.7645451246534933d-16*tv12
     $  + 3.0721858829145982d-17*tv13
     $  + 5.5504108664821579d-02*HY2(0,1)
      HY5(0,1,-1,-1,1) = 
     $  + 4.4239483860061491d-02
     $  - 3.2493754537326350d-01*tu01
     $  + 4.4462186515298750d-02*tu02
     $  - 3.8108642536456408d-03*tu03
     $  + 2.2293231798666157d-04*tu04
     $  - 8.9901447113222695d-06*tu05
     $  + 2.2959225139995976d-07*tu06
     $  - 2.1985163564934731d-09*tu07
     $  - 7.0460025669936385d-11*tu08
     $  + 2.1618428565949139d-12*tu09
     $  + 3.7961815108691857d-14*tu10
     $  - 2.5062265861081495d-15*tu11
     $  + 5.7009070532142637d-01*HY2(0,1)
     $  - 5.8224052646501250d-01*HY3(0,1,-1)
     $  - 6.9314718055994530d-01*HY4(0,1,-1,-1)
      HY5(0,1,-1,1,-1) = 
     $  + 4.1824168833514614d-02
     $  + 2.9878381943835587d-01*tv01
     $  + 4.2063530444775636d-02*tv02
     $  + 3.8145573036522436d-03*tv03
     $  + 2.3963305636044003d-04*tv04
     $  + 1.0345567145401369d-05*tv05
     $  + 2.7130517497620301d-07*tv06
     $  + 1.7273255620176561d-09*tv07
     $  - 1.3978876407646153d-10*tv08
     $  - 2.7492226699257445d-12*tv09
     $  + 1.3553332679630927d-13*tv10
     $  + 4.7656087734784391d-15*tv11
     $  - 1.4423054189885051d-16*tv12
     $  - 5.2253066044918489d-01*HY2(0,1)
     $  + 4.8045301391820142d-01*HY3(0,1,-1)
     $  + 6.9314718055994530d-01*HY4(0,1,-1,1)
      HY5(0,1,-1,1,1) = 
     $  - 1.5970421220241350d-02
     $  + 9.8098905332736657d-02*tu01
     $  - 1.6167874610049324d-02*tu02
     $  + 2.1909422567422571d-03*tu03
     $  - 1.9780330842002595d-04*tu04
     $  + 1.1143616753061323d-05*tu05
     $  - 3.4970788285258598d-07*tu06
     $  + 2.4737298345971262d-09*tu07
     $  + 2.1049905271699688d-10*tu08
     $  - 3.6829751452188057d-12*tu09
     $  - 2.2984898297836793d-13*tu10
     $  + 5.8041435355729247d-15*tu11
     $  + 2.9753327209204345d-16*tu12
     $  - 1.6651232599446473d-01*HY2(0,1)
     $  + 2.4022650695910071d-01*HY3(0,1,-1)
      HY5(0,1,0,1,-1) = 
     $  + 4.6155088778537896d-02
     $  + 3.5912469692896398d-01*tv01
     $  + 4.6380859407450647d-02*tv02
     $  + 3.8633681101991486d-03*tv03
     $  + 2.2602267773980855d-04*tv04
     $  + 9.3668504249552797d-06*tv05
     $  + 2.5195369443104781d-07*tv06
     $  + 2.5727032623578792d-09*tv07
     $  - 9.5079191917240018d-11*tv08
     $  - 3.3913895763562877d-12*tv09
     $  + 5.3430708391974222d-14*tv10
     $  + 4.8668530674676388d-15*tv11
     $  - 6.3196619783816790d-01*HY2(0,1)
     $  + 3.4657359027997265d-01*HY2(0,1)*HY2(0,1)
     $  - 1.3862943611198906d+00*HY4(0,0,1,1)
      HY5(0,1,1,-1,-1) = 
     $  + 2.6517569052959506d-02
     $  + 1.7779971519631702d-01*tv01
     $  + 2.6707704581331518d-02*tv02
     $  + 2.6648397927870009d-03*tv03
     $  + 1.9046558373551644d-04*tv04
     $  + 9.7589088614312887d-06*tv05
     $  + 3.2991992381966999d-07*tv06
     $  + 4.8561819501772588d-09*tv07
     $  - 1.3539320891044926d-10*tv08
     $  - 7.2902837287752652d-12*tv09
     $  + 4.6565266727231378d-14*tv10
     $  + 1.0737158384822919d-14*tv11
     $  + 8.8429222559013066d-17*tv12
     $  - 3.0882537509683393d-01*HY2(0,1)
     $  + 2.4022650695910071d-01*HY3(0,1,1)
      HY5(0,1,1,-1,1) = 
     $  - 3.4235225044469286d-02
     $  + 2.3328399380319120d-01*tu01
     $  - 3.4504895043303540d-02*tu02
     $  + 3.7912859711035303d-03*tu03
     $  - 2.6991423023257518d-04*tu04
     $  + 1.1593278610997365d-05*tu05
     $  - 2.4410984754189600d-07*tu06
     $  - 9.0404989403384568d-10*tu07
     $  + 1.2141001239198240d-10*tu08
     $  + 2.0560625998471416d-12*tu09
     $  - 1.4056423603087806d-13*tu10
     $  - 3.5816535397430763d-15*tu11
     $  + 2.0254744342336250d-16*tu12
     $  - 4.0357837932696163d-01*HY2(0,1)
     $  + 4.8045301391820142d-01*HY3(0,1,1)
     $  - 6.9314718055994530d-01*HY4(0,1,1,-1)
      HY5(0,1,1,1,-1) = 
     $  - 3.9520972886442092d-02
     $  - 3.0546490933012253d-01*tv01
     $  - 3.9720322777238976d-02*tv02
     $  - 3.3467933739052856d-03*tv03
     $  - 1.9959549584313573d-04*tv04
     $  - 8.5499586955813381d-06*tv05
     $  - 2.4553310777499893d-07*tv06
     $  - 3.2119489955423060d-09*tv07
     $  + 7.1927096749843198d-11*tv08
     $  + 3.7520798598557195d-12*tv09
     $  - 1.1440674588088539d-14*tv10
     $  - 4.6108325614199602d-15*tv11
     $  - 6.0832641844878483d-17*tv12
     $  + 5.3721319360804020d-01*HY2(0,1)
     $  - 5.8224052646501250d-01*HY3(0,1,1)
     $  + 6.9314718055994530d-01*HY4(0,1,1,1)
      endif 
** nw > 4 endif 
      endif 
** (n1,n2) = (-1,1) -- completion endif 
      return 
      end 
************************************************************************ 
      subroutine fillirr1dhplat1(r,nw,HR1,HR2,HR3,HR4,HR5, 
     $                                HY1,HY2,HY3,HY4,HY5, 
     $                                Hi1,Hi2,Hi3,Hi4,Hi5,n1,n2) 
** evaluates the HPL for r2m1 < y < r2p1
** fillirr1dhplat1 is called by eval1dhplat1 after calling 
** fillirr1dhplat0 with argument r=(1-y)/(1+y) 
** it is guaranteed that nw is in the range 2:4, and that (n1,n2) 
** take one of the pairs of values (0,1), (-1,0) or (-1,1) 
      implicit double precision (a-h,o-z) 
      dimension HR1(-1:1),HR2(-1:1,-1:1),HR3(-1:1,-1:1,-1:1), 
     $          HR4(-1:1,-1:1,-1:1,-1:1),
     $          HR5(-1:1,-1:1,-1:1,-1:1,-1:1) 
      dimension HY1(n1:n2),HY2(n1:n2,n1:n2),HY3(n1:n2,n1:n2,n1:n2), 
     $          HY4(n1:n2,n1:n2,n1:n2,n1:n2), 
     $          HY5(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2), 
     $          Hi5(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2) 
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
      if ( nw.gt.4 ) then 
      HY5(0,0,0,0,1) = 
     $  + 1.0369277551433699d+00
     $  - 1.0823232337111381d+00*HR1(-1)
     $  + 6.0102845157979714d-01*HR1(-1)*HR1(-1)
     $  - 2.7415567780803773d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 2.8881132523331054d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 8.3333333333333333d-03*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 8.3333333333333333d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR2(0,1)
     $  - 2.5000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)*HR1(1)
     $  - 8.2246703342411321d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - 1.7328679513998632d-01*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  + 8.3333333333333333d-02*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR2(0,-1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR2(0,1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(0,-1,-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(0,-1,1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(0,1,-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(0,1,1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(0)*HR1(1)*HR1(1)*HR1(1)
     $  + 1.2020569031595942d+00*HR1(-1)*HR1(1)
     $  - 8.2246703342411321d-01*HR1(-1)*HR1(1)*HR1(1)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(1)*HR1(1)*HR1(1)
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(1)*HR1(1)*HR1(1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(1)*HR1(1)*HR2(0,-1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(1)*HR1(1)*HR2(0,1)
     $  - HR1( -1)*HR1(1)*HR3(0,-1,-1)
     $  - HR1( -1)*HR1(1)*HR3(0,-1,1)
     $  - HR1( -1)*HR1(1)*HR3(0,1,-1)
     $  - HR1( -1)*HR1(1)*HR3(0,1,1)
     $  + HR1( -1)*HR4(0,-1,-1,-1)
     $  + HR1( -1)*HR4(0,-1,-1,1)
     $  + HR1( -1)*HR4(0,-1,1,-1)
     $  + HR1( -1)*HR4(0,-1,1,1)
     $  + HR1( -1)*HR4(0,1,-1,-1)
     $  + HR1( -1)*HR4(0,1,-1,1)
     $  + HR1( -1)*HR4(0,1,1,-1)
     $  + HR1( -1)*HR4(0,1,1,1)
     $  - 4.1666666666666666d-02*HR1(0)*HR1(1)*HR1(1)*HR1(1)*HR1(1)
     $  - 1.0823232337111381d+00*HR1(1)
     $  + 6.0102845157979714d-01*HR1(1)*HR1(1)
     $  - 2.7415567780803773d-01*HR1(1)*HR1(1)*HR1(1)
     $  - 2.8881132523331054d-02*HR1(1)*HR1(1)*HR1(1)*HR1(1)
     $  - 1.6666666666666666d-01*HR1(1)*HR1(1)*HR1(1)*HR2(-1,1)
     $  + 1.6666666666666666d-01*HR1(1)*HR1(1)*HR1(1)*HR2(0,-1)
     $  + 1.6666666666666666d-01*HR1(1)*HR1(1)*HR1(1)*HR2(0,1)
     $  - 5.0000000000000000d-01*HR1(1)*HR1(1)*HR3(-1,-1,1)
     $  + 5.0000000000000000d-01*HR1(1)*HR1(1)*HR3(-1,1,1)
     $  - 5.0000000000000000d-01*HR1(1)*HR1(1)*HR3(0,-1,-1)
     $  - 5.0000000000000000d-01*HR1(1)*HR1(1)*HR3(0,-1,1)
     $  - 5.0000000000000000d-01*HR1(1)*HR1(1)*HR3(0,1,-1)
     $  - 5.0000000000000000d-01*HR1(1)*HR1(1)*HR3(0,1,1)
     $  - HR1(1) *HR4(-1,-1,-1,1)
     $  + HR1(1) *HR4(-1,-1,1,1)
     $  - HR1(1) *HR4(-1,1,1,1)
     $  + HR1(1) *HR4(0,-1,-1,-1)
     $  + HR1(1) *HR4(0,-1,-1,1)
     $  + HR1(1) *HR4(0,-1,1,-1)
     $  + HR1(1) *HR4(0,-1,1,1)
     $  + HR1(1) *HR4(0,1,-1,-1)
     $  + HR1(1) *HR4(0,1,-1,1)
     $  + HR1(1) *HR4(0,1,1,-1)
     $  + HR1(1) *HR4(0,1,1,1)
     $  - HR5( -1,-1,-1,-1,1)
     $  + HR5( -1,-1,-1,1,1)
     $  - HR5( -1,-1,1,1,1)
     $  + HR5( -1,1,1,1,1)
     $  - HR5(0, -1,-1,-1,-1)
     $  - HR5(0, -1,-1,-1,1)
     $  - HR5(0, -1,-1,1,-1)
     $  - HR5(0, -1,-1,1,1)
     $  - HR5(0, -1,1,-1,-1)
     $  - HR5(0, -1,1,-1,1)
     $  - HR5(0, -1,1,1,-1)
     $  - HR5(0, -1,1,1,1)
     $  - HR5(0,1, -1,-1,-1)
     $  - HR5(0,1, -1,-1,1)
     $  - HR5(0,1, -1,1,-1)
     $  - HR5(0,1, -1,1,1)
     $  - HR5(0,1,1, -1,-1)
     $  - HR5(0,1,1, -1,1)
     $  - HR5(0,1,1,1, -1)
     $  - HR5(0,1,1,1,1) 
      HY5(0,0,0,1,1) = 
     $  + 9.6551159989443734d-02
     $  - 2.7058080842778454d-01*HR1(-1)
     $  + 6.0102845157979714d-01*HR1(-1)*HR1(-1)
     $  - 4.0037751159850118d-02*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 2.8881132523331054d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 8.3333333333333333d-03*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  - 8.3333333333333333d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)*HR1(0)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 8.3333333333333333d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  - 2.5000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(0)*HR1(1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  + 2.5000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR2(0,-1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR2(0,1)
     $  - 1.2011325347955035d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  + 1.7328679513998632d-01*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  - 8.3333333333333333d-02*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)*HR1(1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR2(0,1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(0,-1,-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(0,0,-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(0,0,1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(0,1,-1)
     $  - 2.5000000000000000d-01*HR1(-1)*HR1(0)*HR1(0)*HR1(1)*HR1(1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(0)*HR1(1)*HR1(1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(0)*HR1(1)*HR1(1)*HR1(1)
     $  + HR1( -1)*HR1(0)*HR1(1)*HR2(0,-1)
     $  + HR1( -1)*HR1(0)*HR1(1)*HR2(0,1)
     $  - HR1( -1)*HR1(0)*HR3(0,-1,-1)
     $  - HR1( -1)*HR1(0)*HR3(0,-1,1)
     $  - HR1( -1)*HR1(0)*HR3(0,1,-1)
     $  - HR1( -1)*HR1(0)*HR3(0,1,1)
     $  + 1.2020569031595942d+00*HR1(-1)*HR1(1)
     $  - 1.2011325347955035d-01*HR1(-1)*HR1(1)*HR1(1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(1)*HR1(1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(1)*HR1(1)*HR2(-1,1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(1)*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(1)*HR2(0,1)
     $  + HR1( -1)*HR1(1)*HR3(-1,-1,1)
     $  - HR1( -1)*HR1(1)*HR3(-1,1,1)
     $  - HR1( -1)*HR1(1)*HR3(0,-1,-1)
     $  - HR1( -1)*HR1(1)*HR3(0,0,-1)
     $  - HR1( -1)*HR1(1)*HR3(0,0,1)
     $  - HR1( -1)*HR1(1)*HR3(0,1,-1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR3(0,-1,-1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR3(0,-1,1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR3(0,1,-1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR3(0,1,1)
     $  + HR1( -1)*HR4(-1,-1,-1,1)
     $  - HR1( -1)*HR4(-1,-1,1,1)
     $  + HR1( -1)*HR4(-1,1,1,1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR4(0,-1,-1,-1)
     $  + HR1( -1)*HR4(0,-1,-1,1)
     $  + HR1( -1)*HR4(0,-1,1,-1)
     $  + HR1( -1)*HR4(0,0,-1,-1)
     $  + HR1( -1)*HR4(0,0,-1,1)
     $  + HR1( -1)*HR4(0,0,1,-1)
     $  + HR1( -1)*HR4(0,0,1,1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR4(0,1,-1,-1)
     $  + HR1( -1)*HR4(0,1,-1,1)
     $  + HR1( -1)*HR4(0,1,1,-1)
     $  - 8.3333333333333333d-02*HR1(0)*HR1(0)*HR1(1)*HR1(1)*HR1(1)
     $  - 1.1552453009332421d-01*HR1(0)*HR1(1)*HR1(1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(0)*HR1(1)*HR1(1)*HR2(-1,1)
     $  + 5.0000000000000000d-01*HR1(0)*HR1(1)*HR1(1)*HR2(0,-1)
     $  + 5.0000000000000000d-01*HR1(0)*HR1(1)*HR1(1)*HR2(0,1)
     $  - HR1(0) *HR1(1)*HR3(-1,-1,1)
     $  + HR1(0) *HR1(1)*HR3(-1,1,1)
     $  - HR1(0) *HR1(1)*HR3(0,-1,-1)
     $  - HR1(0) *HR1(1)*HR3(0,-1,1)
     $  - HR1(0) *HR1(1)*HR3(0,1,-1)
     $  - HR1(0) *HR1(1)*HR3(0,1,1)
     $  - HR1(0) *HR4(-1,-1,-1,1)
     $  + HR1(0) *HR4(-1,-1,1,1)
     $  - HR1(0) *HR4(-1,1,1,1)
     $  + HR1(0) *HR4(0,-1,-1,-1)
     $  + HR1(0) *HR4(0,-1,-1,1)
     $  + HR1(0) *HR4(0,-1,1,-1)
     $  + HR1(0) *HR4(0,-1,1,1)
     $  + HR1(0) *HR4(0,1,-1,-1)
     $  + HR1(0) *HR4(0,1,-1,1)
     $  + HR1(0) *HR4(0,1,1,-1)
     $  + HR1(0) *HR4(0,1,1,1)
     $  - 2.7058080842778454d-01*HR1(1)
     $  + 6.0102845157979714d-01*HR1(1)*HR1(1)
     $  - 4.0037751159850118d-02*HR1(1)*HR1(1)*HR1(1)
     $  - 3.4657359027997265d-01*HR1(1)*HR1(1)*HR2(-1,1)
     $  + 3.4657359027997265d-01*HR1(1)*HR1(1)*HR2(0,-1)
     $  + 3.4657359027997265d-01*HR1(1)*HR1(1)*HR2(0,1)
     $  - 5.0000000000000000d-01*HR1(1)*HR1(1)*HR3(-1,-1,1)
     $  - 5.0000000000000000d-01*HR1(1)*HR1(1)*HR3(0,-1,-1)
     $  - 5.0000000000000000d-01*HR1(1)*HR1(1)*HR3(0,0,-1)
     $  - 5.0000000000000000d-01*HR1(1)*HR1(1)*HR3(0,0,1)
     $  - 5.0000000000000000d-01*HR1(1)*HR1(1)*HR3(0,1,-1)
     $  - 6.9314718055994530d-01*HR1(1)*HR3(-1,-1,1)
     $  + 6.9314718055994530d-01*HR1(1)*HR3(-1,1,1)
     $  - 6.9314718055994530d-01*HR1(1)*HR3(0,-1,-1)
     $  - 6.9314718055994530d-01*HR1(1)*HR3(0,-1,1)
     $  - 6.9314718055994530d-01*HR1(1)*HR3(0,1,-1)
     $  - 6.9314718055994530d-01*HR1(1)*HR3(0,1,1)
     $  - 2.0000000000000000d+00*HR1(1)*HR4(-1,-1,-1,1)
     $  + HR1(1) *HR4(-1,-1,1,1)
     $  + 2.0000000000000000d+00*HR1(1)*HR4(0,-1,-1,-1)
     $  + HR1(1) *HR4(0,-1,-1,1)
     $  + HR1(1) *HR4(0,-1,1,-1)
     $  + HR1(1) *HR4(0,0,-1,-1)
     $  + HR1(1) *HR4(0,0,-1,1)
     $  + HR1(1) *HR4(0,0,1,-1)
     $  + HR1(1) *HR4(0,0,1,1)
     $  + 2.0000000000000000d+00*HR1(1)*HR4(0,1,-1,-1)
     $  + HR1(1) *HR4(0,1,-1,1)
     $  + HR1(1) *HR4(0,1,1,-1)
     $  - 6.9314718055994530d-01*HR4(-1,-1,-1,1)
     $  + 6.9314718055994530d-01*HR4(-1,-1,1,1)
     $  - 6.9314718055994530d-01*HR4(-1,1,1,1)
     $  + 6.9314718055994530d-01*HR4(0,-1,-1,-1)
     $  + 6.9314718055994530d-01*HR4(0,-1,-1,1)
     $  + 6.9314718055994530d-01*HR4(0,-1,1,-1)
     $  + 6.9314718055994530d-01*HR4(0,-1,1,1)
     $  + 6.9314718055994530d-01*HR4(0,1,-1,-1)
     $  + 6.9314718055994530d-01*HR4(0,1,-1,1)
     $  + 6.9314718055994530d-01*HR4(0,1,1,-1)
     $  + 6.9314718055994530d-01*HR4(0,1,1,1)
     $  - 3.0000000000000000d+00*HR5(-1,-1,-1,-1,1)
     $  + 2.0000000000000000d+00*HR5(-1,-1,-1,1,1)
     $  - HR5( -1,-1,1,1,1)
     $  - 3.0000000000000000d+00*HR5(0,-1,-1,-1,-1)
     $  - 2.0000000000000000d+00*HR5(0,-1,-1,-1,1)
     $  - 2.0000000000000000d+00*HR5(0,-1,-1,1,-1)
     $  - HR5(0, -1,-1,1,1)
     $  - 2.0000000000000000d+00*HR5(0,-1,1,-1,-1)
     $  - HR5(0, -1,1,-1,1)
     $  - HR5(0, -1,1,1,-1)
     $  - HR5(0,0, -1,-1,-1)
     $  - HR5(0,0, -1,-1,1)
     $  - HR5(0,0, -1,1,-1)
     $  - HR5(0,0, -1,1,1)
     $  - HR5(0,0,1, -1,-1)
     $  - HR5(0,0,1, -1,1)
     $  - HR5(0,0,1,1, -1)
     $  - HR5(0,0,1,1,1) 
     $  - 3.0000000000000000d+00*HR5(0,1,-1,-1,-1)
     $  - 2.0000000000000000d+00*HR5(0,1,-1,-1,1)
     $  - 2.0000000000000000d+00*HR5(0,1,-1,1,-1)
     $  - HR5(0,1, -1,1,1)
     $  - 2.0000000000000000d+00*HR5(0,1,1,-1,-1)
     $  - HR5(0,1,1, -1,1)
     $  - HR5(0,1,1,1, -1)
      HY5(0,0,1,0,1) = 
     $  + 2.2881039760335375d-01
     $  - 8.1174242528335364d-01*HR1(-1)
     $  - 1.7721476084810206d+00*HR1(-1)*HR1(-1)
     $  + 2.7415567780803773d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 2.8881132523331054d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 8.3333333333333333d-03*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 8.3333333333333333d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR2(0,1)
     $  - 8.2246703342411321d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  + 2.5000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR2(0,-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR2(0,1)
     $  + 8.2246703342411321d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  + 1.7328679513998632d-01*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR2(-1,1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR2(0,-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR2(0,1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR2(0,1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(-1,-1,1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(-1,1,1)
     $  + HR1( -1)*HR1(-1)*HR3(0,-1,-1)
     $  + HR1( -1)*HR1(-1)*HR3(0,0,-1)
     $  + HR1( -1)*HR1(-1)*HR3(0,0,1)
     $  + HR1( -1)*HR1(-1)*HR3(0,1,-1)
     $  - 1.6449340668482264d+00*HR1(-1)*HR1(0)*HR1(1)
     $  - HR1( -1)*HR1(0)*HR1(1)*HR2(-1,1)
     $  - HR1( -1)*HR1(0)*HR1(1)*HR2(0,-1)
     $  - HR1( -1)*HR1(0)*HR1(1)*HR2(0,1)
     $  - HR1( -1)*HR1(0)*HR3(-1,-1,1)
     $  + HR1( -1)*HR1(0)*HR3(-1,1,1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR1(0)*HR3(0,-1,-1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR1(0)*HR3(0,-1,1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR1(0)*HR3(0,1,-1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR1(0)*HR3(0,1,1)
     $  - 3.5442952169620413d+00*HR1(-1)*HR1(1)
     $  + 8.2246703342411321d-01*HR1(-1)*HR1(1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(1)*HR1(1)*HR2(-1,1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(1)*HR1(1)*HR2(0,-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(1)*HR1(1)*HR2(0,1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(1)*HR2(-1,1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(1)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(1)*HR2(0,1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR1(1)*HR3(-1,-1,1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR1(1)*HR3(-1,1,1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR1(1)*HR3(0,-1,-1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR1(1)*HR3(0,0,-1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR1(1)*HR3(0,0,1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR1(1)*HR3(0,1,-1)
     $  + 1.6449340668482264d+00*HR1(-1)*HR2(0,-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR2(0,-1)*HR2(0,-1)
     $  - HR1( -1)*HR2(0,-1)*HR2(0,1)
     $  + 1.6449340668482264d+00*HR1(-1)*HR2(0,1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR2(0,1)*HR2(0,1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR3(-1,-1,1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR3(-1,1,1)
     $  + 1.3862943611198906d+00*HR1(-1)*HR3(0,-1,-1)
     $  + 1.3862943611198906d+00*HR1(-1)*HR3(0,-1,1)
     $  + 1.3862943611198906d+00*HR1(-1)*HR3(0,1,-1)
     $  + 1.3862943611198906d+00*HR1(-1)*HR3(0,1,1)
     $  - 3.0000000000000000d+00*HR1(-1)*HR4(-1,-1,-1,1)
     $  + 3.0000000000000000d+00*HR1(-1)*HR4(-1,-1,1,1)
     $  - 3.0000000000000000d+00*HR1(-1)*HR4(-1,1,1,1)
     $  - 4.0000000000000000d+00*HR1(-1)*HR4(0,-1,-1,-1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR4(0,-1,-1,1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR4(0,-1,1,-1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR4(0,0,-1,-1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR4(0,0,-1,1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR4(0,0,1,-1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR4(0,0,1,1)
     $  - 4.0000000000000000d+00*HR1(-1)*HR4(0,1,-1,-1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR4(0,1,-1,1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR4(0,1,1,-1)
     $  - 8.2246703342411321d-01*HR1(0)*HR1(1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(0)*HR1(1)*HR1(1)*HR2(-1,1)
     $  - 5.0000000000000000d-01*HR1(0)*HR1(1)*HR1(1)*HR2(0,-1)
     $  - 5.0000000000000000d-01*HR1(0)*HR1(1)*HR1(1)*HR2(0,1)
     $  + 2.0000000000000000d+00*HR1(0)*HR1(1)*HR3(-1,-1,1)
     $  - 2.0000000000000000d+00*HR1(0)*HR1(1)*HR3(-1,1,1)
     $  + 2.0000000000000000d+00*HR1(0)*HR1(1)*HR3(0,-1,-1)
     $  + 2.0000000000000000d+00*HR1(0)*HR1(1)*HR3(0,-1,1)
     $  + 2.0000000000000000d+00*HR1(0)*HR1(1)*HR3(0,1,-1)
     $  + 2.0000000000000000d+00*HR1(0)*HR1(1)*HR3(0,1,1)
     $  + 3.0000000000000000d+00*HR1(0)*HR4(-1,-1,-1,1)
     $  - 3.0000000000000000d+00*HR1(0)*HR4(-1,-1,1,1)
     $  + 3.0000000000000000d+00*HR1(0)*HR4(-1,1,1,1)
     $  - 3.0000000000000000d+00*HR1(0)*HR4(0,-1,-1,-1)
     $  - 3.0000000000000000d+00*HR1(0)*HR4(0,-1,-1,1)
     $  - 3.0000000000000000d+00*HR1(0)*HR4(0,-1,1,-1)
     $  - 3.0000000000000000d+00*HR1(0)*HR4(0,-1,1,1)
     $  - 3.0000000000000000d+00*HR1(0)*HR4(0,1,-1,-1)
     $  - 3.0000000000000000d+00*HR1(0)*HR4(0,1,-1,1)
     $  - 3.0000000000000000d+00*HR1(0)*HR4(0,1,1,-1)
     $  - 3.0000000000000000d+00*HR1(0)*HR4(0,1,1,1)
     $  - 8.1174242528335364d-01*HR1(1)
     $  - 1.7721476084810206d+00*HR1(1)*HR1(1)
     $  + 3.4657359027997265d-01*HR1(1)*HR1(1)*HR2(-1,1)
     $  - 3.4657359027997265d-01*HR1(1)*HR1(1)*HR2(0,-1)
     $  - 3.4657359027997265d-01*HR1(1)*HR1(1)*HR2(0,1)
     $  + HR1(1) *HR1(1)*HR3(-1,-1,1)
     $  + HR1(1) *HR1(1)*HR3(0,-1,-1)
     $  + HR1(1) *HR1(1)*HR3(0,0,-1)
     $  + HR1(1) *HR1(1)*HR3(0,0,1)
     $  + HR1(1) *HR1(1)*HR3(0,1,-1)
     $  - 1.6449340668482264d+00*HR1(1)*HR2(-1,1)
     $  - 5.0000000000000000d-01*HR1(1)*HR2(-1,1)*HR2(-1,1)
     $  + HR1(1) *HR2(-1,1)*HR2(0,-1)
     $  + HR1(1) *HR2(-1,1)*HR2(0,1)
     $  + 1.6449340668482264d+00*HR1(1)*HR2(0,-1)
     $  - 5.0000000000000000d-01*HR1(1)*HR2(0,-1)*HR2(0,-1)
     $  - HR1(1) *HR2(0,-1)*HR2(0,1)
     $  + 1.6449340668482264d+00*HR1(1)*HR2(0,1)
     $  - 5.0000000000000000d-01*HR1(1)*HR2(0,1)*HR2(0,1)
     $  + 1.3862943611198906d+00*HR1(1)*HR3(-1,-1,1)
     $  - 1.3862943611198906d+00*HR1(1)*HR3(-1,1,1)
     $  + 1.3862943611198906d+00*HR1(1)*HR3(0,-1,-1)
     $  + 1.3862943611198906d+00*HR1(1)*HR3(0,-1,1)
     $  + 1.3862943611198906d+00*HR1(1)*HR3(0,1,-1)
     $  + 1.3862943611198906d+00*HR1(1)*HR3(0,1,1)
     $  + 4.0000000000000000d+00*HR1(1)*HR4(-1,-1,-1,1)
     $  - 2.0000000000000000d+00*HR1(1)*HR4(-1,-1,1,1)
     $  - 4.0000000000000000d+00*HR1(1)*HR4(0,-1,-1,-1)
     $  - 2.0000000000000000d+00*HR1(1)*HR4(0,-1,-1,1)
     $  - 2.0000000000000000d+00*HR1(1)*HR4(0,-1,1,-1)
     $  - 2.0000000000000000d+00*HR1(1)*HR4(0,0,-1,-1)
     $  - 2.0000000000000000d+00*HR1(1)*HR4(0,0,-1,1)
     $  - 2.0000000000000000d+00*HR1(1)*HR4(0,0,1,-1)
     $  - 2.0000000000000000d+00*HR1(1)*HR4(0,0,1,1)
     $  - 4.0000000000000000d+00*HR1(1)*HR4(0,1,-1,-1)
     $  - 2.0000000000000000d+00*HR1(1)*HR4(0,1,-1,1)
     $  - 2.0000000000000000d+00*HR1(1)*HR4(0,1,1,-1)
     $  + HR2( -1,1)*HR3(-1,1,1)
     $  + HR2(0, -1)*HR3(-1,-1,1)
     $  - HR2(0, -1)*HR3(-1,1,1)
     $  + HR2(0, -1)*HR3(0,-1,-1)
     $  + HR2(0, -1)*HR3(0,-1,1)
     $  + HR2(0, -1)*HR3(0,1,-1)
     $  + HR2(0, -1)*HR3(0,1,1)
     $  + HR2(0,1) *HR3(-1,-1,1)
     $  - HR2(0,1) *HR3(-1,1,1)
     $  + HR2(0,1) *HR3(0,1,-1)
     $  + HR2(0,1) *HR3(0,1,1)
     $  - 1.6449340668482264d+00*HR3(-1,-1,1)
     $  + 1.6449340668482264d+00*HR3(-1,1,1)
     $  - 1.6449340668482264d+00*HR3(0,-1,-1)
     $  - 1.6449340668482264d+00*HR3(0,-1,1)
     $  - 1.6449340668482264d+00*HR3(0,1,-1)
     $  - 1.6449340668482264d+00*HR3(0,1,1)
     $  + 2.0794415416798359d+00*HR4(-1,-1,-1,1)
     $  - 2.0794415416798359d+00*HR4(-1,-1,1,1)
     $  + 2.0794415416798359d+00*HR4(-1,1,1,1)
     $  - 2.0794415416798359d+00*HR4(0,-1,-1,-1)
     $  - 2.0794415416798359d+00*HR4(0,-1,-1,1)
     $  - 2.0794415416798359d+00*HR4(0,-1,1,-1)
     $  - 2.0794415416798359d+00*HR4(0,-1,1,1)
     $  - 2.0794415416798359d+00*HR4(0,1,-1,-1)
     $  - 2.0794415416798359d+00*HR4(0,1,-1,1)
     $  - 2.0794415416798359d+00*HR4(0,1,1,-1)
     $  - 2.0794415416798359d+00*HR4(0,1,1,1)
     $  + 7.0000000000000000d+00*HR5(-1,-1,-1,-1,1)
     $  - 7.0000000000000000d+00*HR5(-1,-1,-1,1,1)
     $  - HR5( -1,-1,1,-1,1)
     $  - HR5( -1,1,-1,1,1)
     $  + 7.0000000000000000d+00*HR5(0,-1,-1,-1,-1)
     $  + 4.0000000000000000d+00*HR5(0,-1,-1,-1,1)
     $  + HR5(0, -1,-1,0,1)
     $  + 5.0000000000000000d+00*HR5(0,-1,-1,1,-1)
     $  + 2.0000000000000000d+00*HR5(0,-1,-1,1,1)
     $  - HR5(0, -1,0,-1,-1)
     $  + HR5(0, -1,0,1,1)
     $  + 5.0000000000000000d+00*HR5(0,-1,1,-1,-1)
     $  + 2.0000000000000000d+00*HR5(0,-1,1,-1,1)
     $  + HR5(0, -1,1,0,1)
     $  + 3.0000000000000000d+00*HR5(0,-1,1,1,-1)
     $  + 2.0000000000000000d+00*HR5(0,0,-1,-1,1)
     $  + 2.0000000000000000d+00*HR5(0,0,-1,1,-1)
     $  + 4.0000000000000000d+00*HR5(0,0,-1,1,1)
     $  + 2.0000000000000000d+00*HR5(0,0,1,-1,-1)
     $  + 2.0000000000000000d+00*HR5(0,0,1,-1,1)
     $  + 7.0000000000000000d+00*HR5(0,1,-1,-1,-1)
     $  + 4.0000000000000000d+00*HR5(0,1,-1,-1,1)
     $  + 5.0000000000000000d+00*HR5(0,1,-1,1,-1)
     $  + 2.0000000000000000d+00*HR5(0,1,-1,1,1)
     $  - HR5(0,1,0,1, -1)
     $  - HR5(0,1,0,1,1) 
     $  + 5.0000000000000000d+00*HR5(0,1,1,-1,-1)
     $  + 2.0000000000000000d+00*HR5(0,1,1,-1,1)
     $  + 3.0000000000000000d+00*HR5(0,1,1,1,-1)
      HY5(0,0,1,1,1) = 
     $  + 9.6551159989443734d-02
     $  - 1.0823232337111381d+00*HR1(-1)
     $  - 2.7752054332410789d-02*HR1(-1)*HR1(-1)
     $  + 4.0037751159850118d-02*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 2.8881132523331054d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 8.3333333333333333d-03*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  + 8.3333333333333333d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)*HR1(0)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 8.3333333333333333d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  - 1.2011325347955035d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  - 1.7328679513998632d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(0)
     $  - 8.3333333333333333d-02*HR1(-1)*HR1(-1)*HR1(0)*HR1(0)*HR1(0)
     $  + 2.5000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(0)*HR1(1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  - 2.5000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)*HR1(1)
     $  + 1.2011325347955035d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - 1.7328679513998632d-01*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR2(-1,1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(-1,-1,1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(-1,1,1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(0)*HR1(0)*HR1(0)*HR1(1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(0)*HR1(0)*HR1(1)
     $  + 2.5000000000000000d-01*HR1(-1)*HR1(0)*HR1(0)*HR1(1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(0)*HR1(0)*HR2(0,-1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(0)*HR1(0)*HR2(0,1)
     $  - 2.4022650695910071d-01*HR1(-1)*HR1(0)*HR1(1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(0)*HR1(1)*HR1(1)
     $  + HR1( -1)*HR1(0)*HR1(1)*HR2(-1,1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(0)*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(0)*HR2(0,1)
     $  + HR1( -1)*HR1(0)*HR3(-1,-1,1)
     $  - HR1( -1)*HR1(0)*HR3(-1,1,1)
     $  - HR1( -1)*HR1(0)*HR3(0,-1,-1)
     $  - HR1( -1)*HR1(0)*HR3(0,0,-1)
     $  - HR1( -1)*HR1(0)*HR3(0,0,1)
     $  - HR1( -1)*HR1(0)*HR3(0,1,-1)
     $  - 5.5504108664821579d-02*HR1(-1)*HR1(1)
     $  + 1.2011325347955035d-01*HR1(-1)*HR1(1)*HR1(1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(1)*HR2(-1,1)
     $  + HR1( -1)*HR1(1)*HR3(-1,-1,1)
     $  + 2.4022650695910071d-01*HR1(-1)*HR2(0,-1)
     $  + 2.4022650695910071d-01*HR1(-1)*HR2(0,1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR3(-1,-1,1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR3(-1,1,1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR3(0,-1,-1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR3(0,0,-1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR3(0,0,1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR3(0,1,-1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR4(-1,-1,-1,1)
     $  - HR1( -1)*HR4(-1,-1,1,1)
     $  + HR1( -1)*HR4(0,-1,-1,-1)
     $  + HR1( -1)*HR4(0,0,-1,-1)
     $  + HR1( -1)*HR4(0,0,0,-1)
     $  + HR1( -1)*HR4(0,0,0,1)
     $  + HR1( -1)*HR4(0,0,1,-1)
     $  + HR1( -1)*HR4(0,1,-1,-1)
     $  - 8.3333333333333333d-02*HR1(0)*HR1(0)*HR1(0)*HR1(1)*HR1(1)
     $  - 1.7328679513998632d-01*HR1(0)*HR1(0)*HR1(1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(0)*HR1(0)*HR1(1)*HR2(-1,1)
     $  + 5.0000000000000000d-01*HR1(0)*HR1(0)*HR1(1)*HR2(0,-1)
     $  + 5.0000000000000000d-01*HR1(0)*HR1(0)*HR1(1)*HR2(0,1)
     $  - 5.0000000000000000d-01*HR1(0)*HR1(0)*HR3(-1,-1,1)
     $  + 5.0000000000000000d-01*HR1(0)*HR1(0)*HR3(-1,1,1)
     $  - 5.0000000000000000d-01*HR1(0)*HR1(0)*HR3(0,-1,-1)
     $  - 5.0000000000000000d-01*HR1(0)*HR1(0)*HR3(0,-1,1)
     $  - 5.0000000000000000d-01*HR1(0)*HR1(0)*HR3(0,1,-1)
     $  - 5.0000000000000000d-01*HR1(0)*HR1(0)*HR3(0,1,1)
     $  - 1.2011325347955035d-01*HR1(0)*HR1(1)*HR1(1)
     $  - 6.9314718055994530d-01*HR1(0)*HR1(1)*HR2(-1,1)
     $  + 6.9314718055994530d-01*HR1(0)*HR1(1)*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR1(0)*HR1(1)*HR2(0,1)
     $  - HR1(0) *HR1(1)*HR3(-1,-1,1)
     $  - HR1(0) *HR1(1)*HR3(0,-1,-1)
     $  - HR1(0) *HR1(1)*HR3(0,0,-1)
     $  - HR1(0) *HR1(1)*HR3(0,0,1)
     $  - HR1(0) *HR1(1)*HR3(0,1,-1)
     $  - 6.9314718055994530d-01*HR1(0)*HR3(-1,-1,1)
     $  + 6.9314718055994530d-01*HR1(0)*HR3(-1,1,1)
     $  - 6.9314718055994530d-01*HR1(0)*HR3(0,-1,-1)
     $  - 6.9314718055994530d-01*HR1(0)*HR3(0,-1,1)
     $  - 6.9314718055994530d-01*HR1(0)*HR3(0,1,-1)
     $  - 6.9314718055994530d-01*HR1(0)*HR3(0,1,1)
     $  - 2.0000000000000000d+00*HR1(0)*HR4(-1,-1,-1,1)
     $  + HR1(0) *HR4(-1,-1,1,1)
     $  + 2.0000000000000000d+00*HR1(0)*HR4(0,-1,-1,-1)
     $  + HR1(0) *HR4(0,-1,-1,1)
     $  + HR1(0) *HR4(0,-1,1,-1)
     $  + HR1(0) *HR4(0,0,-1,-1)
     $  + HR1(0) *HR4(0,0,-1,1)
     $  + HR1(0) *HR4(0,0,1,-1)
     $  + HR1(0) *HR4(0,0,1,1)
     $  + 2.0000000000000000d+00*HR1(0)*HR4(0,1,-1,-1)
     $  + HR1(0) *HR4(0,1,-1,1)
     $  + HR1(0) *HR4(0,1,1,-1)
     $  - 1.0823232337111381d+00*HR1(1)
     $  - 2.7752054332410789d-02*HR1(1)*HR1(1)
     $  - 2.4022650695910071d-01*HR1(1)*HR2(-1,1)
     $  + 2.4022650695910071d-01*HR1(1)*HR2(0,-1)
     $  + 2.4022650695910071d-01*HR1(1)*HR2(0,1)
     $  - 6.9314718055994530d-01*HR1(1)*HR3(-1,-1,1)
     $  - 6.9314718055994530d-01*HR1(1)*HR3(0,-1,-1)
     $  - 6.9314718055994530d-01*HR1(1)*HR3(0,0,-1)
     $  - 6.9314718055994530d-01*HR1(1)*HR3(0,0,1)
     $  - 6.9314718055994530d-01*HR1(1)*HR3(0,1,-1)
     $  - HR1(1) *HR4(-1,-1,-1,1)
     $  + HR1(1) *HR4(0,-1,-1,-1)
     $  + HR1(1) *HR4(0,0,-1,-1)
     $  + HR1(1) *HR4(0,0,0,-1)
     $  + HR1(1) *HR4(0,0,0,1)
     $  + HR1(1) *HR4(0,0,1,-1)
     $  + HR1(1) *HR4(0,1,-1,-1)
     $  - 2.4022650695910071d-01*HR3(-1,-1,1)
     $  + 2.4022650695910071d-01*HR3(-1,1,1)
     $  - 2.4022650695910071d-01*HR3(0,-1,-1)
     $  - 2.4022650695910071d-01*HR3(0,-1,1)
     $  - 2.4022650695910071d-01*HR3(0,1,-1)
     $  - 2.4022650695910071d-01*HR3(0,1,1)
     $  - 1.3862943611198906d+00*HR4(-1,-1,-1,1)
     $  + 6.9314718055994530d-01*HR4(-1,-1,1,1)
     $  + 1.3862943611198906d+00*HR4(0,-1,-1,-1)
     $  + 6.9314718055994530d-01*HR4(0,-1,-1,1)
     $  + 6.9314718055994530d-01*HR4(0,-1,1,-1)
     $  + 6.9314718055994530d-01*HR4(0,0,-1,-1)
     $  + 6.9314718055994530d-01*HR4(0,0,-1,1)
     $  + 6.9314718055994530d-01*HR4(0,0,1,-1)
     $  + 6.9314718055994530d-01*HR4(0,0,1,1)
     $  + 1.3862943611198906d+00*HR4(0,1,-1,-1)
     $  + 6.9314718055994530d-01*HR4(0,1,-1,1)
     $  + 6.9314718055994530d-01*HR4(0,1,1,-1)
     $  - 3.0000000000000000d+00*HR5(-1,-1,-1,-1,1)
     $  + HR5( -1,-1,-1,1,1)
     $  - 3.0000000000000000d+00*HR5(0,-1,-1,-1,-1)
     $  - HR5(0, -1,-1,-1,1)
     $  - HR5(0, -1,-1,1,-1)
     $  - HR5(0, -1,1,-1,-1)
     $  - 2.0000000000000000d+00*HR5(0,0,-1,-1,-1)
     $  - HR5(0,0, -1,-1,1)
     $  - HR5(0,0, -1,1,-1)
     $  - HR5(0,0,0, -1,-1)
     $  - HR5(0,0,0, -1,1)
     $  - HR5(0,0,0,1, -1)
     $  - HR5(0,0,0,1,1) 
     $  - 2.0000000000000000d+00*HR5(0,0,1,-1,-1)
     $  - HR5(0,0,1, -1,1)
     $  - HR5(0,0,1,1, -1)
     $  - 3.0000000000000000d+00*HR5(0,1,-1,-1,-1)
     $  - HR5(0,1, -1,-1,1)
     $  - HR5(0,1, -1,1,-1)
     $  - HR5(0,1,1, -1,-1)
      HY5(0,1,0,1,1) = 
     $  + 2.2881039760335375d-01
     $  + 4.0801720544311065d+00*HR1(-1)
     $  - 6.0102845157979714d-01*HR1(-1)*HR1(-1)
     $  + 4.0037751159850118d-02*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 2.8881132523331054d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 8.3333333333333333d-03*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  + 8.3333333333333333d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)*HR1(0)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  + 2.5000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(0)*HR1(1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR2(-1,1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR2(0,-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR2(0,1)
     $  + 1.2011325347955035d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR2(-1,1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR2(0,1)
     $  + HR1( -1)*HR1(-1)*HR3(-1,-1,1)
     $  - HR1( -1)*HR1(-1)*HR3(-1,1,1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(0,-1,-1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(0,0,-1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(0,0,1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(0,1,-1)
     $  + 1.2020569031595942d+00*HR1(-1)*HR1(0)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(0)*HR1(0)*HR2(-1,1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(0)*HR1(0)*HR2(0,-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(0)*HR1(0)*HR2(0,1)
     $  - HR1( -1)*HR1(0)*HR1(1)*HR2(-1,1)
     $  - HR1( -1)*HR1(0)*HR1(1)*HR2(0,-1)
     $  - HR1( -1)*HR1(0)*HR1(1)*HR2(0,1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(0)*HR2(-1,1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(0)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(0)*HR2(0,1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR1(0)*HR3(-1,-1,1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR1(0)*HR3(-1,1,1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR1(0)*HR3(0,-1,-1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR1(0)*HR3(0,0,-1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR1(0)*HR3(0,0,1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR1(0)*HR3(0,1,-1)
     $  - 1.2020569031595942d+00*HR1(-1)*HR1(1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(1)*HR2(-1,1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(1)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(1)*HR2(0,1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR1(1)*HR3(-1,-1,1)
     $  + HR1( -1)*HR1(1)*HR3(0,-1,-1)
     $  + HR1( -1)*HR1(1)*HR3(0,0,-1)
     $  + HR1( -1)*HR1(1)*HR3(0,0,1)
     $  + HR1( -1)*HR1(1)*HR3(0,1,-1)
     $  - 2.4022650695910071d-01*HR1(-1)*HR2(-1,1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR2(-1,1)*HR2(-1,1)
     $  - 2.4022650695910071d-01*HR1(-1)*HR2(0,-1)
     $  - 2.4022650695910071d-01*HR1(-1)*HR2(0,1)
     $  - 1.3862943611198906d+00*HR1(-1)*HR3(-1,-1,1)
     $  + 1.3862943611198906d+00*HR1(-1)*HR3(-1,1,1)
     $  + 1.3862943611198906d+00*HR1(-1)*HR3(0,-1,-1)
     $  + 1.3862943611198906d+00*HR1(-1)*HR3(0,0,-1)
     $  + 1.3862943611198906d+00*HR1(-1)*HR3(0,0,1)
     $  + 1.3862943611198906d+00*HR1(-1)*HR3(0,1,-1)
     $  - 4.0000000000000000d+00*HR1(-1)*HR4(-1,-1,-1,1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR4(-1,-1,1,1)
     $  - 3.0000000000000000d+00*HR1(-1)*HR4(0,-1,-1,-1)
     $  - 3.0000000000000000d+00*HR1(-1)*HR4(0,0,-1,-1)
     $  - 3.0000000000000000d+00*HR1(-1)*HR4(0,0,0,-1)
     $  - 3.0000000000000000d+00*HR1(-1)*HR4(0,0,0,1)
     $  - 3.0000000000000000d+00*HR1(-1)*HR4(0,0,1,-1)
     $  - 3.0000000000000000d+00*HR1(-1)*HR4(0,1,-1,-1)
     $  + 5.0000000000000000d-01*HR1(0)*HR1(0)*HR1(1)*HR2(-1,1)
     $  - 5.0000000000000000d-01*HR1(0)*HR1(0)*HR1(1)*HR2(0,-1)
     $  - 5.0000000000000000d-01*HR1(0)*HR1(0)*HR1(1)*HR2(0,1)
     $  + HR1(0) *HR1(0)*HR3(-1,-1,1)
     $  - HR1(0) *HR1(0)*HR3(-1,1,1)
     $  + HR1(0) *HR1(0)*HR3(0,-1,-1)
     $  + HR1(0) *HR1(0)*HR3(0,-1,1)
     $  + HR1(0) *HR1(0)*HR3(0,1,-1)
     $  + HR1(0) *HR1(0)*HR3(0,1,1)
     $  + 1.2020569031595942d+00*HR1(0)*HR1(1)
     $  + 6.9314718055994530d-01*HR1(0)*HR1(1)*HR2(-1,1)
     $  - 6.9314718055994530d-01*HR1(0)*HR1(1)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(0)*HR1(1)*HR2(0,1)
     $  + 2.0000000000000000d+00*HR1(0)*HR1(1)*HR3(-1,-1,1)
     $  + 2.0000000000000000d+00*HR1(0)*HR1(1)*HR3(0,-1,-1)
     $  + 2.0000000000000000d+00*HR1(0)*HR1(1)*HR3(0,0,-1)
     $  + 2.0000000000000000d+00*HR1(0)*HR1(1)*HR3(0,0,1)
     $  + 2.0000000000000000d+00*HR1(0)*HR1(1)*HR3(0,1,-1)
     $  - 5.0000000000000000d-01*HR1(0)*HR2(-1,1)*HR2(-1,1)
     $  + HR1(0) *HR2(-1,1)*HR2(0,-1)
     $  + HR1(0) *HR2(-1,1)*HR2(0,1)
     $  - 5.0000000000000000d-01*HR1(0)*HR2(0,-1)*HR2(0,-1)
     $  - HR1(0) *HR2(0,-1)*HR2(0,1)
     $  - 5.0000000000000000d-01*HR1(0)*HR2(0,1)*HR2(0,1)
     $  + 1.3862943611198906d+00*HR1(0)*HR3(-1,-1,1)
     $  - 1.3862943611198906d+00*HR1(0)*HR3(-1,1,1)
     $  + 1.3862943611198906d+00*HR1(0)*HR3(0,-1,-1)
     $  + 1.3862943611198906d+00*HR1(0)*HR3(0,-1,1)
     $  + 1.3862943611198906d+00*HR1(0)*HR3(0,1,-1)
     $  + 1.3862943611198906d+00*HR1(0)*HR3(0,1,1)
     $  + 4.0000000000000000d+00*HR1(0)*HR4(-1,-1,-1,1)
     $  - 2.0000000000000000d+00*HR1(0)*HR4(-1,-1,1,1)
     $  - 4.0000000000000000d+00*HR1(0)*HR4(0,-1,-1,-1)
     $  - 2.0000000000000000d+00*HR1(0)*HR4(0,-1,-1,1)
     $  - 2.0000000000000000d+00*HR1(0)*HR4(0,-1,1,-1)
     $  - 2.0000000000000000d+00*HR1(0)*HR4(0,0,-1,-1)
     $  - 2.0000000000000000d+00*HR1(0)*HR4(0,0,-1,1)
     $  - 2.0000000000000000d+00*HR1(0)*HR4(0,0,1,-1)
     $  - 2.0000000000000000d+00*HR1(0)*HR4(0,0,1,1)
     $  - 4.0000000000000000d+00*HR1(0)*HR4(0,1,-1,-1)
     $  - 2.0000000000000000d+00*HR1(0)*HR4(0,1,-1,1)
     $  - 2.0000000000000000d+00*HR1(0)*HR4(0,1,1,-1)
     $  + 4.0801720544311065d+00*HR1(1)
     $  + 2.4022650695910071d-01*HR1(1)*HR2(-1,1)
     $  - 2.4022650695910071d-01*HR1(1)*HR2(0,-1)
     $  - 2.4022650695910071d-01*HR1(1)*HR2(0,1)
     $  + 1.3862943611198906d+00*HR1(1)*HR3(-1,-1,1)
     $  + 1.3862943611198906d+00*HR1(1)*HR3(0,-1,-1)
     $  + 1.3862943611198906d+00*HR1(1)*HR3(0,0,-1)
     $  + 1.3862943611198906d+00*HR1(1)*HR3(0,0,1)
     $  + 1.3862943611198906d+00*HR1(1)*HR3(0,1,-1)
     $  + 3.0000000000000000d+00*HR1(1)*HR4(-1,-1,-1,1)
     $  - 3.0000000000000000d+00*HR1(1)*HR4(0,-1,-1,-1)
     $  - 3.0000000000000000d+00*HR1(1)*HR4(0,0,-1,-1)
     $  - 3.0000000000000000d+00*HR1(1)*HR4(0,0,0,-1)
     $  - 3.0000000000000000d+00*HR1(1)*HR4(0,0,0,1)
     $  - 3.0000000000000000d+00*HR1(1)*HR4(0,0,1,-1)
     $  - 3.0000000000000000d+00*HR1(1)*HR4(0,1,-1,-1)
     $  + 1.2020569031595942d+00*HR2(-1,1)
     $  - 3.4657359027997265d-01*HR2(-1,1)*HR2(-1,1)
     $  + 6.9314718055994530d-01*HR2(-1,1)*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR2(-1,1)*HR2(0,1)
     $  - HR2( -1,1)*HR3(-1,-1,1)
     $  - HR2( -1,1)*HR3(0,-1,-1)
     $  - HR2( -1,1)*HR3(0,0,-1)
     $  - HR2( -1,1)*HR3(0,0,1)
     $  - HR2( -1,1)*HR3(0,1,-1)
     $  - 1.2020569031595942d+00*HR2(0,-1)
     $  - 3.4657359027997265d-01*HR2(0,-1)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR2(0,-1)*HR2(0,1)
     $  + HR2(0, -1)*HR3(0,0,-1)
     $  + HR2(0, -1)*HR3(0,0,1)
     $  - 1.2020569031595942d+00*HR2(0,1)
     $  - 3.4657359027997265d-01*HR2(0,1)*HR2(0,1)
     $  + HR2(0,1) *HR3(0,-1,-1)
     $  + HR2(0,1) *HR3(0,0,-1)
     $  + HR2(0,1) *HR3(0,0,1)
     $  + 4.8045301391820142d-01*HR3(-1,-1,1)
     $  - 4.8045301391820142d-01*HR3(-1,1,1)
     $  + 4.8045301391820142d-01*HR3(0,-1,-1)
     $  + 4.8045301391820142d-01*HR3(0,-1,1)
     $  + 4.8045301391820142d-01*HR3(0,1,-1)
     $  + 4.8045301391820142d-01*HR3(0,1,1)
     $  + 2.7725887222397812d+00*HR4(-1,-1,-1,1)
     $  - 1.3862943611198906d+00*HR4(-1,-1,1,1)
     $  - 2.7725887222397812d+00*HR4(0,-1,-1,-1)
     $  - 1.3862943611198906d+00*HR4(0,-1,-1,1)
     $  - 1.3862943611198906d+00*HR4(0,-1,1,-1)
     $  - 1.3862943611198906d+00*HR4(0,0,-1,-1)
     $  - 1.3862943611198906d+00*HR4(0,0,-1,1)
     $  - 1.3862943611198906d+00*HR4(0,0,1,-1)
     $  - 1.3862943611198906d+00*HR4(0,0,1,1)
     $  - 2.7725887222397812d+00*HR4(0,1,-1,-1)
     $  - 1.3862943611198906d+00*HR4(0,1,-1,1)
     $  - 1.3862943611198906d+00*HR4(0,1,1,-1)
     $  + 7.0000000000000000d+00*HR5(-1,-1,-1,-1,1)
     $  + HR5( -1,-1,1,-1,1)
     $  + 7.0000000000000000d+00*HR5(0,-1,-1,-1,-1)
     $  + 3.0000000000000000d+00*HR5(0,-1,-1,-1,1)
     $  - HR5(0, -1,-1,0,1)
     $  + 2.0000000000000000d+00*HR5(0,-1,-1,1,-1)
     $  + HR5(0, -1,0,-1,-1)
     $  - HR5(0, -1,0,-1,1)
     $  + 2.0000000000000000d+00*HR5(0,-1,1,-1,-1)
     $  + 7.0000000000000000d+00*HR5(0,0,-1,-1,-1)
     $  + HR5(0,0, -1,-1,1)
     $  - HR5(0,0, -1,0,-1)
     $  - HR5(0,0, -1,0,1)
     $  + 2.0000000000000000d+00*HR5(0,0,-1,1,-1)
     $  + 5.0000000000000000d+00*HR5(0,0,1,-1,-1)
     $  + 3.0000000000000000d+00*HR5(0,0,1,-1,1)
     $  - HR5(0,0,1,0, -1)
     $  - HR5(0,0,1,0,1) 
     $  + 4.0000000000000000d+00*HR5(0,0,1,1,-1)
     $  + 7.0000000000000000d+00*HR5(0,1,-1,-1,-1)
     $  + 3.0000000000000000d+00*HR5(0,1,-1,-1,1)
     $  + 2.0000000000000000d+00*HR5(0,1,-1,1,-1)
     $  + HR5(0,1,0,1, -1)
     $  + 2.0000000000000000d+00*HR5(0,1,1,-1,-1)
      HY5(0,1,1,1,1) = 
     $  + 1.0369277551433699d+00
     $  - 9.6181291076284771d-03*HR1(-1)
     $  + 2.7752054332410789d-02*HR1(-1)*HR1(-1)
     $  - 4.0037751159850118d-02*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 2.8881132523331054d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 8.3333333333333333d-03*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  - 8.3333333333333333d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)*HR1(0)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  + 1.2011325347955035d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  + 1.7328679513998632d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(0)
     $  + 8.3333333333333333d-02*HR1(-1)*HR1(-1)*HR1(0)*HR1(0)*HR1(0)
     $  - 2.5000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(0)*HR1(1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR2(-1,1)
     $  - 1.2011325347955035d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(-1,-1,1)
     $  - 5.5504108664821579d-02*HR1(-1)*HR1(0)
     $  - 1.2011325347955035d-01*HR1(-1)*HR1(0)*HR1(0)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(0)*HR1(0)*HR1(0)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(0)*HR1(0)*HR1(0)*HR1(0)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(0)*HR1(0)*HR1(0)*HR1(1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(0)*HR1(0)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(0)*HR1(0)*HR2(-1,1)
     $  + 2.4022650695910071d-01*HR1(-1)*HR1(0)*HR1(1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(0)*HR2(-1,1)
     $  + HR1( -1)*HR1(0)*HR3(-1,-1,1)
     $  + 5.5504108664821579d-02*HR1(-1)*HR1(1)
     $  + 2.4022650695910071d-01*HR1(-1)*HR2(-1,1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR3(-1,-1,1)
     $  + HR1( -1)*HR4(-1,-1,-1,1)
     $  - 4.1666666666666666d-02*HR1(0)*HR1(0)*HR1(0)*HR1(0)*HR1(1)
     $  - 1.1552453009332421d-01*HR1(0)*HR1(0)*HR1(0)*HR1(1)
     $  - 1.6666666666666666d-01*HR1(0)*HR1(0)*HR1(0)*HR2(-1,1)
     $  + 1.6666666666666666d-01*HR1(0)*HR1(0)*HR1(0)*HR2(0,-1)
     $  + 1.6666666666666666d-01*HR1(0)*HR1(0)*HR1(0)*HR2(0,1)
     $  - 1.2011325347955035d-01*HR1(0)*HR1(0)*HR1(1)
     $  - 3.4657359027997265d-01*HR1(0)*HR1(0)*HR2(-1,1)
     $  + 3.4657359027997265d-01*HR1(0)*HR1(0)*HR2(0,-1)
     $  + 3.4657359027997265d-01*HR1(0)*HR1(0)*HR2(0,1)
     $  - 5.0000000000000000d-01*HR1(0)*HR1(0)*HR3(-1,-1,1)
     $  - 5.0000000000000000d-01*HR1(0)*HR1(0)*HR3(0,-1,-1)
     $  - 5.0000000000000000d-01*HR1(0)*HR1(0)*HR3(0,0,-1)
     $  - 5.0000000000000000d-01*HR1(0)*HR1(0)*HR3(0,0,1)
     $  - 5.0000000000000000d-01*HR1(0)*HR1(0)*HR3(0,1,-1)
     $  - 5.5504108664821579d-02*HR1(0)*HR1(1)
     $  - 2.4022650695910071d-01*HR1(0)*HR2(-1,1)
     $  + 2.4022650695910071d-01*HR1(0)*HR2(0,-1)
     $  + 2.4022650695910071d-01*HR1(0)*HR2(0,1)
     $  - 6.9314718055994530d-01*HR1(0)*HR3(-1,-1,1)
     $  - 6.9314718055994530d-01*HR1(0)*HR3(0,-1,-1)
     $  - 6.9314718055994530d-01*HR1(0)*HR3(0,0,-1)
     $  - 6.9314718055994530d-01*HR1(0)*HR3(0,0,1)
     $  - 6.9314718055994530d-01*HR1(0)*HR3(0,1,-1)
     $  - HR1(0) *HR4(-1,-1,-1,1)
     $  + HR1(0) *HR4(0,-1,-1,-1)
     $  + HR1(0) *HR4(0,0,-1,-1)
     $  + HR1(0) *HR4(0,0,0,-1)
     $  + HR1(0) *HR4(0,0,0,1)
     $  + HR1(0) *HR4(0,0,1,-1)
     $  + HR1(0) *HR4(0,1,-1,-1)
     $  - 9.6181291076284771d-03*HR1(1)
     $  - 5.5504108664821579d-02*HR2(-1,1)
     $  + 5.5504108664821579d-02*HR2(0,-1)
     $  + 5.5504108664821579d-02*HR2(0,1)
     $  - 2.4022650695910071d-01*HR3(-1,-1,1)
     $  - 2.4022650695910071d-01*HR3(0,-1,-1)
     $  - 2.4022650695910071d-01*HR3(0,0,-1)
     $  - 2.4022650695910071d-01*HR3(0,0,1)
     $  - 2.4022650695910071d-01*HR3(0,1,-1)
     $  - 6.9314718055994530d-01*HR4(-1,-1,-1,1)
     $  + 6.9314718055994530d-01*HR4(0,-1,-1,-1)
     $  + 6.9314718055994530d-01*HR4(0,0,-1,-1)
     $  + 6.9314718055994530d-01*HR4(0,0,0,-1)
     $  + 6.9314718055994530d-01*HR4(0,0,0,1)
     $  + 6.9314718055994530d-01*HR4(0,0,1,-1)
     $  + 6.9314718055994530d-01*HR4(0,1,-1,-1)
     $  - HR5( -1,-1,-1,-1,1)
     $  - HR5(0, -1,-1,-1,-1)
     $  - HR5(0,0, -1,-1,-1)
     $  - HR5(0,0,0, -1,-1)
     $  - HR5(0,0,0,0, -1)
     $  - HR5(0,0,0,0,1) 
     $  - HR5(0,0,0,1, -1)
     $  - HR5(0,0,1, -1,-1)
     $  - HR5(0,1, -1,-1,-1)
      if (r.lt.0d0) then 
      HY5(0,0,0,1,1) = HY5(0,0,0,1,1) 
     $  + 8.2246703342411321d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 2.4674011002723396d+00*HR1(-1)*HR1(-1)*HR1(1)
     $  + 2.4674011002723396d+00*HR1(-1)*HR1(1)*HR1(1)
     $  + 8.2246703342411321d-01*HR1(1)*HR1(1)*HR1(1)
      HY5(0,0,1,1,1) = HY5(0,0,1,1,1) 
     $  + 1.7102721159642791d+00*HR1(-1)*HR1(-1)
     $  - 8.2246703342411321d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 2.4674011002723396d+00*HR1(-1)*HR1(-1)*HR1(0)
     $  - 2.4674011002723396d+00*HR1(-1)*HR1(-1)*HR1(1)
     $  + 4.9348022005446793d+00*HR1(-1)*HR1(0)*HR1(1)
     $  + 3.4205442319285582d+00*HR1(-1)*HR1(1)
     $  - 2.4674011002723396d+00*HR1(-1)*HR1(1)*HR1(1)
     $  - 4.9348022005446793d+00*HR1(-1)*HR2(0,-1)
     $  - 4.9348022005446793d+00*HR1(-1)*HR2(0,1)
     $  + 2.4674011002723396d+00*HR1(0)*HR1(1)*HR1(1)
     $  + 1.7102721159642791d+00*HR1(1)*HR1(1)
     $  + 4.9348022005446793d+00*HR1(1)*HR2(-1,1)
     $  - 4.9348022005446793d+00*HR1(1)*HR2(0,-1)
     $  - 4.9348022005446793d+00*HR1(1)*HR2(0,1)
     $  + 4.9348022005446793d+00*HR3(-1,-1,1)
     $  - 4.9348022005446793d+00*HR3(-1,1,1)
     $  + 4.9348022005446793d+00*HR3(0,-1,-1)
     $  + 4.9348022005446793d+00*HR3(0,-1,1)
     $  + 4.9348022005446793d+00*HR3(0,1,-1)
     $  + 4.9348022005446793d+00*HR3(0,1,1)
      HY5(0,1,0,1,1) = HY5(0,1,0,1,1) 
     $  - 8.2246703342411321d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 2.4674011002723396d+00*HR1(-1)*HR1(-1)*HR1(1)
     $  + 4.9348022005446793d+00*HR1(-1)*HR2(-1,1)
     $  + 4.9348022005446793d+00*HR1(-1)*HR2(0,-1)
     $  + 4.9348022005446793d+00*HR1(-1)*HR2(0,1)
     $  - 4.9348022005446793d+00*HR1(1)*HR2(-1,1)
     $  + 4.9348022005446793d+00*HR1(1)*HR2(0,-1)
     $  + 4.9348022005446793d+00*HR1(1)*HR2(0,1)
     $  - 9.8696044010893586d+00*HR3(-1,-1,1)
     $  + 9.8696044010893586d+00*HR3(-1,1,1)
     $  - 9.8696044010893586d+00*HR3(0,-1,-1)
     $  - 9.8696044010893586d+00*HR3(0,-1,1)
     $  - 9.8696044010893586d+00*HR3(0,1,-1)
     $  - 9.8696044010893586d+00*HR3(0,1,1)
      HY5(0,1,1,1,1) = HY5(0,1,1,1,1) 
     $  - 2.8732418312458363d+00*HR1(-1)
     $  - 1.7102721159642791d+00*HR1(-1)*HR1(-1)
     $  + 8.2246703342411321d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 2.4674011002723396d+00*HR1(-1)*HR1(-1)*HR1(0)
     $  + 2.4674011002723396d+00*HR1(-1)*HR1(-1)*HR1(1)
     $  + 3.4205442319285582d+00*HR1(-1)*HR1(0)
     $  + 2.4674011002723396d+00*HR1(-1)*HR1(0)*HR1(0)
     $  - 4.9348022005446793d+00*HR1(-1)*HR1(0)*HR1(1)
     $  - 3.4205442319285582d+00*HR1(-1)*HR1(1)
     $  - 4.9348022005446793d+00*HR1(-1)*HR2(-1,1)
     $  + 2.4674011002723396d+00*HR1(0)*HR1(0)*HR1(1)
     $  + 3.4205442319285582d+00*HR1(0)*HR1(1)
     $  + 4.9348022005446793d+00*HR1(0)*HR2(-1,1)
     $  - 4.9348022005446793d+00*HR1(0)*HR2(0,-1)
     $  - 4.9348022005446793d+00*HR1(0)*HR2(0,1)
     $  - 2.8732418312458363d+00*HR1(1)
     $  + 3.4205442319285582d+00*HR2(-1,1)
     $  - 3.4205442319285582d+00*HR2(0,-1)
     $  - 3.4205442319285582d+00*HR2(0,1)
     $  + 4.9348022005446793d+00*HR3(-1,-1,1)
     $  + 4.9348022005446793d+00*HR3(0,-1,-1)
     $  + 4.9348022005446793d+00*HR3(0,0,-1)
     $  + 4.9348022005446793d+00*HR3(0,0,1)
     $  + 4.9348022005446793d+00*HR3(0,1,-1)
      Hi5(0,0,0,0,1) = 
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 2.5000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(1)*HR1(1)*HR1(1)
     $  + 4.1666666666666666d-02*HR1(1)*HR1(1)*HR1(1)*HR1(1)
      Hi5(0,0,0,1,1) = 
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
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(1)*HR1(1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(1)*HR1(1)*HR1(1)
     $  - HR1( -1)*HR1(1)*HR2(0,-1)
     $  - HR1( -1)*HR1(1)*HR2(0,1)
     $  + HR1( -1)*HR3(0,-1,-1)
     $  + HR1( -1)*HR3(0,-1,1)
     $  + HR1( -1)*HR3(0,1,-1)
     $  + HR1( -1)*HR3(0,1,1)
     $  + 1.6666666666666666d-01*HR1(0)*HR1(1)*HR1(1)*HR1(1)
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
      Hi5(0,0,1,0,1) = 
     $  + 8.2246703342411321d-01*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 2.5000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(0,1)
     $  + 1.6449340668482264d+00*HR1(-1)*HR1(1)
     $  + HR1( -1)*HR1(1)*HR2(-1,1)
     $  + HR1( -1)*HR1(1)*HR2(0,-1)
     $  + HR1( -1)*HR1(1)*HR2(0,1)
     $  + HR1( -1)*HR3(-1,-1,1)
     $  - HR1( -1)*HR3(-1,1,1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR3(0,-1,-1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR3(0,-1,1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR3(0,1,-1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR3(0,1,1)
     $  + 8.2246703342411321d-01*HR1(1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(1)*HR1(1)*HR2(-1,1)
     $  + 5.0000000000000000d-01*HR1(1)*HR1(1)*HR2(0,-1)
     $  + 5.0000000000000000d-01*HR1(1)*HR1(1)*HR2(0,1)
     $  - 2.0000000000000000d+00*HR1(1)*HR3(-1,-1,1)
     $  + 2.0000000000000000d+00*HR1(1)*HR3(-1,1,1)
     $  - 2.0000000000000000d+00*HR1(1)*HR3(0,-1,-1)
     $  - 2.0000000000000000d+00*HR1(1)*HR3(0,-1,1)
     $  - 2.0000000000000000d+00*HR1(1)*HR3(0,1,-1)
     $  - 2.0000000000000000d+00*HR1(1)*HR3(0,1,1)
     $  - 3.0000000000000000d+00*HR4(-1,-1,-1,1)
     $  + 3.0000000000000000d+00*HR4(-1,-1,1,1)
     $  - 3.0000000000000000d+00*HR4(-1,1,1,1)
     $  + 3.0000000000000000d+00*HR4(0,-1,-1,-1)
     $  + 3.0000000000000000d+00*HR4(0,-1,-1,1)
     $  + 3.0000000000000000d+00*HR4(0,-1,1,-1)
     $  + 3.0000000000000000d+00*HR4(0,-1,1,1)
     $  + 3.0000000000000000d+00*HR4(0,1,-1,-1)
     $  + 3.0000000000000000d+00*HR4(0,1,-1,1)
     $  + 3.0000000000000000d+00*HR4(0,1,1,-1)
     $  + 3.0000000000000000d+00*HR4(0,1,1,1)
      Hi5(0,0,1,1,1) = 
     $  - 7.0235377994456286d-01*HR1(-1)*HR1(-1)
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
     $  - 1.4047075598891257d+00*HR1(-1)*HR1(1)
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
     $  - 7.0235377994456286d-01*HR1(1)*HR1(1)
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
      Hi5(0,1,0,1,1) = 
     $  - 1.2020569031595942d+00*HR1(-1)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(0,1)
     $  + HR1( -1)*HR1(0)*HR2(-1,1)
     $  + HR1( -1)*HR1(0)*HR2(0,-1)
     $  + HR1( -1)*HR1(0)*HR2(0,1)
     $  + HR1( -1)*HR1(1)*HR2(-1,1)
     $  + HR1( -1)*HR1(1)*HR2(0,-1)
     $  + HR1( -1)*HR1(1)*HR2(0,1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR2(-1,1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR2(0,1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR3(-1,-1,1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR3(-1,1,1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR3(0,-1,-1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR3(0,0,-1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR3(0,0,1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR3(0,1,-1)
     $  - HR1(0) *HR1(1)*HR2(-1,1)
     $  + HR1(0) *HR1(1)*HR2(0,-1)
     $  + HR1(0) *HR1(1)*HR2(0,1)
     $  - 2.0000000000000000d+00*HR1(0)*HR3(-1,-1,1)
     $  + 2.0000000000000000d+00*HR1(0)*HR3(-1,1,1)
     $  - 2.0000000000000000d+00*HR1(0)*HR3(0,-1,-1)
     $  - 2.0000000000000000d+00*HR1(0)*HR3(0,-1,1)
     $  - 2.0000000000000000d+00*HR1(0)*HR3(0,1,-1)
     $  - 2.0000000000000000d+00*HR1(0)*HR3(0,1,1)
     $  - 1.2020569031595942d+00*HR1(1)
     $  - 6.9314718055994530d-01*HR1(1)*HR2(-1,1)
     $  + 6.9314718055994530d-01*HR1(1)*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR1(1)*HR2(0,1)
     $  - 2.0000000000000000d+00*HR1(1)*HR3(-1,-1,1)
     $  - 2.0000000000000000d+00*HR1(1)*HR3(0,-1,-1)
     $  - 2.0000000000000000d+00*HR1(1)*HR3(0,0,-1)
     $  - 2.0000000000000000d+00*HR1(1)*HR3(0,0,1)
     $  - 2.0000000000000000d+00*HR1(1)*HR3(0,1,-1)
     $  + 5.0000000000000000d-01*HR2(-1,1)*HR2(-1,1)
     $  - HR2( -1,1)*HR2(0,-1)
     $  - HR2( -1,1)*HR2(0,1)
     $  + 5.0000000000000000d-01*HR2(0,-1)*HR2(0,-1)
     $  + HR2(0, -1)*HR2(0,1)
     $  + 5.0000000000000000d-01*HR2(0,1)*HR2(0,1)
     $  - 1.3862943611198906d+00*HR3(-1,-1,1)
     $  + 1.3862943611198906d+00*HR3(-1,1,1)
     $  - 1.3862943611198906d+00*HR3(0,-1,-1)
     $  - 1.3862943611198906d+00*HR3(0,-1,1)
     $  - 1.3862943611198906d+00*HR3(0,1,-1)
     $  - 1.3862943611198906d+00*HR3(0,1,1)
     $  - 4.0000000000000000d+00*HR4(-1,-1,-1,1)
     $  + 2.0000000000000000d+00*HR4(-1,-1,1,1)
     $  + 4.0000000000000000d+00*HR4(0,-1,-1,-1)
     $  + 2.0000000000000000d+00*HR4(0,-1,-1,1)
     $  + 2.0000000000000000d+00*HR4(0,-1,1,-1)
     $  + 2.0000000000000000d+00*HR4(0,0,-1,-1)
     $  + 2.0000000000000000d+00*HR4(0,0,-1,1)
     $  + 2.0000000000000000d+00*HR4(0,0,1,-1)
     $  + 2.0000000000000000d+00*HR4(0,0,1,1)
     $  + 4.0000000000000000d+00*HR4(0,1,-1,-1)
     $  + 2.0000000000000000d+00*HR4(0,1,-1,1)
     $  + 2.0000000000000000d+00*HR4(0,1,1,-1)
      Hi5(0,1,1,1,1) = 
     $  - 1.0846773019780311d+00*HR1(-1)
     $  + 7.0235377994456286d-01*HR1(-1)*HR1(-1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  - 2.5000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(0)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  - 1.4047075598891257d+00*HR1(-1)*HR1(0)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(0)*HR1(0)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(0)*HR1(0)*HR1(0)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(0)*HR1(0)*HR1(1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(0)*HR1(1)
     $  - HR1( -1)*HR1(0)*HR2(-1,1)
     $  + 1.4047075598891257d+00*HR1(-1)*HR1(1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR2(-1,1)
     $  - HR1( -1)*HR3(-1,-1,1)
     $  + 1.6666666666666666d-01*HR1(0)*HR1(0)*HR1(0)*HR1(1)
     $  + 3.4657359027997265d-01*HR1(0)*HR1(0)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(0)*HR1(0)*HR2(-1,1)
     $  - 5.0000000000000000d-01*HR1(0)*HR1(0)*HR2(0,-1)
     $  - 5.0000000000000000d-01*HR1(0)*HR1(0)*HR2(0,1)
     $  - 1.4047075598891257d+00*HR1(0)*HR1(1)
     $  + 6.9314718055994530d-01*HR1(0)*HR2(-1,1)
     $  - 6.9314718055994530d-01*HR1(0)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(0)*HR2(0,1)
     $  + HR1(0) *HR3(-1,-1,1)
     $  + HR1(0) *HR3(0,-1,-1)
     $  + HR1(0) *HR3(0,0,-1)
     $  + HR1(0) *HR3(0,0,1)
     $  + HR1(0) *HR3(0,1,-1)
     $  - 1.0846773019780311d+00*HR1(1)
     $  - 1.4047075598891257d+00*HR2(-1,1)
     $  + 1.4047075598891257d+00*HR2(0,-1)
     $  + 1.4047075598891257d+00*HR2(0,1)
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
      endif 
      endif 
** nw > 4 endif 
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
      if ( nw.gt.4 ) then
      HY5(0,0,0,0,-1) = 
     $  + 9.7211977044690930d-01
     $  - 9.4703282949724591d-01*HR1(-1)
     $  + 4.5077133868484785d-01*HR1(-1)*HR1(-1)
     $  - 1.3707783890401886d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 2.8881132523331054d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 8.3333333333333333d-03*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 8.3333333333333333d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  - 4.1123351671205660d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  + 1.7328679513998632d-01*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  - 8.3333333333333333d-02*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)*HR1(1)
     $  + 9.0154267736969571d-01*HR1(-1)*HR1(1)
     $  - 4.1123351671205660d-01*HR1(-1)*HR1(1)*HR1(1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(1)*HR1(1)*HR1(1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(1)*HR1(1)*HR1(1)*HR1(1)
     $  - 9.4703282949724591d-01*HR1(1)
     $  + 4.5077133868484785d-01*HR1(1)*HR1(1)
     $  - 1.3707783890401886d-01*HR1(1)*HR1(1)*HR1(1)
     $  + 2.8881132523331054d-02*HR1(1)*HR1(1)*HR1(1)*HR1(1)
     $  + 1.6666666666666666d-01*HR1(1)*HR1(1)*HR1(1)*HR2(-1,1)
     $  + 5.0000000000000000d-01*HR1(1)*HR1(1)*HR3(-1,-1,1)
     $  - 5.0000000000000000d-01*HR1(1)*HR1(1)*HR3(-1,1,1)
     $  + HR1(1) *HR4(-1,-1,-1,1)
     $  - HR1(1) *HR4(-1,-1,1,1)
     $  + HR1(1) *HR4(-1,1,1,1)
     $  + HR5( -1,-1,-1,-1,1)
     $  - HR5( -1,-1,-1,1,1)
     $  + HR5( -1,-1,1,1,1)
     $  - HR5( -1,1,1,1,1)
      HY5(0,0,0,-1,-1) = 
     $  + 4.8936397049969063d-02
     $  - 8.7785671568655302d-02*HR1(-1)
     $  + 7.5128556447474642d-02*HR1(-1)*HR1(-1)
     $  - 4.0037751159850118d-02*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 2.8881132523331054d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 8.3333333333333333d-03*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 8.3333333333333333d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  - 1.2011325347955035d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  + 1.7328679513998632d-01*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  - 8.3333333333333333d-02*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)*HR1(1)
     $  + 1.5025711289494928d-01*HR1(-1)*HR1(1)
     $  - 1.2011325347955035d-01*HR1(-1)*HR1(1)*HR1(1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(1)*HR1(1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(1)*HR1(1)*HR2(-1,1)
     $  + HR1( -1)*HR1(1)*HR3(-1,-1,1)
     $  - HR1( -1)*HR1(1)*HR3(-1,1,1)
     $  + HR1( -1)*HR4(-1,-1,-1,1)
     $  - HR1( -1)*HR4(-1,-1,1,1)
     $  + HR1( -1)*HR4(-1,1,1,1)
     $  - 8.7785671568655302d-02*HR1(1)
     $  + 7.5128556447474642d-02*HR1(1)*HR1(1)
     $  - 4.0037751159850118d-02*HR1(1)*HR1(1)*HR1(1)
     $  - 3.4657359027997265d-01*HR1(1)*HR1(1)*HR2(-1,1)
     $  - 5.0000000000000000d-01*HR1(1)*HR1(1)*HR3(-1,-1,1)
     $  - 6.9314718055994530d-01*HR1(1)*HR3(-1,-1,1)
     $  + 6.9314718055994530d-01*HR1(1)*HR3(-1,1,1)
     $  - 2.0000000000000000d+00*HR1(1)*HR4(-1,-1,-1,1)
     $  + HR1(1) *HR4(-1,-1,1,1)
     $  - 6.9314718055994530d-01*HR4(-1,-1,-1,1)
     $  + 6.9314718055994530d-01*HR4(-1,-1,1,1)
     $  - 6.9314718055994530d-01*HR4(-1,1,1,1)
     $  - 3.0000000000000000d+00*HR5(-1,-1,-1,-1,1)
     $  + 2.0000000000000000d+00*HR5(-1,-1,-1,1,1)
     $  - HR5( -1,-1,1,1,1)
      HY5(0,0,-1,0,-1) = 
     $  + 9.2748467341632644d-02
     $  - 1.6265466739742008d-01*HR1(-1)
     $  + 1.3478823976576390d-01*HR1(-1)*HR1(-1)
     $  - 1.3707783890401886d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 2.8881132523331054d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 8.3333333333333333d-03*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 8.3333333333333333d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  - 4.1123351671205660d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  + 1.7328679513998632d-01*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR2(-1,1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(-1,-1,1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(-1,1,1)
     $  + 2.6957647953152780d-01*HR1(-1)*HR1(1)
     $  - 4.1123351671205660d-01*HR1(-1)*HR1(1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(1)*HR1(1)*HR2(-1,1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(1)*HR2(-1,1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR1(1)*HR3(-1,-1,1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR1(1)*HR3(-1,1,1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR3(-1,-1,1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR3(-1,1,1)
     $  - 3.0000000000000000d+00*HR1(-1)*HR4(-1,-1,-1,1)
     $  + 3.0000000000000000d+00*HR1(-1)*HR4(-1,-1,1,1)
     $  - 3.0000000000000000d+00*HR1(-1)*HR4(-1,1,1,1)
     $  - 1.6265466739742008d-01*HR1(1)
     $  + 1.3478823976576390d-01*HR1(1)*HR1(1)
     $  + 3.4657359027997265d-01*HR1(1)*HR1(1)*HR2(-1,1)
     $  + HR1(1) *HR1(1)*HR3(-1,-1,1)
     $  + 8.2246703342411321d-01*HR1(1)*HR2(-1,1)
     $  - 5.0000000000000000d-01*HR1(1)*HR2(-1,1)*HR2(-1,1)
     $  + 1.3862943611198906d+00*HR1(1)*HR3(-1,-1,1)
     $  - 1.3862943611198906d+00*HR1(1)*HR3(-1,1,1)
     $  + 4.0000000000000000d+00*HR1(1)*HR4(-1,-1,-1,1)
     $  - 2.0000000000000000d+00*HR1(1)*HR4(-1,-1,1,1)
     $  + HR2( -1,1)*HR3(-1,1,1)
     $  + 8.2246703342411321d-01*HR3(-1,-1,1)
     $  - 8.2246703342411321d-01*HR3(-1,1,1)
     $  + 2.0794415416798359d+00*HR4(-1,-1,-1,1)
     $  - 2.0794415416798359d+00*HR4(-1,-1,1,1)
     $  + 2.0794415416798359d+00*HR4(-1,1,1,1)
     $  + 7.0000000000000000d+00*HR5(-1,-1,-1,-1,1)
     $  - 7.0000000000000000d+00*HR5(-1,-1,-1,1,1)
     $  - HR5( -1,-1,1,-1,1)
     $  - HR5( -1,1,-1,1,1)
      HY5(0,0,-1,-1,-1) = 
     $  + 9.6015684431298325d-03
     $  - 2.3752366322618485d-02*HR1(-1)
     $  + 2.7752054332410789d-02*HR1(-1)*HR1(-1)
     $  - 4.0037751159850118d-02*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 2.8881132523331054d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 8.3333333333333333d-03*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 8.3333333333333333d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  - 1.2011325347955035d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  + 1.7328679513998632d-01*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR2(-1,1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(-1,-1,1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(-1,1,1)
     $  + 5.5504108664821579d-02*HR1(-1)*HR1(1)
     $  - 1.2011325347955035d-01*HR1(-1)*HR1(1)*HR1(1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(1)*HR2(-1,1)
     $  - HR1( -1)*HR1(1)*HR3(-1,-1,1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR3(-1,-1,1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR3(-1,1,1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR4(-1,-1,-1,1)
     $  + HR1( -1)*HR4(-1,-1,1,1)
     $  - 2.3752366322618485d-02*HR1(1)
     $  + 2.7752054332410789d-02*HR1(1)*HR1(1)
     $  + 2.4022650695910071d-01*HR1(1)*HR2(-1,1)
     $  + 6.9314718055994530d-01*HR1(1)*HR3(-1,-1,1)
     $  + HR1(1) *HR4(-1,-1,-1,1)
     $  + 2.4022650695910071d-01*HR3(-1,-1,1)
     $  - 2.4022650695910071d-01*HR3(-1,1,1)
     $  + 1.3862943611198906d+00*HR4(-1,-1,-1,1)
     $  - 6.9314718055994530d-01*HR4(-1,-1,1,1)
     $  + 3.0000000000000000d+00*HR5(-1,-1,-1,-1,1)
     $  - HR5( -1,-1,-1,1,1)
      HY5(0,-1,0,-1,-1) = 
     $  + 1.3531263989594243d-02
     $  - 3.2893195194356041d-02*HR1(-1)
     $  + 7.5128556447474642d-02*HR1(-1)*HR1(-1)
     $  - 4.0037751159850118d-02*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 2.8881132523331054d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 8.3333333333333333d-03*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  - 1.2011325347955035d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR2(-1,1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  - HR1( -1)*HR1(-1)*HR3(-1,-1,1)
     $  + HR1( -1)*HR1(-1)*HR3(-1,1,1)
     $  + 1.5025711289494928d-01*HR1(-1)*HR1(1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(1)*HR2(-1,1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR1(1)*HR3(-1,-1,1)
     $  + 2.4022650695910071d-01*HR1(-1)*HR2(-1,1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR2(-1,1)*HR2(-1,1)
     $  + 1.3862943611198906d+00*HR1(-1)*HR3(-1,-1,1)
     $  - 1.3862943611198906d+00*HR1(-1)*HR3(-1,1,1)
     $  + 4.0000000000000000d+00*HR1(-1)*HR4(-1,-1,-1,1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR4(-1,-1,1,1)
     $  - 3.2893195194356041d-02*HR1(1)
     $  - 2.4022650695910071d-01*HR1(1)*HR2(-1,1)
     $  - 1.3862943611198906d+00*HR1(1)*HR3(-1,-1,1)
     $  - 3.0000000000000000d+00*HR1(1)*HR4(-1,-1,-1,1)
     $  - 1.5025711289494928d-01*HR2(-1,1)
     $  + 3.4657359027997265d-01*HR2(-1,1)*HR2(-1,1)
     $  + HR2( -1,1)*HR3(-1,-1,1)
     $  - 4.8045301391820142d-01*HR3(-1,-1,1)
     $  + 4.8045301391820142d-01*HR3(-1,1,1)
     $  - 2.7725887222397812d+00*HR4(-1,-1,-1,1)
     $  + 1.3862943611198906d+00*HR4(-1,-1,1,1)
     $  - 7.0000000000000000d+00*HR5(-1,-1,-1,-1,1)
     $  - HR5( -1,-1,1,-1,1)
       HY5(0,-1,-1,-1,-1) = 
     $  + 3.1350096016808622d-03
     $  - 9.6181291076284771d-03*HR1(-1)
     $  + 2.7752054332410789d-02*HR1(-1)*HR1(-1)
     $  - 4.0037751159850118d-02*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 2.8881132523331054d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 8.3333333333333333d-03*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  - 1.2011325347955035d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(-1,-1,1)
     $  + 5.5504108664821579d-02*HR1(-1)*HR1(1)
     $  + 2.4022650695910071d-01*HR1(-1)*HR2(-1,1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR3(-1,-1,1)
     $  + HR1( -1)*HR4(-1,-1,-1,1)
     $  - 9.6181291076284771d-03*HR1(1)
     $  - 5.5504108664821579d-02*HR2(-1,1)
     $  - 2.4022650695910071d-01*HR3(-1,-1,1)
     $  - 6.9314718055994530d-01*HR4(-1,-1,-1,1)
     $  - HR5( -1,-1,-1,-1,1)
      endif 
** nw > 4 endif 
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
      if ( nw.gt.4 ) then 
      HY5(-1,-1,-1,-1,1) = 
     $  + 1.8016537870380179d-03
     $  - 1.4134237214990008d-02*HR1(-1)
     $  + 4.7376502115063852d-02*HR1(-1)*HR1(-1)
     $  - 9.7040087744168750d-02*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 2.8881132523331054d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 8.3333333333333333d-03*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(0,-1,-1)
     $  + HR1( -1)*HR4(0,-1,-1,-1)
     $  - HR5(0, -1,-1,-1,-1)
      HY5(-1,-1,-1,1,1) = 
     $  + 3.8760673146652637d-03
     $  - 4.0758239159309251d-02*HR1(-1)
     $  + 2.6860659680402010d-01*HR1(-1)*HR1(-1)
     $  - 4.0037751159850118d-02*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 2.8881132523331054d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 8.3333333333333333d-03*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  - 8.3333333333333333d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)*HR1(0)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR2(0,-1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(0,-1,-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(0,0,-1)
     $  - HR1( -1)*HR1(0)*HR3(0,-1,-1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR3(0,-1,-1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR4(0,-1,-1,-1)
     $  + HR1( -1)*HR4(0,0,-1,-1)
     $  + HR1(0) *HR4(0,-1,-1,-1)
     $  + 6.9314718055994530d-01*HR4(0,-1,-1,-1)
     $  - 3.0000000000000000d+00*HR5(0,-1,-1,-1,-1)
     $  - HR5(0,0, -1,-1,-1)
      HY5(-1,-1,1,-1,1) = 
     $  + 6.2154684604081354d-03
     $  - 8.7985537010508960d-02*HR1(-1)
     $  - 7.3900238327152102d-01*HR1(-1)*HR1(-1)
     $  + 9.7040087744168750d-02*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 2.8881132523331054d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 8.3333333333333333d-03*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  - 2.9112026323250625d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR2(0,-1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  + HR1( -1)*HR1(-1)*HR3(0,-1,-1)
     $  + HR1( -1)*HR1(-1)*HR3(0,0,-1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR1(0)*HR3(0,-1,-1)
     $  + 5.8224052646501250d-01*HR1(-1)*HR2(0,-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR2(0,-1)*HR2(0,-1)
     $  + 1.3862943611198906d+00*HR1(-1)*HR3(0,-1,-1)
     $  - 4.0000000000000000d+00*HR1(-1)*HR4(0,-1,-1,-1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR4(0,0,-1,-1)
     $  - 3.0000000000000000d+00*HR1(0)*HR4(0,-1,-1,-1)
     $  + HR2(0, -1)*HR3(0,-1,-1)
     $  - 5.8224052646501250d-01*HR3(0,-1,-1)
     $  - 2.0794415416798359d+00*HR4(0,-1,-1,-1)
     $  + 7.0000000000000000d+00*HR5(0,-1,-1,-1,-1)
     $  - HR5(0, -1,0,-1,-1)
      HY5(-1,-1,1,1,1) = 
     $  + 1.8530786065466613d-02
     $  - 5.1747906167389938d-01*HR1(-1)
     $  - 2.7752054332410789d-02*HR1(-1)*HR1(-1)
     $  + 4.0037751159850118d-02*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 2.8881132523331054d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 8.3333333333333333d-03*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  + 8.3333333333333333d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)*HR1(0)
     $  - 1.2011325347955035d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  - 1.7328679513998632d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(0)
     $  - 8.3333333333333333d-02*HR1(-1)*HR1(-1)*HR1(0)*HR1(0)*HR1(0)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(0)*HR1(0)*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(0)*HR2(0,-1)
     $  - HR1( -1)*HR1(0)*HR3(0,-1,-1)
     $  - HR1( -1)*HR1(0)*HR3(0,0,-1)
     $  + 2.4022650695910071d-01*HR1(-1)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR3(0,-1,-1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR3(0,0,-1)
     $  + HR1( -1)*HR4(0,-1,-1,-1)
     $  + HR1( -1)*HR4(0,0,-1,-1)
     $  + HR1( -1)*HR4(0,0,0,-1)
     $  - 5.0000000000000000d-01*HR1(0)*HR1(0)*HR3(0,-1,-1)
     $  - 6.9314718055994530d-01*HR1(0)*HR3(0,-1,-1)
     $  + 2.0000000000000000d+00*HR1(0)*HR4(0,-1,-1,-1)
     $  + HR1(0) *HR4(0,0,-1,-1)
     $  - 2.4022650695910071d-01*HR3(0,-1,-1)
     $  + 1.3862943611198906d+00*HR4(0,-1,-1,-1)
     $  + 6.9314718055994530d-01*HR4(0,0,-1,-1)
     $  - 3.0000000000000000d+00*HR5(0,-1,-1,-1,-1)
     $  - 2.0000000000000000d+00*HR5(0,0,-1,-1,-1)
     $  - HR5(0,0,0, -1,-1)
      HY5(-1,1,-1,1,1) = 
     $  + 3.8880058841843904d-02
     $  + 1.9248049955307152d+00*HR1(-1)
     $  - 2.6860659680402010d-01*HR1(-1)*HR1(-1)
     $  + 4.0037751159850118d-02*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 2.8881132523331054d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 8.3333333333333333d-03*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  + 8.3333333333333333d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)*HR1(0)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR2(0,-1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(0,-1,-1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(0,0,-1)
     $  + 5.3721319360804020d-01*HR1(-1)*HR1(0)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(0)*HR1(0)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(0)*HR2(0,-1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR1(0)*HR3(0,-1,-1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR1(0)*HR3(0,0,-1)
     $  - 2.4022650695910071d-01*HR1(-1)*HR2(0,-1)
     $  + 1.3862943611198906d+00*HR1(-1)*HR3(0,-1,-1)
     $  + 1.3862943611198906d+00*HR1(-1)*HR3(0,0,-1)
     $  - 3.0000000000000000d+00*HR1(-1)*HR4(0,-1,-1,-1)
     $  - 3.0000000000000000d+00*HR1(-1)*HR4(0,0,-1,-1)
     $  - 3.0000000000000000d+00*HR1(-1)*HR4(0,0,0,-1)
     $  + HR1(0) *HR1(0)*HR3(0,-1,-1)
     $  - 5.0000000000000000d-01*HR1(0)*HR2(0,-1)*HR2(0,-1)
     $  + 1.3862943611198906d+00*HR1(0)*HR3(0,-1,-1)
     $  - 4.0000000000000000d+00*HR1(0)*HR4(0,-1,-1,-1)
     $  - 2.0000000000000000d+00*HR1(0)*HR4(0,0,-1,-1)
     $  - 5.3721319360804020d-01*HR2(0,-1)
     $  - 3.4657359027997265d-01*HR2(0,-1)*HR2(0,-1)
     $  + HR2(0, -1)*HR3(0,0,-1)
     $  + 4.8045301391820142d-01*HR3(0,-1,-1)
     $  - 2.7725887222397812d+00*HR4(0,-1,-1,-1)
     $  - 1.3862943611198906d+00*HR4(0,0,-1,-1)
     $  + 7.0000000000000000d+00*HR5(0,-1,-1,-1,-1)
     $  + HR5(0, -1,0,-1,-1)
     $  + 7.0000000000000000d+00*HR5(0,0,-1,-1,-1)
     $  - HR5(0,0, -1,0,-1)
      HY5(-1,1,1,1,1) = 
     $  + 5.0840057924226870d-01
     $  - 9.6181291076284771d-03*HR1(-1)
     $  + 2.7752054332410789d-02*HR1(-1)*HR1(-1)
     $  - 4.0037751159850118d-02*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 2.8881132523331054d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 8.3333333333333333d-03*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  - 8.3333333333333333d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)*HR1(0)
     $  + 1.2011325347955035d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  + 1.7328679513998632d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(0)
     $  + 8.3333333333333333d-02*HR1(-1)*HR1(-1)*HR1(0)*HR1(0)*HR1(0)
     $  - 5.5504108664821579d-02*HR1(-1)*HR1(0)
     $  - 1.2011325347955035d-01*HR1(-1)*HR1(0)*HR1(0)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(0)*HR1(0)*HR1(0)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(0)*HR1(0)*HR1(0)*HR1(0)
     $  + 1.6666666666666666d-01*HR1(0)*HR1(0)*HR1(0)*HR2(0,-1)
     $  + 3.4657359027997265d-01*HR1(0)*HR1(0)*HR2(0,-1)
     $  - 5.0000000000000000d-01*HR1(0)*HR1(0)*HR3(0,-1,-1)
     $  - 5.0000000000000000d-01*HR1(0)*HR1(0)*HR3(0,0,-1)
     $  + 2.4022650695910071d-01*HR1(0)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(0)*HR3(0,-1,-1)
     $  - 6.9314718055994530d-01*HR1(0)*HR3(0,0,-1)
     $  + HR1(0) *HR4(0,-1,-1,-1)
     $  + HR1(0) *HR4(0,0,-1,-1)
     $  + HR1(0) *HR4(0,0,0,-1)
     $  + 5.5504108664821579d-02*HR2(0,-1)
     $  - 2.4022650695910071d-01*HR3(0,-1,-1)
     $  - 2.4022650695910071d-01*HR3(0,0,-1)
     $  + 6.9314718055994530d-01*HR4(0,-1,-1,-1)
     $  + 6.9314718055994530d-01*HR4(0,0,-1,-1)
     $  + 6.9314718055994530d-01*HR4(0,0,0,-1)
     $  - HR5(0, -1,-1,-1,-1)
     $  - HR5(0,0, -1,-1,-1)
     $  - HR5(0,0,0, -1,-1)
     $  - HR5(0,0,0,0, -1)
      HY5(0,-1,-1,-1,1) = 
     $  + 4.1914400448554060d-03
     $  - 1.4134237214990008d-02*HR1(-1)
     $  + 4.7376502115063852d-02*HR1(-1)*HR1(-1)
     $  - 9.7040087744168750d-02*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 2.8881132523331054d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 8.3333333333333333d-03*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR2(-1,1)
     $  - 2.9112026323250625d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR2(0,-1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(-1,-1,1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(0,-1,-1)
     $  - HR1( -1)*HR1(0)*HR3(-1,-1,1)
     $  + 9.4753004230127705d-02*HR1(-1)*HR1(1)
     $  - HR1( -1)*HR1(1)*HR3(0,-1,-1)
     $  + 5.8224052646501250d-01*HR1(-1)*HR2(-1,1)
     $  - HR1( -1)*HR2(-1,1)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR3(-1,-1,1)
     $  - HR1( -1)*HR4(-1,-1,-1,1)
     $  + HR1( -1)*HR4(0,-1,-1,-1)
     $  + HR1(0) *HR4(-1,-1,-1,1)
     $  - 1.4134237214990008d-02*HR1(1)
     $  + HR1(1) *HR4(0,-1,-1,-1)
     $  - 9.4753004230127705d-02*HR2(-1,1)
     $  + HR2( -1,1)*HR3(0,-1,-1)
     $  + HR2(0, -1)*HR3(-1,-1,1)
     $  - 5.8224052646501250d-01*HR3(-1,-1,1)
     $  + 6.9314718055994530d-01*HR4(-1,-1,-1,1)
     $  + HR5( -1,-1,-1,-1,1)
     $  - HR5(0, -1,-1,-1,-1)
     $  - HR5(0, -1,-1,-1,1)
      HY5(0,-1,-1,0,1) = 
     $  + 3.0172237496701167d-02
     $  - 7.7340900566758219d-02*HR1(-1)
     $  + 1.9444792308405316d-01*HR1(-1)*HR1(-1)
     $  - 2.7415567780803773d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 2.8881132523331054d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 8.3333333333333333d-03*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR2(0,1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR2(-1,1)
     $  - 8.2246703342411321d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR2(0,-1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR2(0,1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(-1,-1,1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(0,-1,-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(0,1,-1)
     $  - HR1( -1)*HR1(0)*HR3(-1,-1,1)
     $  + 3.8889584616810632d-01*HR1(-1)*HR1(1)
     $  + HR1( -1)*HR1(1)*HR3(-1,-1,1)
     $  - HR1( -1)*HR1(1)*HR3(0,-1,-1)
     $  - HR1( -1)*HR1(1)*HR3(0,1,-1)
     $  + 1.6449340668482264d+00*HR1(-1)*HR2(-1,1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR2(-1,1)*HR2(-1,1)
     $  - HR1( -1)*HR2(-1,1)*HR2(0,-1)
     $  - HR1( -1)*HR2(-1,1)*HR2(0,1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR3(-1,-1,1)
     $  + HR1( -1)*HR4(0,-1,-1,-1)
     $  + HR1( -1)*HR4(0,1,-1,-1)
     $  - HR1(0) *HR1(1)*HR3(-1,-1,1)
     $  + 5.0000000000000000d-01*HR1(0)*HR2(-1,1)*HR2(-1,1)
     $  - 7.7340900566758219d-02*HR1(1)
     $  - 6.9314718055994530d-01*HR1(1)*HR3(-1,-1,1)
     $  - 3.0000000000000000d+00*HR1(1)*HR4(-1,-1,-1,1)
     $  + HR1(1) *HR4(0,-1,-1,-1)
     $  + HR1(1) *HR4(0,1,-1,-1)
     $  - 3.8889584616810632d-01*HR2(-1,1)
     $  + 3.4657359027997265d-01*HR2(-1,1)*HR2(-1,1)
     $  + 2.0000000000000000d+00*HR2(-1,1)*HR3(-1,-1,1)
     $  + HR2( -1,1)*HR3(0,-1,-1)
     $  + HR2( -1,1)*HR3(0,1,-1)
     $  + HR2(0, -1)*HR3(-1,-1,1)
     $  + HR2(0,1) *HR3(-1,-1,1)
     $  - 1.6449340668482264d+00*HR3(-1,-1,1)
     $  - 3.0000000000000000d+00*HR5(-1,-1,-1,-1,1)
     $  - 6.0000000000000000d+00*HR5(-1,-1,-1,1,1)
     $  - 3.0000000000000000d+00*HR5(-1,-1,1,-1,1)
     $  - HR5(0, -1,-1,-1,-1)
     $  - HR5(0, -1,-1,-1,1)
     $  - HR5(0,1, -1,-1,-1)
     $  - HR5(0,1, -1,-1,1)
      HY5(0,-1,-1,1,-1) = 
     $  + 5.9459097989450212d-03
     $  - 2.3275066086727564d-02*HR1(-1)
     $  + 1.0703618543335311d-01*HR1(-1)*HR1(-1)
     $  + 1.7711559006386898d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 2.8881132523331054d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 8.3333333333333333d-03*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  + 5.3134677019160696d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR2(0,-1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(-1,-1,1)
     $  + HR1( -1)*HR1(-1)*HR3(0,-1,-1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(0)*HR2(-1,1)
     $  + 2.1407237086670622d-01*HR1(-1)*HR1(1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(1)*HR2(0,-1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR1(1)*HR3(0,-1,-1)
     $  - 1.0626935403832139d+00*HR1(-1)*HR2(-1,1)
     $  + HR1( -1)*HR2(-1,1)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR3(-1,-1,1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR3(0,-1,-1)
     $  - HR1( -1)*HR4(-1,-1,-1,1)
     $  - 3.0000000000000000d+00*HR1(-1)*HR4(0,-1,-1,-1)
     $  + 6.9314718055994530d-01*HR1(0)*HR3(-1,-1,1)
     $  - 2.3275066086727564d-02*HR1(1)
     $  + 6.9314718055994530d-01*HR1(1)*HR3(0,-1,-1)
     $  - 3.0000000000000000d+00*HR1(1)*HR4(0,-1,-1,-1)
     $  - 2.1407237086670622d-01*HR2(-1,1)
     $  + 6.9314718055994530d-01*HR2(-1,1)*HR2(0,-1)
     $  - 2.0000000000000000d+00*HR2(-1,1)*HR3(0,-1,-1)
     $  - HR2(0, -1)*HR3(-1,-1,1)
     $  + 1.0626935403832139d+00*HR3(-1,-1,1)
     $  + 6.9314718055994530d-01*HR4(-1,-1,-1,1)
     $  - 6.9314718055994530d-01*HR4(0,-1,-1,-1)
     $  - 6.9314718055994530d-01*HR4(0,-1,-1,1)
     $  + HR5( -1,-1,-1,-1,1)
     $  + 4.0000000000000000d+00*HR5(0,-1,-1,-1,-1)
     $  + 3.0000000000000000d+00*HR5(0,-1,-1,-1,1)
     $  + HR5(0, -1,-1,1,-1)
      HY5(0,-1,-1,1,1) = 
     $  + 8.7734377821481916d-03
     $  - 4.0758239159309251d-02*HR1(-1)
     $  + 2.6860659680402010d-01*HR1(-1)*HR1(-1)
     $  - 4.0037751159850118d-02*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 2.8881132523331054d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 8.3333333333333333d-03*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  - 8.3333333333333333d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)*HR1(0)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  - 2.5000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(0)*HR1(1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR2(-1,1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR2(0,-1)
     $  - 1.2011325347955035d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(-1,-1,1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(0,-1,-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(0,0,-1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(0)*HR1(0)*HR2(-1,1)
     $  + HR1( -1)*HR1(0)*HR1(1)*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(0)*HR2(-1,1)
     $  + HR1( -1)*HR1(0)*HR3(-1,-1,1)
     $  - HR1( -1)*HR1(0)*HR3(0,-1,-1)
     $  + 5.3721319360804020d-01*HR1(-1)*HR1(1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(1)*HR2(0,-1)
     $  - HR1( -1)*HR1(1)*HR3(0,-1,-1)
     $  - HR1( -1)*HR1(1)*HR3(0,0,-1)
     $  + 2.4022650695910071d-01*HR1(-1)*HR2(-1,1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR3(-1,-1,1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR3(0,-1,-1)
     $  + HR1( -1)*HR4(-1,-1,-1,1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR4(0,-1,-1,-1)
     $  + HR1( -1)*HR4(0,0,-1,-1)
     $  - 5.0000000000000000d-01*HR1(0)*HR1(0)*HR3(-1,-1,1)
     $  - HR1(0) *HR1(1)*HR3(0,-1,-1)
     $  - HR1(0) *HR2(-1,1)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(0)*HR3(-1,-1,1)
     $  - HR1(0) *HR4(-1,-1,-1,1)
     $  + HR1(0) *HR4(0,-1,-1,-1)
     $  + HR1(0) *HR4(0,-1,-1,1)
     $  - 4.0758239159309251d-02*HR1(1)
     $  - 6.9314718055994530d-01*HR1(1)*HR3(0,-1,-1)
     $  + 2.0000000000000000d+00*HR1(1)*HR4(0,-1,-1,-1)
     $  + HR1(1) *HR4(0,0,-1,-1)
     $  - 5.3721319360804020d-01*HR2(-1,1)
     $  - 6.9314718055994530d-01*HR2(-1,1)*HR2(0,-1)
     $  + HR2( -1,1)*HR3(0,-1,-1)
     $  + HR2( -1,1)*HR3(0,0,-1)
     $  - 2.4022650695910071d-01*HR3(-1,-1,1)
     $  - 6.9314718055994530d-01*HR4(-1,-1,-1,1)
     $  + 6.9314718055994530d-01*HR4(0,-1,-1,-1)
     $  + 6.9314718055994530d-01*HR4(0,-1,-1,1)
     $  - HR5( -1,-1,-1,-1,1)
     $  - 3.0000000000000000d+00*HR5(0,-1,-1,-1,-1)
     $  - 2.0000000000000000d+00*HR5(0,-1,-1,-1,1)
     $  - HR5(0, -1,-1,1,-1)
     $  - HR5(0,0, -1,-1,-1)
     $  - HR5(0,0, -1,-1,1)
      HY5(0,-1,0,-1,1) = 
     $  + 1.7042475614121991d-02
     $  - 4.5512223866526096d-02*HR1(-1)
     $  + 1.2153517583503078d-01*HR1(-1)*HR1(-1)
     $  - 9.7040087744168750d-02*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 2.8881132523331054d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 8.3333333333333333d-03*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR2(-1,1)
     $  - 2.9112026323250625d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR2(-1,1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR2(0,-1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  + HR1( -1)*HR1(-1)*HR3(-1,-1,1)
     $  - HR1( -1)*HR1(-1)*HR3(-1,1,1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(0,-1,-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(0,-1,1)
     $  - HR1( -1)*HR1(0)*HR1(1)*HR2(-1,1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR1(0)*HR3(-1,-1,1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR1(0)*HR3(-1,1,1)
     $  + 2.4307035167006157d-01*HR1(-1)*HR1(1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(1)*HR2(-1,1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR1(1)*HR3(-1,-1,1)
     $  - HR1( -1)*HR1(1)*HR3(0,-1,-1)
     $  - HR1( -1)*HR1(1)*HR3(0,-1,1)
     $  + 5.8224052646501250d-01*HR1(-1)*HR2(-1,1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR2(-1,1)*HR2(-1,1)
     $  - HR1( -1)*HR2(-1,1)*HR2(0,-1)
     $  - 1.3862943611198906d+00*HR1(-1)*HR3(-1,-1,1)
     $  + 1.3862943611198906d+00*HR1(-1)*HR3(-1,1,1)
     $  - 4.0000000000000000d+00*HR1(-1)*HR4(-1,-1,-1,1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR4(-1,-1,1,1)
     $  + HR1( -1)*HR4(0,-1,-1,-1)
     $  + HR1( -1)*HR4(0,-1,1,-1)
     $  + 2.0000000000000000d+00*HR1(0)*HR1(1)*HR3(-1,-1,1)
     $  - 5.0000000000000000d-01*HR1(0)*HR2(-1,1)*HR2(-1,1)
     $  + 4.0000000000000000d+00*HR1(0)*HR4(-1,-1,-1,1)
     $  - 2.0000000000000000d+00*HR1(0)*HR4(-1,-1,1,1)
     $  - 4.5512223866526096d-02*HR1(1)
     $  - 5.8224052646501250d-01*HR1(1)*HR2(-1,1)
     $  + HR1(1) *HR2(-1,1)*HR2(0,-1)
     $  + 1.3862943611198906d+00*HR1(1)*HR3(-1,-1,1)
     $  + 3.0000000000000000d+00*HR1(1)*HR4(-1,-1,-1,1)
     $  + HR1(1) *HR4(0,-1,-1,-1)
     $  + HR1(1) *HR4(0,-1,1,-1)
     $  - 2.4307035167006157d-01*HR2(-1,1)
     $  - 3.4657359027997265d-01*HR2(-1,1)*HR2(-1,1)
     $  - HR2( -1,1)*HR3(-1,-1,1)
     $  + HR2( -1,1)*HR3(0,-1,-1)
     $  + HR2( -1,1)*HR3(0,-1,1)
     $  + 2.0000000000000000d+00*HR2(0,-1)*HR3(-1,-1,1)
     $  - 2.0000000000000000d+00*HR2(0,-1)*HR3(-1,1,1)
     $  - 1.1644810529300250d+00*HR3(-1,-1,1)
     $  + 1.1644810529300250d+00*HR3(-1,1,1)
     $  + 2.7725887222397812d+00*HR4(-1,-1,-1,1)
     $  - 1.3862943611198906d+00*HR4(-1,-1,1,1)
     $  + 7.0000000000000000d+00*HR5(-1,-1,-1,-1,1)
     $  + HR5( -1,-1,1,-1,1)
     $  - HR5(0, -1,-1,-1,-1)
     $  - HR5(0, -1,-1,-1,1)
     $  - HR5(0, -1,1,-1,-1)
     $  - HR5(0, -1,1,-1,1)
      HY5(0,-1,0,1,-1) = 
     $  + 2.2495758621687517d-02
     $  - 6.9368034302854577d-02*HR1(-1)
     $  + 2.5410760640234242d-01*HR1(-1)*HR1(-1)
     $  + 1.7711559006386898d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 2.8881132523331054d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 8.3333333333333333d-03*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  + 5.3134677019160696d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR2(-1,1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR2(0,-1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR2(0,1)
     $  + HR1( -1)*HR1(-1)*HR3(-1,-1,1)
     $  - HR1( -1)*HR1(-1)*HR3(-1,1,1)
     $  + HR1( -1)*HR1(-1)*HR3(0,-1,-1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(0,-1,1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(0,1,-1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(0)*HR2(-1,1)
     $  + 5.0821521280468485d-01*HR1(-1)*HR1(1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(1)*HR2(-1,1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(1)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(1)*HR2(0,1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR1(1)*HR3(-1,-1,1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR1(1)*HR3(0,-1,-1)
     $  + HR1( -1)*HR1(1)*HR3(0,-1,1)
     $  + HR1( -1)*HR1(1)*HR3(0,1,-1)
     $  - 1.0626935403832139d+00*HR1(-1)*HR2(-1,1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR2(-1,1)*HR2(-1,1)
     $  + HR1( -1)*HR2(-1,1)*HR2(0,-1)
     $  - 1.3862943611198906d+00*HR1(-1)*HR3(-1,-1,1)
     $  + 1.3862943611198906d+00*HR1(-1)*HR3(-1,1,1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR3(0,-1,-1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR3(0,1,-1)
     $  - 4.0000000000000000d+00*HR1(-1)*HR4(-1,-1,-1,1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR4(-1,-1,1,1)
     $  - 3.0000000000000000d+00*HR1(-1)*HR4(0,-1,-1,-1)
     $  - HR1( -1)*HR4(0,-1,1,-1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR4(0,1,-1,-1)
     $  + 6.9314718055994530d-01*HR1(0)*HR1(1)*HR2(-1,1)
     $  + 1.3862943611198906d+00*HR1(0)*HR3(-1,-1,1)
     $  - 1.3862943611198906d+00*HR1(0)*HR3(-1,1,1)
     $  - 6.9368034302854577d-02*HR1(1)
     $  + 1.0626935403832139d+00*HR1(1)*HR2(-1,1)
     $  - HR1(1) *HR2(-1,1)*HR2(0,-1)
     $  + 1.3862943611198906d+00*HR1(1)*HR3(-1,-1,1)
     $  + 6.9314718055994530d-01*HR1(1)*HR3(0,-1,-1)
     $  + 6.9314718055994530d-01*HR1(1)*HR3(0,1,-1)
     $  + 3.0000000000000000d+00*HR1(1)*HR4(-1,-1,-1,1)
     $  - 3.0000000000000000d+00*HR1(1)*HR4(0,-1,-1,-1)
     $  - HR1(1) *HR4(0,-1,1,-1)
     $  - 2.0000000000000000d+00*HR1(1)*HR4(0,1,-1,-1)
     $  - 5.0821521280468485d-01*HR2(-1,1)
     $  - 3.4657359027997265d-01*HR2(-1,1)*HR2(-1,1)
     $  + 6.9314718055994530d-01*HR2(-1,1)*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR2(-1,1)*HR2(0,1)
     $  - HR2( -1,1)*HR3(-1,-1,1)
     $  - 2.0000000000000000d+00*HR2(-1,1)*HR3(0,-1,-1)
     $  - HR2( -1,1)*HR3(0,-1,1)
     $  - HR2( -1,1)*HR3(0,1,-1)
     $  - 2.0000000000000000d+00*HR2(0,-1)*HR3(-1,-1,1)
     $  + 2.0000000000000000d+00*HR2(0,-1)*HR3(-1,1,1)
     $  + 2.1253870807664278d+00*HR3(-1,-1,1)
     $  - 2.1253870807664278d+00*HR3(-1,1,1)
     $  + 2.7725887222397812d+00*HR4(-1,-1,-1,1)
     $  - 1.3862943611198906d+00*HR4(-1,-1,1,1)
     $  - 6.9314718055994530d-01*HR4(0,-1,-1,-1)
     $  - 6.9314718055994530d-01*HR4(0,-1,-1,1)
     $  - 6.9314718055994530d-01*HR4(0,1,-1,-1)
     $  - 6.9314718055994530d-01*HR4(0,1,-1,1)
     $  + 7.0000000000000000d+00*HR5(-1,-1,-1,-1,1)
     $  + HR5( -1,-1,1,-1,1)
     $  + 4.0000000000000000d+00*HR5(0,-1,-1,-1,-1)
     $  + 3.0000000000000000d+00*HR5(0,-1,-1,-1,1)
     $  + HR5(0, -1,-1,1,-1)
     $  + HR5(0, -1,1,-1,-1)
     $  + HR5(0, -1,1,-1,1)
     $  + 3.0000000000000000d+00*HR5(0,1,-1,-1,-1)
     $  + 2.0000000000000000d+00*HR5(0,1,-1,-1,1)
     $  + HR5(0,1, -1,1,-1)
      HY5(0,-1,0,1,1) = 
     $  + 3.0833054551948363d-02
     $  - 1.1285749644390297d-01*HR1(-1)
     $  + 6.0102845157979714d-01*HR1(-1)*HR1(-1)
     $  - 4.0037751159850118d-02*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 2.8881132523331054d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 8.3333333333333333d-03*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  - 8.3333333333333333d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)*HR1(0)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  - 2.5000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(0)*HR1(1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR2(-1,1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR2(0,-1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR2(0,1)
     $  - 1.2011325347955035d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR2(-1,1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR2(0,1)
     $  - HR1( -1)*HR1(-1)*HR3(-1,-1,1)
     $  + HR1( -1)*HR1(-1)*HR3(-1,1,1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(0,-1,-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(0,0,-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(0,0,1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(0,1,-1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(0)*HR1(0)*HR2(-1,1)
     $  + HR1( -1)*HR1(0)*HR1(1)*HR2(-1,1)
     $  + HR1( -1)*HR1(0)*HR1(1)*HR2(0,-1)
     $  + HR1( -1)*HR1(0)*HR1(1)*HR2(0,1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(0)*HR2(-1,1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR1(0)*HR3(-1,-1,1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR1(0)*HR3(-1,1,1)
     $  - HR1( -1)*HR1(0)*HR3(0,-1,-1)
     $  - HR1( -1)*HR1(0)*HR3(0,1,-1)
     $  + 1.2020569031595942d+00*HR1(-1)*HR1(1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(1)*HR2(-1,1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(1)*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(1)*HR2(0,1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR1(1)*HR3(-1,-1,1)
     $  - HR1( -1)*HR1(1)*HR3(0,-1,-1)
     $  - HR1( -1)*HR1(1)*HR3(0,0,-1)
     $  - HR1( -1)*HR1(1)*HR3(0,0,1)
     $  - HR1( -1)*HR1(1)*HR3(0,1,-1)
     $  + 2.4022650695910071d-01*HR1(-1)*HR2(-1,1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR2(-1,1)*HR2(-1,1)
     $  + 1.3862943611198906d+00*HR1(-1)*HR3(-1,-1,1)
     $  - 1.3862943611198906d+00*HR1(-1)*HR3(-1,1,1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR3(0,-1,-1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR3(0,1,-1)
     $  + 4.0000000000000000d+00*HR1(-1)*HR4(-1,-1,-1,1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR4(-1,-1,1,1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR4(0,-1,-1,-1)
     $  + HR1( -1)*HR4(0,0,-1,-1)
     $  + HR1( -1)*HR4(0,0,1,-1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR4(0,1,-1,-1)
     $  - 5.0000000000000000d-01*HR1(0)*HR1(0)*HR1(1)*HR2(-1,1)
     $  - HR1(0) *HR1(0)*HR3(-1,-1,1)
     $  + HR1(0) *HR1(0)*HR3(-1,1,1)
     $  - 6.9314718055994530d-01*HR1(0)*HR1(1)*HR2(-1,1)
     $  - 2.0000000000000000d+00*HR1(0)*HR1(1)*HR3(-1,-1,1)
     $  - HR1(0) *HR1(1)*HR3(0,-1,-1)
     $  - HR1(0) *HR1(1)*HR3(0,1,-1)
     $  + 5.0000000000000000d-01*HR1(0)*HR2(-1,1)*HR2(-1,1)
     $  - HR1(0) *HR2(-1,1)*HR2(0,-1)
     $  - HR1(0) *HR2(-1,1)*HR2(0,1)
     $  - 1.3862943611198906d+00*HR1(0)*HR3(-1,-1,1)
     $  + 1.3862943611198906d+00*HR1(0)*HR3(-1,1,1)
     $  - 4.0000000000000000d+00*HR1(0)*HR4(-1,-1,-1,1)
     $  + 2.0000000000000000d+00*HR1(0)*HR4(-1,-1,1,1)
     $  + HR1(0) *HR4(0,-1,-1,-1)
     $  + HR1(0) *HR4(0,-1,-1,1)
     $  + HR1(0) *HR4(0,1,-1,-1)
     $  + HR1(0) *HR4(0,1,-1,1)
     $  - 1.1285749644390297d-01*HR1(1)
     $  - 2.4022650695910071d-01*HR1(1)*HR2(-1,1)
     $  - 1.3862943611198906d+00*HR1(1)*HR3(-1,-1,1)
     $  - 6.9314718055994530d-01*HR1(1)*HR3(0,-1,-1)
     $  - 6.9314718055994530d-01*HR1(1)*HR3(0,1,-1)
     $  - 3.0000000000000000d+00*HR1(1)*HR4(-1,-1,-1,1)
     $  + 2.0000000000000000d+00*HR1(1)*HR4(0,-1,-1,-1)
     $  + HR1(1) *HR4(0,0,-1,-1)
     $  + HR1(1) *HR4(0,0,1,-1)
     $  + 2.0000000000000000d+00*HR1(1)*HR4(0,1,-1,-1)
     $  - 1.2020569031595942d+00*HR2(-1,1)
     $  + 3.4657359027997265d-01*HR2(-1,1)*HR2(-1,1)
     $  - 6.9314718055994530d-01*HR2(-1,1)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR2(-1,1)*HR2(0,1)
     $  + HR2( -1,1)*HR3(-1,-1,1)
     $  + HR2( -1,1)*HR3(0,-1,-1)
     $  + HR2( -1,1)*HR3(0,0,-1)
     $  + HR2( -1,1)*HR3(0,0,1)
     $  + HR2( -1,1)*HR3(0,1,-1)
     $  - 4.8045301391820142d-01*HR3(-1,-1,1)
     $  + 4.8045301391820142d-01*HR3(-1,1,1)
     $  - 2.7725887222397812d+00*HR4(-1,-1,-1,1)
     $  + 1.3862943611198906d+00*HR4(-1,-1,1,1)
     $  + 6.9314718055994530d-01*HR4(0,-1,-1,-1)
     $  + 6.9314718055994530d-01*HR4(0,-1,-1,1)
     $  + 6.9314718055994530d-01*HR4(0,1,-1,-1)
     $  + 6.9314718055994530d-01*HR4(0,1,-1,1)
     $  - 7.0000000000000000d+00*HR5(-1,-1,-1,-1,1)
     $  - HR5( -1,-1,1,-1,1)
     $  - 3.0000000000000000d+00*HR5(0,-1,-1,-1,-1)
     $  - 2.0000000000000000d+00*HR5(0,-1,-1,-1,1)
     $  - HR5(0, -1,-1,1,-1)
     $  - HR5(0,0, -1,-1,-1)
     $  - HR5(0,0, -1,-1,1)
     $  - HR5(0,0,1, -1,-1)
     $  - HR5(0,0,1, -1,1)
     $  - 3.0000000000000000d+00*HR5(0,1,-1,-1,-1)
     $  - 2.0000000000000000d+00*HR5(0,1,-1,-1,1)
     $  - HR5(0,1, -1,1,-1)
      HY5(0,-1,1,-1,-1) = 
     $  + 9.4133341974174110d-03
     $  - 5.0916764064292634d-02*HR1(-1)
     $  - 2.3766885054564933d-01*HR1(-1)*HR1(-1)
     $  + 4.0037751159850118d-02*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 2.8881132523331054d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 8.3333333333333333d-03*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  - 1.2011325347955035d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  + 1.2011325347955035d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(-1,-1,1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(0,-1,-1)
     $  - 2.4022650695910071d-01*HR1(-1)*HR1(0)*HR1(1)
     $  - 4.7533770109129867d-01*HR1(-1)*HR1(1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(1)*HR2(0,-1)
     $  - HR1( -1)*HR1(1)*HR3(0,-1,-1)
     $  - 2.4022650695910071d-01*HR1(-1)*HR2(-1,1)
     $  + 2.4022650695910071d-01*HR1(-1)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR3(-1,-1,1)
     $  - 1.3862943611198906d+00*HR1(-1)*HR3(0,-1,-1)
     $  - HR1( -1)*HR4(-1,-1,-1,1)
     $  + 3.0000000000000000d+00*HR1(-1)*HR4(0,-1,-1,-1)
     $  + 2.4022650695910071d-01*HR1(0)*HR2(-1,1)
     $  - 5.0916764064292634d-02*HR1(1)
     $  + 2.4022650695910071d-01*HR1(1)*HR2(0,-1)
     $  - 1.3862943611198906d+00*HR1(1)*HR3(0,-1,-1)
     $  + 3.0000000000000000d+00*HR1(1)*HR4(0,-1,-1,-1)
     $  + 4.7533770109129867d-01*HR2(-1,1)
     $  - 6.9314718055994530d-01*HR2(-1,1)*HR2(0,-1)
     $  + HR2( -1,1)*HR3(0,-1,-1)
     $  + 2.4022650695910071d-01*HR3(-1,-1,1)
     $  - 2.4022650695910071d-01*HR3(0,-1,-1)
     $  - 2.4022650695910071d-01*HR3(0,-1,1)
     $  + 6.9314718055994530d-01*HR4(-1,-1,-1,1)
     $  + 2.0794415416798359d+00*HR4(0,-1,-1,-1)
     $  + 1.3862943611198906d+00*HR4(0,-1,-1,1)
     $  + 6.9314718055994530d-01*HR4(0,-1,1,-1)
     $  + HR5( -1,-1,-1,-1,1)
     $  - 6.0000000000000000d+00*HR5(0,-1,-1,-1,-1)
     $  - 3.0000000000000000d+00*HR5(0,-1,-1,-1,1)
     $  - 2.0000000000000000d+00*HR5(0,-1,-1,1,-1)
     $  - HR5(0, -1,1,-1,-1)
      HY5(0,-1,1,-1,1) = 
     $  + 1.3833955759762555d-02
     $  - 8.7985537010508960d-02*HR1(-1)
     $  - 7.3900238327152102d-01*HR1(-1)*HR1(-1)
     $  + 9.7040087744168750d-02*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 2.8881132523331054d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 8.3333333333333333d-03*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  - 2.9112026323250625d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR2(-1,1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR2(0,-1)
     $  + 2.9112026323250625d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR2(0,-1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(-1,-1,1)
     $  + HR1( -1)*HR1(-1)*HR3(0,-1,-1)
     $  + HR1( -1)*HR1(-1)*HR3(0,0,-1)
     $  - 5.8224052646501250d-01*HR1(-1)*HR1(0)*HR1(1)
     $  - HR1( -1)*HR1(0)*HR1(1)*HR2(0,-1)
     $  + HR1( -1)*HR1(0)*HR3(-1,-1,1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR1(0)*HR3(0,-1,-1)
     $  - 1.4780047665430420d+00*HR1(-1)*HR1(1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(1)*HR2(0,-1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR1(1)*HR3(0,-1,-1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR1(1)*HR3(0,0,-1)
     $  - 5.8224052646501250d-01*HR1(-1)*HR2(-1,1)
     $  + HR1( -1)*HR2(-1,1)*HR2(0,-1)
     $  + 5.8224052646501250d-01*HR1(-1)*HR2(0,-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR2(0,-1)*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR3(-1,-1,1)
     $  + 1.3862943611198906d+00*HR1(-1)*HR3(0,-1,-1)
     $  + HR1( -1)*HR4(-1,-1,-1,1)
     $  - 4.0000000000000000d+00*HR1(-1)*HR4(0,-1,-1,-1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR4(0,0,-1,-1)
     $  + 2.0000000000000000d+00*HR1(0)*HR1(1)*HR3(0,-1,-1)
     $  + 5.8224052646501250d-01*HR1(0)*HR2(-1,1)
     $  + HR1(0) *HR2(-1,1)*HR2(0,-1)
     $  - HR1(0) *HR4(-1,-1,-1,1)
     $  - 3.0000000000000000d+00*HR1(0)*HR4(0,-1,-1,-1)
     $  - 2.0000000000000000d+00*HR1(0)*HR4(0,-1,-1,1)
     $  - HR1(0) *HR4(0,-1,1,-1)
     $  - 8.7985537010508960d-02*HR1(1)
     $  + 5.8224052646501250d-01*HR1(1)*HR2(0,-1)
     $  - 5.0000000000000000d-01*HR1(1)*HR2(0,-1)*HR2(0,-1)
     $  + 1.3862943611198906d+00*HR1(1)*HR3(0,-1,-1)
     $  - 4.0000000000000000d+00*HR1(1)*HR4(0,-1,-1,-1)
     $  - 2.0000000000000000d+00*HR1(1)*HR4(0,0,-1,-1)
     $  + 1.4780047665430420d+00*HR2(-1,1)
     $  + 6.9314718055994530d-01*HR2(-1,1)*HR2(0,-1)
     $  - 2.0000000000000000d+00*HR2(-1,1)*HR3(0,-1,-1)
     $  - 2.0000000000000000d+00*HR2(-1,1)*HR3(0,0,-1)
     $  - HR2(0, -1)*HR3(-1,-1,1)
     $  + HR2(0, -1)*HR3(0,-1,-1)
     $  + HR2(0, -1)*HR3(0,-1,1)
     $  + 5.8224052646501250d-01*HR3(-1,-1,1)
     $  - 5.8224052646501250d-01*HR3(0,-1,-1)
     $  - 5.8224052646501250d-01*HR3(0,-1,1)
     $  - 6.9314718055994530d-01*HR4(-1,-1,-1,1)
     $  - 2.0794415416798359d+00*HR4(0,-1,-1,-1)
     $  - 1.3862943611198906d+00*HR4(0,-1,-1,1)
     $  - 6.9314718055994530d-01*HR4(0,-1,1,-1)
     $  - HR5( -1,-1,-1,-1,1)
     $  + 7.0000000000000000d+00*HR5(0,-1,-1,-1,-1)
     $  + 4.0000000000000000d+00*HR5(0,-1,-1,-1,1)
     $  + 2.0000000000000000d+00*HR5(0,-1,-1,1,-1)
     $  - HR5(0, -1,0,-1,-1)
     $  - HR5(0, -1,0,-1,1)
     $  + HR5(0, -1,1,-1,-1)
      HY5(0,-1,1,0,1) = 
     $  + 7.6026642213084631d-02
     $  - 3.5228267839753708d-01*HR1(-1)
     $  - 1.7721476084810206d+00*HR1(-1)*HR1(-1)
     $  + 2.7415567780803773d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 2.8881132523331054d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 8.3333333333333333d-03*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR2(0,1)
     $  - 8.2246703342411321d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR2(-1,1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR2(0,-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR2(0,1)
     $  + 8.2246703342411321d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR2(0,-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR2(0,1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR2(0,1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(-1,-1,1)
     $  + HR1( -1)*HR1(-1)*HR3(0,-1,-1)
     $  + HR1( -1)*HR1(-1)*HR3(0,0,-1)
     $  + HR1( -1)*HR1(-1)*HR3(0,0,1)
     $  + HR1( -1)*HR1(-1)*HR3(0,1,-1)
     $  - 1.6449340668482264d+00*HR1(-1)*HR1(0)*HR1(1)
     $  - HR1( -1)*HR1(0)*HR1(1)*HR2(0,-1)
     $  - HR1( -1)*HR1(0)*HR1(1)*HR2(0,1)
     $  + HR1( -1)*HR1(0)*HR3(-1,-1,1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR1(0)*HR3(0,-1,-1)
     $  + HR1( -1)*HR1(0)*HR3(0,-1,1)
     $  + HR1( -1)*HR1(0)*HR3(0,1,-1)
     $  - 3.5442952169620413d+00*HR1(-1)*HR1(1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(1)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(1)*HR2(0,1)
     $  - HR1( -1)*HR1(1)*HR3(-1,-1,1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR1(1)*HR3(0,-1,-1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR1(1)*HR3(0,0,-1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR1(1)*HR3(0,0,1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR1(1)*HR3(0,1,-1)
     $  - 1.6449340668482264d+00*HR1(-1)*HR2(-1,1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR2(-1,1)*HR2(-1,1)
     $  + HR1( -1)*HR2(-1,1)*HR2(0,-1)
     $  + HR1( -1)*HR2(-1,1)*HR2(0,1)
     $  + 1.6449340668482264d+00*HR1(-1)*HR2(0,-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR2(0,-1)*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR3(-1,-1,1)
     $  + 1.3862943611198906d+00*HR1(-1)*HR3(0,-1,-1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR3(0,-1,1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR3(0,1,-1)
     $  - 4.0000000000000000d+00*HR1(-1)*HR4(0,-1,-1,-1)
     $  - HR1( -1)*HR4(0,-1,0,1)
     $  - HR1( -1)*HR4(0,-1,1,-1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR4(0,0,-1,-1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR4(0,0,-1,1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR4(0,0,1,-1)
     $  - 3.0000000000000000d+00*HR1(-1)*HR4(0,1,-1,-1)
     $  + HR1(0) *HR1(1)*HR3(-1,-1,1)
     $  + 2.0000000000000000d+00*HR1(0)*HR1(1)*HR3(0,-1,-1)
     $  + HR1(0) *HR1(1)*HR3(0,-1,1)
     $  + HR1(0) *HR1(1)*HR3(0,1,-1)
     $  + 1.6449340668482264d+00*HR1(0)*HR2(-1,1)
     $  - 5.0000000000000000d-01*HR1(0)*HR2(-1,1)*HR2(-1,1)
     $  + HR1(0) *HR2(-1,1)*HR2(0,-1)
     $  + HR1(0) *HR2(-1,1)*HR2(0,1)
     $  - 3.0000000000000000d+00*HR1(0)*HR4(0,-1,-1,-1)
     $  - 3.0000000000000000d+00*HR1(0)*HR4(0,-1,-1,1)
     $  - 2.0000000000000000d+00*HR1(0)*HR4(0,-1,1,-1)
     $  - 2.0000000000000000d+00*HR1(0)*HR4(0,-1,1,1)
     $  - HR1(0) *HR4(0,1,-1,-1)
     $  - HR1(0) *HR4(0,1,-1,1)
     $  - 3.5228267839753708d-01*HR1(1)
     $  + 1.6449340668482264d+00*HR1(1)*HR2(0,-1)
     $  - 5.0000000000000000d-01*HR1(1)*HR2(0,-1)*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR1(1)*HR3(-1,-1,1)
     $  + 1.3862943611198906d+00*HR1(1)*HR3(0,-1,-1)
     $  + 6.9314718055994530d-01*HR1(1)*HR3(0,-1,1)
     $  + 6.9314718055994530d-01*HR1(1)*HR3(0,1,-1)
     $  + 3.0000000000000000d+00*HR1(1)*HR4(-1,-1,-1,1)
     $  - 4.0000000000000000d+00*HR1(1)*HR4(0,-1,-1,-1)
     $  - HR1(1) *HR4(0,-1,0,1)
     $  - HR1(1) *HR4(0,-1,1,-1)
     $  - 2.0000000000000000d+00*HR1(1)*HR4(0,0,-1,-1)
     $  - 2.0000000000000000d+00*HR1(1)*HR4(0,0,-1,1)
     $  - 2.0000000000000000d+00*HR1(1)*HR4(0,0,1,-1)
     $  - 3.0000000000000000d+00*HR1(1)*HR4(0,1,-1,-1)
     $  + 3.5442952169620413d+00*HR2(-1,1)
     $  - 3.4657359027997265d-01*HR2(-1,1)*HR2(-1,1)
     $  + 6.9314718055994530d-01*HR2(-1,1)*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR2(-1,1)*HR2(0,1)
     $  - 2.0000000000000000d+00*HR2(-1,1)*HR3(-1,-1,1)
     $  - 2.0000000000000000d+00*HR2(-1,1)*HR3(0,-1,-1)
     $  - 2.0000000000000000d+00*HR2(-1,1)*HR3(0,0,-1)
     $  - 2.0000000000000000d+00*HR2(-1,1)*HR3(0,0,1)
     $  - 2.0000000000000000d+00*HR2(-1,1)*HR3(0,1,-1)
     $  - HR2(0, -1)*HR3(-1,-1,1)
     $  + HR2(0, -1)*HR3(0,-1,-1)
     $  + HR2(0, -1)*HR3(0,-1,1)
     $  - HR2(0,1) *HR3(-1,-1,1)
     $  + 1.6449340668482264d+00*HR3(-1,-1,1)
     $  - 1.6449340668482264d+00*HR3(0,-1,-1)
     $  - 1.6449340668482264d+00*HR3(0,-1,1)
     $  - 2.0794415416798359d+00*HR4(0,-1,-1,-1)
     $  - 2.0794415416798359d+00*HR4(0,-1,-1,1)
     $  - 1.3862943611198906d+00*HR4(0,-1,1,-1)
     $  - 1.3862943611198906d+00*HR4(0,-1,1,1)
     $  - 6.9314718055994530d-01*HR4(0,1,-1,-1)
     $  - 6.9314718055994530d-01*HR4(0,1,-1,1)
     $  + 3.0000000000000000d+00*HR5(-1,-1,-1,-1,1)
     $  + 6.0000000000000000d+00*HR5(-1,-1,-1,1,1)
     $  + 3.0000000000000000d+00*HR5(-1,-1,1,-1,1)
     $  + 7.0000000000000000d+00*HR5(0,-1,-1,-1,-1)
     $  + 4.0000000000000000d+00*HR5(0,-1,-1,-1,1)
     $  + HR5(0, -1,-1,0,1)
     $  + 3.0000000000000000d+00*HR5(0,-1,-1,1,-1)
     $  - HR5(0, -1,0,-1,-1)
     $  + HR5(0, -1,0,1,-1)
     $  + 2.0000000000000000d+00*HR5(0,-1,0,1,1)
     $  + 3.0000000000000000d+00*HR5(0,-1,1,-1,-1)
     $  + HR5(0, -1,1,-1,1)
     $  + HR5(0, -1,1,0,1)
     $  + 2.0000000000000000d+00*HR5(0,-1,1,1,-1)
     $  + 2.0000000000000000d+00*HR5(0,0,-1,-1,1)
     $  + 2.0000000000000000d+00*HR5(0,0,-1,1,-1)
     $  + 4.0000000000000000d+00*HR5(0,0,-1,1,1)
     $  + 2.0000000000000000d+00*HR5(0,0,1,-1,-1)
     $  + 2.0000000000000000d+00*HR5(0,0,1,-1,1)
     $  + 4.0000000000000000d+00*HR5(0,1,-1,-1,-1)
     $  + 3.0000000000000000d+00*HR5(0,1,-1,-1,1)
     $  + HR5(0,1, -1,1,-1)
      HY5(0,-1,1,1,-1) = 
     $  + 2.2801059128486651d-02
     $  - 2.0286579517988963d-01*HR1(-1)
     $  + 5.5365194946473328d-01*HR1(-1)*HR1(-1)
     $  - 1.7711559006386898d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 2.8881132523331054d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 8.3333333333333333d-03*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  + 5.3134677019160696d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  + 1.7328679513998632d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(0)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  - 5.3134677019160696d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR2(0,-1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(-1,-1,1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(0,-1,-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(0,0,-1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(0)*HR1(0)*HR1(1)
     $  + 1.0626935403832139d+00*HR1(-1)*HR1(0)*HR1(1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(0)*HR2(-1,1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(0)*HR2(0,-1)
     $  + 1.1073038989294665d+00*HR1(-1)*HR1(1)
     $  - HR1( -1)*HR1(1)*HR3(0,-1,-1)
     $  - HR1( -1)*HR1(1)*HR3(0,0,-1)
     $  + 1.0626935403832139d+00*HR1(-1)*HR2(-1,1)
     $  - HR1( -1)*HR2(-1,1)*HR2(0,-1)
     $  - 1.0626935403832139d+00*HR1(-1)*HR2(0,-1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR2(0,-1)*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR3(-1,-1,1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR3(0,-1,-1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR3(0,0,-1)
     $  + HR1( -1)*HR4(-1,-1,-1,1)
     $  - 3.4657359027997265d-01*HR1(0)*HR1(0)*HR2(-1,1)
     $  - 6.9314718055994530d-01*HR1(0)*HR1(1)*HR2(0,-1)
     $  - 1.0626935403832139d+00*HR1(0)*HR2(-1,1)
     $  - 6.9314718055994530d-01*HR1(0)*HR3(-1,-1,1)
     $  + 6.9314718055994530d-01*HR1(0)*HR3(0,-1,-1)
     $  + 6.9314718055994530d-01*HR1(0)*HR3(0,-1,1)
     $  - 2.0286579517988963d-01*HR1(1)
     $  - 1.0626935403832139d+00*HR1(1)*HR2(0,-1)
     $  + 5.0000000000000000d-01*HR1(1)*HR2(0,-1)*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR1(1)*HR3(0,-1,-1)
     $  + 6.9314718055994530d-01*HR1(1)*HR3(0,0,-1)
     $  - 1.1073038989294665d+00*HR2(-1,1)
     $  + HR2( -1,1)*HR3(0,-1,-1)
     $  + HR2( -1,1)*HR3(0,0,-1)
     $  + HR2(0, -1)*HR3(-1,-1,1)
     $  - HR2(0, -1)*HR3(0,-1,-1)
     $  - HR2(0, -1)*HR3(0,-1,1)
     $  - 1.0626935403832139d+00*HR3(-1,-1,1)
     $  + 1.0626935403832139d+00*HR3(0,-1,-1)
     $  + 1.0626935403832139d+00*HR3(0,-1,1)
     $  - 6.9314718055994530d-01*HR4(-1,-1,-1,1)
     $  - 1.3862943611198906d+00*HR4(0,-1,-1,-1)
     $  - 6.9314718055994530d-01*HR4(0,-1,-1,1)
     $  - 6.9314718055994530d-01*HR4(0,-1,1,-1)
     $  - 6.9314718055994530d-01*HR4(0,0,-1,-1)
     $  - 6.9314718055994530d-01*HR4(0,0,-1,1)
     $  - HR5( -1,-1,-1,-1,1)
     $  + 2.0000000000000000d+00*HR5(0,-1,-1,-1,-1)
     $  + HR5(0, -1,-1,1,-1)
     $  + HR5(0, -1,0,-1,-1)
     $  + HR5(0, -1,0,-1,1)
     $  + HR5(0, -1,1,-1,-1)
     $  + 3.0000000000000000d+00*HR5(0,0,-1,-1,-1)
     $  + 2.0000000000000000d+00*HR5(0,0,-1,-1,1)
     $  + HR5(0,0, -1,1,-1)
      HY5(0,-1,1,1,1) = 
     $  + 3.9984858137537496d-02
     $  - 5.1747906167389938d-01*HR1(-1)
     $  - 2.7752054332410789d-02*HR1(-1)*HR1(-1)
     $  + 4.0037751159850118d-02*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 2.8881132523331054d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 8.3333333333333333d-03*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  + 8.3333333333333333d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)*HR1(0)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  - 1.2011325347955035d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  - 1.7328679513998632d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(0)
     $  - 8.3333333333333333d-02*HR1(-1)*HR1(-1)*HR1(0)*HR1(0)*HR1(0)
     $  + 2.5000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(0)*HR1(1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR2(-1,1)
     $  + 1.2011325347955035d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(-1,-1,1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(0)*HR1(0)*HR1(0)*HR1(1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(0)*HR1(0)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(0)*HR1(0)*HR2(-1,1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(0)*HR1(0)*HR2(0,-1)
     $  - 2.4022650695910071d-01*HR1(-1)*HR1(0)*HR1(1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(0)*HR2(-1,1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(0)*HR2(0,-1)
     $  - HR1( -1)*HR1(0)*HR3(-1,-1,1)
     $  - HR1( -1)*HR1(0)*HR3(0,-1,-1)
     $  - HR1( -1)*HR1(0)*HR3(0,0,-1)
     $  - 5.5504108664821579d-02*HR1(-1)*HR1(1)
     $  - 2.4022650695910071d-01*HR1(-1)*HR2(-1,1)
     $  + 2.4022650695910071d-01*HR1(-1)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR3(-1,-1,1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR3(0,-1,-1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR3(0,0,-1)
     $  - HR1( -1)*HR4(-1,-1,-1,1)
     $  + HR1( -1)*HR4(0,-1,-1,-1)
     $  + HR1( -1)*HR4(0,0,-1,-1)
     $  + HR1( -1)*HR4(0,0,0,-1)
     $  + 1.6666666666666666d-01*HR1(0)*HR1(0)*HR1(0)*HR2(-1,1)
     $  + 5.0000000000000000d-01*HR1(0)*HR1(0)*HR1(1)*HR2(0,-1)
     $  + 3.4657359027997265d-01*HR1(0)*HR1(0)*HR2(-1,1)
     $  + 5.0000000000000000d-01*HR1(0)*HR1(0)*HR3(-1,-1,1)
     $  - 5.0000000000000000d-01*HR1(0)*HR1(0)*HR3(0,-1,-1)
     $  - 5.0000000000000000d-01*HR1(0)*HR1(0)*HR3(0,-1,1)
     $  + 6.9314718055994530d-01*HR1(0)*HR1(1)*HR2(0,-1)
     $  - HR1(0) *HR1(1)*HR3(0,-1,-1)
     $  - HR1(0) *HR1(1)*HR3(0,0,-1)
     $  + 2.4022650695910071d-01*HR1(0)*HR2(-1,1)
     $  + 6.9314718055994530d-01*HR1(0)*HR3(-1,-1,1)
     $  - 6.9314718055994530d-01*HR1(0)*HR3(0,-1,-1)
     $  - 6.9314718055994530d-01*HR1(0)*HR3(0,-1,1)
     $  + HR1(0) *HR4(-1,-1,-1,1)
     $  + 2.0000000000000000d+00*HR1(0)*HR4(0,-1,-1,-1)
     $  + HR1(0) *HR4(0,-1,-1,1)
     $  + HR1(0) *HR4(0,-1,1,-1)
     $  + HR1(0) *HR4(0,0,-1,-1)
     $  + HR1(0) *HR4(0,0,-1,1)
     $  - 5.1747906167389938d-01*HR1(1)
     $  + 2.4022650695910071d-01*HR1(1)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(1)*HR3(0,-1,-1)
     $  - 6.9314718055994530d-01*HR1(1)*HR3(0,0,-1)
     $  + HR1(1) *HR4(0,-1,-1,-1)
     $  + HR1(1) *HR4(0,0,-1,-1)
     $  + HR1(1) *HR4(0,0,0,-1)
     $  + 5.5504108664821579d-02*HR2(-1,1)
     $  + 2.4022650695910071d-01*HR3(-1,-1,1)
     $  - 2.4022650695910071d-01*HR3(0,-1,-1)
     $  - 2.4022650695910071d-01*HR3(0,-1,1)
     $  + 6.9314718055994530d-01*HR4(-1,-1,-1,1)
     $  + 1.3862943611198906d+00*HR4(0,-1,-1,-1)
     $  + 6.9314718055994530d-01*HR4(0,-1,-1,1)
     $  + 6.9314718055994530d-01*HR4(0,-1,1,-1)
     $  + 6.9314718055994530d-01*HR4(0,0,-1,-1)
     $  + 6.9314718055994530d-01*HR4(0,0,-1,1)
     $  + HR5( -1,-1,-1,-1,1)
     $  - 3.0000000000000000d+00*HR5(0,-1,-1,-1,-1)
     $  - HR5(0, -1,-1,-1,1)
     $  - HR5(0, -1,-1,1,-1)
     $  - HR5(0, -1,1,-1,-1)
     $  - 2.0000000000000000d+00*HR5(0,0,-1,-1,-1)
     $  - HR5(0,0, -1,-1,1)
     $  - HR5(0,0, -1,1,-1)
     $  - HR5(0,0,0, -1,-1)
     $  - HR5(0,0,0, -1,1)
      HY5(0,0,-1,-1,1) = 
     $  + 1.2444228784499648d-02
     $  - 3.4159126166513913d-02*HR1(-1)
     $  + 4.7376502115063852d-02*HR1(-1)*HR1(-1)
     $  - 9.7040087744168750d-02*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 2.8881132523331054d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 8.3333333333333333d-03*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 8.3333333333333333d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  - 2.5000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)*HR1(1)
     $  - 2.9112026323250625d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - 1.7328679513998632d-01*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR2(-1,1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR2(0,-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(-1,-1,1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(-1,1,1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(0,-1,-1)
     $  + HR1( -1)*HR1(0)*HR1(1)*HR2(-1,1)
     $  + HR1( -1)*HR1(0)*HR3(-1,-1,1)
     $  - HR1( -1)*HR1(0)*HR3(-1,1,1)
     $  + 9.4753004230127705d-02*HR1(-1)*HR1(1)
     $  - 2.9112026323250625d-01*HR1(-1)*HR1(1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(1)*HR1(1)*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(1)*HR2(-1,1)
     $  + HR1( -1)*HR1(1)*HR3(-1,-1,1)
     $  - HR1( -1)*HR1(1)*HR3(0,-1,-1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR3(-1,-1,1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR3(-1,1,1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR4(-1,-1,-1,1)
     $  - HR1( -1)*HR4(-1,-1,1,1)
     $  + HR1( -1)*HR4(0,-1,-1,-1)
     $  + HR1( -1)*HR4(0,-1,-1,1)
     $  - HR1(0) *HR1(1)*HR3(-1,-1,1)
     $  - 2.0000000000000000d+00*HR1(0)*HR4(-1,-1,-1,1)
     $  + HR1(0) *HR4(-1,-1,1,1)
     $  - 3.4159126166513913d-02*HR1(1)
     $  + 4.7376502115063852d-02*HR1(1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(1)*HR1(1)*HR3(0,-1,-1)
     $  + 5.8224052646501250d-01*HR1(1)*HR2(-1,1)
     $  - HR1(1) *HR2(-1,1)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(1)*HR3(-1,-1,1)
     $  - HR1(1) *HR4(-1,-1,-1,1)
     $  + HR1(1) *HR4(0,-1,-1,-1)
     $  + HR1(1) *HR4(0,-1,-1,1)
     $  - HR2(0, -1)*HR3(-1,-1,1)
     $  + HR2(0, -1)*HR3(-1,1,1)
     $  + 5.8224052646501250d-01*HR3(-1,-1,1)
     $  - 5.8224052646501250d-01*HR3(-1,1,1)
     $  - 1.3862943611198906d+00*HR4(-1,-1,-1,1)
     $  + 6.9314718055994530d-01*HR4(-1,-1,1,1)
     $  - 3.0000000000000000d+00*HR5(-1,-1,-1,-1,1)
     $  + HR5( -1,-1,-1,1,1)
     $  - HR5(0, -1,-1,-1,-1)
     $  - HR5(0, -1,-1,-1,1)
     $  - HR5(0, -1,-1,1,-1)
     $  - HR5(0, -1,-1,1,1)
      HY5(0,0,-1,0,1) = 
     $  + 1.0679981350605469d-01
     $  - 2.0293560632083841d-01*HR1(-1)
     $  + 1.9444792308405316d-01*HR1(-1)*HR1(-1)
     $  - 2.7415567780803773d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 2.8881132523331054d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 8.3333333333333333d-03*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 8.3333333333333333d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR2(0,1)
     $  - 2.5000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)*HR1(1)
     $  - 8.2246703342411321d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - 1.7328679513998632d-01*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR2(-1,1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR2(0,-1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR2(0,1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(-1,-1,1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(-1,1,1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(0,-1,-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(0,1,-1)
     $  + HR1( -1)*HR1(0)*HR1(1)*HR2(-1,1)
     $  + HR1( -1)*HR1(0)*HR3(-1,-1,1)
     $  - HR1( -1)*HR1(0)*HR3(-1,1,1)
     $  + 3.8889584616810632d-01*HR1(-1)*HR1(1)
     $  - 8.2246703342411321d-01*HR1(-1)*HR1(1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(1)*HR1(1)*HR2(-1,1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(1)*HR1(1)*HR2(0,-1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(1)*HR1(1)*HR2(0,1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(1)*HR2(-1,1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR1(1)*HR3(-1,-1,1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR1(1)*HR3(-1,1,1)
     $  - HR1( -1)*HR1(1)*HR3(0,-1,-1)
     $  - HR1( -1)*HR1(1)*HR3(0,1,-1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR3(-1,-1,1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR3(-1,1,1)
     $  + 3.0000000000000000d+00*HR1(-1)*HR4(-1,-1,-1,1)
     $  - 3.0000000000000000d+00*HR1(-1)*HR4(-1,-1,1,1)
     $  + 3.0000000000000000d+00*HR1(-1)*HR4(-1,1,1,1)
     $  + HR1( -1)*HR4(0,-1,-1,-1)
     $  + HR1( -1)*HR4(0,-1,-1,1)
     $  + HR1( -1)*HR4(0,1,-1,-1)
     $  + HR1( -1)*HR4(0,1,-1,1)
     $  - 5.0000000000000000d-01*HR1(0)*HR1(1)*HR1(1)*HR2(-1,1)
     $  - 2.0000000000000000d+00*HR1(0)*HR1(1)*HR3(-1,-1,1)
     $  + 2.0000000000000000d+00*HR1(0)*HR1(1)*HR3(-1,1,1)
     $  - 3.0000000000000000d+00*HR1(0)*HR4(-1,-1,-1,1)
     $  + 3.0000000000000000d+00*HR1(0)*HR4(-1,-1,1,1)
     $  - 3.0000000000000000d+00*HR1(0)*HR4(-1,1,1,1)
     $  - 2.0293560632083841d-01*HR1(1)
     $  + 1.9444792308405316d-01*HR1(1)*HR1(1)
     $  - 3.4657359027997265d-01*HR1(1)*HR1(1)*HR2(-1,1)
     $  - HR1(1) *HR1(1)*HR3(-1,-1,1)
     $  - 5.0000000000000000d-01*HR1(1)*HR1(1)*HR3(0,-1,-1)
     $  - 5.0000000000000000d-01*HR1(1)*HR1(1)*HR3(0,1,-1)
     $  + 1.6449340668482264d+00*HR1(1)*HR2(-1,1)
     $  + 5.0000000000000000d-01*HR1(1)*HR2(-1,1)*HR2(-1,1)
     $  - HR1(1) *HR2(-1,1)*HR2(0,-1)
     $  - HR1(1) *HR2(-1,1)*HR2(0,1)
     $  - 1.3862943611198906d+00*HR1(1)*HR3(-1,-1,1)
     $  + 1.3862943611198906d+00*HR1(1)*HR3(-1,1,1)
     $  - 4.0000000000000000d+00*HR1(1)*HR4(-1,-1,-1,1)
     $  + 2.0000000000000000d+00*HR1(1)*HR4(-1,-1,1,1)
     $  + HR1(1) *HR4(0,-1,-1,-1)
     $  + HR1(1) *HR4(0,-1,-1,1)
     $  + HR1(1) *HR4(0,1,-1,-1)
     $  + HR1(1) *HR4(0,1,-1,1)
     $  - HR2( -1,1)*HR3(-1,1,1)
     $  - HR2(0, -1)*HR3(-1,-1,1)
     $  + HR2(0, -1)*HR3(-1,1,1)
     $  - HR2(0,1) *HR3(-1,-1,1)
     $  + HR2(0,1) *HR3(-1,1,1)
     $  + 1.6449340668482264d+00*HR3(-1,-1,1)
     $  - 1.6449340668482264d+00*HR3(-1,1,1)
     $  - 2.0794415416798359d+00*HR4(-1,-1,-1,1)
     $  + 2.0794415416798359d+00*HR4(-1,-1,1,1)
     $  - 2.0794415416798359d+00*HR4(-1,1,1,1)
     $  - 7.0000000000000000d+00*HR5(-1,-1,-1,-1,1)
     $  + 7.0000000000000000d+00*HR5(-1,-1,-1,1,1)
     $  + HR5( -1,-1,1,-1,1)
     $  + HR5( -1,1,-1,1,1)
     $  - HR5(0, -1,-1,-1,-1)
     $  - HR5(0, -1,-1,-1,1)
     $  - HR5(0, -1,-1,1,-1)
     $  - HR5(0, -1,-1,1,1)
     $  - HR5(0,1, -1,-1,-1)
     $  - HR5(0,1, -1,-1,1)
     $  - HR5(0,1, -1,1,-1)
     $  - HR5(0,1, -1,1,1)
      HY5(0,0,-1,1,-1) = 
     $  + 1.6991592326175436d-02
     $  - 5.4653052738263652d-02*HR1(-1)
     $  + 1.0703618543335311d-01*HR1(-1)*HR1(-1)
     $  + 1.7711559006386898d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 2.8881132523331054d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 8.3333333333333333d-03*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 8.3333333333333333d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  + 5.3134677019160696d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - 1.7328679513998632d-01*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR2(-1,1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR2(0,-1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(-1,-1,1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(-1,1,1)
     $  + HR1( -1)*HR1(-1)*HR3(0,-1,-1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(0)*HR1(1)*HR1(1)
     $  + 2.1407237086670622d-01*HR1(-1)*HR1(1)
     $  + 5.3134677019160696d-01*HR1(-1)*HR1(1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(1)*HR1(1)*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(1)*HR2(-1,1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(1)*HR2(0,-1)
     $  + HR1( -1)*HR1(1)*HR3(-1,-1,1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR1(1)*HR3(0,-1,-1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR3(-1,-1,1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR3(-1,1,1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR3(0,-1,-1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR3(0,-1,1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR4(-1,-1,-1,1)
     $  - HR1( -1)*HR4(-1,-1,1,1)
     $  - 3.0000000000000000d+00*HR1(-1)*HR4(0,-1,-1,-1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR4(0,-1,-1,1)
     $  - HR1( -1)*HR4(0,-1,1,-1)
     $  - 6.9314718055994530d-01*HR1(0)*HR1(1)*HR2(-1,1)
     $  - 6.9314718055994530d-01*HR1(0)*HR3(-1,-1,1)
     $  + 6.9314718055994530d-01*HR1(0)*HR3(-1,1,1)
     $  - 5.4653052738263652d-02*HR1(1)
     $  + 1.0703618543335311d-01*HR1(1)*HR1(1)
     $  - 3.4657359027997265d-01*HR1(1)*HR1(1)*HR2(0,-1)
     $  + HR1(1) *HR1(1)*HR3(0,-1,-1)
     $  - 1.0626935403832139d+00*HR1(1)*HR2(-1,1)
     $  + HR1(1) *HR2(-1,1)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(1)*HR3(-1,-1,1)
     $  + 6.9314718055994530d-01*HR1(1)*HR3(0,-1,-1)
     $  + 6.9314718055994530d-01*HR1(1)*HR3(0,-1,1)
     $  - HR1(1) *HR4(-1,-1,-1,1)
     $  - 3.0000000000000000d+00*HR1(1)*HR4(0,-1,-1,-1)
     $  - 2.0000000000000000d+00*HR1(1)*HR4(0,-1,-1,1)
     $  - HR1(1) *HR4(0,-1,1,-1)
     $  + HR2(0, -1)*HR3(-1,-1,1)
     $  - HR2(0, -1)*HR3(-1,1,1)
     $  - 1.0626935403832139d+00*HR3(-1,-1,1)
     $  + 1.0626935403832139d+00*HR3(-1,1,1)
     $  - 1.3862943611198906d+00*HR4(-1,-1,-1,1)
     $  + 6.9314718055994530d-01*HR4(-1,-1,1,1)
     $  - 6.9314718055994530d-01*HR4(0,-1,-1,-1)
     $  - 6.9314718055994530d-01*HR4(0,-1,-1,1)
     $  - 6.9314718055994530d-01*HR4(0,-1,1,-1)
     $  - 6.9314718055994530d-01*HR4(0,-1,1,1)
     $  - 3.0000000000000000d+00*HR5(-1,-1,-1,-1,1)
     $  + HR5( -1,-1,-1,1,1)
     $  + 4.0000000000000000d+00*HR5(0,-1,-1,-1,-1)
     $  + 3.0000000000000000d+00*HR5(0,-1,-1,-1,1)
     $  + 3.0000000000000000d+00*HR5(0,-1,-1,1,-1)
     $  + 2.0000000000000000d+00*HR5(0,-1,-1,1,1)
     $  + 2.0000000000000000d+00*HR5(0,-1,1,-1,-1)
     $  + HR5(0, -1,1,-1,1)
     $  + HR5(0, -1,1,1,-1)
      HY5(0,0,-1,1,1) = 
     $  + 2.4107342184124538d-02
     $  - 9.3097125991768577d-02*HR1(-1)
     $  + 2.6860659680402010d-01*HR1(-1)*HR1(-1)
     $  - 4.0037751159850118d-02*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 2.8881132523331054d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 8.3333333333333333d-03*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  - 8.3333333333333333d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)*HR1(0)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 8.3333333333333333d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  - 2.5000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(0)*HR1(1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  + 2.5000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR2(0,-1)
     $  - 1.2011325347955035d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  + 1.7328679513998632d-01*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR2(-1,1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(-1,-1,1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(-1,1,1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(0,-1,-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(0,0,-1)
     $  - 2.5000000000000000d-01*HR1(-1)*HR1(0)*HR1(0)*HR1(1)*HR1(1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(0)*HR1(1)*HR1(1)
     $  - HR1( -1)*HR1(0)*HR1(1)*HR2(-1,1)
     $  + HR1( -1)*HR1(0)*HR1(1)*HR2(0,-1)
     $  - HR1( -1)*HR1(0)*HR3(-1,-1,1)
     $  + HR1( -1)*HR1(0)*HR3(-1,1,1)
     $  - HR1( -1)*HR1(0)*HR3(0,-1,-1)
     $  - HR1( -1)*HR1(0)*HR3(0,-1,1)
     $  + 5.3721319360804020d-01*HR1(-1)*HR1(1)
     $  - 1.2011325347955035d-01*HR1(-1)*HR1(1)*HR1(1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(1)*HR2(-1,1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(1)*HR2(0,-1)
     $  - HR1( -1)*HR1(1)*HR3(-1,-1,1)
     $  - HR1( -1)*HR1(1)*HR3(0,-1,-1)
     $  - HR1( -1)*HR1(1)*HR3(0,0,-1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR3(-1,-1,1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR3(-1,1,1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR3(0,-1,-1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR3(0,-1,1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR4(-1,-1,-1,1)
     $  + HR1( -1)*HR4(-1,-1,1,1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR4(0,-1,-1,-1)
     $  + HR1( -1)*HR4(0,-1,-1,1)
     $  + HR1( -1)*HR4(0,-1,1,-1)
     $  + HR1( -1)*HR4(0,0,-1,-1)
     $  + HR1( -1)*HR4(0,0,-1,1)
     $  + 5.0000000000000000d-01*HR1(0)*HR1(0)*HR1(1)*HR2(-1,1)
     $  + 5.0000000000000000d-01*HR1(0)*HR1(0)*HR3(-1,-1,1)
     $  - 5.0000000000000000d-01*HR1(0)*HR1(0)*HR3(-1,1,1)
     $  + 5.0000000000000000d-01*HR1(0)*HR1(1)*HR1(1)*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR1(0)*HR1(1)*HR2(-1,1)
     $  + HR1(0) *HR1(1)*HR3(-1,-1,1)
     $  - HR1(0) *HR1(1)*HR3(0,-1,-1)
     $  - HR1(0) *HR1(1)*HR3(0,-1,1)
     $  + 6.9314718055994530d-01*HR1(0)*HR3(-1,-1,1)
     $  - 6.9314718055994530d-01*HR1(0)*HR3(-1,1,1)
     $  + 2.0000000000000000d+00*HR1(0)*HR4(-1,-1,-1,1)
     $  - HR1(0) *HR4(-1,-1,1,1)
     $  + HR1(0) *HR4(0,-1,-1,-1)
     $  + HR1(0) *HR4(0,-1,-1,1)
     $  + HR1(0) *HR4(0,-1,1,-1)
     $  + HR1(0) *HR4(0,-1,1,1)
     $  - 9.3097125991768577d-02*HR1(1)
     $  + 2.6860659680402010d-01*HR1(1)*HR1(1)
     $  + 3.4657359027997265d-01*HR1(1)*HR1(1)*HR2(0,-1)
     $  - 5.0000000000000000d-01*HR1(1)*HR1(1)*HR3(0,-1,-1)
     $  - 5.0000000000000000d-01*HR1(1)*HR1(1)*HR3(0,0,-1)
     $  + 2.4022650695910071d-01*HR1(1)*HR2(-1,1)
     $  + 6.9314718055994530d-01*HR1(1)*HR3(-1,-1,1)
     $  - 6.9314718055994530d-01*HR1(1)*HR3(0,-1,-1)
     $  - 6.9314718055994530d-01*HR1(1)*HR3(0,-1,1)
     $  + HR1(1) *HR4(-1,-1,-1,1)
     $  + 2.0000000000000000d+00*HR1(1)*HR4(0,-1,-1,-1)
     $  + HR1(1) *HR4(0,-1,-1,1)
     $  + HR1(1) *HR4(0,-1,1,-1)
     $  + HR1(1) *HR4(0,0,-1,-1)
     $  + HR1(1) *HR4(0,0,-1,1)
     $  + 2.4022650695910071d-01*HR3(-1,-1,1)
     $  - 2.4022650695910071d-01*HR3(-1,1,1)
     $  + 1.3862943611198906d+00*HR4(-1,-1,-1,1)
     $  - 6.9314718055994530d-01*HR4(-1,-1,1,1)
     $  + 6.9314718055994530d-01*HR4(0,-1,-1,-1)
     $  + 6.9314718055994530d-01*HR4(0,-1,-1,1)
     $  + 6.9314718055994530d-01*HR4(0,-1,1,-1)
     $  + 6.9314718055994530d-01*HR4(0,-1,1,1)
     $  + 3.0000000000000000d+00*HR5(-1,-1,-1,-1,1)
     $  - HR5( -1,-1,-1,1,1)
     $  - 3.0000000000000000d+00*HR5(0,-1,-1,-1,-1)
     $  - 2.0000000000000000d+00*HR5(0,-1,-1,-1,1)
     $  - 2.0000000000000000d+00*HR5(0,-1,-1,1,-1)
     $  - HR5(0, -1,-1,1,1)
     $  - 2.0000000000000000d+00*HR5(0,-1,1,-1,-1)
     $  - HR5(0, -1,1,-1,1)
     $  - HR5(0, -1,1,1,-1)
     $  - HR5(0,0, -1,-1,-1)
     $  - HR5(0,0, -1,-1,1)
     $  - HR5(0,0, -1,1,-1)
     $  - HR5(0,0, -1,1,1)
      HY5(0,0,0,-1,1) = 
     $  + 5.9142607400864533d-02
     $  - 1.1787599965050932d-01*HR1(-1)
     $  + 1.2153517583503078d-01*HR1(-1)*HR1(-1)
     $  - 9.7040087744168750d-02*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 2.8881132523331054d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 8.3333333333333333d-03*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 8.3333333333333333d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  - 2.5000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)*HR1(1)
     $  - 2.9112026323250625d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - 1.7328679513998632d-01*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  + 8.3333333333333333d-02*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR2(0,-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(0,-1,-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(0,-1,1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(0)*HR1(1)*HR1(1)*HR1(1)
     $  + 2.4307035167006157d-01*HR1(-1)*HR1(1)
     $  - 2.9112026323250625d-01*HR1(-1)*HR1(1)*HR1(1)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(1)*HR1(1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(1)*HR1(1)*HR2(-1,1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(1)*HR1(1)*HR2(0,-1)
     $  - HR1( -1)*HR1(1)*HR3(-1,-1,1)
     $  + HR1( -1)*HR1(1)*HR3(-1,1,1)
     $  - HR1( -1)*HR1(1)*HR3(0,-1,-1)
     $  - HR1( -1)*HR1(1)*HR3(0,-1,1)
     $  - HR1( -1)*HR4(-1,-1,-1,1)
     $  + HR1( -1)*HR4(-1,-1,1,1)
     $  - HR1( -1)*HR4(-1,1,1,1)
     $  + HR1( -1)*HR4(0,-1,-1,-1)
     $  + HR1( -1)*HR4(0,-1,-1,1)
     $  + HR1( -1)*HR4(0,-1,1,-1)
     $  + HR1( -1)*HR4(0,-1,1,1)
     $  + 5.0000000000000000d-01*HR1(0)*HR1(1)*HR1(1)*HR2(-1,1)
     $  + HR1(0) *HR1(1)*HR3(-1,-1,1)
     $  - HR1(0) *HR1(1)*HR3(-1,1,1)
     $  + HR1(0) *HR4(-1,-1,-1,1)
     $  - HR1(0) *HR4(-1,-1,1,1)
     $  + HR1(0) *HR4(-1,1,1,1)
     $  - 1.1787599965050932d-01*HR1(1)
     $  + 1.2153517583503078d-01*HR1(1)*HR1(1)
     $  - 9.7040087744168750d-02*HR1(1)*HR1(1)*HR1(1)
     $  + 1.6666666666666666d-01*HR1(1)*HR1(1)*HR1(1)*HR2(0,-1)
     $  + 3.4657359027997265d-01*HR1(1)*HR1(1)*HR2(-1,1)
     $  + 5.0000000000000000d-01*HR1(1)*HR1(1)*HR3(-1,-1,1)
     $  - 5.0000000000000000d-01*HR1(1)*HR1(1)*HR3(0,-1,-1)
     $  - 5.0000000000000000d-01*HR1(1)*HR1(1)*HR3(0,-1,1)
     $  + 6.9314718055994530d-01*HR1(1)*HR3(-1,-1,1)
     $  - 6.9314718055994530d-01*HR1(1)*HR3(-1,1,1)
     $  + 2.0000000000000000d+00*HR1(1)*HR4(-1,-1,-1,1)
     $  - HR1(1) *HR4(-1,-1,1,1)
     $  + HR1(1) *HR4(0,-1,-1,-1)
     $  + HR1(1) *HR4(0,-1,-1,1)
     $  + HR1(1) *HR4(0,-1,1,-1)
     $  + HR1(1) *HR4(0,-1,1,1)
     $  + 6.9314718055994530d-01*HR4(-1,-1,-1,1)
     $  - 6.9314718055994530d-01*HR4(-1,-1,1,1)
     $  + 6.9314718055994530d-01*HR4(-1,1,1,1)
     $  + 3.0000000000000000d+00*HR5(-1,-1,-1,-1,1)
     $  - 2.0000000000000000d+00*HR5(-1,-1,-1,1,1)
     $  + HR5( -1,-1,1,1,1)
     $  - HR5(0, -1,-1,-1,-1)
     $  - HR5(0, -1,-1,-1,1)
     $  - HR5(0, -1,-1,1,-1)
     $  - HR5(0, -1,-1,1,1)
     $  - HR5(0, -1,1,-1,-1)
     $  - HR5(0, -1,1,-1,1)
     $  - HR5(0, -1,1,1,-1)
     $  - HR5(0, -1,1,1,1)
      HY5(0,0,0,1,-1) = 
     $  + 7.4276054639867797d-02
     $  - 1.7284527823898438d-01*HR1(-1)
     $  + 2.5410760640234242d-01*HR1(-1)*HR1(-1)
     $  + 1.7711559006386898d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 2.8881132523331054d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 8.3333333333333333d-03*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 8.3333333333333333d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  + 5.3134677019160696d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - 1.7328679513998632d-01*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  + 8.3333333333333333d-02*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR2(0,-1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR2(0,1)
     $  + HR1( -1)*HR1(-1)*HR3(0,-1,-1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(0,-1,1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(0,1,-1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(0)*HR1(1)*HR1(1)
     $  + 5.0821521280468485d-01*HR1(-1)*HR1(1)
     $  + 5.3134677019160696d-01*HR1(-1)*HR1(1)*HR1(1)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(1)*HR1(1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(1)*HR1(1)*HR2(-1,1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(1)*HR1(1)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(1)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(1)*HR2(0,1)
     $  - HR1( -1)*HR1(1)*HR3(-1,-1,1)
     $  + HR1( -1)*HR1(1)*HR3(-1,1,1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR1(1)*HR3(0,-1,-1)
     $  + HR1( -1)*HR1(1)*HR3(0,-1,1)
     $  + HR1( -1)*HR1(1)*HR3(0,1,-1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR3(0,-1,-1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR3(0,-1,1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR3(0,1,-1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR3(0,1,1)
     $  - HR1( -1)*HR4(-1,-1,-1,1)
     $  + HR1( -1)*HR4(-1,-1,1,1)
     $  - HR1( -1)*HR4(-1,1,1,1)
     $  - 3.0000000000000000d+00*HR1(-1)*HR4(0,-1,-1,-1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR4(0,-1,-1,1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR4(0,-1,1,-1)
     $  - HR1( -1)*HR4(0,-1,1,1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR4(0,1,-1,-1)
     $  - HR1( -1)*HR4(0,1,-1,1)
     $  - HR1( -1)*HR4(0,1,1,-1)
     $  + 1.1552453009332421d-01*HR1(0)*HR1(1)*HR1(1)*HR1(1)
     $  - 1.7284527823898438d-01*HR1(1)
     $  + 2.5410760640234242d-01*HR1(1)*HR1(1)
     $  + 1.7711559006386898d-01*HR1(1)*HR1(1)*HR1(1)
     $  - 1.6666666666666666d-01*HR1(1)*HR1(1)*HR1(1)*HR2(0,-1)
     $  + 3.4657359027997265d-01*HR1(1)*HR1(1)*HR2(-1,1)
     $  - 3.4657359027997265d-01*HR1(1)*HR1(1)*HR2(0,-1)
     $  - 3.4657359027997265d-01*HR1(1)*HR1(1)*HR2(0,1)
     $  + 5.0000000000000000d-01*HR1(1)*HR1(1)*HR3(-1,-1,1)
     $  + HR1(1) *HR1(1)*HR3(0,-1,-1)
     $  + 5.0000000000000000d-01*HR1(1)*HR1(1)*HR3(0,-1,1)
     $  + 5.0000000000000000d-01*HR1(1)*HR1(1)*HR3(0,1,-1)
     $  + 6.9314718055994530d-01*HR1(1)*HR3(-1,-1,1)
     $  - 6.9314718055994530d-01*HR1(1)*HR3(-1,1,1)
     $  + 6.9314718055994530d-01*HR1(1)*HR3(0,-1,-1)
     $  + 6.9314718055994530d-01*HR1(1)*HR3(0,-1,1)
     $  + 6.9314718055994530d-01*HR1(1)*HR3(0,1,-1)
     $  + 6.9314718055994530d-01*HR1(1)*HR3(0,1,1)
     $  + 2.0000000000000000d+00*HR1(1)*HR4(-1,-1,-1,1)
     $  - HR1(1) *HR4(-1,-1,1,1)
     $  - 3.0000000000000000d+00*HR1(1)*HR4(0,-1,-1,-1)
     $  - 2.0000000000000000d+00*HR1(1)*HR4(0,-1,-1,1)
     $  - 2.0000000000000000d+00*HR1(1)*HR4(0,-1,1,-1)
     $  - HR1(1) *HR4(0,-1,1,1)
     $  - 2.0000000000000000d+00*HR1(1)*HR4(0,1,-1,-1)
     $  - HR1(1) *HR4(0,1,-1,1)
     $  - HR1(1) *HR4(0,1,1,-1)
     $  + 6.9314718055994530d-01*HR4(-1,-1,-1,1)
     $  - 6.9314718055994530d-01*HR4(-1,-1,1,1)
     $  + 6.9314718055994530d-01*HR4(-1,1,1,1)
     $  - 6.9314718055994530d-01*HR4(0,-1,-1,-1)
     $  - 6.9314718055994530d-01*HR4(0,-1,-1,1)
     $  - 6.9314718055994530d-01*HR4(0,-1,1,-1)
     $  - 6.9314718055994530d-01*HR4(0,-1,1,1)
     $  - 6.9314718055994530d-01*HR4(0,1,-1,-1)
     $  - 6.9314718055994530d-01*HR4(0,1,-1,1)
     $  - 6.9314718055994530d-01*HR4(0,1,1,-1)
     $  - 6.9314718055994530d-01*HR4(0,1,1,1)
     $  + 3.0000000000000000d+00*HR5(-1,-1,-1,-1,1)
     $  - 2.0000000000000000d+00*HR5(-1,-1,-1,1,1)
     $  + HR5( -1,-1,1,1,1)
     $  + 4.0000000000000000d+00*HR5(0,-1,-1,-1,-1)
     $  + 3.0000000000000000d+00*HR5(0,-1,-1,-1,1)
     $  + 3.0000000000000000d+00*HR5(0,-1,-1,1,-1)
     $  + 2.0000000000000000d+00*HR5(0,-1,-1,1,1)
     $  + 3.0000000000000000d+00*HR5(0,-1,1,-1,-1)
     $  + 2.0000000000000000d+00*HR5(0,-1,1,-1,1)
     $  + 2.0000000000000000d+00*HR5(0,-1,1,1,-1)
     $  + HR5(0, -1,1,1,1)
     $  + 3.0000000000000000d+00*HR5(0,1,-1,-1,-1)
     $  + 2.0000000000000000d+00*HR5(0,1,-1,-1,1)
     $  + 2.0000000000000000d+00*HR5(0,1,-1,1,-1)
     $  + HR5(0,1, -1,1,1)
     $  + 2.0000000000000000d+00*HR5(0,1,1,-1,-1)
     $  + HR5(0,1,1, -1,1)
     $  + HR5(0,1,1,1, -1)
      HY5(0,0,1,-1,-1) = 
     $  + 2.5535023438634211d-02
     $  - 1.1412342741606084d-01*HR1(-1)
     $  - 2.3766885054564933d-01*HR1(-1)*HR1(-1)
     $  + 4.0037751159850118d-02*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 2.8881132523331054d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 8.3333333333333333d-03*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 8.3333333333333333d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  - 1.2011325347955035d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  + 1.2011325347955035d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - 1.7328679513998632d-01*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR2(-1,1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(-1,-1,1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(-1,1,1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(0,-1,-1)
     $  - 2.4022650695910071d-01*HR1(-1)*HR1(0)*HR1(1)
     $  - 4.7533770109129867d-01*HR1(-1)*HR1(1)
     $  + 1.2011325347955035d-01*HR1(-1)*HR1(1)*HR1(1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(1)*HR2(-1,1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(1)*HR2(0,-1)
     $  + HR1( -1)*HR1(1)*HR3(-1,-1,1)
     $  - HR1( -1)*HR1(1)*HR3(0,-1,-1)
     $  + 2.4022650695910071d-01*HR1(-1)*HR2(0,-1)
     $  + 2.4022650695910071d-01*HR1(-1)*HR2(0,1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR3(-1,-1,1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR3(-1,1,1)
     $  - 1.3862943611198906d+00*HR1(-1)*HR3(0,-1,-1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR3(0,-1,1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR3(0,1,-1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR4(-1,-1,-1,1)
     $  - HR1( -1)*HR4(-1,-1,1,1)
     $  + 3.0000000000000000d+00*HR1(-1)*HR4(0,-1,-1,-1)
     $  + HR1( -1)*HR4(0,-1,-1,1)
     $  + HR1( -1)*HR4(0,-1,1,-1)
     $  + HR1( -1)*HR4(0,1,-1,-1)
     $  - 1.2011325347955035d-01*HR1(0)*HR1(1)*HR1(1)
     $  - 1.1412342741606084d-01*HR1(1)
     $  - 2.3766885054564933d-01*HR1(1)*HR1(1)
     $  + 3.4657359027997265d-01*HR1(1)*HR1(1)*HR2(0,-1)
     $  - 5.0000000000000000d-01*HR1(1)*HR1(1)*HR3(0,-1,-1)
     $  - 2.4022650695910071d-01*HR1(1)*HR2(-1,1)
     $  + 2.4022650695910071d-01*HR1(1)*HR2(0,-1)
     $  + 2.4022650695910071d-01*HR1(1)*HR2(0,1)
     $  - 6.9314718055994530d-01*HR1(1)*HR3(-1,-1,1)
     $  - 1.3862943611198906d+00*HR1(1)*HR3(0,-1,-1)
     $  - 6.9314718055994530d-01*HR1(1)*HR3(0,-1,1)
     $  - 6.9314718055994530d-01*HR1(1)*HR3(0,1,-1)
     $  - HR1(1) *HR4(-1,-1,-1,1)
     $  + 3.0000000000000000d+00*HR1(1)*HR4(0,-1,-1,-1)
     $  + HR1(1) *HR4(0,-1,-1,1)
     $  + HR1(1) *HR4(0,-1,1,-1)
     $  + HR1(1) *HR4(0,1,-1,-1)
     $  - 2.4022650695910071d-01*HR3(-1,-1,1)
     $  + 2.4022650695910071d-01*HR3(-1,1,1)
     $  - 2.4022650695910071d-01*HR3(0,-1,-1)
     $  - 2.4022650695910071d-01*HR3(0,-1,1)
     $  - 2.4022650695910071d-01*HR3(0,1,-1)
     $  - 2.4022650695910071d-01*HR3(0,1,1)
     $  - 1.3862943611198906d+00*HR4(-1,-1,-1,1)
     $  + 6.9314718055994530d-01*HR4(-1,-1,1,1)
     $  + 2.0794415416798359d+00*HR4(0,-1,-1,-1)
     $  + 1.3862943611198906d+00*HR4(0,-1,-1,1)
     $  + 1.3862943611198906d+00*HR4(0,-1,1,-1)
     $  + 6.9314718055994530d-01*HR4(0,-1,1,1)
     $  + 1.3862943611198906d+00*HR4(0,1,-1,-1)
     $  + 6.9314718055994530d-01*HR4(0,1,-1,1)
     $  + 6.9314718055994530d-01*HR4(0,1,1,-1)
     $  - 3.0000000000000000d+00*HR5(-1,-1,-1,-1,1)
     $  + HR5( -1,-1,-1,1,1)
     $  - 6.0000000000000000d+00*HR5(0,-1,-1,-1,-1)
     $  - 3.0000000000000000d+00*HR5(0,-1,-1,-1,1)
     $  - 3.0000000000000000d+00*HR5(0,-1,-1,1,-1)
     $  - HR5(0, -1,-1,1,1)
     $  - 3.0000000000000000d+00*HR5(0,-1,1,-1,-1)
     $  - HR5(0, -1,1,-1,1)
     $  - HR5(0, -1,1,1,-1)
     $  - 3.0000000000000000d+00*HR5(0,1,-1,-1,-1)
     $  - HR5(0,1, -1,-1,1)
     $  - HR5(0,1, -1,1,-1)
     $  - HR5(0,1,1, -1,-1)
      HY5(0,0,1,-1,1) = 
     $  + 3.6321732111088421d-02
     $  - 1.9355535381306524d-01*HR1(-1)
     $  - 7.3900238327152102d-01*HR1(-1)*HR1(-1)
     $  + 9.7040087744168750d-02*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 2.8881132523331054d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 8.3333333333333333d-03*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 8.3333333333333333d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  - 2.9112026323250625d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  + 2.5000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR2(0,-1)
     $  + 2.9112026323250625d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  + 1.7328679513998632d-01*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR2(-1,1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR2(0,-1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(-1,-1,1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(-1,1,1)
     $  + HR1( -1)*HR1(-1)*HR3(0,-1,-1)
     $  + HR1( -1)*HR1(-1)*HR3(0,0,-1)
     $  - 5.8224052646501250d-01*HR1(-1)*HR1(0)*HR1(1)
     $  - HR1( -1)*HR1(0)*HR1(1)*HR2(-1,1)
     $  - HR1( -1)*HR1(0)*HR1(1)*HR2(0,-1)
     $  - HR1( -1)*HR1(0)*HR3(-1,-1,1)
     $  + HR1( -1)*HR1(0)*HR3(-1,1,1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR1(0)*HR3(0,-1,-1)
     $  + HR1( -1)*HR1(0)*HR3(0,-1,1)
     $  + HR1( -1)*HR1(0)*HR3(0,1,-1)
     $  - 1.4780047665430420d+00*HR1(-1)*HR1(1)
     $  + 2.9112026323250625d-01*HR1(-1)*HR1(1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(1)*HR1(1)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(1)*HR2(-1,1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(1)*HR2(0,-1)
     $  - HR1( -1)*HR1(1)*HR3(-1,-1,1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR1(1)*HR3(0,-1,-1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR1(1)*HR3(0,0,-1)
     $  + 5.8224052646501250d-01*HR1(-1)*HR2(0,-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR2(0,-1)*HR2(0,-1)
     $  - HR1( -1)*HR2(0,-1)*HR2(0,1)
     $  + 5.8224052646501250d-01*HR1(-1)*HR2(0,1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR3(-1,-1,1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR3(-1,1,1)
     $  + 1.3862943611198906d+00*HR1(-1)*HR3(0,-1,-1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR3(0,-1,1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR3(0,1,-1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR4(-1,-1,-1,1)
     $  + HR1( -1)*HR4(-1,-1,1,1)
     $  - 4.0000000000000000d+00*HR1(-1)*HR4(0,-1,-1,-1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR4(0,-1,-1,1)
     $  + HR1( -1)*HR4(0,-1,0,1)
     $  - HR1( -1)*HR4(0,-1,1,-1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR4(0,0,-1,-1)
     $  - HR1( -1)*HR4(0,1,-1,-1)
     $  - 2.9112026323250625d-01*HR1(0)*HR1(1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(0)*HR1(1)*HR1(1)*HR2(0,-1)
     $  + HR1(0) *HR1(1)*HR3(-1,-1,1)
     $  + 2.0000000000000000d+00*HR1(0)*HR1(1)*HR3(0,-1,-1)
     $  + HR1(0) *HR1(1)*HR3(0,-1,1)
     $  + HR1(0) *HR1(1)*HR3(0,1,-1)
     $  + 2.0000000000000000d+00*HR1(0)*HR4(-1,-1,-1,1)
     $  - HR1(0) *HR4(-1,-1,1,1)
     $  - 3.0000000000000000d+00*HR1(0)*HR4(0,-1,-1,-1)
     $  - 2.0000000000000000d+00*HR1(0)*HR4(0,-1,-1,1)
     $  - 2.0000000000000000d+00*HR1(0)*HR4(0,-1,1,-1)
     $  - HR1(0) *HR4(0,-1,1,1)
     $  - 2.0000000000000000d+00*HR1(0)*HR4(0,1,-1,-1)
     $  - HR1(0) *HR4(0,1,-1,1)
     $  - HR1(0) *HR4(0,1,1,-1)
     $  - 1.9355535381306524d-01*HR1(1)
     $  - 7.3900238327152102d-01*HR1(1)*HR1(1)
     $  - 3.4657359027997265d-01*HR1(1)*HR1(1)*HR2(0,-1)
     $  + HR1(1) *HR1(1)*HR3(0,-1,-1)
     $  + HR1(1) *HR1(1)*HR3(0,0,-1)
     $  - 5.8224052646501250d-01*HR1(1)*HR2(-1,1)
     $  + HR1(1) *HR2(-1,1)*HR2(0,-1)
     $  + 5.8224052646501250d-01*HR1(1)*HR2(0,-1)
     $  - 5.0000000000000000d-01*HR1(1)*HR2(0,-1)*HR2(0,-1)
     $  - HR1(1) *HR2(0,-1)*HR2(0,1)
     $  + 5.8224052646501250d-01*HR1(1)*HR2(0,1)
     $  + 6.9314718055994530d-01*HR1(1)*HR3(-1,-1,1)
     $  + 1.3862943611198906d+00*HR1(1)*HR3(0,-1,-1)
     $  + 6.9314718055994530d-01*HR1(1)*HR3(0,-1,1)
     $  + 6.9314718055994530d-01*HR1(1)*HR3(0,1,-1)
     $  + HR1(1) *HR4(-1,-1,-1,1)
     $  - 4.0000000000000000d+00*HR1(1)*HR4(0,-1,-1,-1)
     $  - 2.0000000000000000d+00*HR1(1)*HR4(0,-1,-1,1)
     $  + HR1(1) *HR4(0,-1,0,1)
     $  - HR1(1) *HR4(0,-1,1,-1)
     $  - 2.0000000000000000d+00*HR1(1)*HR4(0,0,-1,-1)
     $  - HR1(1) *HR4(0,1,-1,-1)
     $  + HR2(0, -1)*HR3(-1,-1,1)
     $  - HR2(0, -1)*HR3(-1,1,1)
     $  + HR2(0, -1)*HR3(0,-1,-1)
     $  + HR2(0, -1)*HR3(0,-1,1)
     $  + HR2(0, -1)*HR3(0,1,-1)
     $  + HR2(0, -1)*HR3(0,1,1)
     $  - 5.8224052646501250d-01*HR3(-1,-1,1)
     $  + 5.8224052646501250d-01*HR3(-1,1,1)
     $  - 5.8224052646501250d-01*HR3(0,-1,-1)
     $  - 5.8224052646501250d-01*HR3(0,-1,1)
     $  - 5.8224052646501250d-01*HR3(0,1,-1)
     $  - 5.8224052646501250d-01*HR3(0,1,1)
     $  + 1.3862943611198906d+00*HR4(-1,-1,-1,1)
     $  - 6.9314718055994530d-01*HR4(-1,-1,1,1)
     $  - 2.0794415416798359d+00*HR4(0,-1,-1,-1)
     $  - 1.3862943611198906d+00*HR4(0,-1,-1,1)
     $  - 1.3862943611198906d+00*HR4(0,-1,1,-1)
     $  - 6.9314718055994530d-01*HR4(0,-1,1,1)
     $  - 1.3862943611198906d+00*HR4(0,1,-1,-1)
     $  - 6.9314718055994530d-01*HR4(0,1,-1,1)
     $  - 6.9314718055994530d-01*HR4(0,1,1,-1)
     $  + 3.0000000000000000d+00*HR5(-1,-1,-1,-1,1)
     $  - HR5( -1,-1,-1,1,1)
     $  + 7.0000000000000000d+00*HR5(0,-1,-1,-1,-1)
     $  + 4.0000000000000000d+00*HR5(0,-1,-1,-1,1)
     $  + 4.0000000000000000d+00*HR5(0,-1,-1,1,-1)
     $  + 2.0000000000000000d+00*HR5(0,-1,-1,1,1)
     $  - HR5(0, -1,0,-1,-1)
     $  - HR5(0, -1,0,-1,1)
     $  - HR5(0, -1,0,1,-1)
     $  - HR5(0, -1,0,1,1)
     $  + 3.0000000000000000d+00*HR5(0,-1,1,-1,-1)
     $  + HR5(0, -1,1,-1,1)
     $  + HR5(0, -1,1,1,-1)
     $  + 3.0000000000000000d+00*HR5(0,1,-1,-1,-1)
     $  + HR5(0,1, -1,-1,1)
     $  + HR5(0,1, -1,1,-1)
     $  + HR5(0,1,1, -1,-1)
      HY5(0,0,1,0,-1) = 
     $  + 1.8615775173851248d-01
     $  - 5.6852588003909690d-01*HR1(-1)
     $  - 6.6068813489808640d-01*HR1(-1)*HR1(-1)
     $  + 1.3707783890401886d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 2.8881132523331054d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 8.3333333333333333d-03*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 8.3333333333333333d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  - 4.1123351671205660d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  + 4.1123351671205660d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - 1.7328679513998632d-01*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR2(-1,1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR2(0,1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(-1,-1,1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(-1,1,1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(0,-1,-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(0,1,-1)
     $  - 8.2246703342411321d-01*HR1(-1)*HR1(0)*HR1(1)
     $  - 1.3213762697961728d+00*HR1(-1)*HR1(1)
     $  + 4.1123351671205660d-01*HR1(-1)*HR1(1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(1)*HR1(1)*HR2(-1,1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(1)*HR2(-1,1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(1)*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(1)*HR2(0,1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR1(1)*HR3(-1,-1,1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR1(1)*HR3(-1,1,1)
     $  - HR1( -1)*HR1(1)*HR3(0,-1,-1)
     $  - HR1( -1)*HR1(1)*HR3(0,1,-1)
     $  + 8.2246703342411321d-01*HR1(-1)*HR2(0,-1)
     $  + 8.2246703342411321d-01*HR1(-1)*HR2(0,1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR3(-1,-1,1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR3(-1,1,1)
     $  - 1.3862943611198906d+00*HR1(-1)*HR3(0,-1,-1)
     $  - 1.3862943611198906d+00*HR1(-1)*HR3(0,-1,1)
     $  - 1.3862943611198906d+00*HR1(-1)*HR3(0,1,-1)
     $  - 1.3862943611198906d+00*HR1(-1)*HR3(0,1,1)
     $  + 3.0000000000000000d+00*HR1(-1)*HR4(-1,-1,-1,1)
     $  - 3.0000000000000000d+00*HR1(-1)*HR4(-1,-1,1,1)
     $  + 3.0000000000000000d+00*HR1(-1)*HR4(-1,1,1,1)
     $  + 3.0000000000000000d+00*HR1(-1)*HR4(0,-1,-1,-1)
     $  + HR1( -1)*HR4(0,-1,-1,1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR4(0,-1,1,-1)
     $  + 3.0000000000000000d+00*HR1(-1)*HR4(0,1,-1,-1)
     $  + HR1( -1)*HR4(0,1,-1,1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR4(0,1,1,-1)
     $  - 4.1123351671205660d-01*HR1(0)*HR1(1)*HR1(1)
     $  - 5.6852588003909690d-01*HR1(1)
     $  - 6.6068813489808640d-01*HR1(1)*HR1(1)
     $  - 3.4657359027997265d-01*HR1(1)*HR1(1)*HR2(-1,1)
     $  + 3.4657359027997265d-01*HR1(1)*HR1(1)*HR2(0,-1)
     $  + 3.4657359027997265d-01*HR1(1)*HR1(1)*HR2(0,1)
     $  - HR1(1) *HR1(1)*HR3(-1,-1,1)
     $  - 5.0000000000000000d-01*HR1(1)*HR1(1)*HR3(0,-1,-1)
     $  - 5.0000000000000000d-01*HR1(1)*HR1(1)*HR3(0,1,-1)
     $  - 8.2246703342411321d-01*HR1(1)*HR2(-1,1)
     $  + 5.0000000000000000d-01*HR1(1)*HR2(-1,1)*HR2(-1,1)
     $  + 8.2246703342411321d-01*HR1(1)*HR2(0,-1)
     $  + 8.2246703342411321d-01*HR1(1)*HR2(0,1)
     $  - 1.3862943611198906d+00*HR1(1)*HR3(-1,-1,1)
     $  + 1.3862943611198906d+00*HR1(1)*HR3(-1,1,1)
     $  - 1.3862943611198906d+00*HR1(1)*HR3(0,-1,-1)
     $  - 1.3862943611198906d+00*HR1(1)*HR3(0,-1,1)
     $  - 1.3862943611198906d+00*HR1(1)*HR3(0,1,-1)
     $  - 1.3862943611198906d+00*HR1(1)*HR3(0,1,1)
     $  - 4.0000000000000000d+00*HR1(1)*HR4(-1,-1,-1,1)
     $  + 2.0000000000000000d+00*HR1(1)*HR4(-1,-1,1,1)
     $  + 3.0000000000000000d+00*HR1(1)*HR4(0,-1,-1,-1)
     $  + HR1(1) *HR4(0,-1,-1,1)
     $  + 2.0000000000000000d+00*HR1(1)*HR4(0,-1,1,-1)
     $  + 3.0000000000000000d+00*HR1(1)*HR4(0,1,-1,-1)
     $  + HR1(1) *HR4(0,1,-1,1)
     $  + 2.0000000000000000d+00*HR1(1)*HR4(0,1,1,-1)
     $  - HR2( -1,1)*HR3(-1,1,1)
     $  - 8.2246703342411321d-01*HR3(-1,-1,1)
     $  + 8.2246703342411321d-01*HR3(-1,1,1)
     $  - 8.2246703342411321d-01*HR3(0,-1,-1)
     $  - 8.2246703342411321d-01*HR3(0,-1,1)
     $  - 8.2246703342411321d-01*HR3(0,1,-1)
     $  - 8.2246703342411321d-01*HR3(0,1,1)
     $  - 2.0794415416798359d+00*HR4(-1,-1,-1,1)
     $  + 2.0794415416798359d+00*HR4(-1,-1,1,1)
     $  - 2.0794415416798359d+00*HR4(-1,1,1,1)
     $  + 2.0794415416798359d+00*HR4(0,-1,-1,-1)
     $  + 2.0794415416798359d+00*HR4(0,-1,-1,1)
     $  + 2.0794415416798359d+00*HR4(0,-1,1,-1)
     $  + 2.0794415416798359d+00*HR4(0,-1,1,1)
     $  + 2.0794415416798359d+00*HR4(0,1,-1,-1)
     $  + 2.0794415416798359d+00*HR4(0,1,-1,1)
     $  + 2.0794415416798359d+00*HR4(0,1,1,-1)
     $  + 2.0794415416798359d+00*HR4(0,1,1,1)
     $  - 7.0000000000000000d+00*HR5(-1,-1,-1,-1,1)
     $  + 7.0000000000000000d+00*HR5(-1,-1,-1,1,1)
     $  + HR5( -1,-1,1,-1,1)
     $  + HR5( -1,1,-1,1,1)
     $  - 6.0000000000000000d+00*HR5(0,-1,-1,-1,-1)
     $  - 3.0000000000000000d+00*HR5(0,-1,-1,-1,1)
     $  - 4.0000000000000000d+00*HR5(0,-1,-1,1,-1)
     $  - HR5(0, -1,-1,1,1)
     $  - 5.0000000000000000d+00*HR5(0,-1,1,-1,-1)
     $  - 2.0000000000000000d+00*HR5(0,-1,1,-1,1)
     $  - 3.0000000000000000d+00*HR5(0,-1,1,1,-1)
     $  - 6.0000000000000000d+00*HR5(0,1,-1,-1,-1)
     $  - 3.0000000000000000d+00*HR5(0,1,-1,-1,1)
     $  - 4.0000000000000000d+00*HR5(0,1,-1,1,-1)
     $  - HR5(0,1, -1,1,1)
     $  - 5.0000000000000000d+00*HR5(0,1,1,-1,-1)
     $  - 2.0000000000000000d+00*HR5(0,1,1,-1,1)
     $  - 3.0000000000000000d+00*HR5(0,1,1,1,-1)
      HY5(0,0,1,1,-1) = 
     $  + 5.7353571803049304d-02
     $  - 4.3369237704895519d-01*HR1(-1)
     $  + 5.5365194946473328d-01*HR1(-1)*HR1(-1)
     $  - 1.7711559006386898d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 2.8881132523331054d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 8.3333333333333333d-03*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 8.3333333333333333d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  + 5.3134677019160696d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  + 1.7328679513998632d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(0)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  - 5.3134677019160696d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  + 1.7328679513998632d-01*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR2(-1,1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR2(0,-1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(-1,-1,1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(-1,1,1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(0,-1,-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(0,0,-1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(0)*HR1(0)*HR1(1)
     $  + 1.0626935403832139d+00*HR1(-1)*HR1(0)*HR1(1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(0)*HR1(1)*HR1(1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(0)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(0)*HR2(0,1)
     $  + 1.1073038989294665d+00*HR1(-1)*HR1(1)
     $  - 5.3134677019160696d-01*HR1(-1)*HR1(1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(1)*HR1(1)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(1)*HR2(-1,1)
     $  - HR1( -1)*HR1(1)*HR3(-1,-1,1)
     $  - HR1( -1)*HR1(1)*HR3(0,-1,-1)
     $  - HR1( -1)*HR1(1)*HR3(0,0,-1)
     $  - 1.0626935403832139d+00*HR1(-1)*HR2(0,-1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR2(0,-1)*HR2(0,-1)
     $  + HR1( -1)*HR2(0,-1)*HR2(0,1)
     $  - 1.0626935403832139d+00*HR1(-1)*HR2(0,1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR3(-1,-1,1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR3(-1,1,1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR3(0,-1,-1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR3(0,0,-1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR3(0,0,1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR3(0,1,-1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR4(-1,-1,-1,1)
     $  + HR1( -1)*HR4(-1,-1,1,1)
     $  + HR1( -1)*HR4(0,-1,-1,1)
     $  - HR1( -1)*HR4(0,-1,0,1)
     $  - HR1( -1)*HR4(0,0,-1,1)
     $  - HR1( -1)*HR4(0,0,1,-1)
     $  - HR1( -1)*HR4(0,1,-1,-1)
     $  + 1.7328679513998632d-01*HR1(0)*HR1(0)*HR1(1)*HR1(1)
     $  + 5.3134677019160696d-01*HR1(0)*HR1(1)*HR1(1)
     $  + 6.9314718055994530d-01*HR1(0)*HR1(1)*HR2(-1,1)
     $  - 6.9314718055994530d-01*HR1(0)*HR1(1)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(0)*HR1(1)*HR2(0,1)
     $  + 6.9314718055994530d-01*HR1(0)*HR3(-1,-1,1)
     $  - 6.9314718055994530d-01*HR1(0)*HR3(-1,1,1)
     $  + 6.9314718055994530d-01*HR1(0)*HR3(0,-1,-1)
     $  + 6.9314718055994530d-01*HR1(0)*HR3(0,-1,1)
     $  + 6.9314718055994530d-01*HR1(0)*HR3(0,1,-1)
     $  + 6.9314718055994530d-01*HR1(0)*HR3(0,1,1)
     $  - 4.3369237704895519d-01*HR1(1)
     $  + 5.5365194946473328d-01*HR1(1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(1)*HR1(1)*HR3(0,-1,-1)
     $  - 5.0000000000000000d-01*HR1(1)*HR1(1)*HR3(0,0,-1)
     $  + 1.0626935403832139d+00*HR1(1)*HR2(-1,1)
     $  - HR1(1) *HR2(-1,1)*HR2(0,-1)
     $  - 1.0626935403832139d+00*HR1(1)*HR2(0,-1)
     $  + 5.0000000000000000d-01*HR1(1)*HR2(0,-1)*HR2(0,-1)
     $  + HR1(1) *HR2(0,-1)*HR2(0,1)
     $  - 1.0626935403832139d+00*HR1(1)*HR2(0,1)
     $  + 6.9314718055994530d-01*HR1(1)*HR3(-1,-1,1)
     $  + 6.9314718055994530d-01*HR1(1)*HR3(0,-1,-1)
     $  + 6.9314718055994530d-01*HR1(1)*HR3(0,0,-1)
     $  + 6.9314718055994530d-01*HR1(1)*HR3(0,0,1)
     $  + 6.9314718055994530d-01*HR1(1)*HR3(0,1,-1)
     $  + HR1(1) *HR4(-1,-1,-1,1)
     $  + HR1(1) *HR4(0,-1,-1,1)
     $  - HR1(1) *HR4(0,-1,0,1)
     $  - HR1(1) *HR4(0,0,-1,1)
     $  - HR1(1) *HR4(0,0,1,-1)
     $  - HR1(1) *HR4(0,1,-1,-1)
     $  - HR2(0, -1)*HR3(-1,-1,1)
     $  + HR2(0, -1)*HR3(-1,1,1)
     $  - HR2(0, -1)*HR3(0,-1,-1)
     $  - HR2(0, -1)*HR3(0,-1,1)
     $  - HR2(0, -1)*HR3(0,1,-1)
     $  - HR2(0, -1)*HR3(0,1,1)
     $  + 1.0626935403832139d+00*HR3(-1,-1,1)
     $  - 1.0626935403832139d+00*HR3(-1,1,1)
     $  + 1.0626935403832139d+00*HR3(0,-1,-1)
     $  + 1.0626935403832139d+00*HR3(0,-1,1)
     $  + 1.0626935403832139d+00*HR3(0,1,-1)
     $  + 1.0626935403832139d+00*HR3(0,1,1)
     $  + 1.3862943611198906d+00*HR4(-1,-1,-1,1)
     $  - 6.9314718055994530d-01*HR4(-1,-1,1,1)
     $  - 1.3862943611198906d+00*HR4(0,-1,-1,-1)
     $  - 6.9314718055994530d-01*HR4(0,-1,-1,1)
     $  - 6.9314718055994530d-01*HR4(0,-1,1,-1)
     $  - 6.9314718055994530d-01*HR4(0,0,-1,-1)
     $  - 6.9314718055994530d-01*HR4(0,0,-1,1)
     $  - 6.9314718055994530d-01*HR4(0,0,1,-1)
     $  - 6.9314718055994530d-01*HR4(0,0,1,1)
     $  - 1.3862943611198906d+00*HR4(0,1,-1,-1)
     $  - 6.9314718055994530d-01*HR4(0,1,-1,1)
     $  - 6.9314718055994530d-01*HR4(0,1,1,-1)
     $  + 3.0000000000000000d+00*HR5(-1,-1,-1,-1,1)
     $  - HR5( -1,-1,-1,1,1)
     $  + 2.0000000000000000d+00*HR5(0,-1,-1,-1,-1)
     $  - HR5(0, -1,-1,1,1)
     $  + HR5(0, -1,0,-1,-1)
     $  + HR5(0, -1,0,-1,1)
     $  + HR5(0, -1,0,1,-1)
     $  + HR5(0, -1,0,1,1)
     $  + HR5(0, -1,1,-1,-1)
     $  + 3.0000000000000000d+00*HR5(0,0,-1,-1,-1)
     $  + 2.0000000000000000d+00*HR5(0,0,-1,-1,1)
     $  + 2.0000000000000000d+00*HR5(0,0,-1,1,-1)
     $  + HR5(0,0, -1,1,1)
     $  + 2.0000000000000000d+00*HR5(0,0,1,-1,-1)
     $  + HR5(0,0,1, -1,1)
     $  + HR5(0,0,1,1, -1)
     $  + 3.0000000000000000d+00*HR5(0,1,-1,-1,-1)
     $  + HR5(0,1, -1,-1,1)
     $  + HR5(0,1, -1,1,-1)
     $  + HR5(0,1,1, -1,-1)
      HY5(0,1,-1,-1,-1) = 
     $  + 1.9555438852482933d-02
     $  + 1.2679858379652411d-01*HR1(-1)
     $  - 2.7752054332410789d-02*HR1(-1)*HR1(-1)
     $  + 4.0037751159850118d-02*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 2.8881132523331054d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 8.3333333333333333d-03*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  + 1.2011325347955035d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(-1,-1,1)
     $  + 5.5504108664821579d-02*HR1(-1)*HR1(0)
     $  - 5.5504108664821579d-02*HR1(-1)*HR1(1)
     $  - 2.4022650695910071d-01*HR1(-1)*HR2(-1,1)
     $  - 2.4022650695910071d-01*HR1(-1)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR3(-1,-1,1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR3(0,-1,-1)
     $  - HR1( -1)*HR4(-1,-1,-1,1)
     $  - HR1( -1)*HR4(0,-1,-1,-1)
     $  + 5.5504108664821579d-02*HR1(0)*HR1(1)
     $  + 1.2679858379652411d-01*HR1(1)
     $  - 2.4022650695910071d-01*HR1(1)*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR1(1)*HR3(0,-1,-1)
     $  - HR1(1) *HR4(0,-1,-1,-1)
     $  + 5.5504108664821579d-02*HR2(-1,1)
     $  - 5.5504108664821579d-02*HR2(0,-1)
     $  - 5.5504108664821579d-02*HR2(0,1)
     $  + 2.4022650695910071d-01*HR3(-1,-1,1)
     $  + 4.8045301391820142d-01*HR3(0,-1,-1)
     $  + 2.4022650695910071d-01*HR3(0,-1,1)
     $  + 2.4022650695910071d-01*HR3(0,1,-1)
     $  + 6.9314718055994530d-01*HR4(-1,-1,-1,1)
     $  - 2.0794415416798359d+00*HR4(0,-1,-1,-1)
     $  - 6.9314718055994530d-01*HR4(0,-1,-1,1)
     $  - 6.9314718055994530d-01*HR4(0,-1,1,-1)
     $  - 6.9314718055994530d-01*HR4(0,1,-1,-1)
     $  + HR5( -1,-1,-1,-1,1)
     $  + 4.0000000000000000d+00*HR5(0,-1,-1,-1,-1)
     $  + HR5(0, -1,-1,-1,1)
     $  + HR5(0, -1,-1,1,-1)
     $  + HR5(0, -1,1,-1,-1)
     $  + HR5(0,1, -1,-1,-1)
      HY5(0,1,-1,-1,1) = 
     $  + 2.8668668263701248d-02
     $  + 2.3517979306082505d-01*HR1(-1)
     $  - 4.7376502115063852d-02*HR1(-1)*HR1(-1)
     $  + 9.7040087744168750d-02*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 2.8881132523331054d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 8.3333333333333333d-03*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR2(-1,1)
     $  + 2.9112026323250625d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR2(0,-1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(-1,-1,1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(0,-1,-1)
     $  + 9.4753004230127705d-02*HR1(-1)*HR1(0)
     $  + HR1( -1)*HR1(0)*HR3(-1,-1,1)
     $  - HR1( -1)*HR1(0)*HR3(0,-1,-1)
     $  - 9.4753004230127705d-02*HR1(-1)*HR1(1)
     $  + HR1( -1)*HR1(1)*HR3(0,-1,-1)
     $  - 5.8224052646501250d-01*HR1(-1)*HR2(-1,1)
     $  + HR1( -1)*HR2(-1,1)*HR2(0,-1)
     $  - 5.8224052646501250d-01*HR1(-1)*HR2(0,-1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR2(0,-1)*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR3(-1,-1,1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR3(0,-1,-1)
     $  + HR1( -1)*HR4(-1,-1,-1,1)
     $  + 9.4753004230127705d-02*HR1(0)*HR1(1)
     $  - HR1(0) *HR1(1)*HR3(0,-1,-1)
     $  - HR1(0) *HR4(-1,-1,-1,1)
     $  + 3.0000000000000000d+00*HR1(0)*HR4(0,-1,-1,-1)
     $  + HR1(0) *HR4(0,-1,-1,1)
     $  + HR1(0) *HR4(0,-1,1,-1)
     $  + HR1(0) *HR4(0,1,-1,-1)
     $  + 2.3517979306082505d-01*HR1(1)
     $  - 5.8224052646501250d-01*HR1(1)*HR2(0,-1)
     $  + 5.0000000000000000d-01*HR1(1)*HR2(0,-1)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(1)*HR3(0,-1,-1)
     $  + 9.4753004230127705d-02*HR2(-1,1)
     $  - HR2( -1,1)*HR3(0,-1,-1)
     $  - 9.4753004230127705d-02*HR2(0,-1)
     $  - HR2(0, -1)*HR3(-1,-1,1)
     $  - 2.0000000000000000d+00*HR2(0,-1)*HR3(0,-1,-1)
     $  - HR2(0, -1)*HR3(0,-1,1)
     $  - HR2(0, -1)*HR3(0,1,-1)
     $  - 9.4753004230127705d-02*HR2(0,1)
     $  + HR2(0,1) *HR3(0,-1,-1)
     $  + 5.8224052646501250d-01*HR3(-1,-1,1)
     $  + 1.1644810529300250d+00*HR3(0,-1,-1)
     $  + 5.8224052646501250d-01*HR3(0,-1,1)
     $  + 5.8224052646501250d-01*HR3(0,1,-1)
     $  - 6.9314718055994530d-01*HR4(-1,-1,-1,1)
     $  + 2.0794415416798359d+00*HR4(0,-1,-1,-1)
     $  + 6.9314718055994530d-01*HR4(0,-1,-1,1)
     $  + 6.9314718055994530d-01*HR4(0,-1,1,-1)
     $  + 6.9314718055994530d-01*HR4(0,1,-1,-1)
     $  - HR5( -1,-1,-1,-1,1)
     $  - 3.0000000000000000d+00*HR5(0,-1,-1,-1,-1)
     $  - HR5(0, -1,-1,0,1)
     $  - HR5(0, -1,-1,1,-1)
     $  + 3.0000000000000000d+00*HR5(0,-1,0,-1,-1)
     $  - HR5(0, -1,1,-1,-1)
     $  + 6.0000000000000000d+00*HR5(0,0,-1,-1,-1)
     $  - HR5(0,1, -1,-1,-1)
      HY5(0,1,-1,1,-1) = 
     $  + 4.7069633474401836d-02
     $  + 6.4210078767232862d-01*HR1(-1)
     $  - 1.0703618543335311d-01*HR1(-1)*HR1(-1)
     $  - 1.7711559006386898d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 2.8881132523331054d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 8.3333333333333333d-03*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  - 5.3134677019160696d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR2(0,-1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(-1,-1,1)
     $  - HR1( -1)*HR1(-1)*HR3(0,-1,-1)
     $  + 2.1407237086670622d-01*HR1(-1)*HR1(0)
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(0)*HR2(-1,1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(0)*HR2(0,-1)
     $  - 2.1407237086670622d-01*HR1(-1)*HR1(1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(1)*HR2(0,-1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR1(1)*HR3(0,-1,-1)
     $  + 1.0626935403832139d+00*HR1(-1)*HR2(-1,1)
     $  - HR1( -1)*HR2(-1,1)*HR2(0,-1)
     $  + 1.0626935403832139d+00*HR1(-1)*HR2(0,-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR2(0,-1)*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR3(-1,-1,1)
     $  - 1.3862943611198906d+00*HR1(-1)*HR3(0,-1,-1)
     $  - 1.3862943611198906d+00*HR1(-1)*HR3(0,0,-1)
     $  + HR1( -1)*HR4(-1,-1,-1,1)
     $  + 4.0000000000000000d+00*HR1(-1)*HR4(0,-1,-1,-1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR4(0,0,-1,-1)
     $  + 2.1407237086670622d-01*HR1(0)*HR1(1)
     $  + 6.9314718055994530d-01*HR1(0)*HR1(1)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(0)*HR3(-1,-1,1)
     $  - 1.3862943611198906d+00*HR1(0)*HR3(0,-1,-1)
     $  - 6.9314718055994530d-01*HR1(0)*HR3(0,-1,1)
     $  - 6.9314718055994530d-01*HR1(0)*HR3(0,1,-1)
     $  + 6.4210078767232862d-01*HR1(1)
     $  + 1.0626935403832139d+00*HR1(1)*HR2(0,-1)
     $  - 5.0000000000000000d-01*HR1(1)*HR2(0,-1)*HR2(0,-1)
     $  - 1.3862943611198906d+00*HR1(1)*HR3(0,-1,-1)
     $  - 1.3862943611198906d+00*HR1(1)*HR3(0,0,-1)
     $  + 4.0000000000000000d+00*HR1(1)*HR4(0,-1,-1,-1)
     $  + 2.0000000000000000d+00*HR1(1)*HR4(0,0,-1,-1)
     $  + 2.1407237086670622d-01*HR2(-1,1)
     $  - 6.9314718055994530d-01*HR2(-1,1)*HR2(0,-1)
     $  + 2.0000000000000000d+00*HR2(-1,1)*HR3(0,-1,-1)
     $  - 2.1407237086670622d-01*HR2(0,-1)
     $  + 3.4657359027997265d-01*HR2(0,-1)*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR2(0,-1)*HR2(0,1)
     $  + HR2(0, -1)*HR3(-1,-1,1)
     $  + 2.0000000000000000d+00*HR2(0,-1)*HR3(0,-1,-1)
     $  + HR2(0, -1)*HR3(0,-1,1)
     $  + HR2(0, -1)*HR3(0,1,-1)
     $  - 2.1407237086670622d-01*HR2(0,1)
     $  - 2.0000000000000000d+00*HR2(0,1)*HR3(0,-1,-1)
     $  - 1.0626935403832139d+00*HR3(-1,-1,1)
     $  - 2.1253870807664278d+00*HR3(0,-1,-1)
     $  - 1.0626935403832139d+00*HR3(0,-1,1)
     $  - 1.0626935403832139d+00*HR3(0,1,-1)
     $  - 6.9314718055994530d-01*HR4(-1,-1,-1,1)
     $  + 2.7725887222397812d+00*HR4(0,-1,-1,-1)
     $  + 1.3862943611198906d+00*HR4(0,-1,-1,1)
     $  - 6.9314718055994530d-01*HR4(0,-1,0,1)
     $  + 6.9314718055994530d-01*HR4(0,-1,1,-1)
     $  + 1.3862943611198906d+00*HR4(0,0,-1,-1)
     $  + 6.9314718055994530d-01*HR4(0,1,-1,-1)
     $  - HR5( -1,-1,-1,-1,1)
     $  - 8.0000000000000000d+00*HR5(0,-1,-1,-1,-1)
     $  - 4.0000000000000000d+00*HR5(0,-1,-1,-1,1)
     $  + 2.0000000000000000d+00*HR5(0,-1,-1,0,1)
     $  - 2.0000000000000000d+00*HR5(0,-1,-1,1,-1)
     $  - 4.0000000000000000d+00*HR5(0,-1,0,-1,-1)
     $  + HR5(0, -1,0,-1,1)
     $  + HR5(0, -1,0,1,-1)
     $  - HR5(0, -1,1,-1,-1)
     $  - 1.2000000000000000d+01*HR5(0,0,-1,-1,-1)
     $  - HR5(0,1, -1,-1,-1)
      HY5(0,1,-1,1,1) = 
     $  + 8.2208743029471844d-02
     $  + 1.9248049955307152d+00*HR1(-1)
     $  - 2.6860659680402010d-01*HR1(-1)*HR1(-1)
     $  + 4.0037751159850118d-02*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 2.8881132523331054d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 8.3333333333333333d-03*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  + 8.3333333333333333d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)*HR1(0)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  + 2.5000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(0)*HR1(1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR2(-1,1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR2(0,-1)
     $  + 1.2011325347955035d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(-1,-1,1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(0,-1,-1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(0,0,-1)
     $  + 5.3721319360804020d-01*HR1(-1)*HR1(0)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(0)*HR1(0)*HR2(-1,1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(0)*HR1(0)*HR2(0,-1)
     $  - HR1( -1)*HR1(0)*HR1(1)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(0)*HR2(-1,1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(0)*HR2(0,-1)
     $  - HR1( -1)*HR1(0)*HR3(-1,-1,1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR1(0)*HR3(0,-1,-1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR1(0)*HR3(0,0,-1)
     $  - 5.3721319360804020d-01*HR1(-1)*HR1(1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(1)*HR2(0,-1)
     $  + HR1( -1)*HR1(1)*HR3(0,-1,-1)
     $  + HR1( -1)*HR1(1)*HR3(0,0,-1)
     $  - 2.4022650695910071d-01*HR1(-1)*HR2(-1,1)
     $  - 2.4022650695910071d-01*HR1(-1)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR3(-1,-1,1)
     $  + 1.3862943611198906d+00*HR1(-1)*HR3(0,-1,-1)
     $  + 1.3862943611198906d+00*HR1(-1)*HR3(0,0,-1)
     $  - HR1( -1)*HR4(-1,-1,-1,1)
     $  - 3.0000000000000000d+00*HR1(-1)*HR4(0,-1,-1,-1)
     $  - 3.0000000000000000d+00*HR1(-1)*HR4(0,0,-1,-1)
     $  - 3.0000000000000000d+00*HR1(-1)*HR4(0,0,0,-1)
     $  - 5.0000000000000000d-01*HR1(0)*HR1(0)*HR1(1)*HR2(0,-1)
     $  + 5.0000000000000000d-01*HR1(0)*HR1(0)*HR3(-1,-1,1)
     $  + HR1(0) *HR1(0)*HR3(0,-1,-1)
     $  + 5.0000000000000000d-01*HR1(0)*HR1(0)*HR3(0,-1,1)
     $  + 5.0000000000000000d-01*HR1(0)*HR1(0)*HR3(0,1,-1)
     $  + 5.3721319360804020d-01*HR1(0)*HR1(1)
     $  - 6.9314718055994530d-01*HR1(0)*HR1(1)*HR2(0,-1)
     $  + 2.0000000000000000d+00*HR1(0)*HR1(1)*HR3(0,-1,-1)
     $  + 2.0000000000000000d+00*HR1(0)*HR1(1)*HR3(0,0,-1)
     $  + HR1(0) *HR2(-1,1)*HR2(0,-1)
     $  - 5.0000000000000000d-01*HR1(0)*HR2(0,-1)*HR2(0,-1)
     $  - HR1(0) *HR2(0,-1)*HR2(0,1)
     $  + 6.9314718055994530d-01*HR1(0)*HR3(-1,-1,1)
     $  + 1.3862943611198906d+00*HR1(0)*HR3(0,-1,-1)
     $  + 6.9314718055994530d-01*HR1(0)*HR3(0,-1,1)
     $  + 6.9314718055994530d-01*HR1(0)*HR3(0,1,-1)
     $  + HR1(0) *HR4(-1,-1,-1,1)
     $  - 4.0000000000000000d+00*HR1(0)*HR4(0,-1,-1,-1)
     $  - 2.0000000000000000d+00*HR1(0)*HR4(0,-1,-1,1)
     $  + HR1(0) *HR4(0,-1,0,1)
     $  - HR1(0) *HR4(0,-1,1,-1)
     $  - 2.0000000000000000d+00*HR1(0)*HR4(0,0,-1,-1)
     $  - HR1(0) *HR4(0,1,-1,-1)
     $  + 1.9248049955307152d+00*HR1(1)
     $  - 2.4022650695910071d-01*HR1(1)*HR2(0,-1)
     $  + 1.3862943611198906d+00*HR1(1)*HR3(0,-1,-1)
     $  + 1.3862943611198906d+00*HR1(1)*HR3(0,0,-1)
     $  - 3.0000000000000000d+00*HR1(1)*HR4(0,-1,-1,-1)
     $  - 3.0000000000000000d+00*HR1(1)*HR4(0,0,-1,-1)
     $  - 3.0000000000000000d+00*HR1(1)*HR4(0,0,0,-1)
     $  + 5.3721319360804020d-01*HR2(-1,1)
     $  + 6.9314718055994530d-01*HR2(-1,1)*HR2(0,-1)
     $  - HR2( -1,1)*HR3(0,-1,-1)
     $  - HR2( -1,1)*HR3(0,0,-1)
     $  - 5.3721319360804020d-01*HR2(0,-1)
     $  - 3.4657359027997265d-01*HR2(0,-1)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR2(0,-1)*HR2(0,1)
     $  + HR2(0, -1)*HR3(0,0,-1)
     $  - 5.3721319360804020d-01*HR2(0,1)
     $  + HR2(0,1) *HR3(0,-1,-1)
     $  + HR2(0,1) *HR3(0,0,-1)
     $  + 2.4022650695910071d-01*HR3(-1,-1,1)
     $  + 4.8045301391820142d-01*HR3(0,-1,-1)
     $  + 2.4022650695910071d-01*HR3(0,-1,1)
     $  + 2.4022650695910071d-01*HR3(0,1,-1)
     $  + 6.9314718055994530d-01*HR4(-1,-1,-1,1)
     $  - 2.7725887222397812d+00*HR4(0,-1,-1,-1)
     $  - 1.3862943611198906d+00*HR4(0,-1,-1,1)
     $  + 6.9314718055994530d-01*HR4(0,-1,0,1)
     $  - 6.9314718055994530d-01*HR4(0,-1,1,-1)
     $  - 1.3862943611198906d+00*HR4(0,0,-1,-1)
     $  - 6.9314718055994530d-01*HR4(0,1,-1,-1)
     $  + HR5( -1,-1,-1,-1,1)
     $  + 7.0000000000000000d+00*HR5(0,-1,-1,-1,-1)
     $  + 3.0000000000000000d+00*HR5(0,-1,-1,-1,1)
     $  - HR5(0, -1,-1,0,1)
     $  + 2.0000000000000000d+00*HR5(0,-1,-1,1,-1)
     $  + HR5(0, -1,0,-1,-1)
     $  - HR5(0, -1,0,-1,1)
     $  - HR5(0, -1,0,1,-1)
     $  + HR5(0, -1,1,-1,-1)
     $  + 7.0000000000000000d+00*HR5(0,0,-1,-1,-1)
     $  + HR5(0,0, -1,-1,1)
     $  - HR5(0,0, -1,0,-1)
     $  - HR5(0,0, -1,0,1)
     $  + HR5(0,1, -1,-1,-1)
      HY5(0,1,0,1,-1) = 
     $  + 1.4122347902560834d-01
     $  + 1.4132080497842155d+00*HR1(-1)
     $  - 2.5410760640234242d-01*HR1(-1)*HR1(-1)
     $  - 1.7711559006386898d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 2.8881132523331054d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 8.3333333333333333d-03*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  - 5.3134677019160696d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR2(-1,1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR2(0,-1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR2(0,1)
     $  - HR1( -1)*HR1(-1)*HR3(-1,-1,1)
     $  + HR1( -1)*HR1(-1)*HR3(-1,1,1)
     $  - HR1( -1)*HR1(-1)*HR3(0,-1,-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(0,-1,1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(0,1,-1)
     $  + 5.0821521280468485d-01*HR1(-1)*HR1(0)
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(0)*HR2(-1,1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(0)*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(0)*HR2(0,1)
     $  - 5.0821521280468485d-01*HR1(-1)*HR1(1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(1)*HR2(-1,1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(1)*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(1)*HR2(0,1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR1(1)*HR3(-1,-1,1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR1(1)*HR3(0,-1,-1)
     $  - HR1( -1)*HR1(1)*HR3(0,-1,1)
     $  - HR1( -1)*HR1(1)*HR3(0,1,-1)
     $  + 1.0626935403832139d+00*HR1(-1)*HR2(-1,1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR2(-1,1)*HR2(-1,1)
     $  - HR1( -1)*HR2(-1,1)*HR2(0,-1)
     $  + 1.0626935403832139d+00*HR1(-1)*HR2(0,-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR2(0,-1)*HR2(0,-1)
     $  - HR1( -1)*HR2(0,-1)*HR2(0,1)
     $  + 1.0626935403832139d+00*HR1(-1)*HR2(0,1)
     $  + 1.3862943611198906d+00*HR1(-1)*HR3(-1,-1,1)
     $  - 1.3862943611198906d+00*HR1(-1)*HR3(-1,1,1)
     $  - 1.3862943611198906d+00*HR1(-1)*HR3(0,-1,-1)
     $  - 1.3862943611198906d+00*HR1(-1)*HR3(0,0,-1)
     $  - 1.3862943611198906d+00*HR1(-1)*HR3(0,0,1)
     $  - 1.3862943611198906d+00*HR1(-1)*HR3(0,1,-1)
     $  + 4.0000000000000000d+00*HR1(-1)*HR4(-1,-1,-1,1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR4(-1,-1,1,1)
     $  + 4.0000000000000000d+00*HR1(-1)*HR4(0,-1,-1,-1)
     $  + HR1( -1)*HR4(0,-1,0,1)
     $  + HR1( -1)*HR4(0,-1,1,-1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR4(0,0,-1,-1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR4(0,0,-1,1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR4(0,0,1,-1)
     $  + 3.0000000000000000d+00*HR1(-1)*HR4(0,1,-1,-1)
     $  + 5.0821521280468485d-01*HR1(0)*HR1(1)
     $  - 6.9314718055994530d-01*HR1(0)*HR1(1)*HR2(-1,1)
     $  + 6.9314718055994530d-01*HR1(0)*HR1(1)*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR1(0)*HR1(1)*HR2(0,1)
     $  - 1.3862943611198906d+00*HR1(0)*HR3(-1,-1,1)
     $  + 1.3862943611198906d+00*HR1(0)*HR3(-1,1,1)
     $  - 1.3862943611198906d+00*HR1(0)*HR3(0,-1,-1)
     $  - 1.3862943611198906d+00*HR1(0)*HR3(0,-1,1)
     $  - 1.3862943611198906d+00*HR1(0)*HR3(0,1,-1)
     $  - 1.3862943611198906d+00*HR1(0)*HR3(0,1,1)
     $  + 1.4132080497842155d+00*HR1(1)
     $  - 1.0626935403832139d+00*HR1(1)*HR2(-1,1)
     $  + HR1(1) *HR2(-1,1)*HR2(0,-1)
     $  + 1.0626935403832139d+00*HR1(1)*HR2(0,-1)
     $  - 5.0000000000000000d-01*HR1(1)*HR2(0,-1)*HR2(0,-1)
     $  - HR1(1) *HR2(0,-1)*HR2(0,1)
     $  + 1.0626935403832139d+00*HR1(1)*HR2(0,1)
     $  - 1.3862943611198906d+00*HR1(1)*HR3(-1,-1,1)
     $  - 1.3862943611198906d+00*HR1(1)*HR3(0,-1,-1)
     $  - 1.3862943611198906d+00*HR1(1)*HR3(0,0,-1)
     $  - 1.3862943611198906d+00*HR1(1)*HR3(0,0,1)
     $  - 1.3862943611198906d+00*HR1(1)*HR3(0,1,-1)
     $  - 3.0000000000000000d+00*HR1(1)*HR4(-1,-1,-1,1)
     $  + 4.0000000000000000d+00*HR1(1)*HR4(0,-1,-1,-1)
     $  + HR1(1) *HR4(0,-1,0,1)
     $  + HR1(1) *HR4(0,-1,1,-1)
     $  + 2.0000000000000000d+00*HR1(1)*HR4(0,0,-1,-1)
     $  + 2.0000000000000000d+00*HR1(1)*HR4(0,0,-1,1)
     $  + 2.0000000000000000d+00*HR1(1)*HR4(0,0,1,-1)
     $  + 3.0000000000000000d+00*HR1(1)*HR4(0,1,-1,-1)
     $  + 5.0821521280468485d-01*HR2(-1,1)
     $  + 3.4657359027997265d-01*HR2(-1,1)*HR2(-1,1)
     $  - 6.9314718055994530d-01*HR2(-1,1)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR2(-1,1)*HR2(0,1)
     $  + HR2( -1,1)*HR3(-1,-1,1)
     $  + 2.0000000000000000d+00*HR2(-1,1)*HR3(0,-1,-1)
     $  + HR2( -1,1)*HR3(0,-1,1)
     $  + HR2( -1,1)*HR3(0,1,-1)
     $  - 5.0821521280468485d-01*HR2(0,-1)
     $  + 3.4657359027997265d-01*HR2(0,-1)*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR2(0,-1)*HR2(0,1)
     $  + 2.0000000000000000d+00*HR2(0,-1)*HR3(-1,-1,1)
     $  - 2.0000000000000000d+00*HR2(0,-1)*HR3(-1,1,1)
     $  + 2.0000000000000000d+00*HR2(0,-1)*HR3(0,-1,-1)
     $  + 2.0000000000000000d+00*HR2(0,-1)*HR3(0,-1,1)
     $  + 2.0000000000000000d+00*HR2(0,-1)*HR3(0,1,-1)
     $  + 2.0000000000000000d+00*HR2(0,-1)*HR3(0,1,1)
     $  - 5.0821521280468485d-01*HR2(0,1)
     $  + 3.4657359027997265d-01*HR2(0,1)*HR2(0,1)
     $  - 2.0000000000000000d+00*HR2(0,1)*HR3(0,-1,-1)
     $  - HR2(0,1) *HR3(0,-1,1)
     $  - 2.1253870807664278d+00*HR3(-1,-1,1)
     $  + 2.1253870807664278d+00*HR3(-1,1,1)
     $  - 2.1253870807664278d+00*HR3(0,-1,-1)
     $  - 2.1253870807664278d+00*HR3(0,-1,1)
     $  - 2.1253870807664278d+00*HR3(0,1,-1)
     $  - 2.1253870807664278d+00*HR3(0,1,1)
     $  - 2.7725887222397812d+00*HR4(-1,-1,-1,1)
     $  + 1.3862943611198906d+00*HR4(-1,-1,1,1)
     $  + 2.7725887222397812d+00*HR4(0,-1,-1,-1)
     $  + 1.3862943611198906d+00*HR4(0,-1,-1,1)
     $  + 1.3862943611198906d+00*HR4(0,-1,1,-1)
     $  + 1.3862943611198906d+00*HR4(0,0,-1,-1)
     $  + 1.3862943611198906d+00*HR4(0,0,-1,1)
     $  + 1.3862943611198906d+00*HR4(0,0,1,-1)
     $  + 1.3862943611198906d+00*HR4(0,0,1,1)
     $  + 2.7725887222397812d+00*HR4(0,1,-1,-1)
     $  + 1.3862943611198906d+00*HR4(0,1,-1,1)
     $  + 1.3862943611198906d+00*HR4(0,1,1,-1)
     $  - 7.0000000000000000d+00*HR5(-1,-1,-1,-1,1)
     $  - HR5( -1,-1,1,-1,1)
     $  - 8.0000000000000000d+00*HR5(0,-1,-1,-1,-1)
     $  - 4.0000000000000000d+00*HR5(0,-1,-1,-1,1)
     $  + 2.0000000000000000d+00*HR5(0,-1,-1,0,1)
     $  - 2.0000000000000000d+00*HR5(0,-1,-1,1,-1)
     $  - 4.0000000000000000d+00*HR5(0,-1,0,-1,-1)
     $  - HR5(0, -1,0,-1,1)
     $  - HR5(0, -1,0,1,-1)
     $  - 3.0000000000000000d+00*HR5(0,-1,1,-1,-1)
     $  - HR5(0, -1,1,-1,1)
     $  + HR5(0, -1,1,0,1)
     $  - 1.2000000000000000d+01*HR5(0,0,-1,-1,-1)
     $  - 4.0000000000000000d+00*HR5(0,0,-1,-1,1)
     $  - 4.0000000000000000d+00*HR5(0,0,-1,1,-1)
     $  - 4.0000000000000000d+00*HR5(0,0,1,-1,-1)
     $  - 2.0000000000000000d+00*HR5(0,0,1,-1,1)
     $  - 4.0000000000000000d+00*HR5(0,0,1,1,-1)
     $  - 7.0000000000000000d+00*HR5(0,1,-1,-1,-1)
     $  - 3.0000000000000000d+00*HR5(0,1,-1,-1,1)
     $  - 2.0000000000000000d+00*HR5(0,1,-1,1,-1)
     $  - HR5(0,1,0,1, -1)
     $  - 2.0000000000000000d+00*HR5(0,1,1,-1,-1)
      HY5(0,1,1,-1,-1) = 
     $  + 1.0254618242743082d-01
     $  - 6.0337978402921669d-01*HR1(-1)
     $  + 2.3766885054564933d-01*HR1(-1)*HR1(-1)
     $  - 4.0037751159850118d-02*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 2.8881132523331054d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 8.3333333333333333d-03*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  + 1.2011325347955035d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  - 1.2011325347955035d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(-1,-1,1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(0,-1,-1)
     $  - 4.7533770109129867d-01*HR1(-1)*HR1(0)
     $  - 1.2011325347955035d-01*HR1(-1)*HR1(0)*HR1(0)
     $  + 2.4022650695910071d-01*HR1(-1)*HR1(0)*HR1(1)
     $  + 4.7533770109129867d-01*HR1(-1)*HR1(1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(1)*HR2(0,-1)
     $  + HR1( -1)*HR1(1)*HR3(0,-1,-1)
     $  + 2.4022650695910071d-01*HR1(-1)*HR2(-1,1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR3(-1,-1,1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR3(0,-1,-1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR3(0,0,-1)
     $  + HR1( -1)*HR4(-1,-1,-1,1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR4(0,-1,-1,-1)
     $  - HR1( -1)*HR4(0,0,-1,-1)
     $  - 1.2011325347955035d-01*HR1(0)*HR1(0)*HR1(1)
     $  - 4.7533770109129867d-01*HR1(0)*HR1(1)
     $  - 2.4022650695910071d-01*HR1(0)*HR2(-1,1)
     $  + 2.4022650695910071d-01*HR1(0)*HR2(0,-1)
     $  + 2.4022650695910071d-01*HR1(0)*HR2(0,1)
     $  - 6.0337978402921669d-01*HR1(1)
     $  + 6.9314718055994530d-01*HR1(1)*HR3(0,-1,-1)
     $  + 6.9314718055994530d-01*HR1(1)*HR3(0,0,-1)
     $  - 2.0000000000000000d+00*HR1(1)*HR4(0,-1,-1,-1)
     $  - HR1(1) *HR4(0,0,-1,-1)
     $  - 4.7533770109129867d-01*HR2(-1,1)
     $  + 6.9314718055994530d-01*HR2(-1,1)*HR2(0,-1)
     $  - HR2( -1,1)*HR3(0,-1,-1)
     $  + 4.7533770109129867d-01*HR2(0,-1)
     $  - 3.4657359027997265d-01*HR2(0,-1)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR2(0,-1)*HR2(0,1)
     $  + 4.7533770109129867d-01*HR2(0,1)
     $  + HR2(0,1) *HR3(0,-1,-1)
     $  - 2.4022650695910071d-01*HR3(-1,-1,1)
     $  - 2.4022650695910071d-01*HR3(0,-1,-1)
     $  - 2.4022650695910071d-01*HR3(0,0,-1)
     $  - 2.4022650695910071d-01*HR3(0,0,1)
     $  - 2.4022650695910071d-01*HR3(0,1,-1)
     $  - 6.9314718055994530d-01*HR4(-1,-1,-1,1)
     $  - 6.9314718055994530d-01*HR4(0,-1,-1,1)
     $  + 6.9314718055994530d-01*HR4(0,-1,0,1)
     $  + 6.9314718055994530d-01*HR4(0,0,-1,1)
     $  + 6.9314718055994530d-01*HR4(0,0,1,-1)
     $  + 6.9314718055994530d-01*HR4(0,1,-1,-1)
     $  - HR5( -1,-1,-1,-1,1)
     $  + 2.0000000000000000d+00*HR5(0,-1,-1,-1,-1)
     $  + 2.0000000000000000d+00*HR5(0,-1,-1,-1,1)
     $  - HR5(0, -1,-1,0,1)
     $  + HR5(0, -1,-1,1,-1)
     $  + HR5(0, -1,0,-1,-1)
     $  - HR5(0, -1,0,-1,1)
     $  - HR5(0, -1,0,1,-1)
     $  + 3.0000000000000000d+00*HR5(0,0,-1,-1,-1)
     $  - HR5(0,0, -1,-1,1)
     $  - HR5(0,0, -1,1,-1)
     $  - HR5(0,0,1, -1,-1)
     $  - HR5(0,1, -1,-1,-1)
      HY5(0,1,1,-1,1) = 
     $  + 1.7692012816527167d-01
     $  - 2.4370424139224501d+00*HR1(-1)
     $  + 7.3900238327152102d-01*HR1(-1)*HR1(-1)
     $  - 9.7040087744168750d-02*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 2.8881132523331054d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 8.3333333333333333d-03*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  + 2.9112026323250625d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR2(-1,1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR2(0,-1)
     $  - 2.9112026323250625d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR2(0,-1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(-1,-1,1)
     $  - HR1( -1)*HR1(-1)*HR3(0,-1,-1)
     $  - HR1( -1)*HR1(-1)*HR3(0,0,-1)
     $  - 1.4780047665430420d+00*HR1(-1)*HR1(0)
     $  - 2.9112026323250625d-01*HR1(-1)*HR1(0)*HR1(0)
     $  + 5.8224052646501250d-01*HR1(-1)*HR1(0)*HR1(1)
     $  + HR1( -1)*HR1(0)*HR1(1)*HR2(0,-1)
     $  - HR1( -1)*HR1(0)*HR3(-1,-1,1)
     $  - HR1( -1)*HR1(0)*HR3(0,-1,-1)
     $  - HR1( -1)*HR1(0)*HR3(0,0,-1)
     $  + 1.4780047665430420d+00*HR1(-1)*HR1(1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(1)*HR2(0,-1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR1(1)*HR3(0,-1,-1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR1(1)*HR3(0,0,-1)
     $  + 5.8224052646501250d-01*HR1(-1)*HR2(-1,1)
     $  - HR1( -1)*HR2(-1,1)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR3(-1,-1,1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR3(0,-1,-1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR3(0,0,-1)
     $  - HR1( -1)*HR4(-1,-1,-1,1)
     $  + 3.0000000000000000d+00*HR1(-1)*HR4(0,-1,-1,-1)
     $  + 3.0000000000000000d+00*HR1(-1)*HR4(0,0,-1,-1)
     $  + 3.0000000000000000d+00*HR1(-1)*HR4(0,0,0,-1)
     $  - 2.9112026323250625d-01*HR1(0)*HR1(0)*HR1(1)
     $  - 1.4780047665430420d+00*HR1(0)*HR1(1)
     $  - HR1(0) *HR1(1)*HR3(0,-1,-1)
     $  - HR1(0) *HR1(1)*HR3(0,0,-1)
     $  - 5.8224052646501250d-01*HR1(0)*HR2(-1,1)
     $  - HR1(0) *HR2(-1,1)*HR2(0,-1)
     $  + 5.8224052646501250d-01*HR1(0)*HR2(0,-1)
     $  + 5.0000000000000000d-01*HR1(0)*HR2(0,-1)*HR2(0,-1)
     $  + HR1(0) *HR2(0,-1)*HR2(0,1)
     $  + 5.8224052646501250d-01*HR1(0)*HR2(0,1)
     $  + HR1(0) *HR4(-1,-1,-1,1)
     $  + HR1(0) *HR4(0,-1,-1,1)
     $  - HR1(0) *HR4(0,-1,0,1)
     $  - HR1(0) *HR4(0,0,-1,1)
     $  - HR1(0) *HR4(0,0,1,-1)
     $  - HR1(0) *HR4(0,1,-1,-1)
     $  - 2.4370424139224501d+00*HR1(1)
     $  - 6.9314718055994530d-01*HR1(1)*HR3(0,-1,-1)
     $  - 6.9314718055994530d-01*HR1(1)*HR3(0,0,-1)
     $  + 3.0000000000000000d+00*HR1(1)*HR4(0,-1,-1,-1)
     $  + 3.0000000000000000d+00*HR1(1)*HR4(0,0,-1,-1)
     $  + 3.0000000000000000d+00*HR1(1)*HR4(0,0,0,-1)
     $  - 1.4780047665430420d+00*HR2(-1,1)
     $  - 6.9314718055994530d-01*HR2(-1,1)*HR2(0,-1)
     $  + 2.0000000000000000d+00*HR2(-1,1)*HR3(0,-1,-1)
     $  + 2.0000000000000000d+00*HR2(-1,1)*HR3(0,0,-1)
     $  + 1.4780047665430420d+00*HR2(0,-1)
     $  + 3.4657359027997265d-01*HR2(0,-1)*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR2(0,-1)*HR2(0,1)
     $  + HR2(0, -1)*HR3(-1,-1,1)
     $  + HR2(0, -1)*HR3(0,-1,-1)
     $  - 2.0000000000000000d+00*HR2(0,-1)*HR3(0,0,-1)
     $  + HR2(0, -1)*HR3(0,1,-1)
     $  + 1.4780047665430420d+00*HR2(0,1)
     $  - 2.0000000000000000d+00*HR2(0,1)*HR3(0,-1,-1)
     $  - 2.0000000000000000d+00*HR2(0,1)*HR3(0,0,-1)
     $  - 5.8224052646501250d-01*HR3(-1,-1,1)
     $  - 5.8224052646501250d-01*HR3(0,-1,-1)
     $  - 5.8224052646501250d-01*HR3(0,0,-1)
     $  - 5.8224052646501250d-01*HR3(0,0,1)
     $  - 5.8224052646501250d-01*HR3(0,1,-1)
     $  + 6.9314718055994530d-01*HR4(-1,-1,-1,1)
     $  + 6.9314718055994530d-01*HR4(0,-1,-1,1)
     $  - 6.9314718055994530d-01*HR4(0,-1,0,1)
     $  - 6.9314718055994530d-01*HR4(0,0,-1,1)
     $  - 6.9314718055994530d-01*HR4(0,0,1,-1)
     $  - 6.9314718055994530d-01*HR4(0,1,-1,-1)
     $  + HR5( -1,-1,-1,-1,1)
     $  - 3.0000000000000000d+00*HR5(0,-1,-1,-1,-1)
     $  - 3.0000000000000000d+00*HR5(0,-1,-1,-1,1)
     $  + 2.0000000000000000d+00*HR5(0,-1,-1,0,1)
     $  - HR5(0, -1,-1,1,-1)
     $  - 3.0000000000000000d+00*HR5(0,-1,0,-1,-1)
     $  + 2.0000000000000000d+00*HR5(0,-1,0,-1,1)
     $  + HR5(0, -1,0,1,-1)
     $  - 9.0000000000000000d+00*HR5(0,0,-1,-1,-1)
     $  + HR5(0,0, -1,-1,1)
     $  + 3.0000000000000000d+00*HR5(0,0,-1,0,-1)
     $  + 2.0000000000000000d+00*HR5(0,0,-1,0,1)
     $  + HR5(0,0, -1,1,-1)
     $  + 6.0000000000000000d+00*HR5(0,0,0,-1,-1)
     $  + 3.0000000000000000d+00*HR5(0,0,0,-1,1)
     $  + 3.0000000000000000d+00*HR5(0,0,0,1,-1)
     $  + HR5(0,0,1, -1,-1)
     $  + HR5(0,0,1,0, -1)
     $  + HR5(0,1, -1,-1,-1)
      HY5(0,1,1,1,-1) = 
     $  + 4.0707197178331534d-01
     $  + 1.0681889964961481d+00*HR1(-1)
     $  - 5.5365194946473328d-01*HR1(-1)*HR1(-1)
     $  + 1.7711559006386898d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 2.8881132523331054d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 8.3333333333333333d-03*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  - 5.3134677019160696d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  - 1.7328679513998632d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(0)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  + 5.3134677019160696d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR2(0,-1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(-1,-1,1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(0,-1,-1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR3(0,0,-1)
     $  + 1.1073038989294665d+00*HR1(-1)*HR1(0)
     $  + 5.3134677019160696d-01*HR1(-1)*HR1(0)*HR1(0)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(0)*HR1(0)*HR1(0)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(0)*HR1(0)*HR1(1)
     $  - 1.0626935403832139d+00*HR1(-1)*HR1(0)*HR1(1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(0)*HR2(-1,1)
     $  - 1.1073038989294665d+00*HR1(-1)*HR1(1)
     $  + HR1( -1)*HR1(1)*HR3(0,-1,-1)
     $  + HR1( -1)*HR1(1)*HR3(0,0,-1)
     $  - 1.0626935403832139d+00*HR1(-1)*HR2(-1,1)
     $  + HR1( -1)*HR2(-1,1)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR3(-1,-1,1)
     $  - HR1( -1)*HR4(-1,-1,-1,1)
     $  - HR1( -1)*HR4(0,-1,-1,-1)
     $  - HR1( -1)*HR4(0,0,-1,-1)
     $  - HR1( -1)*HR4(0,0,0,-1)
     $  + 1.1552453009332421d-01*HR1(0)*HR1(0)*HR1(0)*HR1(1)
     $  + 5.3134677019160696d-01*HR1(0)*HR1(0)*HR1(1)
     $  + 3.4657359027997265d-01*HR1(0)*HR1(0)*HR2(-1,1)
     $  - 3.4657359027997265d-01*HR1(0)*HR1(0)*HR2(0,-1)
     $  - 3.4657359027997265d-01*HR1(0)*HR1(0)*HR2(0,1)
     $  + 1.1073038989294665d+00*HR1(0)*HR1(1)
     $  + 1.0626935403832139d+00*HR1(0)*HR2(-1,1)
     $  - 1.0626935403832139d+00*HR1(0)*HR2(0,-1)
     $  - 1.0626935403832139d+00*HR1(0)*HR2(0,1)
     $  + 6.9314718055994530d-01*HR1(0)*HR3(-1,-1,1)
     $  + 6.9314718055994530d-01*HR1(0)*HR3(0,-1,-1)
     $  + 6.9314718055994530d-01*HR1(0)*HR3(0,0,-1)
     $  + 6.9314718055994530d-01*HR1(0)*HR3(0,0,1)
     $  + 6.9314718055994530d-01*HR1(0)*HR3(0,1,-1)
     $  + 1.0681889964961481d+00*HR1(1)
     $  - HR1(1) *HR4(0,-1,-1,-1)
     $  - HR1(1) *HR4(0,0,-1,-1)
     $  - HR1(1) *HR4(0,0,0,-1)
     $  + 1.1073038989294665d+00*HR2(-1,1)
     $  - HR2( -1,1)*HR3(0,-1,-1)
     $  - HR2( -1,1)*HR3(0,0,-1)
     $  - 1.1073038989294665d+00*HR2(0,-1)
     $  - HR2(0, -1)*HR3(-1,-1,1)
     $  - HR2(0, -1)*HR3(0,-1,-1)
     $  + HR2(0, -1)*HR3(0,0,-1)
     $  - HR2(0, -1)*HR3(0,1,-1)
     $  - 1.1073038989294665d+00*HR2(0,1)
     $  + HR2(0,1) *HR3(0,-1,-1)
     $  + HR2(0,1) *HR3(0,0,-1)
     $  + 1.0626935403832139d+00*HR3(-1,-1,1)
     $  + 1.0626935403832139d+00*HR3(0,-1,-1)
     $  + 1.0626935403832139d+00*HR3(0,0,-1)
     $  + 1.0626935403832139d+00*HR3(0,0,1)
     $  + 1.0626935403832139d+00*HR3(0,1,-1)
     $  + 6.9314718055994530d-01*HR4(-1,-1,-1,1)
     $  - 6.9314718055994530d-01*HR4(0,-1,-1,-1)
     $  - 6.9314718055994530d-01*HR4(0,0,-1,-1)
     $  - 6.9314718055994530d-01*HR4(0,0,0,-1)
     $  - 6.9314718055994530d-01*HR4(0,0,0,1)
     $  - 6.9314718055994530d-01*HR4(0,0,1,-1)
     $  - 6.9314718055994530d-01*HR4(0,1,-1,-1)
     $  + HR5( -1,-1,-1,-1,1)
     $  + 2.0000000000000000d+00*HR5(0,-1,-1,-1,-1)
     $  + HR5(0, -1,-1,-1,1)
     $  - HR5(0, -1,-1,0,1)
     $  + 2.0000000000000000d+00*HR5(0,-1,0,-1,-1)
     $  - HR5(0, -1,0,-1,1)
     $  + 6.0000000000000000d+00*HR5(0,0,-1,-1,-1)
     $  - HR5(0,0, -1,-1,1)
     $  - 2.0000000000000000d+00*HR5(0,0,-1,0,-1)
     $  - HR5(0,0, -1,0,1)
     $  - 4.0000000000000000d+00*HR5(0,0,0,-1,-1)
     $  - 2.0000000000000000d+00*HR5(0,0,0,-1,1)
     $  - 2.0000000000000000d+00*HR5(0,0,0,1,-1)
     $  + HR5(0,0,1, -1,-1)
     $  - HR5(0,0,1,0, -1)
     $  + HR5(0,1, -1,-1,-1)
      if (r.lt.0d0) then 
      HY5(-1,-1,-1,1,1) = HY5(-1,-1,-1,1,1) 
     $  + 8.2246703342411321d-01*HR1(-1)*HR1(-1)*HR1(-1)
      HY5(-1,-1,1,1,1) = HY5(-1,-1,1,1,1) 
     $  + 1.7102721159642791d+00*HR1(-1)*HR1(-1)
     $  - 8.2246703342411321d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 2.4674011002723396d+00*HR1(-1)*HR1(-1)*HR1(0)
     $  - 4.9348022005446793d+00*HR1(-1)*HR2(0,-1)
     $  + 4.9348022005446793d+00*HR3(0,-1,-1)
      HY5(-1,1,-1,1,1) = HY5(-1,1,-1,1,1) 
     $  - 8.2246703342411321d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 4.9348022005446793d+00*HR1(-1)*HR2(0,-1)
     $  - 9.8696044010893586d+00*HR3(0,-1,-1)
      HY5(-1,1,1,1,1) = HY5(-1,1,1,1,1) 
     $  - 2.8732418312458363d+00*HR1(-1)
     $  - 1.7102721159642791d+00*HR1(-1)*HR1(-1)
     $  + 8.2246703342411321d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 2.4674011002723396d+00*HR1(-1)*HR1(-1)*HR1(0)
     $  + 3.4205442319285582d+00*HR1(-1)*HR1(0)
     $  + 2.4674011002723396d+00*HR1(-1)*HR1(0)*HR1(0)
     $  - 4.9348022005446793d+00*HR1(0)*HR2(0,-1)
     $  - 3.4205442319285582d+00*HR2(0,-1)
     $  + 4.9348022005446793d+00*HR3(0,-1,-1)
     $  + 4.9348022005446793d+00*HR3(0,0,-1)
      HY5(0,-1,-1,1,1) = HY5(0,-1,-1,1,1) 
     $  + 8.2246703342411321d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 2.4674011002723396d+00*HR1(-1)*HR1(-1)*HR1(1)
     $  - 4.9348022005446793d+00*HR1(-1)*HR2(-1,1)
     $  + 4.9348022005446793d+00*HR3(-1,-1,1)
      HY5(0,-1,0,1,1) = HY5(0,-1,0,1,1) 
     $  + 8.2246703342411321d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 2.4674011002723396d+00*HR1(-1)*HR1(-1)*HR1(1)
     $  - 4.9348022005446793d+00*HR1(-1)*HR2(-1,1)
     $  + 4.9348022005446793d+00*HR1(1)*HR2(-1,1)
     $  + 9.8696044010893586d+00*HR3(-1,-1,1)
     $  - 9.8696044010893586d+00*HR3(-1,1,1)
      HY5(0,-1,1,1,-1) = HY5(0,-1,1,1,-1) 
     $  - 1.7102721159642791d+00*HR1(-1)*HR1(-1)
     $  - 3.4205442319285582d+00*HR1(-1)*HR1(1)
     $  + 3.4205442319285582d+00*HR2(-1,1)
      HY5(0,-1,1,1,1) = HY5(0,-1,1,1,1) 
     $  + 1.7102721159642791d+00*HR1(-1)*HR1(-1)
     $  - 8.2246703342411321d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 2.4674011002723396d+00*HR1(-1)*HR1(-1)*HR1(0)
     $  - 2.4674011002723396d+00*HR1(-1)*HR1(-1)*HR1(1)
     $  + 4.9348022005446793d+00*HR1(-1)*HR1(0)*HR1(1)
     $  + 3.4205442319285582d+00*HR1(-1)*HR1(1)
     $  + 4.9348022005446793d+00*HR1(-1)*HR2(-1,1)
     $  - 4.9348022005446793d+00*HR1(-1)*HR2(0,-1)
     $  - 4.9348022005446793d+00*HR1(0)*HR2(-1,1)
     $  - 4.9348022005446793d+00*HR1(1)*HR2(0,-1)
     $  - 3.4205442319285582d+00*HR2(-1,1)
     $  - 4.9348022005446793d+00*HR3(-1,-1,1)
     $  + 4.9348022005446793d+00*HR3(0,-1,-1)
     $  + 4.9348022005446793d+00*HR3(0,-1,1)
      HY5(0,0,-1,1,1) = HY5(0,0,-1,1,1) 
     $  + 8.2246703342411321d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 2.4674011002723396d+00*HR1(-1)*HR1(-1)*HR1(1)
     $  + 2.4674011002723396d+00*HR1(-1)*HR1(1)*HR1(1)
     $  - 4.9348022005446793d+00*HR1(1)*HR2(-1,1)
     $  - 4.9348022005446793d+00*HR3(-1,-1,1)
     $  + 4.9348022005446793d+00*HR3(-1,1,1)
      HY5(0,0,1,1,-1) = HY5(0,0,1,1,-1) 
     $  - 1.7102721159642791d+00*HR1(-1)*HR1(-1)
     $  - 3.4205442319285582d+00*HR1(-1)*HR1(1)
     $  - 1.7102721159642791d+00*HR1(1)*HR1(1)
      HY5(0,1,-1,1,1) = HY5(0,1,-1,1,1) 
     $  - 8.2246703342411321d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 2.4674011002723396d+00*HR1(-1)*HR1(-1)*HR1(1)
     $  + 4.9348022005446793d+00*HR1(-1)*HR2(-1,1)
     $  + 4.9348022005446793d+00*HR1(-1)*HR2(0,-1)
     $  + 4.9348022005446793d+00*HR1(1)*HR2(0,-1)
     $  - 4.9348022005446793d+00*HR3(-1,-1,1)
     $  - 9.8696044010893586d+00*HR3(0,-1,-1)
     $  - 4.9348022005446793d+00*HR3(0,-1,1)
     $  - 4.9348022005446793d+00*HR3(0,1,-1)
      HY5(0,1,1,-1,-1) = HY5(0,1,1,-1,-1) 
     $  + 1.1854702951709319d+00*HR1(-1)
     $  + 1.1854702951709319d+00*HR1(1)
      HY5(0,1,1,-1,1) = HY5(0,1,1,-1,1) 
     $  + 2.8732418312458363d+00*HR1(-1)
     $  + 2.8732418312458363d+00*HR1(1)
      HY5(0,1,1,1,-1) = HY5(0,1,1,1,-1) 
     $  - 5.2441824215877001d+00*HR1(-1)
     $  + 1.7102721159642791d+00*HR1(-1)*HR1(-1)
     $  - 3.4205442319285582d+00*HR1(-1)*HR1(0)
     $  + 3.4205442319285582d+00*HR1(-1)*HR1(1)
     $  - 3.4205442319285582d+00*HR1(0)*HR1(1)
     $  - 5.2441824215877001d+00*HR1(1)
     $  - 3.4205442319285582d+00*HR2(-1,1)
     $  + 3.4205442319285582d+00*HR2(0,-1)
     $  + 3.4205442319285582d+00*HR2(0,1)
      Hi5(-1,-1,-1,-1,1) = 
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
      Hi5(-1,-1,-1,1,1) = 
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  + HR1( -1)*HR3(0,-1,-1)
     $  - HR4(0, -1,-1,-1)
      Hi5(-1,-1,1,-1,1) = 
     $  + 2.9112026323250625d-01*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR3(0,-1,-1)
     $  + 3.0000000000000000d+00*HR4(0,-1,-1,-1)
      Hi5(-1,-1,1,1,1) = 
     $  - 7.0235377994456286d-01*HR1(-1)*HR1(-1)
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
      Hi5(-1,1,-1,1,1) = 
     $  - 5.3721319360804020d-01*HR1(-1)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  + HR1( -1)*HR1(0)*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR2(0,-1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR3(0,-1,-1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR3(0,0,-1)
     $  - 2.0000000000000000d+00*HR1(0)*HR3(0,-1,-1)
     $  + 5.0000000000000000d-01*HR2(0,-1)*HR2(0,-1)
     $  - 1.3862943611198906d+00*HR3(0,-1,-1)
     $  + 4.0000000000000000d+00*HR4(0,-1,-1,-1)
     $  + 2.0000000000000000d+00*HR4(0,0,-1,-1)
      Hi5(-1,1,1,1,1) = 
     $  - 1.0846773019780311d+00*HR1(-1)
     $  + 7.0235377994456286d-01*HR1(-1)*HR1(-1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  - 2.5000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(0)
     $  - 1.4047075598891257d+00*HR1(-1)*HR1(0)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(0)*HR1(0)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(0)*HR1(0)*HR1(0)
     $  - 5.0000000000000000d-01*HR1(0)*HR1(0)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(0)*HR2(0,-1)
     $  + HR1(0) *HR3(0,-1,-1)
     $  + HR1(0) *HR3(0,0,-1)
     $  + 1.4047075598891257d+00*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR3(0,-1,-1)
     $  + 6.9314718055994530d-01*HR3(0,0,-1)
     $  - HR4(0, -1,-1,-1)
     $  - HR4(0,0, -1,-1)
     $  - HR4(0,0,0, -1)
      Hi5(0,-1,-1,-1,1) = 
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  + HR1( -1)*HR3(-1,-1,1)
     $  - HR4( -1,-1,-1,1)
      Hi5(0,-1,-1,0,1) = 
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  + HR1( -1)*HR3(-1,-1,1)
     $  + HR1(1) *HR3(-1,-1,1)
     $  - 5.0000000000000000d-01*HR2(-1,1)*HR2(-1,1)
      Hi5(0,-1,-1,1,-1) = 
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR2(-1,1)
     $  - 6.9314718055994530d-01*HR3(-1,-1,1)
      Hi5(0,-1,-1,1,1) = 
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  - HR1( -1)*HR1(0)*HR2(-1,1)
     $  - HR1( -1)*HR1(1)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR2(-1,1)
     $  - HR1( -1)*HR3(-1,-1,1)
     $  + HR1( -1)*HR3(0,-1,-1)
     $  + HR1(0) *HR3(-1,-1,1)
     $  + HR1(1) *HR3(0,-1,-1)
     $  + HR2( -1,1)*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR3(-1,-1,1)
     $  + HR4( -1,-1,-1,1)
     $  - HR4(0, -1,-1,-1)
     $  - HR4(0, -1,-1,1)
      Hi5(0,-1,0,-1,1) = 
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  + HR1( -1)*HR1(1)*HR2(-1,1)
     $  + 2.0000000000000000d+00*HR1(-1)*HR3(-1,-1,1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR3(-1,1,1)
     $  - 2.0000000000000000d+00*HR1(1)*HR3(-1,-1,1)
     $  + 5.0000000000000000d-01*HR2(-1,1)*HR2(-1,1)
     $  - 4.0000000000000000d+00*HR4(-1,-1,-1,1)
     $  + 2.0000000000000000d+00*HR4(-1,-1,1,1)
      Hi5(0,-1,0,1,-1) = 
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR2(-1,1)
     $  - 6.9314718055994530d-01*HR1(1)*HR2(-1,1)
     $  - 1.3862943611198906d+00*HR3(-1,-1,1)
     $  + 1.3862943611198906d+00*HR3(-1,1,1)
      Hi5(0,-1,0,1,1) = 
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
     $  + 6.9314718055994530d-01*HR1(1)*HR2(-1,1)
     $  + 2.0000000000000000d+00*HR1(1)*HR3(-1,-1,1)
     $  + HR1(1) *HR3(0,-1,-1)
     $  + HR1(1) *HR3(0,1,-1)
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
      Hi5(0,-1,1,-1,-1) = 
     $  + 1.2011325347955035d-01*HR1(-1)*HR1(-1)
     $  + 2.4022650695910071d-01*HR1(-1)*HR1(1)
     $  - 2.4022650695910071d-01*HR2(-1,1)
      Hi5(0,-1,1,-1,1) = 
     $  + 2.9112026323250625d-01*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  + 5.8224052646501250d-01*HR1(-1)*HR1(1)
     $  + HR1( -1)*HR1(1)*HR2(0,-1)
     $  - HR1( -1)*HR3(-1,-1,1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR3(0,-1,-1)
     $  - 2.0000000000000000d+00*HR1(1)*HR3(0,-1,-1)
     $  - 5.8224052646501250d-01*HR2(-1,1)
     $  - HR2( -1,1)*HR2(0,-1)
     $  + HR4( -1,-1,-1,1)
     $  + 3.0000000000000000d+00*HR4(0,-1,-1,-1)
     $  + 2.0000000000000000d+00*HR4(0,-1,-1,1)
     $  + HR4(0, -1,1,-1)
      Hi5(0,-1,1,0,1) = 
     $  + 8.2246703342411321d-01*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(0,1)
     $  + 1.6449340668482264d+00*HR1(-1)*HR1(1)
     $  + HR1( -1)*HR1(1)*HR2(0,-1)
     $  + HR1( -1)*HR1(1)*HR2(0,1)
     $  - HR1( -1)*HR3(-1,-1,1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR3(0,-1,-1)
     $  - HR1( -1)*HR3(0,-1,1)
     $  - HR1( -1)*HR3(0,1,-1)
     $  - HR1(1) *HR3(-1,-1,1)
     $  - 2.0000000000000000d+00*HR1(1)*HR3(0,-1,-1)
     $  - HR1(1) *HR3(0,-1,1)
     $  - HR1(1) *HR3(0,1,-1)
     $  - 1.6449340668482264d+00*HR2(-1,1)
     $  + 5.0000000000000000d-01*HR2(-1,1)*HR2(-1,1)
     $  - HR2( -1,1)*HR2(0,-1)
     $  - HR2( -1,1)*HR2(0,1)
     $  + 3.0000000000000000d+00*HR4(0,-1,-1,-1)
     $  + 3.0000000000000000d+00*HR4(0,-1,-1,1)
     $  + 2.0000000000000000d+00*HR4(0,-1,1,-1)
     $  + 2.0000000000000000d+00*HR4(0,-1,1,1)
     $  + HR4(0,1, -1,-1)
     $  + HR4(0,1, -1,1)
      Hi5(0,-1,1,1,-1) = 
     $  - 5.3134677019160696d-01*HR1(-1)*HR1(-1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(0)*HR1(1)
     $  - 1.0626935403832139d+00*HR1(-1)*HR1(1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR2(-1,1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR1(0)*HR2(-1,1)
     $  + 6.9314718055994530d-01*HR1(1)*HR2(0,-1)
     $  + 1.0626935403832139d+00*HR2(-1,1)
     $  + 6.9314718055994530d-01*HR3(-1,-1,1)
     $  - 6.9314718055994530d-01*HR3(0,-1,-1)
     $  - 6.9314718055994530d-01*HR3(0,-1,1)
      Hi5(0,-1,1,1,1) = 
     $  - 7.0235377994456286d-01*HR1(-1)*HR1(-1)
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
     $  - 1.4047075598891257d+00*HR1(-1)*HR1(1)
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
     $  - 6.9314718055994530d-01*HR1(1)*HR2(0,-1)
     $  + HR1(1) *HR3(0,-1,-1)
     $  + HR1(1) *HR3(0,0,-1)
     $  + 1.4047075598891257d+00*HR2(-1,1)
     $  - 6.9314718055994530d-01*HR3(-1,-1,1)
     $  + 6.9314718055994530d-01*HR3(0,-1,-1)
     $  + 6.9314718055994530d-01*HR3(0,-1,1)
     $  - HR4( -1,-1,-1,1)
     $  - 2.0000000000000000d+00*HR4(0,-1,-1,-1)
     $  - HR4(0, -1,-1,1)
     $  - HR4(0, -1,1,-1)
     $  - HR4(0,0, -1,-1)
     $  - HR4(0,0, -1,1)
      Hi5(0,0,-1,-1,1) = 
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 2.5000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  - HR1( -1)*HR1(1)*HR2(-1,1)
     $  - HR1( -1)*HR3(-1,-1,1)
     $  + HR1( -1)*HR3(-1,1,1)
     $  + HR1(1) *HR3(-1,-1,1)
     $  + 2.0000000000000000d+00*HR4(-1,-1,-1,1)
     $  - HR4( -1,-1,1,1)
      Hi5(0,0,-1,0,1) = 
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 2.5000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  - HR1( -1)*HR1(1)*HR2(-1,1)
     $  - HR1( -1)*HR3(-1,-1,1)
     $  + HR1( -1)*HR3(-1,1,1)
     $  + 5.0000000000000000d-01*HR1(1)*HR1(1)*HR2(-1,1)
     $  + 2.0000000000000000d+00*HR1(1)*HR3(-1,-1,1)
     $  - 2.0000000000000000d+00*HR1(1)*HR3(-1,1,1)
     $  + 3.0000000000000000d+00*HR4(-1,-1,-1,1)
     $  - 3.0000000000000000d+00*HR4(-1,-1,1,1)
     $  + 3.0000000000000000d+00*HR4(-1,1,1,1)
      Hi5(0,0,-1,1,-1) = 
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(1)*HR1(1)
     $  + 6.9314718055994530d-01*HR1(1)*HR2(-1,1)
     $  + 6.9314718055994530d-01*HR3(-1,-1,1)
     $  - 6.9314718055994530d-01*HR3(-1,1,1)
      Hi5(0,0,-1,1,1) = 
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - 2.5000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(0)*HR1(1)*HR1(1)
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
      Hi5(0,0,0,-1,1) = 
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 2.5000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(1)*HR1(1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(1)*HR1(1)*HR2(-1,1)
     $  - HR1(1) *HR3(-1,-1,1)
     $  + HR1(1) *HR3(-1,1,1)
     $  - HR4( -1,-1,-1,1)
     $  + HR4( -1,-1,1,1)
     $  - HR4( -1,1,1,1)
      Hi5(0,0,0,1,-1) = 
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(1)*HR1(1)
     $  - 1.1552453009332421d-01*HR1(1)*HR1(1)*HR1(1)
      Hi5(0,0,1,-1,-1) = 
     $  + 1.2011325347955035d-01*HR1(-1)*HR1(-1)
     $  + 2.4022650695910071d-01*HR1(-1)*HR1(1)
     $  + 1.2011325347955035d-01*HR1(1)*HR1(1)
      Hi5(0,0,1,-1,1) = 
     $  + 2.9112026323250625d-01*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 2.5000000000000000d-01*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  + 5.8224052646501250d-01*HR1(-1)*HR1(1)
     $  + HR1( -1)*HR1(1)*HR2(-1,1)
     $  + HR1( -1)*HR1(1)*HR2(0,-1)
     $  + HR1( -1)*HR3(-1,-1,1)
     $  - HR1( -1)*HR3(-1,1,1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR3(0,-1,-1)
     $  - HR1( -1)*HR3(0,-1,1)
     $  - HR1( -1)*HR3(0,1,-1)
     $  + 2.9112026323250625d-01*HR1(1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(1)*HR1(1)*HR2(0,-1)
     $  - HR1(1) *HR3(-1,-1,1)
     $  - 2.0000000000000000d+00*HR1(1)*HR3(0,-1,-1)
     $  - HR1(1) *HR3(0,-1,1)
     $  - HR1(1) *HR3(0,1,-1)
     $  - 2.0000000000000000d+00*HR4(-1,-1,-1,1)
     $  + HR4( -1,-1,1,1)
     $  + 3.0000000000000000d+00*HR4(0,-1,-1,-1)
     $  + 2.0000000000000000d+00*HR4(0,-1,-1,1)
     $  + 2.0000000000000000d+00*HR4(0,-1,1,-1)
     $  + HR4(0, -1,1,1)
     $  + 2.0000000000000000d+00*HR4(0,1,-1,-1)
     $  + HR4(0,1, -1,1)
     $  + HR4(0,1,1, -1)
      Hi5(0,0,1,0,-1) = 
     $  + 4.1123351671205660d-01*HR1(-1)*HR1(-1)
     $  + 8.2246703342411321d-01*HR1(-1)*HR1(1)
     $  + 4.1123351671205660d-01*HR1(1)*HR1(1)
      Hi5(0,0,1,1,-1) = 
     $  - 5.3134677019160696d-01*HR1(-1)*HR1(-1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR1(0)*HR1(1)
     $  - 1.0626935403832139d+00*HR1(-1)*HR1(1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(1)*HR1(1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR2(0,1)
     $  - 3.4657359027997265d-01*HR1(0)*HR1(1)*HR1(1)
     $  - 5.3134677019160696d-01*HR1(1)*HR1(1)
     $  - 6.9314718055994530d-01*HR1(1)*HR2(-1,1)
     $  + 6.9314718055994530d-01*HR1(1)*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR1(1)*HR2(0,1)
     $  - 6.9314718055994530d-01*HR3(-1,-1,1)
     $  + 6.9314718055994530d-01*HR3(-1,1,1)
     $  - 6.9314718055994530d-01*HR3(0,-1,-1)
     $  - 6.9314718055994530d-01*HR3(0,-1,1)
     $  - 6.9314718055994530d-01*HR3(0,1,-1)
     $  - 6.9314718055994530d-01*HR3(0,1,1)
      Hi5(0,1,-1,-1,-1) = 
     $  - 5.5504108664821579d-02*HR1(-1)
     $  - 5.5504108664821579d-02*HR1(1)
      Hi5(0,1,-1,-1,1) = 
     $  - 9.4753004230127705d-02*HR1(-1)
     $  - 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  - HR1( -1)*HR3(-1,-1,1)
     $  + HR1( -1)*HR3(0,-1,-1)
     $  - 9.4753004230127705d-02*HR1(1)
     $  + HR1(1) *HR3(0,-1,-1)
     $  + HR4( -1,-1,-1,1)
     $  - 3.0000000000000000d+00*HR4(0,-1,-1,-1)
     $  - HR4(0, -1,-1,1)
     $  - HR4(0, -1,1,-1)
     $  - HR4(0,1, -1,-1)
      Hi5(0,1,-1,1,-1) = 
     $  - 2.1407237086670622d-01*HR1(-1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR2(-1,1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR2(0,-1)
     $  - 2.1407237086670622d-01*HR1(1)
     $  - 6.9314718055994530d-01*HR1(1)*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR3(-1,-1,1)
     $  + 1.3862943611198906d+00*HR3(0,-1,-1)
     $  + 6.9314718055994530d-01*HR3(0,-1,1)
     $  + 6.9314718055994530d-01*HR3(0,1,-1)
      Hi5(0,1,-1,1,1) = 
     $  - 5.3721319360804020d-01*HR1(-1)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  + 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  + HR1( -1)*HR1(0)*HR2(-1,1)
     $  + HR1( -1)*HR1(0)*HR2(0,-1)
     $  + HR1( -1)*HR1(1)*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR2(-1,1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR2(0,-1)
     $  + HR1( -1)*HR3(-1,-1,1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR3(0,-1,-1)
     $  - 2.0000000000000000d+00*HR1(-1)*HR3(0,0,-1)
     $  + HR1(0) *HR1(1)*HR2(0,-1)
     $  - HR1(0) *HR3(-1,-1,1)
     $  - 2.0000000000000000d+00*HR1(0)*HR3(0,-1,-1)
     $  - HR1(0) *HR3(0,-1,1)
     $  - HR1(0) *HR3(0,1,-1)
     $  - 5.3721319360804020d-01*HR1(1)
     $  + 6.9314718055994530d-01*HR1(1)*HR2(0,-1)
     $  - 2.0000000000000000d+00*HR1(1)*HR3(0,-1,-1)
     $  - 2.0000000000000000d+00*HR1(1)*HR3(0,0,-1)
     $  - HR2( -1,1)*HR2(0,-1)
     $  + 5.0000000000000000d-01*HR2(0,-1)*HR2(0,-1)
     $  + HR2(0, -1)*HR2(0,1)
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
      Hi5(0,1,0,1,-1) = 
     $  - 5.0821521280468485d-01*HR1(-1)
     $  + 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR2(-1,1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(-1)*HR2(0,1)
     $  - 5.0821521280468485d-01*HR1(1)
     $  + 6.9314718055994530d-01*HR1(1)*HR2(-1,1)
     $  - 6.9314718055994530d-01*HR1(1)*HR2(0,-1)
     $  - 6.9314718055994530d-01*HR1(1)*HR2(0,1)
     $  + 1.3862943611198906d+00*HR3(-1,-1,1)
     $  - 1.3862943611198906d+00*HR3(-1,1,1)
     $  + 1.3862943611198906d+00*HR3(0,-1,-1)
     $  + 1.3862943611198906d+00*HR3(0,-1,1)
     $  + 1.3862943611198906d+00*HR3(0,1,-1)
     $  + 1.3862943611198906d+00*HR3(0,1,1)
      Hi5(0,1,1,-1,-1) = 
     $  + 4.7533770109129867d-01*HR1(-1)
     $  - 1.2011325347955035d-01*HR1(-1)*HR1(-1)
     $  + 2.4022650695910071d-01*HR1(-1)*HR1(0)
     $  - 2.4022650695910071d-01*HR1(-1)*HR1(1)
     $  + 2.4022650695910071d-01*HR1(0)*HR1(1)
     $  + 4.7533770109129867d-01*HR1(1)
     $  + 2.4022650695910071d-01*HR2(-1,1)
     $  - 2.4022650695910071d-01*HR2(0,-1)
     $  - 2.4022650695910071d-01*HR2(0,1)
      Hi5(0,1,1,-1,1) = 
     $  + 1.4780047665430420d+00*HR1(-1)
     $  - 2.9112026323250625d-01*HR1(-1)*HR1(-1)
     $  + 4.1666666666666666d-02*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666d-01*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  - 5.0000000000000000d-01*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  + 5.8224052646501250d-01*HR1(-1)*HR1(0)
     $  - 5.8224052646501250d-01*HR1(-1)*HR1(1)
     $  - HR1( -1)*HR1(1)*HR2(0,-1)
     $  + HR1( -1)*HR3(-1,-1,1)
     $  + HR1( -1)*HR3(0,-1,-1)
     $  + HR1( -1)*HR3(0,0,-1)
     $  + 5.8224052646501250d-01*HR1(0)*HR1(1)
     $  + 1.4780047665430420d+00*HR1(1)
     $  + HR1(1) *HR3(0,-1,-1)
     $  + HR1(1) *HR3(0,0,-1)
     $  + 5.8224052646501250d-01*HR2(-1,1)
     $  + HR2( -1,1)*HR2(0,-1)
     $  - 5.8224052646501250d-01*HR2(0,-1)
     $  - 5.0000000000000000d-01*HR2(0,-1)*HR2(0,-1)
     $  - HR2(0, -1)*HR2(0,1)
     $  - 5.8224052646501250d-01*HR2(0,1)
     $  - HR4( -1,-1,-1,1)
     $  - HR4(0, -1,-1,1)
     $  + HR4(0, -1,0,1)
     $  + HR4(0,0, -1,1)
     $  + HR4(0,0,1, -1)
     $  + HR4(0,1, -1,-1)
      Hi5(0,1,1,1,-1) = 
     $  + 3.2877511713386177d-02*HR1(-1)
     $  + 5.3134677019160696d-01*HR1(-1)*HR1(-1)
     $  - 1.1552453009332421d-01*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(0)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(-1)*HR1(1)
     $  - 1.0626935403832139d+00*HR1(-1)*HR1(0)
     $  - 3.4657359027997265d-01*HR1(-1)*HR1(0)*HR1(0)
     $  + 6.9314718055994530d-01*HR1(-1)*HR1(0)*HR1(1)
     $  + 1.0626935403832139d+00*HR1(-1)*HR1(1)
     $  + 6.9314718055994530d-01*HR1(-1)*HR2(-1,1)
     $  - 3.4657359027997265d-01*HR1(0)*HR1(0)*HR1(1)
     $  - 1.0626935403832139d+00*HR1(0)*HR1(1)
     $  - 6.9314718055994530d-01*HR1(0)*HR2(-1,1)
     $  + 6.9314718055994530d-01*HR1(0)*HR2(0,-1)
     $  + 6.9314718055994530d-01*HR1(0)*HR2(0,1)
     $  + 3.2877511713386177d-02*HR1(1)
     $  - 1.0626935403832139d+00*HR2(-1,1)
     $  + 1.0626935403832139d+00*HR2(0,-1)
     $  + 1.0626935403832139d+00*HR2(0,1)
     $  - 6.9314718055994530d-01*HR3(-1,-1,1)
     $  - 6.9314718055994530d-01*HR3(0,-1,-1)
     $  - 6.9314718055994530d-01*HR3(0,0,-1)
     $  - 6.9314718055994530d-01*HR3(0,0,1)
     $  - 6.9314718055994530d-01*HR3(0,1,-1)
      endif  
      endif 
** nw > 4 endif 
      endif 
** (n1,n2) = (-1,1) -- completion endif 
      return 
      end 
************************************************************************ 
      subroutine fillirr1dhplatinf(x,nw,HX1,HX2,HX3,HX4,HX5, 
     $                                HY1,HY2,HY3,HY4,HY5, 
     $                                Hi1,Hi2,Hi3,Hi4,Hi5,n1,n2) 
** evaluates the HPL for y > r2p1
** fillirr1dhplatinf is called by eval1dhplatinf after calling 
** fillirr1dhplat0 with argument r=1/y 
** it is guaranteed that nw is in the range 2:4, and that (n1,n2) 
** take one of the pairs of values (0,1), (-1,0) or (-1,1) 
      implicit double precision (a-h,o-z) 
      dimension HX1(n1:n2),HX2(n1:n2,n1:n2),HX3(n1:n2,n1:n2,n1:n2), 
     $          HX4(n1:n2,n1:n2,n1:n2,n1:n2), 
     $          HX5(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension HY1(n1:n2),HY2(n1:n2,n1:n2),HY3(n1:n2,n1:n2,n1:n2), 
     $          HY4(n1:n2,n1:n2,n1:n2,n1:n2), 
     $          HY5(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2), 
     $          Hi5(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2) 
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
      if ( nw.gt.4 ) then 
      HY5(0,0,0,0,1) = 
     $  - 2.1646464674222763d+00*HX1(0)
     $  - 5.4831135561607547d-01*HX1(0)*HX1(0)*HX1(0)
     $  + 8.3333333333333333d-03*HX1(0)*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + HX5(0,0,0,0,1) 
      HY5(0,0,0,1,1) = 
     $  - 2.9176809454512223d+00
     $  - 2.1646464674222763d+00*HX1(0)
     $  + 6.0102845157979714d-01*HX1(0)*HX1(0)
     $  + 8.2246703342411321d-01*HX1(0)*HX1(0)*HX1(0)
     $  - 8.3333333333333333d-03*HX1(0)*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  - HX1(0) *HX4(0,0,0,1)
     $  + 3.0000000000000000d+00*HX5(0,0,0,0,1)
     $  - HX5(0,0,0,1,1) 
      HY5(0,0,1,0,1) = 
     $  + 3.7615063806157047d+00
     $  - 1.0823232337111381d+00*HX1(0)
     $  - 1.2020569031595942d+00*HX1(0)*HX1(0)
     $  + 5.4831135561607547d-01*HX1(0)*HX1(0)*HX1(0)
     $  - 8.3333333333333333d-03*HX1(0)*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0)*HX3(0,0,1)
     $  + 3.0000000000000000d+00*HX1(0)*HX4(0,0,0,1)
     $  + 3.2898681336964528d+00*HX3(0,0,1)
     $  - 7.0000000000000000d+00*HX5(0,0,0,0,1)
     $  - HX5(0,0,1,0,1) 
      HY5(0,0,1,1,1) = 
     $  + 3.0142321054406660d+00
     $  + 5.1410353601279064d+00*HX1(0)
     $  - 8.2246703342411321d-01*HX1(0)*HX1(0)*HX1(0)
     $  + 8.3333333333333333d-03*HX1(0)*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0)*HX3(0,0,1)
     $  - 2.0000000000000000d+00*HX1(0)*HX4(0,0,0,1)
     $  + HX1(0) *HX4(0,0,1,1)
     $  - 4.9348022005446793d+00*HX3(0,0,1)
     $  + 3.0000000000000000d+00*HX5(0,0,0,0,1)
     $  - 2.0000000000000000d+00*HX5(0,0,0,1,1)
     $  + HX5(0,0,1,1,1) 
      HY5(0,1,0,1,1) = 
     $  - 1.5553916327150548d+00
     $  + 8.1174242528335364d-01*HX1(0)
     $  - 6.0102845157979714d-01*HX1(0)*HX1(0)
     $  - 8.2246703342411321d-01*HX1(0)*HX1(0)*HX1(0)
     $  + 8.3333333333333333d-03*HX1(0)*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0)*HX2(0,1)
     $  - HX1(0) *HX1(0)*HX3(0,0,1)
     $  - 4.9348022005446793d+00*HX1(0)*HX2(0,1)
     $  + 5.0000000000000000d-01*HX1(0)*HX2(0,1)*HX2(0,1)
     $  + 4.0000000000000000d+00*HX1(0)*HX4(0,0,0,1)
     $  - 2.0000000000000000d+00*HX1(0)*HX4(0,0,1,1)
     $  - 1.2020569031595942d+00*HX2(0,1)
     $  - HX2(0,1) *HX3(0,0,1)
     $  + 9.8696044010893586d+00*HX3(0,0,1)
     $  - 7.0000000000000000d+00*HX5(0,0,0,0,1)
     $  + 7.0000000000000000d+00*HX5(0,0,0,1,1)
     $  + HX5(0,0,1,0,1) 
     $  + HX5(0,1,0,1,1) 
      HY5(0,1,1,1,1) = 
     $  + 1.0369277551433699d+00
     $  - 4.0587121264167682d+00*HX1(0)
     $  + 8.2246703342411321d-01*HX1(0)*HX1(0)*HX1(0)
     $  - 8.3333333333333333d-03*HX1(0)*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  - 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0)*HX2(0,1)
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0)*HX3(0,0,1)
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0)*HX3(0,1,1)
     $  + 4.9348022005446793d+00*HX1(0)*HX2(0,1)
     $  - HX1(0) *HX4(0,0,0,1)
     $  + HX1(0) *HX4(0,0,1,1)
     $  - HX1(0) *HX4(0,1,1,1)
     $  - 4.9348022005446793d+00*HX3(0,0,1)
     $  + 4.9348022005446793d+00*HX3(0,1,1)
     $  + HX5(0,0,0,0,1) 
     $  - HX5(0,0,0,1,1) 
     $  + HX5(0,0,1,1,1) 
     $  - HX5(0,1,1,1,1) 
      Hi5(0,0,0,0,1) = 
     $  + 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0)
      Hi5(0,0,0,1,1) = 
     $  + 1.0823232337111381d+00
     $  + 1.2020569031595942d+00*HX1(0)
     $  + 8.2246703342411321d-01*HX1(0)*HX1(0)
     $  - 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  - HX4(0,0,0,1) 
      Hi5(0,0,1,0,1) = 
     $  - 3.2469697011334145d+00
     $  - 2.4041138063191885d+00*HX1(0)
     $  - 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  - HX1(0) *HX3(0,0,1)
     $  + 3.0000000000000000d+00*HX4(0,0,0,1)
      Hi5(0,0,1,1,1) = 
     $  + 1.8940656589944918d+00
     $  - 8.2246703342411321d-01*HX1(0)*HX1(0)
     $  + 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + HX1(0) *HX3(0,0,1)
     $  - 2.0000000000000000d+00*HX4(0,0,0,1)
     $  + HX4(0,0,1,1) 
      Hi5(0,1,0,1,1) = 
     $  - 2.4352272758500609d+00
     $  - 1.2020569031595942d+00*HX1(0)
     $  - 8.2246703342411321d-01*HX1(0)*HX1(0)
     $  + 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0)*HX2(0,1)
     $  - 2.0000000000000000d+00*HX1(0)*HX3(0,0,1)
     $  - 1.6449340668482264d+00*HX2(0,1)
     $  + 5.0000000000000000d-01*HX2(0,1)*HX2(0,1)
     $  + 4.0000000000000000d+00*HX4(0,0,0,1)
     $  - 2.0000000000000000d+00*HX4(0,0,1,1)
      Hi5(0,1,1,1,1) = 
     $  - 8.1174242528335364d-01
     $  + 8.2246703342411321d-01*HX1(0)*HX1(0)
     $  - 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0)*HX2(0,1)
     $  + HX1(0) *HX3(0,0,1)
     $  - HX1(0) *HX3(0,1,1)
     $  + 1.6449340668482264d+00*HX2(0,1)
     $  - HX4(0,0,0,1) 
     $  + HX4(0,0,1,1) 
     $  - HX4(0,1,1,1) 
      endif 
** nw > 4 endif 
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
      if ( nw.gt.4 ) then
      HY5(0,0,0,0,-1) = 
     $  - 1.8940656589944918d+00*HX1(0)
     $  - 2.7415567780803773d-01*HX1(0)*HX1(0)*HX1(0)
     $  - 8.3333333333333333d-03*HX1(0)*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + HX5(0,0,0,0, -1)
      HY5(0,0,0,-1,-1) = 
     $  + 3.0142321054406660d+00
     $  + 1.8940656589944918d+00*HX1(0)
     $  + 6.0102845157979714d-01*HX1(0)*HX1(0)
     $  - 8.3333333333333333d-03*HX1(0)*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + HX1(0) *HX4(0,0,0,-1)
     $  - HX5(0,0,0, -1,-1)
     $  - 3.0000000000000000d+00*HX5(0,0,0,0,-1)
      HY5(0,0,-1,0,-1) = 
     $  - 8.1023197211680719d+00
     $  - 5.1410353601279064d+00*HX1(0)
     $  - 1.2020569031595942d+00*HX1(0)*HX1(0)
     $  - 2.7415567780803773d-01*HX1(0)*HX1(0)*HX1(0)
     $  - 8.3333333333333333d-03*HX1(0)*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0)*HX3(0,0,-1)
     $  - 3.0000000000000000d+00*HX1(0)*HX4(0,0,0,-1)
     $  + 1.6449340668482264d+00*HX3(0,0,-1)
     $  - HX5(0,0, -1,0,-1)
     $  + 7.0000000000000000d+00*HX5(0,0,0,0,-1)
      HY5(0,0,-1,-1,-1) = 
     $  - 3.0142321054406660d+00
     $  - 1.0823232337111381d+00*HX1(0)
     $  - 8.3333333333333333d-03*HX1(0)*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0)*HX3(0,0,-1)
     $  - HX1(0) *HX4(0,0,-1,-1)
     $  - 2.0000000000000000d+00*HX1(0)*HX4(0,0,0,-1)
     $  + HX5(0,0, -1,-1,-1)
     $  + 2.0000000000000000d+00*HX5(0,0,0,-1,-1)
     $  + 3.0000000000000000d+00*HX5(0,0,0,0,-1)
      HY5(0,-1,0,-1,-1) = 
     $  + 7.4873046836069432d+00
     $  + 3.2469697011334145d+00*HX1(0)
     $  + 6.0102845157979714d-01*HX1(0)*HX1(0)
     $  - 8.3333333333333333d-03*HX1(0)*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0)*HX2(0,-1)
     $  - HX1(0) *HX1(0)*HX3(0,0,-1)
     $  - 5.0000000000000000d-01*HX1(0)*HX2(0,-1)*HX2(0,-1)
     $  + 2.0000000000000000d+00*HX1(0)*HX4(0,0,-1,-1)
     $  + 4.0000000000000000d+00*HX1(0)*HX4(0,0,0,-1)
     $  - 1.2020569031595942d+00*HX2(0,-1)
     $  + HX2(0, -1)*HX3(0,0,-1)
     $  + HX5(0, -1,0,-1,-1)
     $  - HX5(0,0, -1,0,-1)
     $  - 7.0000000000000000d+00*HX5(0,0,0,-1,-1)
     $  - 7.0000000000000000d+00*HX5(0,0,0,0,-1)
      HY5(0,-1,-1,-1,-1) = 
     $  + 1.0369277551433699d+00
     $  - 8.3333333333333333d-03*HX1(0)*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0)*HX2(0,-1)
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0)*HX3(0,-1,-1)
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0)*HX3(0,0,-1)
     $  + HX1(0) *HX4(0,-1,-1,-1)
     $  + HX1(0) *HX4(0,0,-1,-1)
     $  + HX1(0) *HX4(0,0,0,-1)
     $  - HX5(0, -1,-1,-1,-1)
     $  - HX5(0,0, -1,-1,-1)
     $  - HX5(0,0,0, -1,-1)
     $  - HX5(0,0,0,0, -1)
      endif 
** nw > 4 endif 
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
      if ( nw.gt.4 ) then 
      HY5(-1,-1,-1,-1,1) = 
     $  - 2.1900870176160439d+00
     $  + 2.4278628067547031d+00*HX1(-1)
     $  - 1.3810359531144620d+00*HX1(-1)*HX1(-1)
     $  + 4.1123351671205660d-01*HX1(-1)*HX1(-1)*HX1(-1)
     $  + 4.1666666666666666d-02*HX1(-1)*HX1(-1)*HX1(-1)*HX1(-1)*HX1(0)
     $  - 8.3333333333333333d-02*HX1(-1)*HX1(-1)*HX1(-1)*HX1(0)*HX1(0)
     $  - 1.6666666666666666d-01*HX1(-1)*HX1(-1)*HX1(-1)*HX2(0,-1)
     $  - 1.6666666666666666d-01*HX1(-1)*HX1(-1)*HX1(-1)*HX2(0,1)
     $  - 1.2337005501361698d+00*HX1(-1)*HX1(-1)*HX1(0)
     $  + 8.3333333333333333d-02*HX1(-1)*HX1(-1)*HX1(0)*HX1(0)*HX1(0)
     $  + 5.0000000000000000d-01*HX1(-1)*HX1(-1)*HX3(0,-1,-1)
     $  + 5.0000000000000000d-01*HX1(-1)*HX1(-1)*HX3(0,0,-1)
     $  + 5.0000000000000000d-01*HX1(-1)*HX1(-1)*HX3(0,0,1)
     $  + 5.0000000000000000d-01*HX1(-1)*HX1(-1)*HX3(0,1,-1)
     $  + 2.7620719062289241d+00*HX1(-1)*HX1(0)
     $  + 1.2337005501361698d+00*HX1(-1)*HX1(0)*HX1(0)
     $  - 4.1666666666666666d-02*HX1(-1)*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  - HX1( -1)*HX4(0,-1,-1,-1)
     $  - HX1( -1)*HX4(0,0,-1,-1)
     $  - HX1( -1)*HX4(0,0,0,-1)
     $  - HX1( -1)*HX4(0,0,0,1)
     $  - HX1( -1)*HX4(0,0,1,-1)
     $  - HX1( -1)*HX4(0,1,-1,-1)
     $  - 2.4278628067547031d+00*HX1(0)
     $  - 1.3810359531144620d+00*HX1(0)*HX1(0)
     $  - 4.1123351671205660d-01*HX1(0)*HX1(0)*HX1(0)
     $  + 8.3333333333333333d-03*HX1(0)*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + HX5( -1,-1,-1,-1,1)
     $  + HX5(0, -1,-1,-1,-1)
     $  + HX5(0,0, -1,-1,-1)
     $  + HX5(0,0,0, -1,-1)
     $  + HX5(0,0,0,0, -1)
     $  + HX5(0,0,0,0,1) 
     $  + HX5(0,0,0,1, -1)
     $  + HX5(0,0,1, -1,-1)
     $  + HX5(0,1, -1,-1,-1)
      HY5(-1,-1,-1,1,1) = 
     $  - 5.5581622138319701d+00
     $  + 2.0293560632083841d+00*HX1(-1)
     $  + 1.3810359531144620d+00*HX1(-1)*HX1(-1)
     $  - 8.2246703342411321d-01*HX1(-1)*HX1(-1)*HX1(-1)
     $  + 8.3333333333333333d-02*HX1(-1)*HX1(-1)*HX1(-1)*HX1(0)*HX1(0)
     $  + 2.4674011002723396d+00*HX1(-1)*HX1(-1)*HX1(0)
     $  - 8.3333333333333333d-02*HX1(-1)*HX1(-1)*HX1(0)*HX1(0)*HX1(0)
     $  - 5.0000000000000000d-01*HX1(-1)*HX1(-1)*HX1(0)*HX2(0,-1)
     $  - 5.0000000000000000d-01*HX1(-1)*HX1(-1)*HX1(0)*HX2(0,1)
     $  - 5.0000000000000000d-01*HX1(-1)*HX1(-1)*HX3(0,-1,1)
     $  + 5.0000000000000000d-01*HX1(-1)*HX1(-1)*HX3(0,0,-1)
     $  + 5.0000000000000000d-01*HX1(-1)*HX1(-1)*HX3(0,0,1)
     $  - 5.0000000000000000d-01*HX1(-1)*HX1(-1)*HX3(0,1,1)
     $  - 2.7620719062289241d+00*HX1(-1)*HX1(0)
     $  - 2.4674011002723396d+00*HX1(-1)*HX1(0)*HX1(0)
     $  + 4.1666666666666666d-02*HX1(-1)*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + HX1( -1)*HX1(0)*HX3(0,-1,-1)
     $  + HX1( -1)*HX1(0)*HX3(0,0,-1)
     $  + HX1( -1)*HX1(0)*HX3(0,0,1)
     $  + HX1( -1)*HX1(0)*HX3(0,1,-1)
     $  + HX1( -1)*HX4(0,-1,-1,1)
     $  + HX1( -1)*HX4(0,-1,1,-1)
     $  - HX1( -1)*HX4(0,0,-1,-1)
     $  + HX1( -1)*HX4(0,0,-1,1)
     $  - 2.0000000000000000d+00*HX1(-1)*HX4(0,0,0,-1)
     $  - 2.0000000000000000d+00*HX1(-1)*HX4(0,0,0,1)
     $  - HX1( -1)*HX4(0,0,1,-1)
     $  + HX1( -1)*HX4(0,0,1,1)
     $  + HX1( -1)*HX4(0,1,-1,1)
     $  + HX1( -1)*HX4(0,1,1,-1)
     $  - 2.0293560632083841d+00*HX1(0)
     $  + 1.3810359531144620d+00*HX1(0)*HX1(0)
     $  + 8.2246703342411321d-01*HX1(0)*HX1(0)*HX1(0)
     $  - 8.3333333333333333d-03*HX1(0)*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + HX1(0) *HX4(-1,-1,-1,1)
     $  - HX1(0) *HX4(0,-1,-1,-1)
     $  - HX1(0) *HX4(0,0,-1,-1)
     $  - HX1(0) *HX4(0,0,0,-1)
     $  - HX1(0) *HX4(0,0,0,1)
     $  - HX1(0) *HX4(0,0,1,-1)
     $  - HX1(0) *HX4(0,1,-1,-1)
     $  + HX5( -1,-1,-1,1,1)
     $  - HX5(0, -1,-1,-1,1)
     $  - HX5(0, -1,-1,1,-1)
     $  - HX5(0, -1,1,-1,-1)
     $  + HX5(0,0, -1,-1,-1)
     $  - HX5(0,0, -1,-1,1)
     $  - HX5(0,0, -1,1,-1)
     $  + 2.0000000000000000d+00*HX5(0,0,0,-1,-1)
     $  - HX5(0,0,0, -1,1)
     $  + 3.0000000000000000d+00*HX5(0,0,0,0,-1)
     $  + 3.0000000000000000d+00*HX5(0,0,0,0,1)
     $  + 2.0000000000000000d+00*HX5(0,0,0,1,-1)
     $  - HX5(0,0,0,1,1) 
     $  + HX5(0,0,1, -1,-1)
     $  - HX5(0,0,1, -1,1)
     $  - HX5(0,0,1,1, -1)
     $  - HX5(0,1, -1,-1,1)
     $  - HX5(0,1, -1,1,-1)
     $  - HX5(0,1,1, -1,-1)
      HY5(-1,-1,1,-1,1) = 
     $  + 9.2329419831013177d+00
     $  - 3.3856186219460558d+00*HX1(-1)
     $  + 6.5847232569963413d-01*HX1(-1)*HX1(-1)
     $  + 1.2337005501361698d+00*HX1(-1)*HX1(-1)*HX1(0)
     $  - 8.3333333333333333d-02*HX1(-1)*HX1(-1)*HX1(0)*HX1(0)*HX1(0)
     $  + 5.0000000000000000d-01*HX1(-1)*HX1(-1)*HX1(0)*HX2(0,-1)
     $  + 5.0000000000000000d-01*HX1(-1)*HX1(-1)*HX1(0)*HX2(0,1)
     $  + HX1( -1)*HX1(-1)*HX3(0,-1,1)
     $  - HX1( -1)*HX1(-1)*HX3(0,0,-1)
     $  - HX1( -1)*HX1(-1)*HX3(0,0,1)
     $  + HX1( -1)*HX1(-1)*HX3(0,1,1)
     $  - 1.3169446513992682d+00*HX1(-1)*HX1(0)
     $  - 1.2337005501361698d+00*HX1(-1)*HX1(0)*HX1(0)
     $  + 4.1666666666666666d-02*HX1(-1)*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + 5.0000000000000000d-01*HX1(-1)*HX1(0)*HX1(0)*HX2(0,-1)
     $  + 5.0000000000000000d-01*HX1(-1)*HX1(0)*HX1(0)*HX2(0,1)
     $  + HX1( -1)*HX1(0)*HX3(-1,-1,1)
     $  - 2.0000000000000000d+00*HX1(-1)*HX1(0)*HX3(0,-1,-1)
     $  - 2.0000000000000000d+00*HX1(-1)*HX1(0)*HX3(0,0,-1)
     $  - 2.0000000000000000d+00*HX1(-1)*HX1(0)*HX3(0,0,1)
     $  - 2.0000000000000000d+00*HX1(-1)*HX1(0)*HX3(0,1,-1)
     $  - 2.4674011002723396d+00*HX1(-1)*HX2(0,-1)
     $  + 5.0000000000000000d-01*HX1(-1)*HX2(0,-1)*HX2(0,-1)
     $  + HX1( -1)*HX2(0,-1)*HX2(0,1)
     $  - 2.4674011002723396d+00*HX1(-1)*HX2(0,1)
     $  + 5.0000000000000000d-01*HX1(-1)*HX2(0,1)*HX2(0,1)
     $  - 2.0000000000000000d+00*HX1(-1)*HX4(0,-1,-1,1)
     $  - 2.0000000000000000d+00*HX1(-1)*HX4(0,-1,1,-1)
     $  + 2.0000000000000000d+00*HX1(-1)*HX4(0,0,-1,-1)
     $  - 2.0000000000000000d+00*HX1(-1)*HX4(0,0,-1,1)
     $  + 4.0000000000000000d+00*HX1(-1)*HX4(0,0,0,-1)
     $  + 4.0000000000000000d+00*HX1(-1)*HX4(0,0,0,1)
     $  + 2.0000000000000000d+00*HX1(-1)*HX4(0,0,1,-1)
     $  - 2.0000000000000000d+00*HX1(-1)*HX4(0,0,1,1)
     $  - 2.0000000000000000d+00*HX1(-1)*HX4(0,1,-1,1)
     $  - 2.0000000000000000d+00*HX1(-1)*HX4(0,1,1,-1)
     $  + 3.3856186219460558d+00*HX1(0)
     $  + 6.5847232569963413d-01*HX1(0)*HX1(0)
     $  + 4.1123351671205660d-01*HX1(0)*HX1(0)*HX1(0)
     $  - 8.3333333333333333d-03*HX1(0)*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0)*HX3(-1,-1,1)
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0)*HX3(0,-1,-1)
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0)*HX3(0,0,-1)
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0)*HX3(0,0,1)
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0)*HX3(0,1,-1)
     $  - 3.0000000000000000d+00*HX1(0)*HX4(-1,-1,-1,1)
     $  + 3.0000000000000000d+00*HX1(0)*HX4(0,-1,-1,-1)
     $  + 3.0000000000000000d+00*HX1(0)*HX4(0,0,-1,-1)
     $  + 3.0000000000000000d+00*HX1(0)*HX4(0,0,0,-1)
     $  + 3.0000000000000000d+00*HX1(0)*HX4(0,0,0,1)
     $  + 3.0000000000000000d+00*HX1(0)*HX4(0,0,1,-1)
     $  + 3.0000000000000000d+00*HX1(0)*HX4(0,1,-1,-1)
     $  - HX2(0, -1)*HX3(-1,-1,1)
     $  - HX2(0, -1)*HX3(0,-1,-1)
     $  - HX2(0, -1)*HX3(0,1,-1)
     $  - HX2(0,1) *HX3(-1,-1,1)
     $  - HX2(0,1) *HX3(0,1,-1)
     $  + 2.4674011002723396d+00*HX3(-1,-1,1)
     $  + 2.4674011002723396d+00*HX3(0,-1,-1)
     $  + 2.4674011002723396d+00*HX3(0,0,-1)
     $  + 2.4674011002723396d+00*HX3(0,0,1)
     $  + 2.4674011002723396d+00*HX3(0,1,-1)
     $  + HX5( -1,-1,1,-1,1)
     $  + 3.0000000000000000d+00*HX5(0,-1,-1,-1,1)
     $  - HX5(0, -1,-1,0,1)
     $  + 2.0000000000000000d+00*HX5(0,-1,-1,1,-1)
     $  + HX5(0, -1,0,-1,-1)
     $  - HX5(0, -1,0,-1,1)
     $  + 2.0000000000000000d+00*HX5(0,-1,1,-1,-1)
     $  + HX5(0,0, -1,-1,1)
     $  - HX5(0,0, -1,0,-1)
     $  - HX5(0,0, -1,0,1)
     $  + 2.0000000000000000d+00*HX5(0,0,-1,1,-1)
     $  - 7.0000000000000000d+00*HX5(0,0,0,-1,-1)
     $  - 7.0000000000000000d+00*HX5(0,0,0,0,-1)
     $  - 7.0000000000000000d+00*HX5(0,0,0,0,1)
     $  - 7.0000000000000000d+00*HX5(0,0,0,1,-1)
     $  - 2.0000000000000000d+00*HX5(0,0,1,-1,-1)
     $  + 3.0000000000000000d+00*HX5(0,0,1,-1,1)
     $  - HX5(0,0,1,0, -1)
     $  - HX5(0,0,1,0,1) 
     $  + 4.0000000000000000d+00*HX5(0,0,1,1,-1)
     $  + 3.0000000000000000d+00*HX5(0,1,-1,-1,1)
     $  + 2.0000000000000000d+00*HX5(0,1,-1,1,-1)
     $  + HX5(0,1,0,1, -1)
     $  + 2.0000000000000000d+00*HX5(0,1,1,-1,-1)
      HY5(-1,-1,1,1,1) = 
     $  + 5.5581622138319701d+00
     $  - 6.4865749331714713d+00*HX1(-1)
     $  - 2.4674011002723396d+00*HX1(-1)*HX1(-1)*HX1(0)
     $  + 8.3333333333333333d-02*HX1(-1)*HX1(-1)*HX1(0)*HX1(0)*HX1(0)
     $  + 2.4674011002723396d+00*HX1(-1)*HX1(0)*HX1(0)
     $  - 4.1666666666666666d-02*HX1(-1)*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  - 5.0000000000000000d-01*HX1(-1)*HX1(0)*HX1(0)*HX2(0,-1)
     $  - 5.0000000000000000d-01*HX1(-1)*HX1(0)*HX1(0)*HX2(0,1)
     $  - HX1( -1)*HX1(0)*HX3(0,-1,1)
     $  + HX1( -1)*HX1(0)*HX3(0,0,-1)
     $  + HX1( -1)*HX1(0)*HX3(0,0,1)
     $  - HX1( -1)*HX1(0)*HX3(0,1,1)
     $  + 4.9348022005446793d+00*HX1(-1)*HX2(0,-1)
     $  + 4.9348022005446793d+00*HX1(-1)*HX2(0,1)
     $  - HX1( -1)*HX4(0,-1,1,1)
     $  + HX1( -1)*HX4(0,0,-1,1)
     $  - HX1( -1)*HX4(0,0,0,-1)
     $  - HX1( -1)*HX4(0,0,0,1)
     $  + HX1( -1)*HX4(0,0,1,1)
     $  - HX1( -1)*HX4(0,1,1,1)
     $  + 6.4865749331714713d+00*HX1(0)
     $  - 8.2246703342411321d-01*HX1(0)*HX1(0)*HX1(0)
     $  + 8.3333333333333333d-03*HX1(0)*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0)*HX3(-1,-1,1)
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0)*HX3(0,-1,-1)
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0)*HX3(0,0,-1)
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0)*HX3(0,0,1)
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0)*HX3(0,1,-1)
     $  + HX1(0) *HX4(-1,-1,1,1)
     $  + HX1(0) *HX4(0,-1,-1,1)
     $  + HX1(0) *HX4(0,-1,1,-1)
     $  - HX1(0) *HX4(0,0,-1,-1)
     $  + HX1(0) *HX4(0,0,-1,1)
     $  - 2.0000000000000000d+00*HX1(0)*HX4(0,0,0,-1)
     $  - 2.0000000000000000d+00*HX1(0)*HX4(0,0,0,1)
     $  - HX1(0) *HX4(0,0,1,-1)
     $  + HX1(0) *HX4(0,0,1,1)
     $  + HX1(0) *HX4(0,1,-1,1)
     $  + HX1(0) *HX4(0,1,1,-1)
     $  - 4.9348022005446793d+00*HX3(-1,-1,1)
     $  - 4.9348022005446793d+00*HX3(0,-1,-1)
     $  - 4.9348022005446793d+00*HX3(0,0,-1)
     $  - 4.9348022005446793d+00*HX3(0,0,1)
     $  - 4.9348022005446793d+00*HX3(0,1,-1)
     $  + HX5( -1,-1,1,1,1)
     $  + HX5(0, -1,-1,1,1)
     $  + HX5(0, -1,1,-1,1)
     $  + HX5(0, -1,1,1,-1)
     $  - HX5(0,0, -1,-1,1)
     $  - HX5(0,0, -1,1,-1)
     $  + HX5(0,0, -1,1,1)
     $  + HX5(0,0,0, -1,-1)
     $  - 2.0000000000000000d+00*HX5(0,0,0,-1,1)
     $  + 3.0000000000000000d+00*HX5(0,0,0,0,-1)
     $  + 3.0000000000000000d+00*HX5(0,0,0,0,1)
     $  + HX5(0,0,0,1, -1)
     $  - 2.0000000000000000d+00*HX5(0,0,0,1,1)
     $  - HX5(0,0,1, -1,1)
     $  - HX5(0,0,1,1, -1)
     $  + HX5(0,0,1,1,1) 
     $  + HX5(0,1, -1,1,1)
     $  + HX5(0,1,1, -1,1)
     $  + HX5(0,1,1,1, -1)
      HY5(-1,1,-1,1,1) = 
     $  - 7.7439193717015657d-01
     $  + 8.5393570350547734d-01*HX1(-1)
     $  + 2.7620719062289241d+00*HX1(-1)*HX1(0)
     $  + 2.4674011002723396d+00*HX1(-1)*HX1(0)*HX1(0)
     $  - 4.1666666666666666d-02*HX1(-1)*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + 5.0000000000000000d-01*HX1(-1)*HX1(0)*HX1(0)*HX2(-1,1)
     $  + 5.0000000000000000d-01*HX1(-1)*HX1(0)*HX1(0)*HX2(0,-1)
     $  + 5.0000000000000000d-01*HX1(-1)*HX1(0)*HX1(0)*HX2(0,1)
     $  + 2.0000000000000000d+00*HX1(-1)*HX1(0)*HX3(0,-1,1)
     $  - 2.0000000000000000d+00*HX1(-1)*HX1(0)*HX3(0,0,-1)
     $  - 2.0000000000000000d+00*HX1(-1)*HX1(0)*HX3(0,0,1)
     $  + 2.0000000000000000d+00*HX1(-1)*HX1(0)*HX3(0,1,1)
     $  - 4.9348022005446793d+00*HX1(-1)*HX2(-1,1)
     $  - 4.9348022005446793d+00*HX1(-1)*HX2(0,-1)
     $  - 4.9348022005446793d+00*HX1(-1)*HX2(0,1)
     $  + 3.0000000000000000d+00*HX1(-1)*HX4(0,-1,1,1)
     $  - 3.0000000000000000d+00*HX1(-1)*HX4(0,0,-1,1)
     $  + 3.0000000000000000d+00*HX1(-1)*HX4(0,0,0,-1)
     $  + 3.0000000000000000d+00*HX1(-1)*HX4(0,0,0,1)
     $  - 3.0000000000000000d+00*HX1(-1)*HX4(0,0,1,1)
     $  + 3.0000000000000000d+00*HX1(-1)*HX4(0,1,1,1)
     $  - 8.5393570350547734d-01*HX1(0)
     $  - 1.3810359531144620d+00*HX1(0)*HX1(0)
     $  - 8.2246703342411321d-01*HX1(0)*HX1(0)*HX1(0)
     $  + 8.3333333333333333d-03*HX1(0)*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  - 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0)*HX2(-1,1)
     $  + 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0)*HX2(0,-1)
     $  + 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0)*HX2(0,1)
     $  - HX1(0) *HX1(0)*HX3(-1,-1,1)
     $  - HX1(0) *HX1(0)*HX3(0,-1,-1)
     $  - HX1(0) *HX1(0)*HX3(0,0,-1)
     $  - HX1(0) *HX1(0)*HX3(0,0,1)
     $  - HX1(0) *HX1(0)*HX3(0,1,-1)
     $  + 4.9348022005446793d+00*HX1(0)*HX2(-1,1)
     $  + 5.0000000000000000d-01*HX1(0)*HX2(-1,1)*HX2(-1,1)
     $  - HX1(0) *HX2(-1,1)*HX2(0,-1)
     $  - HX1(0) *HX2(-1,1)*HX2(0,1)
     $  - 4.9348022005446793d+00*HX1(0)*HX2(0,-1)
     $  + 5.0000000000000000d-01*HX1(0)*HX2(0,-1)*HX2(0,-1)
     $  + HX1(0) *HX2(0,-1)*HX2(0,1)
     $  - 4.9348022005446793d+00*HX1(0)*HX2(0,1)
     $  + 5.0000000000000000d-01*HX1(0)*HX2(0,1)*HX2(0,1)
     $  - 2.0000000000000000d+00*HX1(0)*HX4(-1,-1,1,1)
     $  - 2.0000000000000000d+00*HX1(0)*HX4(0,-1,-1,1)
     $  - 2.0000000000000000d+00*HX1(0)*HX4(0,-1,1,-1)
     $  + 2.0000000000000000d+00*HX1(0)*HX4(0,0,-1,-1)
     $  - 2.0000000000000000d+00*HX1(0)*HX4(0,0,-1,1)
     $  + 4.0000000000000000d+00*HX1(0)*HX4(0,0,0,-1)
     $  + 4.0000000000000000d+00*HX1(0)*HX4(0,0,0,1)
     $  + 2.0000000000000000d+00*HX1(0)*HX4(0,0,1,-1)
     $  - 2.0000000000000000d+00*HX1(0)*HX4(0,0,1,1)
     $  - 2.0000000000000000d+00*HX1(0)*HX4(0,1,-1,1)
     $  - 2.0000000000000000d+00*HX1(0)*HX4(0,1,1,-1)
     $  + 2.7620719062289241d+00*HX2(-1,1)
     $  - HX2( -1,1)*HX3(0,-1,1)
     $  + HX2( -1,1)*HX3(0,0,-1)
     $  + HX2( -1,1)*HX3(0,0,1)
     $  - HX2( -1,1)*HX3(0,1,1)
     $  - 2.7620719062289241d+00*HX2(0,-1)
     $  - HX2(0, -1)*HX3(0,0,-1)
     $  - HX2(0, -1)*HX3(0,0,1)
     $  - 2.7620719062289241d+00*HX2(0,1)
     $  + HX2(0,1) *HX3(0,-1,1)
     $  - HX2(0,1) *HX3(0,0,-1)
     $  - HX2(0,1) *HX3(0,0,1)
     $  + 9.8696044010893586d+00*HX3(-1,-1,1)
     $  + 9.8696044010893586d+00*HX3(0,-1,-1)
     $  + 9.8696044010893586d+00*HX3(0,0,-1)
     $  + 9.8696044010893586d+00*HX3(0,0,1)
     $  + 9.8696044010893586d+00*HX3(0,1,-1)
     $  + HX5( -1,1,-1,1,1)
     $  - 2.0000000000000000d+00*HX5(0,-1,-1,1,1)
     $  + HX5(0, -1,0,-1,1)
     $  - HX5(0, -1,0,1,1)
     $  - 2.0000000000000000d+00*HX5(0,-1,1,-1,1)
     $  - HX5(0, -1,1,0,1)
     $  - 3.0000000000000000d+00*HX5(0,-1,1,1,-1)
     $  + 4.0000000000000000d+00*HX5(0,0,-1,-1,1)
     $  + HX5(0,0, -1,0,-1)
     $  + HX5(0,0, -1,0,1)
     $  + 3.0000000000000000d+00*HX5(0,0,-1,1,-1)
     $  - 4.0000000000000000d+00*HX5(0,0,-1,1,1)
     $  + 7.0000000000000000d+00*HX5(0,0,0,-1,1)
     $  - 7.0000000000000000d+00*HX5(0,0,0,0,-1)
     $  - 7.0000000000000000d+00*HX5(0,0,0,0,1)
     $  + 7.0000000000000000d+00*HX5(0,0,0,1,1)
     $  + 2.0000000000000000d+00*HX5(0,0,1,-1,1)
     $  + HX5(0,0,1,0, -1)
     $  + HX5(0,0,1,0,1) 
     $  + 3.0000000000000000d+00*HX5(0,0,1,1,-1)
     $  - 2.0000000000000000d+00*HX5(0,1,-1,1,1)
     $  + HX5(0,1,0,1,1) 
     $  - 2.0000000000000000d+00*HX5(0,1,1,-1,1)
     $  - 3.0000000000000000d+00*HX5(0,1,1,1,-1)
      HY5(-1,1,1,1,1) = 
     $  + 2.1900870176160439d+00
     $  + 4.0587121264167682d+00*HX1(-1)
     $  - 2.4674011002723396d+00*HX1(-1)*HX1(0)*HX1(0)
     $  + 4.1666666666666666d-02*HX1(-1)*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  - 4.0587121264167682d+00*HX1(0)
     $  + 8.2246703342411321d-01*HX1(0)*HX1(0)*HX1(0)
     $  - 8.3333333333333333d-03*HX1(0)*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0)*HX2(-1,1)
     $  - 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0)*HX2(0,-1)
     $  - 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0)*HX2(0,1)
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0)*HX3(-1,1,1)
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0)*HX3(0,-1,1)
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0)*HX3(0,0,-1)
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0)*HX3(0,0,1)
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0)*HX3(0,1,1)
     $  - 4.9348022005446793d+00*HX1(0)*HX2(-1,1)
     $  + 4.9348022005446793d+00*HX1(0)*HX2(0,-1)
     $  + 4.9348022005446793d+00*HX1(0)*HX2(0,1)
     $  + HX1(0) *HX4(-1,1,1,1)
     $  - HX1(0) *HX4(0,-1,1,1)
     $  + HX1(0) *HX4(0,0,-1,1)
     $  - HX1(0) *HX4(0,0,0,-1)
     $  - HX1(0) *HX4(0,0,0,1)
     $  + HX1(0) *HX4(0,0,1,1)
     $  - HX1(0) *HX4(0,1,1,1)
     $  - 4.9348022005446793d+00*HX3(-1,1,1)
     $  + 4.9348022005446793d+00*HX3(0,-1,1)
     $  - 4.9348022005446793d+00*HX3(0,0,-1)
     $  - 4.9348022005446793d+00*HX3(0,0,1)
     $  + 4.9348022005446793d+00*HX3(0,1,1)
     $  + HX5( -1,1,1,1,1)
     $  - HX5(0, -1,1,1,1)
     $  + HX5(0,0, -1,1,1)
     $  - HX5(0,0,0, -1,1)
     $  + HX5(0,0,0,0, -1)
     $  + HX5(0,0,0,0,1) 
     $  - HX5(0,0,0,1,1) 
     $  + HX5(0,0,1,1,1) 
     $  - HX5(0,1,1,1,1) 
      HY5(0,-1,-1,-1,1) = 
     $  - 2.1501021594785064d+00
     $  - 2.4278628067547031d+00*HX1(0)
     $  - 1.3810359531144620d+00*HX1(0)*HX1(0)
     $  - 4.1123351671205660d-01*HX1(0)*HX1(0)*HX1(0)
     $  + 8.3333333333333333d-03*HX1(0)*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  - 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0)*HX2(0,-1)
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0)*HX3(0,-1,-1)
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0)*HX3(0,0,-1)
     $  + 2.4674011002723396d+00*HX1(0)*HX2(0,-1)
     $  - HX1(0) *HX4(0,-1,-1,-1)
     $  - HX1(0) *HX4(0,0,-1,-1)
     $  - HX1(0) *HX4(0,0,0,-1)
     $  + 2.7620719062289241d+00*HX2(0,-1)
     $  + HX2(0, -1)*HX3(0,-1,-1)
     $  - HX2(0, -1)*HX3(0,0,-1)
     $  - HX2(0, -1)*HX3(0,0,1)
     $  - 2.4674011002723396d+00*HX3(0,-1,-1)
     $  - 2.4674011002723396d+00*HX3(0,0,-1)
     $  - HX5(0, -1,-1,-1,1)
     $  + HX5(0, -1,-1,0,1)
     $  - 2.0000000000000000d+00*HX5(0,-1,0,-1,-1)
     $  + HX5(0, -1,0,-1,1)
     $  - 4.0000000000000000d+00*HX5(0,0,-1,-1,-1)
     $  + HX5(0,0, -1,-1,1)
     $  + 2.0000000000000000d+00*HX5(0,0,-1,0,-1)
     $  + HX5(0,0, -1,0,1)
     $  + 6.0000000000000000d+00*HX5(0,0,0,-1,-1)
     $  + 2.0000000000000000d+00*HX5(0,0,0,-1,1)
     $  + 2.0000000000000000d+00*HX5(0,0,0,0,-1)
     $  + HX5(0,0,0,0,1) 
     $  + 3.0000000000000000d+00*HX5(0,0,0,1,-1)
     $  + HX5(0,0,1,0, -1)
      HY5(0,-1,-1,0,1) = 
     $  + 2.6781232869596824d+00
     $  + 2.0293560632083841d-01*HX1(0)
     $  - 9.0154267736969571d-01*HX1(0)*HX1(0)
     $  - 5.4831135561607547d-01*HX1(0)*HX1(0)*HX1(0)
     $  + 8.3333333333333333d-03*HX1(0)*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  - 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0)*HX2(0,-1)
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0)*HX3(0,-1,-1)
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0)*HX3(0,0,-1)
     $  + 3.2898681336964528d+00*HX1(0)*HX2(0,-1)
     $  - 5.0000000000000000d-01*HX1(0)*HX2(0,-1)*HX2(0,-1)
     $  + 1.8030853547393914d+00*HX2(0,-1)
     $  + HX2(0, -1)*HX3(0,0,-1)
     $  - HX2(0, -1)*HX3(0,0,1)
     $  - 3.2898681336964528d+00*HX3(0,-1,-1)
     $  - 3.2898681336964528d+00*HX3(0,0,-1)
     $  + HX5(0, -1,-1,0,1)
     $  - HX5(0,0, -1,0,-1)
     $  + HX5(0,0, -1,0,1)
     $  - 3.0000000000000000d+00*HX5(0,0,0,-1,-1)
     $  + 3.0000000000000000d+00*HX5(0,0,0,-1,1)
     $  - 2.0000000000000000d+00*HX5(0,0,0,0,-1)
     $  + HX5(0,0,0,0,1) 
     $  + 3.0000000000000000d+00*HX5(0,0,0,1,-1)
     $  + HX5(0,0,1,0, -1)
      HY5(0,-1,-1,1,-1) = 
     $  + 7.5170363885043517d+00
     $  + 7.2835884202641093d+00*HX1(0)
     $  + 2.7620719062289241d+00*HX1(0)*HX1(0)
     $  + 4.1123351671205660d-01*HX1(0)*HX1(0)*HX1(0)
     $  + 8.3333333333333333d-03*HX1(0)*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  - 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0)*HX2(0,-1)
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0)*HX3(0,-1,-1)
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0)*HX3(0,0,-1)
     $  - 2.4674011002723396d+00*HX1(0)*HX2(0,-1)
     $  - 5.0000000000000000d-01*HX1(0)*HX2(0,-1)*HX2(0,-1)
     $  + HX1(0) *HX4(0,-1,-1,1)
     $  - HX1(0) *HX4(0,-1,0,1)
     $  - HX1(0) *HX4(0,0,-1,1)
     $  + HX1(0) *HX4(0,0,0,1)
     $  - 5.5241438124578482d+00*HX2(0,-1)
     $  - HX2(0, -1)*HX3(0,-1,-1)
     $  + 2.0000000000000000d+00*HX2(0,-1)*HX3(0,0,-1)
     $  + 2.0000000000000000d+00*HX2(0,-1)*HX3(0,0,1)
     $  + 2.4674011002723396d+00*HX3(0,-1,-1)
     $  + 2.4674011002723396d+00*HX3(0,0,-1)
     $  - HX5(0, -1,-1,0,1)
     $  - HX5(0, -1,-1,1,-1)
     $  + 3.0000000000000000d+00*HX5(0,-1,0,-1,-1)
     $  - HX5(0, -1,0,-1,1)
     $  + HX5(0, -1,0,1,-1)
     $  + 6.0000000000000000d+00*HX5(0,0,-1,-1,-1)
     $  - 2.0000000000000000d+00*HX5(0,0,-1,-1,1)
     $  - 3.0000000000000000d+00*HX5(0,0,-1,0,-1)
     $  - HX5(0,0, -1,0,1)
     $  + HX5(0,0, -1,1,-1)
     $  - 9.0000000000000000d+00*HX5(0,0,0,-1,-1)
     $  - 3.0000000000000000d+00*HX5(0,0,0,-1,1)
     $  - 3.0000000000000000d+00*HX5(0,0,0,0,-1)
     $  - 4.0000000000000000d+00*HX5(0,0,0,0,1)
     $  - 7.0000000000000000d+00*HX5(0,0,0,1,-1)
     $  - 2.0000000000000000d+00*HX5(0,0,1,0,-1)
      HY5(0,-1,-1,1,1) = 
     $  - 6.8749092010391177d+00
     $  - 2.0293560632083841d+00*HX1(0)
     $  + 1.3810359531144620d+00*HX1(0)*HX1(0)
     $  + 8.2246703342411321d-01*HX1(0)*HX1(0)*HX1(0)
     $  - 8.3333333333333333d-03*HX1(0)*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0)*HX2(0,-1)
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0)*HX3(0,-1,-1)
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0)*HX3(0,0,-1)
     $  - 4.9348022005446793d+00*HX1(0)*HX2(0,-1)
     $  + 5.0000000000000000d-01*HX1(0)*HX2(0,-1)*HX2(0,-1)
     $  - HX1(0) *HX4(0,-1,-1,1)
     $  + HX1(0) *HX4(0,-1,0,1)
     $  + HX1(0) *HX4(0,0,-1,1)
     $  - HX1(0) *HX4(0,0,0,1)
     $  - 2.7620719062289241d+00*HX2(0,-1)
     $  - HX2(0, -1)*HX3(0,0,-1)
     $  - HX2(0, -1)*HX3(0,0,1)
     $  + 4.9348022005446793d+00*HX3(0,-1,-1)
     $  + 4.9348022005446793d+00*HX3(0,0,-1)
     $  - HX5(0, -1,-1,1,1)
     $  + HX5(0, -1,0,-1,1)
     $  + HX5(0, -1,0,1,1)
     $  + 2.0000000000000000d+00*HX5(0,0,-1,-1,1)
     $  + HX5(0,0, -1,0,-1)
     $  + HX5(0,0, -1,1,1)
     $  + 3.0000000000000000d+00*HX5(0,0,0,-1,-1)
     $  + 2.0000000000000000d+00*HX5(0,0,0,0,-1)
     $  + 3.0000000000000000d+00*HX5(0,0,0,0,1)
     $  + 3.0000000000000000d+00*HX5(0,0,0,1,-1)
     $  - HX5(0,0,0,1,1) 
     $  + HX5(0,0,1,0, -1)
      HY5(0,-1,0,-1,1) = 
     $  - 5.6642178198849788d+00
     $  - 3.4847341165810188d+00*HX1(0)
     $  - 1.2595007772794312d+00*HX1(0)*HX1(0)
     $  - 4.1123351671205660d-01*HX1(0)*HX1(0)*HX1(0)
     $  + 8.3333333333333333d-03*HX1(0)*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  - 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0)*HX2(0,-1)
     $  + HX1(0) *HX1(0)*HX3(0,0,-1)
     $  + 2.4674011002723396d+00*HX1(0)*HX2(0,-1)
     $  + 5.0000000000000000d-01*HX1(0)*HX2(0,-1)*HX2(0,-1)
     $  - 2.0000000000000000d+00*HX1(0)*HX4(0,0,-1,-1)
     $  - 4.0000000000000000d+00*HX1(0)*HX4(0,0,0,-1)
     $  + 2.5190015545588625d+00*HX2(0,-1)
     $  - 2.0000000000000000d+00*HX2(0,-1)*HX3(0,0,-1)
     $  - HX2(0, -1)*HX3(0,0,1)
     $  - 4.9348022005446793d+00*HX3(0,0,-1)
     $  + HX5(0, -1,0,-1,1)
     $  + 4.0000000000000000d+00*HX5(0,0,-1,0,-1)
     $  + 2.0000000000000000d+00*HX5(0,0,-1,0,1)
     $  + 1.2000000000000000d+01*HX5(0,0,0,-1,-1)
     $  + 2.0000000000000000d+00*HX5(0,0,0,-1,1)
     $  + 8.0000000000000000d+00*HX5(0,0,0,0,-1)
     $  + HX5(0,0,0,0,1) 
     $  + 3.0000000000000000d+00*HX5(0,0,0,1,-1)
     $  + HX5(0,0,1,0, -1)
      HY5(0,-1,0,1,-1) = 
     $  - 3.4903562040510891d-01
     $  + 3.0788629039393420d+00*HX1(0)
     $  + 2.1610434546491269d+00*HX1(0)*HX1(0)
     $  + 4.1123351671205660d-01*HX1(0)*HX1(0)*HX1(0)
     $  + 8.3333333333333333d-03*HX1(0)*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  - 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0)*HX2(0,-1)
     $  + HX1(0) *HX1(0)*HX3(0,0,-1)
     $  - 2.4674011002723396d+00*HX1(0)*HX2(0,-1)
     $  - HX1(0) *HX4(0,-1,0,1)
     $  - 3.0000000000000000d+00*HX1(0)*HX4(0,0,0,-1)
     $  + HX1(0) *HX4(0,0,0,1)
     $  - 4.3220869092982539d+00*HX2(0,-1)
     $  + HX2(0, -1)*HX3(0,0,-1)
     $  + 2.0000000000000000d+00*HX2(0,-1)*HX3(0,0,1)
     $  + 4.9348022005446793d+00*HX3(0,0,-1)
     $  + HX5(0, -1,0,1,-1)
     $  - 3.0000000000000000d+00*HX5(0,0,-1,0,-1)
     $  - 2.0000000000000000d+00*HX5(0,0,-1,0,1)
     $  - 6.0000000000000000d+00*HX5(0,0,0,-1,-1)
     $  - 6.0000000000000000d+00*HX5(0,0,0,-1,1)
     $  + 3.0000000000000000d+00*HX5(0,0,0,0,-1)
     $  - 4.0000000000000000d+00*HX5(0,0,0,0,1)
     $  - 7.0000000000000000d+00*HX5(0,0,0,1,-1)
     $  - 2.0000000000000000d+00*HX5(0,0,1,0,-1)
      HY5(0,-1,0,1,1) = 
     $  - 6.6327319875542747d+00
     $  - 3.8557765200959298d+00*HX1(0)
     $  + 6.0102845157979714d-01*HX1(0)*HX1(0)
     $  + 8.2246703342411321d-01*HX1(0)*HX1(0)*HX1(0)
     $  - 8.3333333333333333d-03*HX1(0)*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0)*HX2(0,-1)
     $  - HX1(0) *HX1(0)*HX3(0,0,-1)
     $  - 4.9348022005446793d+00*HX1(0)*HX2(0,-1)
     $  + HX1(0) *HX4(0,-1,0,1)
     $  + 3.0000000000000000d+00*HX1(0)*HX4(0,0,0,-1)
     $  - HX1(0) *HX4(0,0,0,1)
     $  - 1.2020569031595942d+00*HX2(0,-1)
     $  - HX2(0, -1)*HX3(0,0,1)
     $  + 9.8696044010893586d+00*HX3(0,0,-1)
     $  + HX5(0, -1,0,1,1)
     $  + 3.0000000000000000d+00*HX5(0,0,0,-1,1)
     $  - 4.0000000000000000d+00*HX5(0,0,0,0,-1)
     $  + 3.0000000000000000d+00*HX5(0,0,0,0,1)
     $  + 3.0000000000000000d+00*HX5(0,0,0,1,-1)
     $  - HX5(0,0,0,1,1) 
     $  + HX5(0,0,1,0, -1)
      HY5(0,-1,1,-1,-1) = 
     $  - 1.4571792178347941d+01
     $  - 7.2835884202641093d+00*HX1(0)
     $  - 1.3810359531144620d+00*HX1(0)*HX1(0)
     $  + 8.3333333333333333d-03*HX1(0)*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  - 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0)*HX2(0,-1)
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0)*HX3(0,-1,1)
     $  + HX1(0) *HX1(0)*HX3(0,0,-1)
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0)*HX3(0,0,1)
     $  + 5.0000000000000000d-01*HX1(0)*HX2(0,-1)*HX2(0,-1)
     $  + HX1(0) *HX4(0,-1,0,1)
     $  + HX1(0) *HX4(0,-1,1,-1)
     $  - 2.0000000000000000d+00*HX1(0)*HX4(0,0,-1,-1)
     $  + 2.0000000000000000d+00*HX1(0)*HX4(0,0,-1,1)
     $  - 4.0000000000000000d+00*HX1(0)*HX4(0,0,0,-1)
     $  - 3.0000000000000000d+00*HX1(0)*HX4(0,0,0,1)
     $  - HX1(0) *HX4(0,0,1,-1)
     $  + 2.7620719062289241d+00*HX2(0,-1)
     $  - HX2(0, -1)*HX3(0,0,-1)
     $  - HX2(0, -1)*HX3(0,0,1)
     $  - HX5(0, -1,0,-1,-1)
     $  - HX5(0, -1,0,1,-1)
     $  - HX5(0, -1,1,-1,-1)
     $  + HX5(0,0, -1,0,-1)
     $  - 2.0000000000000000d+00*HX5(0,0,-1,1,-1)
     $  + 7.0000000000000000d+00*HX5(0,0,0,-1,-1)
     $  + 7.0000000000000000d+00*HX5(0,0,0,0,-1)
     $  + 6.0000000000000000d+00*HX5(0,0,0,0,1)
     $  + 6.0000000000000000d+00*HX5(0,0,0,1,-1)
     $  + HX5(0,0,1, -1,-1)
     $  + HX5(0,0,1,0, -1)
      HY5(0,-1,1,-1,1) = 
     $  + 8.4053642852001231d+00
     $  + 3.3856186219460558d+00*HX1(0)
     $  + 6.5847232569963413d-01*HX1(0)*HX1(0)
     $  + 4.1123351671205660d-01*HX1(0)*HX1(0)*HX1(0)
     $  - 8.3333333333333333d-03*HX1(0)*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0)*HX2(0,-1)
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0)*HX3(0,-1,1)
     $  - HX1(0) *HX1(0)*HX3(0,0,-1)
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0)*HX3(0,0,1)
     $  - 2.4674011002723396d+00*HX1(0)*HX2(0,-1)
     $  - 5.0000000000000000d-01*HX1(0)*HX2(0,-1)*HX2(0,-1)
     $  - HX1(0) *HX4(0,-1,0,1)
     $  - HX1(0) *HX4(0,-1,1,-1)
     $  + 2.0000000000000000d+00*HX1(0)*HX4(0,0,-1,-1)
     $  - 2.0000000000000000d+00*HX1(0)*HX4(0,0,-1,1)
     $  + 4.0000000000000000d+00*HX1(0)*HX4(0,0,0,-1)
     $  + 3.0000000000000000d+00*HX1(0)*HX4(0,0,0,1)
     $  + HX1(0) *HX4(0,0,1,-1)
     $  - 1.3169446513992682d+00*HX2(0,-1)
     $  + HX2(0, -1)*HX3(0,-1,1)
     $  + 2.0000000000000000d+00*HX2(0,-1)*HX3(0,0,-1)
     $  + 2.0000000000000000d+00*HX2(0,-1)*HX3(0,0,1)
     $  - 2.4674011002723396d+00*HX3(0,-1,1)
     $  + 4.9348022005446793d+00*HX3(0,0,-1)
     $  + 2.4674011002723396d+00*HX3(0,0,1)
     $  - 3.0000000000000000d+00*HX5(0,-1,0,-1,1)
     $  - HX5(0, -1,1,-1,1)
     $  + HX5(0, -1,1,0,1)
     $  - 4.0000000000000000d+00*HX5(0,0,-1,-1,1)
     $  - 4.0000000000000000d+00*HX5(0,0,-1,0,-1)
     $  - 2.0000000000000000d+00*HX5(0,0,-1,0,1)
     $  - 1.2000000000000000d+01*HX5(0,0,0,-1,-1)
     $  - 2.0000000000000000d+00*HX5(0,0,0,-1,1)
     $  - 8.0000000000000000d+00*HX5(0,0,0,0,-1)
     $  - 7.0000000000000000d+00*HX5(0,0,0,0,1)
     $  - 9.0000000000000000d+00*HX5(0,0,0,1,-1)
     $  + HX5(0,0,1, -1,1)
     $  - 3.0000000000000000d+00*HX5(0,0,1,0,-1)
     $  - HX5(0,0,1,0,1) 
      HY5(0,-1,1,0,1) = 
     $  - 4.5667800194983164d-01
     $  - 3.4847341165810188d+00*HX1(0)
     $  - 1.2020569031595942d+00*HX1(0)*HX1(0)
     $  + 5.4831135561607547d-01*HX1(0)*HX1(0)*HX1(0)
     $  - 8.3333333333333333d-03*HX1(0)*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0)*HX2(0,-1)
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0)*HX3(0,-1,1)
     $  - HX1(0) *HX1(0)*HX3(0,0,-1)
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0)*HX3(0,0,1)
     $  - 3.2898681336964528d+00*HX1(0)*HX2(0,-1)
     $  - HX1(0) *HX4(0,-1,0,1)
     $  - 2.0000000000000000d+00*HX1(0)*HX4(0,0,-1,1)
     $  + 3.0000000000000000d+00*HX1(0)*HX4(0,0,0,-1)
     $  + 3.0000000000000000d+00*HX1(0)*HX4(0,0,0,1)
     $  + 2.4041138063191885d+00*HX2(0,-1)
     $  + 2.0000000000000000d+00*HX2(0,-1)*HX3(0,0,1)
     $  - 3.2898681336964528d+00*HX3(0,-1,1)
     $  + 6.5797362673929057d+00*HX3(0,0,-1)
     $  + 3.2898681336964528d+00*HX3(0,0,1)
     $  + HX5(0, -1,1,0,1)
     $  - 2.0000000000000000d+00*HX5(0,0,-1,0,1)
     $  - 3.0000000000000000d+00*HX5(0,0,0,-1,1)
     $  - 4.0000000000000000d+00*HX5(0,0,0,0,-1)
     $  - 7.0000000000000000d+00*HX5(0,0,0,0,1)
     $  - 6.0000000000000000d+00*HX5(0,0,0,1,-1)
     $  - 2.0000000000000000d+00*HX5(0,0,1,0,-1)
     $  - HX5(0,0,1,0,1) 
      HY5(0,-1,1,1,-1) = 
     $  + 8.0778313134277241d+00
     $  + 6.7309350447071233d-01*HX1(0)
     $  - 2.0395082788140962d+00*HX1(0)*HX1(0)
     $  - 4.1123351671205660d-01*HX1(0)*HX1(0)*HX1(0)
     $  - 8.3333333333333333d-03*HX1(0)*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0)*HX2(0,-1)
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0)*HX3(0,-1,1)
     $  - HX1(0) *HX1(0)*HX3(0,0,-1)
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0)*HX3(0,0,1)
     $  + 2.4674011002723396d+00*HX1(0)*HX2(0,-1)
     $  + HX1(0) *HX4(0,-1,1,1)
     $  - 2.0000000000000000d+00*HX1(0)*HX4(0,0,-1,1)
     $  + 3.0000000000000000d+00*HX1(0)*HX4(0,0,0,-1)
     $  + 2.0000000000000000d+00*HX1(0)*HX4(0,0,0,1)
     $  - HX1(0) *HX4(0,0,1,1)
     $  + 4.0790165576281924d+00*HX2(0,-1)
     $  - HX2(0, -1)*HX3(0,-1,1)
     $  - HX2(0, -1)*HX3(0,0,-1)
     $  - HX2(0, -1)*HX3(0,0,1)
     $  + 2.4674011002723396d+00*HX3(0,-1,1)
     $  - 4.9348022005446793d+00*HX3(0,0,-1)
     $  - 2.4674011002723396d+00*HX3(0,0,1)
     $  + 2.0000000000000000d+00*HX5(0,-1,0,-1,1)
     $  - HX5(0, -1,0,1,1)
     $  - HX5(0, -1,1,0,1)
     $  - HX5(0, -1,1,1,-1)
     $  + 4.0000000000000000d+00*HX5(0,0,-1,-1,1)
     $  + 3.0000000000000000d+00*HX5(0,0,-1,0,-1)
     $  + 2.0000000000000000d+00*HX5(0,0,-1,0,1)
     $  + 2.0000000000000000d+00*HX5(0,0,-1,1,-1)
     $  - 2.0000000000000000d+00*HX5(0,0,-1,1,1)
     $  + 6.0000000000000000d+00*HX5(0,0,0,-1,-1)
     $  + 6.0000000000000000d+00*HX5(0,0,0,-1,1)
     $  - 3.0000000000000000d+00*HX5(0,0,0,0,-1)
     $  - 2.0000000000000000d+00*HX5(0,0,0,0,1)
     $  + 4.0000000000000000d+00*HX5(0,0,0,1,-1)
     $  + 3.0000000000000000d+00*HX5(0,0,0,1,1)
     $  + 2.0000000000000000d+00*HX5(0,0,1,0,-1)
     $  + HX5(0,0,1,0,1) 
     $  + HX5(0,0,1,1, -1)
      HY5(0,-1,1,1,1) = 
     $  + 6.9367501878022481d+00
     $  + 6.4865749331714713d+00*HX1(0)
     $  - 8.2246703342411321d-01*HX1(0)*HX1(0)*HX1(0)
     $  + 8.3333333333333333d-03*HX1(0)*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  - 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0)*HX2(0,-1)
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0)*HX3(0,-1,1)
     $  + HX1(0) *HX1(0)*HX3(0,0,-1)
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0)*HX3(0,0,1)
     $  + 4.9348022005446793d+00*HX1(0)*HX2(0,-1)
     $  - HX1(0) *HX4(0,-1,1,1)
     $  + 2.0000000000000000d+00*HX1(0)*HX4(0,0,-1,1)
     $  - 3.0000000000000000d+00*HX1(0)*HX4(0,0,0,-1)
     $  - 2.0000000000000000d+00*HX1(0)*HX4(0,0,0,1)
     $  + HX1(0) *HX4(0,0,1,1)
     $  + 4.9348022005446793d+00*HX3(0,-1,1)
     $  - 9.8696044010893586d+00*HX3(0,0,-1)
     $  - 4.9348022005446793d+00*HX3(0,0,1)
     $  - HX5(0, -1,1,1,1)
     $  + 2.0000000000000000d+00*HX5(0,0,-1,1,1)
     $  - 3.0000000000000000d+00*HX5(0,0,0,-1,1)
     $  + 4.0000000000000000d+00*HX5(0,0,0,0,-1)
     $  + 3.0000000000000000d+00*HX5(0,0,0,0,1)
     $  - 2.0000000000000000d+00*HX5(0,0,0,1,1)
     $  + HX5(0,0,1,1,1) 
      HY5(0,0,-1,-1,1) = 
     $  - 9.1165552893555038d-01
     $  - 2.5209599327464717d+00*HX1(0)
     $  - 1.3810359531144620d+00*HX1(0)*HX1(0)
     $  - 4.1123351671205660d-01*HX1(0)*HX1(0)*HX1(0)
     $  + 8.3333333333333333d-03*HX1(0)*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0)*HX3(0,0,-1)
     $  + HX1(0) *HX4(0,0,-1,-1)
     $  + 2.0000000000000000d+00*HX1(0)*HX4(0,0,0,-1)
     $  + 2.4674011002723396d+00*HX3(0,0,-1)
     $  + HX5(0,0, -1,-1,1)
     $  - HX5(0,0, -1,0,-1)
     $  - HX5(0,0, -1,0,1)
     $  - 3.0000000000000000d+00*HX5(0,0,0,-1,-1)
     $  - HX5(0,0,0, -1,1)
     $  - 2.0000000000000000d+00*HX5(0,0,0,0,-1)
     $  + HX5(0,0,0,0,1) 
      HY5(0,0,-1,0,-1) = 
     $  - 8.1023197211680719d+00
     $  - 5.1410353601279064d+00*HX1(0)
     $  - 1.2020569031595942d+00*HX1(0)*HX1(0)
     $  - 2.7415567780803773d-01*HX1(0)*HX1(0)*HX1(0)
     $  - 8.3333333333333333d-03*HX1(0)*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0)*HX3(0,0,-1)
     $  - 3.0000000000000000d+00*HX1(0)*HX4(0,0,0,-1)
     $  + 1.6449340668482264d+00*HX3(0,0,-1)
     $  - HX5(0,0, -1,0,-1)
     $  + 7.0000000000000000d+00*HX5(0,0,0,0,-1)
      HY5(0,0,-1,0,1) = 
     $  + 2.0434339691042511d+00
     $  - 9.4703282949724591d-01*HX1(0)
     $  - 9.0154267736969571d-01*HX1(0)*HX1(0)
     $  - 5.4831135561607547d-01*HX1(0)*HX1(0)*HX1(0)
     $  + 8.3333333333333333d-03*HX1(0)*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0)*HX3(0,0,-1)
     $  + 3.0000000000000000d+00*HX1(0)*HX4(0,0,0,-1)
     $  + 3.2898681336964528d+00*HX3(0,0,-1)
     $  - HX5(0,0, -1,0,1)
     $  - 6.0000000000000000d+00*HX5(0,0,0,0,-1)
     $  + HX5(0,0,0,0,1) 
      HY5(0,0,-1,1,-1) = 
     $  + 1.2874316759375889d+01
     $  + 8.5266539820739622d+00*HX1(0)
     $  + 2.7620719062289241d+00*HX1(0)*HX1(0)
     $  + 4.1123351671205660d-01*HX1(0)*HX1(0)*HX1(0)
     $  + 8.3333333333333333d-03*HX1(0)*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0)*HX3(0,0,-1)
     $  - HX1(0) *HX4(0,0,-1,1)
     $  + 3.0000000000000000d+00*HX1(0)*HX4(0,0,0,-1)
     $  + HX1(0) *HX4(0,0,0,1)
     $  - 2.4674011002723396d+00*HX3(0,0,-1)
     $  + HX5(0,0, -1,0,-1)
     $  + HX5(0,0, -1,0,1)
     $  + HX5(0,0, -1,1,-1)
     $  + 3.0000000000000000d+00*HX5(0,0,0,-1,1)
     $  - 7.0000000000000000d+00*HX5(0,0,0,0,-1)
     $  - 4.0000000000000000d+00*HX5(0,0,0,0,1)
     $  - HX5(0,0,0,1, -1)
      HY5(0,0,-1,1,1) = 
     $  - 4.2205881177506120d+00
     $  - 6.2689427375197987d-01*HX1(0)
     $  + 1.3810359531144620d+00*HX1(0)*HX1(0)
     $  + 8.2246703342411321d-01*HX1(0)*HX1(0)*HX1(0)
     $  - 8.3333333333333333d-03*HX1(0)*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0)*HX3(0,0,-1)
     $  + HX1(0) *HX4(0,0,-1,1)
     $  - 3.0000000000000000d+00*HX1(0)*HX4(0,0,0,-1)
     $  - HX1(0) *HX4(0,0,0,1)
     $  - 4.9348022005446793d+00*HX3(0,0,-1)
     $  + HX5(0,0, -1,1,1)
     $  - 3.0000000000000000d+00*HX5(0,0,0,-1,1)
     $  + 6.0000000000000000d+00*HX5(0,0,0,0,-1)
     $  + 3.0000000000000000d+00*HX5(0,0,0,0,1)
     $  - HX5(0,0,0,1,1) 
      HY5(0,0,0,-1,1) = 
     $  - 4.8071216221292780d+00
     $  - 3.9234217222028759d+00*HX1(0)
     $  - 1.2595007772794312d+00*HX1(0)*HX1(0)
     $  - 4.1123351671205660d-01*HX1(0)*HX1(0)*HX1(0)
     $  + 8.3333333333333333d-03*HX1(0)*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  - HX1(0) *HX4(0,0,0,-1)
     $  - HX5(0,0,0, -1,1)
     $  + 4.0000000000000000d+00*HX5(0,0,0,0,-1)
     $  + HX5(0,0,0,0,1) 
      HY5(0,0,0,1,-1) = 
     $  + 5.2683829003001246d+00
     $  + 4.1940025306306604d+00*HX1(0)
     $  + 2.1610434546491269d+00*HX1(0)*HX1(0)
     $  + 4.1123351671205660d-01*HX1(0)*HX1(0)*HX1(0)
     $  + 8.3333333333333333d-03*HX1(0)*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + HX1(0) *HX4(0,0,0,1)
     $  - HX5(0,0,0,0, -1)
     $  - 4.0000000000000000d+00*HX5(0,0,0,0,1)
     $  - HX5(0,0,0,1, -1)
      HY5(0,0,1,-1,-1) = 
     $  - 9.4096904031705199d+00
     $  - 5.8027584430066521d+00*HX1(0)
     $  - 1.3810359531144620d+00*HX1(0)*HX1(0)
     $  + 8.3333333333333333d-03*HX1(0)*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0)*HX3(0,0,1)
     $  - HX1(0) *HX4(0,0,0,-1)
     $  - 3.0000000000000000d+00*HX1(0)*HX4(0,0,0,1)
     $  - HX1(0) *HX4(0,0,1,-1)
     $  + HX5(0,0,0, -1,-1)
     $  + 3.0000000000000000d+00*HX5(0,0,0,0,-1)
     $  + 6.0000000000000000d+00*HX5(0,0,0,0,1)
     $  + 3.0000000000000000d+00*HX5(0,0,0,1,-1)
     $  + HX5(0,0,1, -1,-1)
      HY5(0,0,1,-1,1) = 
     $  + 8.7596705482058876d+00
     $  + 4.3326514514433017d+00*HX1(0)
     $  + 6.5847232569963413d-01*HX1(0)*HX1(0)
     $  + 4.1123351671205660d-01*HX1(0)*HX1(0)*HX1(0)
     $  - 8.3333333333333333d-03*HX1(0)*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0)*HX3(0,0,1)
     $  + HX1(0) *HX4(0,0,0,-1)
     $  + 3.0000000000000000d+00*HX1(0)*HX4(0,0,0,1)
     $  + HX1(0) *HX4(0,0,1,-1)
     $  + 2.4674011002723396d+00*HX3(0,0,1)
     $  + HX5(0,0,0, -1,1)
     $  - 4.0000000000000000d+00*HX5(0,0,0,0,-1)
     $  - 7.0000000000000000d+00*HX5(0,0,0,0,1)
     $  - 3.0000000000000000d+00*HX5(0,0,0,1,-1)
     $  + HX5(0,0,1, -1,1)
     $  - HX5(0,0,1,0, -1)
     $  - HX5(0,0,1,0,1) 
      HY5(0,0,1,0,-1) = 
     $  - 6.8544356072335814d+00
     $  - 5.0057449559140141d+00*HX1(0)
     $  - 9.0154267736969571d-01*HX1(0)*HX1(0)
     $  + 2.7415567780803773d-01*HX1(0)*HX1(0)*HX1(0)
     $  + 8.3333333333333333d-03*HX1(0)*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0)*HX3(0,0,1)
     $  - 3.0000000000000000d+00*HX1(0)*HX4(0,0,0,1)
     $  + 1.6449340668482264d+00*HX3(0,0,1)
     $  - HX5(0,0,0,0, -1)
     $  + 6.0000000000000000d+00*HX5(0,0,0,0,1)
     $  - HX5(0,0,1,0, -1)
      HY5(0,0,1,1,-1) = 
     $  + 3.2887945813357094d+00
     $  + 1.5001934240460787d-01*HX1(0)
     $  - 2.0395082788140962d+00*HX1(0)*HX1(0)
     $  - 4.1123351671205660d-01*HX1(0)*HX1(0)*HX1(0)
     $  - 8.3333333333333333d-03*HX1(0)*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0)*HX3(0,0,1)
     $  + 2.0000000000000000d+00*HX1(0)*HX4(0,0,0,1)
     $  - HX1(0) *HX4(0,0,1,1)
     $  - 2.4674011002723396d+00*HX3(0,0,1)
     $  + HX5(0,0,0,0, -1)
     $  - 2.0000000000000000d+00*HX5(0,0,0,0,1)
     $  + HX5(0,0,0,1, -1)
     $  + 3.0000000000000000d+00*HX5(0,0,0,1,1)
     $  + HX5(0,0,1,0, -1)
     $  + HX5(0,0,1,0,1) 
     $  + HX5(0,0,1,1, -1)
      HY5(0,1,-1,-1,-1) = 
     $  + 7.4845537829591132d+00
     $  + 2.4278628067547031d+00*HX1(0)
     $  + 8.3333333333333333d-03*HX1(0)*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0)*HX2(0,1)
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0)*HX3(0,0,-1)
     $  - HX1(0) *HX1(0)*HX3(0,0,1)
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0)*HX3(0,1,-1)
     $  + HX1(0) *HX4(0,0,-1,-1)
     $  + 2.0000000000000000d+00*HX1(0)*HX4(0,0,0,-1)
     $  + 3.0000000000000000d+00*HX1(0)*HX4(0,0,0,1)
     $  + 2.0000000000000000d+00*HX1(0)*HX4(0,0,1,-1)
     $  + HX1(0) *HX4(0,1,-1,-1)
     $  - HX5(0,0, -1,-1,-1)
     $  - 2.0000000000000000d+00*HX5(0,0,0,-1,-1)
     $  - 3.0000000000000000d+00*HX5(0,0,0,0,-1)
     $  - 4.0000000000000000d+00*HX5(0,0,0,0,1)
     $  - 3.0000000000000000d+00*HX5(0,0,0,1,-1)
     $  - 2.0000000000000000d+00*HX5(0,0,1,-1,-1)
     $  - HX5(0,1, -1,-1,-1)
      HY5(0,1,-1,-1,1) = 
     $  + 7.4122526767078634d-02
     $  + 3.0440340948125761d+00*HX1(0)
     $  + 1.3810359531144620d+00*HX1(0)*HX1(0)
     $  + 4.1123351671205660d-01*HX1(0)*HX1(0)*HX1(0)
     $  - 8.3333333333333333d-03*HX1(0)*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  - 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0)*HX2(0,1)
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0)*HX3(0,0,-1)
     $  + HX1(0) *HX1(0)*HX3(0,0,1)
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0)*HX3(0,1,-1)
     $  + 2.4674011002723396d+00*HX1(0)*HX2(0,1)
     $  - HX1(0) *HX4(0,0,-1,-1)
     $  - 2.0000000000000000d+00*HX1(0)*HX4(0,0,0,-1)
     $  - 3.0000000000000000d+00*HX1(0)*HX4(0,0,0,1)
     $  - 2.0000000000000000d+00*HX1(0)*HX4(0,0,1,-1)
     $  - HX1(0) *HX4(0,1,-1,-1)
     $  + HX2(0, -1)*HX3(0,1,-1)
     $  + 2.7620719062289241d+00*HX2(0,1)
     $  - HX2(0,1) *HX3(0,-1,-1)
     $  - HX2(0,1) *HX3(0,0,-1)
     $  - HX2(0,1) *HX3(0,0,1)
     $  + HX2(0,1) *HX3(0,1,-1)
     $  - 2.4674011002723396d+00*HX3(0,0,-1)
     $  - 4.9348022005446793d+00*HX3(0,0,1)
     $  - 2.4674011002723396d+00*HX3(0,1,-1)
     $  + HX5(0, -1,-1,0,1)
     $  + HX5(0, -1,0,-1,1)
     $  + HX5(0,0, -1,-1,1)
     $  + HX5(0,0, -1,0,-1)
     $  + 2.0000000000000000d+00*HX5(0,0,-1,0,1)
     $  + 3.0000000000000000d+00*HX5(0,0,0,-1,-1)
     $  + 4.0000000000000000d+00*HX5(0,0,0,-1,1)
     $  + 2.0000000000000000d+00*HX5(0,0,0,0,-1)
     $  + 3.0000000000000000d+00*HX5(0,0,0,0,1)
     $  + 6.0000000000000000d+00*HX5(0,0,0,1,-1)
     $  + 6.0000000000000000d+00*HX5(0,0,0,1,1)
     $  - 2.0000000000000000d+00*HX5(0,0,1,-1,1)
     $  + 2.0000000000000000d+00*HX5(0,0,1,0,-1)
     $  + 3.0000000000000000d+00*HX5(0,0,1,0,1)
     $  - 4.0000000000000000d+00*HX5(0,0,1,1,-1)
     $  - HX5(0,1, -1,-1,1)
     $  - 2.0000000000000000d+00*HX5(0,1,0,1,-1)
      HY5(0,1,-1,1,-1) = 
     $  - 1.8385793928455058d+01
     $  - 9.4736868115712082d+00*HX1(0)
     $  - 2.7620719062289241d+00*HX1(0)*HX1(0)
     $  - 4.1123351671205660d-01*HX1(0)*HX1(0)*HX1(0)
     $  - 8.3333333333333333d-03*HX1(0)*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  - 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0)*HX2(0,1)
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0)*HX3(0,0,-1)
     $  + HX1(0) *HX1(0)*HX3(0,0,1)
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0)*HX3(0,1,-1)
     $  - HX1(0) *HX2(0,-1)*HX2(0,1)
     $  - 2.4674011002723396d+00*HX1(0)*HX2(0,1)
     $  - 5.0000000000000000d-01*HX1(0)*HX2(0,1)*HX2(0,1)
     $  + HX1(0) *HX4(0,-1,0,1)
     $  + 3.0000000000000000d+00*HX1(0)*HX4(0,0,-1,1)
     $  - 3.0000000000000000d+00*HX1(0)*HX4(0,0,0,-1)
     $  - 4.0000000000000000d+00*HX1(0)*HX4(0,0,0,1)
     $  + 2.0000000000000000d+00*HX1(0)*HX4(0,0,1,1)
     $  + HX1(0) *HX4(0,1,-1,1)
     $  - HX2(0, -1)*HX3(0,1,-1)
     $  - 5.5241438124578482d+00*HX2(0,1)
     $  + 2.0000000000000000d+00*HX2(0,1)*HX3(0,-1,-1)
     $  + 2.0000000000000000d+00*HX2(0,1)*HX3(0,0,-1)
     $  + 2.0000000000000000d+00*HX2(0,1)*HX3(0,0,1)
     $  - HX2(0,1) *HX3(0,1,-1)
     $  + 2.4674011002723396d+00*HX3(0,0,-1)
     $  + 4.9348022005446793d+00*HX3(0,0,1)
     $  + 2.4674011002723396d+00*HX3(0,1,-1)
     $  - 2.0000000000000000d+00*HX5(0,-1,-1,0,1)
     $  - 2.0000000000000000d+00*HX5(0,-1,0,-1,1)
     $  - HX5(0, -1,0,1,-1)
     $  - 4.0000000000000000d+00*HX5(0,0,-1,-1,1)
     $  - HX5(0,0, -1,0,-1)
     $  - 3.0000000000000000d+00*HX5(0,0,-1,0,1)
     $  - 3.0000000000000000d+00*HX5(0,0,-1,1,-1)
     $  - 9.0000000000000000d+00*HX5(0,0,0,-1,1)
     $  + 7.0000000000000000d+00*HX5(0,0,0,0,-1)
     $  + 8.0000000000000000d+00*HX5(0,0,0,0,1)
     $  - 2.0000000000000000d+00*HX5(0,0,0,1,-1)
     $  - 1.2000000000000000d+01*HX5(0,0,0,1,1)
     $  - 2.0000000000000000d+00*HX5(0,0,1,0,-1)
     $  - 4.0000000000000000d+00*HX5(0,0,1,0,1)
     $  + 4.0000000000000000d+00*HX5(0,0,1,1,-1)
     $  - HX5(0,1, -1,1,-1)
     $  + 3.0000000000000000d+00*HX5(0,1,0,1,-1)
      HY5(0,1,-1,1,1) = 
     $  + 7.0189712804378578d-01
     $  - 8.5393570350547734d-01*HX1(0)
     $  - 1.3810359531144620d+00*HX1(0)*HX1(0)
     $  - 8.2246703342411321d-01*HX1(0)*HX1(0)*HX1(0)
     $  + 8.3333333333333333d-03*HX1(0)*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0)*HX2(0,1)
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0)*HX3(0,0,-1)
     $  - HX1(0) *HX1(0)*HX3(0,0,1)
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0)*HX3(0,1,-1)
     $  + HX1(0) *HX2(0,-1)*HX2(0,1)
     $  - 4.9348022005446793d+00*HX1(0)*HX2(0,1)
     $  + 5.0000000000000000d-01*HX1(0)*HX2(0,1)*HX2(0,1)
     $  - HX1(0) *HX4(0,-1,0,1)
     $  - 3.0000000000000000d+00*HX1(0)*HX4(0,0,-1,1)
     $  + 3.0000000000000000d+00*HX1(0)*HX4(0,0,0,-1)
     $  + 4.0000000000000000d+00*HX1(0)*HX4(0,0,0,1)
     $  - 2.0000000000000000d+00*HX1(0)*HX4(0,0,1,1)
     $  - HX1(0) *HX4(0,1,-1,1)
     $  - 2.7620719062289241d+00*HX2(0,1)
     $  + HX2(0,1) *HX3(0,-1,1)
     $  - HX2(0,1) *HX3(0,0,-1)
     $  - HX2(0,1) *HX3(0,0,1)
     $  + 4.9348022005446793d+00*HX3(0,0,-1)
     $  + 9.8696044010893586d+00*HX3(0,0,1)
     $  + 4.9348022005446793d+00*HX3(0,1,-1)
     $  - 2.0000000000000000d+00*HX5(0,-1,0,1,1)
     $  - HX5(0, -1,1,0,1)
     $  + HX5(0,0, -1,0,1)
     $  - 5.0000000000000000d+00*HX5(0,0,-1,1,1)
     $  + 6.0000000000000000d+00*HX5(0,0,0,-1,1)
     $  - 6.0000000000000000d+00*HX5(0,0,0,0,-1)
     $  - 7.0000000000000000d+00*HX5(0,0,0,0,1)
     $  + 7.0000000000000000d+00*HX5(0,0,0,1,1)
     $  + HX5(0,0,1,0,1) 
     $  - HX5(0,1, -1,1,1)
     $  + HX5(0,1,0,1,1) 
      HY5(0,1,0,1,-1) = 
     $  - 1.0648543070067714d+01
     $  - 4.6326901362525175d+00*HX1(0)
     $  - 2.1610434546491269d+00*HX1(0)*HX1(0)
     $  - 4.1123351671205660d-01*HX1(0)*HX1(0)*HX1(0)
     $  - 8.3333333333333333d-03*HX1(0)*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  - 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0)*HX2(0,1)
     $  + HX1(0) *HX1(0)*HX3(0,0,1)
     $  - 2.4674011002723396d+00*HX1(0)*HX2(0,1)
     $  - 5.0000000000000000d-01*HX1(0)*HX2(0,1)*HX2(0,1)
     $  - 4.0000000000000000d+00*HX1(0)*HX4(0,0,0,1)
     $  + 2.0000000000000000d+00*HX1(0)*HX4(0,0,1,1)
     $  - 4.3220869092982539d+00*HX2(0,1)
     $  + HX2(0,1) *HX3(0,0,-1)
     $  + 2.0000000000000000d+00*HX2(0,1)*HX3(0,0,1)
     $  + 4.9348022005446793d+00*HX3(0,0,1)
     $  - HX5(0,0, -1,0,1)
     $  - 3.0000000000000000d+00*HX5(0,0,0,-1,1)
     $  + HX5(0,0,0,0, -1)
     $  + 8.0000000000000000d+00*HX5(0,0,0,0,1)
     $  - 2.0000000000000000d+00*HX5(0,0,0,1,-1)
     $  - 1.2000000000000000d+01*HX5(0,0,0,1,1)
     $  - 2.0000000000000000d+00*HX5(0,0,1,0,-1)
     $  - 4.0000000000000000d+00*HX5(0,0,1,0,1)
     $  + HX5(0,1,0,1, -1)
      HY5(0,1,1,-1,-1) = 
     $  + 3.7073850386582300d+00
     $  + 4.4002966535502479d+00*HX1(0)
     $  + 1.3810359531144620d+00*HX1(0)*HX1(0)
     $  - 8.3333333333333333d-03*HX1(0)*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  - 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0)*HX2(0,1)
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0)*HX3(0,0,1)
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0)*HX3(0,1,1)
     $  + HX1(0) *HX2(0,-1)*HX2(0,1)
     $  + 5.0000000000000000d-01*HX1(0)*HX2(0,1)*HX2(0,1)
     $  - HX1(0) *HX4(0,-1,0,1)
     $  - 2.0000000000000000d+00*HX1(0)*HX4(0,0,-1,1)
     $  + HX1(0) *HX4(0,0,0,-1)
     $  - HX1(0) *HX4(0,0,1,-1)
     $  + HX1(0) *HX4(0,1,1,-1)
     $  + 2.7620719062289241d+00*HX2(0,1)
     $  - HX2(0,1) *HX3(0,-1,-1)
     $  - HX2(0,1) *HX3(0,0,-1)
     $  - HX2(0,1) *HX3(0,0,1)
     $  + HX5(0, -1,-1,0,1)
     $  + HX5(0, -1,0,-1,1)
     $  + HX5(0, -1,0,1,-1)
     $  + 2.0000000000000000d+00*HX5(0,0,-1,-1,1)
     $  + HX5(0,0, -1,0,1)
     $  + 2.0000000000000000d+00*HX5(0,0,-1,1,-1)
     $  - HX5(0,0,0, -1,-1)
     $  + 3.0000000000000000d+00*HX5(0,0,0,-1,1)
     $  - 3.0000000000000000d+00*HX5(0,0,0,0,-1)
     $  - 2.0000000000000000d+00*HX5(0,0,0,0,1)
     $  + 3.0000000000000000d+00*HX5(0,0,0,1,1)
     $  + HX5(0,0,1, -1,-1)
     $  + HX5(0,0,1,0,1) 
     $  - 2.0000000000000000d+00*HX5(0,0,1,1,-1)
     $  - HX5(0,1,0,1, -1)
     $  - HX5(0,1,1, -1,-1)
      HY5(0,1,1,-1,1) = 
     $  - 1.0133791970287534d+01
     $  - 5.5757170132531547d+00*HX1(0)
     $  - 6.5847232569963413d-01*HX1(0)*HX1(0)
     $  - 4.1123351671205660d-01*HX1(0)*HX1(0)*HX1(0)
     $  + 8.3333333333333333d-03*HX1(0)*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0)*HX2(0,1)
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0)*HX3(0,0,1)
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0)*HX3(0,1,1)
     $  - HX1(0) *HX2(0,-1)*HX2(0,1)
     $  - 2.4674011002723396d+00*HX1(0)*HX2(0,1)
     $  - 5.0000000000000000d-01*HX1(0)*HX2(0,1)*HX2(0,1)
     $  + HX1(0) *HX4(0,-1,0,1)
     $  + 2.0000000000000000d+00*HX1(0)*HX4(0,0,-1,1)
     $  - HX1(0) *HX4(0,0,0,-1)
     $  + HX1(0) *HX4(0,0,1,-1)
     $  - HX1(0) *HX4(0,1,1,-1)
     $  + HX2(0, -1)*HX3(0,1,1)
     $  - 1.3169446513992682d+00*HX2(0,1)
     $  - 2.0000000000000000d+00*HX2(0,1)*HX3(0,-1,1)
     $  + 2.0000000000000000d+00*HX2(0,1)*HX3(0,0,-1)
     $  + 2.0000000000000000d+00*HX2(0,1)*HX3(0,0,1)
     $  + HX2(0,1) *HX3(0,1,1)
     $  + 2.4674011002723396d+00*HX3(0,0,1)
     $  - 2.4674011002723396d+00*HX3(0,1,1)
     $  + 3.0000000000000000d+00*HX5(0,-1,0,1,1)
     $  + 2.0000000000000000d+00*HX5(0,-1,1,0,1)
     $  - 2.0000000000000000d+00*HX5(0,0,-1,0,1)
     $  + 6.0000000000000000d+00*HX5(0,0,-1,1,1)
     $  - 7.0000000000000000d+00*HX5(0,0,0,-1,1)
     $  + 4.0000000000000000d+00*HX5(0,0,0,0,-1)
     $  + 3.0000000000000000d+00*HX5(0,0,0,0,1)
     $  - 3.0000000000000000d+00*HX5(0,0,0,1,-1)
     $  - 9.0000000000000000d+00*HX5(0,0,0,1,1)
     $  + HX5(0,0,1, -1,1)
     $  - HX5(0,0,1,0, -1)
     $  - 3.0000000000000000d+00*HX5(0,0,1,0,1)
     $  - 6.0000000000000000d+00*HX5(0,0,1,1,1)
     $  - 3.0000000000000000d+00*HX5(0,1,0,1,1)
     $  - HX5(0,1,1, -1,1)
      HY5(0,1,1,1,-1) = 
     $  + 7.7484048807851773d-01
     $  - 5.6922216412839301d-02*HX1(0)
     $  + 2.0395082788140962d+00*HX1(0)*HX1(0)
     $  + 4.1123351671205660d-01*HX1(0)*HX1(0)*HX1(0)
     $  + 8.3333333333333333d-03*HX1(0)*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + 1.6666666666666666d-01*HX1(0)*HX1(0)*HX1(0)*HX2(0,1)
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0)*HX3(0,0,1)
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0)*HX3(0,1,1)
     $  + 2.4674011002723396d+00*HX1(0)*HX2(0,1)
     $  + HX1(0) *HX4(0,0,0,1)
     $  - HX1(0) *HX4(0,0,1,1)
     $  + HX1(0) *HX4(0,1,1,1)
     $  - HX2(0, -1)*HX3(0,1,1)
     $  + 4.0790165576281924d+00*HX2(0,1)
     $  + HX2(0,1) *HX3(0,-1,1)
     $  - HX2(0,1) *HX3(0,0,-1)
     $  - HX2(0,1) *HX3(0,0,1)
     $  - HX2(0,1) *HX3(0,1,1)
     $  - 2.4674011002723396d+00*HX3(0,0,1)
     $  + 2.4674011002723396d+00*HX3(0,1,1)
     $  - HX5(0, -1,0,1,1)
     $  - HX5(0, -1,1,0,1)
     $  + HX5(0,0, -1,0,1)
     $  - 2.0000000000000000d+00*HX5(0,0,-1,1,1)
     $  + 3.0000000000000000d+00*HX5(0,0,0,-1,1)
     $  - HX5(0,0,0,0, -1)
     $  - 2.0000000000000000d+00*HX5(0,0,0,0,1)
     $  + 2.0000000000000000d+00*HX5(0,0,0,1,-1)
     $  + 6.0000000000000000d+00*HX5(0,0,0,1,1)
     $  + HX5(0,0,1,0, -1)
     $  + 2.0000000000000000d+00*HX5(0,0,1,0,1)
     $  + HX5(0,0,1,1, -1)
     $  + 4.0000000000000000d+00*HX5(0,0,1,1,1)
     $  + 2.0000000000000000d+00*HX5(0,1,0,1,1)
     $  - HX5(0,1,1,1, -1)
      Hi5(-1,-1,-1,-1,1) = 
     $  + 9.6181291076284771d-03
     $  - 5.5504108664821579d-02*HX1(-1)
     $  + 1.2011325347955035d-01*HX1(-1)*HX1(-1)
     $  - 1.1552453009332421d-01*HX1(-1)*HX1(-1)*HX1(-1)
     $  + 4.1666666666666666d-02*HX1(-1)*HX1(-1)*HX1(-1)*HX1(-1)
     $  - 1.6666666666666666d-01*HX1(-1)*HX1(-1)*HX1(-1)*HX1(0)
     $  + 3.4657359027997265d-01*HX1(-1)*HX1(-1)*HX1(0)
     $  + 2.5000000000000000d-01*HX1(-1)*HX1(-1)*HX1(0)*HX1(0)
     $  - 2.4022650695910071d-01*HX1(-1)*HX1(0)
     $  - 3.4657359027997265d-01*HX1(-1)*HX1(0)*HX1(0)
     $  - 1.6666666666666666d-01*HX1(-1)*HX1(0)*HX1(0)*HX1(0)
     $  + 5.5504108664821579d-02*HX1(0)
     $  + 1.2011325347955035d-01*HX1(0)*HX1(0)
     $  + 1.1552453009332421d-01*HX1(0)*HX1(0)*HX1(0)
     $  + 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0)
      Hi5(-1,-1,-1,1,1) = 
     $  + 2.3395367393886929d+00
     $  - 2.4532465311320902d+00*HX1(-1)
     $  + 9.4258028690366357d-01*HX1(-1)*HX1(-1)
     $  + 1.6666666666666666d-01*HX1(-1)*HX1(-1)*HX1(-1)*HX1(0)
     $  - 2.5000000000000000d-01*HX1(-1)*HX1(-1)*HX1(0)*HX1(0)
     $  - 5.0000000000000000d-01*HX1(-1)*HX1(-1)*HX2(0,-1)
     $  - 5.0000000000000000d-01*HX1(-1)*HX1(-1)*HX2(0,1)
     $  - 1.8851605738073271d+00*HX1(-1)*HX1(0)
     $  + 1.6666666666666666d-01*HX1(-1)*HX1(0)*HX1(0)*HX1(0)
     $  + HX1( -1)*HX3(0,-1,-1)
     $  + HX1( -1)*HX3(0,0,-1)
     $  + HX1( -1)*HX3(0,0,1)
     $  + HX1( -1)*HX3(0,1,-1)
     $  + 2.4532465311320902d+00*HX1(0)
     $  + 9.4258028690366357d-01*HX1(0)*HX1(0)
     $  - 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + HX4( -1,-1,-1,1)
     $  - HX4(0, -1,-1,-1)
     $  - HX4(0,0, -1,-1)
     $  - HX4(0,0,0, -1)
     $  - HX4(0,0,0,1) 
     $  - HX4(0,0,1, -1)
     $  - HX4(0,1, -1,-1)
      Hi5(-1,-1,1,-1,1) = 
     $  - 5.1782796940106865d+00
     $  + 3.1962209462999012d+00*HX1(-1)
     $  - 6.5146002367115732d-01*HX1(-1)*HX1(-1)
     $  - 3.4657359027997265d-01*HX1(-1)*HX1(-1)*HX1(0)
     $  - 2.5000000000000000d-01*HX1(-1)*HX1(-1)*HX1(0)*HX1(0)
     $  + 5.0000000000000000d-01*HX1(-1)*HX1(-1)*HX2(0,-1)
     $  + 5.0000000000000000d-01*HX1(-1)*HX1(-1)*HX2(0,1)
     $  + 1.3029200473423146d+00*HX1(-1)*HX1(0)
     $  + 3.4657359027997265d-01*HX1(-1)*HX1(0)*HX1(0)
     $  + 1.6666666666666666d-01*HX1(-1)*HX1(0)*HX1(0)*HX1(0)
     $  + HX1( -1)*HX1(0)*HX2(0,-1)
     $  + HX1( -1)*HX1(0)*HX2(0,1)
     $  + 6.9314718055994530d-01*HX1(-1)*HX2(0,-1)
     $  + 6.9314718055994530d-01*HX1(-1)*HX2(0,1)
     $  + HX1( -1)*HX3(-1,-1,1)
     $  - 2.0000000000000000d+00*HX1(-1)*HX3(0,-1,-1)
     $  - 2.0000000000000000d+00*HX1(-1)*HX3(0,0,-1)
     $  - 2.0000000000000000d+00*HX1(-1)*HX3(0,0,1)
     $  - 2.0000000000000000d+00*HX1(-1)*HX3(0,1,-1)
     $  - 3.1962209462999012d+00*HX1(0)
     $  - 6.5146002367115732d-01*HX1(0)*HX1(0)
     $  - 1.1552453009332421d-01*HX1(0)*HX1(0)*HX1(0)
     $  - 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  - HX1(0) *HX3(-1,-1,1)
     $  - HX1(0) *HX3(0,-1,-1)
     $  - HX1(0) *HX3(0,0,-1)
     $  - HX1(0) *HX3(0,0,1)
     $  - HX1(0) *HX3(0,1,-1)
     $  - 6.9314718055994530d-01*HX3(-1,-1,1)
     $  - 6.9314718055994530d-01*HX3(0,-1,-1)
     $  - 6.9314718055994530d-01*HX3(0,0,-1)
     $  - 6.9314718055994530d-01*HX3(0,0,1)
     $  - 6.9314718055994530d-01*HX3(0,1,-1)
     $  - 3.0000000000000000d+00*HX4(-1,-1,-1,1)
     $  + 3.0000000000000000d+00*HX4(0,-1,-1,-1)
     $  + 3.0000000000000000d+00*HX4(0,0,-1,-1)
     $  + 3.0000000000000000d+00*HX4(0,0,0,-1)
     $  + 3.0000000000000000d+00*HX4(0,0,0,1)
     $  + 3.0000000000000000d+00*HX4(0,0,1,-1)
     $  + 3.0000000000000000d+00*HX4(0,1,-1,-1)
      Hi5(-1,-1,1,1,1) = 
     $  + 3.1512791646720465d+00
     $  - 5.5504108664821579d-02*HX1(-1)
     $  - 8.2246703342411321d-01*HX1(-1)*HX1(-1)
     $  + 2.5000000000000000d-01*HX1(-1)*HX1(-1)*HX1(0)*HX1(0)
     $  + 1.6449340668482264d+00*HX1(-1)*HX1(0)
     $  - 1.6666666666666666d-01*HX1(-1)*HX1(0)*HX1(0)*HX1(0)
     $  - HX1( -1)*HX1(0)*HX2(0,-1)
     $  - HX1( -1)*HX1(0)*HX2(0,1)
     $  - HX1( -1)*HX3(0,-1,1)
     $  + HX1( -1)*HX3(0,0,-1)
     $  + HX1( -1)*HX3(0,0,1)
     $  - HX1( -1)*HX3(0,1,1)
     $  + 5.5504108664821579d-02*HX1(0)
     $  - 8.2246703342411321d-01*HX1(0)*HX1(0)
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
      Hi5(-1,1,-1,1,1) = 
     $  - 5.6883244754482793d+00
     $  + 2.9285842322233888d+00*HX1(-1)
     $  + 1.8851605738073271d+00*HX1(-1)*HX1(0)
     $  - 1.6666666666666666d-01*HX1(-1)*HX1(0)*HX1(0)*HX1(0)
     $  + HX1( -1)*HX1(0)*HX2(-1,1)
     $  + HX1( -1)*HX1(0)*HX2(0,-1)
     $  + HX1( -1)*HX1(0)*HX2(0,1)
     $  + 2.0000000000000000d+00*HX1(-1)*HX3(0,-1,1)
     $  - 2.0000000000000000d+00*HX1(-1)*HX3(0,0,-1)
     $  - 2.0000000000000000d+00*HX1(-1)*HX3(0,0,1)
     $  + 2.0000000000000000d+00*HX1(-1)*HX3(0,1,1)
     $  - 2.9285842322233888d+00*HX1(0)
     $  - 9.4258028690366357d-01*HX1(0)*HX1(0)
     $  + 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0)*HX2(-1,1)
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0)*HX2(0,-1)
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0)*HX2(0,1)
     $  - 2.0000000000000000d+00*HX1(0)*HX3(-1,-1,1)
     $  - 2.0000000000000000d+00*HX1(0)*HX3(0,-1,-1)
     $  - 2.0000000000000000d+00*HX1(0)*HX3(0,0,-1)
     $  - 2.0000000000000000d+00*HX1(0)*HX3(0,0,1)
     $  - 2.0000000000000000d+00*HX1(0)*HX3(0,1,-1)
     $  + 1.8851605738073271d+00*HX2(-1,1)
     $  + 5.0000000000000000d-01*HX2(-1,1)*HX2(-1,1)
     $  - HX2( -1,1)*HX2(0,-1)
     $  - HX2( -1,1)*HX2(0,1)
     $  - 1.8851605738073271d+00*HX2(0,-1)
     $  + 5.0000000000000000d-01*HX2(0,-1)*HX2(0,-1)
     $  + HX2(0, -1)*HX2(0,1)
     $  - 1.8851605738073271d+00*HX2(0,1)
     $  + 5.0000000000000000d-01*HX2(0,1)*HX2(0,1)
     $  - 2.0000000000000000d+00*HX4(-1,-1,1,1)
     $  - 2.0000000000000000d+00*HX4(0,-1,-1,1)
     $  - 2.0000000000000000d+00*HX4(0,-1,1,-1)
     $  + 2.0000000000000000d+00*HX4(0,0,-1,-1)
     $  - 2.0000000000000000d+00*HX4(0,0,-1,1)
     $  + 4.0000000000000000d+00*HX4(0,0,0,-1)
     $  + 4.0000000000000000d+00*HX4(0,0,0,1)
     $  + 2.0000000000000000d+00*HX4(0,0,1,-1)
     $  - 2.0000000000000000d+00*HX4(0,0,1,1)
     $  - 2.0000000000000000d+00*HX4(0,1,-1,1)
     $  - 2.0000000000000000d+00*HX4(0,1,1,-1)
      Hi5(-1,1,1,1,1) = 
     $  - 8.0212429617572516d-01
     $  - 1.6449340668482264d+00*HX1(-1)*HX1(0)
     $  + 1.6666666666666666d-01*HX1(-1)*HX1(0)*HX1(0)*HX1(0)
     $  + 8.2246703342411321d-01*HX1(0)*HX1(0)
     $  - 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0)*HX2(-1,1)
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0)*HX2(0,-1)
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0)*HX2(0,1)
     $  + HX1(0) *HX3(-1,1,1)
     $  - HX1(0) *HX3(0,-1,1)
     $  + HX1(0) *HX3(0,0,-1)
     $  + HX1(0) *HX3(0,0,1)
     $  - HX1(0) *HX3(0,1,1)
     $  - 1.6449340668482264d+00*HX2(-1,1)
     $  + 1.6449340668482264d+00*HX2(0,-1)
     $  + 1.6449340668482264d+00*HX2(0,1)
     $  + HX4( -1,1,1,1)
     $  - HX4(0, -1,1,1)
     $  + HX4(0,0, -1,1)
     $  - HX4(0,0,0, -1)
     $  - HX4(0,0,0,1) 
     $  + HX4(0,0,1,1) 
     $  - HX4(0,1,1,1) 
      Hi5(0,-1,-1,-1,1) = 
     $  + 5.2709719078152786d-01
     $  + 5.5504108664821579d-02*HX1(0)
     $  + 1.2011325347955035d-01*HX1(0)*HX1(0)
     $  + 1.1552453009332421d-01*HX1(0)*HX1(0)*HX1(0)
     $  + 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0)*HX2(0,-1)
     $  - 6.9314718055994530d-01*HX1(0)*HX2(0,-1)
     $  + HX1(0) *HX3(0,-1,-1)
     $  + HX1(0) *HX3(0,0,-1)
     $  - 2.4022650695910071d-01*HX2(0,-1)
     $  + 6.9314718055994530d-01*HX3(0,-1,-1)
     $  + 6.9314718055994530d-01*HX3(0,0,-1)
     $  - HX4(0, -1,-1,-1)
     $  - HX4(0,0, -1,-1)
     $  - HX4(0,0,0, -1)
      Hi5(0,-1,-1,0,1) = 
     $  - 3.3822601053473068d-01
     $  - 1.0517997902646449d+00*HX1(0)
     $  - 4.1123351671205660d-01*HX1(0)*HX1(0)
     $  + 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0)*HX2(0,-1)
     $  + HX1(0) *HX3(0,-1,-1)
     $  + HX1(0) *HX3(0,0,-1)
     $  + 8.2246703342411321d-01*HX2(0,-1)
     $  - 5.0000000000000000d-01*HX2(0,-1)*HX2(0,-1)
      Hi5(0,-1,-1,1,-1) = 
     $  + 3.3389529407850318d-01
     $  - 1.6651232599446473d-01*HX1(0)
     $  - 2.4022650695910071d-01*HX1(0)*HX1(0)
     $  - 1.1552453009332421d-01*HX1(0)*HX1(0)*HX1(0)
     $  + 6.9314718055994530d-01*HX1(0)*HX2(0,-1)
     $  + 4.8045301391820142d-01*HX2(0,-1)
     $  - 6.9314718055994530d-01*HX3(0,-1,-1)
     $  - 6.9314718055994530d-01*HX3(0,0,-1)
      Hi5(0,-1,-1,1,1) = 
     $  + 2.0079271680389850d+00
     $  + 2.4532465311320902d+00*HX1(0)
     $  + 9.4258028690366357d-01*HX1(0)*HX1(0)
     $  - 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0)*HX2(0,-1)
     $  - HX1(0) *HX3(0,-1,-1)
     $  - HX1(0) *HX3(0,0,-1)
     $  - 1.8851605738073271d+00*HX2(0,-1)
     $  + 5.0000000000000000d-01*HX2(0,-1)*HX2(0,-1)
     $  - HX4(0, -1,-1,1)
     $  + HX4(0, -1,0,1)
     $  + HX4(0,0, -1,1)
     $  - HX4(0,0,0,1) 
      Hi5(0,-1,0,-1,1) = 
     $  + 3.0521251417144869d+00
     $  + 1.5335088752078636d+00*HX1(0)
     $  + 4.1123351671205660d-01*HX1(0)*HX1(0)
     $  + 1.1552453009332421d-01*HX1(0)*HX1(0)*HX1(0)
     $  + 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0)*HX2(0,-1)
     $  - 6.9314718055994530d-01*HX1(0)*HX2(0,-1)
     $  + 2.0000000000000000d+00*HX1(0)*HX3(0,0,-1)
     $  - 8.2246703342411321d-01*HX2(0,-1)
     $  + 5.0000000000000000d-01*HX2(0,-1)*HX2(0,-1)
     $  + 1.3862943611198906d+00*HX3(0,0,-1)
     $  - 2.0000000000000000d+00*HX4(0,0,-1,-1)
     $  - 4.0000000000000000d+00*HX4(0,0,0,-1)
      Hi5(0,-1,0,1,-1) = 
     $  + 1.2498035299465379d+00
     $  + 5.7009070532142637d-01*HX1(0)
     $  - 1.1552453009332421d-01*HX1(0)*HX1(0)*HX1(0)
     $  + 6.9314718055994530d-01*HX1(0)*HX2(0,-1)
     $  - 1.3862943611198906d+00*HX3(0,0,-1)
      Hi5(0,-1,0,1,1) = 
     $  - 6.0880681896251523d-01
     $  + 1.0517997902646449d+00*HX1(0)
     $  + 8.2246703342411321d-01*HX1(0)*HX1(0)
     $  - 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0)*HX2(0,-1)
     $  - 2.0000000000000000d+00*HX1(0)*HX3(0,0,-1)
     $  - 1.6449340668482264d+00*HX2(0,-1)
     $  + HX4(0, -1,0,1)
     $  + 3.0000000000000000d+00*HX4(0,0,0,-1)
     $  - HX4(0,0,0,1) 
      Hi5(0,-1,1,-1,-1) = 
     $  + 1.9757838252848865d-01
     $  + 1.6651232599446473d-01*HX1(0)
     $  + 1.2011325347955035d-01*HX1(0)*HX1(0)
     $  - 2.4022650695910071d-01*HX2(0,-1)
      Hi5(0,-1,1,-1,1) = 
     $  - 5.7090374034238681d+00
     $  - 3.1962209462999012d+00*HX1(0)
     $  - 6.5146002367115732d-01*HX1(0)*HX1(0)
     $  - 1.1552453009332421d-01*HX1(0)*HX1(0)*HX1(0)
     $  - 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0)*HX2(0,-1)
     $  + 6.9314718055994530d-01*HX1(0)*HX2(0,-1)
     $  + HX1(0) *HX3(0,-1,1)
     $  - 2.0000000000000000d+00*HX1(0)*HX3(0,0,-1)
     $  - HX1(0) *HX3(0,0,1)
     $  + 1.3029200473423146d+00*HX2(0,-1)
     $  - 5.0000000000000000d-01*HX2(0,-1)*HX2(0,-1)
     $  + 6.9314718055994530d-01*HX3(0,-1,1)
     $  - 1.3862943611198906d+00*HX3(0,0,-1)
     $  - 6.9314718055994530d-01*HX3(0,0,1)
     $  - HX4(0, -1,0,1)
     $  - HX4(0, -1,1,-1)
     $  + 2.0000000000000000d+00*HX4(0,0,-1,-1)
     $  - 2.0000000000000000d+00*HX4(0,0,-1,1)
     $  + 4.0000000000000000d+00*HX4(0,0,0,-1)
     $  + 3.0000000000000000d+00*HX4(0,0,0,1)
     $  + HX4(0,0,1, -1)
      Hi5(0,-1,1,0,1) = 
     $  - 5.6493805840032952d+00
     $  - 2.6736902858507163d+00*HX1(0)
     $  - 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0)*HX2(0,-1)
     $  + HX1(0) *HX3(0,-1,1)
     $  - 2.0000000000000000d+00*HX1(0)*HX3(0,0,-1)
     $  - HX1(0) *HX3(0,0,1)
     $  - HX4(0, -1,0,1)
     $  - 2.0000000000000000d+00*HX4(0,0,-1,1)
     $  + 3.0000000000000000d+00*HX4(0,0,0,-1)
     $  + 3.0000000000000000d+00*HX4(0,0,0,1)
      Hi5(0,-1,1,1,-1) = 
     $  - 2.3933959928473851d+00
     $  - 1.7102721159642791d+00*HX1(0)
     $  - 2.9112026323250625d-01*HX1(0)*HX1(0)
     $  + 1.1552453009332421d-01*HX1(0)*HX1(0)*HX1(0)
     $  - 6.9314718055994530d-01*HX1(0)*HX2(0,-1)
     $  + 5.8224052646501250d-01*HX2(0,-1)
     $  - 6.9314718055994530d-01*HX3(0,-1,1)
     $  + 1.3862943611198906d+00*HX3(0,0,-1)
     $  + 6.9314718055994530d-01*HX3(0,0,1)
      Hi5(0,-1,1,1,1) = 
     $  + 3.2396052320380567d+00
     $  + 5.5504108664821579d-02*HX1(0)
     $  - 8.2246703342411321d-01*HX1(0)*HX1(0)
     $  + 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0)*HX2(0,-1)
     $  - HX1(0) *HX3(0,-1,1)
     $  + 2.0000000000000000d+00*HX1(0)*HX3(0,0,-1)
     $  + HX1(0) *HX3(0,0,1)
     $  + 1.6449340668482264d+00*HX2(0,-1)
     $  - HX4(0, -1,1,1)
     $  + 2.0000000000000000d+00*HX4(0,0,-1,1)
     $  - 3.0000000000000000d+00*HX4(0,0,0,-1)
     $  - 2.0000000000000000d+00*HX4(0,0,0,1)
     $  + HX4(0,0,1,1) 
      Hi5(0,0,-1,-1,1) = 
     $  - 1.3569495655898781d+00
     $  - 4.8170908494321862d-01*HX1(0)
     $  + 1.2011325347955035d-01*HX1(0)*HX1(0)
     $  + 1.1552453009332421d-01*HX1(0)*HX1(0)*HX1(0)
     $  + 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  - HX1(0) *HX3(0,0,-1)
     $  - 6.9314718055994530d-01*HX3(0,0,-1)
     $  + HX4(0,0, -1,-1)
     $  + 2.0000000000000000d+00*HX4(0,0,0,-1)
      Hi5(0,0,-1,0,1) = 
     $  - 2.8410984884917377d+00
     $  - 1.8030853547393914d+00*HX1(0)
     $  - 4.1123351671205660d-01*HX1(0)*HX1(0)
     $  + 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  - HX1(0) *HX3(0,0,-1)
     $  + 3.0000000000000000d+00*HX4(0,0,0,-1)
      Hi5(0,0,-1,1,-1) = 
     $  - 6.2490176497326899d-01
     $  - 5.7009070532142637d-01*HX1(0)
     $  - 2.4022650695910071d-01*HX1(0)*HX1(0)
     $  - 1.1552453009332421d-01*HX1(0)*HX1(0)*HX1(0)
     $  + 6.9314718055994530d-01*HX3(0,0,-1)
      Hi5(0,0,-1,1,1) = 
     $  + 3.8055457225523666d+00
     $  + 2.7620719062289241d+00*HX1(0)
     $  + 9.4258028690366357d-01*HX1(0)*HX1(0)
     $  - 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + HX1(0) *HX3(0,0,-1)
     $  + HX4(0,0, -1,1)
     $  - 3.0000000000000000d+00*HX4(0,0,0,-1)
     $  - HX4(0,0,0,1) 
      Hi5(0,0,0,-1,1) = 
     $  + 9.4703282949724591d-01
     $  + 9.0154267736969571d-01*HX1(0)
     $  + 4.1123351671205660d-01*HX1(0)*HX1(0)
     $  + 1.1552453009332421d-01*HX1(0)*HX1(0)*HX1(0)
     $  + 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  - HX4(0,0,0, -1)
      Hi5(0,0,0,1,-1) = 
     $  - 1.1552453009332421d-01*HX1(0)*HX1(0)*HX1(0)
      Hi5(0,0,1,-1,-1) = 
     $  + 1.2011325347955035d-01*HX1(0)*HX1(0)
      Hi5(0,0,1,-1,1) = 
     $  - 3.5336454555719528d+00
     $  - 2.6736902858507163d+00*HX1(0)
     $  - 6.5146002367115732d-01*HX1(0)*HX1(0)
     $  - 1.1552453009332421d-01*HX1(0)*HX1(0)*HX1(0)
     $  - 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  - HX1(0) *HX3(0,0,1)
     $  - 6.9314718055994530d-01*HX3(0,0,1)
     $  + HX4(0,0,0, -1)
     $  + 3.0000000000000000d+00*HX4(0,0,0,1)
     $  + HX4(0,0,1, -1)
      Hi5(0,0,1,0,-1) = 
     $  + 4.1123351671205660d-01*HX1(0)*HX1(0)
      Hi5(0,0,1,1,-1) = 
     $  - 8.3320235329769199d-01
     $  - 1.1401814106428527d+00*HX1(0)
     $  - 2.9112026323250625d-01*HX1(0)*HX1(0)
     $  + 1.1552453009332421d-01*HX1(0)*HX1(0)*HX1(0)
     $  + 6.9314718055994530d-01*HX3(0,0,1)
      Hi5(0,1,-1,-1,-1) = 
     $  - 5.5504108664821579d-02*HX1(0)
      Hi5(0,1,-1,-1,1) = 
     $  + 3.4402173672056757d+00
     $  + 1.0517997902646449d+00*HX1(0)
     $  - 1.2011325347955035d-01*HX1(0)*HX1(0)
     $  - 1.1552453009332421d-01*HX1(0)*HX1(0)*HX1(0)
     $  - 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  - 5.0000000000000000d-01*HX1(0)*HX1(0)*HX2(0,1)
     $  - 6.9314718055994530d-01*HX1(0)*HX2(0,1)
     $  + HX1(0) *HX3(0,0,-1)
     $  + 2.0000000000000000d+00*HX1(0)*HX3(0,0,1)
     $  + HX1(0) *HX3(0,1,-1)
     $  - 2.4022650695910071d-01*HX2(0,1)
     $  + 6.9314718055994530d-01*HX3(0,0,-1)
     $  + 1.3862943611198906d+00*HX3(0,0,1)
     $  + 6.9314718055994530d-01*HX3(0,1,-1)
     $  - HX4(0,0, -1,-1)
     $  - 2.0000000000000000d+00*HX4(0,0,0,-1)
     $  - 3.0000000000000000d+00*HX4(0,0,0,1)
     $  - 2.0000000000000000d+00*HX4(0,0,1,-1)
     $  - HX4(0,1, -1,-1)
      Hi5(0,1,-1,1,-1) = 
     $  + 1.8532608833279382d+00
     $  + 1.0926213657706112d+00*HX1(0)
     $  + 2.4022650695910071d-01*HX1(0)*HX1(0)
     $  + 1.1552453009332421d-01*HX1(0)*HX1(0)*HX1(0)
     $  + 6.9314718055994530d-01*HX1(0)*HX2(0,1)
     $  + 4.8045301391820142d-01*HX2(0,1)
     $  - 6.9314718055994530d-01*HX3(0,0,-1)
     $  - 1.3862943611198906d+00*HX3(0,0,1)
     $  - 6.9314718055994530d-01*HX3(0,1,-1)
      Hi5(0,1,-1,1,1) = 
     $  - 5.4839540823383124d+00
     $  - 2.9285842322233888d+00*HX1(0)
     $  - 9.4258028690366357d-01*HX1(0)*HX1(0)
     $  + 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0)*HX2(0,1)
     $  - HX1(0) *HX3(0,0,-1)
     $  - 2.0000000000000000d+00*HX1(0)*HX3(0,0,1)
     $  - HX1(0) *HX3(0,1,-1)
     $  + HX2(0, -1)*HX2(0,1)
     $  - 1.8851605738073271d+00*HX2(0,1)
     $  + 5.0000000000000000d-01*HX2(0,1)*HX2(0,1)
     $  - HX4(0, -1,0,1)
     $  - 3.0000000000000000d+00*HX4(0,0,-1,1)
     $  + 3.0000000000000000d+00*HX4(0,0,0,-1)
     $  + 4.0000000000000000d+00*HX4(0,0,0,1)
     $  - 2.0000000000000000d+00*HX4(0,0,1,1)
     $  - HX4(0,1, -1,1)
      Hi5(0,1,0,1,-1) = 
     $  + 1.6664047065953839d+00
     $  + 6.3196619783816790d-01*HX1(0)
     $  + 1.1552453009332421d-01*HX1(0)*HX1(0)*HX1(0)
     $  + 6.9314718055994530d-01*HX1(0)*HX2(0,1)
     $  - 1.3862943611198906d+00*HX3(0,0,1)
      Hi5(0,1,1,-1,-1) = 
     $  + 3.9515676505697730d-01
     $  + 3.0882537509683393d-01*HX1(0)
     $  - 1.2011325347955035d-01*HX1(0)*HX1(0)
     $  - 2.4022650695910071d-01*HX2(0,1)
      Hi5(0,1,1,-1,1) = 
     $  + 1.3317828347103277d+00
     $  + 2.2701119065237547d+00*HX1(0)
     $  + 6.5146002367115732d-01*HX1(0)*HX1(0)
     $  + 1.1552453009332421d-01*HX1(0)*HX1(0)*HX1(0)
     $  + 4.1666666666666666d-02*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + 5.0000000000000000d-01*HX1(0)*HX1(0)*HX2(0,1)
     $  + 6.9314718055994530d-01*HX1(0)*HX2(0,1)
     $  - HX1(0) *HX3(0,0,1)
     $  + HX1(0) *HX3(0,1,1)
     $  - HX2(0, -1)*HX2(0,1)
     $  + 1.3029200473423146d+00*HX2(0,1)
     $  - 5.0000000000000000d-01*HX2(0,1)*HX2(0,1)
     $  - 6.9314718055994530d-01*HX3(0,0,1)
     $  + 6.9314718055994530d-01*HX3(0,1,1)
     $  + HX4(0, -1,0,1)
     $  + 2.0000000000000000d+00*HX4(0,0,-1,1)
     $  - HX4(0,0,0, -1)
     $  + HX4(0,0,1, -1)
     $  - HX4(0,1,1, -1)
      Hi5(0,1,1,1,-1) = 
     $  - 9.5774727708194543d-01
     $  + 6.0296821703481255d-01*HX1(0)
     $  + 2.9112026323250625d-01*HX1(0)*HX1(0)
     $  - 1.1552453009332421d-01*HX1(0)*HX1(0)*HX1(0)
     $  - 6.9314718055994530d-01*HX1(0)*HX2(0,1)
     $  + 5.8224052646501250d-01*HX2(0,1)
     $  + 6.9314718055994530d-01*HX3(0,0,1)
     $  - 6.9314718055994530d-01*HX3(0,1,1)
      endif 
** nw > 4 endif 
      endif 
** (n1,n2) = (-1,1) -- completion endif 
      return 
      end 
************************************************************************ 
      subroutine fillirr1dhplin1(y,nw,HY1,HY2,HY3,HY4,HY5,n1,n2) 
** evaluates the irreducible HPL for y =1
** it is guaranteed that nw is in the range 2:4, and that (n1,n2) 
** take one of the pairs of values (0,1), (-1,0) or (-1,1) 
      implicit double precision (a-h,o-z) 
      dimension HY1(n1:n2),HY2(n1:n2,n1:n2),HY3(n1:n2,n1:n2,n1:n2), 
     $          HY4(n1:n2,n1:n2,n1:n2,n1:n2), 
     $          HY5(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2) 
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
      if (nw.gt.4) then
      HY5(0,0,0,0,1) = 
     $  + 1.0369277551433699d+00
      HY5(0,0,0,1,1) = 
     $  + 9.6551159989443734d-02
      HY5(0,0,1,0,1) = 
     $  + 2.2881039760335375d-01
      HY5(0,0,1,1,1) = 
     $  + 9.6551159989443734d-02
      HY5(0,1,0,1,1) = 
     $  + 2.2881039760335375d-01
      HY5(0,1,1,1,1) = 
     $  + 1.0369277551433699d+00
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
      if (nw.gt.4) then
      HY5(0,0,0,0,-1) = 
     $  + 9.7211977044690930d-01
      HY5(0,0,0,-1,-1) = 
     $  + 4.8936397049969063d-02
      HY5(0,0,-1,0,-1) = 
     $  + 9.2748467341632644d-02
      HY5(0,0,-1,-1,-1) = 
     $  + 9.6015684431298325d-03
      HY5(0,-1,0,-1,-1) = 
     $  + 1.3531263989594243d-02
      HY5(0,-1,-1,-1,-1) = 
     $  + 3.1350096016808622d-03
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
      if (nw.gt.4) then 
      HY5(-1,-1,-1,-1,1) = 
     $  + 1.8016537870380179d-03
      HY5(-1,-1,-1,1,1) = 
     $  + 3.8760673146652637d-03
      HY5(-1,-1,1,-1,1) = 
     $  + 6.2154684604081354d-03
      HY5(-1,-1,1,1,1) = 
     $  + 1.8530786065466613d-02
      HY5(-1,1,-1,1,1) = 
     $  + 3.8880058841843904d-02
      HY5(-1,1,1,1,1) = 
     $  + 5.0840057924226870d-01
      HY5(0,-1,-1,-1,1) = 
     $  + 4.1914400448554060d-03
      HY5(0,-1,-1,0,1) = 
     $  + 3.0172237496701167d-02
      HY5(0,-1,-1,1,-1) = 
     $  + 5.9459097989450212d-03
      HY5(0,-1,-1,1,1) = 
     $  + 8.7734377821481916d-03
      HY5(0,-1,0,-1,1) = 
     $  + 1.7042475614121991d-02
      HY5(0,-1,0,1,-1) = 
     $  + 2.2495758621687517d-02
      HY5(0,-1,0,1,1) = 
     $  + 3.0833054551948363d-02
      HY5(0,-1,1,-1,-1) = 
     $  + 9.4133341974174110d-03
      HY5(0,-1,1,-1,1) = 
     $  + 1.3833955759762555d-02
      HY5(0,-1,1,0,1) = 
     $  + 7.6026642213084631d-02
      HY5(0,-1,1,1,-1) = 
     $  + 2.2801059128486651d-02
      HY5(0,-1,1,1,1) = 
     $  + 3.9984858137537496d-02
      HY5(0,0,-1,-1,1) = 
     $  + 1.2444228784499648d-02
      HY5(0,0,-1,0,1) = 
     $  + 1.0679981350605469d-01
      HY5(0,0,-1,1,-1) = 
     $  + 1.6991592326175436d-02
      HY5(0,0,-1,1,1) = 
     $  + 2.4107342184124538d-02
      HY5(0,0,0,-1,1) = 
     $  + 5.9142607400864533d-02
      HY5(0,0,0,1,-1) = 
     $  + 7.4276054639867797d-02
      HY5(0,0,1,-1,-1) = 
     $  + 2.5535023438634211d-02
      HY5(0,0,1,-1,1) = 
     $  + 3.6321732111088421d-02
      HY5(0,0,1,0,-1) = 
     $  + 1.8615775173851248d-01
      HY5(0,0,1,1,-1) = 
     $  + 5.7353571803049304d-02
      HY5(0,1,-1,-1,-1) = 
     $  + 1.9555438852482933d-02
      HY5(0,1,-1,-1,1) = 
     $  + 2.8668668263701248d-02
      HY5(0,1,-1,1,-1) = 
     $  + 4.7069633474401836d-02
      HY5(0,1,-1,1,1) = 
     $  + 8.2208743029471844d-02
      HY5(0,1,0,1,-1) = 
     $  + 1.4122347902560834d-01
      HY5(0,1,1,-1,-1) = 
     $  + 1.0254618242743082d-01
      HY5(0,1,1,-1,1) = 
     $  + 1.7692012816527167d-01
      HY5(0,1,1,1,-1) = 
     $  + 4.0707197178331534d-01
      endif
      endif
** (n1,n2) = (-1,1) -- completion endif 
      return
      end
