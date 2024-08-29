      program moments
*
      implicit none
*
      double precision x
      integer m,j
      double precision r,rmax
      double precision theta,sigma,t
      double precision WGPLG
      double precision S121mx, Li21mx, lnx, lnx2, lnx3
      double precision x2, ln1mx, ln1mx2, ln1mx3, ln1px, ln1px2
      double precision Li2mx, Li31mx, Li3mx, S12mx
      double complex N,tmp
      double complex fN
      double precision fx, fxe
      double precision pi
      double complex PSI, DPSI
      double complex NM, N1, N2, NI, NMI, N1I, N2I
      double complex S1, S2, S3, S4, S5, S6
      double complex S1M, S2M, S3M, S11, S21, S31, S22
      double complex G0QG_HAT, G0GQ
      double complex POL1, POL2
      double complex V1, V2, V3
      double complex ACG3

      parameter(pi = 3.1415926535897932385d0)
      double precision CF
      parameter(CF = 4d0 / 3d0)
      double precision CA
      parameter(CA = 3d0)
      double precision TR
      parameter(TR = 1d0 / 2d0)
      double precision EMC, ZETA2, ZETA3, ZETA4, ZETA5, ZETA6
      parameter(EMC   = 0.57721566490153D0)
      parameter(ZETA2 = 1.644934066848226D0)
      parameter(ZETA3 = 1.202056903159594D0)
      parameter(ZETA4 = 1.082323233711138D0)
      parameter(ZETA5 = 1.036927755143370D0)
      parameter(ZETA6 = 1.017343061984449D0)
*
      x = 0.1d0
*
      t = - dlog(x)
      m = 103                  ! Must be odd
      r = 2d0 * m / 5d0 / t
      rmax = 10d0
      if(r.gt.rmax) r = rmax
      tmp = (0d0,0d0)
      do j=1,m-1
         theta = - pi + dble(j) * ( 2d0 * pi / dble(m) )
         sigma = theta + (theta/tan(theta)-1d0)/tan(theta)
         N     = r * theta * dcmplx(1d0/tan(theta),1d0)
*
         NM = N - 1.
         N1 = N + 1.
         N2 = N + 2.
         NI = 1./N
         NMI = 1./NM
         N1I = 1./N1
         N2I = 1./N2
*
         S1 = PSI(N1) + EMC
         S2 = ZETA2 - DPSI (N1,1)
         S3 = ZETA3 + 0.5 * DPSI (N1,2)
         S4 = ZETA4 - 1./6.D0 * DPSI (N1,3)
         S5 = ZETA5 + 1./24.D0 * DPSI (N1,4)
         S6 = ZETA6 - 1./120.D0 * DPSI (N1,5)
*
         S1M = S1 - NI
         S2M = S2 - NI*NI
         S3M = S3 - NI**3
         S11 = S1 + N1I
         S21 = S2 + N1I*N1I
         S31 = S3 + N1I**3
         S22 = S21 + N2I*N2I
*
         G0QG_HAT = - 8.D0 * TR * NM/N/N1
         G0GQ = - 4.D0 * CF * N2/N/N1
         POL1 = 12.D0*N**8 + 52.D0*N**7 + 60.D0*N**6 -25.D0*N**4 
     1        - 2.D0*N**3 + 3.D0*N**2 + 8.D0*N + 4.D0  
         POL2 = 2.*N**8 + 10.*N**7 + 22.*N**6 +36.*N**5 + 29.*N**4 
     1        + 4.*N**3 + 33.*N**2 + 12.*N + 4.D0
*
         V1 = (PSI((N1+1)/2.d0) - PSI(N1/2.d0))/2.d0
         V2 = (DPSI((N1+1)/2.d0,1) - DPSI(N1/2.d0,1))/4.d0
         V3 = (DPSI((N1+1)/2.d0,2) - DPSI(N1/2.d0,2))/8.d0       
*
*     APS2Hq(N) [THIS WORKS]
*
c         fN = - 4.D0 * CF * TR * N2/N**2/N1**2 * (NM * (2.D0*S2+ZETA2)
c     1        + (4.D0*N**3 - 4.D0*N**2 - 3.D0*N - 1.D0)/N**2/N1**2)
*
*     AS2polHg(N) [THIS WORKS]
*
         fN = CF*TR * (
     1         4./3.D0*NM/N/N1*(-4.*S3 + S1**3 +3.*S1*S2+ 6.*S1*ZETA2) 
     2        - 4.D0*(N**4 + 17.*N**3 + 43.*N**2 +33.*N
     3        + 2.D0)*S2/N**2/N1**2/N2 
     4        - 4.D0*(3.*N**2+3.*N-2.D0)*S1**2/N**2/N1/N2 
     5        - 2.D0*NM*(3.*N**2+3.*N+2.D0)*ZETA2/N**2/N1**2
     6        - 4.D0*(N**3-2.*N**2-22.*N-36.D0)*S1/N**2/N1/N2
     7        - 2.D0*POL1/N**4/N1**4/N2)
     8        + CA * TR * (
     9           4.D0*(N**2+4.*N+5.D0)*S1**2/N/N1**2/N2  
     1        + 4.D0*(7.*N**3+24.*N**2+15.*N-16.D0)*S2/N**2/N1**2/N2
     2        + 8.D0*NM*N2*ZETA2/N**2/N1**2 
     3        + 4.D0*(N**4+4.*N**3-N**2-10.*N+2.d0)*S1/N/N1**3/N2
     4        - 4.D0*POL2/N**4/N1**4/N2
     5        - 16.D0*NM/N/N1**2*V2 
     6        + 4.D0*NM/3./N/N1*(12.*ACG3(N)+3.*V3-8.*S3-S1**3
     7        - 9.*S1*S2 -12.*S1*V2 -12.*V1*ZETA2-3.*ZETA3))

         tmp = tmp + exp(t * N) * dcmplx(1d0, sigma) * fN
      enddo
*
      fx = dreal(r * tmp / m)
*
*     Definitions
*
      S121mx = wgplg(1, 2, 1d0 - x);
      Li21mx = wgplg(1, 1, 1d0 - x);
      lnx    = dlog(x);
      lnx2   = lnx * lnx;
      lnx3   = lnx * lnx2;
      x2     = x * x
      ln1mx  = dlog(1d0 - x)
      ln1mx2 = ln1mx * ln1mx
      ln1mx3 = ln1mx * ln1mx2
      ln1px  = dlog(1d0 + x)
      ln1px2 = ln1px * ln1px
      Li2mx  = wgplg(1, 1, - x)
      Li31mx = wgplg(2, 1, 1d0 - x)
      Li3mx  = wgplg(2, 1, - x)
      S121mx = wgplg(1, 2, 1d0 - x)
      S12mx  = wgplg(1, 2, - x)
*
*     APS2Hq(x) [THIS WORKS]
*
c      fxe = CF * TR * ( ( 1d0 + x ) * ( 32d0 * S121mx
c     1     + 16d0 * lnx * Li21mx - 24d0 * zeta2 * lnx
c     2     - 4d0 * lnx3 / 3d0 )
c     3     + 20d0 * ( 1d0 - x ) * ( 2d0 * Li21mx - 3d0 * zeta2 )
c     4     - ( 2d0 - 6d0 * x ) * lnx2 
c     5     - ( 12d0 + 60d0 * x ) * lnx - 72d0 * ( 1d0 - x ) )
*
*     AS2Hg(x) [THIS WORKS]
*
      fxe = - ( CF * TR * ( ( - 1 + 2 * x ) * ( 8 * zeta3
     1     + 8 * zeta2 * ln1mx + 4 * ln1mx3 / 3 - 8 * ln1mx * Li21mx
     2     + 4 * zeta2 * lnx - 4 * lnx * ln1mx2 + 2 * lnx3 / 3
     3     - 8 * lnx * Li21mx + 8 * Li31mx - 24 * S121mx )
     4     - ( 116 - 48 * x - 16 * x2 ) * Li21mx
     5     + ( 50 - 32 * x - 8 * x2 ) * zeta2
     6     - ( 72 - 16 * x - 8 * x2 ) * lnx * ln1mx 
     7     + ( 12 - 8 * x - 4 * x2 ) * ln1mx2
     8     - ( 5 - 8 * x - 4 * x2 ) * lnx2 - ( 64 - 60 * x ) * ln1mx
     9     - ( 16 + 50 * x ) * lnx - 22 + 46 * x )
     1     + CA * TR * ( ( - 1 + 2 * x ) * ( - 8 * zeta2 * ln1mx
     2     - 4 * ln1mx3 / 3 + 8 * ln1mx * Li21mx - 8 * Li31mx )
     3     + ( 1 + 2 * x) * ( - 4 * lnx3 / 3 - 8 * zeta2 * ln1px
     4     - 16 * ln1px * Li2mx - 8 * lnx * ln1px2 + 4 * lnx2 * ln1px
     5     + 8 * lnx * Li2mx - 8 * Li3mx - 16 * S12mx )
     6     + 16 * ( 1 + x ) * ( 4 * S121mx + 2 * lnx * Li21mx
     7     - 3 * zeta2 * lnx + Li2mx + lnx * ln1px )
     8     - 16 * ( 1 - x ) * zeta3
     9     + ( 100 - 112 * x - 8 * x2 ) * Li21mx
     1     - ( 132 - 144 * x - 4 * x2 ) * zeta2
     2     - 4 * x * ( 4 +  x ) * lnx * ln1mx 
     3     - ( 10 - 8 * x - 2 * x2 ) * ln1mx2
     4     - ( 6 + 2 * x2 ) * lnx2 + 4 * ln1mx
     5     - ( 56 + 148 * x ) * lnx - 204 + 212 * x )
     6     )
*
      write(6,*) x, fx, fxe, fx / fxe
*
      end
*
* =================================================================av==
*
* ..File psifcts.f
*
*
* ..The complex psi function,  PZI(Z),  and its m'th derivatives, 
*    DPZI(Z,M),  calculated from the asymtotic expansions. The 
*    functional equations are used for  |Im(Z)| < 10  to shift the 
*    argument to  Re(Z) >= 10  before applying the expansions.
*
* =====================================================================
*
*
       FUNCTION PSI (Z)
*
       IMPLICIT DOUBLE COMPLEX (A-Z)
       SUB = 0.D0
       ZZ = Z
*
* ---------------------------------------------------------------------
*
* ..Shift of the argument using the functional equation
*
       IF (DABS(DIMAG(ZZ)) .LT. 10.D0) THEN
*
  1      CONTINUE
         IF (DBLE(ZZ) .LT. 10.D0) THEN
           SUB = SUB - 1./ ZZ
           ZZ = ZZ + 1.
           GOTO 1
         END IF
*
       END IF
*
* ..Use of the asymtotic expansion (at the shifted argument)
*
       RZ = 1./ ZZ
       DZ = RZ * RZ
       PSI = SUB + LOG(ZZ) - 0.5 * RZ - DZ/5040.D0 * ( 420.+ DZ *
     1       ( - 42. + DZ * (20. - 21. * DZ) ) )
*
* ---------------------------------------------------------------------
*
       RETURN
       END
*
* =====================================================================
*
*
       FUNCTION DPSI (Z,M)
*
       IMPLICIT DOUBLE COMPLEX (A-Z)
       INTEGER M, K1, K2
       SUB = 0.D0
       ZZ = Z
*
* ---------------------------------------------------------------------
*
* ..Shift of the argument using the functional equations
*
       IF (DABS(DIMAG(ZZ)) .LT. 10.D0) THEN
*
  1      CONTINUE
         SUBM = -1./ZZ
         DO 10 K1 = 1, M
           SUBM = - SUBM * K1 / ZZ
 10      CONTINUE
*
         IF (DBLE(ZZ) .LT. 10.D0) THEN
           SUB = SUB + SUBM
           ZZ = ZZ + 1.
           GOTO 1
         END IF
*
       END IF
*
* ---------------------------------------------------------------------
*
* ..Expansion coefficients for the first derivative
*
       A1 =  1.D0
       A2 =  1./2.D0
       A3 =  1./6.D0
       A4 = -1./30.D0
       A5 =  1./42.D0
       A6 = -1./30.D0
       A7 =  5./66.D0
*
* ..Expansion coefficients for the higher derivatives
*
       IF (M .EQ. 1) GO TO 2
       DO 11 K2 = 2, M
         A1 = A1 * (K2-1.)
         A2 = A2 *  K2
         A3 = A3 * (K2+1.)
         A4 = A4 * (K2+3.)
         A5 = A5 * (K2+5.)
         A6 = A6 * (K2+7.)
         A7 = A7 * (K2+9.)
  11   CONTINUE
  2    CONTINUE 
*
* ..Use of the asymtotic expansion (at the shifted argument)
*
       RZ = 1./ ZZ
       DZ = RZ * RZ
       DPSI = SUB + (-1)**(M+1) * RZ**M * ( A1 + RZ * (A2 + RZ * 
     1        (A3 + DZ * (A4 + DZ * (A5 + DZ * (A6 + A7 * DZ ))))) )
*
* ---------------------------------------------------------------------
*
       RETURN
       END
*
* =================================================================av==
C   24/08/89 101231638  MEMBER NAME  WGPLG    (ZWPROD.S)    F77
      double precision FUNCTION WGPLG(N,P,X)
 
      INTEGER P,P1,NC(10),INDEX(31)
      DOUBLE PRECISION FCT(0:4),SGN(0:4),U(0:4),S1(4,4),C(4,4)
      DOUBLE PRECISION A(0:30,10)
      DOUBLE PRECISION X,X1,H,ALFA,R,Q,C1,C2,B0,B1,B2,ZERO,HALF
 
      COMPLEX*16 V(0:5),SK,SM
 
      DATA FCT /1.0D0,1.0D0,2.0D0,6.0D0,24.0D0/
      DATA SGN /1.0D0,-1.0D0,1.0D0,-1.0D0,1.0D0/
      DATA ZERO /0.0D0/, HALF /0.5D0/
      DATA C1 /1.33333 33333 333D0/, C2 /0.33333 33333 3333D0/
 
      DATA S1(1,1) /1.64493 40668 482D0/
      DATA S1(1,2) /1.20205 69031 596D0/
      DATA S1(1,3) /1.08232 32337 111D0/
      DATA S1(1,4) /1.03692 77551 434D0/
      DATA S1(2,1) /1.20205 69031 596D0/
      DATA S1(2,2) /2.70580 80842 778D-1/
      DATA S1(2,3) /9.65511 59989 444D-2/
      DATA S1(3,1) /1.08232 32337 111D0/
      DATA S1(3,2) /9.65511 59989 444D-2/
      DATA S1(4,1) /1.03692 77551 434D0/
 
      DATA C(1,1) / 1.64493 40668 482D0/
      DATA C(1,2) / 1.20205 69031 596D0/
      DATA C(1,3) / 1.08232 32337 111D0/
      DATA C(1,4) / 1.03692 77551 434D0/
      DATA C(2,1) / 0.00000 00000 000D0/
      DATA C(2,2) /-1.89406 56589 945D0/
      DATA C(2,3) /-3.01423 21054 407D0/
      DATA C(3,1) / 1.89406 56589 945D0/
      DATA C(3,2) / 3.01423 21054 407D0/
      DATA C(4,1) / 0.00000 00000 000D0/
 
      DATA INDEX /1,2,3,4,6*0,5,6,7,7*0,8,9,8*0,10/
 
      DATA NC /24,26,28,30,22,24,26,19,22,17/
 
      DATA A( 0,1) / .96753 21504 3498D0/
      DATA A( 1,1) / .16607 30329 2785D0/
      DATA A( 2,1) / .02487 93229 2423D0/
      DATA A( 3,1) / .00468 63619 5945D0/
      DATA A( 4,1) / .00100 16274 9616D0/
      DATA A( 5,1) / .00023 20021 9609D0/
      DATA A( 6,1) / .00005 68178 2272D0/
      DATA A( 7,1) / .00001 44963 0056D0/
      DATA A( 8,1) / .00000 38163 2946D0/
      DATA A( 9,1) / .00000 10299 0426D0/
      DATA A(10,1) / .00000 02835 7538D0/
      DATA A(11,1) / .00000 00793 8705D0/
      DATA A(12,1) / .00000 00225 3670D0/
      DATA A(13,1) / .00000 00064 7434D0/
      DATA A(14,1) / .00000 00018 7912D0/
      DATA A(15,1) / .00000 00005 5029D0/
      DATA A(16,1) / .00000 00001 6242D0/
      DATA A(17,1) / .00000 00000 4827D0/
      DATA A(18,1) / .00000 00000 1444D0/
      DATA A(19,1) / .00000 00000 0434D0/
      DATA A(20,1) / .00000 00000 0131D0/
      DATA A(21,1) / .00000 00000 0040D0/
      DATA A(22,1) / .00000 00000 0012D0/
      DATA A(23,1) / .00000 00000 0004D0/
      DATA A(24,1) / .00000 00000 0001D0/
 
      DATA A( 0,2) / .95180 88912 7832D0/
      DATA A( 1,2) / .43131 13184 6532D0/
      DATA A( 2,2) / .10002 25071 4905D0/
      DATA A( 3,2) / .02442 41559 5220D0/
      DATA A( 4,2) / .00622 51246 3724D0/
      DATA A( 5,2) / .00164 07883 1235D0/
      DATA A( 6,2) / .00044 40792 0265D0/
      DATA A( 7,2) / .00012 27749 4168D0/
      DATA A( 8,2) / .00003 45398 1284D0/
      DATA A( 9,2) / .00000 98586 9565D0/
      DATA A(10,2) / .00000 28485 6995D0/
      DATA A(11,2) / .00000 08317 0847D0/
      DATA A(12,2) / .00000 02450 3950D0/
      DATA A(13,2) / .00000 00727 6496D0/
      DATA A(14,2) / .00000 00217 5802D0/
      DATA A(15,2) / .00000 00065 4616D0/
      DATA A(16,2) / .00000 00019 8033D0/
      DATA A(17,2) / .00000 00006 0204D0/
      DATA A(18,2) / .00000 00001 8385D0/
      DATA A(19,2) / .00000 00000 5637D0/
      DATA A(20,2) / .00000 00000 1735D0/
      DATA A(21,2) / .00000 00000 0536D0/
      DATA A(22,2) / .00000 00000 0166D0/
      DATA A(23,2) / .00000 00000 0052D0/
      DATA A(24,2) / .00000 00000 0016D0/
      DATA A(25,2) / .00000 00000 0005D0/
      DATA A(26,2) / .00000 00000 0002D0/
 
      DATA A( 0,3) / .98161 02799 1365D0/
      DATA A( 1,3) / .72926 80632 0726D0/
      DATA A( 2,3) / .22774 71490 9321D0/
      DATA A( 3,3) / .06809 08329 6197D0/
      DATA A( 4,3) / .02013 70118 3064D0/
      DATA A( 5,3) / .00595 47848 0197D0/
      DATA A( 6,3) / .00176 76901 3959D0/
      DATA A( 7,3) / .00052 74821 8502D0/
      DATA A( 8,3) / .00015 82746 1460D0/
      DATA A( 9,3) / .00004 77492 2076D0/
      DATA A(10,3) / .00001 44792 0408D0/
      DATA A(11,3) / .00000 44115 4886D0/
      DATA A(12,3) / .00000 13500 3870D0/
      DATA A(13,3) / .00000 04148 1779D0/
      DATA A(14,3) / .00000 01279 3307D0/
      DATA A(15,3) / .00000 00395 9070D0/
      DATA A(16,3) / .00000 00122 9055D0/
      DATA A(17,3) / .00000 00038 2658D0/
      DATA A(18,3) / .00000 00011 9459D0/
      DATA A(19,3) / .00000 00003 7386D0/
      DATA A(20,3) / .00000 00001 1727D0/
      DATA A(21,3) / .00000 00000 3687D0/
      DATA A(22,3) / .00000 00000 1161D0/
      DATA A(23,3) / .00000 00000 0366D0/
      DATA A(24,3) / .00000 00000 0116D0/
      DATA A(25,3) / .00000 00000 0037D0/
      DATA A(26,3) / .00000 00000 0012D0/
      DATA A(27,3) / .00000 00000 0004D0/
      DATA A(28,3) / .00000 00000 0001D0/
 
      DATA A( 0,4) /1.06405 21184 614 D0/
      DATA A( 1,4) /1.06917 20744 981 D0/
      DATA A( 2,4) / .41527 19325 1768D0/
      DATA A( 3,4) / .14610 33293 6222D0/
      DATA A( 4,4) / .04904 73264 8784D0/
      DATA A( 5,4) / .01606 34086 0396D0/
      DATA A( 6,4) / .00518 88935 0790D0/
      DATA A( 7,4) / .00166 29871 7324D0/
      DATA A( 8,4) / .00053 05827 9969D0/
      DATA A( 9,4) / .00016 88702 9251D0/
      DATA A(10,4) / .00005 36832 8059D0/
      DATA A(11,4) / .00001 70592 3313D0/
      DATA A(12,4) / .00000 54217 4374D0/
      DATA A(13,4) / .00000 17239 4082D0/
      DATA A(14,4) / .00000 05485 3275D0/
      DATA A(15,4) / .00000 01746 7795D0/
      DATA A(16,4) / .00000 00556 7550D0/
      DATA A(17,4) / .00000 00177 6234D0/
      DATA A(18,4) / .00000 00056 7224D0/
      DATA A(19,4) / .00000 00018 1313D0/
      DATA A(20,4) / .00000 00005 8012D0/
      DATA A(21,4) / .00000 00001 8579D0/
      DATA A(22,4) / .00000 00000 5955D0/
      DATA A(23,4) / .00000 00000 1911D0/
      DATA A(24,4) / .00000 00000 0614D0/
      DATA A(25,4) / .00000 00000 0197D0/
      DATA A(26,4) / .00000 00000 0063D0/
      DATA A(27,4) / .00000 00000 0020D0/
      DATA A(28,4) / .00000 00000 0007D0/
      DATA A(29,4) / .00000 00000 0002D0/
      DATA A(30,4) / .00000 00000 0001D0/
 
      DATA A( 0,5) / .97920 86066 9175D0/
      DATA A( 1,5) / .08518 81314 8683D0/
      DATA A( 2,5) / .00855 98522 2013D0/
      DATA A( 3,5) / .00121 17721 4413D0/
      DATA A( 4,5) / .00020 72276 8531D0/
      DATA A( 5,5) / .00003 99695 8691D0/
      DATA A( 6,5) / .00000 83806 4065D0/
      DATA A( 7,5) / .00000 18684 8945D0/
      DATA A( 8,5) / .00000 04366 6087D0/
      DATA A( 9,5) / .00000 01059 1733D0/
      DATA A(10,5) / .00000 00264 7892D0/
      DATA A(11,5) / .00000 00067 8700D0/
      DATA A(12,5) / .00000 00017 7654D0/
      DATA A(13,5) / .00000 00004 7342D0/
      DATA A(14,5) / .00000 00001 2812D0/
      DATA A(15,5) / .00000 00000 3514D0/
      DATA A(16,5) / .00000 00000 0975D0/
      DATA A(17,5) / .00000 00000 0274D0/
      DATA A(18,5) / .00000 00000 0077D0/
      DATA A(19,5) / .00000 00000 0022D0/
      DATA A(20,5) / .00000 00000 0006D0/
      DATA A(21,5) / .00000 00000 0002D0/
      DATA A(22,5) / .00000 00000 0001D0/
 
      DATA A( 0,6) / .95021 85196 3952D0/
      DATA A( 1,6) / .29052 52916 1433D0/
      DATA A( 2,6) / .05081 77406 1716D0/
      DATA A( 3,6) / .00995 54376 7280D0/
      DATA A( 4,6) / .00211 73389 5031D0/
      DATA A( 5,6) / .00047 85947 0550D0/
      DATA A( 6,6) / .00011 33432 1308D0/
      DATA A( 7,6) / .00002 78473 3104D0/
      DATA A( 8,6) / .00000 70478 8108D0/
      DATA A( 9,6) / .00000 18278 8740D0/
      DATA A(10,6) / .00000 04838 7492D0/
      DATA A(11,6) / .00000 01303 3842D0/
      DATA A(12,6) / .00000 00356 3769D0/
      DATA A(13,6) / .00000 00098 7174D0/
      DATA A(14,6) / .00000 00027 6586D0/
      DATA A(15,6) / .00000 00007 8279D0/
      DATA A(16,6) / .00000 00002 2354D0/
      DATA A(17,6) / .00000 00000 6435D0/
      DATA A(18,6) / .00000 00000 1866D0/
      DATA A(19,6) / .00000 00000 0545D0/
      DATA A(20,6) / .00000 00000 0160D0/
      DATA A(21,6) / .00000 00000 0047D0/
      DATA A(22,6) / .00000 00000 0014D0/
      DATA A(23,6) / .00000 00000 0004D0/
      DATA A(24,6) / .00000 00000 0001D0/
 
      DATA A( 0,7) / .95064 03218 6777D0/
      DATA A( 1,7) / .54138 28546 5171D0/
      DATA A( 2,7) / .13649 97959 0321D0/
      DATA A( 3,7) / .03417 94232 8207D0/
      DATA A( 4,7) / .00869 02788 3583D0/
      DATA A( 5,7) / .00225 28408 4155D0/
      DATA A( 6,7) / .00059 51608 9806D0/
      DATA A( 7,7) / .00015 99561 7766D0/
      DATA A( 8,7) / .00004 36521 3096D0/
      DATA A( 9,7) / .00001 20747 4688D0/
      DATA A(10,7) / .00000 33801 8176D0/
      DATA A(11,7) / .00000 09563 2476D0/
      DATA A(12,7) / .00000 02731 3129D0/
      DATA A(13,7) / .00000 00786 6968D0/
      DATA A(14,7) / .00000 00228 3195D0/
      DATA A(15,7) / .00000 00066 7205D0/
      DATA A(16,7) / .00000 00019 6191D0/
      DATA A(17,7) / .00000 00005 8018D0/
      DATA A(18,7) / .00000 00001 7246D0/
      DATA A(19,7) / .00000 00000 5151D0/
      DATA A(20,7) / .00000 00000 1545D0/
      DATA A(21,7) / .00000 00000 0465D0/
      DATA A(22,7) / .00000 00000 0141D0/
      DATA A(23,7) / .00000 00000 0043D0/
      DATA A(24,7) / .00000 00000 0013D0/
      DATA A(25,7) / .00000 00000 0004D0/
      DATA A(26,7) / .00000 00000 0001D0/
 
      DATA A( 0,8) / .98800 01167 2229D0/
      DATA A( 1,8) / .04364 06760 9601D0/
      DATA A( 2,8) / .00295 09117 8278D0/
      DATA A( 3,8) / .00031 47780 9720D0/
      DATA A( 4,8) / .00004 31484 6029D0/
      DATA A( 5,8) / .00000 69381 8230D0/
      DATA A( 6,8) / .00000 12464 0350D0/
      DATA A( 7,8) / .00000 02429 3628D0/
      DATA A( 8,8) / .00000 00504 0827D0/
      DATA A( 9,8) / .00000 00109 9075D0/
      DATA A(10,8) / .00000 00024 9467D0/
      DATA A(11,8) / .00000 00005 8540D0/
      DATA A(12,8) / .00000 00001 4127D0/
      DATA A(13,8) / .00000 00000 3492D0/
      DATA A(14,8) / .00000 00000 0881D0/
      DATA A(15,8) / .00000 00000 0226D0/
      DATA A(16,8) / .00000 00000 0059D0/
      DATA A(17,8) / .00000 00000 0016D0/
      DATA A(18,8) / .00000 00000 0004D0/
      DATA A(19,8) / .00000 00000 0001D0/
 
      DATA A( 0,9) / .95768 50654 6350D0/
      DATA A( 1,9) / .19725 24967 9534D0/
      DATA A( 2,9) / .02603 37031 3918D0/
      DATA A( 3,9) / .00409 38216 8261D0/
      DATA A( 4,9) / .00072 68170 7110D0/
      DATA A( 5,9) / .00014 09187 9261D0/
      DATA A( 6,9) / .00002 92045 8914D0/
      DATA A( 7,9) / .00000 63763 1144D0/
      DATA A( 8,9) / .00000 14516 7850D0/
      DATA A( 9,9) / .00000 03420 5281D0/
      DATA A(10,9) / .00000 00829 4302D0/
      DATA A(11,9) / .00000 00206 0784D0/
      DATA A(12,9) / .00000 00052 2823D0/
      DATA A(13,9) / .00000 00013 5066D0/
      DATA A(14,9) / .00000 00003 5451D0/
      DATA A(15,9) / .00000 00000 9436D0/
      DATA A(16,9) / .00000 00000 2543D0/
      DATA A(17,9) / .00000 00000 0693D0/
      DATA A(18,9) / .00000 00000 0191D0/
      DATA A(19,9) / .00000 00000 0053D0/
      DATA A(20,9) / .00000 00000 0015D0/
      DATA A(21,9) / .00000 00000 0004D0/
      DATA A(22,9) / .00000 00000 0001D0/
 
      DATA A( 0,10) / .99343 65167 1347D0/
      DATA A( 1,10) / .02225 77012 6826D0/
      DATA A( 2,10) / .00101 47557 4703D0/
      DATA A( 3,10) / .00008 17515 6250D0/
      DATA A( 4,10) / .00000 89997 3547D0/
      DATA A( 5,10) / .00000 12082 3987D0/
      DATA A( 6,10) / .00000 01861 6913D0/
      DATA A( 7,10) / .00000 00317 4723D0/
      DATA A( 8,10) / .00000 00058 5215D0/
      DATA A( 9,10) / .00000 00011 4739D0/
      DATA A(10,10) / .00000 00002 3652D0/
      DATA A(11,10) / .00000 00000 5082D0/
      DATA A(12,10) / .00000 00000 1131D0/
      DATA A(13,10) / .00000 00000 0259D0/
      DATA A(14,10) / .00000 00000 0061D0/
      DATA A(15,10) / .00000 00000 0015D0/
      DATA A(16,10) / .00000 00000 0004D0/
      DATA A(17,10) / .00000 00000 0001D0/
 
      IF(N .LT. 1 .OR. N .GT. 4 .OR. P .LT. 1 .OR. P .GT. 4 .OR.
     1   N+P .GT. 5) THEN
       WGPLG=ZERO
       PRINT 1000, N,P
       RETURN
      END IF
      IF(X .EQ. SGN(0)) THEN
       WGPLG=S1(N,P)
       RETURN
      END IF
 
      IF(X .GT. FCT(2) .OR. X .LT. SGN(1)) THEN
       X1=SGN(0)/X
       H=C1*X1+C2
       ALFA=H+H
       V(0)=SGN(0)
       V(1)=LOG(DCMPLX(-X,ZERO))
       DO 33 L = 2,N+P
   33  V(L)=V(1)*V(L-1)/L
       SK=ZERO
       DO 34 K = 0,P-1
       P1=P-K
       R=X1**P1/(FCT(P1)*FCT(N-1))
       SM=ZERO
       DO 35 M = 0,K
       N1=N+K-M
       L=INDEX(10*N1+P1-10)
       B1=ZERO
       B2=ZERO
       DO 31 I = NC(L),0,-1
       B0=A(I,L)+ALFA*B1-B2
       B2=B1
   31  B1=B0
       Q=(FCT(N1-1)/FCT(K-M))*(B0-H*B2)*R/P1**N1
   35  SM=SM+V(M)*Q
   34  SK=SK+SGN(K)*SM
       SM=ZERO
       DO 36 M = 0,N-1
   36  SM=SM+V(M)*C(N-M,P)
       WGPLG=SGN(N)*SK+SGN(P)*(SM+V(N+P))
       RETURN
      END IF
 
      IF(X .GT. HALF) THEN
       X1=SGN(0)-X
       H=C1*X1+C2
       ALFA=H+H
       V(0)=SGN(0)
       U(0)=SGN(0)
       V(1)=LOG(DCMPLX(X1,ZERO))
       U(1)=LOG(X)
       DO 23 L = 2,P
   23  V(L)=V(1)*V(L-1)/L
       DO 26 L = 2,N
   26  U(L)=U(1)*U(L-1)/L
       SK=ZERO
       DO 24 K = 0,N-1
       P1=N-K
       R=X1**P1/FCT(P1)
       SM=ZERO
       DO 25 M = 0,P-1
       N1=P-M
       L=INDEX(10*N1+P1-10)
       B1=ZERO
       B2=ZERO
       DO 12 I = NC(L),0,-1
       B0=A(I,L)+ALFA*B1-B2
       B2=B1
   12  B1=B0
       Q=SGN(M)*(B0-H*B2)*R/P1**N1
   25  SM=SM+V(M)*Q
   24  SK=SK+U(K)*(S1(P1,P)-SM)
       WGPLG=SK+SGN(P)*U(N)*V(P)
       RETURN
      END IF
 
      L=INDEX(10*N+P-10)
      H=C1*X+C2
      ALFA=H+H
      B1=ZERO
      B2=ZERO
      DO 11 I = NC(L),0,-1
      B0=A(I,L)+ALFA*B1-B2
      B2=B1
   11 B1=B0
      WGPLG=(B0-H*B2)*X**P/(FCT(P)*P**N)
      RETURN
 1000 FORMAT(/' ***** CERN SUBROUTINE WGPLG ... ILLEGAL VALUES',
     1        '   N = ',I3,'   P = ',I3)
      END
C     ----------------------------
C************************************************************************
      COMPLEX*16 FUNCTION ACG3(ZN1)
C     ----------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: LI2(X)/(1+X)
C---  FOR COMPLEX ARGUMENT
C
C************************************************************************
C
      IMPLICIT COMPLEX*16 (A-H,O-Z)
      DOUBLE PRECISION ZET2
      parameter(ZET2 = 1.644934066848226D0)
      DOUBLE PRECISION DL, G3
      PARAMETER(DL = DLOG(2.0D0))
      PARAMETER(GE = 0.57721566490153D+0)
      DOUBLE PRECISION AK1(9)
      DATA AK1/0.999999974532240D+0,
     &        -0.499995525890027D+0,
     &         0.333203435554182D+0,
     &        -0.248529457735332D+0,
     &         0.191451164493502D+0,
     &        -0.137466222203386D+0,
     &         0.792107405737825D-1,
     &        -0.301109652783781D-1,
     &         0.538406198111749D-2/
C
      ZN=ZN1+1D0
C
10    T=DCMPLX(DL*ZET2,0D0)
      Z=ZN1
      Z1=Z+1D0
      DO 1 L=1,9
      ZL =DCMPLX(DBLE(L),0D0)
      ZL1=Z+ZL+1D0
      PS = PSI(ZL1)
      S1=PS+GE
      T=T-AK1(L)*(ZET2*Z/(Z+ZL)+ZL/(Z+ZL)**2*S1)
1     CONTINUE
      GOTO 100
C
20    T=1.01/(ZN+1D0)-0.846/(ZN+1D0*2)+1.155/(ZN+1D0*3)
     &   -1.074/(ZN+1D0*4)+0.55/(ZN+1D0*5)
      GOTO 100
30    T=1.004/(ZN+1D0)-0.846/(ZN+1D0*2)+1.342/(ZN+1D0*3)
     &   -1.532/(ZN+1D0*4)+0.839/(ZN+1D0*5)
100   CONTINUE
      ACG3=T
C
      RETURN
      END
