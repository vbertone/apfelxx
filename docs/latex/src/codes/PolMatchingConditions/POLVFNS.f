      PROGRAM POLVNS
*     --------------
*     ---------------------------------------------------------------
*
*     J. Bluemlein Jan 03 2024: The two loop two mass VFNS transision
*     matrix emlements in the polarized case and for transversity
*
*     To be quoted at any use:
*     \bibitem{Bierenbaum:2022biv}
*     I.~Bierenbaum, J.~Bl\"umlein, A.~De Freitas, A.~Goedicke, S.~Klein 
*     and K.~Sch\"onwald,
*     {\it O(\ensuremath{\alpha}s2) polarized heavy flavor corrections to 
*     deep-inelastic scattering at Q2???\ensuremath{\gg}???m2},
*     Nucl. Phys. B \textbf{988} (2023), 116114
*     doi:10.1016/j.nuclphysb.2023.116114
*     [arXiv:2211.15337 [hep-ph]].
*
*     For transversity one has to quote:
*     \bibitem{Blumlein:2009rg}
*     J.~Bl\"umlein, S.~Klein and B.~T\"odtli,
*     {\it $O(alpha_s^2)$ and $O(\alpha_s^3)$ Heavy Flavor Contributions to Transversity 
*     at $Q^2 \gg m^2$},
*     Phys. Rev. D \textbf{80} (2009), 094010
*     doi:10.1103/PhysRevD.80.094010
*     [arXiv:0909.1547 [hep-ph]].

*     Furthermore, the code uses routines from:
*     R. Piessens, Angew. Informatik 9 (1973) 399???401. [daind,daind1]
*
*     \bibitem{Blumlein:2000hw}
*      J.~Bl\"umlein,
*      {\it Analytic continuation of Mellin transforms up to two loop order},
*      Comput. Phys. Commun. \textbf{133} (2000), 76-104
*      doi:10.1016/S0010-4655(00)00156-9
*      [arXiv:hep-ph/0003100 [hep-ph]].
*
*     shall also be quoted.
*
*     The main programme allows for a test of its parts.
*     The structure of the VFNS has been given in 
*
*      \bibitem{Blumlein:2018jfm}
*      J.~Bl\"umlein, A.~De Freitas, C.~Schneider and K.~Sch\"onwald,
*      {\it The Variable Flavor Number Scheme at Next-to-Leading Order},
*      Phys. Lett. B \textbf{782} (2018), 362-366
*      doi:10.1016/j.physletb.2018.05.054
*      [arXiv:1804.03129 [hep-ph]].
*
*     The Mellin convolution for the NS and gluonic contributions is split into three
*     parts as described in Blumlein:2018jfm, otherwise only the regular or delta-parts
*     contribute.  
*
*     The quark masses are in the OMS, i.e.:
*     m_c = 1.59 GeV, m_b = 4.78 GeV
*
*     Compilation by gfortran.
*     ---------------------------------------------------------------
*
      IMPLICIT NONE
*
      REAL*8 Q2,XM2,XC2,XB2,CF,CA,TF,Z3,Z2
      REAL*8 NSDEL,NSPLU,NSREG,PSREG,QG1REG,QG2REG,GQREG
      REAL*8 GG1DEL,GG2DEL,GG2PLU,GG2REG,TWOMQG2,TWOMGG2
      REAL*8 X,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13
      REAL*8 T14,T15,T16,TRNSDEL,TRNSPLU,TRNSREG
*
      COMMON/VAR/ Q2
      COMMON/MASS/ XM2,XC2,XB2
      COMMON/COLOR/ CF,CA,TF
      COMMON/CONST/ Z3
      COMMON/CONST1/ Z2
*
      Z3=1.2020569031595942854D0
      Z2=1.6449340668482264365D0
*
      CF=4.0D0/3.0D0
      CA=3.0D0
      TF=1.0D0/2.0D0
      X  =1.0D0/3.0D0
      Q2 = 5.0D1 
      XC2 = (1.59D0)**2
      XB2 = (4.78D0)**2
      XM2 = XC2
      T1=NSDEL(X)   ! ANS2qqH
      T2=NSPLU(X)
      T3=NSREG(X)
      T4=PSREG(X)   ! APS2Hq
      T5=QG1REG(X)  ! AS1Hg
      T6=QG2REG(X)  ! AS2Hg
      T7=GQREG(X)   ! AS2gqH
      T8=GG1DEL(X)  ! AS2ggH
      T9=GG2DEL(X)
      T10=GG2PLU(X)
      T11=GG2REG(X)

C     Two-mass matching
      T12=TWOMQG2(X)
      T13=TWOMGG2(X)

C     Transversity matching
      T14=TRNSDEL(X)
      T15=TRNSPLU(X)
      T16=TRNSREG(X)
*
      WRITE(6,*) 'X,Q2,NSDEL=',X,Q2,T1
      WRITE(6,*) 'X,Q2,NSPLU=',X,Q2,T2
      WRITE(6,*) 'X,Q2,NSREG=',X,Q2,T3
      WRITE(6,*) 'X,Q2,PSREG=',X,Q2,T4
      WRITE(6,*) 'X,Q2,QG1REG=',X,Q2,T5
      WRITE(6,*) 'X,Q2,QG2REG=',X,Q2,T6
      WRITE(6,*) 'X,Q2,GQREG=',X,Q2,T7
      WRITE(6,*) 'X,Q2,GG1DEL=',X,Q2,T8
      WRITE(6,*) 'X,Q2,GG2DEL=',X,Q2,T9
      WRITE(6,*) 'X,Q2,GG2PLU=',X,Q2,T10
      WRITE(6,*) 'X,Q2,GG2REG=',X,Q2,T11
      WRITE(6,*) 'X,Q2,TWOMQG2=',X,Q2,T12
      WRITE(6,*) 'X,Q2,TWOMGG2=',X,Q2,T13
      WRITE(6,*) 'X,Q2,TRNSDEL=',X,Q2,T14
      WRITE(6,*) 'X,Q2,TRNSPLU=',X,Q2,T15
      WRITE(6,*) 'X,Q2,TRNSREG=',X,Q2,T16
*
      STOP
      END
      REAL*8 FUNCTION NSDEL(X)
*     ------------------------
*     ----------------------------------------------------------------
*
*     J. Bluemlein: 02.01.2024 pol VFNS Delta-part NS 2 loop
*     (This is a constant.)
*     ----------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 X,LM,XM2,XC2,XB2,Q2,CF,CA,TF,Z2,Z3
*
      COMMON/VAR/ Q2
      COMMON/MASS/ XM2,XC2,XB2
      COMMON/COLOR/ CF,CA,TF
*
      Z2=1.6449340668482264365D0
      Z3=1.2020569031595942854D0 
      LM=LOG(XM2/Q2)
*
      NSDEL=CF*TF/18.0D0*(73.0D0 + 12.0D0*LM + 36.0D0*LM**2 
     &     + 80.0D0*Z2+ 96.0D0*LM*Z2 - 48.0D0*Z3)
*
      RETURN
      END
      REAL*8 FUNCTION NSPLU(X)
*     ------------------------
*     ----------------------------------------------------------------
*
*     J. Bluemlein: 02.01.2024 pol VFNS plus-part NS 2 loop
*     ----------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 X,LM,XM2,XC2,XB2,Q2,CF,CA,TF
*
      COMMON/VAR/ Q2
      COMMON/MASS/ XM2,XC2,XB2
      COMMON/COLOR/ CF,CA,TF
*
      LM=LOG(XM2/Q2)
      NSPLU=CF*TF*(224.0D0/27.0D0 + 80.0D0*LM/9.0D0 
     &     + 8.0D0*LM**2/3.0D0)/(1.0D0 - X)
*
      RETURN
      END
      REAL*8 FUNCTION NSREG(X)
*     ------------------------
*     ----------------------------------------------------------------
*
*     J. Bluemlein: 02.01.2024 pol VFNS regular part NS 2 loop
*     ----------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 X,LM,XM2,XC2,XB2,Q2,CF,CA,TF
       REAL*8 w(7) 
*
      COMMON/VAR/ Q2
      COMMON/MASS/ XM2,XC2,XB2
      COMMON/COLOR/ CF,CA,TF
*
      LM=LOG(XM2/Q2)
*
      w(1)=Log(x)
      w(2)=1/( - 9.0D0 + 9*x)
      w(3)=1/( - 3.0D0 + 3*x)
      w(4)=x + 1.0D0
      w(4)=1.D0/3.D0*w(4)
      w(5)=w(4) + 2*w(3)
      w(5)=w(5)*w(1)
      w(6)= - 1 + 11*x
      w(7)= - 20*w(2) - 2.D0/9.D0*w(6) - w(5)
      w(7)=w(1)*w(7)
      w(5)= - 1.D0/9.D0*w(6) - w(5)
      w(4)= - LM*w(4)
      w(4)=2*w(5) + w(4)
      w(4)=LM*w(4)
      w(5)=11 - 67*x
      w(4)=2*w(4) + 2.D0/27.D0*w(5) + w(7)
      w(4)=2*w(4)
*
      NSREG = TF*CF*w(4)
*
      RETURN
      END
      REAL*8 FUNCTION TRNSDEL(X)
*     ------------------------
*     ----------------------------------------------------------------
*
*     J. Bluemlein: 02.01.2024 pol VFNS Delta-part Transversity 2 loop
*     (This is a constant.)
*     ----------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 X,LM,XM2,XC2,XB2,Q2,CF,CA,TF,Z2,Z3
*
      COMMON/VAR/ Q2
      COMMON/MASS/ XM2,XC2,XB2
      COMMON/COLOR/ CF,CA,TF
*
      Z2=1.6449340668482264365D0
      Z3=1.2020569031595942854D0 
      LM=LOG(XM2/Q2)
*
      TRNSDEL=CF*TF/18.0D0*(73.0D0 + 12.0D0*LM + 36.0D0*LM**2 
     &     + 80.0D0*Z2+ 96.0D0*LM*Z2 - 48.0D0*Z3)
*
      RETURN
      END
      REAL*8 FUNCTION TRNSPLU(X)
*     --------------------------
*     ----------------------------------------------------------------
*
*     J. Bluemlein: 02.01.2024 pol VFNS plus-parttransversity  NS 2 loop
*     ----------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 X,LM,XM2,XC2,XB2,Q2,CF,CA,TF
*
      COMMON/VAR/ Q2
      COMMON/MASS/ XM2,XC2,XB2
      COMMON/COLOR/ CF,CA,TF
*
      LM=LOG(XM2/Q2)
      TRNSPLU=CF*TF*8.0D0/27.0D0*(28.0D0 + 30.0D0*LM 
     &     + 9.0D0*LM**2)/(1.0D0 - X)
*
      RETURN
      END
      REAL*8 FUNCTION TRNSREG(X)
*     --------------------------
*     ----------------------------------------------------------------
*
*     J. Bluemlein: 02.01.2024 pol VFNS transversirty-NS regular part at 2 loop
*     ----------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 X,LM,XM2,XC2,XB2,Q2,CF,CA,TF
*
      COMMON/VAR/ Q2
      COMMON/MASS/ XM2,XC2,XB2
      COMMON/COLOR/ CF,CA,TF
*
      LM=LOG(XM2/Q2)
      TRNSREG=CF*TF*(
     &-8.0D0*LM**2/3.0D0 - 4.0D0/27.0D0*(47.0D0+9.0D0*X)
     &+40.0D0/9.0D0*X*LOG(X)/(1.0D0-X)+4.0D0/3.0D0*X*LOG(X)**2/(1.0D0-X)
     &+LM*(-80.0D0/9.0D0 + 16.0D0/3.0D0*X*LOG(X)/(1.0D0-X)))
*
      RETURN
      END
      REAL*8 FUNCTION PSREG(X)
*     ------------------------
*     ----------------------------------------------------------------
*
*     J. Bluemlein: 02.01.2024 pol VFNS PS 2 loop
*     ----------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 X,LM,XM2,XC2,XB2,Q2,CF,CA,TF,Z3,FLI2,FLI3
      REAL*8 w(11)
*
      COMMON/VAR/ Q2
      COMMON/MASS/ XM2,XC2,XB2
      COMMON/COLOR/ CF,CA,TF
      COMMON/CONST/ Z3
*
      LM=LOG(XM2/Q2)
*
      w(1)=Log(x)
      w(2)=Log(1.0D0 - x)
      w(3)=FLI2(x)
      w(4)=FLI3(x)
      w(5)=x - 1
      w(6)=x + 1
      w(7)=w(6)*w(1)
      w(8)=w(7) + 1 - 3*x
      w(8)=w(1)*w(8)
      w(8)=w(8) + w(5)
      w(9)=2*w(1)
      w(9)=w(9)*w(6)
      w(10)=5*w(5)
      w(11)=w(10) - w(9)
      w(11)=LM*w(11)
      w(8)=2*w(8) + w(11)
      w(8)=LM*w(8)
      w(11)=w(4) - z3
      w(11)= - 16*w(11)
      w(6)=w(6)*w(11)
      w(11)=w(2)*w(5)
      w(11)=10*w(11) + 3 - 17*x
      w(7)= - 2.D0/3.D0*w(7) + 3 + 5*x
      w(7)=w(1)*w(7)
      w(7)=2*w(11) + w(7)
      w(7)=w(1)*w(7)
      w(9)=w(10) + w(9)
      w(9)=w(3)*w(9)
      w(5)=4*w(9) + 2*w(8) + 28*w(5) + w(7) + w(6)
      w(5)=2*w(5)
*
      PSREG = TF*CF*w(5)
*
      RETURN
      END
      REAL*8 FUNCTION GQREG(X)
*     ------------------------
*     ----------------------------------------------------------------
*
*     J. Bluemlein: 02.01.2024 pol VFNS GQ 2 loop
*     ----------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 X,LM,XM2,XC2,XB2,Q2,CF,CA,TF,Z3,FLI2,FLI3
*
      COMMON/VAR/ Q2
      COMMON/MASS/ XM2,XC2,XB2
      COMMON/COLOR/ CF,CA,TF
      COMMON/CONST/ Z3
*
      LM=LOG(XM2/Q2)
      GQREG= CF*TF*((-8.0D0*LM**2*(-2.0D0 + X))/3.0D0 - (32.0D0*(-11.0D0 
     &     + 4.0D0*X))/27.0D0 + (8.0D0*(4.0D0 + X)*Log(1 - X))/9.0D0 - 
     &     (4.0D0*(-2.0D0 + X)*Log(1.0D0 - X)**2)/3.0D0 + 
     &     LM*((16.0D0*(4.0D0 + X))/9.0D0 - (16.0D0*(-2.0D0 + X)
     &     *Log(1 - X))/3.0D0))
*
      RETURN
      END
      REAL*8 FUNCTION GG1DEL(X)
*     ------------------------
*     ----------------------------------------------------------------
*
*     J. Bluemlein: 02.01.2024 pol VFNS Delta-part gluon 1 loop
*     (This is a constant.)
*     ----------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 X,LM,XM2,XC2,XB2,Q2,CF,CA,TF
*
      COMMON/VAR/ Q2
      COMMON/MASS/ XM2,XC2,XB2
      COMMON/COLOR/ CF,CA,TF
*
      LM=LOG(XM2/Q2)
      GG1DEL = 4.0D0/3.0D0*LM*TF
*
      RETURN
      END
      REAL*8 FUNCTION GG2DEL(X)
*     ------------------------
*     ----------------------------------------------------------------
*
*     J. Bluemlein: 02.01.2024 pol VFNS Delta-part gluon 2 loop
*     (This is a constant.)
*     ----------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 X,LM,XM2,XC2,XB2,Q2,CF,CA,TF
*
      COMMON/VAR/ Q2
      COMMON/MASS/ XM2,XC2,XB2
      COMMON/COLOR/ CF,CA,TF
*
      LM=LOG(XM2/Q2)
      GG2DEL=CF*(-15.0D0 + 4.0D0*LM)*TF 
     &      + 2.0D0/9.0D0*CA*(5.0D0 + 24.0D0*LM)*TF 
     &      + 16.0D0*LM**2*TF**2/9.0D0
*
      RETURN
      END
      REAL*8 FUNCTION GG2PLU(X)
*     ------------------------
*     ----------------------------------------------------------------
*
*     J. Bluemlein: 02.01.2024 pol VFNS plus-part gluon 2 loop
*     ----------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 X,LM,XM2,XC2,XB2,Q2,CF,CA,TF
*
      COMMON/VAR/ Q2
      COMMON/MASS/ XM2,XC2,XB2
      COMMON/COLOR/ CF,CA,TF
*
      LM=LOG(XM2/Q2)
      GG2PLU=(CA*(224.0D0/27.0D0 + (80.0D0*LM)/9.0D0 
     &      + (8.0D0*LM**2)/3.0D0)*TF)/(1 - X)
*
      RETURN
      END
      REAL*8 FUNCTION GG2REG(X)
*     ------------------------
*     ----------------------------------------------------------------
*
*     J. Bluemlein: 02.01.2024 pol VFNS regular part gluon 2 loop
*     ----------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 X,LM,XM2,XC2,XB2,Q2,CF,CA,TF
*
      COMMON/VAR/ Q2
      COMMON/MASS/ XM2,XC2,XB2
      COMMON/COLOR/ CF,CA,TF
*
      LM=LOG(XM2/Q2)
      GG2REG=CA*TF*((-8.0D0*LM**2*(-1.0D0 + 2.0D0*X))/3.0D0 - 
     & (2.0D0*(-337.0D0 + 449.0D0*X))/27.0D0 - (4.0D0*X
     & *Log(1.0D0-X))/3.0D0 + (4.0D0*(22.0D0+X)*Log(X))/9.0D0 
     & +(4.0D0*(1.0D0+X)*Log(X)**2)/3.0D0 + LM*((-16.0D0*(-14.0D0 
     & + 19.0D0*X))/9.0D0 + (16.0D0*(1.0D0 + X)*Log(X))/3.0D0))
     & + CF*TF*(-56.0D0*(-1.0D0 + X) + 12.0D0*(3.0D0 + X)*Log(X) - 
     & 2.0D0*(-5.0D0 + X)*Log(X)**2 + (4.0D0*(1.0D0 + X)
     & *Log(X)**3)/3.0D0 + LM**2*(-20.0D0*(-1.0D0 + X) + 8.0D0
     & *(1.0D0 + X)*Log(X)) +  LM*(-40.0D0*(-1.0D0 + X) 
     & - 8.0D0*(-5.0D0 + X)*Log(X) + 8.0D0*(1.0D0 + X)*Log(X)**2))
*
      RETURN
      END
      REAL*8 FUNCTION QG2REG(X)
*     ------------------------
*     ----------------------------------------------------------------
*
*     J. Bluemlein: 02.01.2024 pol VFNS regular part QG 1 loop
*     ----------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 X,LM,XM2,XC2,XB2,Q2,CF,CA,TF,FLI2,FLI3,S12,Z3,Z2
      REAL*8 TA
      REAL*8 w
      DIMENSION w(31)
*
      COMMON/VAR/ Q2
      COMMON/MASS/ XM2,XC2,XB2
      COMMON/COLOR/ CF,CA,TF
      COMMON/CONST/ Z3
      COMMON/CONST1/ Z2
*
      LM=LOG(XM2/Q2)
*
      w(1)=Log(1.D0 - X)
      w(2)=FLI2(1.D0 - X)
      w(3)=Log(1.D0 - X)
      w(4)=FLI2(1.D0 - X)
      w(5)=Log(1.D0 - X)
      w(6)=Log(1.D0 - X)
      w(7)=Log(1.D0 + X)
      w(8)=FLI2( - X)
      w(9)=Log(X)
      w(10)=FLI2(X)
      w(11)=S12( - X)
      w(12)=FLI3(1.D0 - X)
      w(13)=FLI3( - X)
      w(14)=FLI3(X)
      w(15)=2.D0*X
      w(16)=w(15) + 1.D0
      w(17)=2.D0*w(7)
      w(18)=w(16)*w(17)
      w(19)=X + 4.D0
      w(20)=w(19)*X
      w(21)=4.D0*LM
      w(22)= - w(18) + w(21) + 4.D0 - w(20)
      w(22)=Z2*w(22)
      w(23)=w(11)*w(16)
      w(24)=3.D0 + 5.D0*X
      w(24)=z3*w(24)
      w(23)=w(23) - w(24)
      w(24)=14.D0 + X
      w(24)=w(24)*w(15)
      w(24)= - 25.D0 + w(24)
      w(24)=w(10)*w(24)
      w(25)=2.D0*LM
      w(26)=X - 1.D0
      w(27)=LM*w(26)
      w(27)=6.D0*w(27) - 12.D0 + 11.D0*X
      w(27)=w(27)*w(25)
      w(28)=w(15) - 1.D0
      w(29)=w(12)*w(28)
      w(30)=X + 1.D0
      w(31)=w(14)*w(30)
      w(22)=w(22) + w(27) - 16.D0*w(31) + w(24) - 2.D0*w(29) - 51.D0 + 
     & 53.D0*X - 4.D0*w(23)
      w(23)=8.D0*w(10)
      w(24)=w(30)*w(23)
      w(27)= - w(30)*w(25)
      w(27)=w(27) - 1.D0 - 8.D0*X
      w(27)=w(27)*w(25)
      w(29)=24.D0 + X
      w(29)=X*w(29)
      w(29)= - 25.D0 + w(29)
      w(29)=w(6)*w(29)
      w(24)=w(27) + w(29) + w(24) - 14.D0 - 37.D0*X
      w(27)=w(16)*LM
      w(27)=w(27) + w(30)
      w(29)=w(9)*w(16)
      w(18)= - w(18) + 4.D0*w(27) + w(29)
      w(18)=w(17)*w(18)
      w(29)= - w(21) + 2.D0/3.D0*w(9)
      w(30)= - w(16)*w(29)
      w(31)=X**2
      w(30)= - 3.D0 - w(31) + w(30)
      w(30)=w(9)*w(30)
      w(18)=w(18) + 2.D0*w(24) + w(30)
      w(18)=w(9)*w(18)
      w(24)=w(28)*LM
      w(24)=w(24) + 4.D0*w(26)
      w(24)=w(24)*w(25)
      w(26)=1.D0 - w(24)
      w(21)= - 2.D0/3.D0*w(1) + w(21)
      w(21)=w(28)*w(21)
      w(20)= - 5.D0 + w(20) + w(21)
      w(20)=w(1)*w(20)
      w(20)= - 4.D0*w(2) + 2.D0*w(26) + w(20)
      w(20)=w(1)*w(20)
      w(21)=w(5)*w(4)*X
      w(18)=8.D0*w(21) + w(20) + 2.D0*w(22) + w(18)
      w(17)= - w(17) + w(9)
      w(17)=w(16)*w(17)
      w(17)=2.D0*w(27) + w(17)
      w(17)=w(8)*w(17)
      w(16)=w(13)*w(16)
      w(16)=w(17) - w(16)
      w(16)=2.D0*w(18) + 8.D0*w(16)
      w(16)=CA*w(16)
      w(17)= - w(28)*w(23)
      w(15)= - w(19)*w(15)
      w(18)= - w(6)*w(28)
      w(15)=w(18) + 11.D0 + w(15)
      w(15)=w(6)*w(15)
      w(18)=4.D0*w(6)
      w(19)= - LM + w(18)
      w(19)=w(28)*w(19)
      w(20)=16.D0*X
      w(19)=1.D0 - w(20) + w(19)
      w(19)=w(19)*w(25)
      w(15)=w(19) + 2.D0*w(15) + w(17) - 8.D0 - 25.D0*X
      w(17)= - w(18) + w(29)
      w(17)=w(28)*w(17)
      w(18)=X + 2.D0
      w(18)=w(18)*X
      w(17)= - 5.D0 + 4.D0*w(18) + w(17)
      w(17)=w(9)*w(17)
      w(15)=2.D0*w(15) + w(17)
      w(15)=w(9)*w(15)
      w(17)=1.D0/3.D0*w(1) - w(25)
      w(17)=w(28)*w(17)
      w(17)=3.D0 - w(18) + w(17)
      w(17)=w(1)*w(17)
      w(17)=w(17) + w(24) - 16.D0 + 15.D0*X
      w(17)=w(1)*w(17)
      w(19)=4.D0*X
      w(21)= - 3.D0 - X
      w(19)=w(21)*w(19)
      w(19)=29.D0 + w(19)
      w(19)=w(10)*w(19)
      w(17)=w(19) + w(17)
      w(19)= - 8.D0*z3 + 12.D0*w(14) + 4.D0*w(12)
      w(19)=w(28)*w(19)
      w(21)= - 2.D0 - 3.D0*X
      w(21)=2.D0*w(21) + 3.D0*LM
      w(21)=LM*w(21)
      w(17)=w(21) - 11.D0 + 23.D0*X + w(19) + 2.D0*w(17)
      w(19)=w(28)*w(25)
      w(18)=w(19) - 9.D0 + w(18)
      w(18)=Z2*w(18)
      w(19)= - w(6)*w(20)
      w(19)=8.D0*w(3) + w(19)
      w(19)=w(19)*w(4)
      w(15)=w(19) + 8.D0*w(18) + 2.D0*w(17) + w(15)
      w(15)=CF*w(15)
      w(17)=TF*w(28)*LM**2
      w(15)= - 16.D0/3.D0*w(17) + w(16) + w(15)

      QG2REG = TF*w(15)
*
      RETURN
      END
      REAL*8 FUNCTION QG1REG(X)
*     -------------------------
*     ----------------------------------------------------------------
*
*     J. Bluemlein: 02.01.2024 pol VFNS regular part QG 1 loop
*     ----------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 X,LM,XM2,XC2,XB2,Q2,CF,CA,TF
*
      COMMON/VAR/ Q2
      COMMON/MASS/ XM2,XC2,XB2
      COMMON/COLOR/ CF,CA,TF
*
      LM=LOG(XM2/Q2)
      QG1REG=-4.0D0*LM*TF*(-1.0D0 + 2.0D0*X)
*
      RETURN
      END
      REAL*8 FUNCTION TWOMQG2(X)
*     ---------------------------
*     ----------------------------------------------------------------*
*     J. Bluemlein: 02.01.2024 pol VFNS 2-mass QG 2 loop
*     ----------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 X,LC,LB,XC2,XB2,Q2,CF,CA,TF,B0Q,XM2
*
      COMMON/VAR/ Q2
      COMMON/MASS/ XM2,XC2,XB2
      COMMON/COLOR/ CF,CA,TF
*
      B0Q=-4.0D0/3.0D0*TF
      LC=LOG(XC2/Q2)
      LB=LOG(XB2/Q2)
      TWOMQG2=-B0Q*LC*LB*8.0D0*TF*(1.0D0-2.0D0*X)
*
      RETURN
      END
      REAL*8 FUNCTION TWOMGG2(X)
*     --------------------------
*     ----------------------------------------------------------------
*
*     J. Bluemlein: 02.01.2024 pol VFNS 2-mass GG 2 loop
*     ----------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 X,LC,LB,XC2,XB2,Q2,CF,CA,TF,B0Q,XM2
*
      COMMON/VAR/ Q2
      COMMON/MASS/ XM2,XC2,XB2
      COMMON/COLOR/ CF,CA,TF
*
      B0Q=-4.0D0/3.0D0*TF
      LC=LOG(XC2/Q2)
      LB=LOG(XB2/Q2)
      TWOMGG2=2.0D0*B0Q**2*LC*LB
*
      RETURN
      END
      DOUBLE PRECISION FUNCTION FLI2(X)
C     ---------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  LI2(X)
C
C************************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      EXTERNAL FR1
C
      PI = 3.141592653589793238462643D0
      ZETA2 = PI**2/6.0D0
C
      IF(X.LT.-1.0D0.OR.X.GT.1.0D0) GOTO 100
      IF(X.EQ.1.0D0)  GOTO 5
      IF(X.EQ.-1.0D0) GOTO 6
      IF(X.LT.0.0D0)  GOTO 1
      IF(X.EQ.0.0D0)  GOTO 2
      IF(X.GT.0.5D0)  GOTO 3
      T=FR1(X)
      GOTO 200
100   WRITE(6,*) 'FLI2 -> NOT ALLOWED,X=',X,'STOP ***'
      STOP
1     Y=-X
      IF(Y.GT.0.5D0) GOTO 4
      Y2=Y**2
      T=FR1(Y2)/2.0D0-FR1(Y)
      GOTO 200
2     T=0.0D0
      GOTO 200
3     XM=1.0D0-X
      T=-FR1(XM)-LOG(X)*LOG(XM)+ZETA2
      GOTO 200
4     YM=1.0D0-Y
      T1=-FR1(YM)-LOG(Y)*LOG(YM)+ZETA2
      Y2=Y**2
      IF(Y2.LT.0.5) THEN
      T2=FR1(Y2)
      ELSE
      T2=-FR1(1.0D0-Y2)-LOG(Y2)*LOG(1.0D0-Y2)+ZETA2
      ENDIF
      T=T2/2.0D0-T1
      GOTO 200
5     T=ZETA2
      GOTO 200
6     T=-ZETA2/2.0D0
      GOTO 200
C
200   FLI2=T
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FR1(X)
C     --------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  SUBSIDIARTY ROUTINE FOR  LI2(X)
C
C************************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      Y=LOG(1.0D0-X)
      Y2=Y*Y
C
      T = (-1.D0+(-1.D0/4.D0+(-1.D0/36.D0+(1.D0/3600.D0+(
     &-1.D0/211680.D0+(1.D0/10886400.D0+(-1.D0/526901760.D0
     &+(691.D0/16999766784000.D0+(-1.D0/1120863744000.D0
     &+3617.D0/0.18140058832896D18*Y2)*Y2)*Y2)*Y2)*Y2)*Y2)
     &*Y2)*Y)*Y)*Y
C
      FR1=T
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FLI3(X)
C     -----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  LI3(X) FOR -1. LE . X . LE .+1
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      A=1D0
      F=0D0
      AN=0D0
      TCH=1D-16
1     AN=AN+1D0
      A=A*X
      B=A/AN**3
      F=F+B
      IF((ABS(B)-TCH).LE.0.0D0) GOTO2
      IF((ABS(B)-TCH).GT.0.0D0) GOTO1
2     FLI3=F
      END
      DOUBLE PRECISION FUNCTION S12(Z)
C     --------------------------------
C
C     S_12 IN INTEGRAL REPRESENTATION
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      EXTERNAL FINT
C
      EPS = 1.0D-8
      KEY = 2
      MAX = 10000
C
      F = DAIND1(0.0D0,Z,FINT,EPS,KEY,MAX,KOU,EST)
      S12 = F/2.0D0
      RETURN
      END
      DOUBLE PRECISION FUNCTION FINT(Y)
C     ---------------------------------
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      FINT = LOG(1.0D0-Y)**2/Y
      RETURN
      END
      DOUBLE PRECISION FUNCTION DAIND1(A,B,FUN,EPS,KEY,MAX,KOUNT,EST)
C     ---------------------------------------------------------------
C  INPUTPARAMETERS
C  A,B      LIMITS OF THE INTEGRATION INTERVAL
C  FUN      FUNCTION TO BE INTEGRATED (TO BE DECLARED EXTERNAL IN THE MAIN PR.)
C  EPS      ABSOLUTE OR RELATIVE TOLERANCE,DEPENDING OF THE VALUE OF 'KEY'
C  KEY      =1 THEN 'EPS' DENOTES AN ABSOLUTE, =2 THEN A RELATIVE TOLERANCE
C  MAX      UPPER BOUND ON THE NUMBERS OF INTEGRAND EVALUATIONS (MAX.LE.10000)
C
C  OUTPUTPARAMETERS
C  KOUNT    NUMBER OF INTEGRAND EVALUATIONS
C  EST      ESTIMATION OF THE ABSOLUTE ERROR OF THE APPROXIMATION
      IMPLICIT REAL*8 (A-H,O-Z)
      DOUBLE PRECISION MAXIM,MINIM,MODUL1,MODUL2
      INTEGER RANG(130)
      DIMENSION
     &AINIT(250),END(250),EPSIL(250),PART(250),W1(5),W2(5),W3(6),
     &     X1(5),X2(5)
      DATA X1/0.973906528517D+0,0.865063366689D+0,0.679409568299D+0,
     *         0.433395394129D+0,0.148874338981D+0/
      DATA X2/0.995657163026D+0,0.930157491356D+0,0.780817726586D+0,
     *         0.562757134669D+0,0.294392862701D+0/
      DATA W1/0.666713443087D-1,0.149451349151D+0,0.219086362516D+0,
     *         0.269266719310D+0,0.295524224715D+0/
      DATA W2/0.325581623080D-1,0.750396748109D-1,0.109387158802D+0,
     *         0.134709217311D+0,0.147739104901D+0/
      DATA W3/0.116946388674D-1,0.547558965744D-1,0.931254545837D-1,
     *         0.123491976262D+0,0.142775938577D+0,0.149445554003D+0/
      DATA TOL/0.23D-15/
      MAX1 = (MAX+21)/42+1
      MAX2 = MAX1/2+2
      ALFA = A
      BETA = B
      MAAT = 1
C EVALUATION OF GAUSSIAN AND KRONROD FORMULAS
   10 S = 0.5D+0*(BETA-ALFA)
      U = 0.5D+0*(BETA+ALFA)
      RES1 = 0.0D+0
      RES2 = W3(6)*FUN(U)
      DO 20 K = 1,5
        C = S*X1(K)
        C = FUN(C+U)+FUN(U-C)
        RES1 = RES1+W1(K)*C
        RES2 = RES2+W2(K)*C
        C = S*X2(K)
        RES2 =RES2+W3(K)*(FUN(C+U)+FUN(U-C))
20    CONTINUE      
      PAT = RES2*S
      MODUL2 = ABS(PAT-RES1*S)
      IF(MAAT.GT.1) GOTO 50
      EST = MODUL2
      BINT = PAT
      KOUNT =21
      PART(1) = BINT
      GOTO 90
   30 RANG(1) = 1
      AINIT(1) = A
      END(1) = B
      EPSIL(1) = EST
   40 NR = RANG(1)
      BINT = BINT-PART(NR)
      EST =EST-EPSIL(NR)
C THE SUBINTERVAL WITH LARGEST ERROR IS SPLIT UP INTO TWO EQUAL PARTS
      ALFA = AINIT(NR)
      BETA = (AINIT(NR)+END(NR))*0.5
      JJ = 1
      MAAT = MAAT+1
      GOTO 10
   50 EST = EST+MODUL2
      BINT = BINT+PAT
      IF(JJ.EQ.0) GOTO 60
      MODUL1 = MODUL2
      PAT1 = PAT
      ALFA = BETA
      BETA = END(NR)
      JJ = 0
      GOTO 10
   60 MA = MAAT
      IF(MAAT.GT.MAX2) MA = MAX1+3-MAAT
      IF(MODUL1.GT.MODUL2) GOTO 70
      EPSIL(NR) = MODUL2
      EPSIL(MAAT) = MODUL1
      AINIT(MAAT) = AINIT(NR)
      AINIT(NR) = ALFA
      END(MAAT) = ALFA
      MAXIM = MODUL2
      MINIM = MODUL1
      PART(NR) = PAT
      PART(MAAT) = PAT1
      GOTO 80
   70 EPSIL(NR) = MODUL1
      EPSIL(MAAT) = MODUL2
      END(MAAT) = BETA
      END(NR) = ALFA
      AINIT(MAAT) = ALFA
      MAXIM = MODUL1
      MINIM = MODUL2
      PART(NR) = PAT1
      PART(MAAT) = PAT
   80 KOUNT = KOUNT+42
C TEST ON THE NUMBER OF FUNCTION EVALUATIONS
      IF(KOUNT.GE.MAX) GOTO 190
   90 GOTO (100,110),KEY
C TEST ON ABSOLUTE ACCURACY
  100 IF(EST.LE.EPS) GOTO 190
      GOTO 120
C TEST ON RELATIVE ACCURACY
  110 IF(ABS(EPS*BINT).LE.TOL) GOTO 100
      IF(EST.LE.ABS(EPS*BINT)) GOTO 190
  120 IF(MAAT.EQ.1) GOTO 30
      IF(MAAT.GT.2) GOTO 130
      RANG(2) = 2
      GOTO 40
  130 MB = MA-1
C SEARCH FOR THE SUBINTERVAL WITH LARGEST ERROR
      DO 140 I = 2,MB
        IR = RANG(I)
        IF(MAXIM.GE.EPSIL(IR)) GOTO 150
        RANG(I-1) = RANG(I)
140     CONTINUE    
      RANG(MB) = NR
      RANG(MA) = MAAT
      GOTO 40
  150 RANG(I-1) = NR
      DO 160 K = I,MB
        IR = RANG(K)
        IF(MINIM.GE.EPSIL(IR)) GOTO 170
  160   CONTINUE
      RANG(MA) = MAAT
      GOTO 40
  170 DO 180 I = K,MB
        KK = MB-I+K
       RANG(KK+1) = RANG(KK)
180   CONTINUE      
      RANG(K) = MAAT
      GOTO 40
C CALCULATION OF THE INTEGRAL
  190 AIND1 = 0.0D+0
      DO 200 K = 1,MAAT
      AIND1 = AIND1+PART(K)
200   CONTINUE       
      IF(AIND1.EQ.0.0D+0)
     &        WRITE(6,*) '**** AIND=0.**** EST NOT CALCULATED'
      IF(AIND1.NE.0.0D+0) EST=EST/AIND1
      DAIND1=AIND1
      RETURN
      END
      DOUBLE PRECISION FUNCTION DAIND(A,B,FUN,EPS,KEY,MAX,KOUNT,EST)
C     --------------------------------------------------------------
C  INPUTPARAMETERS
C  A,B      LIMITS OF THE INTEGRATION INTERVAL
C  FUN      FUNCTION TO BE INTEGRATED (TO BE DECLARED EXTERNAL IN THE MAIN PR.)
C  EPS      ABSOLUTE OR RELATIVE TOLERANCE,DEPENDING OF THE VALUE OF 'KEY'
C  KEY      =1 THEN 'EPS' DENOTES AN ABSOLUTE, =2 THEN A RELATIVE TOLERANCE
C  MAX      UPPER BOUND ON THE NUMBERS OF INTEGRAND EVALUATIONS (MAX.LE.10000)
C
C  OUTPUTPARAMETERS
C  KOUNT    NUMBER OF INTEGRAND EVALUATIONS
C  EST      ESTIMATION OF THE ABSOLUTE ERROR OF THE APPROXIMATION
      IMPLICIT REAL*8 (A-H,O-Z)
      DOUBLE PRECISION MAXIM,MINIM,MODUL1,MODUL2
      INTEGER RANG(130)
      DIMENSION
     &AINIT(250),END(250),EPSIL(250),PART(250),W1(5),W2(5),W3(6),
     &     X1(5),X2(5)
      DATA X1/0.973906528517D+0,0.865063366689D+0,0.679409568299D+0,
     *         0.433395394129D+0,0.148874338981D+0/
      DATA X2/0.995657163026D+0,0.930157491356D+0,0.780817726586D+0,
     *         0.562757134669D+0,0.294392862701D+0/
      DATA W1/0.666713443087D-1,0.149451349151D+0,0.219086362516D+0,
     *         0.269266719310D+0,0.295524224715D+0/
      DATA W2/0.325581623080D-1,0.750396748109D-1,0.109387158802D+0,
     *         0.134709217311D+0,0.147739104901D+0/
      DATA W3/0.116946388674D-1,0.547558965744D-1,0.931254545837D-1,
     *         0.123491976262D+0,0.142775938577D+0,0.149445554003D+0/
      DATA TOL/0.23D-15/
      MAX1 = (MAX+21)/42+1
      MAX2 = MAX1/2+2
      ALFA = A
      BETA = B
      MAAT = 1
C EVALUATION OF GAUSSIAN AND KRONROD FORMULAS
   10 S = 0.5D+0*(BETA-ALFA)
      U = 0.5D+0*(BETA+ALFA)
      RES1 = 0.0D+0
      RES2 = W3(6)*FUN(U)
      DO 20 K = 1,5
        C = S*X1(K)
        C = FUN(C+U)+FUN(U-C)
        RES1 = RES1+W1(K)*C
        RES2 = RES2+W2(K)*C
        C = S*X2(K)
        RES2 =RES2+W3(K)*(FUN(C+U)+FUN(U-C))
20    CONTINUE      
      PAT = RES2*S
      MODUL2 = ABS(PAT-RES1*S)
      IF(MAAT.GT.1) GOTO 50
      EST = MODUL2
      BINT = PAT
      KOUNT =21
      PART(1) = BINT
      GOTO 90
   30 RANG(1) = 1
      AINIT(1) = A
      END(1) = B
      EPSIL(1) = EST
   40 NR = RANG(1)
      BINT = BINT-PART(NR)
      EST =EST-EPSIL(NR)
C THE SUBINTERVAL WITH LARGEST ERROR IS SPLIT UP INTO TWO EQUAL PARTS
      ALFA = AINIT(NR)
      BETA = (AINIT(NR)+END(NR))*0.5
      JJ = 1
      MAAT = MAAT+1
      GOTO 10
   50 EST = EST+MODUL2
      BINT = BINT+PAT
      IF(JJ.EQ.0) GOTO 60
      MODUL1 = MODUL2
      PAT1 = PAT
      ALFA = BETA
      BETA = END(NR)
      JJ = 0
      GOTO 10
   60 MA = MAAT
      IF(MAAT.GT.MAX2) MA = MAX1+3-MAAT
      IF(MODUL1.GT.MODUL2) GOTO 70
      EPSIL(NR) = MODUL2
      EPSIL(MAAT) = MODUL1
      AINIT(MAAT) = AINIT(NR)
      AINIT(NR) = ALFA
      END(MAAT) = ALFA
      MAXIM = MODUL2
      MINIM = MODUL1
      PART(NR) = PAT
      PART(MAAT) = PAT1
      GOTO 80
   70 EPSIL(NR) = MODUL1
      EPSIL(MAAT) = MODUL2
      END(MAAT) = BETA
      END(NR) = ALFA
      AINIT(MAAT) = ALFA
      MAXIM = MODUL1
      MINIM = MODUL2
      PART(NR) = PAT1
      PART(MAAT) = PAT
   80 KOUNT = KOUNT+42
C TEST ON THE NUMBER OF FUNCTION EVALUATIONS
      IF(KOUNT.GE.MAX) GOTO 190
   90 GOTO (100,110),KEY
C TEST ON ABSOLUTE ACCURACY
  100 IF(EST.LE.EPS) GOTO 190
      GOTO 120
C TEST ON RELATIVE ACCURACY
  110 IF(ABS(EPS*BINT).LE.TOL) GOTO 100
      IF(EST.LE.ABS(EPS*BINT)) GOTO 190
  120 IF(MAAT.EQ.1) GOTO 30
      IF(MAAT.GT.2) GOTO 130
      RANG(2) = 2
      GOTO 40
  130 MB = MA-1
C SEARCH FOR THE SUBINTERVAL WITH LARGEST ERROR
      DO 140 I = 2,MB
        IR = RANG(I)
        IF(MAXIM.GE.EPSIL(IR)) GOTO 150
        RANG(I-1) = RANG(I)
140     CONTINUE    
      RANG(MB) = NR
      RANG(MA) = MAAT
      GOTO 40
  150 RANG(I-1) = NR
      DO 160 K = I,MB
        IR = RANG(K)
        IF(MINIM.GE.EPSIL(IR)) GOTO 170
  160   CONTINUE
      RANG(MA) = MAAT
      GOTO 40
  170 DO 180 I = K,MB
        KK = MB-I+K
       RANG(KK+1) = RANG(KK)
180   CONTINUE      
      RANG(K) = MAAT
      GOTO 40
C CALCULATION OF THE INTEGRAL
  190 AIND1 = 0.0D+0
      DO 200 K = 1,MAAT
      AIND1 = AIND1+PART(K)
200   CONTINUE       
      IF(AIND1.EQ.0.0D+0)
     &        WRITE(6,*) '**** AIND=0.**** EST NOT CALCULATED'
      IF(AIND1.NE.0.0D+0) EST=EST/AIND1
      DAIND1=AIND1
      RETURN
      END

