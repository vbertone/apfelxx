Replacements = {"pow(NC,2)"->"NC2", "pow(NC,-1)"->"NCi", "pow(NC,-2)"->"NCi2", "log(2)"->"l2", "pow(log(2),2)"->"l22", "pow(log(2),3)"->"l23", "pow(pi,2)"->"pi2", "pow(pi,4)"->"pi4", "pow(x,2)"->"x2", "pow(x,3)"->"x3", "pow(x,4)"->"x4", "pow(x,5)"->"x5", "pow(x,6)"->"x6", "pow(x,7)"->"x7", "pow(x,-1)"->"xi", "pow(x,-2)"->"xi2", "pow(z,2)"->"z2", "pow(z,3)"->"z3", "pow(z,-1)"->"zi", "pow(z,-2)"->"zi2", "pow(1 - x,-1)"->"omxi", "pow(1 + x,-1)"->"opxi", "pow(1 - z,-1)"->"omzi", "pow(1 + z,-1)"->"opzi", "pow(-1 + z,-1)"->"mopzi", "log(x)"->"lx", "pow(log(x),2)"->"lx2", "pow(log(x),3)"->"lx3", "log(z)"->"lz", "pow(log(z),2)"->"lz2", "pow(log(z),3)"->"lz3", "log(1 - x)"->"lomx", "pow(log(1 - x),2)"->"lomx2", "pow(log(1 - x),3)"->"lomx3", "log(1 - z)"->"lomz", "pow(log(1 - z),2)"->"lomz2", "pow(log(1 - z),3)"->"lomz3", "log(1 + x)"->"lopx", "pow(log(1 + x),2)"->"lopx2", "pow(log(1 + x),3)"->"lopx3", "log(1 + z)"->"lopz", "pow(log(1 + z),2)"->"lopz2", "dilog(x)"->"Li2x", "dilog(-x)"->"Li2mx", "dilog(z)"->"Li2z", "dilog(-z)"->"Li2mz", "trilog(x)"->"Li3x", "trilog(1 - x)"->"Li3omx", "trilog(0.5 - x/2.)"->"Li3omxh", "trilog((1 + x)/2.)"->"Li3opxh", "trilog(-((-1 + x)*opxi))"->"Li3omxopxi", "trilog(z)"->"Li3z", "trilog(-z)"->"Li3mz", "trilog(1 - z)"->"Li3omz", "trilog(0.5 - z/2.)"->"Li3omzh", "trilog((1 + z)/2.)"->"Li3opzh", "trilog(opzi)"->"Li3opzi", "trilog(-((-1 + z)*opzi))"->"Li3omzopzi", "pow(1 - x - z,-1)"->"omxmzi", "pow(1 - x - z,-2)"->"omxmzi2", "pow(x - z,-1)"->"xmzi", "pow(x - z,-2)"->"xmzi2", "log(1 - x - z)"->"lomxmz", "log(-1 + x + z)"->"lmopxpz", "log(x - z)"->"lxmz", "log(-x + z)"->"lmxpz", "log(x + z)"->"lxpz", "log(1 + x*z)"->"lopxz", "log(1 + x*zi)"->"lopxzi", "dilog(1 - x*zi)"->"li2omxzi", "pow(1 + sqrtxz1 - z,-1)"->"spec1i", "log(1 + sqrtxz1 - z)"->"lspec1", "pow(log(1 + sqrtxz1 - z),2)"->"lspec1_2", "log(1 + sqrtxz1 + z)"->"lspec2", "log(sqrtxz3)"->"lspec3", "log(sqrtxz3*z)"->"lspec4", "log(1 - sqrtxz2 + x)"->"lspec5", "log(1 + sqrtxz2 + x)"->"lspec6", "log(1 - 2*z + 4*x*z + z2)"->"lspec7", "pow(log(1 - 2*z + 4*x*z + z2),2)"->"lspec7_2", "dilog(0.5 - sqrtxz1/2. - z/2.)"->"li2spec1", "dilog(0.5 - sqrtxz1/2. + z/2.)"->"li2spec2", "dilog(0.5 - zi/2. - (sqrtxz1*zi)/2.)"->"li2spec3", "dilog(0.5 + zi/2. - (sqrtxz1*zi)/2.)"->"li2spec4", "dilog(0.5 - sqrtxz2/2. - x/2.)"->"li2spec5", "dilog(0.5 + sqrtxz2/2. - x/2.)"->"li2spec6", "dilog(0.5 - xi/2. - (sqrtxz2*xi)/2.)"->"li2spec7", "dilog(0.5 - xi/2. + (sqrtxz2*xi)/2.)"->"li2spec8", "dilog((1 - z)*omxi)"->"li2spec9", "dilog(x*(1 - z)*omxi*zi)"->"li2spec10", "dilog((1 - x)*omzi)"->"li2spec11", "dilog((1 - x)*xiz*omzi)"->"li2spec12", "dilog(z*omxi)"->"li2spec13", "dilog((1 - z)*z*omxi*xi)"->"li2spec14", "dilog(xz*omxi*omzi)"->"li2spec15", "dilog((1 - x)*zi)"->"li2spec16", "dilog((1 - x)*(1 - z)*xizi)"->"li2spec17", "dilog((1 - x)*x*omzi*zi)"->"li2spec18", "dilog(-(x*z))"->"li2spec19", "dilog(-(x*zi))"->"li2spec20", "dilog((1 - z)*sqrtxz1i)"->"li2spec21", "dilog(-spec1i + sqrtxz1*spec1i + z*spec1i)"->"li2spec22", "dilog(sqrtxz1*mopzi)"->"li2spec23", "dilog((1 - x)*z*xi*omzi)"->"li2spec24","dilog(x*z*omxi*omzi)"->"li2spec25","dilog((1 - x)*(1 - z)*xi*zi)"->"li2spec26","dilog((1 - z)*pow(sqrtxz1,-1))"->"li2spec27","atan(sqrtxz3)"->"atanspec1", "atan(sqrtxz3*z)"->"atanspec2", "InvTanInt(-sqrtxz3)"->"itani1", "InvTanInt(sqrtxz3)"->"itani2", "InvTanInt(sqrtxz3*z)"->"itani3", "T(r1)"->"Tr1", "T(r2)"->"Tr2", "T(t1)"->"Tt1", "T(t2)"->"Tt2", "T(u1)"->"Tu1", "T(u2)"->"Tu2", "T(u3)"->"Tu3", "T(u4)"->"Tu4","pow(sqrtxz1,-1)"->"sqrtxz1i","pow(poly2,-1)"->"poly2i","pow(poly2,-2)"->"poly2i2","pow(sqrtxz2,-1)"->"sqrtxz2i"};

/////////////////////////////////////////////////////////////////////////////////
Declarations = {{"NC2", "    const double NC2 = NC * NC;"}, {"NCi", "    const double NCi = 1. / NC;"}, {"NCi2", "    const double NCi2 = NCi * NCi;"}, {"l2", "    const double l2 = log(2);"}, {"l22", "    const double l22 = l2 * l2;"}, {"l23", "    const double l23 = l2 * l22;"}, {"pi2", "    const double pi2 = M_PI * M_PI;"}, {"pi4", "    const double pi4 = pi2 * pi2;"}, {"x2", "    const double x2 = x * x;"}, {"x3", "    const double x3 = x * x2;"}, {"x4", "    const double x4 = x * x3;"}, {"x5", "    const double x5 = x * x4;"}, {"x6", "    const double x6 = x * x5;"}, {"x7", "    const double x7 = x * x6;"}, {"xi", "    const double xi = 1. / x;"}, {"xi2", "    const double xi2 = xi * xi;"}, {"z2", "    const double z2 = z * z;"}, {"z3", "    const double z3 = z * z2;"}, {"zi", "    const double zi = 1. / z;"}, {"zi2", "    const double zi2 = zi * zi;"}, {"omxi", "    const double omxi = 1. / ( 1 - x );"}, {"opxi", "    const double opxi = 1. / ( 1 + x );"}, {"omzi", "    const double omzi = 1. / ( 1 - z );"}, {"opzi", "    const double opzi = 1. / ( 1 + z );"}, {"mopzi", "    const double mopzi = 1. / ( - 1 + z );"}, {"lx", "    const double lx = log(x);"}, {"lx2", "    const double lx2 = lx * lx;"}, {"lx3", "    const double lx3 = lx * lx2;"}, {"lz", "    const double lz = log(z);"}, {"lz2", "    const double lz2 = lz * lz;"}, {"lz3", "    const double lz3 = lz * lz2;"}, {"lomx", "    const double lomx = log(1 - x);"}, {"lomx2", "    const double lomx2 = lomx * lomx;"}, {"lomx3", "    const double lomx3 = lomx * lomx2;"}, {"lomz", "    const double lomz = log(1 - z);"}, {"lomz2", "    const double lomz2 = lomz * lomz;"}, {"lomz3", "    const double lomz3 = lomz * lomz2;"}, {"lopx", "    const double lopx = log(1 + x);"}, {"lopx2", "    const double lopx2 = lopx * lopx;"}, {"lopx3", "    const double lopx3 = lopx * lopx2;"}, {"lopz", "    const double lopz = log(1 + z);"}, {"lopz2", "    const double lopz2 = lopz * lopz;"}, {"Li2x", "    const double Li2x = dilog(x);"}, {"Li2mx", "    const double Li2mx = dilog(-x);"}, {"Li2z", "    const double Li2z = dilog(z);"}, {"Li2mz", "    const double Li2mz = dilog(-z);"}, {"Li3x", "    const double Li3x = trilog(x);"}, {"Li3omx", "    const double Li3omx = trilog(1 - x);"}, {"Li3omxh", "    const double Li3omxh = trilog(0.5 - x/2.);"}, {"Li3opxh", "    const double Li3opxh = trilog((1 + x)/2.);"}, {"Li3omxopxi", "    const double Li3omxopxi = trilog(-((-1 + x)*opxi));"}, {"Li3z", "    const double Li3z = trilog(z);"}, {"Li3mz", "    const double Li3mz = trilog(-z);"}, {"Li3omz", "    const double Li3omz = trilog(1 - z);"}, {"Li3omzh", "    const double Li3omzh = trilog(0.5 - z/2.);"}, {"Li3opzh", "    const double Li3opzh = trilog((1 + z)/2.);"}, {"Li3opzi", "    const double Li3opzi = trilog(opzi);"}, {"Li3omzopzi", "    const double Li3omzopzi = trilog(-((-1 + z)*opzi));"}, {"omxmzi", "    const double omxmzi = 1. / ( 1 - x - z );"}, {"omxmzi2", "    const double omxmzi2 = omxmzi * omxmzi;"}, {"xmzi", "    const double xmzi = 1. / ( x - z );"}, {"xmzi2", "    const double xmzi2 = xmzi * xmzi;"}, {"lomxmz", "    const double lomxmz = log(1 - x - z);"}, {"lmopxpz", "    const double lmopxpz = log(-1 + x + z);"}, {"lxmz", "    const double lxmz = log(x - z);"}, {"lmxpz", "    const double lmxpz = log(-x + z);"}, {"lxpz", "    const double lxpz = log(x + z);"}, {"lopxz", "    const double lopxz = log(1 + x*z);"}, {"lopxzi", "    const double lopxzi = log(1 + x*zi);"}, {"li2omxzi", "    const double li2omxzi = dilog(1 - x*zi);"}, {"sqrtxz1", "    const double sqrtxz1 = sqrt(1 - 2*z + z2 + 4*x*z);"}, {"sqrtxz1i", "    const double sqrtxz1i = 1. / sqrtxz1;"}, {"poly2", "    const double poly2 = 1 + 2*x + x2 - 4*x*z;"}, {"poly2i", "    const double poly2i = 1. / poly2;"}, {"poly2i2", "    const double poly2i2 = poly2i * poly2i;"}, {"sqrtxz2", "    const double sqrtxz2 = sqrt(poly2);"}, {"sqrtxz2i", "    const double sqrtxz2i = 1. / sqrtxz2;"}, {"sqrtxz3", "    const double sqrtxz3 = sqrt(x/z);"}, {"spec1i", "    const double spec1i = 1. / ( 1 + sqrtxz1 - z );"}, {"lspec1", "    const double lspec1 = log(1 + sqrtxz1 - z);"}, {"lspec1_2", "    const double lspec1_2 = lspec1 * lspec1;"}, {"lspec2", "    const double lspec2 = log(1 + sqrtxz1 + z);"}, {"lspec3", "    const double lspec3 = log(sqrtxz3);"}, {"lspec4", "    const double lspec4 = log(sqrtxz3*z);"}, {"lspec5", "    const double lspec5 = log(1 - sqrtxz2 + x);"}, {"lspec6", "    const double lspec6 = log(1 + sqrtxz2 + x);"}, {"lspec7", "    const double lspec7 = log(1 - 2*z + 4*x*z + z2);"}, {"lspec7_2", "    const double lspec7_2 = lspec7 * lspec7;"}, {"li2spec1", "    const double li2spec1 = dilog(0.5 - sqrtxz1/2. - z/2.);"}, {"li2spec2", "    const double li2spec2 = dilog(0.5 - sqrtxz1/2. + z/2.);"}, {"li2spec3", "    const double li2spec3 = dilog(0.5 - zi/2. - (sqrtxz1*zi)/2.);"}, {"li2spec4", "    const double li2spec4 = dilog(0.5 + zi/2. - (sqrtxz1*zi)/2.);"}, {"li2spec5", "    const double li2spec5 = dilog(0.5 - sqrtxz2/2. - x/2.);"}, {"li2spec6", "    const double li2spec6 = dilog(0.5 + sqrtxz2/2. - x/2.);"}, {"li2spec7", "    const double li2spec7 = dilog(0.5 - xi/2. - (sqrtxz2*xi)/2.);"}, {"li2spec8", "    const double li2spec8 = dilog(0.5 - xi/2. + (sqrtxz2*xi)/2.);"}, {"li2spec9", "    const double li2spec9 = dilog((1 - z)*omxi);"}, {"li2spec10", "    const double li2spec10 = dilog(x*(1 - z)*omxi*zi);"}, {"li2spec11", "    const double li2spec11 = dilog((1 - x)*omzi);"}, {"li2spec12", "    const double li2spec12 = dilog((1 - x)*xiz*omzi);"}, {"li2spec13", "    const double li2spec13 = dilog(z*omxi);"}, {"li2spec14", "    const double li2spec14 = dilog((1 - z)*z*omxi*xi);"}, {"li2spec15", "    const double li2spec15 = dilog(x*z*omxi*omzi);"}, {"li2spec16", "    const double li2spec16 = dilog((1 - x)*zi);"}, {"li2spec17", "    const double li2spec17 = dilog((1 - x)*(1 - z)*xizi);"}, {"li2spec18", "    const double li2spec18 = dilog((1 - x)*x*omzi*zi);"}, {"li2spec19", "    const double li2spec19 = dilog(-x*z);"}, {"li2spec20", "    const double li2spec20 = dilog(-x*zi);"}, {"li2spec21", "    const double li2spec21 = dilog((1 - z)*sqrtxz1i);"}, {"li2spec22", "    const double li2spec22 = dilog(-spec1i + sqrtxz1*spec1i + z*spec1i);"}, {"li2spec23", "    const double li2spec23 = dilog(sqrtxz1*mopzi);"}, {"li2spec24", "    const double li2spec24 = dilog((1 - x)*z*xi*omzi);"}, {"li2spec25", "    const double li2spec25 = dilog(x*z*omxi*omzi);"}, {"li2spec26", "    const double li2spec26 = dilog((1 - x)*(1 - z)*xi*zi);"}, {"li2spec27", "    const double li2spec27 = dilog((1 - z)*pow(sqrtxz1,-1));"}, {"atanspec1", "    const double atanspec1 = atan(sqrtxz3);"}, {"atanspec2", "    const double atanspec2 = atan(sqrtxz3*z);"}, {"itani1", "    const double itani1 = InvTanInt(-sqrtxz3);"}, {"itani2", "    const double itani2 = InvTanInt(sqrtxz3);"}, {"itani3", "    const double itani3 = InvTanInt(sqrtxz3*z);"}, {"Tr1", "    const double Tr1 = (z < 1 - x ? 1 : 0);"}, {"Tr2", "    const double Tr2 = (z > 1 - x ? 1 : 0);"}, {"Tt1", "    const double Tt1 = (z > x ? 1 : 0);"}, {"Tt2", "    const double Tt2 = (z < x ? 1 : 0);"}, {"Tu1", "    const double Tu1 = (z < 1 - x && z < x ? 1 : 0);"}, {"Tu2", "    const double Tu2 = (z > 1 - x && z < x ? 1 : 0);"}, {"Tu3", "    const double Tu3 = (z < 1 - x && z > x ? 1 : 0);"}, {"Tu4", "    const double Tu4 = (z > 1 - x && z > x ? 1 : 0);"}};