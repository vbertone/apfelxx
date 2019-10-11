//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/doubleexponentialquadrature.h"

#include <math.h>

namespace apfel
{
  //_____________________________________________________________________________
  DoubleExponentialQuadrature::DoubleExponentialQuadrature(double const& eps)
  {
    /* ---- adjustable parameter ---- */
    int lmax = 5;
    double efs = 0.1, enoff = 0.40, pqoff = 2.9, ppoff = -0.72;
    /* ------------------------------ */
    int noff0, nk0, noff, k, nk, j, lenaw;
    double pi4, tinyln, epsln, frq4, per2, pp, pq, ehp, ehm, h, t,
           ep, em, tk, xw, wg, xa, tiny;

    lenaw = 8000;
    tiny = 1e-30;
    pi4 = atan(1.0);
    tinyln = -log(tiny);
    epsln = 1 - log(efs * eps);
    frq4 = 1 / (2 * pi4);
    per2 = 4 * pi4;
    pq = pqoff / epsln;
    pp = ppoff - log(pq * pq * frq4);
    ehp = exp(2 * pq);
    ehm = 1 / ehp;
    _aw[3] = lmax;
    _aw[4] = eps;
    _aw[5] = sqrt(efs * eps);
    noff0 = 6;
    nk0 = 1 + (int) (enoff * epsln);
    _aw[1] = nk0;
    noff = 2 * nk0 + noff0;
    wg = 0;
    xw = 1;
    for (k = 1; k <= nk0; k++)
      {
        wg += xw;
        _aw[noff - 2 * k] = wg;
        _aw[noff - 2 * k + 1] = xw;
        xw = xw * (nk0 - k) / k;
      }
    wg = per2 / wg;
    for (k = noff0; k <= noff - 2; k += 2)
      {
        _aw[k] *= wg;
        _aw[k + 1] *= wg;
      }
    xw = exp(pp - 2 * pi4);
    _aw[noff] = sqrt(xw * (per2 * 0.5));
    _aw[noff + 1] = xw * pq;
    _aw[noff + 2] = per2 * 0.5;
    h = 2;
    nk = 0;
    k = noff + 3;
    do
      {
        t = h * 0.5;
        do
          {
            em = exp(2 * pq * t);
            ep = pi4 * em;
            em = pi4 / em;
            tk = t;
            j = k;
            do
              {
                xw = exp(pp - ep - em);
                wg = sqrt(frq4 * xw + tk * tk);
                xa = xw / (tk + wg);
                wg = (pq * xw * (ep - em) + xa) / wg;
                _aw[j] = xa;
                _aw[j + 1] = xw * pq;
                _aw[j + 2] = wg;
                ep *= ehp;
                em *= ehm;
                tk += 1;
                j += 3;
              }
            while (ep < tinyln && j <= lenaw - 3);
            t += h;
            k += nk;
          }
        while (t < 1);
        h *= 0.5;
        if (nk == 0)
          {
            if (j > lenaw - 6) j -= 3;
            nk = j - noff;
            k += nk;
            _aw[2] = nk;
          }
      }
    while (2 * k - noff - 3 <= lenaw);
    _aw[0] = k - 3;
  }

  //_____________________________________________________________________________
  double DoubleExponentialQuadrature::transform(std::function<double(double const&)> const& f, double const& qT) const
  {
    int lenawm, nk0, noff0, nk, noff, lmax, m, k, j, jm, l;
    double eps, per, perw, w02, ir, h, iback, irback, t, tk,
           xa, fm, fp, errh, s0, s1, s2, errd, i, err, a;

    a = 0;
    i = 0;
    errh = 0;
    fm = 0;
    fp = 0;
    lenawm = (int) (_aw[0] + 0.5);
    nk0 = (int) (_aw[1] + 0.5);
    noff0 = 6;
    nk = (int) (_aw[2] + 0.5);
    noff = 2 * nk0 + noff0;
    lmax = (int) (_aw[3] + 0.5);
    eps = _aw[4];
    per = 1 / fabs(qT);
    w02 = 2 * _aw[noff + 2];
    perw = per * w02;
    i = f(a + _aw[noff] * per) * j0(qT * (a + _aw[noff] * per));
    ir = i * _aw[noff + 1];
    i *= _aw[noff + 2];
    err = fabs(i);
    h = 2;
    m = 1;
    k = noff;
    do
      {
        iback = i;
        irback = ir;
        t = h * 0.5;
        do
          {
            if (k == noff)
              {
                tk = 1;
                k += nk;
                j = noff;
                do
                  {
                    j += 3;
                    xa = per * _aw[j];
                    fm = f(a + xa) * j0(qT * (a + xa));
                    fp = f(a + xa + perw * tk) * j0(qT * (a + xa + perw * tk));
                    ir += (fm + fp) * _aw[j + 1];
                    fm *= _aw[j + 2];
                    fp *= w02 - _aw[j + 2];
                    i += fm + fp;
                    err += fabs(fm) + fabs(fp);
                    tk += 1;
                  }
                while (_aw[j] > eps && j < k);
                errh = err * _aw[5];
                err *= eps;
                jm = j - noff;
              }
            else
              {
                tk = t;
                for (j = k + 3; j <= k + jm; j += 3)
                  {
                    xa = per * _aw[j];
                    fm = f(a + xa) * j0(qT * (a + xa));
                    fp = f(a + xa + perw * tk) * j0(qT * (a + xa + perw * tk));
                    ir += (fm + fp) * _aw[j + 1];
                    fm *= _aw[j + 2];
                    fp *= w02 - _aw[j + 2];
                    i += fm + fp;
                    tk += 1;
                  }
                j = k + jm;
                k += nk;
              }
            while (fabs(fm) > err && j < k)
              {
                j += 3;
                fm = f(a + per * _aw[j]) * j0(qT * (a + per * _aw[j]));
                ir += fm * _aw[j + 1];
                fm *= _aw[j + 2];
                i += fm;
              }
            fm = f(a + perw * tk) * j0(qT * (a + perw * tk));
            s2 = w02 * fm;
            i += s2;
            if (fabs(fp) > err || fabs(s2) > err)
              {
                l = 0;
                for (;;)
                  {
                    l++;
                    s0 = 0;
                    s1 = 0;
                    s2 = fm * _aw[noff0 + 1];
                    for (j = noff0 + 2; j <= noff - 2; j += 2)
                      {
                        tk += 1;
                        fm = f(a + perw * tk) * j0(qT * (a + perw * tk));
                        s0 += fm;
                        s1 += fm * _aw[j];
                        s2 += fm * _aw[j + 1];
                      }
                    if (s2 <= err || l >= lmax) break;
                    i += w02 * s0;
                  }
                i += s1;
                if (s2 > err) err = s2;
              }
            t += h;
          }
        while (t < 1);
        if (m == 1)
          {
            errd = 1 + 2 * errh;
          }
        else
          {
            errd = h * (fabs(i - 2 * iback) + fabs(ir - 2 * irback));
          }
        h *= 0.5;
        m *= 2;
      }
    while (errd > errh && 2 * k - noff <= lenawm);
    i *= h * per;
    if (errd > errh)
      {
        err = -errd * per;
      }
    else
      {
        err *= per * m * 0.5;
      }
    return i;
  }
}
