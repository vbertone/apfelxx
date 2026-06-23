//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include <apfel/apfelxx.h>

#include <sstream>
#include <iostream>
#include <cmath>
#include <cstring>
#include <cstdio>

// Non-regression test for the cereal binary (de)serialisation of
// DoubleOperator. A DoubleOperator is built on two overlapping
// multi-subgrid grids (the configuration that exercises the grid
// "locking" performed at construction) and round-tripped both through
// an in-memory stream and through a file (the std::string overloads).
// In each case the whole nested payload is compared bit-for-bit
// (memcmp on each double, not within a tolerance).

// Compare two DoubleOperator payloads bit-for-bit. Returns true if the
// shapes match; reports the element and bit-difference counts.
bool BitCompare(std::vector<std::vector<apfel::matrix<apfel::matrix<double>>>> const& a,
                std::vector<std::vector<apfel::matrix<apfel::matrix<double>>>> const& b,
                std::size_t& nelem, std::size_t& bitdiffs)
{
  nelem = 0;
  bitdiffs = 0;
  bool shape = (a.size() == b.size());
  for (std::size_t i = 0; shape && i < a.size(); i++)
    {
      shape = (a[i].size() == b[i].size());
      for (std::size_t j = 0; shape && j < a[i].size(); j++)
        {
          shape = (a[i][j].size(0) == b[i][j].size(0) && a[i][j].size(1) == b[i][j].size(1));
          for (std::size_t r = 0; shape && r < a[i][j].size(0); r++)
            for (std::size_t c = 0; shape && c < a[i][j].size(1); c++)
              {
                const apfel::matrix<double>& ma = a[i][j](r, c);
                const apfel::matrix<double>& mb = b[i][j](r, c);
                shape = (ma.size(0) == mb.size(0) && ma.size(1) == mb.size(1));
                for (std::size_t p = 0; shape && p < ma.size(0); p++)
                  for (std::size_t q = 0; q < ma.size(1); q++)
                    {
                      const double va = ma(p, q), vb = mb(p, q);
                      nelem++;
                      if (std::memcmp(&va, &vb, sizeof(double)) != 0)
                        bitdiffs++;
                    }
              }
        }
    }
  return shape;
}

int main()
{
  // Two distinct grids, each with overlapping subgrids.
  const apfel::Grid gx{{apfel::SubGrid{20, 1e-3, 3}, apfel::SubGrid{10, 1e-1, 3}}};
  const apfel::Grid gz{{apfel::SubGrid{15, 1e-2, 3}, apfel::SubGrid{8, 5e-1, 3}}};

  // Build a DoubleOperator from a concrete DoubleExpression.
  const apfel::DoubleOperator O{gx, gz, apfel::DoubleIdentity{}};
  const std::vector<std::vector<apfel::matrix<apfel::matrix<double>>>> a = O.GetDoubleOperator();

  // 1) Round-trip through an in-memory portable-binary stream.
  std::stringstream ss(std::ios::in | std::ios::out | std::ios::binary);
  O.EmitDoubleOperatorBinary(ss);
  ss.seekg(0);
  const apfel::DoubleOperator Rs = apfel::DoubleOperator::ReadBinary(ss, gx, gz, apfel::DoubleIdentity{});
  std::size_t nelem_s = 0, bitdiffs_s = 0;
  const bool shape_s = BitCompare(a, Rs.GetDoubleOperator(), nelem_s, bitdiffs_s);

  // 2) Round-trip through a file (external-grid std::string overloads).
  const std::string path = "cereal_roundtrip_test.bin";
  O.EmitDoubleOperatorBinary(path);
  const apfel::DoubleOperator Rf = apfel::DoubleOperator::ReadBinary(path, gx, gz, apfel::DoubleIdentity{});
  std::size_t nelem_f = 0, bitdiffs_f = 0;
  const bool shape_f = BitCompare(a, Rf.GetDoubleOperator(), nelem_f, bitdiffs_f);

  // 3) Self-contained round-trip: read with no external grids and no
  // DoubleExpression. The grids are rebuilt from the file and owned by
  // the returned object; they must value-equal the originals, and the
  // payload must still be bit-exact.
  const apfel::DoubleOperator Rc = apfel::DoubleOperator::ReadBinary(path);
  std::remove(path.c_str());
  std::size_t nelem_c = 0, bitdiffs_c = 0;
  const bool shape_c = BitCompare(a, Rc.GetDoubleOperator(), nelem_c, bitdiffs_c);
  const bool gridsmatch = (Rc.GetFirstGrid() == gx && Rc.GetSecondGrid() == gz);

  // 4) A self-contained operator owns its grids through shared_ptr; the
  // new owning constructor must keep the grid references valid across
  // copy and move. Exercise both and re-check grids and payload.
  apfel::DoubleOperator Rcopy = Rc;                  // copy construction
  std::size_t nelem_cp = 0, bitdiffs_cp = 0;
  const bool shape_cp = BitCompare(a, Rcopy.GetDoubleOperator(), nelem_cp, bitdiffs_cp);
  const bool copygrids = (Rcopy.GetFirstGrid() == gx && Rcopy.GetSecondGrid() == gz);

  apfel::DoubleOperator Rmove = std::move(Rcopy);    // move construction
  std::size_t nelem_mv = 0, bitdiffs_mv = 0;
  const bool shape_mv = BitCompare(a, Rmove.GetDoubleOperator(), nelem_mv, bitdiffs_mv);
  const bool movegrids = (Rmove.GetFirstGrid() == gx && Rmove.GetSecondGrid() == gz);

  const bool copymove_ok = shape_cp && shape_mv && copygrids && movegrids
                           && nelem_cp == nelem_s && nelem_mv == nelem_s
                           && bitdiffs_cp == 0 && bitdiffs_mv == 0;

  const bool epsmatch  = (O.GetIntegrationAccuracy() == Rs.GetIntegrationAccuracy()
                          && O.GetIntegrationAccuracy() == Rf.GetIntegrationAccuracy()
                          && O.GetIntegrationAccuracy() == Rc.GetIntegrationAccuracy());
  const bool namematch = (O.GetDoubleExpressionName() == Rs.GetDoubleExpressionName()
                          && O.GetDoubleExpressionName() == Rf.GetDoubleExpressionName()
                          && O.GetDoubleExpressionName() == Rc.GetDoubleExpressionName());

  std::cout << "shape consistent       : " << ((shape_s && shape_f && shape_c) ? "yes" : "no") << "\n";
  std::cout << "elements               : " << nelem_s << "\n";
  std::cout << "bit differences (mem)  : " << bitdiffs_s << "\n";
  std::cout << "bit differences (file) : " << bitdiffs_f << "\n";
  std::cout << "bit differences (self) : " << bitdiffs_c << "\n";
  std::cout << "self-contained grids   : " << (gridsmatch ? "yes" : "no") << "\n";
  std::cout << "owning copy/move ok    : " << (copymove_ok ? "yes" : "no") << "\n";
  std::cout << "eps match              : " << (epsmatch ? "yes" : "no") << "\n";
  std::cout << "name match             : " << (namematch ? "yes" : "no") << "\n";

  const bool pass = shape_s && shape_f && shape_c && nelem_s > 0
                    && nelem_f == nelem_s && nelem_c == nelem_s
                    && bitdiffs_s == 0 && bitdiffs_f == 0 && bitdiffs_c == 0
                    && gridsmatch && copymove_ok && epsmatch && namematch;
  std::cout << (pass ? "ROUND-TRIP PASS (bit-exact)" : "ROUND-TRIP FAIL") << std::endl;

  return pass ? 0 : 1;
}
