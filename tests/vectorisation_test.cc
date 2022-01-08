//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include <apfel/apfelxx.h>
#if __AVX512F__ == 1
#include <immintrin.h>

// Use one of these commands from shell to list the environment
// variables set by the compiler:
//
// clang++ -mavx512f -dM -E - < /dev/null | egrep "SSE|AVX" | sort
// clang++ -mavx2 -dM -E - < /dev/null | egrep "SSE|AVX" | sort
// clang++ -msse3 -dM -E - < /dev/null | egrep "SSE|AVX" | sort

double DotProduct(const double *a, const double *b, int const& N)
{
  double answer = 0;
  for(int i = 0; i < N; i++)
    answer += a[i] * b[i];
  return answer;
}

double DotProductAVX2_2(const double *a, const double *b, int const& N)
{
  // Initialise
  __m256d sum = _mm256_set_pd(0, 0, 0, 0);

  // Add up partial dot-products in blocks of 256 bits
  for (int i = 0; i < N - 4; i += 4)
    {
      __m256d x = _mm256_loadu_pd(a + i);
      __m256d y = _mm256_loadu_pd(b + i);
      sum = _mm256_add_pd(sum, _mm256_mul_pd(x, y));
    }
  // Grab the first 128 bits, grab the lower 128 bits and then add
  // them together with a 128 bit add instruction.
  __m256d hsm = _mm256_hadd_pd(sum, sum);
  __m128d dot = _mm_add_pd(_mm256_extractf128_pd(hsm, 1), _mm256_castpd256_pd128(hsm));

  // Find the partial dot-product for the remaining elements after
  // dealing with all 256-bit blocks.
  double rest = 0;
  for (int i = N - N%4; i < N; i++)
    rest += a[i] * b[i];

  return ((double*)&dot)[0] + rest;
}

double DotProductAVX2(const double *x, const double *y, int const& N)
{
  // divide cache line into two variables handled by the same processor
  __m256d res1 = _mm256_setzero_pd();
  __m256d res2 = _mm256_setzero_pd();

  for (int i = 0; i < N - 4; i += 8)
    {
      // use partial loop unrolling in order to cover entire cache line
      res1 = _mm256_fmadd_pd(_mm256_loadu_pd(x + i),     _mm256_loadu_pd(y + i),     res1);
      res2 = _mm256_fmadd_pd(_mm256_loadu_pd(x + i + 4), _mm256_loadu_pd(y + i + 4), res2);
    }

  // reduce all intrinsics to single double
  __m256d res = _mm256_add_pd(res1, res2);
  __m256d sum = _mm256_hadd_pd(res, res);
  return ((double*)&sum)[0] + ((double*)&sum)[2];
}

double DotProductAVX512(const double *x, const double *y, int const& N)
{
  __m512d res = _mm512_setzero_pd();
  for (int i = 0; i < N; i += 8)
    res = _mm512_fmadd_pd(_mm512_loadu_pd(x + i), _mm512_loadu_pd(y + i), res);

  // Reduce intrinsic to single double
  return _mm512_reduce_add_pd(res);
}

int main()
{
  const int N = 3095;
  double a[N], b[N];

  for (int i = 0; i < N; i++)
    a[i] = b[i] = i / sqrt(N);

  apfel::Timer t;
  double pr = 0;
  for (int i = 0; i < 1000000; i++)
    pr += DotProduct(a, b, N);
  std::cout << pr << std::endl;
  t.stop();

  t.start();
  double p = 0;
  for (int i = 0; i < 1000000; i++)
    p += DotProductAVX2(a, b, N);
  std::cout << ( p - pr ) / pr << std::endl;
  t.stop();

  t.start();
  p = 0;
  for (int i = 0; i < 1000000; i++)
    p += DotProductAVX512(a, b, N);
  std::cout << ( p - pr ) / pr  << std::endl;
  t.stop();

  t.start();
  p = 0;
  for (int i = 0; i < 1000000; i++)
    p += DotProductAVX2_2(a, b, N);
  std::cout << ( p - pr ) / pr  << std::endl;
  t.stop();

  return 0;
}
#else
int main()
{
  return 0;
}
#endif
