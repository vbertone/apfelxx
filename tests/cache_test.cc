#include <iostream>
#include <apfel/cache.h>
#include <apfel/timer.h>
using namespace apfel;
using namespace std;

/**
 * Example 1 with 1 variable
 */
class Fibonacci
{
public:
  Fibonacci()
  {
    _cache = unique_ptr<Cache<int,int>>(new Cache<int,int>());
  }

  int fib(int n)
  {
    if (n < 2) return n;
    return fib(n-1) + fib(n-2);
  }

  int fib_cache(int n)
  {
    int res = n;
    if (_cache->tryGet(n, res)) return res;

    if (n >= 2)
      res = fib_cache(n-1) + fib_cache(n-2);

    _cache->insert(n,res);
    return res;
  }

private:
  unique_ptr<Cache<int,int>> _cache;
};

/**
 * Example how to use with 2 or more variables
 */
class Test2D
{
public:
  Test2D()
  {
    _cache = unique_ptr<Cache<tuple<double,double>,double>>
        (new Cache<tuple<double,double>,double>());
  }

  double fun(double const& x, double const& y)
  {
    double res = 0;
    auto key = make_tuple(x,y);
    if (_cache->tryGet(key,res)) return res;
    for (int i = 0; i < 1e8; i++)
      {}
    _cache->insert(key,res);
    return res;
  }

private:
  unique_ptr<Cache<tuple<double,double>,double>> _cache;
};


int main()
{
  Fibonacci f;

  Timer t;
  t.start();
  cout << "Computing fibonacci without cache..." << endl;
  cout << "result f(45) = " << f.fib(45) << endl;
  cout << "result f(44) = " << f.fib(44) << endl;
  cout << "result f(43) = " << f.fib(43) << endl;
  cout << "result f(42) = " << f.fib(42) << endl;
  cout << "result f(41) = " << f.fib(41) << endl;
  t.stop();

  t.start();
  cout << "Computing fibonacci with cache..." << endl;
  cout << "result f(45) = " << f.fib_cache(45) << endl;
  cout << "result f(44) = " << f.fib_cache(44) << endl;
  cout << "result f(43) = " << f.fib_cache(43) << endl;
  cout << "result f(42) = " << f.fib_cache(42) << endl;
  cout << "result f(41) = " << f.fib_cache(41) << endl;
  t.stop();

  Test2D t2d;
  t.start();
  cout << "Computing 2d example with cache..." << endl;
  cout << "result f(1,2) = " << t2d.fun(1,2) << endl;
  t.stop();

  return 0;
}
