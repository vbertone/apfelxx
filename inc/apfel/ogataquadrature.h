//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include <vector>
#include <functional>

namespace apfel
{
  /**
   * @brief The OgataQuadrature class implements the Hankel-transform
   * of the input function using the Ogata quadrature method described
   * here:
   * http://www.kurims.kyoto-u.ac.jp/~okamoto/paper/Publ_RIMS_DE/41-4-40.pdf.
   */
  class OgataQuadrature
  {
  public:
    /**
     * @brief The Integrator constructor.
     * @param nu: the order of the Bessel function (default: 0)
     * @param CutOff: the accuracy computed as a cutoff on the size of the last computed term relative to the total (default: 10<SUP>-5</SUP>)
     * @param h: internal variable of the algorithm (default: 0.001)
     * @note Note that the default value of the parameter 'h' (0.001)
     * is based of studies of Drell-Yan transverse-momentum
     * distributions. However, this values could possibly be badly
     * non-optimal in other contexts. A good value of 'h' should take
     * into account the decay rate of the integrand. In particular,
     * the faster the decay the smaller the value of 'h' has to be.
     */
    OgataQuadrature(int    const& nu = 0,
                    double const& CutOff = 1e-5,
                    double const& h = 0.001);

    /**
     * @brief Function that transform the input function.
     * @param func: function to be transformed
     * @param qT: value of qT in which to compute the transform
     * @param nmax: maximum number of terms in the Ogata quadrature (default: 1000)
     * @return the value of the transform
     */
    template<typename T>
    T transform(std::function<T(double const&)> const& func, double const& qT, int const& nmax = 1000) const;

    /**
     * @brief Function that returns the unscaled coordinates used in
     * the Ogata quadrature.
     * @return the unscales coordinates
     */
    std::vector<double> GetCoordinates() const { return _xf; }

    /**
     * @brief Function that returns the weights used in the Ogata
     * quadrature.
     * @return the weights
     */
    std::vector<double> GetWeights() const { return _weights; }

    /**
     * @brief Function that writes on screan the first 1000 zeros of the
     * Bessel function J0. This function essentially generates the
     * std::vector<double> j0Zeros above. This function requires BOOST
     * and thus it is not available by defaults.
     * @param nu: the order of the Bessel function
     */
    void JnuZerosGenerator(int const& nu) const;

  private:
    double              const _CutOff;  //!< The target accuracy parameter
    double              const _h;       //!< The step parameter
    std::vector<double>       _xf;      //!< Unscaled coordinates
    std::vector<double>       _weights; //!< Weights of the quadrature
  };
}
