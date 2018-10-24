//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
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
     * @param CutOff: the accuracy computed as a cutoff on the size of the last computed term relative to the total (default: 10<SUP>-5</SUP>)
     * @param h: internal variable of the algorithm (default: 0.001)
     */
    OgataQuadrature(double const& CutOff = 1e-5,
		    double const& h = 0.001);

    /**
     * @brief Function that integrates the integrand with a given
     * relative accuracy.
     * @param xmin: the lower bound integration bound
     * @param xmax: the upper bound integration bound
     * @param eps: the required relative accuracy
     * @return the value of the integral
     */
    double transform(std::function<double(double const&)> const& func, double const& qT) const;

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
     */
    void J0ZerosGenerator() const;

  private:
    double              const _CutOff;  //!< The target accuracy parameter
    double              const _h;       //!< The step parameter
    std::vector<double>       _xf;      //!< Unscaled coordinates
    std::vector<double>       _weights; //!< Weights of the quadrature
  };
}
