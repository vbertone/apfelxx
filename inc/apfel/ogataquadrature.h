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
     * @param nZeroMax: maximum number of terms in the Ogata quadrature (default: 1000)
     * @note Note that the default value of the parameter 'h' (0.001)
     * is based of studies of Drell-Yan transverse-momentum
     * distributions. However, this values could possibly be badly
     * non-optimal in other contexts. A good value of 'h' should take
     * into account the decay rate of the integrand. In particular,
     * the faster the decay the smaller the value of 'h' has to be.
     */
    OgataQuadrature(int    const& nu = 0,
                    double const& CutOff = 1e-5,
                    double const& h = 0.001,
                    int    const& nZeroMax = 1000);

    /**
     * @brief Function that initialises the coordinates and the
     * weights for the Ogata quadrature.
     * @param nZeroMax: maximum number of terms in the Ogata quadrature
     */
    void InitialiseWeights(int const& nZeroMax);

    /**
     * @brief Function that transform the input function.
     * @param func: function to be transformed
     * @param qT: value of qT in which to compute the transform
     * @param Dynh: switch to compute the step parameter _h dynamically (default: true)
     * @param nmax: maximum number of terms in the Ogata quadrature (default: 1000)
     * @param period: interval across which the integral is checked to determine where the sum is truncated (default: 10)
     * @return the value of the transform
     */
    template<typename T>
    T transform(std::function<T(double const&)> const& func, double const& qT, bool const& Dynh = true, int const& nmax = 1000, int const& period = 10) const;

    /**
     * @brief Function that returns the Bessel order.
     * @return the Bessel order
     */
    double GetBesselOrder() const { return _nu; }

    /**
     * @brief Function that returns the Ogata cut-off parameter.
     * @return the cut-off parameter
     */
    double GetCutOff() const { return _CutOff; }

    /**
     * @brief Function that returns the Ogata step parameter.
     * @return the step parameter
     */
    double GetStepParameter() const { return _h; }

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
     * @brief Function that sets the Ogata step parameter.
     * @param h: step parameter
     */
    void SetStepParameter(double const& h) { _h = h; }

    /**
     * @brief Function that writes on screan the first 1000 zeros of the
     * Bessel function J0. This function essentially generates the
     * std::vector<double> j0Zeros above. This function requires BOOST
     * and thus it is not available by defaults.
     * @param nu: the order of the Bessel function
     */
    void JnuZerosGenerator(int const& nu) const;

  private:
    int                 const _nu;       //!< The Bessel order
    double              const _CutOff;   //!< The target accuracy parameter
    double                    _h;        //!< The step parameter
    int                 const _nZeroMax; //!< The maximum number of zero's initialised
    std::vector<double>       _xf;       //!< Unscaled coordinates
    std::vector<double>       _weights;  //!< Weights of the quadrature
  };
}
