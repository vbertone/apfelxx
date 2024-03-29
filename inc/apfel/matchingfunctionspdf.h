//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include "apfel/expression.h"

namespace apfel
{
  /**
   * @defgroup TMDMatchingFunctions TMD matching functions
   * The perturbative matching functions for PDFs and FFs (references:
   * https://arxiv.org/pdf/2012.03256.pdf,
   * https://arxiv.org/pdf/1604.07869.pdf, and
   * https://arxiv.org/pdf/1706.01473.pdf). Notice that the
   * expressions for FFs are implicitly multiplied by a factor
   * x<SUP>2</SUP>. This is because the convolution with FFs has the
   * form:
   *
   * \f$D(x) = C(x) \otimes d(x) / x^2 = (1/x^2) \int_x^1 (dy/y) [y^2 C(y)] d(x/y)\f$
   *
   * For the implementation of the O(&alpha;<SUB>s</SUB><SUP>2</SUP>)
   * and O(&alpha;<SUB>s</SUB><SUP>3</SUP>) expressions I have used
   * the parameterization provided with
   * https://arxiv.org/pdf/2012.03256.pdf. A notebook with the
   * expressions can be found in docs/latex/src/codes/.
   */
  ///@{
  ///@}
  /**
   * @defgroup SLMatchFunc Space-like matching functions
   * TMD mathing functions for PDFs (space-like hard scales)
   * @ingroup TMDMatchingFunctions
   */
  ///@{
  ///@}
  /**
   * @defgroup NLOmatch NLO matching functions for PDFs
   * @ingroup SLMatchFunc
   */
  ///@{
  /**
   * @brief The O(&alpha;<SUB>s</SUB>) non-singlet matching function
   * for PDFs (references: https://arxiv.org/pdf/1604.07869.pdf and
   * https://arxiv.org/pdf/1706.01473.pdf).
   */
  class C1nspdf: public Expression
  {
  public:
    C1nspdf();
    double Regular(double const& x) const;
    double Local(double const&)     const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB>) quark-gluon matching function
   * for PDFs (references: https://arxiv.org/pdf/1604.07869.pdf and
   * https://arxiv.org/pdf/1706.01473.pdf).
   */
  class C1qgpdf: public Expression
  {
  public:
    C1qgpdf();
    double Regular(double const& x) const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB>) gluon-quark matching function
   * for PDFs (references: https://arxiv.org/pdf/1604.07869.pdf and
   * https://arxiv.org/pdf/1706.01473.pdf).
   */
  class C1gqpdf: public Expression
  {
  public:
    C1gqpdf();
    double Regular(double const& x) const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB>) gluon-gluon matching function
   * for PDFs (references: https://arxiv.org/pdf/1604.07869.pdf and
   * https://arxiv.org/pdf/1706.01473.pdf).
   */
  class C1ggpdf: public Expression
  {
  public:
    C1ggpdf();
    double Local(double const&) const;
  };
  ///@}

  /**
   * @defgroup NNLOmatch NNLO matching functions for PDFs
   * @ingroup SLMatchFunc
   */
  ///@{
  /**
   * @brief The O(&alpha;<SUB>s</SUB><SUP>2</SUP>) valence plus
   * matching function for PDFs (reference:
   * https://arxiv.org/pdf/2012.03256.pdf).
   */
  class C2nsppdf: public Expression
  {
  public:
    C2nsppdf(int const& nf);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  protected:
    int const _nf;
    double    _A2;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB><SUP>2</SUP>) valence minus
   * matching function for PDFs (reference:
   * https://arxiv.org/pdf/2012.03256.pdf).
   */
  class C2nsmpdf: public Expression
  {
  public:
    C2nsmpdf(int const& nf);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  protected:
    int const _nf;
    double    _A2;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB><SUP>2</SUP>) pure-singlet
   * matching function for PDFs (reference:
   * https://arxiv.org/pdf/2012.03256.pdf).
   */
  class C2pspdf: public Expression
  {
  public:
    C2pspdf();
    double Regular(double const& x) const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB><SUP>2</SUP>) quark-gluon
   * matching function for PDFs (reference:
   * https://arxiv.org/pdf/2012.03256.pdf).
   */
  class C2qgpdf: public Expression
  {
  public:
    C2qgpdf();
    double Regular(double const& x) const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB><SUP>2</SUP>) gluon-quark
   * matching function for PDFs (reference:
   * https://arxiv.org/pdf/2012.03256.pdf).
   */
  class C2gqpdf: public Expression
  {
  public:
    C2gqpdf(int const& nf);
    double Regular(double const& x) const;
  private:
    int const _nf;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB><SUP>2</SUP>) gluon-gluon
   * matching function for PDFs (reference:
   * https://arxiv.org/pdf/2012.03256.pdf).
   */
  class C2ggpdf: public Expression
  {
  public:
    C2ggpdf(int const& nf);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    int const _nf;
    double    _A2;
  };
  ///@}

  /**
   * @defgroup NNNLOmatch NNNLO matching functions for PDFs
   * @ingroup SLMatchFunc
   */
  ///@{
  /**
   * @brief The O(&alpha;<SUB>s</SUB><SUP>3</SUP>) valence plus
   * matching function for PDFs (reference:
   * https://arxiv.org/pdf/2012.03256.pdf).
   */
  class C3nsppdf: public Expression
  {
  public:
    C3nsppdf(int const& nf);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  protected:
    int const _nf;
    double    _A2;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB><SUP>3</SUP>) valence minus
   * matching function for PDFs (reference:
   * https://arxiv.org/pdf/2012.03256.pdf).
   */
  class C3nsmpdf: public Expression
  {
  public:
    C3nsmpdf(int const& nf);
    double Regular(double const& x) const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  protected:
    int const _nf;
    double    _A2;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB><SUP>3</SUP>) pure-valence
   * matching function for PDFs (reference:
   * https://arxiv.org/pdf/2012.03256.pdf).
   */
  class C3pvpdf: public Expression
  {
  public:
    C3pvpdf();
    double Regular(double const& x) const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB><SUP>3</SUP>) pure-singlet
   * matching function for PDFs (reference:
   * https://arxiv.org/pdf/2012.03256.pdf).
   */
  class C3pspdf: public Expression
  {
  public:
    C3pspdf(int const& nf);
    double Regular(double const& x) const;
  protected:
    int const _nf;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB><SUP>3</SUP>) quark-gluon
   * matching function for PDFs (reference:
   * https://arxiv.org/pdf/2012.03256.pdf).
   */
  class C3qgpdf: public Expression
  {
  public:
    C3qgpdf(int const& nf);
    double Regular(double const& x) const;
  protected:
    int const _nf;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB><SUP>3</SUP>) gluon-quark
   * matching function for PDFs (reference:
   * https://arxiv.org/pdf/2012.03256.pdf).
   */
  class C3gqpdf: public Expression
  {
  public:
    C3gqpdf(int const& nf);
    double Regular(double const& x) const;
  protected:
    int const _nf;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB><SUP>3</SUP>) gluon-gluon
   * matching function for PDFs (reference:
   * https://arxiv.org/pdf/2012.03256.pdf).
   */
  class C3ggpdf: public Expression
  {
  public:
    C3ggpdf(int const& nf);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  protected:
    int const _nf;
    double    _A2;
  };
  ///@}

  /**
   * @defgroup NLOBM NLO matching functions for Boer-Mulders PDFs
   * NLO matching functions for linearly polarised gluon PDF (Boer-Mulders)
   * @ingroup SLMatchFunc
   */
  ///@{
  /**
   * @brief The O(&alpha;<SUB>s</SUB>) gluon-quark matching function
   * for linearly polarised gluon PDF (reference:
   * https://arxiv.org/pdf/1907.03780.pdf).
   */
  class C1gqpdfBM: public Expression
  {
  public:
    C1gqpdfBM();
    double Regular(double const& x) const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB>) gluon-gluon matching function
   * for linearly polarised gluon PDF (reference:
   * https://arxiv.org/pdf/1907.03780.pdf).
   */
  class C1ggpdfBM: public Expression
  {
  public:
    C1ggpdfBM();
    double Regular(double const& x) const;
  };
  ///@}

  /**
   * @defgroup NNLOBM NNLO matching functions for Boer-Mulders PDFs
   * NNLO matching functions for linearly polarised gluon PDF (Boer-Mulders)
   * @ingroup SLMatchFunc
   */
  ///@{
  /**
   * @brief The O(&alpha;<SUB>s</SUB><SUP>2</SUP>) gluon-quark
   * matching function for linearly polarised gluon PDF (reference:
   * https://arxiv.org/pdf/1907.03780.pdf).
   */
  class C2gqpdfBM: public Expression
  {
  public:
    C2gqpdfBM(int const& nf);
    double Regular(double const& x) const;
  protected:
    int const _nf;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB><SUP>2</SUP>) gluon-gluon
   * matching function for linearly polarised gluon PDFs (reference:
   * https://arxiv.org/pdf/1907.03780.pdf).
   */
  class C2ggpdfBM: public Expression
  {
  public:
    C2ggpdfBM(int const& nf);
    double Regular(double const& x) const;
  protected:
    int const _nf;
  };
  ///@}

  /**
  * @defgroup NLOSivers NLO matching functions for Sivers quark PDFs
  * NLO matching functions for Sivers quark PDFs
  * @ingroup SLMatchFunc
  */
  ///@{
  /**
   * @brief The O(&alpha;<SUB>s</SUB>) non-singlet matching function
   * for Sivers PDFs (see Eq. (A.9) of
   * https://arxiv.org/pdf/2009.10710.pdf).
   */
  class C1nspdfSivers: public Expression
  {
  public:
    C1nspdfSivers();
    double Regular(double const& x) const;
    double Local(double const&)     const;
  };
  ///@}

  /**
   * @defgroup NLOg1 NLO matching functions for helicity PDFs g1
   * NLO matching functions for the helicity PDFs g1
   * @ingroup SLMatchFunc
   */
  ///@{
  /**
   * @brief The O(&alpha;<SUB>s</SUB>) non-singlet matching function
   * for g1 PDFs (reference: https://arxiv.org/pdf/1702.06558.pdf).
   */
  class C1nspdfg1: public Expression
  {
  public:
    C1nspdfg1();
    double Regular(double const& x) const;
    double Local(double const&)     const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB>) quark-gluon matching function
   * for g1 PDFs (reference: https://arxiv.org/pdf/1702.06558.pdf).
   */
  class C1qgpdfg1: public Expression
  {
  public:
    C1qgpdfg1();
    double Regular(double const& x) const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB>) gluon-quark matching function
   * for g1 PDFs (reference: https://arxiv.org/pdf/1702.06558.pdf).
   */
  class C1gqpdfg1: public Expression
  {
  public:
    C1gqpdfg1();
    double Regular(double const& x) const;
  };

  /**
   * @brief The O(&alpha;<SUB>s</SUB>) gluon-gluon matching function
   * for g1 PDFs (reference: https://arxiv.org/pdf/1702.06558.pdf).
   */
  class C1ggpdfg1: public Expression
  {
  public:
    C1ggpdfg1();
    double Regular(double const& x) const;
    double Local(double const&)     const;
  };
  ///@}
}
