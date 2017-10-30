#include <functional>
#include <vector>


#include "potential/electrostatics/Electrostatics.h"
#include "potential/electrostatics/aliases.h"

#define DEBUG // XXX: Probably pass this as a compiler flag


namespace {
  // Utility constants
  const double g34 = std::exp(elec::gammln(3.0/4.0));

  // Pure fmap for vector
  template <typename Fn, typename T>
  auto vecMap(const Fn& f, const Vec<T>& v) {
    Vec<decltype(f(v[0]))> mapped;
    std::transform(v.cbegin(), v.cend(), mapped.begin(), f);
    return mapped;
  }


  // Replicate a vector elementwise, e.g.
  //   replicateElt 2 [1, 2, 3, 4] = [1, 1, 2, 2, 3, 3, 4, 4]
  template <typename It>
  auto replicateElt(
    size_t n, const It& srcBegin, const It& srcEnd) {

    using Elt = typename std::iterator_traits<It>::value_type;

    // Sanity check
    auto srcSize = std::distance(srcEnd, srcBegin);
    if (srcSize < 0) {
      throw std::range_error("Invalid range for source");
    }

    Vec<Elt> replicated(srcSize * n);


    auto destIt = replicated.begin();
    for (auto it = srcBegin; it != srcEnd; ++it) {
      std::fill(destIt, destIt + n, *it);
      destIt += n;
    }

    return replicated;
  }

} // namespace


namespace elec {


  Electrostatics::Electrostatics(
		const Vec<double>& chg,
    const Vec<double>& polfac,
    const Vec<double>& pol,
    const Vec<double>& orig_xyz,
    const Vec<Id>& mon_id,
    const Vec<size_t>& sites,
    const Vec<size_t>& first_ind,
    const AssocVec<Id, size_t>& mon_type_count,
    const double tolerance,
    const size_t maxit,
    const bool do_grads,
    Vec<double>& grad) :

    // Init system params
    _chg(chg),
    _polfac(polfac),
    _pol(pol),
    _orig_xyz(orig_xyz),
    _mon_id(mon_id),
    _sites(sites),
    _first_ind(first_ind),
    _mon_type_count(mon_type_count),
    _tolerance(tolerance),
    _maxit(maxit),
    _do_grads(do_grads),
    _grad(grad),

    // Calculate which pairs to skip
    _notExcluded(getNonExcluded()),

		_nsites(chg.size()),

		// Temporaries
    _sqrtpol(
      // Squareroot all pols
      vecMap(
        [](double x) { return std::sqrt(x); },
        replicateElt(3, pol.cbegin(), pol.cend())
      )
    ),
		_Efq(3*_nsites),
    _Efd(3*_nsites),
    _phi(_nsites),
    _mu(3*_nsites),
    _xyz(3*_nsites),

    _phii(maxMonomerCount(),0.0),
    _phij(maxMonomerCount(),0.0),
    _Efqix(maxMonomerCount(),0.0),
    _Efqjx(maxMonomerCount(),0.0),
    _Efqiy(maxMonomerCount(),0.0),
    _Efqjy(maxMonomerCount(),0.0),
    _Efqiz(maxMonomerCount(),0.0),
    _Efqjz(maxMonomerCount(),0.0)
  {
    // It's in the name
    initTemps();
  }


  void Electrostatics::initTemps() {
    // Zero out everything
    std::fill(_Efq.begin(), _Efq.end(), 0.0);
    std::fill(_Efd.begin(), _Efd.end(), 0.0);
    std::fill(_phi.begin(), _phi.end(), 0.0);
    std::fill(_mu.begin(), _mu.end(), 0.0);

    std::fill(_phii.begin(), _phii.end(), 0.0);
    std::fill(_phij.begin(), _phij.end(), 0.0);
    std::fill(_Efqix.begin(), _Efqix.end(), 0.0);
    std::fill(_Efqjx.begin(), _Efqjx.end(), 0.0);
    std::fill(_Efqiy.begin(), _Efqiy.end(), 0.0);
    std::fill(_Efqjy.begin(), _Efqjy.end(), 0.0);
    std::fill(_Efqiz.begin(), _Efqiz.end(), 0.0);
    std::fill(_Efqjz.begin(), _Efqjz.end(), 0.0);

    // Order coordinates by site and monomer
    organizeCrds();
  }


  void Electrostatics::organizeCrds() {
    size_t fi_mon = 0;
    size_t fi_crd = 0;

    for (size_t mt = 0; mt < _mon_type_count.size(); mt++) {
      size_t ns = _sites[fi_mon];
      size_t nmon = _mon_type_count[mt].second;
      size_t nmon2 = nmon*2;
      for (size_t m = 0; m < nmon; m++) {
        size_t mns3 = m*ns*3;
        for (size_t i = 0; i < ns; i++) {
          size_t inmon3 = 3*i*nmon;
          _xyz[inmon3 + m + fi_crd] =
                 _orig_xyz[fi_crd + mns3 + 3*i];
          _xyz[inmon3 + m + fi_crd + nmon] =
                 _orig_xyz[fi_crd + mns3 + 3*i + 1];
          _xyz[inmon3 + m + fi_crd + nmon2] =
                 _orig_xyz[fi_crd + mns3 + 3*i + 2];
        }
      }
      fi_mon += nmon;
      fi_crd += nmon*ns*3;
    }
  }


  AssocVec<Pair<Id, size_t>, Vec<IndexPair>> Electrostatics::getNonExcluded() const {
    AssocVec<Pair<Id, size_t>, Vec<IndexPair>> notExcluded;

    excluded_set_type exc12;
    excluded_set_type exc13;
    excluded_set_type exc14;

    size_t fi_mon = 0;

    for (size_t mt = 0; mt < _mon_type_count.size(); mt++) {
      size_t ns = _sites[fi_mon];

      // FIXME/TODO: i dunno if the exclusion maps are actualy maps; if they
      // are, it'd be a lot easier just to copy or use them more directy
      systools::GetExcluded(_mon_id[fi_mon], exc12, exc13, exc14);

      // Make an entry for this monomer type
      notExcluded.push_back(
        std::make_pair(_mon_type_count[mt], Vec<IndexPair>{}));

      auto& newEntry = notExcluded.back();

      for (size_t i = 0; i < ns - 1; i++) {
        for (size_t j = i + 1; j < ns; j++) {
          // Skip if i and j are bonded
          // FIXME: verify the above^ is true; i'm just assuming
          bool is12 = systools::IsExcluded(exc12, i, j);
          bool is13 = systools::IsExcluded(exc13, i, j);
          bool is14 = systools::IsExcluded(exc14, i, j);
          if (is12 || is13 || is14) continue;

          // Add this pair to the list
          newEntry.second.push_back(std::make_pair(i, j));
        }
      }
    }

	  return notExcluded;
  }


  void Electrostatics::computePermanent() {
    computeIntraSites();
    computeInterSites();
  }


  void Electrostatics::recomputePermanent() {
    recomputeIntraSites();
    recomputeInterSites();
  }


  void Electrostatics::computeInterSites() {
    // TODO
  }


  void Electrostatics::recomputeIntraSites() {
    // TODO
  }


  void Electrostatics::recomputeInterSites() {
    // TODO
  }

	// Calculate screening coefficients (set as out params)
  void Electrostatics::screenCoeff(const double rsq, const double A,
    double& s0r, double& s1r3) {

    const double scaledA = std::pow(A,1.0/6.0);
    const double Ai = 1/scaledA;
    const double Asqsq = scaledA*scaledA*scaledA*scaledA;

    const double r = std::sqrt(rsq);
    const double ri = 1.0/r;
    const double risq = ri*ri;
    const double rsqsq = rsq * rsq;

    // Some values that will be used in the screening functions
    const double rA4 = rsqsq/Asqsq;
    const double arA4 = aCC*rA4;
    // TODO look at the exponential function intel vec
    const double exp1 = std::exp(-arA4);
    const double gq = gammq(3.0/4.0, arA4);
    //const double gq = 1.0;
    const double a_mrt = std::pow(aCC, 1.0/4.0);
    //const double a4 = aCC * 4.0;

    // Get screening functions
    const double s1r = ri - exp1*ri;

    s0r = s1r + a_mrt * Ai * g34 * gq;
    s1r3 = s1r * risq;
  }


  void Electrostatics::computeIntraSites() {
    size_t  fi_mon = 0, fi_sites = 0, fi_crd = 0;

    // FIXME: there's gotta be a prettier way to bind this (until c++17...)
    for (auto& typeCount_indexPair : _notExcluded) {
      auto& typeCount = typeCount_indexPair.first;
      auto& siteIndexPairs = typeCount_indexPair.second;

      size_t ns	= _sites[fi_mon]; // num sites for this type
      size_t nmon = typeCount.second; // num monomers for this type
      size_t nmon2 = nmon * 2;


      // Iterate through all non-excluded pairs
			for (auto& indexPair : siteIndexPairs) {
        auto i = indexPair.first;
        auto j = indexPair.second;

        // Get a1a2 and check if it's (floating pt) 0
			  // FIXME: if it is, then i'm assuming we need to screen? is that
			  // the only difference?
        double A = _polfac[fi_sites + i] * _polfac[fi_sites + j];
        bool screeningNeeded = (A > constants::EPS);

        size_t inmon  = i*nmon;
        size_t inmon3  = 3*inmon;
        size_t i3 = 3*i;
        size_t jnmon  = j*nmon;
        size_t jnmon3  = 3*jnmon;
        size_t j3 = 3*j;


        // Iterate through this pair (in each monomer)
        for (size_t m = 0; m < nmon; m++) {
          // Reset temporaries
          // FIXME: This is probably overkill, but I'm not 100% sure about the
          // safety of omitting it
          initTemps();


          //  Distances and values that will be reused
          const double rijx = _xyz[fi_crd + inmon3 + m]
                            - _xyz[fi_crd + jnmon3 + m];
          const double rijy = _xyz[fi_crd + inmon3 + m + nmon]
                            - _xyz[fi_crd + jnmon3 + m + nmon];
          const double rijz = _xyz[fi_crd + inmon3 + m + nmon2]
                            - _xyz[fi_crd + jnmon3 + m + nmon2];
          const double rsq = rijx*rijx + rijy*rijy + rijz*rijz;
          const double r = std::sqrt(rsq);
          const double ri = 1.0/r;
          const double risq = ri*ri;

				  // Update the field
      	  const size_t shift = _first_ind[fi_mon + m];
      	  const size_t spi = shift + i;
      	  const size_t spj = shift + j;

      	  _phii[m] = ri * _chg[spj];
      	  _phij[m] = ri * _chg[spi];


      	  // Update Electric field
      	  _Efqix[m] = ri * risq * _chg[spj] * rijx;
      	  _Efqjx[m] = ri * risq * _chg[spi] * rijx;
      	  _Efqiy[m] = ri * risq * _chg[spj] * rijy;
      	  _Efqjy[m] = ri * risq * _chg[spi] * rijy;
      	  _Efqiz[m] = ri * risq * _chg[spj] * rijz;
      	  _Efqjz[m] = ri * risq * _chg[spi] * rijz;

          // Multiply by screening coeffs. if needed
          if (screeningNeeded) {
            double s0r, s1r3;
            screenCoeff(rsq, A, s0r, s1r3);

      	    _phii[m] *= s0r;
      	    _phij[m] *= s0r;

      	    _Efqix[m] *= s1r3;
      	    _Efqjx[m] *= s1r3;
      	    _Efqiy[m] *= s1r3;
      	    _Efqjy[m] *= s1r3;
      	    _Efqiz[m] *= s1r3;
      	    _Efqjz[m] *= s1r3;
          }
        }

        // Almost there...
        for (size_t m = 0; m < nmon; m++) {
          const size_t shift = _first_ind[fi_mon + m];
          const size_t shift3 = 3*shift;
          _phi[i + shift] += _phii[m];
          _phi[j + shift] += _phij[m];
          _Efq[shift3 + i3] += _Efqix[m];
          _Efq[shift3 + j3] -= _Efqjx[m];
          _Efq[shift3 + i3 + 1] += _Efqiy[m];
          _Efq[shift3 + j3 + 1] -= _Efqjy[m];
          _Efq[shift3 + i3 + 2] += _Efqiz[m];
          _Efq[shift3 + j3 + 2] -= _Efqjz[m];
        }
      }


      // Increment first indices
      fi_mon += nmon, fi_sites += nmon * ns, fi_crd += nmon * ns * 3;
    }
  }


  double Electrostatics::eval() {
    // Reset temporaries
    initTemps();

    computePermanent();
    computeDipoles();
    recomputePermanent();

    // TODO

    return 0.0;
  }

} // namespace elec
