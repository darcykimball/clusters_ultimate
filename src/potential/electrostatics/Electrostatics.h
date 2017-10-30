#ifndef ELECTROSTATICS_H
#define ELECTROSTATICS_H

#include <vector>
#include <cstdlib>
#include <string>
#include <cmath>
#include <iostream>

#include "bblock/sys_tools.h"
#include "tools/definitions.h"
#include "tools/constants.h"
#include "potential/electrostatics/gammq.h"
#include "potential/electrostatics/aliases.h"


namespace elec {


  class Electrostatics {
    public:

      // Exceptions
      struct MaxIterExceeded : std::exception {
        const char* what() const noexcept override {
          return "Maximum iterations exceeded; dipoles failed to converge";
        }
      };

      struct DipolesDiverged : std::exception {
        const char* what() const noexcept override {
          return "Dipoles diverged";
        }
      };

      
      // XXX: Are these parameters always fixed? Anyways, factor them out if you want
      // to leave them as params
      static constexpr double aCC = 0.4;
      static constexpr double aDD = 0.055;

    private:
      // These are parameters of a given system/simulation instance, so I'm
      // assuming it makes sense to put them here, i.e there'll only be one,
      // unchanging Electrostatics instance per run.
      const Vec<double>& _chg;
      const Vec<double>& _polfac;
      const Vec<double>& _pol;
      const Vec<double>& _orig_xyz;
      const Vec<Id>& _mon_id;
      const Vec<size_t>& _sites;
      const Vec<size_t>& _first_ind;
      const AssocVec<Id, size_t>& _mon_type_count;
      const double _tolerance;
      const size_t _maxit;
      const bool _do_grads;
      Vec<double>& _grad;

      // Mapping from each monomer type/count to the set of index pairs that are
      // not excluded from site-site calculations
      const AssocVec<Pair<Id, size_t>, Vec<IndexPair>> _notExcluded;

      // These are temporaries used for each computation. Figured they could be
      // alloc'd here once and for all

      // These shouldn't change for a given instance
      const size_t _nsites;
      const Vec<double> _sqrtpol;

      // Electric fields and potential
      Vec<double> _Efq;
      Vec<double> _Efd;
      Vec<double> _phi;
      Vec<double> _mu;
      Vec<double> _xyz;

		  Vec<double> _phii;
		  Vec<double> _phij;
		  Vec<double> _Efqix;
		  Vec<double> _Efqjx;
		  Vec<double> _Efqiy;
		  Vec<double> _Efqjy;
		  Vec<double> _Efqiz;
		  Vec<double> _Efqjz;


      // Initialize/reset temporaries
      void initTemps();


			// Organize xyz so we have x1_1 x1_2 ... y1_1 y1_2...
			// where xN_M is read as coordinate x of site N of monomer M
	    // for the first monomer type. Then follows the second, and so on.
      void organizeCrds();


      // Figure out which pairs are (not) excluded
      AssocVec<Pair<Id, size_t>, Vec<IndexPair>> getNonExcluded() const;


      // Compute the permanent electric field (Efq and friends? FIXME)
      // Postcondition: the temp vectors that hold the components of the perm.
      // elec. field are initialized with, well, the field values.
      void computePermanent();


      // Compute dipole field (iteratively)
      void computeDipoles();


      // Compute intra-monomer site-site contributions
      void computeIntraSites();


      // Compute inter-monomer site-site contributions
      void computeInterSites();


      // Recompute permanent field
      // Precondition: computeDipoles() was called successfully
      void recomputePermanent();


      // Recompute intra-monomer site-site contributions
      // Precondition: computeDipoles() was called successfully
      void recomputeIntraSites();


      // Recompute inter-monomer site-site contributions
      // Precondition: computeDipoles() was called successfully
      void recomputeInterSites();


      // Check if convergence was achieved
      bool converged() {
        // TODO
        return false;
      }


      // XXX: Sigh...
      size_t maxMonomerCount() {
        return std::max_element(_mon_type_count.cbegin(), _mon_type_count.cend(),
            [](auto const& tc1, auto const& tc2) {
              return tc1.second < tc2.second;
            }
          )->second;
      }


			// XXX: This might be worth making public
			// Calculate screening coefficients (set as output params)
      static void screenCoeff(const double rsq, const double A,
				double& s0r, double& s1r3);

    public:

      Electrostatics(
        const Vec<double>& chg,
        const Vec<double>& polfac,
        const Vec<double>& pol,
        const Vec<double>& orig_xyz,
        const Vec<Id>& mon_id,
        const Vec<size_t>& sites,
        const Vec<size_t>& first_ind,
        const AssocVec<Id, size_t>& mon_type_count,
        const double tolerance, const size_t maxit, const bool do_grads,
        Vec<double>& grad);

      // XXX: if it's true that there really should only be a unique instance..
      Electrostatics(Electrostatics&&) = delete;
      Electrostatics(Electrostatics const&) = delete;
      Electrostatics& operator=(Electrostatics const&) = delete;
      Electrostatics& operator=(Electrostatics&&) = delete;


      // Evaluate the energy
      // Side-effects: grad vector is updated if do_grads is true
      double eval();

  };

} //namespace elec

#endif
