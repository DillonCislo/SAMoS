/*!
 * \file rng_vm.hpp
 * \author Dillon Cislo, dilloncislo@gmail.com
 * \date 08-Apr-2022
 * \brief Class RNG_VM provides wrappers for the GSL random number generators.
 */

#ifndef __RNG_VM_HPP__
#define __RNG_VM_HPP__

#include <memory>
#include <math.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "rng.hpp"

using std::shared_ptr;

/*! Class handles random numbers in the system.
 * Has additional handling for drawing angles from a
 * von Mises distribution
 */
class RNG_VM : public RNG
{

  public:

    //! Constructor (initialize random number generator)
    RNG_VM(int seed) : RNG(seed) {};

    //! Destructor
    ~RNG_VM() {};

    //! Evaluate a von Mises distribution
    //! \param x the value at which the function is evaluated
    //! \param mu the mean of the circular distribution
    //! \param k the concentration of the distribution
    //! \param n the number of peaks in the distribution
    //! \return the random angle drawn from a von Mises distribution
    double von_mises(double x, double mu, double k, int n);

    //! Return a von Mises distributed angle using rejection sampling
    //! \param mu the mean of the circular distribution
    //! \param k the concentration of the distribution
    //! \param n the number of peaks in the distribution
    //! \param truncateDomain whether or not to truncate the basic domain
    //! \return a von Mises distributed angle on [-pi, pi]
    double von_mises_rng(double mu, double k, int n, bool truncateDomain);

    //! Remainder after division (modulo operation - matches MATLAB's mod(q,m))
    //! \param q dividend
    //! \param m divisor
    //! \return Remainder after division
    double mmod(double q, double m);

};

typedef shared_ptr<RNG_VM> RNGVMPtr;

#endif
