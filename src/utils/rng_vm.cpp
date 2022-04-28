/*!
 * \file rng_vm.cpp
 * \author Dillon Cislo, dilloncislo@gmail.com
 * \date 08-Apr-2022
 * \brief Class RNG_VM provides wrappers for the GSL random number generators
 */

#include "rng_vm.hpp"


//! Evaluate a von Mises distribution
//! \param x value at which the function is evaluated
//! \param mu mean of the circular distribution
//! \param k concentration of the distribution
//! \param n number of peaks in the distribution
//! \return random angle drawn from a von Mises distribution
double RNG_VM::von_mises(double x, double mu, double k, int n)
{

  double vm = std::exp( k * std::cos(((double) n)*(x-mu)) ) /
    ( 2.0 * M_PI * std::cyl_bessel_i(0.0, k) );

  return vm;

};

//! Return a von Mises distributed angle using rejection sampling
//! \param mu the mean of the circular distribution
//! \param k the concentration of the distribution
//! \param n the number of peaks in the distribution
//! \param truncateDomain whether or not to truncate the basic domain
//! \return a von Mises distributed angle on [-pi, pi]
double RNG_VM::von_mises_rng(double mu, double k, int n, bool truncateDomain)
{

  // Wrap mean to the domain [-pi, pi]
  mu = std::atan2(std::sin(mu), std::cos(mu));

  double vm;
  bool sampleFound = false;

  while ( !sampleFound ) {

    // Draw a number uniformly on the domain [-pi, pi]
    double curX = 2.0 * M_PI * this->drnd() - M_PI;

    // Draw a test sample from the proposed (uniform) distribution
    double curY = this->drnd();

    if (curY < this->von_mises(curX, mu, k, n)) {
      vm = curX;
      sampleFound = true;
      break;
    }

  }

  // Wrap the output to the appropriate domain
  if ( truncateDomain ) {

    // Shift the angles to be centered about the mean
    vm = vm - mu;

    // Divide the domain by the number of peaks
    // The domain should still be symmetric around zero
    double dmin = -M_PI / ((double) n);
    double dmax = M_PI / ((double) n);
    double D = dmax - dmin;

    // Wrap the angles to [-D/2, D/2]
    if ((vm < dmin) || (dmax < vm)) {

      vm = vm + (D / 2.0);
      bool positiveInput = vm > 0;
      
      vm = this->mmod(vm, D);
      if ((vm == 0) && positiveInput)
        vm = D;

       vm = vm - D/2;

    }

    // Shift the angles back
    vm = vm + mu;

  }

  // Wrap angle to the domain [-pi, pi]
  vm = std::atan2(std::sin(vm), std::cos(vm));

  return vm;

};

//! Remainder after division (modulo operation - matches MATLAB's mod(q,m))
//! \param q dividend
//! \param m divisor
//! \return Remainder after division
double RNG_VM::mmod(double q, double m)
{

  if (m == 0.0)
    return q;

  double result = std::fmod(q, m);
  result = ((result >= 0 && m > 0) || (q <= 0 && m < 0)) ? result : (result + m);
  return result;

};
