/*!
 * \file population_sequence_vonmises.hpp
 * \author Dillon Cislo, dilloncislo@gmail.com
 * \date 09-Apr-2022
 * \brief Declaration of PopulationSequenceVonMises class.
 */

#ifndef __POPULATION_SEQUENCE_VON_MISES__
#define __POPULATION_SEQUENCE_VON_MISES__

#include <list>
#include <stdexcept>
#include <cmath>

#include "population.hpp"
#include "rng_vm.hpp"

using std::list;
using std::runtime_error;
using std::atan2;
using std::cos;
using std::sin;

/*! PopulationSequenceVonMises class handles a basic procedure in which cells
 * divide according to a prescribed sequence. Cells will divide in the order specified
 * by their particle ID. Once a cell divides, the progeny are moved to the end of the
 * queue. After a division, the director field of the progeny may either be
 * inherited from the parent or drawn at random from a von Mises distribution.
 *
 * Particles are always divided along the direction of the orientation vector n
 *
 */
class PopulationSequenceVonMises : public Population
{

public:

  //! Construct PopulationSequenceVonMises object
  //! \param sys Reference to the System object
  //! \param msg Reference to the system wide messenger
  //! \param param Reference to the parameters for population control
  PopulationSequenceVonMises(SystemPtr sys, const MessengerPtr msg, pairs_type& param) :
    Population(sys, msg, param)
  {

    if (param.find("seed") == param.end())
    {
      m_msg->msg(Messenger::WARNING,
          "Cell population control. No random number generator seed specified. Using default 0.");
      m_rng_vm = make_shared<RNG_VM>(0);
      m_msg->write_config("population.sequence_von_mises.seed", lexical_cast<string>(0));
    }
    else
    {
      m_msg->msg(Messenger::INFO,
          "Cell population control. Setting random number generator seed to "+param["seed"]+".");
      m_rng_vm = make_shared<RNG_VM>(lexical_cast<int>(param["seed"]));
      m_msg->write_config("population.sequence_von_mises.seed", param["seed"]);
    }

    if (param.find("mu") == param.end())
    {
      m_msg->msg(Messenger::WARNING, "Cell population control. No circular mean specified. Using default 0.0");
      m_mu = 0.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO, "Cell population control. Setting circular mean "+param["mu"]+".");
      m_mu = lexical_cast<double>(param["mu"]);
      m_mu = atan2(sin(m_mu), cos(m_mu));
    }
    m_msg->write_config("population.sequence_von_mises.mu", lexical_cast<string>(m_mu));

    if (param.find("k") == param.end())
    {
      m_msg->msg(Messenger::WARNING, "Cell population control. No distribution concentration specified. Usuing default 1.0");
      m_k = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO, "Cell population control. Setting distribution concentration "+param["k"]+".");
      m_k = lexical_cast<double>(param["k"]);
      if (m_k < 0.0)
      {
        m_msg->msg(Messenger::ERROR, "Cell population control. Distribution concentration must be nonnegative.");
        throw runtime_error("Wrong distribution concentration.");
      }
    }
    m_msg->write_config("population.sequence_von_mises.k", lexical_cast<string>(m_k));

    if (param.find("n") == param.end())
    {
      m_msg->msg(Messenger::WARNING, "Cell population control. No distribution peak number specified. Using default 1");
      m_n = 1;
    }
    else
    {
      m_msg->msg(Messenger::INFO, "Cell population control. Setting distribution peak number "+param["n"]+".");
      m_n = lexical_cast<int>(param["n"]);
      if (m_n <= 0)
      {
        m_msg->msg(Messenger::ERROR, "Cell population control. Distribution peak number must be greater than 0.");
        throw runtime_error("Wrong distribution peak number.");
      }
    }
    m_msg->write_config("population.sequence_von_mises.n", lexical_cast<string>(m_n));

    if (param.find("truncate_domain") == param.end())
    {
      m_msg->msg(Messenger::WARNING, "Cell population control. Distribution domain will be truncated.");
      m_truncate_domain = true;
      m_msg->write_config("population.sequence_von_mises.truncate_domain", "true");
    }
    else
    {
      m_msg->msg(Messenger::INFO, "Cell population control. Distribution domain will not be truncated.");
      m_truncate_domain = false;
      m_msg->write_config("population.sequence_von_mises.truncate_domain", "false");
    }

    if (param.find("inherit_director") == param.end())
    {
      m_msg->msg(Messenger::WARNING, "Cell population control. Divsion children will inherit their parent's orientation vector.");
      m_inherit_director = true;
      m_msg->write_config("population.sequence_von_mises.inherit_director", "true");
    }
    else
    {
      m_msg->msg(Messenger::INFO, "Cell population control. Division children will not inherit their parent's orientation vector.");
      m_inherit_director = false;
      m_msg->write_config("population.sequence_von_mises.inherit_director", "false");
    }

    if (param.find("split_distance") == param.end())
    {
      m_msg->msg(Messenger::WARNING, "Cell population control. No split distance set. Using default 0.5.");
      m_split_distance = 0.5;
    }
    else
    {
      m_msg->msg(Messenger::INFO, "Cell population control. Split distance set to "+param["split_distance"]+".");
      m_split_distance = lexical_cast<double>(param["split_distance"]);
    }
    m_msg->write_config("population.sequence_von_mises.split_distance", lexical_cast<string>(m_split_distance));

    if (param.find("move_split") == param.end())
    {
      m_msg->msg(Messenger::WARNING, "Cell population control. No split of the offspring separation set. Assuming 0.5");
      m_alpha = 0.5;
    }
    else
    {
      m_msg->msg(Messenger::INFO, "Cell population control. Split of the offspring separation set to "+param["move_split"]+".");
      m_alpha = lexical_cast<double>(param["move_split"]);
    }
    m_msg->write_config("population.sequence_von_mises.move_split", lexical_cast<string>(m_alpha));

    if (param.find("death_rate") == param.end())
    {
      m_msg->msg(Messenger::WARNING, "Cell population control. No death rate set. Using default 0.0.");
      m_death_rate = 0.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO, "Cell population control. Setting death rate "+param["death_rate"]+".");
      m_death_rate = lexical_cast<double>(param["death_rate"]);
    }
    m_msg->write_config("population.sequence_von_mises.death_rate", lexical_cast<string>(m_death_rate));

    // Initialize division sequence
    // Division sequence is set by particle ID order
    vector<int> particles = m_system->get_group(m_group_name)->get_particles();
    int N = m_system->get_group(m_group_name)->get_size();

    m_div_queue.reserve(N);
    for (int i = 0; i < particles.size(); i++)
    {
      int pi = particles[i];
      Particle &p = m_system->get_particle(pi);
      if (p.in_tissue && !p.boundary)
        m_div_queue.push_back(p);
    }

  };

  //! Particle division (emulated cell division)
  void divide(int);

  //! Remove particle (emulates cell death)
  void remove(int time) { };

  //! Add particle (has no direct biological application)
  void add(int time) { };

  //! Change cell native area
  void grow(int time) { };

  //! Change particle length (makes no sense here)
  void elongate(int time) { };

private:

  RNGVMPtr m_rng_vm;             //!< Cell number/division angle generator
  double m_mu;                   //!< Mean of von Mises distribution
  double m_k;                    //!< Concentration of von Mises distribution
  int m_n;                       //!< Number of peaks in the distribution
  bool m_truncate_domain;        //!< If true, truncate basic distribution domain

  vector<Particle> m_div_queue;  //!< Queue of cell division events
  bool m_inherit_director;       //!< If true, division children inherit their parent's orientation vector
  double m_split_distance;       //!< Fraction of the particle radius to split after the division
  double m_alpha;                //!< When dividing particles, move new one to alpha*m_split_distance and the old one to (1-alpha)*split_distance

  double m_death_rate;           //!< Rate of cell death

};

typedef shared_ptr<PopulationSequenceVonMises> PopulationSequenceVonMisesPtr;

#endif
