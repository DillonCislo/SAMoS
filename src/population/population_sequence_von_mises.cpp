/*!
 * \file population_sequence_von_mises.cpp
 * \author Dillon Cislo, dilloncislo@gmail.com
 * \date 10-Apr-2022
 * \brief Implementation of PopulationSequenceVonMises class
 */

#include "population_sequence_von_mises.hpp"

/*! Divide cells according to their place in a predetermined sequence
 * Division order is set by the initial ordering of particle ideas supplied
 * by the user. Cells are popped off the front of the queue, divided, and their
 * offspring are then placed at the back of the queue.
 *
 * After the division, each daugther cell has the same native area as
 * the mother cell
 *
 * Position of daughter cells is determined by the orientation vector n
 * \param t current time step
 *
 */
void PopulationSequenceVonMises::divide(int t)
{

  if (m_freq > 0 && t % m_freq == 0) // Attempt division only at certain time steps
  {

    if (!m_system->group_ok(m_group_name))
    {
      cout << "Before cell division: Group info mismatch for group : " << m_group_name << endl;
      throw runtime_error("Group mismatch");
    }

    // Pop the next dividing cell from the queue
    Particle &p = m_div_queue[0];

    if ( !(p.in_tissue && !p.boundary) )
    {
      cout << "Before cell division: Bad cell picked to divide" << endl;
      throw runtime_error("Bad cell picked to divide");
    }

    // Generate the new cell
    Particle p_new(m_system->size(), p.get_type(), p.get_radius());
    p_new.x = p.x + m_alpha*m_split_distance*p.get_radius()*p.nx;
    p_new.y = p.y + m_alpha*m_split_distance*p.get_radius()*p.ny;
    p_new.z = p.z + m_alpha*m_split_distance*p.get_radius()*p.nz;
    p_new.set_parent(p.get_flag());
    m_system->apply_periodic(p_new.x,p_new.y,p_new.z);
    
    p.x -= (1.0-m_alpha)*m_split_distance*p.get_radius()*p.nx;
    p.y -= (1.0-m_alpha)*m_split_distance*p.get_radius()*p.ny;
    p.z -= (1.0-m_alpha)*m_split_distance*p.get_radius()*p.nz;
    p.age = 0.0;
    p.A0 = p.get_A0();
    m_system->apply_periodic(p.x,p.y,p.z);

    p_new.vx = p.vx; p_new.vy = p.vy; p_new.vz = p.vz;
    p_new.Nx = p.Nx; p_new.Ny = p.Ny; p_new.Nz = p.Nz;
    p_new.age = 0.0;
    p_new.A0 = p.get_A0();
    p_new.set_radius(p.get_radius());
    p_new.set_type(p.get_type());
    p_new.set_default_area(p.get_A0());
    p_new.in_tissue = true;

    if (m_inherit_director)
    {
      p_new.nx = p.nx; p_new.ny = p.ny; p_new.nz = p.nz;
    }
    else
    {
      double nL = std::sqrt( p.nx*p.nx + p.ny*p.ny + p.nz*p.nz );
      double phi1 = m_rng_vm->von_mises_rng(m_mu, m_k, m_n, m_truncate_domain);
      double phi2 = m_rng_vm->von_mises_rng(m_mu, m_k, m_n, m_truncate_domain);

      // NOTE: ONLY VALID WITH PLANE CONSTRAINTS
      p.nx = nL * std::cos(phi1); p.ny = nL * std::sin(phi1); p.nz = 0.0;
      p_new.nx = nL * std::cos(phi2); p_new.ny = nL * std::sin(phi2); p_new.nz = 0.0;
    }

    for(list<string>::iterator it_g = p.groups.begin(); it_g != p.groups.end(); it_g++)
      p_new.add_group(*it_g);

    m_system->add_particle(p_new);

    m_div_queue.push_back(p);
    m_div_queue.push_back(p_new);
    m_div_queue.erase(m_div_queue.begin());

    if (!m_system->group_ok(m_group_name))
    {
      cout << "After cell division: Group info mismatch for group : " << m_group_name << endl;
      throw runtime_error("Group mismatch.");
    }

    m_system->set_force_nlist_rebuild(true);

  }

}
