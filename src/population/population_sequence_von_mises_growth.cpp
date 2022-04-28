/*!
 * \file population_sequence_von_mises_growth.cpp
 * \author Dillon Cislo, dilloncislo@gmail.com
 * \date 28-Apr-2022
 * \brief Implementation of PopulationSequenceVonMisesGrowth class
 */

#include "population_sequence_von_mises_growth.hpp"

/*! Divide cells according to their place in a predetermined sequence
 * Division order is set by the initial ordering of particle ideas supplied
 * by the user. Cells are programmed to grow their native area with a delay
 * specified by the ordering. Cells divide once their native area crosses a
 * particular threshold.
 *
 * After the division, each daugther cell has the same native area as
 * the mother cell
 *
 * Position of daughter cells is determined by the orientation vector n
 * \param t current time step
 *
 */
void PopulationSequenceVonMisesGrowth::divide(int t)
{

  if (m_freq > 0 && t % m_freq == 0) // Attempt division only at certain time steps
  {

    if (!m_system->group_ok(m_group_name))
    {
      cout << "Before cell division: Group info mismatch for group : " << m_group_name << endl;
      throw runtime_error("Group mismatch");
    }

    int N = m_system->get_group(m_group_name)->get_size();
    vector<int> particles = m_system->get_group(m_group_name)->get_particles();

    for (int i = 0; i < N; i++ )
    {

      int pi = particles[i];
      Particle &p = m_system->get_particle(pi);

      if (p.in_tissue && !p.boundary && p.A0 > m_max_A0)
      {

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

      }

    }

    if (!m_system->group_ok(m_group_name))
    {
      cout << "After cell division: Group info mismatch for group : " << m_group_name << endl;
      throw runtime_error("Group mismatch.");
    }

    m_system->set_force_nlist_rebuild(true);

  }

};

/*! Remove particle
 *
 * This is essentially a clean-up operation to preserve tissue integrity.
 * An internal cell is removed if it has fewer than a minimum threshold of
 * internal neighbors
 *
 * \param t current time step
 *
 */
void PopulationSequenceVonMisesGrowth::remove(int t)
{
  if (m_freq > 0 && t % m_freq == 0 && m_death_rate > 0.0 && m_min_nn > 0)
  {

    if (!m_system->group_ok(m_group_name))
    {
      cout << "Before Cell Remove: Group info mismatch for group : " << m_group_name << endl;
      throw runtime_error("Group mismatch.");
    }

    int N = m_system->get_group(m_group_name)->get_size();
    vector<int> particles = m_system->get_group(m_group_name)->get_particles();
    vector<int> to_remove;

    // Determine which cells have fewer than the minimum number of neighbors
    for (int i = 0; i < N, i++)
    {

      int pi = particles[i];
      Particle &p = m_system->get_particle(pi);

      int nn = 0;
      for (int j = 0; j < m_nlist->m_contact_list[pi].size(); j++)
      {
        int pj = m_nlist->m_contact_list[pi][j];
        Particle &pn = m_system->get_particle(pj);
        if (pn.in_tissue && !pn.boundary)
          nn++;
      }

      if (p.in_tissue && !p.boundary && nn < m_min_nn)
        to_remove.push_back(p.get_id());

    }

    // Remove the invalid cells from the system
    int offset = 0;
    for (vector<int>::iterator it = to_remove.begin(); it != to_remove.end(); it++)
    {
      m_system->remove_particle((*it)-offset);
      offset++;      
    }

    if (m_system->size() == 0)
    {
      m_msg->msg(Messenger::ERROR,"Cell population control. No cells left in the system. Please reduce that death rate.");
      throw runtime_error("No cells left in the system.");
    }

    if (!m_system->group_ok(m_group_name))
    {
      cout << "After Cell Remove: Group info mismatch for group : " << m_group_name << endl;
      throw runtime_error("Group mismatch.");
    }

    m_system->set_force_nlist_rebuild(true);

  }
};

/*! Grow (rescale) cell native area by a given amount.
 *
 *  \param t current time step
 *
 */
void PopulationSequenceVonMisesGrowth::grow(int t)
{
  if (m_freq > 0 && t % m_freq == 0 && m_growth_rate > 0.0)
  {

    double fact = m_freq*m_system->get_integrator_step()*m_growth_rate;
    int N = m_system->get_group(m_group_name)->get_size();
    vector<int> particles = m_system->get_group(m_group_name)->get_particles();

    // Implement delay by not allowing cells in the queue to grow
    if (!m_div_queue.empty())
    {

      // Remove the next cell in the queue allowing it to begin to grow
      if (m_growth_delay_freq > 0 && t % m_growth_delay_freq == 0)
        m_div_queue.pop_back();

      for (int i = 0; i < m_div_queue.size(); i++)
      {

        Particle &p = m_div_queue[i];
        int pi = p.get_id();
        vector<int>::iterator position = find(particles.begin(), particles.end(), pi);
        if (position != particles.end()) // == particles.end() means it was not found
          particles.erase(position);
      }

    }

    for (int i = 0; i < particles.size(); i++)
    {
      if (m_rng_vm->drnd() < m_growth_prob)
      {
        int pi = particles[i];
        Particle &p = m_system->get_particle(pi);
        if (p.in_tissue)
          p.A0 *= (1.0+fact);
      }
    }

    m_system->set_force_nlist_rebuild(true);
    if (m_rescale_contacts)
      m_system->set_nlist_rescale(sqrt(1.0+fact));

  }

};
