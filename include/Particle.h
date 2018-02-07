#ifndef PARTICLE_H
#define PARTICLE_H

#include <string>
#include <vector>
#include "TVector3.h"

namespace selection{
  
  /**
   * @brief  Particle class
   */
  class Particle{

    public : 

      /**
       * @brief  Constructor for MC particles 
       *
       * @param  pdg of the particle
       * @param  mass mass of the particle
       * @param  energy total energy of the particle
       * @param  vertex start point of the track
       * @param  end end of the track
       * @param  momentum momentum of the track
       *
       */
      Particle(const int pdg, const float mass, const float energy, const TVector3 &vertex, const TVector3 &end, const TVector3 &momentum);

      /**
       * @brief  Constructor for reconstructed tracks 
       *
       * @param  pdg of the particle
       * @param  kinetic_energy total energy of the particle
       * @param  length length of the particle
       * @param  vertex start point of the track
       * @param  end end of the track
       *
       */
      Particle(const int pdg, const float kinetic_energy, const float length, const TVector3 &vertex, const TVector3 &end);

      /**
       * @brief  Constructor for reconstructed showers 
       *
       * @param  pdg of the particle
       * @param  vertex start point of the track
       * @param  end end of the track
       *
       */
      Particle(const int pdg, const TVector3 &vertex, const TVector3 &end);

      /**
       * @brief  Get the mass from the pdg code
       *
       * @param  pdg pdg of the particle
       *
       */
      float GetMassFromPdg(const int pdg) const;

      /**
       * @brief  Get the pdg code
       */
      int GetPdgCode() const;

      /**
       * @brief  Get the mass
       */
      float GetMass() const;

      /**
       * @brief  Get the energy
       */
      float GetEnergy() const;

      /**
       * @brief  Get the length
       */
      float GetLength() const;

      /**
       * @brief  Get the vertex
       */
      TVector3 GetVertex() const;

      /**
       * @brief  Get the end
       */
      TVector3 GetEnd() const;

      /**
       * @brief  Get the momentum
       */
      TVector3 GetMomentum() const;

      /**
       * @brief  Get whether the particle has calorimetry
       */
      bool GetHasCalorimetry() const;


    private : 

      int      m_pdg;             ///< pdg code
      float    m_mass;            ///< mass of the particle
      float    m_energy;          ///< energy of the particle
      float    m_length;          ///< energy of the particle
      bool     m_has_calorimetry; ///< whether or not the particle has calorimetry
      TVector3 m_vertex;          ///< particle start position
      TVector3 m_end;             ///< particle end position
      TVector3 m_momentum;        ///< particle momentum


  }; // Particle
} // Selection
#endif
