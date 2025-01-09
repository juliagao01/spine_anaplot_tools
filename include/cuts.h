/**
 * @file cuts.h
 * @brief Header file for definitions of selection cuts.
 * @author justin.mueller@colostate.edu 
*/

#include <functional>
#include <vector>
#include <TVector3.h>
#include <string>
#include <sstream>
#include <numeric>
#include <iostream>

#ifndef CUTS_H
#define CUTS_H

namespace cuts
{
    /**
     * Apply a cut on whether a match exists.
     * @tparam T the object type (true or reco, interaction or particle).
     * @param obj the object to select on.
     * @return true if the object is matched.
    */
    template<class T>
        bool matched(const T & obj) { return obj.match_ids.size() > 0; }

    /**
     * Apply a cut on the validity of the flash match.
     * @tparam T the type of interaction (true or reco).
     * @param interaction on which to place the flash validity cut.
     * @return true if the interaction is flash matched and the time is valid.
    */
    template<class T>
        bool valid_flashmatch(const T & interaction) { return !std::isnan(interaction.flash_time) && interaction.is_flash_matched == 1; }

    /**
     * Check if the particle meets final state signal requirements.
     * Particles must be primary and have an energy above threshold.
     * Muons must have a length of at least 50 cm (143.425 MeV), protons
     * must have an energy above 50 MeV, and all other particles must have
     * an energy above 25 MeV.
     * @tparam T the type of particle (true or reco).
     * @param particle to check.
     * @return true if the particle is a final state signal particle.
    */
    template<class T>
        bool final_state_signal(const T & p)
        {
            bool passes(false);
            if(p.is_primary)
            {
                double energy(p.pid > 1 ? p.csda_ke : p.calo_ke);
                if constexpr (std::is_same_v<T, caf::SRParticleTruthDLPProxy>)
                    energy = p.energy_deposit;

                if((p.pid == 2 && energy > 143.425) || (p.pid != 2 && p.pid < 4 && energy > 25) || (p.pid == 4 && energy > 50))
                    passes = true;
            }
            return passes;
        }

    /**
     * Count the primaries of the interaction with cuts applied to each particle.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to find the topology of.
     * @return the count of primaries of each particle type within the
     * interaction.
     */
    template<class T>
        std::vector<uint32_t> count_primaries(const T & interaction)
        {
            std::vector<uint32_t> counts(5, 0);
            for(auto &p : interaction.particles)
            {
                if(final_state_signal(p))
                    ++counts[p.pid];
            }
            return counts;
        }

    /**
     * Find the topology of the interaction with cuts applied to each particle.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to find the topology of.
     * @return the topology of the interaction as a string (e.g 0ph0e1mu0pi1p).
     */
    template<class T>
        std::string topology(const T & interaction)
        {
            std::vector<uint32_t> counts(count_primaries(interaction));
            std::stringstream ss;
            ss  << counts[0] << "ph"
                << counts[1] << "e"
                << counts[2] << "mu"
                << counts[3] << "pi"
                << counts[4] << "p";
            return ss.str();
        }

    /**
     * Apply selection for 1muNp topology.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction has 1muNptopology.
     */
    template <class T> bool topological_1muNp_cut(const T & interaction)
      {
	std::vector<uint32_t> c(count_primaries(interaction));
	return c[0] == 0 && c[1] == 0 && c[2] == 1 && c[3] == 0 && c[4] >= 1 && c[5] == 0;
      }

    /**
     * Apply a fiducial volume cut. Interaction vertex must be within 25 cm of
     * x and y detector faces, 50 cm of downstream (+) z face, and 30 cm of
     * upstream (-) z face.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the vertex is in the fiducial volume.
     */
    template<class T>
        bool fiducial_cut(const T & interaction)
        {
            return interaction.is_fiducial && !(interaction.vertex[0] > 210.215 && interaction.vertex[1] > 60 && (interaction.vertex[2] > 290 && interaction.vertex[2] < 390));
        }
    
    /**
     * Apply a containment volume cut. All points within the interaction must be
     * at least 5 cm from the detector boundaries.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction is contained.
     */
    template<class T>
        bool containment_cut(const T & interaction) { return interaction.is_contained; }
    
    /**
     * Apply a flash time cut.  The interaction must be matched to an in-time
     * flash.  The in-time definition is valid for BNB simulation.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if interaction has been matched to an in-time flash.
     */
    template<class T> bool flash_cut_bnb(const T & interaction)
      {
	if(!valid_flashmatch(interaction))
          {
            return false;
          }
        else
	  {
	    return (interaction.flash_time >= 0) && (interaction.flash_time <= 1.6);
	  }
      }

    /**
     * Apply a flash time cut.  The interaction must be matched to an in-time
     * flash.  The in-time definition is valid for NuMI simulation.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction has been matched to an in-time flash.
     */
    template<class T> bool flash_cut_numi(const T & interaction)
      {
	if(!valid_flashmatch(interaction))
	  {
	    return false;
	  }
	else
	  {
	  return (interaction.flash_time >= 0) && (interaction.flash_time <= 9.6);
	  }
      }

    /**
     * Apply fiducial, track containment, topological (1mu + 2gamma + 0pi), and flahs time (NuMI) cut
     * (logical "and" of each).
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction passes above cuts.
     */
    template<class T> bool all_1muNp_cut(const T & interaction)
      {
	return topological_1muNp_cut<T>(interaction) && fiducial_cut<T>(interaction) && containment_cut<T>(interaction) && flash_cut_numi<T>(interaction);
      }
    
    /**
     * Defined the true neutrino interaction classification.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction is a neutrino interaction.
    */
    template<class T>
        bool neutrino(const T & interaction) { return interaction.nu_id > -1; }

    /**
     * Define the true cosmic interaction classification.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction is a cosmic.
     */
    template<class T>
        bool cosmic(const T & interaction) { return interaction.nu_id == -1; }

    /**
     * Define the true neutrino interaction classification.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction is a neutrino interaction.
    */
    template<class T>
        bool matched_neutrino(const T & interaction) { return interaction.match_ids.size() > 0 && neutrino(interaction); }

    /**
     * Define the true neutrino interaction classification.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction is a neutrino interaction.
    */
    template<class T>
        bool matched_cosmic(const T & interaction) { return interaction.match_ids.size() > 0 && cosmic(interaction); }

    /**
     * Define the true 1mu1pi0 interaction classification.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction is a 1mu1p neutrino interaction.
     */
    template<class T> bool signal_1muNp(const T & interaction)
      {
	return topological_1muNp_cut<T>(interaction) && neutrino(interaction);
      }
    
    template<class T> bool other_nu_1muNp(const T & interaction)
      {
	return !topological_1muNp_cut<T>(interaction) && neutrino(interaction);
      }
}
#endif
