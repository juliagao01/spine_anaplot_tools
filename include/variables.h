/**
 * @file variables.h
 * @brief Header file for definitions of selection variables.
 * @author justin.mueller@colostate.edu
*/
#ifndef VARIABLES_H
#define VARIABLES_H

#define ELECTRON_MASS 0.5109989461
#define MUON_MASS 105.6583745
#define PION_MASS 139.57039
#define PROTON_MASS 938.2720813

#include <algorithm>
#include <iostream>
#include <TVector3.h>
#include <string>
#include <iostream>

namespace vars
{

    /**
     * Variable for counting interactions/particles.
     * @tparam T the type of object (true or reco, interaction or particle).
     * @param interaction/particle to apply the variable on.
     * @return 1.0 (always).
     */
    template<class T>
        double count(const T & obj) { return 1.0; }
   
    /**
     * Variable for id (unique identifier for the object).
     * @tparam T the type of object (true or reco, interaction or particle).
     * @param interaction/particle to apply the variable on.
     * @return the id of the interaction/particle.
    */
    template<class T>
        double id(const T & obj) { return obj.id; }

    /**
     * Variable for enumerating interaction categories.  This is a basic
     * categorization using only signal, neutrino background, and cosmic
     * background as the three categories.
     * 0: 1muNp (contained and fiducial)
     * 1: 1muNp (not contained or fiducial)               
     * 2: Other nu                
     * 3: cosmic
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the enumerated category of the interaction.   
     */
    template<class T>
      double category(const T & interaction)
      {
        // Cosmic background               
        double cat(3);

        // Signal
        if(cuts::signal_1muNp(interaction))
          {
	    if(cuts::fiducial_cut(interaction) && cuts::containment_cut(interaction))
	      {
		cat = 0;
	      }
	    else cat = 1;
          }
        // Neutrino Background                         
        else if(cuts::other_nu_1muNp(interaction))
          {
            cat = 2;
          }
        return cat;
      }

    /**
     * Variable for enumerating interaction categories. This classifies the          
     * interactions based on the visible final states.        
     * @tparam T the type of interaction (true or reco).             
     * @param interaction to apply the variable on.      
     * @return the enumerated category of the interaction.                                               
     */
    template<class T>
        double category_topology(const T & interaction)
        {
            uint16_t cat(6);
            if(interaction.nu_id >= 0)
            {
                std::vector<uint32_t> counts(cuts::count_primaries(interaction));
                if(counts[0] == 0 && counts[1] == 0 && counts[2] == 1)
                {
                    if(counts[3] == 0 && counts[4] == 1 && interaction.is_contained && interaction.is_fiducial) cat = 2;
                    else if(counts[3] == 0 && counts[4] == 1) cat = 7;
                    else if(counts[3] == 0 && counts[4] == 0) cat = 4;
                    else if(counts[3] == 0 && counts[4] > 1 && interaction.is_contained && interaction.is_fiducial) cat = 2;
                    else if(counts[3] == 0 && counts[4] > 1) cat = 7;
                    else if(counts[3] == 1 && counts[4] == 1) cat = 4;
                    else if(interaction.current_type == 0) cat = 4;
                }
                else if(interaction.current_type == 0) cat = 4;
                else if(interaction.current_type == 1) cat = 5;
            }
            return cat;
        }

    /**
     * Variable for enumerating interaction categories. This categorization
     * uses the interaction type (generator truth) classify the interactions
     * 0: nu_mu CC QE, 1: nu_mu CC Res, 2: nu_mu CC MEC, 3: nu_mu CC DIS, 4: nu_mu CC Coh, 5: nu_e CC, 6: NC, 7: Cosmic
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the enumerated category of the interaction.
    */
    template<class T>
        double category_interaction_mode(const T & interaction)
        {
            double cat(7);

            if(interaction.nu_id > -1)
            {
                if(interaction.current_type == 0)
                {
                    if(abs(interaction.pdg_code) == 14)
                    {
                        if(interaction.interaction_mode == 0) cat = 0;
                        else if(interaction.interaction_mode == 1) cat = 1;
                        else if(interaction.interaction_mode == 10) cat = 2;
                        else if(interaction.interaction_mode == 2) cat = 3;
                        else if(interaction.interaction_mode == 3) cat = 4;
                        else cat = 8;
                    }
                    else cat = 5;
                }
                else cat = 6;
            }

            return cat;
        }

    /**
     * Variable for total visible energy of interaction.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the total visible energy of the interaction.
     */
    template<class T>
      double visible_energy(const T & interaction)
      {
	double energy(0);
	for(const auto & p : interaction.particles)
	  {
	    if(p.is_primary)
	      {
		if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
			       {
				 energy += p.energy_deposit;
			       }
		else
		  {
		    if(p.pid < 2) energy += p.calo_ke;
		    else energy += p.csda_ke;
		  }
		if(p.pid == 2) energy += MUON_MASS;
		else if(p.pid == 3) energy += PION_MASS;
	      }
	  }
	return energy;
      }


}

#endif
