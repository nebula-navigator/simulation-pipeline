//
// Created by lex on 2/27/23.
//

#include "IO.h"
#include "tsunami.hpp"
#include "simprof.hpp"
#include <filesystem>

namespace fs = std::filesystem;
namespace Nb = Nbodyalgorithms;
using std::cout, std::endl;

int main(int argc, char *argv[]) {
    try { // Catch errors

        TsunamiCode Tsunami; // Default constructor, no parenthesis :(

        double dt_output = 6.28372892847e-1; // output step
        double tfin = 6.28372892847e2; // Final time

        if (fs::exists(fs::path(Tsunami.paramfile))) {
            IO::read_parameter_file(Tsunami.paramfile, Tsunami.N, tfin, Tsunami.Mscale,
                                    Tsunami.Lscale, dt_output, Tsunami.Conf, Tsunami.Conf.dcoll);
        } else {
            cout << Tsunami.paramfile + " not found. Using default values" << endl;
        }

        bool cont = false;
        std::vector<int> pthold;

        // Read command line
        IO::read_command_line(argc, argv, tfin, dt_output, Tsunami.Lscale,
                              Tsunami.Mscale, cont, Tsunami.Conf, Tsunami.fname_in,
                              Tsunami.N, pthold);

        Tsunami.allocate_arrays();
        Tsunami.set_units(Tsunami.Mscale, Tsunami.Lscale);

        IO::get_output_file_names(Tsunami.fname_in, Tsunami.ene_name, Tsunami.out_name, Tsunami.coll_name);

        if(cont) {
            cont = IO::read_last_timestep(Tsunami.N, Tsunami.System.pos, Tsunami.System.vel,
                                          Tsunami.System.mass, Tsunami.System.radius,
                                          Tsunami.System.xdata, Tsunami.System.spin,
                                          Tsunami.out_name, Tsunami.fname_in);
            remove(Tsunami.coll_name.c_str()); // Removes collision file in case it was generated

            if(cont) { // Continue only if there is a valid output file, otherwise will start from initial conditions

                // Read t0 e0 and eoff from last energy.dat output
                IO::get_restart_infos(Tsunami.time, Tsunami.energy, Tsunami.eoff, Tsunami.ene_name,
                                      Tsunami.System.pcom, Tsunami.System.vcom, Tsunami.Conf.wExt);

                // Exits if final time was already reached
                if (Tsunami.time >= tfin) {
                    return EXIT_SUCCESS;
                }
            }
        }
        if (!cont) {
            // Read input file
            IO::read_input_file(Tsunami.N, Tsunami.System.pos, Tsunami.System.vel, Tsunami.System.mass,
                                Tsunami.System.radius, Tsunami.System.xdata, Tsunami.System.spin, Tsunami.fname_in);
            Nb::scale_to_cdm(Tsunami.System.pos, Tsunami.System.vel, Tsunami.System.mass, Tsunami.N);
            // Compute e0 from particles
            Tsunami.energy = Nb::energy_calculation(Tsunami.System.pos, Tsunami.System.vel,
                                                    Tsunami.System.mass,Tsunami.N, Tsunami.pot, Tsunami.kin);
        }

        // Overwrite particle type
        if(!pthold.empty()) {
            if(pthold.size() != Tsunami.N) throw TsuError("Wrong number of particle type arguments (-pt), pass it after -N option");
            for(size_t i = 0; i<Tsunami.N; i++) {
                Tsunami.System.xdata[i].stype = static_cast<ptype>(pthold[i]);
            }
        }

        //if (Tsunami.System.spin) IO::print_spin(Tsunami.N, Tsunami.System.spin, 0.0);

        // Initialize tidal tables
        Tsunami.initialize_particle_tides();
        // Initialize alpha, beta, gamma for regularization
        Tsunami.initialize_regularization(1, 0, 0);

        // Initialize chain
        Tsunami.initialize_chain();
        //Tsunami.System.print_chain();
        //Tsunami.System.print_chained_vectors();

        // Initialize output
        IO Output(cont, Tsunami.N, Tsunami.fname_in, Tsunami.Conf.wExt);
        if(!cont) {
            Output.write_output(Tsunami.System.pos, Tsunami.System.vel, Tsunami.System.mass,
                                Tsunami.time, Tsunami.eoff, Tsunami.energy, Tsunami.System.spin,
                                Tsunami.System.xdata, Tsunami.System.pcom, Tsunami.System.vcom); // First output


        } else if (TsunamiConfig::useTiming) { // Let's print once, so we know the starting percentage
            Output.timing(Tsunami.dtphysical, Tsunami.time, tfin);
        }

        Tsunami.initialize_integrator();

        double time4output = Tsunami.time + dt_output;
        //Tsunami.timestep = 0.1;
        //Tsunami.nstep = 10;

        while (true) {
            Tsunami.do_step_bulirsh();

            Tsunami.System.update_from_chain_to_com();
            Tsunami.System.find_chain();

            if (Tsunami.check_collision()) {
                //Tsunami.iterate_to_collision_bulirsh(); //FIXME
                //Tsunami.revert_step();
                break;
            }

            Tsunami.update_time_coordinates_chain();

            // Check for input output
            if(Tsunami.time >= time4output) {
                Tsunami.energy = Nbodyalgorithms::energy_calculation(Tsunami.System.pos, Tsunami.System.vel,
                                                                     Tsunami.System.mass, Tsunami.N,
                                                                     Tsunami.pot, Tsunami.kin);
                Tsunami.deltaE = fabs(log(fabs((Tsunami.kin + (Tsunami.Eintegrator))/(Tsunami.pot))));
                Output.write_output(Tsunami.System.pos, Tsunami.System.vel, Tsunami.System.mass,
                                    Tsunami.time, Tsunami.deltaE+Tsunami.eoff, Tsunami.energy, Tsunami.System.spin,
                                    Tsunami.System.xdata, Tsunami.System.pcom, Tsunami.System.vcom);
                time4output = Tsunami.time + dt_output;
                if (TsunamiConfig::useTiming) {
                    Output.timing(Tsunami.dtphysical, Tsunami.time, tfin);
                }
            }
            if(Tsunami.time >= tfin) { // If reached target time, break out

                // One last print out if not already performed at this step
                if (time4output != Tsunami.time + dt_output) {
                    Tsunami.energy = Nbodyalgorithms::energy_calculation(Tsunami.System.pos, Tsunami.System.vel,
                                                                         Tsunami.System.mass, Tsunami.N,
                                                                         Tsunami.pot, Tsunami.kin);
                    Tsunami.deltaE = fabs(log(fabs((Tsunami.kin + (Tsunami.Eintegrator))/(Tsunami.pot))));
                    Output.write_output(Tsunami.System.pos, Tsunami.System.vel, Tsunami.System.mass,
                                        Tsunami.time, Tsunami.deltaE+Tsunami.eoff, Tsunami.energy, Tsunami.System.spin,
                                        Tsunami.System.xdata, Tsunami.System.pcom, Tsunami.System.vcom);
                }
                break;
            }
        }
        std::cout << std::endl; // Flush

        if constexpr (TsunamiConfig::useProfiling) {
            Tsunami.print_profiling();
        }

        if(Tsunami.System.CollInfo.collflag or Tsunami.stopcond) {
            // Last output during/after collision & stopping conditions
            Tsunami.energy = Nbodyalgorithms::energy_calculation(Tsunami.System.pos, Tsunami.System.vel,
                                                                 Tsunami.System.mass, Tsunami.N,
                                                                 Tsunami.pot, Tsunami.kin);
            Tsunami.deltaE = fabs(log(fabs((Tsunami.kin + (Tsunami.Eintegrator))/(Tsunami.pot))));
            Output.write_output(Tsunami.System.pos, Tsunami.System.vel, Tsunami.System.mass,
                                Tsunami.time, Tsunami.deltaE+Tsunami.eoff, Tsunami.energy, Tsunami.System.spin,
                                Tsunami.System.xdata, Tsunami.System.pcom, Tsunami.System.vcom);

            if (Tsunami.System.CollInfo.collflag) {
                Output.collision_log(Tsunami.System.CollInfo, Tsunami.System.pos, Tsunami.System.vel,
                                     Tsunami.System.mass, Tsunami.System.radius, Tsunami.time);
            }
        }

    } catch(TsunamiError &t) {} // End catch errors

    return EXIT_SUCCESS;
}