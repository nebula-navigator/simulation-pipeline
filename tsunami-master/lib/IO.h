//
// Created by lex on 22/12/18.
//

#ifndef TSUNAMI_IO_H
#define TSUNAMI_IO_H

#include <array>
#include <fstream>
#include <iomanip>
#include "custom_types.hpp"
#include "errhand.hpp"
#include "Nbodyalgorithms.hpp"
#include "classification.h"

using std::string;

class IO {
public:
	IO(bool cont, size_t N, string &in_fname, bool has_center_of_mass) {

		Npart = N;

        icname = in_fname;
        has_excom = has_center_of_mass;

        get_output_file_names(in_fname, enename, outname, collname);

        // Open output files
        if(cont) {
			output.open(outname.c_str(), std::ios::out | std::ios::app);
			energy.open(enename.c_str(), std::ios::out | std::ios::app);
		} else {
			output.open(outname.c_str(), std::ios::out);
			energy.open(enename.c_str(), std::ios::out);
		}

		output<<std::scientific<<std::setprecision(16);
		energy<<std::scientific<<std::setprecision(16);
		std::cout << std::scientific << std::setprecision(15);
	}

	~IO() {
		output.close();
		energy.close();
	}

	// Static, no need for initialization, for input and print to screen only

	static void print_usage();

	static void read_command_line(int narg, char **varg, double &tfin, double &dt_output, double &Lscale,
                      double &Mscale, bool &cont, TsunamiConfig &Config, string &in_fname, size_t &N,
                      std::vector<int> &pthold);

	static void read_parameter_file(string &paramfile, size_t &N, double &tfin, double &Mscale, double &Rscale,
                                    double &dt_output, TsunamiConfig &Config, double &dcoll);

	static void read_input_file(size_t N, double3 *pos, double3 *vel, double *mass, double *radius, pinfo *xdata,
                                double3 *spin, const string &fname);

	static void get_output_file_names(const string &in_fname, string &ene_name, string &out_name, string &coll_name);

	static bool read_last_timestep(size_t N, double3 *pos, double3 *vel, double *mass,
                                   double *rad, pinfo *xdata, double3 *spin,
                                   const string &outfile, const string &icfile);

    static bool read_last_timestep_as_string(size_t N, std::vector<string> & datastring,
                                             const string &outfile, const string &icfile);

    static void get_restart_infos(double &timelast, double &e0, double &eoff,
                                  const string &energy_out, double3 &pcom, double3 &vcom, bool wExt);

    static void get_restart_infos_as_string(string &timelast, string &e0, string &eoff, const string &energy_out);

    static void read_input_file_as_string(size_t N, std::vector<string> & datastring, const string &fname);

	static void print_particles(size_t N, double3 *pos, double3 *vel, double *mass, double *soft, double *rad, pinfo *xdata);
	static void print_particles(size_t N, double3 *pos, double3 *vel);

	static void write_input_file(size_t N, double3 *p, double3 *v, double *m, double *rad, pinfo *xdata, double *s,
	                             const string &fname);
	static void datastring_to_chreal(size_t N, double3 *pos, double3 *vel, double *mass,
                                     std::vector<string> datastring);

	static void read_tidal_table(std::map<short int , Classification::TideTableElement> &tidetable, string &tidefile);

    static string get_exe_path();

    static bool read_line(std::ifstream &in, std::istringstream &stream);

    static void print_spin(size_t N, double3 *spin, double time);

    // Non-static, need initialization, for output only

	void write_output(double3 *pos, double3 *vel, double *mass, double time, double deltaE,
                      double E, double3 *spin, pinfo *xdata, double3 &pcom, double3 &vcom);

	void write_output_as_string(const std::vector<string> & datastring, string time, string deltaE, string E);

	void collision_log(Nbodyalgorithms::CollisionLog &CollInfo, double3 *pos, double3 *vel, double *mass,
	                       double *rad, double time);

	void timing(const double &dt, const double &time, const double &tfin);

private:

	std::ofstream output, energy; //output file to collect all pos,vel,mass over time and energy variation
	string collname, outname, icname, enename, excomname;
	bool has_excom;
	size_t Npart;

    static constexpr size_t Ncols_input = 9;
    static constexpr size_t Ncols_output = 8;
    static constexpr size_t Ncols_energy = 3;

};


#endif //TSUNAMI_IO_H
