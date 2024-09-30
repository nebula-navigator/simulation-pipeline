/**
 * @file IO.cpp
 * @author Alessandro A. Trani
 * @author Mario Spera
 * @brief IO file of the TSUNAMI code
 * @details Contains the the main functions for Input/Output.
 * @date Tuesday, December 18, 2018
 */

#include <iostream>
#include <iomanip>
#include <algorithm>
#include "errhand.hpp"
#include "IO.h"
#include "Nbodyalgorithms.hpp"
#ifdef __linux__
#include "unistd.h"
#endif

namespace Nb = Nbodyalgorithms;

///////////////////////////////////////////////////////////////////
///						STATIC MEMBERS       				    ///
///////////////////////////////////////////////////////////////////

/**
 * Prints command line usage to terminal
 */
void IO::print_usage() {
    std::cout << "USAGE:" << std::endl;
    std::cout << "  tsunami.x [file.dat] [options...]" << std::endl << std::endl;
    std::cout << "  [file.dat] takes the first optional position," << std::endl;
    std::cout << "             otherwise it is read from input/tsunami_default_input.dat" << std::endl << std::endl;
    std::cout << "  Default options are read from input/tsunami_parameters.txt" << std::endl << std::endl;

    std::cout << "OPTIONS:" << std::endl;
    std::cout << "  file.dat             : initial condition file" << std::endl;
    std::cout << "  -N <value>           : number of particles" << std::endl;
    std::cout << "  -c                   : continue simulation from last timestep" << std::endl;
    std::cout << "  -ft <value>          : total integration time in nbody units" << std::endl;
    std::cout << "  -dt <value>          : timestep for output in nbody units" << std::endl;
    std::cout << "  -fcol <value>        : collision check multiplier" << std::endl;
    std::cout << "  -L <value>           : length scale (only for light speed in PN terms)" << std::endl;
    std::cout << "  -T                   : enable equilibrium tides, if not enabled in parameter file" << std::endl;
    std::cout << "  -TD                  : enable dynamical tides, if not enabled in parameter file" << std::endl;
    std::cout << "  -PN                  : enable post newtonians, if not enabled in parameter file" << std::endl;
    std::cout << "  -E                   : enable post external potentials, if not enabled in parameter file" << std::endl;
    std::cout << "  -pt <pt1> .. <ptN>   : particle type override" << std::endl;
    std::cout << "  -h                   : prints this screen" << std::endl;
}

/**
 * Parses the command line
 * @param[in] narg 			Number of command line argument
 * @param[in] varg 			Command line arguments pointer
 * @param[out] tfin 		Final integration time
 * @param[out] dt_output 	Timestep for output
 * @param[out] dcoll 		Collisions check multiplier
 * @param[out] cont 		Continue simulation flag
 * @param[out] in_fname 	Initial condition filename
 * @param[out] outname 		Output filename
 * @param[out] enename 		Energy output filename
 * @param[out] collname 	Collision log filename
 * @param[out] wPNs			Postnewtonian flag
 * @param[out] wTides		Tides flag
 * @param[out] N			Number of particles
 */

void IO::read_command_line(int narg, char **varg, double &tfin, double &dt_output, double &Lscale,
                           double &Mscale, bool &cont, TsunamiConfig &Config, string &in_fname, size_t &N,
                           std::vector<int> &pthold) {

    for (int i = 1; i < narg; i++) {
        string argv_got = string(&varg[i][0]);
        //std::cout<<"argv_got   "<<argv_got<<std::endl;

        if (argv_got == string("-c")) {
            cont = true;
        } else if (argv_got == string("-ft")) {
            tfin = std::strtod(&varg[i + 1][0], nullptr);
            if (tfin <= 0)
                throw TsuError("Total integration time cannot be negative (-ft)");
            i++;

        } else if (argv_got == string("-N")) {
            N = (size_t) std::strtol(&varg[i + 1][0], nullptr, 10);
            if (N <= 1)
                throw TsuError("Can't simulate less than 2 particles (-N)");
            i++;

        } else if (argv_got == string("-dt")) {
            dt_output = std::strtod(&varg[i + 1][0], nullptr);
            if (dt_output <= 0)
                throw TsuError("Interval time for output cannot be negative (-dt)");
            i++;

        } else if (argv_got == string("-fcol")) {
            Config.dcoll = std::strtod(&varg[i + 1][0], nullptr);
            if (Config.dcoll < 0)
                throw TsuError("Collision check multiplier cannot be negative (-fcol)");
            i++;

        } else if (argv_got == string("-L")) {
            if (string(&varg[i + 1][0]) == "au") {
                Lscale = 1;
                std::cout << "Assuming 1 astronomical unit as length scale" << std::endl;
            } else if (string(&varg[i + 1][0]) == "pc") {
                Lscale = 2.06265e+5;
                std::cout << "Assuming 1 parsec as length scale" << std::endl;
            } else {
                Lscale = std::strtod(&varg[i + 1][0], nullptr);
            }
            if (Lscale <= 0)
                throw TsuError("Length scale cannot be negative (-L)");
            i++;

        } else if (argv_got == string("-M")) {
            if (string(&varg[i + 1][0]) == "Msun") {
                Mscale = 1;
                std::cout << "Assuming 1 solar mass as mass scale" << std::endl;
            } else {
                Lscale = std::strtod(&varg[i + 1][0], nullptr);
            }
            if (Mscale <= 0)
                throw TsuError("Mass scale cannot be negative (-M)");
            i++;

        } else if (argv_got == string("-T")) {
            Config.wEqTides = true;
            std::cout << "Overriding parameter.txt: equilibrium tides on" << std::endl;

        } else if (argv_got == string("-TD")) {
            Config.wDynTides = true;
            std::cout << "Overriding parameter.txt: dynamical tides on" << std::endl;

        } else if (argv_got == string("-PN")) {
            Config.wPNs = true;
            std::cout << "Overriding parameter.txt: PNs on" << std::endl;

        } else if (argv_got == string("-E")) {
            Config.wExt = true;
            std::cout << "Overriding parameter.txt: external potentials on" << std::endl;

        } else if (argv_got == string("-pt")) {
            std::cout << "Overriding particle types: PNs on" << std::endl;
            for (size_t k = 0; k < N; k++) {
                try {
                    pthold.push_back(std::stoi(&varg[i + 1][0]));
                } catch (const std::invalid_argument &) {
                    throw TsuError("Wrong number of particle types (-pt), must be equal to N");
                } catch (const std::logic_error &) {
                    throw TsuError("Wrong number of particle types (-pt), must be same as N");
                }
                i++;
            }

        } else if (argv_got == string("-h")) {
            print_usage();
            exit(EXIT_SUCCESS);

        } else if (i == 1) {
            /* if it's not a recognized option and it is the first argument  *
             * let it be the input file                                     */
            in_fname = &varg[i][0];
            size_t found = in_fname.find_last_of('.');
            if (found == std::string::npos) {
                print_usage();
                throw TsuError(in_fname + " not recognized either as an option or input file (should have an extension)");
            }

        } else {
            IO::print_usage();
            throw TsuError(argv_got + " not recognized");
        }
    }
}

string IO::get_exe_path() {
#ifdef __linux__
    const size_t bufsiz = 200;
    char buf[bufsiz] = ""; // Otherwise valgrind complains
    readlink("/proc/self/exe", buf, bufsiz);
    string exepath(buf);
    size_t found = exepath.find_last_of('/');
    exepath = exepath.substr(0, found) + '/';

    //std::cout << "Found exepath " << exepath << std::endl;
#else
    string exepath = "";
#endif
    return exepath;
}

/**
 * Reads the parameter file
 * @param N				Particle number
 * @param tfin			Final integration time
 * @param Mscale		Mass scale
 * @param Rscale		Length scale
 * @param dt_output		Timestep for output
 * @param wPNs			Postnewtonian flag
 * @param wEqTides		Tides flag
 * @param dcoll			Collisions check multiplier
 * @param nextStepFac	BS next step factor
 */
void IO::read_parameter_file(string &paramfile, size_t &N, double &tfin, double &Mscale, double &Rscale,
                             double &dt_output, TsunamiConfig &Config, double &dcoll) {

    string line;
    string withPNs, withEqTides, withDynTides, withExternal;
    std::istringstream stream;
    std::ifstream in;
    in.open(paramfile, std::ios::in);
    if (!in) throw TsuError("Cannot open file " + paramfile);

    getline(in, line);
    if (in.eof()) throw TsuError("First line of parameters.txt is empty.. please check your file");
    stream.str(line);
    stream >> N;
    if (N <= 1) throw TsuError("Cannot integrate a system with less than 1 body");
    stream.clear();

    getline(in, line);
    if (in.eof()) throw TsuError("Cannot read total integration time.. please check your file");
    stream.str(line);
    stream >> tfin;
    if (tfin <= 0) throw TsuError("Total integration time cannot be negative");
    stream.clear();

    getline(in, line);
    if (in.eof()) throw TsuError("Cannot read scale mass value.. please check your file");
    stream.str(line);
    stream >> Mscale;
    if (Mscale <= 0) throw TsuError("Mass scale cannot be negative");
    stream.clear();

    getline(in, line);
    if (in.eof()) throw TsuError("Cannot read scale length value.. please check your file");
    stream.str(line);
    stream >> Rscale;
    if (Rscale <= 0) throw TsuError("Length scale cannot be negative");
    stream.clear();

    getline(in, line);
    if (in.eof()) throw TsuError("Cannot read interval of time for output.. please check your file");
    stream.str(line);
    stream >> dt_output;
    if (dt_output <= 0) throw TsuError("Interval time for output cannot be negative");
    stream.clear();

    getline(in, line);
    if (in.eof())
        throw TsuError("Cannot read if you want to use PN or not.. please check your file");
    stream.str(line);
    stream >> withPNs;
    if (withPNs != "yes" && withPNs != "no")
        throw TsuError("Please use just 'yes' or 'no' for PN terms");
    stream.clear();

    getline(in, line);
    if (in.eof())
        throw TsuError("Cannot read if you want to use equilibrium tides or not.. please check your file");
    stream.str(line);
    stream >> withEqTides;
    if (withEqTides != "yes" && withEqTides != "no")
        throw TsuError("Please use just 'yes' or 'no' for tidal term");
    stream.clear();

    getline(in, line);
    if (in.eof())
        throw TsuError("Cannot read if you want to use dynamical tides or not.. please check your file");
    stream.str(line);
    stream >> withDynTides;
    if (withDynTides != "yes" && withDynTides != "no")
        throw TsuError("Please use just 'yes' or 'no' for tidal term");
    stream.clear();


    getline(in, line);
    if (in.eof())
        throw TsuError("Cannot read collision check distance factor.. please check your file");
    stream.str(line);
    stream >> dcoll;
    if (dcoll < 0) throw TsuError("Collision check multiplier cannot be negative");
    stream.clear();


    if (withEqTides == "no") Config.wEqTides = false;
    else if (withEqTides == "yes") Config.wEqTides = true;

    if (withDynTides == "no") Config.wDynTides = false;
    else if (withDynTides == "yes") Config.wDynTides = true;

    if (withPNs == "no") Config.wPNs = false;
    else if (withPNs == "yes") Config.wPNs = true;

    getline(in, line);
    if (in.eof())
        throw TsuError("Cannot read if you want to use external potentials or not.. please check your file");
    stream.str(line);
    stream >> withExternal;
    if (withExternal != "yes" && withExternal != "no")
        throw TsuError("Please use just 'yes' or 'no' for external potentials");
    stream.clear();

    if (withExternal == "no") Config.wExt = false;
    else if (withExternal == "yes") Config.wExt = true;

    in.close();
}

bool IO::read_line(std::ifstream &in, std::istringstream &stream) {
    stream.clear();
    std::string line;
    getline(in, line);
    stream.str(line);
    return !in.eof();
}

void IO::print_spin(size_t N, double3 *spin, double time) {
    std::cout << "Spin t=" << time << std::endl;
    for (size_t i = 0; i < N; i++) {
        std::cout << spin[i].x << "   " << spin[i].y << "   " << spin[i].z << std::endl;
    }
    std::cout << std::endl;
}

/**
 * Reads the initial condition file
 * @param[in] N 		Number of particles
 * @param[out] p 		Positions
 * @param[out] v 		Velocities
 * @param[out] m 		Masses
 * @param[out] rad 		Radii
 * @param[out] xdata 	Particles extra data
 * @param[out] s 		softening
 * @param[in] fname 	Input filename
 */
void IO::read_input_file(size_t N, double3 *pos, double3 *vel, double *mass, double *radius, pinfo *xdata,
                         double3 *spin, const string &fname) {
    std::ifstream in;
    in.open(fname.c_str(), std::ios::in);
    if (!in) throw TsuError("Cannot open file " + fname);

    std::istringstream stream;
    //format of input file = x, y, z, vx, vy, vz, mass, radius, id, xspin, yspin, zspin
    size_t Ncols = Ncols_input;
    size_t Nrows = 0;
    size_t colsfound = 0;

    if (spin) {
        Ncols += 3;
    }

    for (size_t i = 0; i < N; i++) {
        std::vector<double> tmp;
        std::string value;
        read_line(in, stream);
        Nrows++;
        while (stream >> value) {
            tmp.push_back(s2n(double, value));
            stream.clear();
        }

        if (tmp.empty() and Nrows <= N) {
            throw TsuError("Wrong number of rows: expected " + n2s(N) + ", found " + n2s(Nrows-1));
        } else if (tmp.size() < Ncols) {

            throw TsuError("Wrong number of columns: expected " + n2s(Ncols) + ", found " + n2s(tmp.size())
            + " (maybe missing spin vectors?)");
        }
        colsfound = tmp.size();

        pos[i] = {tmp[0], tmp[1], tmp[2]};
        vel[i] = {tmp[3], tmp[4], tmp[5]};
        mass[i] = tmp[6];
        radius[i] = tmp[7];
        xdata[i].stype = static_cast<ptype>(tmp[8]);

        if (spin) {
            spin[i] = {tmp[9], tmp[10], tmp[11]};
        }
    }
    if (colsfound > Ncols) {
        std::cout << "Found " << colsfound << " columns, expected " << Ncols << ". Ignoring the remaining ones"
                  << std::endl;
    }
    in.close();
}


/**
 * Reads the initial condition file
 * @param[in] N 		Number of particles
 * @param[out] datastring  Datastring using Brutus format
 * @param[in] fname 	Input filename
 */
void IO::read_input_file_as_string(size_t N, std::vector<string> &datastring, const string &fname) {
    std::ifstream in;
    in.open(fname.c_str(), std::ios::in);
    if (!in) throw TsuError("Cannot open file " + fname);

    std::array<string, 7> coord_tmp;
    string dummy;
    for (size_t i = 0; i < N; i++) {
        /// X Y Z  VX VY VZ M S R PT
        in >> coord_tmp[1] >> coord_tmp[2] >> coord_tmp[3] >> coord_tmp[4] >> coord_tmp[5] >> coord_tmp[6]
           >> coord_tmp[0] >> dummy >> dummy;
        for (auto &tmp : coord_tmp)
            datastring.push_back(tmp);
    }
    in.close();
}

/**
 * Reads the last timestep from output file, and loads the last position, masses and velocities
 * Returns true if the last timestep is found and initial conditions are correctly set from the last output,
 * otherwise return false and it will read back the initial condition file
 * @param[in] N 		Number of particles
 * @param[out] pos		Positions
 * @param[out] vel		Velocities
 * @param[out] mass		Masses
 * @param[out] radius		Radii
 * @param[out] xdata	Particles extra data
 * @param[out] s		Softening
 * @param[in] outfile	Output filename
 * @param[in] icfile	Initial condition filename
 * @return				true if last output was found and loaded, otherwise false
 */
bool IO::read_last_timestep(size_t N, double3 *pos, double3 *vel, double *mass, double *radius, pinfo *xdata,
                            double3 *spin, const string &outfile, const string &icfile) {

    // First read input file
    read_input_file(N, pos, vel, mass, radius, xdata, spin, icfile);

    std::ifstream in;
    // Then read last timestep for particle positions and velocities
    in.open(outfile.c_str(), std::ios::in | std::ios::binary);
    if (!in) {
        std::cout << "No output file, reading from " << icfile << std::endl;
        return false;
    } //throw TsunamiError("Cannot open file "+icfile, __FILE__, __LINE__);

    in.seekg(0, std::ios::end);
    std::streamoff end = in.tellg();
    std::streamoff position = end;

    // Counting '\n' backwards
    char curr;
    size_t count = 0;
    while (position) {
        in.seekg(--position);
        in.get(curr);

        if (curr == '\n') {
            if (count++ == N) {
                break;
            }
        }
    }

    std::streamoff bufsize = in.tellg();
    bufsize = end - bufsize;
    char *buffe = new char[bufsize];
    in.read(buffe, bufsize);

    size_t Ncols = Ncols_output;
    if (spin) {
        Ncols += 3;
    }

    std::stringstream buffer;
    buffer << buffe;
    // override position, velocities and masses with last timestep
    std::istringstream linestream;
    for (size_t i = 0; i < N; i++) {
        std::vector<double> tmp;
        std::string value;

        linestream.clear();
        std::string line;
        getline(buffer, line);
        linestream.str(line);

        while (linestream >> value) {
            tmp.push_back(s2n(double, value));
            linestream.clear();
        }

        if (tmp.size() < Ncols) {
            throw TsuError("Wrong number of columns: expected " + n2s(Ncols) + ", found " + n2s(tmp.size())
                           + " (maybe missing spin vectors?)");
        }

        pos[i] = {tmp[0], tmp[1], tmp[2]};
        vel[i] = {tmp[3], tmp[4], tmp[5]};
        mass[i] = tmp[6];
        xdata[i].eloss = tmp[7];
        if (spin) {
            spin[i] = {tmp[8], tmp[9],tmp[10]};
        }
    }
    delete[] buffe;
    in.close();
    return true;
}

/**
 * Reads the last timestep from output file, and loads the last position, masses and velocities
 * Returns true if the last timestep is found and initial conditions are correctly set from the last output,
 * otherwise return false and it will read back the initial condition file
 * @param[in] N 		Number of particles
 * @param[out] p		Positions
 * @param[out] v		Velocities
 * @param[out] m		Masses
 * @param[out] rad		Radii
 * @param[out] xdata	Particles extra data
 * @param[out] s		Softening
 * @param[in] outfile	Output filename
 * @param[in] icfile	Initial condition filename
 * @return				true if last output was found and loaded, otherwise false
 */
bool IO::read_last_timestep_as_string(size_t N, std::vector<string> &datastring, const string &outfile,
                                      const string &icfile) {

    std::ifstream in;
    // then read last timestep for particle positions and velocities
    in.open(outfile.c_str(), std::ios::in | std::ios::binary);
    if (!in) {
        std::cout << "No output file, reading from " << icfile << std::endl;
        return false;
    } //throw TsunamiError("Cannot open file "+icfile, __FILE__, __LINE__);

    in.seekg(0, std::ios::end);
    std::streamoff end = in.tellg();
    std::streamoff pos = end;

    char curr;
    size_t count = 0;
    while (pos) {
        in.seekg(--pos);
        in.get(curr);

        if (curr == '\n') {
            if (count++ == N) {
                break;
            }
        }
    }
    std::streamoff bufsize = in.tellg();
    bufsize = end - bufsize;
    char *buffe = new char[bufsize];
    in.read(buffe, bufsize);

    std::stringstream buffer;
    buffer << buffe;
    std::array<string, 7> coord_tmp;
    for (size_t i = 0; i < N; i++) {
        /// X Y Z  VX VY VZ  M
        buffer >> coord_tmp[1] >> coord_tmp[2] >> coord_tmp[3] >> coord_tmp[4] >> coord_tmp[5] >> coord_tmp[6]
               >> coord_tmp[0];
        for (auto &tmp : coord_tmp)
            datastring.push_back(tmp);
    }

    delete[] buffe;
    in.close();
    return true;
}


/**
 * Gets energy and time infos for restarting a simulation
 * @param[in] timelast		Time of last output
 * @param[in] e0			Energy of last output
 * @param[in] eoff			DeltaE of last output
 * @param[in] energy_out	Energy filename
 */
void IO::get_restart_infos(double &timelast, double &e0, double &eoff, const string &energy_out,
                           double3 &pcom, double3 &vcom, bool wExt) {

    std::ifstream ene;
    ene.open(energy_out.c_str(), std::ios::in | std::ios::binary);
    if (!ene) throw TsuError("Cannot open file " + energy_out);

    double dummy;
    ene >> dummy >> dummy >> e0;
    //std::cout << "e0: " << e0 << std::endl;

    ene.seekg(0, std::ios::end);
    std::streamoff end = ene.tellg();
    std::streamoff pos = end;

    int count = 0;
    char curr;
    while (pos) {
        ene.seekg(--pos);
        ene.get(curr);
        if (curr == '\n') {
            if (count++ == 1) {
                break;
            }
        }
    }
    std::streamoff bufsize = ene.tellg();
    bufsize = end - bufsize;
    char *buffe = new char[bufsize];
    ene.read(buffe, bufsize);

    std::stringstream buffer;

    buffer << buffe;

    size_t Ncols = Ncols_energy;
    if (wExt) {
        Ncols += 6;
    }

    std::istringstream linestream;
    std::vector<double> tmp;
    std::string value;

    linestream.clear();
    std::string line;
    getline(buffer, line);
    linestream.str(line);

    while (linestream >> value) {
        tmp.push_back(s2n(double, value));
        linestream.clear();
    }

    if (tmp.size() < Ncols) {
        throw TsuError("Wrong number of columns: expected " + n2s(Ncols) + ", found " + n2s(tmp.size())
                       + " (maybe missing CoM position and velocities?)");
    }

    timelast = tmp[0];
    eoff = tmp[1];
    // energy = tmp[2] not needed

    if (wExt) {
        pcom = {tmp[3], tmp[4], tmp[5]};
        vcom = {tmp[6], tmp[7], tmp[8]};
    }

    delete[] buffe;
    ene.close();
}

/**
 * Gets energy and time infos for restarting a simulation
 * @param[out] timelast		Time of last output
 * @param[out] e0			Energy of last output
 * @param[out] eoff			DeltaE of last output
 * @param[in] energy_out	Energy filename
 */
void IO::get_restart_infos_as_string(string &timelast, string &e0, string &eoff, const string &energy_out) {
    std::ifstream ene;
    ene.open(energy_out.c_str(), std::ios::in | std::ios::binary);
    if (!ene) throw TsuError("Cannot open file " + energy_out);

    double dummy;
    ene >> dummy >> dummy >> e0;
    //std::cout << "e0: " << e0 << std::endl;

    ene.seekg(0, std::ios::end);
    std::streamoff end = ene.tellg();
    std::streamoff pos = end;

    int count = 0;
    char curr;
    while (pos) {
        ene.seekg(--pos);
        ene.get(curr);
        if (curr == '\n') {
            if (count++ == 1) {
                break;
            }
        }
    }
    std::streamoff bufsize = ene.tellg();
    bufsize = end - bufsize;
    char *buffe = new char[bufsize];
    ene.read(buffe, bufsize);

    std::stringstream buffer;
    buffer << buffe;
    buffer >> timelast >> eoff >> dummy;
    delete[] buffe;
    ene.close();
}

void IO::get_output_file_names(const string &in_fname, string &ene_name, string &out_name, string &coll_name) {
    size_t found = in_fname.find_last_of('.');
    if (found == std::string::npos) {
        print_usage();
        throw TsuError(in_fname + " not recognized either as an option or input file (should have an extension)");
    }
    string basename = in_fname.substr(0, found);

    coll_name = basename + "_collision.dat";
    out_name = basename + "_output.dat";
    ene_name = basename + "_energy.dat";
}


void IO::datastring_to_chreal(size_t N, double3 *pos, double3 *vel, double *mass, std::vector<string> datastring) {
    for (size_t i = 0; i < N; i++) {
        mass[i] = std::stod(datastring[i * 7 + 0]);
        pos[i].x = std::stod(datastring[i * 7 + 1]);
        pos[i].y = std::stod(datastring[i * 7 + 2]);
        pos[i].z = std::stod(datastring[i * 7 + 3]);
        vel[i].x = std::stod(datastring[i * 7 + 4]);
        vel[i].y = std::stod(datastring[i * 7 + 5]);
        vel[i].z = std::stod(datastring[i * 7 + 6]);
    }
}

/**
 * Print particles to screen
 * @param[in] N Particle number
 * @param[in] pos Positions
 * @param[in] vel Velocities
 * @param[in] mass Masses
 * @param[in] soft Softening
 * @param[in] rad Radii
 * @param[in] xdata Particles extra data
 */
void IO::print_particles(size_t N, double3 *pos, double3 *vel, double *mass, double *soft,
                         double *rad, pinfo *xdata) {
    for (size_t i = 0; i < N; i++) {
        std::cout << pos[i].x << "   " << pos[i].y << "   " << pos[i].z << "   " << vel[i].x << "   " << vel[i].y
                  << "   "
                  << vel[i].z << "   " << mass[i] << "   " << soft[i] << "   " << rad[i] << "   " << xdata[i].stype
                  << std::endl;
    }
    std::cout << std::endl;
}

/**
 * * Print particles to screen, only positions and velocities
 * @param[in] N Particle number
 * @param[in] pos Positions
 * @param[in] vel Velocities
 */
void IO::print_particles(size_t N, double3 *pos, double3 *vel) {
    for (size_t i = 0; i < N; i++) {
        std::cout << pos[i].x << "   " << pos[i].y << "   " << pos[i].z << "   " << vel[i].x << "   " << vel[i].y
                  << "   " << vel[i].z << std::endl;
    }
    std::cout << std::endl;
}

/**
 *  Writes initial condition file
 * @param[in] N Particles number
 * @param[in] p Positions
 * @param[in] v Velocities
 * @param[in] m Masses
 * @param[in] rad Radii
 * @param[in] xdata Particles extra data
 * @param[in] s softenings
 * @param[in] fname file name
 */
void IO::write_input_file(size_t N, double3 *p, double3 *v, double *m, double *rad, pinfo *xdata, double *s,
                          const string &fname) {
    std::ofstream finput;
    finput.open(fname.c_str(), std::ios::out);
    if (!finput) throw TsuError("Cannot open file " + fname);

    for (size_t i = 0; i < N; i++) {
        finput << p[i].x << "    " << p[i].y << "    " << p[i].z << "    " << v[i].x << "    " << v[i].y << "    "
               << v[i].z << "    " << m[i] << "    " << s[i] << "    " << rad[i] << "    " << xdata[i].stype;
    }
    finput.close();
}

void IO::read_tidal_table(std::map<short int, Classification::TideTableElement> &tidetable, string &tidefile) {

    string line;
    std::istringstream stream;
    std::ifstream in;
    in.open(tidefile, std::ios::in);
    if (!in) {
        std::cerr << "WARNING: No tidal table found (" << tidefile << ")" << std::endl;
        std::cerr << "         Tidal parameters will revert to default values" << std::endl;
        return;
    }

    short int ind;
    Classification::TideTableElement tidelem;
    while (getline(in, line)) {
        // Erase comments
        line.erase(std::find(line.begin(), line.end(), '#'), line.end());
        stream.str(line);

        if (!stream.str().empty()) {
            stream >> ind >> tidelem.taulag >> tidelem.kaps >> tidelem.polyt >> tidelem.gyradius;
            //if (tidelem.taulag <= 0 or tidelem.kaps <= 0 or tidelem.polyt <= 0)
            //    throw TsunamiError("Taulag, klove and polytropic index cannot be negative or zero", __FILE__, __LINE__);
            tidetable[ind] = tidelem;
        }
    }
#ifdef TIDELOG
    std::cout << "Tidal table elements:" << std::endl;
    std::cout << " type          taulag        kaps           polyind" << std::endl;
    for (auto elem : tidetable)
        std::cout << std::setw(5) << elem.first << std::setw(15) << elem.second.taulag << std::setw(15)
                  << elem.second.kaps << std::setw(15) << elem.second.polyt << std::endl;
#endif // TIDELOG
}

///////////////////////////////////////////////////////////////////
///						 NON-STATIC MEMBERS       				///
///////////////////////////////////////////////////////////////////

/**
 * Writes the current timestep into the output and energy files
 * @param[in] pos Positions
 * @param[in] vel Velocities
 * @param[in] mass Masses
 * @param[in] time Time
 * @param[in] deltaE Energy error
 * @param[in] E Energy
 */
void IO::write_output(double3 *pos, double3 *vel, double *mass, double time, double deltaE, double E,
                      double3 *spin, pinfo *xdata, double3 &pcom, double3 &vcom) {

    for (size_t i = 0; i < Npart; i++) {
        output << pos[i].x << "   " << pos[i].y << "   " << pos[i].z << "   " << vel[i].x << "   " << vel[i].y << "   "
               << vel[i].z << "   " << mass[i] << "   " << xdata[i].eloss;
        if (spin) output << "   " << spin[i].x << "   " << spin[i].y << "   " << spin[i].z;
        output << "\n";
    }
    energy << time << "   " << deltaE << "   " << E;
    if(has_excom) {
        energy << "     " << pcom.x << "   " << pcom.y << "   " << pcom.z << "    "
        << vcom.x << "   " << vcom.y << "   " << vcom.z;
    }
    energy << "\n";
}

/**
 * Writes the current timestep into the output and energy files
 * @param[in] datastring  Datastring using Brutus format
 * @param[in] time Time
 * @param[in] deltaE Energy error
 * @param[in] E Energy
 */
void IO::write_output_as_string(const std::vector<string> &datastring, string time, string deltaE, string E) {

    for (size_t i = 0; i < Npart; i++)
        output << datastring[i * 7 + 1] << "   "
               << datastring[i * 7 + 2] << "   "
               << datastring[i * 7 + 3] << "   "
               << datastring[i * 7 + 4] << "   "
               << datastring[i * 7 + 5] << "   "
               << datastring[i * 7 + 6] << "   "
               << datastring[i * 7 + 0] << std::endl;

    energy << time << "   " << deltaE << "   " << E << std::endl;
}

/**
 * Writes the collision log in case of collisions
 * @param[in] collision_type Type of collision
 * @param[in] colldata Collision data
 * @param[in] pos Positions
 * @param[in] vel Velocities
 * @param[in] mass Masses
 * @param[in] rad Radii
 * @param[in] time Time
 */
void IO::collision_log(Nb::CollisionLog &CollInfo, double3 *pos, double3 *vel, double *mass,
                       double *rad, double time) {

    size_t k = CollInfo.collind[0];
    size_t j = CollInfo.collind[1];
    time = time + CollInfo.colltime;

    switch (CollInfo.collflag) {
        case Nb::CollType::NOCOLLISION :
            throw TsuError("ERROR: collision_log function called, but collflag is NOCOLLISION");
            break;
        case Nb::CollType::COLLISION_CANDIDATE :
            std::cout << std::scientific << std::setprecision(6) << icname << ": [COLLISION] PREDICTED AT TIME " << time
                      << " (pair " << k << "," << j << ")" << std::endl;
            break;
        case Nb::CollType::COLLISION_REAL :
            std::cout << std::scientific << std::setprecision(6) << icname << ": [COLLISION] DETECTED AT TIME " << time
                      << " (pair " << k << "," << j << ")" << std::endl;
            break;
        case Nb::CollType::COLLISION_TDE :
            std::cout << std::scientific << std::setprecision(6) << icname << ": [TIDAL DISRUPTION] DETECTED AT TIME " << time
            << " (pair " << k << "," << j << ")" << std::endl;
            break;
    }
    /* time
     * particle1_index x y z vx vy vz mass radius
     * particle2_index x y z vx vy vz mass radius      */
    std::ofstream coll;
    coll.open(collname.c_str(), std::ios::out);
    coll << std::scientific << std::setprecision(16);

    coll << time << "   " << CollInfo.collflag << std::endl;
    coll << k << "   " << pos[k].x << "   " << pos[k].y << "   " << pos[k].z << "   " << vel[k].x << "   " << vel[k].y
         << "   " << vel[k].z << "   " << mass[k] << "   " << rad[k] << std::endl;
    coll << j << "   " << pos[j].x << "   " << pos[j].y << "   " << pos[j].z << "   " << vel[j].x << "   " << vel[j].y
         << "   " << vel[j].z << "   " << mass[j] << "   " << rad[j] << std::endl;
    coll.close();
}

/**
 * Logs current timestep, number of particles and time percentage
 * @param[in] dt 	Timestep
 * @param[in] time	Current time
 * @param[in] tfin	Final time
 */
void IO::timing(const double &dt, const double &time, const double &tfin) {
    std::cout << "\r";
    std::cout << std::fixed << std::setprecision(2) << std::setw(3) << "N = " << Npart << "     ";
    std::cout << std::scientific << std::setprecision(6) << std::setw(10) << "t = " << time << "     ";
    std::cout << "dt = " << dt << "      ";
    std::cout << std::fixed << std::setprecision(2) << "t/tfin= " << std::setw(5) << time / tfin * 100.0 << " %"
              << std::flush;
    std::cout << std::scientific << std::setprecision(15);
}
