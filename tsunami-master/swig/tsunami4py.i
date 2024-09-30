%module tsunami

// Begin module
%{
    #include "keplerutils.h"
    #include "tsunami.hpp"
    #include "kozailidov.h"
    #define SWIG_FILE_WITH_INIT
%}
%include "numpy.i"

%init %{
	import_array();
%}
%fragment("NumPy_Fragments");

%include <typemaps.i>
%include <std_array.i>
%include <std_string.i>
%include <exception.i>
%include <std_vector.i>
%include <std_map.i>
%include <std_pair.i>
%feature("autodoc", "3");
%feature("naturalvar");

%rename (CollLogger) Nbodyalgorithms::CollisionLog;
%rename (PTDELogger) Nbodyalgorithms::PTDELog;
%rename (PTDEevent) Nbodyalgorithms::PTDE_event;
%rename (PTDEhistory) Nbodyalgorithms::PTDE_history;
%rename (Tsunami) TsunamiClass;
%template(PTDEevent_vector) std::vector<Nbodyalgorithms::PTDE_event>;
%template(size_t_vector) std::vector<size_t>;
%template(PTDEhistory_map) std::map<size_t, Nbodyalgorithms::PTDE_history>;

%immutable TsunamiClass::System;
%immutable TsunamiClass::BSExtra;
%immutable TsunamiClass::Leapfrog;

// Read only members
%immutable TsunamiClass::deltaE;
%immutable TsunamiClass::energy;
%immutable TsunamiClass::dtphysical;
%immutable TsunamiClass::N;
%immutable TsunamiClass::eoff;
%immutable TsunamiClass::pot;
%immutable TsunamiClass::kin;
%immutable Mscale;
%immutable Lscale;
%immutable Tscale;
%immutable Vscale;

//%ignore Tsunami::eoff;
//%ignore Tsunami::Marx;
//%ignore Tsunami4py::System;
//%ignore Tsunami::BSExtra;
//%ignore Tsunami::Leapfrog;
//%ignore Tsunami::Eintegrator;
//%ignore Tsunami::Conf;

%catches(TsunamiError) TsunamiClass::add_particle_set;
%catches(TsunamiError) TsunamiClass::sync_internal_state;
%catches(TsunamiError) TsunamiClass::override_masses;
%catches(TsunamiError) TsunamiClass::override_position_and_velocities;
%catches(TsunamiError) TsunamiClass::initialize_tidal_parameters;
%catches(TsunamiError) TsunamiClass::evolve_system;
%catches(TsunamiError) TsunamiClass::print_profiling;
%catches(TsunamiError) TsunamiClass::save_restart_file;
%catches(TsunamiError) TsunamiClass::load_restart_file;
%catches(TsunamiError) TsunamiClass::get_chain_vectors;
%catches(TsunamiError) TsunamiClass::get_accelerations_of_particle_pair;
%catches(TsunamiError) TsunamiClass::sync_masses;

%catches(TsunamiError) KeplerUtils::kepl_to_cart;
%catches(TsunamiError) KeplerUtils::cart_to_kepl;

%rename(cart_to_kepl_deprecated) cart_to_kepl(double [6], double [6], double);
%rename(kepl_to_cart_deprecated) kepl_to_cart(double [6], double [6], double, double, double, double, double, double, double);
%feature("kwargs") cart_to_kepl;
%feature("kwargs") kepl_to_cart;
%feature("kwargs") cart_to_kepl2;
%feature("kwargs") kepl_to_cart2;

// add_particle_set
// override_position_and_velocities
%apply (double* IN_ARRAY1, int DIM1) {(double *mass_in, size_t nmass)}
%apply (long* IN_ARRAY1, int DIM1) {(long *stype_in, size_t nstype)}
%apply (double* IN_ARRAY1, int DIM1) {(double *rad_in, size_t nrad)}
%apply (double* IN_ARRAY2, int DIM1, int DIM2) {(double *pos_in, size_t npos, size_t pncoord)}
%apply (double* IN_ARRAY2, int DIM1, int DIM2) {(double *vel_in, size_t nvel, size_t vncoord)}
%apply (double* IN_ARRAY2, int DIM1, int DIM2) {(double *spin_in, size_t nspin, size_t sncoord)}

// sync_internal_state
%apply (double* INPLACE_ARRAY2, int DIM1, int DIM2) {(double *pos_inout, size_t npos, size_t pncoord)}
%apply (double* INPLACE_ARRAY2, int DIM1, int DIM2) {(double *vel_inout, size_t nvel, size_t vncoord)}
%apply (double* INPLACE_ARRAY2, int DIM1, int DIM2) {(double *spin_inout, size_t nspin, size_t sncoord)}

// initialize_tidal_parameters
%apply (double* IN_ARRAY1, int DIM1) {(double *kaps_in, size_t nkaps)}
%apply (double* IN_ARRAY1, int DIM1) {(double *taulag_in, size_t ntaulag)}
%apply (double* IN_ARRAY1, int DIM1) {(double *polyt_in, size_t npolyt)}
%apply (double* IN_ARRAY1, int DIM1) {(double *gyrad_in, size_t ngyrad)}

// get_chain_vectors
%apply (double** ARGOUTVIEWM_ARRAY2, int* DIM1, int* DIM2) {(double **chpos_out, int *npos, int *pncoord)}
%apply (double** ARGOUTVIEWM_ARRAY2, int* DIM1, int* DIM2) {(double **chvel_out, int *nvel, int *vncoord)}

// sync_masses, sync_radii, sync_eloss
%apply (double* INPLACE_ARRAY1, int DIM1) {(double *mass_out, size_t nmass)}
%apply (double* INPLACE_ARRAY1, int DIM1) {(double *rad_out, size_t nrad)}
%apply (double* INPLACE_ARRAY1, int DIM1) {(double *eloss_out, size_t neloss)}

// get_com, set_com
%apply (double ARGOUT_ARRAY1[ANY]) {(double pcom_out[3])}
%apply (double IN_ARRAY1[ANY]) {(double pcom_in[3])}
%apply (double ARGOUT_ARRAY1[ANY]) {(double vcom_in[3])}
%apply (double IN_ARRAY1[ANY]) {(double vcom_out[3])}

%apply size_t &OUTPUT {size_t &id1};
%apply size_t &OUTPUT {size_t &id2};

// get_acceleration
%apply (double ARGOUT_ARRAY1[ANY]) {(double acc_i_out[3])}
%apply (double ARGOUT_ARRAY1[ANY]) {(double acc_j_out[3])}

// waveform
%apply double &OUTPUT {(double &hplus)};
%apply double &OUTPUT {(double &hcross)};

// Member array
%typemap(out) size_t collind [ANY] {
  int i;
  $result = PyList_New($1_dim0);
  for (i = 0; i < $1_dim0; i++) {
    PyObject *o = PyLong_FromSize_t($1[i]);
    PyList_SetItem($result, i, o);
  }
}

/// keplfunc.h
%apply (double ARGOUT_ARRAY1[ANY]) {(double out_posvel[6])}
%apply (double IN_ARRAY1[ANY]) {(double prim_posvel[6])}
%apply (double ARGOUT_ARRAY1[ANY]) {(double out_keplpar[6])}
%apply (double IN_ARRAY1[ANY]) {(double in_posvel[6])}
%apply (double IN_ARRAY1[ANY]) {(double in_pos[3])}
%apply (double IN_ARRAY1[ANY]) {(double in_vel[3])}
%apply double &OUTPUT {double &a_corr};
%apply double &OUTPUT {double &e_corr};


%include "okinami.i"

%include "config.hpp"
%include "custom_types.hpp"
%feature("flatnested");
%include "Nbodyalgorithms.hpp"
%include "chain.hpp"
%import "leapfrog_stepped.hpp"
%import "bulirsch.hpp"
%import "classification.h"
%include "kozailidov.h"
%feature("flatnested", "");
%include "keplerutils.h"
%import "simprof.hpp"
%include "errhand.hpp"
%import "IO.h"
%include "tsunami.hpp"

%template(double3) chvec<double>;
%template(BSTable) BStable<TsunamiConfig::wSpins>;
%template(ChainSys) Chain<TsunamiConfig::wSpins, TsunamiConfig::useTDEtracker, TsunamiConfig::debug_ch>;
%template(LeapfrogStepped) Leapfrog_stepped<TsunamiConfig::wSpins, TsunamiConfig::debug_lf, TsunamiConfig::useProfiling>;
%template(BSExtrap) BSExtrapolator<TsunamiConfig::wSpins, TsunamiConfig::debug_bs>;
%template(TsunamiCode) TsunamiClass<TsunamiConfig::useProfiling, (TsunamiConfig::debug_bs or TsunamiConfig::debug_lf or TsunamiConfig::debug_ch)>;
%template(paired) std::pair<double, double>;
%template(pairlist) std::vector<std::pair<double, double>>;
