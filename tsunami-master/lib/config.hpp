//
// Created by lex on 8/29/21.
//

#ifndef TSUNAMI_CONFIG_HPP
#define TSUNAMI_CONFIG_HPP


struct TsunamiConfig {
    // Runtime options
    bool wPNs = false;
    bool wEqTides = false;
    bool wDynTides = false;
    bool wExt = false;
    bool wExt_vdep = false;
    bool wMassEvol = false;

    double Mscale = 1; ///< Mass unit scale in MSun
    double Lscale = 1; ///< Length unit scale in au

    // Regularization option
    double alpha = 1;
    double beta = 0;
    double gamma = 0;
    double Usafe = 1e-3;
    double dcoll = 1;
    bool TTL = false;

    bool pn1 = true;
    bool pn2 = true;
    bool pn25 = true;
    bool pn3 = true;
    bool pn35 = true;

    bool ss = true; ///< spin-spin coupling
    bool so = true; ///< spin-orbit coupling

    // Compile time options, it is ugly, but I don't know any better syntax
#ifdef ABSORB
    staticconstexpr bool useChaosAbsorb = true;
#else
    static constexpr bool useChaosAbsorb = false;
#endif

#ifdef SPINS
    static constexpr bool wSpins = true;
#else
    static constexpr bool wSpins = false;
#endif

#ifdef DEBUGLF
    static constexpr bool debug_lf = true;
#else
    static constexpr bool debug_lf = false;
#endif

#ifdef DEBUGBS
    static constexpr bool debug_bs = true;
#else
    static constexpr bool debug_bs = false;
#endif

#ifdef DEBUGCH
    static constexpr bool debug_ch = true;
#else
    static constexpr bool debug_ch = false;
#endif

#ifdef PROFILE
    static constexpr bool useProfiling = true;
#else
    static constexpr bool useProfiling = false;
#endif

#ifdef TIMINGU
    static constexpr bool useTiming = true;
#else
    static constexpr bool useTiming = false;
#endif

#ifdef TDETRACKER
    static constexpr bool useTDEtracker = true;
#else
    static constexpr bool useTDEtracker = false;
#endif

};

#endif //TSUNAMI_CONFIG_HPP
