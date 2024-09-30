//
// Created by lex on 6/13/22.
//

#ifndef TSUNAMI_DOP853_HPP
#define TSUNAMI_DOP853_HPP

#include <functional>
#include <array>
#include <cstdlib>
#include <cmath>
#include <iostream>

template<typename T, size_t dim>
class DOP853 {

public:
    DOP853() = default;

    DOP853(std::function<void(const T &, T &, const double)> calc_derivative_func) {
        system = calc_derivative_func;
    }

    void bind_system(std::function<void(const T &, T &, const double)> calc_derivative_func) {
        system = calc_derivative_func;
    }

    void return_derivative(const T & y,T & dydt, const double t) {
        system(y, dydt, t);
    }

    void iterate_step(T &y, const double t, double &dt, double &dt_done) {
        bool reject = true;
        int counter = 0;
        while (reject and counter < counter_max) {
            reject = not try_step(y, t, dt, dt_done);
            counter++;
        }
    }

    bool try_step(T &y, const double t, double &dt, double &dt_done) {
        T y_new;

        system(y, k1, t);

        for (size_t i = 0; i < dim; i++)
            y_new[i] = y[i] + dt * a21 * k1[i];
        system(y_new, k2, t + c2 * dt);

        for (size_t i = 0; i < dim; i++)
            y_new[i] = y[i] + dt * (a31*k1[i] + a32*k2[i]);
        system(y_new, k3, t + c3 * dt);

        for (size_t i = 0; i < dim; i++)
            y_new[i] = y[i] + dt * (a41*k1[i] + a43*k3[i]);
        system(y_new, k4, t + c4 * dt);

        for (size_t i = 0; i < dim; i++)
            y_new[i] = y[i] + dt * (a51*k1[i] + a53*k3[i] + a54*k4[i]);
        system(y_new, k5, t + c5 * dt);

        for (size_t i = 0; i < dim; i++)
            y_new[i] = y[i] + dt * (a61*k1[i] + a64*k4[i] + a65*k5[i]);
        system(y_new, k6, t + c6 * dt);

        for (size_t i = 0; i < dim; i++)
            y_new[i] = y[i] + dt * (a71*k1[i] + a74*k4[i] + a75*k5[i] + a76*k6[i]);
        system(y_new, k7, t + c7 * dt);

        for (size_t i = 0; i < dim; i++)
            y_new[i] = y[i] + dt * (a81*k1[i] + a84*k4[i] + a85*k5[i] + a86*k6[i] +
                                    a87*k7[i]);
        system(y_new, k8, t + c8 * dt);

        for (size_t i = 0; i < dim; i++)
            y_new[i] = y[i] + dt * (a91*k1[i] + a94*k4[i] + a95*k5[i] + a96*k6[i] +
                                    a97*k7[i] + a98*k8[i]);
        system(y_new, k9, t + c9 * dt);

        for (size_t i = 0; i < dim; i++)
            y_new[i] = y[i] + dt * (a101*k1[i] + a104*k4[i] + a105*k5[i] + a106*k6[i] +
                                    a107*k7[i] + a108*k8[i] + a109*k9[i]);
        system(y_new, k10, t + c10 * dt);

        for (size_t i = 0; i < dim; i++)
            y_new[i] = y[i] + dt * (a111*k1[i] + a114*k4[i] + a115*k5[i] + a116*k6[i] +
                                    a117*k7[i] + a118*k8[i] + a119*k9[i] + a1110*k10[i]);
        system(y_new, k2, t + c11 * dt); // Re-using k2

        for (size_t i = 0; i < dim; i++)
            y_new[i] = y[i] + dt * (a121*k1[i] + a124*k4[i] + a125*k5[i] + a126*k6[i] +
                                    a127*k7[i] + a128*k8[i] + a129*k9[i] +
                                    a1210*k10[i] + a1211*k2[i]);
        system(y_new, k3, t + dt); // Re-using k3

        // Final solution, re-using k4
        for (size_t i = 0; i < dim; i++) {
            k4[i] = (b1*k1[i] + b6*k6[i] + b7*k7[i] + b8*k8[i] + b9*k9[i] +
                     b10*k10[i] + b11*k2[i] + b12*k3[i]);
            y_new[i] = y[i] + dt * k4[i];
        }

        //Estimating error, 3rd order and 5th order
        double yerr2_3 = 0.0, yerr2_5 = 0.0;
        for (size_t i = 0; i < dim; i++) {
            double sk = abs_err + rel_err * (ax * std::max(std::fabs(y[i]), std::fabs(y_new[i]))
                                            + adxdt * dt * std::fabs(k1[i]));
            double sqr = (k4[i] - bhh1*k1[i] - bhh2*k9[i] - bhh3*k3[i]) / sk;
            yerr2_3 += sqr*sqr;
            sqr = (er1*k1[i] + er6*k6[i] + er7*k7[i] + er8*k8[i] + er9*k9[i] +
                   er10 * k10[i] + er11*k2[i] + er12*k3[i]) / sk;
            yerr2_5 += sqr*sqr;
        }

        // Estimating 8th order error (Eq. 10.17 Hairer I)
        double deno = yerr2_5 + 0.01 * yerr2_3;
        if (deno <= 0.0) deno = 1.0; // Avoid singularity
        double yerr = std::fabs(dt) * yerr2_5 * sqrt (1.0 / (deno*static_cast<double>(dim)));
        //std::cout << "yerr: " << yerr << std::endl;

        // Compute new timestep dt_new
        double fac11 = pow(yerr, error_exp);
        // PI-control, (Eq. 2.46 Hairer II)
        double fac = fac11 / pow(facold,beta);
        // 0.333 <= dt_new/dt <= 6.0
        fac = std::max(facc2, std::min(facc1, fac/safe));
        //std::cout << "dt_fac: " << fac << std::endl;

        double dt_new = dt / fac;

        if (yerr <= 1.0) {
            // Accept step
            facold = std::max(yerr, 1.0e-4);
            //system(y_new, k4, t + dt);

            //if (std::fabs(dt_new) > dt_max)
            //    dt_new = dt_max;

            // If previously was rejected, do not increase timestep
            if (reject) dt_new = std::min(dt_new, dt);
            reject = false;
            y = y_new;
            dt_done = dt;
            dt = dt_new;
            return true;
        } else {
            // Reject step
            dt_new = dt / std::min(facc1, fac11/safe);

            dt = dt_new;
            reject = true;
            return false;
        }
    }

    void estimate_first_timestep(T &y, const double t, double &dt) {
        system(y, k1, t);

        double dnf = 0.0;
        double dny = 0.0;
        for (size_t i = 0; i < dim; i++) {
            double sk = abs_err + rel_err * std::fabs(y[i]);
            double sqr = k1[i] / sk;
            dnf += sqr*sqr;
            sqr = y[i] / sk;
            dny += sqr*sqr;
        }

        if (dnf <= 1e-10 or dny <= 1e-10) dt = 1e-6;
        else dt = sqrt(dny/dnf) * 0.01;

        // Euler step
        T y1;
        for (size_t i = 0; i < dim; i++)
            y1[i] = y[i] + dt * k1[i];

        // Second step
        system(y1, k2, t+dt);

        // Second derivative
        double der2 = 0.0;
        for (size_t i = 0; i < dim; i++)
        {
            double sk = abs_err + rel_err * std::fabs(y[i]);
            double sqr = (k2[i] - k1[i]) / sk;
            der2 += sqr*sqr;
        }
        der2 = sqrt(der2) / dt;

        // Step size is computed such that dt^8 * max(norm(k1),norm(der2)) = 0.01
        double der12 = std::max(std::fabs(der2), sqrt(dnf));
        double dt1;
        if (der12 <= 1e-15) dt1 = std::max(1e-6, std::fabs(dt)*1e-3);
        else dt1 = pow(0.01/der12, order_exp);

        dt = std::min(100.0 * std::fabs(dt), dt1);
    }

    bool reject = false;

    double rel_err = 1e-7;
    double abs_err = 1e-8;
    double ax = 1.0;
    double adxdt = 1.0;

    double facc1 = 1.0/0.333;
    double facc2 = 1.0/6.0;
    double safe = 0.9;

    void set_beta(double beta_new) {
        beta = beta_new;
        error_exp = order_exp - beta * 0.2;
    }

private:
    //void (*system)(const T &y_const, T &dydt_dev, const double t);
    std::function<void(const T &, T &, const double)> system;
    T k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13;

    T fabs(T vec) {
        T out_vec;
        for (size_t i = 0; i < dim; i++) out_vec[i] = std::fabs(vec[i]);
        return out_vec;
    }

    double max(T vec) {
        double max_value = vec[0];
        for (size_t i = 1; i < dim; i++) max_value = vec[i] = (vec[i] > max_value) ? vec[i] : max_value;
        return max_value;
    }

    static constexpr double order_exp = 1.0 / 8.0;
    double beta = 0.0; // Timestep stabilization
    double error_exp = order_exp;

    double facold = 1e-4;

    int counter_max = 1000;

    static constexpr double c2 = 0.526001519587677318785587544488E-01;
    static constexpr double c3 = 0.789002279381515978178381316732E-01;
    static constexpr double c4  = 0.118350341907227396726757197510E+00;
    static constexpr double c5  = 0.281649658092772603273242802490E+00;
    static constexpr double c6  = 0.333333333333333333333333333333E+00;
    static constexpr double c7  = 0.25E+00;
    static constexpr double c8  = 0.307692307692307692307692307692E+00;
    static constexpr double c9  = 0.651282051282051282051282051282E+00;
    static constexpr double c10 = 0.6E+00;
    static constexpr double c11 = 0.857142857142857142857142857142E+00;
    static constexpr double c14 = 0.1E+00;
    static constexpr double c15 = 0.2E+00;
    static constexpr double c16 = 0.777777777777777777777777777778E+00;

    static constexpr double b1 =   5.42937341165687622380535766363E-2;
    static constexpr double b6 =   4.45031289275240888144113950566E0;
    static constexpr double b7 =   1.89151789931450038304281599044E0;
    static constexpr double b8 =  -5.8012039600105847814672114227E0;
    static constexpr double b9 =   3.1116436695781989440891606237E-1;
    static constexpr double b10 = -1.52160949662516078556178806805E-1;
    static constexpr double b11 =  2.01365400804030348374776537501E-1;
    static constexpr double b12 =  4.47106157277725905176885569043E-2;

    static constexpr double bhh1 = 0.244094488188976377952755905512E+00;
    static constexpr double bhh2 = 0.733846688281611857341361741547E+00;
    static constexpr double bhh3 = 0.220588235294117647058823529412E-01;

    static constexpr double er1  =  0.1312004499419488073250102996E-01;
    static constexpr double er6  = -0.1225156446376204440720569753E+01;
    static constexpr double er7  = -0.4957589496572501915214079952E+00;
    static constexpr double er8  =  0.1664377182454986536961530415E+01;
    static constexpr double er9  = -0.3503288487499736816886487290E+00;
    static constexpr double er10 =  0.3341791187130174790297318841E+00;
    static constexpr double er11 =  0.8192320648511571246570742613E-01;
    static constexpr double er12 = -0.2235530786388629525884427845E-01;

    static constexpr double a21 =    5.26001519587677318785587544488E-2;
    static constexpr double a31 =    1.97250569845378994544595329183E-2;
    static constexpr double a32 =    5.91751709536136983633785987549E-2;
    static constexpr double a41 =    2.95875854768068491816892993775E-2;
    static constexpr double a43 =    8.87627564304205475450678981324E-2;
    static constexpr double a51 =    2.41365134159266685502369798665E-1;
    static constexpr double a53 =   -8.84549479328286085344864962717E-1;
    static constexpr double a54 =    9.24834003261792003115737966543E-1;
    static constexpr double a61 =    3.7037037037037037037037037037E-2;
    static constexpr double a64 =    1.70828608729473871279604482173E-1;
    static constexpr double a65 =    1.25467687566822425016691814123E-1;
    static constexpr double a71 =    3.7109375E-2;
    static constexpr double a74 =    1.70252211019544039314978060272E-1;
    static constexpr double a75 =    6.02165389804559606850219397283E-2;
    static constexpr double a76 =   -1.7578125E-2;

    static constexpr double a81 =    3.70920001185047927108779319836E-2;
    static constexpr double a84 =    1.70383925712239993810214054705E-1;
    static constexpr double a85 =    1.07262030446373284651809199168E-1;
    static constexpr double a86 =   -1.53194377486244017527936158236E-2;
    static constexpr double a87 =    8.27378916381402288758473766002E-3;
    static constexpr double a91 =    6.24110958716075717114429577812E-1;
    static constexpr double a94 =   -3.36089262944694129406857109825E0;
    static constexpr double a95 =   -8.68219346841726006818189891453E-1;
    static constexpr double a96 =    2.75920996994467083049415600797E1;
    static constexpr double a97 =    2.01540675504778934086186788979E1;
    static constexpr double a98 =   -4.34898841810699588477366255144E1;
    static constexpr double a101 =   4.77662536438264365890433908527E-1;
    static constexpr double a104 =  -2.48811461997166764192642586468E0;
    static constexpr double a105 =  -5.90290826836842996371446475743E-1;
    static constexpr double a106 =   2.12300514481811942347288949897E1;
    static constexpr double a107 =   1.52792336328824235832596922938E1;
    static constexpr double a108 =  -3.32882109689848629194453265587E1;
    static constexpr double a109 =  -2.03312017085086261358222928593E-2;

    static constexpr double a111 =  -9.3714243008598732571704021658E-1;
    static constexpr double a114 =   5.18637242884406370830023853209E0;
    static constexpr double a115 =   1.09143734899672957818500254654E0;
    static constexpr double a116 =  -8.14978701074692612513997267357E0;
    static constexpr double a117 =  -1.85200656599969598641566180701E1;
    static constexpr double a118 =   2.27394870993505042818970056734E1;
    static constexpr double a119 =   2.49360555267965238987089396762E0;
    static constexpr double a1110 = -3.0467644718982195003823669022E0;
    static constexpr double a121 =   2.27331014751653820792359768449E0;
    static constexpr double a124 =  -1.05344954667372501984066689879E1;
    static constexpr double a125 =  -2.00087205822486249909675718444E0;
    static constexpr double a126 =  -1.79589318631187989172765950534E1;
    static constexpr double a127 =   2.79488845294199600508499808837E1;
    static constexpr double a128 =  -2.85899827713502369474065508674E0;
    static constexpr double a129 =  -8.87285693353062954433549289258E0;
    static constexpr double a1210 =  1.23605671757943030647266201528E1;
    static constexpr double a1211 =  6.43392746015763530355970484046E-1;

    static constexpr double a141 =  5.61675022830479523392909219681E-2;
    static constexpr double a147 =  2.53500210216624811088794765333E-1;
    static constexpr double a148 = -2.46239037470802489917441475441E-1;
    static constexpr double a149 = -1.24191423263816360469010140626E-1;
    static constexpr double a1410 =  1.5329179827876569731206322685E-1;
    static constexpr double a1411 =  8.20105229563468988491666602057E-3;
    static constexpr double a1412 =  7.56789766054569976138603589584E-3;
    static constexpr double a1413 = -8.298E-3;

    static constexpr double a151 =  3.18346481635021405060768473261E-2;
    static constexpr double a156 =  2.83009096723667755288322961402E-2;
    static constexpr double a157 =  5.35419883074385676223797384372E-2;
    static constexpr double a158 = -5.49237485713909884646569340306E-2;
    static constexpr double a1511 = -1.08347328697249322858509316994E-4;
    static constexpr double a1512 =  3.82571090835658412954920192323E-4;
    static constexpr double a1513 = -3.40465008687404560802977114492E-4;
    static constexpr double a1514 =  1.41312443674632500278074618366E-1;
    static constexpr double a161 = -4.28896301583791923408573538692E-1;
    static constexpr double a166 = -4.69762141536116384314449447206E0;
    static constexpr double a167 =  7.68342119606259904184240953878E0;
    static constexpr double a168 =  4.06898981839711007970213554331E0;
    static constexpr double a169 =  3.56727187455281109270669543021E-1;
    static constexpr double a1613 = -1.39902416515901462129418009734E-3;
    static constexpr double a1614 =  2.9475147891527723389556272149E0;
    static constexpr double a1615 = -9.15095847217987001081870187138E0;

    static constexpr double d41  = -0.84289382761090128651353491142E+01;
    static constexpr double d46  =  0.56671495351937776962531783590E+00;
    static constexpr double d47  = -0.30689499459498916912797304727E+01;
    static constexpr double d48  =  0.23846676565120698287728149680E+01;
    static constexpr double d49  =  0.21170345824450282767155149946E+01;
    static constexpr double d410 = -0.87139158377797299206789907490E+00;
    static constexpr double d411 =  0.22404374302607882758541771650E+01;
    static constexpr double d412 =  0.63157877876946881815570249290E+00;
    static constexpr double d413 = -0.88990336451333310820698117400E-01;
    static constexpr double d414 =  0.18148505520854727256656404962E+02;
    static constexpr double d415 = -0.91946323924783554000451984436E+01;
    static constexpr double d416 = -0.44360363875948939664310572000E+01;

    static constexpr double d51  =  0.10427508642579134603413151009E+02;
    static constexpr double d56  =  0.24228349177525818288430175319E+03;
    static constexpr double d57  =  0.16520045171727028198505394887E+03;
    static constexpr double d58  = -0.37454675472269020279518312152E+03;
    static constexpr double d59  = -0.22113666853125306036270938578E+02;
    static constexpr double d510 =  0.77334326684722638389603898808E+01;
    static constexpr double d511 = -0.30674084731089398182061213626E+02;
    static constexpr double d512 = -0.93321305264302278729567221706E+01;
    static constexpr double d513 =  0.15697238121770843886131091075E+02;
    static constexpr double d514 = -0.31139403219565177677282850411E+02;
    static constexpr double d515 = -0.93529243588444783865713862664E+01;
    static constexpr double d516 =  0.35816841486394083752465898540E+02;

    static constexpr double d61 =  0.19985053242002433820987653617E+02;
    static constexpr double d66 = -0.38703730874935176555105901742E+03;
    static constexpr double d67 = -0.18917813819516756882830838328E+03;
    static constexpr double d68 =  0.52780815920542364900561016686E+03;
    static constexpr double d69 = -0.11573902539959630126141871134E+02;
    static constexpr double d610 =  0.68812326946963000169666922661E+01;
    static constexpr double d611 = -0.10006050966910838403183860980E+01;
    static constexpr double d612 =  0.77771377980534432092869265740E+00;
    static constexpr double d613 = -0.27782057523535084065932004339E+01;
    static constexpr double d614 = -0.60196695231264120758267380846E+02;
    static constexpr double d615 =  0.84320405506677161018159903784E+02;
    static constexpr double d616 =  0.11992291136182789328035130030E+02;

    static constexpr double d71  = -0.25693933462703749003312586129E+02;
    static constexpr double d76  = -0.15418974869023643374053993627E+03;
    static constexpr double d77  = -0.23152937917604549567536039109E+03;
    static constexpr double d78  =  0.35763911791061412378285349910E+03;
    static constexpr double d79  =  0.93405324183624310003907691704E+02;
    static constexpr double d710 = -0.37458323136451633156875139351E+02;
    static constexpr double d711 =  0.10409964950896230045147246184E+03;
    static constexpr double d712 =  0.29840293426660503123344363579E+02;
    static constexpr double d713 = -0.43533456590011143754432175058E+02;
    static constexpr double d714 =  0.96324553959188282948394950600E+02;
    static constexpr double d715 = -0.39177261675615439165231486172E+02;
    static constexpr double d716 = -0.14972683625798562581422125276E+03;

};


#endif //TSUNAMI_DOP853_HPP
