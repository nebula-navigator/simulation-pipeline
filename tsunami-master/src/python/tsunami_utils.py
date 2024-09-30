try:
    from amuse.units import units, constants

    AMUSE_FOUND = True
except ModuleNotFoundError:
    AMUSE_FOUND = False
import numpy as np
from numpy.lib import recfunctions

ic_type = np.dtype({
    'names': ['x', 'y', 'z', 'vx', 'vy', 'vz', 'mass', 'eps', 'radius', 'type'],
    'formats': ['f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'i8'] })

output_type = np.dtype({ 'names': ['x', 'y', 'z', 'vx', 'vy', 'vz', 'mass'],
                         'formats': ['f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8'] })

energy_type = np.dtype({ 'names': ['time', 'dE', 'E'],
                         'formats': ['f8', 'f8', 'f8'] })


def load_tsunami_input(fname):
    ictable = np.genfromtxt(fname, dtype=ic_type)
    return ictable


def load_tsunami_output(fname):
    outtable = np.genfromtxt(fname, dtype=output_type)
    return outtable


def load_tsunami_output_at_time(outname, time_nbody, N=None, as_input_table=False):
    icname = outname.replace("_output.dat", ".dat")
    enename = outname.replace("_output.dat", "_energy.dat")

    t = np.loadtxt(enename, usecols=(0,))
    if time_nbody > t[-1]:
        raise ValueError(
            "Total time is shorter than requested one\nRequested: {:g}\nFinal: {:g}".format(time_nbody, t[-1]))

    Ntime = len(t)
    time_ind = np.argwhere(t > time_nbody)[0][0]

    if N is None:
        ictab = np.loadtxt(icname)
        N = len(ictab)

    with open(outname, 'r') as f:
        outlines = f.readlines()
    Nout = len(outlines)
    if Ntime * N != Nout:
        raise ValueError("{:s} and {:s} do not match\n".format(enename, outname))

    shnap = np.genfromtxt(outlines[time_ind * N:(time_ind + 1) * N], dtype=output_type)

    if as_input_table is True:
        ictab = np.loadtxt(icname, dtype=np.dtype({'names': ['eps', 'radius', 'type'], 'formats': ['f8', 'f8', 'i8']}),
                           usecols=(7, 8, 9))
        for nam, dtyp in ictab.dtype.descr:
            shnap = recfunctions.append_fields(shnap, nam, ictab[nam], dtyp, usemask=False)
    return shnap


def write_tsunami_input(ictable, filename):
    """
    ictable is a named numpy array. named field must be:
    x
    y
    z
    vx
    vy
    vz
    mass
    radius [optional, defaults to zero. not knowing the units, not possible to know the speed of light]
    type [optional, defaults to zero]
    eps [optional, defaults to zero]
    """

    namelist = ictable.dtype.names
    if "eps" not in namelist:
        ictable = recfunctions.append_fields(ictable, "eps", np.zeros(N), 'f8', usemask=False)
    if "type" not in namelist:
        ictable = recfunctions.append_fields(ictable, "type", np.zeros(N), 'i8', usemask=False)
    if "radius" not in namelist:
        ictable = recfunctions.append_fields(ictable, "radius", np.zeros(N), 'f8', usemask=False)

    np.savetxt(filename, np.transpose([ictable["x"],
                                       ictable["y"],
                                       ictable["z"],
                                       ictable["vx"],
                                       ictable["vy"],
                                       ictable["vz"],
                                       ictable["mass"],
                                       ictable["eps"],
                                       ictable["radius"],
                                       ictable["type"]]),
               fmt='%.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.0f')


def consistency_check(icname):
    enename = icname.replace(".dat", "_energy.dat")
    outname = enename.replace("_energy.dat", "_output.dat")

    t = np.loadtxt(enename, usecols=(0,))
    if t.ndim == 1:
        Ntime = len(t)
    else: Ntime = 1

    ictab = np.loadtxt(icname, usecols=(0,))
    N = len(ictab)

    outtab = np.loadtxt(outname, usecols=(0,))
    Nout = len(outtab)

    if Ntime * N == Nout:
        return False
    else:
        return True


if AMUSE_FOUND:

    def write_tsunami_input_from_amuse_set(syst, filename, mass_unit=units.MSun, length_unit=units.parsec,
                                           silent=False):
        """
        syst is AMUSE particle dataset with:
        x
        y
        z
        vx
        vy
        vz
        mass
        radius [optional, defaults to Schwartzchild radius]
        type [optional, defaults to zero]
        eps [optional, defaults to zero]
        """
        conv = nbody_system.nbody_to_si(mass_unit, length_unit)
        velocity_unit = conv.to_si(nbody_system.length / nbody_system.time).as_unit()

        attlist = syst.get_attribute_names_defined_in_store()
        if ("eps" not in attlist):
            syst.eps = 0 | units.au
        if ("type" not in attlist):
            syst.type = 0
        if ("radius" not in attlist):
            for p in syst:
                p.radius = 2 * p.mass * constants.G / constants.c ** 2

        np.savetxt(filename, np.transpose([syst.x.value_in(length_unit),
                                           syst.y.value_in(length_unit),
                                           syst.z.value_in(length_unit),
                                           syst.vx.value_in(velocity_unit),
                                           syst.vy.value_in(velocity_unit),
                                           syst.vz.value_in(velocity_unit),
                                           syst.mass.value_in(mass_unit),
                                           syst.eps.value_in(length_unit),
                                           syst.radius.value_in(length_unit),
                                           syst.type]),
                   fmt='%.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.0f')
        if not silent:
            print("written {:s}\n".format(filename))
            print("1 time unit = {:s}".format(conv.to_si(nbody_system.time).as_string_in(units.Myr)))
            print("1 yr = {:g}".format(conv.to_nbody(1 | units.yr).value_in(nbody_system.time)))
