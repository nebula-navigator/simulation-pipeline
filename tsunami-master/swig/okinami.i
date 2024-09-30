%module tsunami

%pythoncode %{
from enum import IntEnum
stateid = IntEnum("id", "e1 e2 g1 g2 h1 a1 H a2 m1 m2 m3 R1 R2 R3", start=0)
%}


%immutable Okinami::y_stop;
%naturalvar Okinami::y_stop;
%ignore Okinami::initialize_integrator;
%ignore Okinami::id;
%immutable Okinami::speed_of_light;

%catches(TsunamiError) Okinami::evolve_system;

%feature("kwargs") set_integrator_tolerance;

%template(delaunay_state) std::array<double,14>;

%apply double &OUTPUT { double &inc1};
%apply double &OUTPUT { double &inc2};
%apply double &OUTPUT { double &H};
%apply double &INOUT { double &t};

// compute_derivatives
%typemap(argout, optimal="1") Okinami::delaunay_state &dydt_dev {
	size_t size = $1->size();
	npy_intp dims[] = { static_cast<npy_intp>(size) };

	PyObject* out_array = PyArray_SimpleNew(1, dims, NPY_DOUBLE);

	if (!out_array) {
		PyErr_SetString(PyExc_ValueError, "Unable to create the output array.");
		SWIG_fail;
	}
	double* out_data = (double*) array_data(out_array);
	for (size_t i = 0; i < size; i++) {
		out_data[i] = static_cast<double>($1->data()[i]);
	}
	$result = out_array;
}
%typemap(in,numinputs=0) Okinami::delaunay_state &dydt_dev (Okinami::delaunay_state temp_dydt) {
	$1 = &temp_dydt;
}
%typemap(in,numinputs=0) const double t (double temp_t) {
	$1 = temp_t;
}

%apply double &INOUT { double &time};
%apply double &INOUT { double &dt};

%typemap(out, optimal="1") double &time {
    $result = PyFloat_FromDouble($1);
}
%typemap(out, optimal="1") double &dt {
    $result = PyFloat_FromDouble($1);
}

%typemap(in) Okinami::delaunay_state &y {
    PyArrayObject *numpy_array = obj_to_array_no_conversion($input, NPY_DOUBLE);
    npy_intp size[1] = {static_cast<long>($1->size())};
    if (!numpy_array || !require_dimensions(numpy_array, 1) || !require_size(numpy_array, size, 1)) SWIG_fail;
    double *pyarr = (double*) array_data(numpy_array);

    $1 = reinterpret_cast<Okinami::delaunay_state *>(pyarr);
}
