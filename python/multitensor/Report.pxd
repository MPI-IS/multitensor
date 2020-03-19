from libcpp.vector cimport vector


cdef extern from "multitensor/utils.hpp" namespace "multitensor::utils":
    cdef cppclass Report:
        int nof_realizations
        vector[double] vec_L2
        vector[size_t] vec_iter;
        vector[const char *] vec_term_reason
        double duration

        Report() except +
        double max_L2()
