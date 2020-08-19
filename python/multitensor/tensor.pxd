from libcpp.vector cimport vector


cdef extern from "multitensor/tensor.hpp" namespace "multitensor::tensor":
    cdef cppclass Tensor[scalar_t]:
        size_t nrows
        size_t ncols
        size_t ntubes
        vector[scalar_t] data

        Tensor() except +
        Tensor(size_t nrows, size_t ncols, size_t ntubes = 1) except +

        void resize(size_t nrows_, size_t ncols_, size_t ntubes)
        size_t get_nrows()
        size_t get_ncols()
        size_t get_ntubes()

        scalar_t &operator()(const size_t i, const size_t j, const size_t k)

    cdef cppclass Matrix[scalar_t](Tensor[scalar_t]):
        Matrix() except +
        Matrix(size_t nrows, size_t ncols) except +

        void resize(size_t nrows, size_t ncols)
        scalar_t &operator()(const size_t i, const size_t j)

    cdef cppclass SymmetricTensor[scalar_t](Tensor[scalar_t]):
        SymmetricTensor() except +
        SymmetricTensor(size_t nrows, size_t ntubes) except +

        void resize(size_t nrwos, size_t ntubes)

    cdef cppclass DiagonalTensor[scalar_t](Tensor[scalar_t]):
        DiagonalTensor() except +
        DiagonalTensor(size_t nrows, size_t ntubes)
