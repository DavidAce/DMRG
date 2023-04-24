#pragma once
#include "general/sfinae.h"
#include "mpi-logger.h"
#include <cassert>
#include <complex>
#include <h5pp/details/h5ppFilesystem.h>
#include <h5pp/details/h5ppFormat.h>
#include <mpi.h>
#include <vector>
namespace mpi {

    inline bool on = false;

    struct comm {
        int id   = 0;
        int size = 1;

        template<typename T>
        T get_id() {
            return static_cast<T>(id);
        }
        template<typename T>
        T get_size() {
            return static_cast<T>(size);
        }
    };

    inline comm world;
    void        init();
    void        finalize();
    void        barrier();
    template<typename T>
    [[nodiscard]] constexpr MPI_Datatype get_dtype() noexcept {
        MPI_Datatype mpi_type = MPI_DATATYPE_NULL;
        using D               = typename std::decay<T>::type;
        if constexpr(std::is_same_v<T, char>) mpi_type = MPI_CHAR;
        if constexpr(std::is_same_v<D, signed char>) mpi_type = MPI_SIGNED_CHAR;
        if constexpr(std::is_same_v<D, unsigned char>) mpi_type = MPI_UNSIGNED_CHAR;
        if constexpr(std::is_same_v<D, wchar_t>) mpi_type = MPI_WCHAR;
        if constexpr(std::is_same_v<D, signed short>) mpi_type = MPI_SHORT;
        if constexpr(std::is_same_v<D, unsigned short>) mpi_type = MPI_UNSIGNED_SHORT;
        if constexpr(std::is_same_v<D, signed int>) mpi_type = MPI_INT;
        if constexpr(std::is_same_v<D, unsigned int>) mpi_type = MPI_UNSIGNED;
        if constexpr(std::is_same_v<D, signed long int>) mpi_type = MPI_LONG;
        if constexpr(std::is_same_v<D, unsigned long int>) mpi_type = MPI_UNSIGNED_LONG;
        if constexpr(std::is_same_v<D, signed long long int>) mpi_type = MPI_LONG_LONG;
        if constexpr(std::is_same_v<D, unsigned long long int>) mpi_type = MPI_UNSIGNED_LONG_LONG;
        if constexpr(std::is_same_v<D, float>) mpi_type = MPI_FLOAT;
        if constexpr(std::is_same_v<D, double>) mpi_type = MPI_DOUBLE;
        if constexpr(std::is_same_v<D, long double>) mpi_type = MPI_LONG_DOUBLE;
        if constexpr(std::is_same_v<D, std::int8_t>) mpi_type = MPI_INT8_T;
        if constexpr(std::is_same_v<D, std::int16_t>) mpi_type = MPI_INT16_T;
        if constexpr(std::is_same_v<D, std::int32_t>) mpi_type = MPI_INT32_T;
        if constexpr(std::is_same_v<D, std::int64_t>) mpi_type = MPI_INT64_T;
        if constexpr(std::is_same_v<D, std::uint8_t>) mpi_type = MPI_UINT8_T;
        if constexpr(std::is_same_v<D, std::uint16_t>) mpi_type = MPI_UINT16_T;
        if constexpr(std::is_same_v<D, std::uint32_t>) mpi_type = MPI_UINT32_T;
        if constexpr(std::is_same_v<D, std::uint64_t>) mpi_type = MPI_UINT64_T;
        if constexpr(std::is_same_v<D, bool>) mpi_type = MPI_C_BOOL;
        if constexpr(std::is_same_v<D, std::complex<float>>) mpi_type = MPI_C_COMPLEX;
        if constexpr(std::is_same_v<D, std::complex<double>>) mpi_type = MPI_C_DOUBLE_COMPLEX;
        if constexpr(std::is_same_v<D, std::complex<long double>>) mpi_type = MPI_C_LONG_DOUBLE_COMPLEX;
        if constexpr(sfinae::has_Scalar_v<D>) return get_dtype<typename D::Scalar>();
        if constexpr(sfinae::has_value_type_v<D>) return get_dtype<typename D::value_type>();
        assert(mpi_type != MPI_DATATYPE_NULL);
        return mpi_type;
    }
    template<typename T>
    [[nodiscard]] void *get_buffer(T &data) {
        if constexpr(sfinae::has_data_v<T>)
            return static_cast<void *>(data.data());
        else if constexpr(std::is_pointer_v<T> or std::is_array_v<T>)
            return static_cast<void *>(data);
        else
            return static_cast<void *>(&data);
    }

    template<typename T>
    [[nodiscard]] const void *get_cbuffer(const T &data) {
        if constexpr(sfinae::has_data_v<T>)
            return static_cast<const void *>(data.data());
        else if constexpr(std::is_pointer_v<T> or std::is_array_v<T>)
            return static_cast<const void *>(data);
        else
            return static_cast<const void *>(&data);
    }

    template<typename T>
    [[nodiscard]] int get_count(T &data) {
        if constexpr(sfinae::has_size_v<T> or std::is_array_v<T>)
            return static_cast<int>(std::size(data));
        else
            return 1;
    }

    template<typename T>
    void send(const T &data, int dst, int tag) {
        if constexpr(sfinae::has_size_v<T>) {
            size_t count = data.size();
            MPI_Send(&count, 1, mpi::get_dtype<size_t>(), dst, 0, MPI_COMM_WORLD);
        }

        MPI_Send(mpi::get_cbuffer(data), mpi::get_count(data), mpi::get_dtype<T>(), dst, tag, MPI_COMM_WORLD);
    }

    template<typename T>
    void recv(T &data, int src, int tag) {
        if constexpr(sfinae::has_size_v<T>) {
            size_t count;
            MPI_Recv(&count, 1, mpi::get_dtype<size_t>(), src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            if constexpr(sfinae::has_resize_v<T>) data.resize(count);
            if(data.size() < count) throw std::runtime_error(fmt::format("mpi::recv: cointainer size {} < count {}", data.size(), count));
        }
        MPI_Recv(mpi::get_buffer(data), mpi::get_count(data), mpi::get_dtype<T>(), src, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    template<typename T>
    void sendrecv(const T &send, const T &recv, int src, int dst, int tag) {
        // Start by sending the data size, so we can resize the receiving buffer accordingly
        int count = mpi::get_count(send);
        MPI_Sendrecv_replace(&count, 1, MPI_INT, dst, tag, src, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if constexpr(sfinae::has_resize_v<T>) {
            recv.resize(count); // Both containers are now ready to receive
        }
        MPI_Sendrecv(mpi::get_buffer(send), mpi::get_count(send), mpi::get_dtype<T>(), dst, tag, mpi::get_buffer(recv), mpi::get_count(recv),
                     mpi::get_dtype<T>(), src, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    template<typename T>
    void sendrecv_replace(const T &data, int src, int dst, int tag) {
        // Start by sending the data size, so we can resize the receiving buffer accordingly
        int count = mpi::get_count(data);
        int err1  = MPI_Sendrecv_replace(&count, 1, MPI_INT, dst, tag, src, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if(err1 != 0) throw std::runtime_error(fmt::format("MPI_Sendrecv_replace error 1: {}", err1));

        if constexpr(sfinae::has_resize_v<T>) {
            data.resize(count); // Should not modify the src container
        }

        int err2 =
            MPI_Sendrecv_replace(mpi::get_buffer(data), mpi::get_count(data), mpi::get_dtype<T>(), dst, tag, src, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if(err2 != 0) throw std::runtime_error(fmt::format("MPI_Sendrecv_replace error 2: {}", err2));
    }

    template<typename T>
    void bcast(const T &data, int src) {
        // Start by sending the data size, so we can resize the receiving buffer accordingly
        int count = mpi::get_count(data);
        int err1  = MPI_Bcast(&count, 1, MPI_INT, src, MPI_COMM_WORLD);
        if(err1 != 0) throw std::runtime_error(fmt::format("MPI_Bcast error: {}", err1));

        if constexpr(sfinae::has_resize_v<T>) {
            data.resize(count); // Should not modify the src container
        }

        int err2 = MPI_Bcast(mpi::get_buffer(data), mpi::get_count(data), mpi::get_dtype<T>(), src, MPI_COMM_WORLD);
        if(err2 != 0) throw std::runtime_error(fmt::format("MPI_Bcast error: {}", err2));
    }

    void scatter(std::vector<h5pp::fs::path> &data, int srcId);
    void scatter_r(std::vector<h5pp::fs::path> &data, int srcId); // Roundrobin
    void broadcast(std::vector<h5pp::fs::path> &data, int srcId);
}