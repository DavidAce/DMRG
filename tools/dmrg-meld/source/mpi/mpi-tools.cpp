
#include "mpi-tools.h"
#include "general/iter.h"
#include "math/num.h"
#include "meld-io/logger.h"
#include <mpi/mpi.h>
void mpi::init() {
    // Initialize the MPI environment
    MPI_Init(nullptr, nullptr);
    MPI_Comm_rank(MPI_COMM_WORLD, &world.id);   // Establish thread number of this worker
    MPI_Comm_size(MPI_COMM_WORLD, &world.size); // Get total number of threads

    if(world.id == 0) tools::logger::log->info("MPI initialized with {} processes", world.size);
    if(world.size > 1) {
        mpi::on = true;

        auto logname       = tools::logger::log->name();
        auto width         = fmt::format("{}", world.size - 1).size();
        auto mpiname       = fmt::format("{}-{:>{}}", tools::logger::log->name(), world.id, width);
        tools::logger::log = tools::logger::setLogger(mpiname, static_cast<size_t>(tools::logger::log->level()));
    }
}

void mpi::finalize() {
    if(world.size > 1 or mpi::on) MPI_Finalize();
}

void mpi::barrier() {
    if(world.size > 1 or mpi::on) MPI_Barrier(MPI_COMM_WORLD);
}

void mpi::scatter(std::vector<h5pp::fs::path> &data, int src) {
    if(world.size == 1) return; // No need to scatter

    // Since the paths are not plain data, we need to mimic the logic of scatter
    // Every path is constructible from a string, so we can simply send the string data around.

    // First we need to decide how to split the data into even chunks
    std::vector<size_t> counts(world.get_size<size_t>(), 0);
    size_t              sum_counts = 0;
    if(world.id == src) {
        for(const auto &[i, d] : iter::enumerate(data)) counts[num::mod<size_t>(i, world.get_size<size_t>())] += 1;
        sum_counts = num::sum(counts);
        if(sum_counts != data.size()) throw std::runtime_error(h5pp::format("mpi::scatter: counts does not add up: {} = {}", counts, num::sum(counts)));
    }
    // Communicate to everyone how many elements they should expect
    MPI_Bcast(mpi::get_buffer(counts), mpi::get_count(counts), mpi::get_dtype<size_t>(), src, MPI_COMM_WORLD);
    sum_counts = num::sum(counts);

    // Reserve space to receive data
    std::vector<h5pp::fs::path> srcData; // Elements to keep in src
    if(world.id != src)
        data.reserve(counts[world.get_id<size_t>()]);
    else
        srcData.reserve(counts[world.get_id<size_t>()]);

    // Now we can start sending data

    size_t srcidx = 0;
    for(int dst = 0; dst < world.size; dst++) {
        if(world.id == src or world.id == dst) {
            size_t c = counts[static_cast<unsigned long>(dst)];
            for(size_t dstidx = 0; dstidx < c; dstidx++) {
                if(world.id == dst and src != dst) {
                    std::string tmp;
                    mpi::recv(tmp, src, 0);
                    data.emplace_back(tmp);
                } else if(world.id == src) {
                    if(dst == src) {
                        srcData.emplace_back(data[srcidx]);
                    } else
                        mpi::send(data[srcidx].string(), dst, 0);
                    srcidx++;
                }
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }
    if(world.id == src) data = srcData; // Now src can let go of the sent data and keep only what corresponds to src

    MPI_Barrier(MPI_COMM_WORLD);
    if(data.size() != counts[world.get_id<size_t>()])
        throw std::logic_error(h5pp::format("mpi::scatter on id {} failed: counts {} | size {}", world.id, counts, data.size()));
}

void mpi::scatter_r(std::vector<h5pp::fs::path> &data, int src) {
    if(world.size == 1) return; // No need to scatter

    // Scatter in roundrobin mode

    // Since the paths are not plain data, we need to mimic the logic of scatter
    // Every path is constructible from a string, so we can simply send the string data around.

    // First we need to decide how to split the data into even chunks
    std::vector<size_t> counts(world.get_size<size_t>(), 0);
    size_t              sum_counts = 0;
    if(world.id == src) {
        for(const auto &[i, d] : iter::enumerate(data)) counts[num::mod<size_t>(i, world.get_size<size_t>())] += 1;
        sum_counts = num::sum(counts);
        if(sum_counts != data.size()) throw std::runtime_error(h5pp::format("mpi::scatter: counts does not add up: {} = {}", counts, num::sum(counts)));
    }

    // Communicate to everyone how many elements they should expect
    MPI_Bcast(mpi::get_buffer(counts), mpi::get_count(counts), mpi::get_dtype<size_t>(), src, MPI_COMM_WORLD);
    sum_counts = num::sum(counts);

    // Reserve space to receive data
    std::vector<h5pp::fs::path> srcData; // Elements to keep in src
    if(world.id != src)
        data.reserve(counts[world.get_id<size_t>()]);
    else
        srcData.reserve(counts[world.get_id<size_t>()]);

    // Now we can start sending data

    for(size_t srcidx = 0; srcidx < sum_counts; srcidx++) {
        int dst = num::mod<int>(static_cast<int>(srcidx), world.size); // Roundrobin destination
        if(world.id == src or world.id == dst) {
            if(world.id == dst and src != dst) {
                std::string tmp;
                mpi::recv(tmp, src, 0);
                data.emplace_back(tmp);
            } else if(world.id == src) {
                if(dst == src) {
                    srcData.emplace_back(data[srcidx]);
                } else
                    mpi::send(data[srcidx].string(), dst, 0);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    if(world.id == src) data = srcData; // Now src can let go of the sent data and keep only what corresponds to src

    MPI_Barrier(MPI_COMM_WORLD);
    if(data.size() != counts[world.get_id<size_t>()])
        throw std::logic_error(h5pp::format("mpi::scatter on id {} failed: counts {} | size {}", world.id, counts, data.size()));
}

void mpi::broadcast(std::vector<h5pp::fs::path> &data, int src) {
    if(world.size == 1) return; // No need to scatter

    // Broadcast so all ranks get a copy of data

    // Since the paths are not plain data, we need to mimic the logic of scatter
    // Every path is constructible from a string, so we can simply send the string data around.

    // First we need to communitate how many elements everyone should allocate before receiving them
    std::vector<size_t> counts(world.get_size<size_t>(), data.size()); // The correct data size lives in src

    // Communicate to everyone how many elements they should expect: Everyone gets the correct counts from src
    MPI_Bcast(mpi::get_buffer(counts), mpi::get_count(counts), mpi::get_dtype<size_t>(), src, MPI_COMM_WORLD);

    // Reserve space to receive data
    if(world.id != src) data.resize(counts.at(world.get_id<size_t>()));

    // Now we can start broadcasting the actual data
    for (const auto & elem : data) {
        mpi::bcast
    }

    for(size_t srcidx = 0; srcidx < sum_counts; srcidx++) {
        int dst = num::mod<int>(static_cast<int>(srcidx), world.size); // Roundrobin destination
        if(world.id == src or world.id == dst) {
            if(world.id == dst and src != dst) {
                std::string tmp;
                mpi::recv(tmp, src, 0);
                data.emplace_back(tmp);
            } else if(world.id == src) {
                if(dst == src) {
                    srcData.emplace_back(data[srcidx]);
                } else
                    mpi::send(data[srcidx].string(), dst, 0);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    if(world.id == src) data = srcData; // Now src can let go of the sent data and keep only what corresponds to src

    MPI_Barrier(MPI_COMM_WORLD);
    if(data.size() != counts[world.get_id<size_t>()])
        throw std::logic_error(h5pp::format("mpi::scatter on id {} failed: counts {} | size {}", world.id, counts, data.size()));
}

