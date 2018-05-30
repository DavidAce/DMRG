
#ifndef DMRG_CLASS_STORAGE_H
#define DMRG_CLASS_STORAGE_H
#include <map>
#include <list>
#include <memory>
#include <variant>
#include <complex>
#include <general/nmspc_tensor_extra.h>
#include <mps_routines/class_environment.h>
#include <mps_routines/class_hamiltonian.h>
#include <sim_parameters/nmspc_model.h>
#include <sim_parameters/nmspc_sim_settings.h>
#include <IO/class_hdf5_file.h>
#include <iostream>

class class_superblock;


/*!
 \class class_environment_storage
 \brief Storage class for environments during finite DMRG
 During the infinite-DMRG part of the algorithm, the current state for each chain-length
 is stored in a `class_finite_chain_storage` object, for later use in the finite-DMRG stage.

 There are as many elements in storage as there are particles in the current environment, i.e. superblock->environment_size.



 The storage is always partitioned into two lists, corresponding to the two halves about the current position.
 The back and front of these lists correspond to the current state of the superblock's left and right parts, respectively,
 while "LA" is the central \f$ \Lambda^A \f$.

 Once full, the particle content of the chain should match the length given in the constant max_length.
 This can be checked with the following conditions:
 ENV_L.back().size + ENV_R.front().size == max_length - 2
 ENV_L.back().size + ENV_R.front().size == superblock->environment_size
 ENV_L.back().size == superblock->Lblock->size
 ENV_R.back().size == superblock->Rblock->size


*/



class class_finite_chain_sweeper {
public:
    using Scalar = std::complex<double>;


private:
    std::list<std::tuple<Textra::Tensor<Scalar,1>,Textra::Tensor<Scalar,3>>>  MPS_L;  /*!< A list of stored \f$ \Lambda^B \Gamma^A...  \f$-tensors. */
    std::list<std::tuple<Textra::Tensor<Scalar,3>,Textra::Tensor<Scalar,1>>>  MPS_R;  /*!< A list of stored \f$ \Gamma^B \Lambda^B...  \f$-tensors. */
    Textra::Tensor<Scalar,1> MPS_C;  //Current center bond matrix;
    std::list<class_environment> ENV_L;
    std::list<class_environment> ENV_R;
    std::list<class_environment_var> ENV2_L;
    std::list<class_environment_var> ENV2_R;
    std::list<class_hamiltonian> MPO_L;                                           /*!< A list of stored Hamiltonian MPO tensors,indexed by chain position. */
    std::list<class_hamiltonian> MPO_R;                                           /*!< A list of stored Hamiltonian MPO tensors,indexed by chain position. */



    std::shared_ptr<class_superblock> superblock;
    std::shared_ptr<class_hdf5_file>  hdf5;
    SimulationType sim_type;
    std::string    sim_name;

    bool max_length_is_set = false;
    bool superblock_is_set = false;
    bool hdf5_file_is_set  = false;

    int direction = 1;
    int sweeps    = 0;
    unsigned long max_length = 0;                                                 /*!< The maximum length of the chain */
    unsigned long current_length = 0;

    template<int tup=-1, typename T>
    void write_list_to_file(const std::list<T> &obj, std::string dataset_name, unsigned long &counter){
        for(auto &it: obj) {
            if constexpr(tup == -1) {
                if constexpr(std::is_same<std::decay_t<decltype(it)>, class_hamiltonian>::value){
                    hdf5->write_dataset(it.MPO, dataset_name + "_" + std::to_string(counter));
                    hdf5->write_attribute_to_dataset(dataset_name + "_" + std::to_string(counter), it.get_site_coupling(), "coupling");
                    hdf5->write_attribute_to_dataset(dataset_name + "_" + std::to_string(counter), it.get_site_field(), "field");
                    hdf5->write_attribute_to_dataset(dataset_name + "_" + std::to_string(counter), it.get_site_random_field(), "random_field");
                    hdf5->write_attribute_to_dataset(dataset_name + "_" + std::to_string(counter), it.get_site_energy(), "energy");
                    counter++;
                }
                else if constexpr(std::is_same<std::decay_t<decltype(it)>, class_environment>::value or
                                  std::is_same<std::decay_t<decltype(it)>, class_environment_var>::value  ){
                    hdf5->write_dataset(it.block, dataset_name + "_" + std::to_string(counter));
                    hdf5->write_attribute_to_dataset(dataset_name + "_" + std::to_string(counter), it.size, "sites");
                    counter++;
                }else{
                    hdf5->write_dataset(it, dataset_name + "_" + std::to_string(counter++));
                }
            } else {
                hdf5->write_dataset(std::get<tup>(it), dataset_name + "_" + std::to_string(counter++));
            }
        }
    }

public:

    void print_error_and_exit(int error_type);
    class_finite_chain_sweeper()=default;
    explicit class_finite_chain_sweeper(int max_length_,                         /*!< The maximum length of the chain. */
                                 std::shared_ptr<class_superblock> superblock_,
                                 std::shared_ptr<class_hdf5_file>  hdf5_,
                                 SimulationType sim_type,
                                 std::string    sim_name
    );

    void set_max_length(int max_length_);                                        /*!< Sets the maximum length of the chain. */
    void set_superblock(std::shared_ptr<class_superblock> superblock_);          /*!< Sets the pointer to a superblock */
    void set_hdf5_file (std::shared_ptr<class_hdf5_file>  hdf5_);                /*!< Sets the pointer to an hdf5-file for storage */
    void update_current_length();
    int  insert();                                                               /*!< Store current MPS and environments indexed by their respective positions on the chain. */
    int  load();                                                                 /*!< Load MPS and environments according to current position. */
    void overwrite_MPS();                                                        /*!< Update the MPS stored at current position.*/
    int  move();                                                                 /*!< Move current position to the left (`direction=1`) or right (`direction=-1`), and store the **newly enlarged** environment. Turn direction around if the edge is reached. */
    void write_chain_to_file();

    const auto & get_MPS_L() const {return std::as_const(MPS_L);}
    const auto & get_MPS_R() const {return std::as_const(MPS_R);}
    const auto & get_MPS_C() const {return std::as_const(MPS_C);}
    const auto & get_MPO_L() const {return std::as_const(MPO_L);}
    const auto & get_MPO_R() const {return std::as_const(MPO_R);}
    const auto & get_ENV_L() const {return std::as_const(ENV_L);}
    const auto & get_ENV_R() const {return std::as_const(ENV_R);}
    const auto & get_ENV2_L()const {return std::as_const(ENV2_L);}
    const auto & get_ENV2_R()const {return std::as_const(ENV2_R);}

    void print_storage();                                                        /*!< Print the tensor dimensions for all \f$\Gamma\f$-tensors. */
    void print_storage_compact();                                                /*!< Print the tensor dimensions for all \f$\Gamma\f$-tensors. */
    void print_hamiltonians();
    int get_direction() const;
    int get_position() const;
    int get_length() const;
    int get_sweeps() const;
    bool position_is_the_middle();
    bool position_is_the_left_edge();
    bool position_is_the_right_edge();


};




#endif //DMRG_CLASS_STORAGE_H
