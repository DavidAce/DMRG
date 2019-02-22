//
// Created by david on 2017-12-02.
//

#ifndef DMRG_CLASS_HDF5_DATASET_BUFFER_H
#define DMRG_CLASS_HDF5_DATASET_BUFFER_H

class class_hdf5_file;
template<typename DataType, typename AttrType = int, typename IterType = int>
class class_hdf5_dataset_buffer : public std::vector<DataType>{
private:
    std::shared_ptr<class_hdf5_file> hdf5_out;
    std::string group_name      = "default_group";
    IterType iteration          = 0;
    std::string dataset_name    = "default_data";
    AttrType attribute          = 0;
    std::string attribute_name  = "default_attribute";
    bool attribute_set          = false;
    int pad_to_width            = 3;
    int max_elements            = 1000;
    std::string dataset_relative_name;
    std::string left_pad(const std::string &temp);
    bool data_has_been_written_to_file = false;

public:

    explicit class_hdf5_dataset_buffer()=default;

    class_hdf5_dataset_buffer(std::shared_ptr<class_hdf5_file> hdf5_out_,
                              const std::string &group_name_,
                              const IterType &iteration_,
                              const std::string &dataset_name_);

    class_hdf5_dataset_buffer(const std::string &group_name_,
                              const IterType &iteration_,
                              const std::string &dataset_name_);

    class_hdf5_dataset_buffer(std::shared_ptr<class_hdf5_file> hdf5_out_,
                              const std::string &group_name_,
                              const IterType &iteration_,
                              const std::string &dataset_name_,
                              const AttrType &attribute_,
                              const std::string &attribute_name_);


    class_hdf5_dataset_buffer(const std::string &group_name_,
                              const IterType &iteration_,
                              const std::string &dataset_name_,
                              const AttrType &attribute_,
                              const std::string &attribute_name_);

    ~class_hdf5_dataset_buffer(){
        if (hdf5_out != nullptr){
            write_buffer_to_file();
        }else if (!data_has_been_written_to_file){
            std::cerr << "Warning: Output data has not been saved to file, yet it is being discarded!\n" << std::endl;
        }
    }

    void hdf5_set_group_name    (const std::string &group_name_);
    void hdf5_set_dataset_name  (const std::string &dataset_name_);
    void hdf5_set_iteration     (const IterType &iteration_);
    void hdf5_set_width     (const int &pad_to_width_);
    void hdf5_set_all(const std::string &group_name_, const IterType iteration_, const std::string &dataset_name_);
    void hdf5_set_attribute(const AttrType &attribute_, const std::string &attribute_name_);
    void hdf5_refresh_relative_name();
    void write_buffer_to_file();


    void push_back(const DataType& element) {
        if ((int)std::vector<DataType>::size() >= max_elements) {
            write_buffer_to_file();
        }
        std::vector<DataType>::push_back(element);
        data_has_been_written_to_file = false;
    }
};


#endif //DMRG_CLASS_HDF5_DATASET_BUFFER_H
