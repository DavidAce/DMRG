//
// Created by david on 2018-12-16.
//


#include <IO/class_hdf5_file.h>
#include <general/nmspc_tensor_extra.h>
#include <iomanip>
#include <iostream>


/*! \brief Prints the content of a vector nicely */
template<typename T>
std::ostream &operator<<(std::ostream &out, const std::vector<T> &v) {
    if (!v.empty()) {
        out << "[ ";
        std::copy(v.begin(), v.end(), std::ostream_iterator<T>(out, " "));
        out << "]";
    }
    return out;
}

using namespace std::complex_literals;



// Store some dummy data to an hdf5 file



int main()
{
    using Scalar = std::complex<double>;

    static_assert(TypeCheck::has_member_data<std::vector<double>>() and "Compile time type-checker failed. Could not properly detect class member data. Check that you are using a supported compiler!");
    static_assert(TypeCheck::has_member_data<std::vector<double>>() and "Compile time type-checker failed. Could not properly detect class member data. Check that you are using a supported compiler!");
    static_assert(TypeCheck::has_member_data<std::vector<double>>() and "Compile time type-checker failed. Could not properly detect class member data. Check that you are using a supported compiler!");

    // Generate dummy data
    std::vector<double> vectorD = { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                               0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                               -1.0, 0.0, 1.0, 0.0,-1.0, 0.0, 0.0,-1.0, 0.0, 1.0, 0.0, 1.0};
    std::vector<Scalar> vectorC = {-0.191154, 0.964728, -0.0351791,0.177544};

    Eigen::Tensor<Scalar,4> tensorC(3,3,3,3);
    tensorC.setRandom();

    std::string stringX = "Dummy string with spaces";


    std::string output_filename         = "hdf5_test.h5";
    std::string output_folder           = "Testing/hdf5";
    bool        create_dir_if_not_found = true;
    bool        overwrite_file_if_found = true;
    bool        resume_from_file        = false;


    class_hdf5_file hdf5(output_filename,
                         output_folder,
                         create_dir_if_not_found,
                         overwrite_file_if_found,
                         resume_from_file );

    std::cout << "Writing vectorD: " <<  vectorD  << std::endl;
    std::cout << "Writing vectorC: " <<  vectorC  << std::endl;
    std::cout << "Writing tensorC: " <<  tensorC  << std::endl;
    std::cout << "Writing stringX: " <<  stringX  << std::endl;

    hdf5.write_dataset(vectorD,"vectorD");
    hdf5.write_dataset(vectorC,"vectorC");
    hdf5.write_dataset(tensorC,"tensorC");
    hdf5.write_dataset(stringX,"stringX");
    hdf5.write_attribute_to_dataset("vectorD" ,std::string("This is an attribute"), "TestAttr");
    hdf5.write_attribute_to_dataset("vectorC" ,std::string("This is an attribute"), "TestAttr");
    hdf5.write_attribute_to_dataset("tensorC" ,std::string("This is an attribute"), "TestAttr");
    hdf5.write_attribute_to_dataset("stringX" ,std::string("This is an attribute"), "TestAttr");


    // Read the data back
    std::vector<double>     vectorD_read;
    std::vector<Scalar>     vectorC_read;
    Eigen::Tensor<Scalar,4> tensorC_read;
    std::string             stringX_read;


    hdf5.read_dataset(vectorD_read,"vectorD");
    hdf5.read_dataset(vectorC_read,"vectorC");
    hdf5.read_dataset(tensorC_read,"tensorC");
    hdf5.read_dataset(stringX_read,"stringX");
    std::cout << "Reading vectorD: " <<  vectorD_read  << std::endl;
    std::cout << "Reading vectorC: " <<  vectorC_read  << std::endl;
    std::cout << "Reading tensorC: " <<  tensorC_read  << std::endl;
    std::cout << "Reading stringX: " <<  stringX_read  << std::endl;

    if (vectorD != vectorD_read){return 1;}
    if (vectorC != vectorC_read){return 1;}
    if (stringX != stringX_read){return 1;}

    //Tensor comparison isn't as straightforward
    Eigen::Map<Eigen::VectorXcd> tensorMap(tensorC.data(), tensorC.size());
    Eigen::Map<Eigen::VectorXcd> tensorMapRead(tensorC_read.data(), tensorC_read.size());
    if (tensorMap != tensorMapRead){return 1;}

    return 0;
}