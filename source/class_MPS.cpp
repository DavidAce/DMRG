#include <class_MPS.h>
using namespace std;


void class_MPS::initialize(const long local_dimension_){
    local_dimension = local_dimension_;
    GA.resize(array3{local_dimension,1,1});
    GA.setZero();
    GA(0,0,0) = 1;
    LA.resize(array1{1});
    LA.setConstant(1.0);

    GB.resize(array3{local_dimension,1,1});
    GB.setZero();
    GB(0,0,0) = 1;
    LB.resize(array1{1});
    LB.setConstant(1.0);
    swapped = false;
    L_tail = LB; /*! \todo , check?  was LA before */

    std::cout << "Initial transfer matrix LL: " << GA.contract(GA, idxlist2{idx2(0,0), idx2(1,1)}) << endl;
    std::cout << "Initial transfer matrix LR: " << GA.contract(GA, idxlist2{idx2(0,0), idx2(2,2)}) << endl;
    std::cout << "Initial transfer matrix RL: " << GB.contract(GB, idxlist2{idx2(0,0), idx2(1,1)}) << endl;
    std::cout << "Initial transfer matrix RR: " << GB.contract(GB, idxlist2{idx2(0,0), idx2(2,2)}) << endl;
}


Tensor4d class_MPS::get_theta() const {
    return asDiagonal(L_tail) //whatever L_A was in the previous step
            .contract(GA, idxlist1{idx2(1,1)})
            .contract(asDiagonal(LA), idxlist1{idx2(2,0)})
            .contract(GB, idxlist1{idx2(2,1)})
            .contract(asDiagonal(LB), idxlist1{idx2(3,0)});
    //Outputs:
    //      1  2
    //   0__|__|__3
}

void class_MPS::swap_AB(){
    swapped    = !swapped;

    //Swap Gamma
    tmp3    = GA;
    GA      = GB;
    GB      = tmp3;

    //Swap Lambda
    tmp1    = LA;
    LA      = LB;
    LB      = tmp1;
    L_tail  = LB;

}

double class_MPS::get_energy(const Tensor4d &Hamiltonian) {
    Tensor4d theta = get_theta();//.shuffle(array4{1,0,2,3});
    Tensor0d E1 =  Hamiltonian.contract(theta.conjugate(), idxlist2{idx2(0,1), idx2(1,2)})
            .contract(theta, idxlist4{idx2(0,1),idx2(1,2),idx2(2,0),idx2(3,3)});

//    swap_AB();
//    theta = get_theta(); //.shuffle(array4{1,0,2,3});
//    Tensor0d E2  =  Hamiltonian.contract(theta.conjugate(), idxlist2{idx2(0,1), idx2(1,2)})
//            .contract(theta, idxlist4{idx2(0,1),idx2(1,2),idx2(2,0),idx2(3,3)});
//    swap_AB();
//    cout << GA <<endl<<endl;
//    cout << GB <<endl<<endl;
//    cout << LA <<endl<<endl;
//    cout << LB <<endl<<endl;


//    cout << "energies: "<< E1 << " " <<  E2 << endl;
    return (E1(0));
//    return (E1(0) + E2(0))/2.0;
//    return  Hamiltonian.contract(LGLGL.conjugate(), idxlist2{idx2(0,0), idx2(1,1)})
//            .contract(LGLGL, idxlist4{idx2(0,0),idx2(1,1),idx2(2,2),idx2(3,3)});

}

double class_MPS::get_variance(const Tensor4d &Hamiltonian) {
    Tensor4d theta = get_theta();//.shuffle(array4{1,0,2,3});
    Tensor4d H_sq = Hamiltonian.square();
    Tensor0d Var =  H_sq.contract(theta.conjugate(), idxlist2{idx2(0,1), idx2(1,2)})
            .contract(theta, idxlist4{idx2(0,1),idx2(1,2),idx2(2,0),idx2(3,3)});

    return (Var(0) - std::pow(get_energy(Hamiltonian), 2));
//    return (Var(0));
}


double class_MPS::get_entropy() {
    Tensor0d result1  = -LB.square().contract(LB.square().log().eval(), idxlist1{idx2(0,0)});
    return result1(0) ;
}