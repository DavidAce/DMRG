//
// Created by david on 2019-01-29.
//


#include <mps_routines/nmspc_mps_tools.h>
#include <mps_routines/class_finite_chain_state.h>
#include <mps_routines/class_superblock.h>
#include <mps_routines/class_mps_2site.h>
#include <mps_routines/class_environment.h>
#include <model/class_hamiltonian_base.h>


using Scalar         = std::complex<double>;



void MPS_Tools::Finite::Chain::copy_superblock_to_chain(class_finite_chain_state &  state, const class_superblock & superblock) {
    MPS_Tools::Finite::Chain::copy_superblock_mps_to_chain(state, superblock);
    MPS_Tools::Finite::Chain::copy_superblock_mpo_to_chain(state, superblock);
    MPS_Tools::Finite::Chain::copy_superblock_env_to_chain(state, superblock);
}

void MPS_Tools::Finite::Chain::copy_superblock_mps_to_chain(class_finite_chain_state &  state, const class_superblock & superblock) {
    auto & MPS_L  = state.get_MPS_L();
    auto & MPS_R  = state.get_MPS_R();
    auto & MPS_C  = state.get_MPS_C();

    assert(!MPS_L.empty() and !MPS_R.empty());
    assert(MPS_L.size() + MPS_R.size() <= state.get_length());
    assert(MPS_L.back().get_position()   == superblock.MPS->MPS_A->get_position());
    assert(MPS_R.front().get_position()  == superblock.MPS->MPS_B->get_position());


    MPS_L.back()    = *superblock.MPS->MPS_A;
    MPS_C           = superblock.MPS->LC;
    MPS_R.front()   = *superblock.MPS->MPS_B;
    state.all_mps_have_been_written_to_hdf5 = false;
}


void MPS_Tools::Finite::Chain::copy_superblock_mpo_to_chain(class_finite_chain_state &  state, const class_superblock & superblock){
//    std::cout << "Current state -- Overwrite: " << std::endl;
//    std::cout << "HA: " << superblock.HA->get_position() << " MPO_L back : " << MPO_L.back()->get_position() << std::endl;
//    std::cout << "HB: " << superblock.HB->get_position() << " MPO_R front: " << MPO_R.front()->get_position() << std::endl;
    auto & MPO_L  = state.get_MPO_L();
    auto & MPO_R  = state.get_MPO_R();

    assert(superblock.HA->get_position() == MPO_L.back()->get_position());
    assert(superblock.HB->get_position() == MPO_R.front()->get_position());
    MPO_L.back()    = superblock.HA->clone();
    MPO_R.front()   = superblock.HB->clone();
    assert(MPO_L.size() + MPO_R.size() == state.get_length());
    state.all_mpo_have_been_written_to_hdf5 = false;
}

void MPS_Tools::Finite::Chain::copy_superblock_env_to_chain(class_finite_chain_state &  state, const class_superblock & superblock){
    auto & ENV_L  = state.get_ENV_L();
    auto & ENV_R  = state.get_ENV_R();
    auto & ENV2_L = state.get_ENV2_L();
    auto & ENV2_R = state.get_ENV2_R();
    assert(superblock.Lblock->get_position() == ENV_L.back().get_position());
    assert(superblock.Rblock->get_position() == ENV_R.front().get_position());
    assert(superblock.Lblock2->get_position() == ENV2_L.back().get_position());
    assert(superblock.Rblock2->get_position() == ENV2_R.front().get_position());

    ENV_L.back()    = *superblock.Lblock;
    ENV2_L.back()   = *superblock.Lblock2;

    ENV_R.front()   = *superblock.Rblock;
    ENV2_R.front()  = *superblock.Rblock2;
    assert(ENV_L.size() + ENV_R.size() == state.get_length());
    state.all_env_have_been_written_to_hdf5 = false;
}


int MPS_Tools::Finite::Chain::insert_superblock_to_chain(class_finite_chain_state &  state, class_superblock & superblock){
    auto & MPS_L  = state.get_MPS_L();
    auto & MPS_R  = state.get_MPS_R();
    auto & MPS_C  = state.get_MPS_C();
    auto & MPO_L  = state.get_MPO_L();
    auto & MPO_R  = state.get_MPO_R();
    auto & ENV_L  = state.get_ENV_L();
    auto & ENV_R  = state.get_ENV_R();
    auto & ENV2_L = state.get_ENV2_L();
    auto & ENV2_R = state.get_ENV2_R();


    assert(ENV_L .size() + ENV_L.size() <= state.get_length());
    assert(MPS_L .size() + MPS_R.size() <= state.get_length());
    assert(ENV_L .size()       == superblock.Lblock->size);
    assert(ENV_R .size()       == superblock.Rblock->size);
    assert(ENV2_L.size()       == superblock.Lblock2->size);
    assert(ENV2_R.size()       == superblock.Rblock2->size);


    MPS_L.emplace_back (*superblock.MPS->MPS_A);
    MPS_R.emplace_front(*superblock.MPS->MPS_B);
    MPS_C = superblock.MPS->LC;

    ENV_L.emplace_back(*superblock.Lblock);
    ENV_R.emplace_front(*superblock.Rblock);
    ENV2_L.emplace_back(*superblock.Lblock2);
    ENV2_R.emplace_front(*superblock.Rblock2);
    MPO_L.emplace_back      (superblock.HA->clone());
    MPO_R.emplace_front     (superblock.HB->clone());


    int pos = 0;
    for (auto &MPS: MPS_L){MPS.set_position(pos++);}
    for (auto &MPS: MPS_R){MPS.set_position(pos++);}
    pos = 0;
    for (auto &ENV: ENV_L){ENV.set_position(pos++);}
    for (auto &ENV: ENV_R){ENV.set_position(pos++);}
    pos = 0;
    for (auto &ENV2: ENV2_L){ENV2.set_position(pos++);}
    for (auto &ENV2: ENV2_R){ENV2.set_position(pos++);}
    pos = 0;
    for (auto &MPO : MPO_L){MPO->set_position(pos++);}
    for (auto &MPO : MPO_R){MPO->set_position(pos++);}

    superblock.MPS->MPS_A->set_position(MPS_L.back().get_position());
    superblock.MPS->MPS_B->set_position(MPS_R.front().get_position());
    superblock.Lblock->set_position(ENV_L.back().get_position());
    superblock.Rblock->set_position(ENV_R.front().get_position());
    superblock.Lblock2->set_position(ENV2_L.back().get_position());
    superblock.Rblock2->set_position(ENV2_R.front().get_position());
    superblock.HA->set_position(MPO_L.back()->get_position());
    superblock.HB->set_position(MPO_R.front()->get_position());

    assert(ENV_L.back().size + ENV_R.front().size == superblock.environment_size);
    assert(ENV_L.back().size   == superblock.Lblock->size);
    assert(ENV_R.front().size  == superblock.Rblock->size);
    assert(ENV2_L.back().size  == superblock.Lblock2->size);
    assert(ENV2_R.front().size == superblock.Rblock2->size);

    state.all_mps_have_been_written_to_hdf5 = false;
    state.all_mpo_have_been_written_to_hdf5 = false;
    state.all_env_have_been_written_to_hdf5 = false;
    return (int)MPS_L.size();
}




int MPS_Tools::Finite::Chain::move_center_point(class_finite_chain_state &  state, class_superblock & superblock){
    //Take current MPS and generate an Lblock one larger and store it in list for later loading
//    std::cout << "Current state -- Direction: " << direction << std::endl;
//    std::cout << "HA: " << superblock.HA->get_position() << " MPO_L back : " << MPO_L.back()->get_position() << std::endl;
//    std::cout << "HB: " << superblock.HB->get_position() << " MPO_R front: " << MPO_R.front()->get_position() << std::endl;
//
    auto & MPS_L  = state.get_MPS_L();
    auto & MPS_R  = state.get_MPS_R();
    auto & MPS_C  = state.get_MPS_C();
    auto & MPO_L  = state.get_MPO_L();
    auto & MPO_R  = state.get_MPO_R();
    auto & ENV_L  = state.get_ENV_L();
    auto & ENV_R  = state.get_ENV_R();
    auto & ENV2_L = state.get_ENV2_L();
    auto & ENV2_R = state.get_ENV2_R();
    assert(!MPS_L.empty() and !MPS_R.empty());
    assert(MPS_L.size() + MPS_R.size() == state.get_length());
    assert(ENV_L.size() + ENV_R.size() == state.get_length());
    assert(ENV_L.back().size + ENV_R.front().size == state.get_length() - 2);
    assert(ENV_L.back().size + ENV_R.front().size == superblock.environment_size);

    assert(MPS_L.back().get_position()   == superblock.MPS->MPS_A->get_position());
    assert(MPS_R.front().get_position()  == superblock.MPS->MPS_B->get_position());


    assert(MPO_L.back()->get_position()  == superblock.HA->get_position() );
    assert(MPO_R.front()->get_position() == superblock.HB->get_position() );

    if (state.get_direction() == 1){
        //Note that Lblock must just have grown!!
//        assert(MPS_R.size() > 1);
        assert(ENV_L.back().size   + 1  == superblock.Lblock->size);
        assert(ENV2_L.back().size  + 1  == superblock.Lblock2->size);
        assert(ENV_R.front().size       == superblock.Rblock->size);
        assert(ENV2_R.front().size      == superblock.Rblock2->size);
        assert(ENV_L.back().get_position()  + 1 == superblock.Lblock->get_position());
        assert(ENV2_L.back().get_position() + 1 == superblock.Lblock2->get_position());
        assert(ENV_R.front().get_position()     == superblock.Rblock->get_position());
        assert(ENV2_R.front().get_position()    == superblock.Rblock2->get_position());


        MPS_L.emplace_back(*superblock.MPS->MPS_B);
        MPO_L.emplace_back (superblock.HB->clone());
        ENV_L.emplace_back (*superblock.Lblock);
        ENV2_L.emplace_back (*superblock.Lblock2);


        MPS_R.pop_front();
        MPO_R.pop_front();
        ENV_R.pop_front();
        ENV2_R.pop_front();


        superblock.MPS->MPS_A->set_L(superblock.MPS->LC);
        superblock.MPS->MPS_A->set_G(superblock.MPS->MPS_B->get_G());
        superblock.MPS->MPS_A->set_position(MPS_L.back().get_position());
        superblock.MPS->LC = superblock.MPS->MPS_B->get_L();
        superblock.MPS->MPS_B->set_G(MPS_R.front().get_G());
        superblock.MPS->MPS_B->set_L(MPS_R.front().get_L());
        superblock.MPS->MPS_B->set_position(MPS_R.front().get_position());
        MPS_C  = superblock.MPS->LC;

        superblock.HA = MPO_L.back()->clone();
        superblock.HB = MPO_R.front()->clone();

        *superblock.Rblock  = ENV_R.front();
        *superblock.Rblock2 = ENV2_R.front();

    }else{
        //Note that Rblock must just have grown!!
//        assert(MPS_L.size() > 1);
        assert(ENV_R.front().size  + 1  == superblock.Rblock->size);
        assert(ENV2_R.front().size + 1  == superblock.Rblock2->size);
        assert(ENV_L.back().size        == superblock.Lblock->size);
        assert(ENV2_L.back().size       == superblock.Lblock2->size);
        assert(ENV_L.back().get_position()      == superblock.Lblock->get_position());
        assert(ENV2_L.back().get_position()     == superblock.Lblock2->get_position());
        assert(ENV_R.front().get_position()  -1 == superblock.Rblock->get_position());
        assert(ENV2_R.front().get_position() -1 == superblock.Rblock2->get_position());
        MPS_R.emplace_front (*superblock.MPS->MPS_A);
        MPO_R.emplace_front (superblock.HA->clone());
        ENV_R.emplace_front (*superblock.Rblock);
        ENV2_R.emplace_front(*superblock.Rblock2);

        MPS_L.pop_back();
        MPO_L.pop_back();
        ENV_L.pop_back();
        ENV2_L.pop_back();

        superblock.MPS->MPS_B->set_L(superblock.MPS->LC);
        superblock.MPS->MPS_B->set_G(superblock.MPS->MPS_A->get_G());
        superblock.MPS->MPS_B->set_position(MPS_R.front().get_position());
        superblock.MPS->LC = superblock.MPS->MPS_A->get_L();
        superblock.MPS->MPS_A->set_G(MPS_L.back().get_G());
        superblock.MPS->MPS_A->set_L(MPS_L.back().get_L());
        superblock.MPS->MPS_A->set_position(MPS_L.back().get_position());

        MPS_C                       = superblock.MPS->LC;

        superblock.HA = MPO_L.back()->clone();
        superblock.HB = MPO_R.front()->clone();

        *superblock.Lblock  = ENV_L.back();
        *superblock.Lblock2 = ENV2_L.back();
    }

    assert(superblock.MPS->MPS_A->get_G().dimension(2) == superblock.MPS->MPS_B->get_G().dimension(1));
    assert(MPO_L.size() + MPO_R.size() == state.get_length());
    assert(superblock.HA->get_position() + 1 == MPO_L.size());
    assert(superblock.HA->get_position() + 1 == state.get_length() - MPO_R.size());
    assert(superblock.HB->get_position() + 1 == MPO_L.size() + 1);
    assert(superblock.HB->get_position() + 1 == state.get_length() - MPO_R.size() + 1);
    assert(MPS_L.back().get_position()   == superblock.MPS->MPS_A->get_position());
    assert(MPS_R.front().get_position()  == superblock.MPS->MPS_B->get_position());

    assert(ENV_L.back().get_position()   == superblock.Lblock->get_position());
    assert(ENV_R.front().get_position()  == superblock.Rblock->get_position());
    assert(ENV2_L.back().get_position()  == superblock.Lblock2->get_position());
    assert(ENV2_R.front().get_position() == superblock.Rblock2->get_position());

    assert(MPO_L.back()->get_position()  == superblock.HA->get_position() );
    assert(MPO_R.front()->get_position() == superblock.HB->get_position() );

    //    Check edge
    if (state.position_is_the_left_edge() or state.position_is_the_right_edge()) {
        state.flip_direction();
    }
    if (state.position_is_the_left_edge()){
        state.increment_sweeps();
    }
    state.all_mps_have_been_written_to_hdf5 = false;
    state.all_mpo_have_been_written_to_hdf5 = false;
    state.all_env_have_been_written_to_hdf5 = false;
    return state.get_sweeps();
}




