#pragma once

namespace Textra{
    extern template idxlistpair<Eigen::Index,1> idx (const Eigen::Index (&list1)[1], const Eigen::Index (&list2)[1]);
    extern template idxlistpair<Eigen::Index,2> idx (const Eigen::Index (&list1)[2], const Eigen::Index (&list2)[2]);
    extern template idxlistpair<Eigen::Index,3> idx (const Eigen::Index (&list1)[3], const Eigen::Index (&list2)[3]);
    extern template idxlistpair<Eigen::Index,4> idx (const Eigen::Index (&list1)[4], const Eigen::Index (&list2)[4]);
    extern template idxlistpair<Eigen::Index,5> idx (const Eigen::Index (&list1)[5], const Eigen::Index (&list2)[5]);
}