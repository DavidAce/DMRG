#pragma once
#include <iterator>
namespace iter {
    template<class T>
    class reverse {
        private:
        T &rng;

        public:
        reverse(T &r) noexcept : rng(r) {}

        auto begin() const noexcept {
            using std::end;
            return std::make_reverse_iterator(end(rng));
        }
        auto end() const noexcept {
            using std::begin;
            return std::make_reverse_iterator(begin(rng));
        }
    };

    namespace internal {

        template<class Iterator, bool reverse = false>
        struct enumerate_iterator {
            using iterator   = Iterator;
            using index_type = std::size_t;
            using reference  = typename std::iterator_traits<iterator>::reference;
            using pointer    = typename std::iterator_traits<iterator>::pointer;

            enumerate_iterator(index_type index, iterator iterator) : index(index), iter(iterator) {}
            enumerate_iterator &operator++() {
                if constexpr(reverse) --index;
                else
                    ++index;
                ++iter;
                return *this;
            }

            bool operator!=(const enumerate_iterator &other) const { return iter != other.iter; }

            std::pair<std::reference_wrapper<index_type>, reference> operator*() { return {index, *iter}; }

            private:
            index_type index;
            iterator   iter;
        };

        template<class Iterator,bool reverse=false>
        struct enumerate_range {
            using index_type = std::size_t;
            using iterator   = enumerate_iterator<Iterator,reverse>;

            enumerate_range(Iterator first, Iterator last, index_type initial, index_type final)
                : first(first), last(last), initial(initial),final(final) {}

            iterator begin() const { return iterator(initial, first); }
            iterator end() const { return iterator(final, last); }

            private:
            Iterator   first;
            Iterator   last;
            index_type initial;
            index_type final;
        };
    }

    template<class Iterator>
    decltype(auto) enumerate(Iterator first, Iterator last, std::size_t initial) {
        return internal::enumerate_range<Iterator,false>(first, last, initial);
    }

    template<class Container>
    decltype(auto) enumerate(Container &content) {
        return internal::enumerate_range<typename Container::iterator,false>(std::begin(content), std::end(content), 0, content.size()-1);
    }

    template<class Container>
    decltype(auto) enumerate(const Container &content) {
        return internal::enumerate_range<typename Container::const_iterator,false>(std::begin(content), std::end(content), 0, content.size()-1);
    }

    template<class Iterator>
    decltype(auto) enumerate_reverse(Iterator first, Iterator last, std::size_t initial,std::size_t final) {
        return internal::enumerate_range<Iterator,true>(first, last, initial, final);
    }

    template<class Container>
    decltype(auto) enumerate_reverse(Container &content) {
        return internal::enumerate_range<typename Container::reverse_iterator,true>(std::rbegin(content), std::rend(content), content.size()-1, 0);
    }

//    template<typename T>
//    struct reverse_enumerate {
//        private:
//        T &container;
//        struct iterator {
//            typename T::reverse_iterator iter;
//            size_t                       index;
//
//            bool operator==(const iterator &other) const { return iter == other.iter; }
//            bool operator!=(const iterator &other) const { return iter != other.iter; }
//            auto operator*() { return std::pair(index, *iter); }
//            void operator++() {
//                --index;
//                ++iter;
//            }
//        };
//
//        public:
//        reverse_enumerate(T &container) : container(container) {}
//        iterator begin() { return {container.size() - 1, container.rbegin()}; }
//        iterator end() { return {0, container.rend()}; }
//    };
}