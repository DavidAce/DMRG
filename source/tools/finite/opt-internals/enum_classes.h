//
// Created by david on 2019-10-03.
//

#include <iostream>
#include <sstream>

namespace tools::finite::opt{
    namespace internal{
        enum class EnumType  {REAL,CPLX};
        enum class EnumMode  {OVERLAP,VARIANCE};
        enum class EnumSpace {SUBSPACE, DIRECT};
        template <typename EnumClass> class EnumBase;
    }
    using TYPE  = internal::EnumType;
    using MODE  = internal::EnumMode;
    using SPACE = internal::EnumSpace;

    class OptSpace;
    class OptType ;
    class OptMode ;
}


template <typename EnumClass>
class tools::finite::opt::internal::EnumBase{
protected:
    virtual void print(std::ostream& str) const = 0;
public:
    EnumClass option;
    explicit EnumBase(const EnumClass &enumVal):option(enumVal){};

    void operator= (EnumClass enumVal){
        option = enumVal;
    }

    bool operator== (const EnumClass &enumVal) const{
        return option == enumVal;
    }

    friend std::ostream& operator<<(std::ostream& str, const EnumBase<EnumClass> &base){
        base.print(str);
        return str;
    }
    [[nodiscard]] std::string str() const{
        std::stringstream ss;
        ss << *this;
        return ss.str();
    }

    explicit operator std::string() const {
        return str();
    }
};


class tools::finite::opt::OptType : public tools::finite::opt::internal::EnumBase<opt::TYPE>{
public:
    using EnumBase::EnumBase;
    using EnumBase::operator=;
    static constexpr opt::TYPE REAL = opt::TYPE::REAL ;
    static constexpr opt::TYPE CPLX = opt::TYPE::CPLX;
private:
    void print(std::ostream& str) const final;
};


class tools::finite::opt::OptMode : public tools::finite::opt::internal::EnumBase<opt::MODE>{
public:
    using EnumBase::EnumBase;
    using EnumBase::operator=;
    static constexpr opt::MODE OVERLAP    = opt::MODE::OVERLAP ;
    static constexpr opt::MODE VARIANCE   = opt::MODE::VARIANCE;
private:
    void print(std::ostream& str) const final;
};



class tools::finite::opt::OptSpace : public tools::finite::opt::internal::EnumBase<opt::SPACE>{
public:
    using EnumBase::EnumBase;
    using EnumBase::operator=;
    static constexpr opt::SPACE SUBSPACE = opt::SPACE::SUBSPACE;
    static constexpr opt::SPACE DIRECT   = opt::SPACE::DIRECT;
private:
    void print(std::ostream& str) const final;
};

