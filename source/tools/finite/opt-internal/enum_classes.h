#pragma once
#include <iostream>
#include <sstream>



namespace tools::finite::opt{
    class OptType{
        private:
        enum class EnumType  {REAL,CPLX};
        EnumType enum_;
        public:
        static constexpr auto REAL = EnumType::REAL;
        static constexpr auto CPLX = EnumType::CPLX;
        OptType(const EnumType & otherEnum): enum_(otherEnum){}
        OptType(const OptType & otherOpt) = default;
        operator EnumType(){return enum_;}
        [[nodiscard]] std::string string() const {
            switch (enum_){
                case EnumType::REAL : return "REAL";
                case EnumType::CPLX : return "CPLX";
            }
            return "ERROR";
        }

        explicit operator std::string() const {
            return this->string();
        }
//        explicit             operator std::string_view() const { return this->string(); }

        friend std::ostream& operator<<(std::ostream& str, const OptType & optType){
            return str << optType.string();
        }


        OptType & operator= (const EnumType & other){
            this->enum_ =other;
            return *this;
        }
        bool operator== (const EnumType &other) const{
            return enum_ == other;
        }
        bool operator== (const OptType & other) const{
            return enum_ == other.enum_;
        }
    };

    class OptMode{
        private:
        enum class EnumMode  {OVERLAP,VARIANCE};
        EnumMode enum_;
        public:
        static constexpr auto OVERLAP  = EnumMode::OVERLAP;
        static constexpr auto VARIANCE = EnumMode::VARIANCE;
        OptMode(const EnumMode & otherEnum): enum_(otherEnum){}
        OptMode(const OptMode & otherOpt) = default;
        operator EnumMode(){return enum_;}
        [[nodiscard]] std::string string() const {
            switch (enum_){
                case EnumMode::OVERLAP  : return "OVERLAP";
                case EnumMode::VARIANCE : return "VARIANCE";
            }
            return "ERROR";
        }
        explicit operator std::string() const {
            return this->string();
        }
//        explicit             operator std::string_view() const { return this->string(); }

        friend std::ostream& operator<<(std::ostream& str, const OptMode & optMode){
            return str << optMode.string();
        }

        OptMode & operator= (const EnumMode & other){
            this->enum_ =other;
            return *this;
        }

        bool operator== (const EnumMode &other) const{
            return enum_ == other;
        }
        bool operator== (const OptMode & other) const{
            return enum_ == other.enum_;
        }

    };

    class OptSpace{
        private:
        enum class EnumSpace : int {SUBSPACE_ONLY, DIRECT, SUBSPACE_AND_DIRECT};
        EnumSpace enum_;
        public:
        static constexpr auto SUBSPACE_ONLY  = EnumSpace::SUBSPACE_ONLY;
        static constexpr auto SUBSPACE_AND_DIRECT = EnumSpace::SUBSPACE_AND_DIRECT;
        static constexpr auto DIRECT = EnumSpace::DIRECT;
        OptSpace(const EnumSpace &other): enum_(other){}
        operator EnumSpace(){return enum_;}
        OptSpace(const OptSpace & otherOpt) = default;
        [[nodiscard]] std::string string() const {
            switch (enum_){
                case EnumSpace::SUBSPACE_ONLY : return "SUBSPACE_ONLY";
                case EnumSpace::SUBSPACE_AND_DIRECT : return "SUBSPACE_AND_DIRECT";
                case EnumSpace::DIRECT : return "DIRECT";
            }
            return "ERROR";
        }

        explicit operator std::string() const {
            return this->string();
        }
//        explicit             operator std::string_view() const { return this->string(); }
        friend std::ostream& operator<<(std::ostream& str, const OptSpace & optSpace){
            return str << optSpace.string();
        }

        OptSpace & operator= (const EnumSpace &other){
            this->enum_ = other;
            return *this;
        }
        bool operator== (const EnumSpace &other) const{
            return enum_ == other;
        }
        bool operator== (const OptSpace & other) const{
            return enum_ == other.enum_;
        }
    };
}


//
//
//namespace tools::finite::opt{
//    namespace internal{
//        enum class EnumType  {REAL,CPLX};
//        enum class EnumMode  {OVERLAP,VARIANCE};
//        enum class EnumSpace {SUBSPACE_ONLY, DIRECT, SUBSPACE_DIRECT};
//        template <typename EnumClass> class EnumBase;
//    }
//    using TYPE  = internal::EnumType;
//    using MODE  = internal::EnumMode;
//    using SPACE = internal::EnumSpace;
//
//    class OptSpace;
//    class OptType ;
//    class OptMode ;
//}

//
//template <typename EnumClass>
//class tools::finite::opt::internal::EnumBase{
//protected:
//    virtual void print(std::ostream& str) const = 0;
//public:
//    EnumClass option;
//    explicit EnumBase(const EnumClass &enumVal):option(enumVal){};
//
//    void operator= (EnumClass enumVal){
//        option = enumVal;
//    }
//
//    bool operator== (const EnumClass &enumVal) const{
//        return option == enumVal;
//    }
//
//    friend std::ostream& operator<<(std::ostream& str, const EnumBase<EnumClass> &base){
//        base.print(str);
//        return str;
//    }
//    [[nodiscard]] std::string str() const{
//        std::stringstream ss;
//        ss << *this;
//        return ss.str();
//    }
//    operator EnumClass() {return option;}
//    explicit operator std::string() const {
//        return str();
//    }
//};
//
//
//class tools::finite::opt::OptType : public tools::finite::opt::internal::EnumBase<opt::TYPE>{
//public:
//    using EnumBase::EnumBase;
//    using EnumBase::operator=;
//    static constexpr opt::TYPE REAL = opt::TYPE::REAL ;
//    static constexpr opt::TYPE CPLX = opt::TYPE::CPLX;
//private:
//    void print(std::ostream& str) const final;
//};
//
//
//class tools::finite::opt::OptMode : public tools::finite::opt::internal::EnumBase<opt::MODE>{
//public:
//    using EnumBase::EnumBase;
//    using EnumBase::operator=;
//    static constexpr opt::MODE OVERLAP    = opt::MODE::OVERLAP ;
//    static constexpr opt::MODE VARIANCE   = opt::MODE::VARIANCE;
//private:
//    void print(std::ostream& str) const final;
//};
//
//
//
//class tools::finite::opt::OptSpace : public tools::finite::opt::internal::EnumBase<opt::SPACE>{
//public:
//    using EnumBase::EnumBase;
//    using EnumBase::operator=;
//    static constexpr opt::SPACE SUBSPACE = opt::SPACE::SUBSPACE_ONLY;
//    static constexpr opt::SPACE DIRECT   = opt::SPACE::DIRECT;
//private:
//    void print(std::ostream& str) const final;
//};

