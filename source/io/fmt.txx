FMT_EXTERN template void fmt::detail::vformat_to(
    buffer<char>& buf,
    basic_string_view<char> fmt,
    basic_format_args<FMT_BUFFER_CONTEXT(type_identity_t<char>)> args,
    locale_ref loc = {});

FMT_EXTERN template std::string fmt::format(format_string<std::string_view>, std::string_view &&);

FMT_EXTERN template std::string fmt::format(format_string<int>, int &&);
FMT_EXTERN template std::string fmt::format(format_string<long>, long &&);
FMT_EXTERN template std::string fmt::format(format_string<size_t>, size_t &&);
FMT_EXTERN template std::string fmt::format(format_string<double>, double &&);

FMT_EXTERN template std::string fmt::format(format_string<int,int>, int &&, int &&);
FMT_EXTERN template std::string fmt::format(format_string<long,long>, long &&, long &&);
FMT_EXTERN template std::string fmt::format(format_string<size_t,size_t>, size_t &&, size_t &&);
FMT_EXTERN template std::string fmt::format(format_string<double,double>, double &&, double &&);

FMT_EXTERN template std::string fmt::format(format_string<int, int, int>, int &&, int &&, int &&);
FMT_EXTERN template std::string fmt::format(format_string<long, long, long>, long &&, long &&, long &&);
FMT_EXTERN template std::string fmt::format(format_string<size_t, size_t, size_t>, size_t &&, size_t &&, size_t &&);
FMT_EXTERN template std::string fmt::format(format_string<double, double, double>, double &&, double &&, double &&);

