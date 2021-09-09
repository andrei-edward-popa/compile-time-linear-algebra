#pragma once

#include <array>
#include <complex>
#include <fmt/format.h>

#include "Definition.hpp"

template <FloatingPoint T>
struct fmt::formatter<std::complex<T>> {

    template <typename ParseContext>
    constexpr typename ParseContext::iterator parse(ParseContext& ctx) {
        return ctx.begin(); 
    }

    template <typename FormatContext>
    constexpr auto format(const std::complex<T>& item, FormatContext& ctx) -> decltype(ctx.out()) {
        if (item.imag() != T{0}) {
            return fmt::format_to(ctx.out(), "{:+.4f}{:+.4f}*i", item.real(), item.imag());
        } else {
            return fmt::format_to(ctx.out(), "{:+.4f}", item.real());
        }
    }
};

template <typename T, std::size_t N>
struct fmt::formatter<std::array<T, N>> {

    static constexpr const char* sExtension = std::is_floating_point_v<T> ? ", {:.4f}" : ", {}";
    static constexpr const char* sFormatter = std::is_floating_point_v<T> ? "{:.4f}" : "{}";

    template <typename ParseContext>
    constexpr typename ParseContext::iterator parse(ParseContext& ctx) {
        return ctx.begin(); 
    }

    template <typename FormatContext, std::size_t... I>
    constexpr auto format_impl(const std::array<T, N>& item, FormatContext& ctx, std::index_sequence<I...>) -> decltype(ctx.out()) {
        std::string extension = std::string(sExtension);
        std::string formatter_string = std::string(sFormatter);
        std::size_t index = 0;
        while(index++ < item.size() - 1) {
            formatter_string += extension;
        }
        return format_to(ctx.out(), fmt::format("[{}]", formatter_string), item[I]...);
    }

    template <typename FormatContext, typename Indeces = std::make_index_sequence<N>>
    constexpr auto format(const std::array<T, N>& item, FormatContext& ctx) -> decltype(ctx.out()) {
        if (item.empty()) {
            return fmt::format_to(ctx.out(), "[]");
        } else {
            return format_impl(item, ctx, Indeces{});
        }
    }
};

template <FloatingPoint T, std::size_t R, std::size_t C>
struct fmt::formatter<cte::mat::Matrix<T, R, C>> {

    static constexpr const char* sExtension = ",\n {}";
    static constexpr const char* sFormatter = "{}";

    template <typename ParseContext>
    typename ParseContext::iterator parse(ParseContext& ctx) {
        return ctx.begin(); 
    }

    template <typename FormatContext, std::size_t... I>
    constexpr auto format_impl(const cte::mat::Matrix<T, R, C>& item, FormatContext& ctx, std::index_sequence<I...>) -> decltype(ctx.out()) {
        std::string extension = std::string(sExtension);
        std::string formatter_string = std::string(sFormatter);
        std::size_t index = 0;
        while(index++ < R - 1) {
            formatter_string += extension;
        }
        return format_to(ctx.out(), fmt::format("[{}]", formatter_string), item[I]...);
    }

    template <typename FormatContext, typename Indeces = std::make_index_sequence<R>>
    constexpr auto format(const cte::mat::Matrix<T, R, C>& item, FormatContext& ctx) -> decltype(ctx.out()) {
        if (R == 0 || C == 0) {
            return fmt::format_to(ctx.out(), "[]");
        } else {
            return format_impl(item, ctx, Indeces{});
        }
    }
};

template <FloatingPoint T, std::size_t N>
struct fmt::formatter<cte::poly::Polynomial<T, N>> {

    static constexpr const char* sExtension = std::is_floating_point_v<T> ? "{:+.4f}" : "{:+}";

    template <typename ParseContext>
    constexpr typename ParseContext::iterator parse(ParseContext& ctx) {
        return ctx.begin(); 
    }

    template <typename FormatContext, std::size_t... I>
    constexpr auto format_impl(const cte::poly::Polynomial<T, N>& item, FormatContext& ctx, std::index_sequence<I...>) -> decltype(ctx.out()) {
        std::string extension = std::string(sExtension);
        std::string formatter_string{""};
        std::size_t index = 0;
        while(index <= item.getDegree()) {
            formatter_string += fmt::format("{}*x^{}", extension, N - index++);
        }
        return fmt::format_to(ctx.out(), fmt::format("{}", formatter_string), item[I]...);
    }

    template <typename FormatContext, typename Indeces = std::make_index_sequence<N + 1>>
    constexpr auto format(const cte::poly::Polynomial<T, N>& item, FormatContext& ctx) -> decltype(ctx.out()) {
        return format_impl(item, ctx, Indeces{});
    }
};

