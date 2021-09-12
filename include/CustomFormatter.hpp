#pragma once

#include <array>
#include <complex>
#include <fmt/format.h>

#include "Definition.hpp"

constexpr auto precisionParser = []<typename ParseContext>(ParseContext& ctx, std::size_t& mPrecision) {
    mPrecision = 0;
    auto it = ctx.begin(), end = ctx.end();
    if (it == end) return it;
    std::size_t index = 0;
    char precision[8];
    if (*it++ == '.') {
        while (*it != 'f') {
            precision[index++] = *it++;
        }
    }
    precision[index] = '\0';
    for (int8_t i = static_cast<int8_t>(strlen(precision)) - 1; i >= 0; i--) {
        mPrecision = 10 * mPrecision + (precision[i] - 48);
    }
    return ++it; 
};

template <FloatingPoint T>
struct fmt::formatter<std::complex<T>> {

    std::size_t mPrecision{};

    template <typename ParseContext>
    constexpr typename ParseContext::iterator parse(ParseContext& ctx) {
        return precisionParser(ctx, mPrecision);
    }

    template <typename FormatContext>
    constexpr auto format(const std::complex<T>& item, FormatContext& ctx) -> decltype(ctx.out()) {
        if (mPrecision == 0) {
            fmt::format_to(ctx.out(), "{}", item.real());
            if (item.imag() != T{0}) {
                return fmt::format_to(ctx.out(), "{:+}i", item.imag());
            }
        } else {
            fmt::format_to(ctx.out(), "{:.{}f}", item.real(), mPrecision);
            if (item.imag() != T{0}) {
                return fmt::format_to(ctx.out(), "{:+.{}f}i", item.imag(), mPrecision);
            }
        }
        return ctx.out();
    }
};

template <typename T, std::size_t N>
struct fmt::formatter<std::array<T, N>> {

    std::size_t mPrecision{};

    template <typename ParseContext>
    constexpr typename ParseContext::iterator parse(ParseContext& ctx) {
        return precisionParser(ctx, mPrecision);
    }

    template <typename FormatContext, std::size_t... I>
    constexpr auto format_impl(const std::array<T, N>& item, FormatContext& ctx, std::index_sequence<I...>) -> decltype(ctx.out()) {
        std::string sExtension{};
        std::string sFormatter{};
        if (mPrecision == 0) {
            sExtension = ", {}";
            sFormatter = "{}";
        } else {
            sExtension = ", {:." + std::to_string(mPrecision) + "f}";
            sFormatter = "{:." + std::to_string(mPrecision) + "f}";
        }
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

    std::size_t mPrecision{};

    template <typename ParseContext>
    constexpr typename ParseContext::iterator parse(ParseContext& ctx) {
        return precisionParser(ctx, mPrecision);
    }

    template <typename FormatContext, std::size_t... I>
    constexpr auto format_impl(const cte::mat::Matrix<T, R, C>& item, FormatContext& ctx, std::index_sequence<I...>) -> decltype(ctx.out()) {
        std::string sExtension{};
        std::string sFormatter{};
        if (mPrecision == 0) {
            sExtension = ",\n {}";
            sFormatter = "{}";
        } else {
            sExtension = ",\n {:." + std::to_string(mPrecision) + "f}";
            sFormatter = "{:." + std::to_string(mPrecision) + "f}";
        }
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

    std::size_t mPrecision{};

    template <typename ParseContext>
    constexpr typename ParseContext::iterator parse(ParseContext& ctx) {
        return precisionParser(ctx, mPrecision);
    }

    template <typename FormatContext, std::size_t... I>
    constexpr auto format_impl(const cte::poly::Polynomial<T, N>& item, FormatContext& ctx, std::index_sequence<I...>) -> decltype(ctx.out()) {
        std::string sExtension{};
        std::string sFormatter{};
        if (mPrecision == 0) {
            sExtension = "{:+}";
            sFormatter = "{}";
        } else {
            sExtension = "{:+." + std::to_string(mPrecision) + "f}";
            sFormatter = "{:." + std::to_string(mPrecision) + "f}";
        }
        std::string extension = std::string(sExtension);
        std::string formatter_string{""};
        std::size_t index = 0;
        formatter_string += fmt::format("{}*x^{}", sFormatter, N - index++);
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

