#pragma once

#include <string>
#include <exception>

#include "Definition.hpp"

namespace cte {

template<Integer IntegerType>
struct DimensionMismatchException : std::exception {

    DimensionMismatchException(const std::string& message, const IntegerType dim1, const IntegerType dim2) 
        : eMessage(message), eDim1(dim1), eDim2(dim2) {}

    virtual const char* what() const noexcept override {
        return (eMessage + ": " + std::to_string(eDim1) + " and " + std::to_string(eDim2)).c_str();
    }

private:
    std::string eMessage;
    IntegerType eDim1, eDim2;
};

template<Integer IntegerType>
struct OutOfRangeException : std::exception {

    OutOfRangeException(const IntegerType index) 
        : eIndex(index) {}

    virtual const char* what() const noexcept override {
        return (eMessage + std::to_string(eIndex)).c_str();
    }

private:
    IntegerType eIndex;
    inline static std::string eMessage = "Index out of range: ";
};

struct CustomMessageException : std::exception {

    CustomMessageException(const std::string& message) 
        : eMessage(message) {}

    virtual const char* what() const noexcept override {
        return eMessage.c_str();
    }

private:
    std::string eMessage;
};

}

