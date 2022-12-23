#pragma once
//The below is only temporary until NVIDIA gets their act together with
//c++20
#define _c20 1
#if(_c20)
#include <concepts>
#include <filesystem>
#define _c20_concept(mycode) mycode
#else
#include <experimental/filesystem>
#define _c20_concept(mycode) typename
#endif