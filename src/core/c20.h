#pragma once

//If we are using nvcc, we need to disable c++20, sad times
#ifdef __NVCC__
#define _c20 0
#else
#define _c20 1
#endif


//The below is only temporary until NVIDIA gets their act together with
//c++20
#if(_c20)
#include <concepts>
#include <filesystem>
#define _c20_concept(mycode) mycode
#else
#include <experimental/filesystem>
#define _c20_concept(mycode) typename
#endif