#pragma once

#include <filesystem>

namespace spade::io
{
    static void mkdir(const std::string& dir)
    {
        if (!std::filesystem::is_directory(dir)) std::filesystem::create_directory(dir);
    };
}