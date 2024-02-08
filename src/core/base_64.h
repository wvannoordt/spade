#pragma once

#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <cstring>
#include "core/print.h"

namespace spade::detail
{
    template <class data_t> static inline void stream_base_64(std::ostream& strm, const data_t* data, const size_t& data_size)
    {
        //This is a very slow implentation for now.
        const char  table[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
        char term = '=';
        char b0, b1, b2, c0, c1, c2, c3;
        std::size_t size = data_size*sizeof(data_t);
        size_t numWriteTotal = (size - size%3)/3;
        char* bytes = (char*)data;
        for (size_t pos = 0; pos < numWriteTotal; pos++)
        {
            b0 = bytes[3*pos];
            b1 = bytes[3*pos+1];
            b2 = bytes[3*pos+2];
            c0 = ((b0 & 0b11111100) >> 2);
            c1 = ((b0 & 0b00000011) << 4) | ((b1 & 0b11110000) >> 4);
            c2 = ((b1 & 0b00001111) << 2) | ((b2 & 0b11000000) >> 6);
            c3 = ((b2 & 0b00111111));
            strm << table[c0] << table[c1] << table[c2] << table[c3];
        }
        int numLeft = size%3;
        if (numLeft==0) return;
        char last3[3] = {0};
        char last4[4] = {0};
        int num2write[3] = {0, 2, 3};
        int numTerms[3] = {0, 2, 1};
        for (int i = 0; i < numLeft; i++)
        {
            last3[i] = bytes[size-(size%3)+i];
        }
        b0 = last3[0];
        b1 = last3[1];
        b2 = last3[2];
        c0 = ((b0 & 0b11111100) >> 2);
        c1 = ((b0 & 0b00000011) << 4) | ((b1 & 0b11110000) >> 4);
        c2 = ((b1 & 0b00001111) << 2) | ((b2 & 0b11000000) >> 6);
        c3 = ((b2 & 0b00111111));
        last4[0] = c0;
        last4[1] = c1;
        last4[2] = c2;
        last4[3] = c3;
        for (int i = 0; i < num2write[numLeft]; i++)
        {
            strm << table[last4[i]];
        }
        for (int i = 0; i < numTerms[numLeft]; i++)
        {
            strm << term;
        }
    }
    
    template<typename data_t, typename alloc_t>
    static inline void stream_base_64(std::ostream& strm, const std::vector<data_t, alloc_t>& data)
    {
        unsigned int bytes = data.size() * sizeof(data_t);
        stream_base_64(strm, &bytes, 1);
        stream_base_64(strm, data.data(), data.size());
    }
}