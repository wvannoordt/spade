#pragma once
#include <type_traits>

template <const int start, const int finish, typename loopCall>
static inline void static_for(const loopCall& loopObj)
{
	if constexpr (start < finish)
	{
		loopObj(std::integral_constant<int, start>{});
		static_for<start+1, finish>(loopObj);
	}
};