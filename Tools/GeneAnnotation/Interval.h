#pragma once

#include <cstddef>

namespace Tools
{
	namespace GeneAnnotation
	{
		class Interval
		{
		public:
			using coord_t = size_t;

		private:
			coord_t _start_pos;
			coord_t _end_pos;

		public:
			Interval(coord_t start_pos, coord_t end_pos);

			coord_t start_pos() const;

			coord_t end_pos() const;

			coord_t length() const;

			inline bool operator<(const Interval &other) const;

			bool is_intercept(const Interval &other) const;

			virtual void merge(const Interval &other);
		};
	}
}
