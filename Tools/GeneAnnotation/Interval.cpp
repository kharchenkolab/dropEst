#include "Interval.h"

#include <algorithm>

namespace Tools
{
	namespace GeneAnnotation
	{
		Interval::Interval(coord_t start_pos, coord_t end_pos)
				: _start_pos(start_pos)
				, _end_pos(end_pos)
		{}

		bool Interval::is_intercept(const Interval &other) const
		{
			return this->start_pos() <= other.end_pos() && this->end_pos() > other.start_pos();
		}

		void Interval::merge(const Interval &other)
		{
			this->_start_pos = std::min(this->_start_pos, other._start_pos);
			this->_end_pos = std::max(this->_end_pos, other._end_pos);
		}

		bool Interval::operator<(const Interval &other) const
		{
			if (this->start_pos() == other.start_pos())
				return this->end_pos() < other.end_pos();

			return this->start_pos() < other.start_pos();
		}

		Interval::coord_t Interval::length() const
		{
			return this->end_pos() - this->start_pos();
		}

		Interval::coord_t Interval::start_pos() const
		{
			return this->_start_pos;
		}

		Interval::coord_t Interval::end_pos() const
		{
			return this->_end_pos;
		}
	}
}