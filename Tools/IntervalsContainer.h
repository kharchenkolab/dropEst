#pragma once

#include <algorithm>
#include <list>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
#include <map>
#include "Interval.h"

namespace TestTools
{
	struct testGeneMerge;
	struct testParseBed;
	struct testInitGtf;
}

namespace Tools
{
template<class IntervalLabelT> class IntervalsContainer
{
	friend struct TestTools::testGeneMerge;
	friend struct TestTools::testParseBed;
	friend struct TestTools::testInitGtf;

public:
	typedef size_t coord_t;
	typedef std::set<IntervalLabelT> interval_labels_t;

private:
	/// Intervals are stored in 0-based coordinate system
	class BaseIntervalInfo : public Interval
	{
	private:
		const IntervalLabelT _label;

	public:
		BaseIntervalInfo(coord_t start_pos, coord_t end_pos, IntervalLabelT label)
			: Interval(start_pos, end_pos)
			, _label(label)
		{}

		const IntervalLabelT& label() const
		{
			return this->_label;
		}

		inline bool operator<(const BaseIntervalInfo& other) const
		{
			if (Interval::operator<(other))
				return true;

			if (static_cast<Interval>(other).operator<(*this))
				return false;

			return this->label() < other.label();
		}
	};

	class QueryInterval : public Interval
	{
	public:
		const interval_labels_t base_interval_labels;

		QueryInterval(coord_t start_pos, coord_t end_pos, const interval_labels_t &base_interval_labels)
			: Interval(start_pos, end_pos)
			, base_interval_labels(base_interval_labels)
		{}
	};

	struct Event
	{
		enum EventType
		{
			OPEN,
			CLOSE
		};
		const IntervalLabelT* label;
		EventType type;

		Event(const IntervalLabelT *label, IntervalsContainer::Event::EventType type)
			: label(label)
			, type(type)
		{}
	};

	typedef std::multimap<coord_t, Event> events_t;

private:
	bool _initialized;

	const bool _allow_intercepts;
	const unsigned _min_interval_len;

	std::vector<QueryInterval> _homogenous_intervals; // Continuous parts with the same labels composition
	std::map<IntervalLabelT, std::list<BaseIntervalInfo>> _base_intervals;

private:
	events_t convert_intervals_to_events(const std::vector<BaseIntervalInfo> &intervals) const
	{
		events_t events;
		for (auto const &interval_info : intervals)
		{
			events.emplace(interval_info.start_pos(), Event(&interval_info.label(), Event::OPEN));
			events.emplace(interval_info.end_pos(), Event(&interval_info.label(), Event::CLOSE));
		}

		return events;
	}

	void extract_homogeneous_intervals(const events_t &events)
	{
		this->_homogenous_intervals.clear();

		coord_t start_pos = 0, end_pos;
		interval_labels_t cur_labels;
		for (auto const &event : events)
		{
			end_pos = event.first;
			if (!cur_labels.empty() && end_pos - start_pos >= this->_min_interval_len)
			{
				if (!this->_allow_intercepts && cur_labels.size() > 1)
					throw std::runtime_error("Intervals intersection at (" + std::to_string(start_pos) + ", " + std::to_string(end_pos) + ")");

				this->_homogenous_intervals.push_back(QueryInterval(start_pos, end_pos, cur_labels));
			}

			if (event.second.type == Event::OPEN)
			{
				cur_labels.insert(*event.second.label);
			}
			else
			{
				cur_labels.erase(*event.second.label);
			}

			start_pos = end_pos;
		}
	}

public:
	IntervalsContainer(bool allow_intercepts = true, unsigned min_interval_length = 1)
		: _initialized(false)
		, _allow_intercepts(allow_intercepts)
		, _min_interval_len(min_interval_length)
	{}

	void add_interval(coord_t start, coord_t end, IntervalLabelT interval_label, bool force=false)
	{
		if (!force && this->_initialized)
			throw std::runtime_error("IntervalsContainer is already initialized");

		BaseIntervalInfo query(start, end, interval_label);

		auto &cur_intervals = this->_base_intervals[interval_label];
		auto cur_iterator = cur_intervals.begin();
		while (cur_iterator != cur_intervals.end() && !query.is_intercept(*cur_iterator))
		{
			if (cur_iterator->start_pos() > query.end_pos())
			{
				cur_intervals.insert(cur_iterator, query);
				return;
			}
			++cur_iterator;
		}

		if (cur_iterator == cur_intervals.end())
		{
			cur_intervals.push_back(query);
			return;
		}

		auto end_iterator = cur_iterator;
		end_iterator++;

		while (end_iterator != cur_intervals.end() && query.is_intercept(*end_iterator))
		{
			query.merge(*end_iterator);
			++end_iterator;
		}

		cur_iterator->merge(query);
		cur_iterator++;
		cur_intervals.erase(cur_iterator, end_iterator);
	}

	void set_initialized(bool clear=true)
	{
		this->_initialized = true;

		std::vector<BaseIntervalInfo> intervals;
		auto back_inserter = std::back_inserter(intervals);
		for (auto &cur_label : this->_base_intervals)
		{
			std::move(cur_label.second.begin(), cur_label.second.end(), back_inserter);
		}

		auto events = this->convert_intervals_to_events(intervals);
		this->extract_homogeneous_intervals(events);

		if (clear)
		{
			this->_base_intervals.clear();
		}
	}

	/// Return intervals, which intersect [start_pos; end_pos)
	/// \param start_pos 0-based start position (inclusive)
	/// \param end_pos 0-based end position (exclusive)
	/// \return labels of intersected intervals
	interval_labels_t get_intervals(coord_t start_pos, coord_t end_pos) const
	{
		if (!this->_initialized)
			throw std::runtime_error("Interval must be initialized");

		auto intercept_it = std::lower_bound(this->_homogenous_intervals.begin(), this->_homogenous_intervals.end(), start_pos,
		                                     [](const QueryInterval &interval, coord_t pos){ return interval.end_pos() <= pos;});

		if (intercept_it == this->_homogenous_intervals.end() || intercept_it->start_pos() >= end_pos)
			return interval_labels_t();

		interval_labels_t res_labels;
		while (intercept_it != this->_homogenous_intervals.end() && intercept_it->start_pos() < end_pos)
		{
			res_labels.insert(intercept_it->base_interval_labels.begin(), intercept_it->base_interval_labels.end());
			++intercept_it;
		}

		return res_labels;
	}
};
}