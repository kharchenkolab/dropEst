#include "UMI.h"

#include "ReadInfo.h"

#include <stdexcept>

namespace Estimation
{
	UMI::UMI(size_t quality_length, size_t read_count)
		: _read_count(read_count)
		, _mark()
		, _sum_quality(quality_length, 0)
	{}

	void UMI::merge(const UMI &umi)
	{
		this->_read_count += umi._read_count;
		this->_mark.add(umi._mark);
	}

	void UMI::add_read(const ReadInfo &read_info)
	{
		this->_read_count++;
		this->_mark.add(read_info.umi_mark);

		if (read_info.params.umi_quality().length() != this->_sum_quality.size())
			throw std::runtime_error("Wrong quality length: " + std::to_string(read_info.params.umi_quality().length()) +
					                         ", expected: " + std::to_string(this->_sum_quality.size()));

		for (size_t i = 0; i < this->_sum_quality.size(); ++i)
		{
			this->_sum_quality[i] += read_info.params.umi_quality()[i];
		}
	}

	size_t UMI::read_count() const
	{
		return this->_read_count;
	}

	const UMI::Mark &UMI::mark() const
	{
		return this->_mark;
	}

	std::vector<double> UMI::mean_quality() const
	{
		std::vector<double> res(this->_sum_quality.size());
		for (int i = 0; i < res.size(); ++i)
		{
			res[i] = (this->_sum_quality[i] - ReadInfo::quality_offset) / this->_read_count;
		}

		return res;
	}

	UMI::Mark::Mark(MarkType type)
		: _mark(type)
	{}

	void UMI::Mark::add(Mark::MarkType type)
	{
		this->_mark |= type;
	}

	bool UMI::Mark::check(Mark::MarkType type) const
	{
		return this->_mark & type;
	}

	void UMI::Mark::add(const UMI::Mark &mark)
	{
		this->_mark |= mark._mark;
	}

	bool UMI::Mark::match(const std::vector<Mark>match_levels) const
	{
		for (auto const &match_level : match_levels)
		{
			if (this->_mark == match_level._mark)
				return true;
		}

		return false;
	}

	void UMI::Mark::add(Tools::GtfRecord::RecordType type)
	{
		switch (type)
		{
			case Tools::GtfRecord::EXON:
				this->add(HAS_EXONS);
				break;
			case Tools::GtfRecord::INTRON:
				this->add(HAS_INTRONS);
				break;
			default:
				throw std::runtime_error("Unexpected GtfRecord type: " + std::to_string(type));
		}
	}

	bool UMI::Mark::operator==(const UMI::Mark::MarkType &other) const
	{
		return this->_mark == other;
	}

	bool UMI::Mark::operator==(const UMI::Mark &other) const
	{
		return this->_mark == other._mark;
	}

	std::vector<UMI::Mark> UMI::Mark::get_by_code(const std::string &code)
	{
		std::vector<Mark> match_levels;
		for (char c : code)
		{
			match_levels.push_back(Mark::get_by_code(c));
		}

		return match_levels;
	}

	UMI::Mark UMI::Mark::get_by_code(char code)
	{
		Mark mark;
		switch (code)
		{
			case 'e':
				mark.add(HAS_EXONS);
				return mark;
			case 'i':
				mark.add(HAS_INTRONS);
				return mark;
			case 'E':
				mark.add(HAS_EXONS);
				mark.add(HAS_NOT_ANNOTATED);
				return mark;
			case 'I':
				mark.add(HAS_INTRONS);
				mark.add(HAS_NOT_ANNOTATED);
				return mark;
			case 'B':
				mark.add(HAS_EXONS);
				mark.add(HAS_INTRONS);
				return mark;
			case 'A':
				mark.add(HAS_EXONS);
				mark.add(HAS_INTRONS);
				mark.add(HAS_NOT_ANNOTATED);
				return mark;
			default:
				throw std::runtime_error(std::string("Unexpected gene match levels: ") + code);
		}
	}
}
