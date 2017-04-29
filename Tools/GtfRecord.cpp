#include "GtfRecord.h"

namespace Tools
{
	GtfRecord::GtfRecord(const std::string &chr_name, const std::string &gene_id, const std::string &gene_name,
	                   coord_t start_pos, coord_t end_pos, RecordType type, const std::string &transcript_id)
		: Interval(start_pos, end_pos)
		, _chr_name(chr_name)
		, _gene_id(gene_id)
		, _gene_name(gene_name == gene_id ? "" : gene_name)
		, _transcript_id(transcript_id)
		, _type(type)
	{}

	bool GtfRecord::is_valid() const
	{
		return this->gene_id() != "";
	}

	bool GtfRecord::operator<(const GtfRecord &other) const
	{
		return this->gene_id() < other.gene_id();
	}

	const std::string& GtfRecord::chr_name() const
	{
		return this->_chr_name;
	}

	const std::string& GtfRecord::gene_id() const
	{
		return this->_gene_id;
	}

	const std::string& GtfRecord::gene_name() const
	{
		return this->_gene_name == "" ? this->_gene_id : this->_gene_name;
	}

	GtfRecord::GtfRecord()
		: Interval(0, 0)
	{}

	const std::string &GtfRecord::transcript_id() const
	{
		return this->_transcript_id == "" ? this->_gene_id : this->_transcript_id;
	}

	GtfRecord::RecordType GtfRecord::type() const
	{
		return this->_type;
	}
}