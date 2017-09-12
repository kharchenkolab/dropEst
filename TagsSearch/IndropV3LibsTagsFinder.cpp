#include <Tools/UtilFunctions.h>
#include "IndropV3LibsTagsFinder.h"

namespace TagsSearch
{

	IndropV3LibsTagsFinder::IndropV3LibsTagsFinder(const std::string &barcode1_fastq_name, const std::string &barcode2_fastq_name,
	                                               const std::string &gene_fastq_name, const std::string &library_fastq_name,
	                                               const std::string &library_tag,
		                                           const boost::property_tree::ptree &barcodes_config,
		                                           const boost::property_tree::ptree &config, bool save_stats)
		: IndropV3TagsFinder(barcode1_fastq_name, barcode2_fastq_name, gene_fastq_name, barcodes_config, config, save_stats)
		, library_tag(library_tag)
		, max_lib_tag_ed(barcodes_config.get<unsigned>("max_libtag_ed", 2))
		, _lib_tag_reader(library_fastq_name)
	{}

	bool IndropV3LibsTagsFinder::parse_fastq_record(FastQReader::FastQRecord &record, Tools::ReadParameters &read_params)
	{
		auto lib_record = this->_lib_tag_reader.get_next_record();
		if (Tools::edit_distance(lib_record.sequence.c_str(), this->library_tag.c_str(), false) > this->max_lib_tag_ed)
		{
			this->skip_records_row();
			read_params = Tools::ReadParameters();
			return true;
		}

		return IndropV3TagsFinder::parse_fastq_record(record, read_params);
	}
}