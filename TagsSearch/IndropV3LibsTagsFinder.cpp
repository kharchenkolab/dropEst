#include <Tools/UtilFunctions.h>
#include "IndropV3LibsTagsFinder.h"

namespace TagsSearch
{

	IndropV3LibsTagsFinder::IndropV3LibsTagsFinder(const std::shared_ptr<FilesProcessor> &files_processor, const std::string &library_tag,
		                                               unsigned max_lib_tag_ed, const boost::property_tree::ptree &barcodes_config,
		                                               const boost::property_tree::ptree &config)
		: IndropV3TagsFinder(files_processor, barcodes_config, config)
		, library_tag(library_tag)
		, max_lib_tag_ed(max_lib_tag_ed)
	{}

	bool
	IndropV3LibsTagsFinder::parse_fastq_record(FilesProcessor::FastQRecord &record, Tools::ReadParameters &read_params)
	{
		auto lib_record = this->_files_processor->get_fastq_record(3);
		if (Tools::edit_distance(lib_record.sequence.c_str(), this->library_tag.c_str(), false) > this->max_lib_tag_ed)
		{
			for (size_t i = 0; i < 3; ++i)
			{
				this->_files_processor->get_fastq_record(i);
			}
			read_params = Tools::ReadParameters();
			return true;
		}

		return IndropV3TagsFinder::parse_fastq_record(record, read_params);
	}
}