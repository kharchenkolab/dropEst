#include <Tools/UtilFunctions.h>
#include "IndropV3LibsTagsFinder.h"

namespace TagsSearch
{

	IndropV3LibsTagsFinder::IndropV3LibsTagsFinder(const std::vector<std::string> &fastq_filenames,
	                                               const std::string &library_tag,
		                                           const boost::property_tree::ptree &barcodes_config,
		                                           const boost::property_tree::ptree &config,
		                                           const std::shared_ptr<ConcurrentGzWriter> &writer,
		                                           bool save_stats, bool save_read_params)
		: IndropV3TagsFinder(fastq_filenames, barcodes_config, config, writer, save_stats, save_read_params)
		, library_tag(library_tag)
		, max_lib_tag_ed(barcodes_config.get<unsigned>("max_libtag_ed", 2))
	{}

	bool IndropV3LibsTagsFinder::parse_fastq_records(FastQReader::FastQRecord &record,
													 Tools::ReadParameters &read_params)
	{
		FastQReader::FastQRecord lib_record;
		if (!this->fastq_reader(3).get_next_record(lib_record))
			return false;

		if (Tools::edit_distance(lib_record.sequence.c_str(), this->library_tag.c_str(), false) > this->max_lib_tag_ed)
		{
			for (size_t file_id = 0; file_id < 3; ++file_id)
			{
				this->fastq_reader(file_id).get_next_record(lib_record);
			}
			read_params = Tools::ReadParameters();
			return true;
		}

		return IndropV3TagsFinder::parse_fastq_records(record, read_params);
	}
}