#pragma once
#include "BamProcessor.h"

namespace Estimation
{
	class ReadMapBamProcessor : public BamProcessor
	{
	private:
		Tools::reads_params_map_t _reads_params;

	private:
		void fill_names_map(const std::string &reads_params_names_str);

	protected:
		virtual bool get_read_params(const BamTools::BamAlignment &alignment, Tools::ReadParameters &read_params) const override;

	public:
		ReadMapBamProcessor(size_t read_prefix_length, const std::string & reads_params_names_str,
		                    const std::string & gtf_path);
	};
}