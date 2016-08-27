#pragma once
#include "BamProcessor.h"
#include <Tools/ReadParameters.h>

namespace Estimation
{
	class ReadMapBamProcessor : public BamProcessor
	{
	private:
		mutable Tools::reads_params_map_t _reads_params;
		const std::string _read_param_names;

	protected:
		virtual bool get_read_params(const BamTools::BamAlignment &alignment, Tools::ReadParameters &read_params) const override;
		virtual void init_temporaries_before_parsing(bool save_read_name) const override;

	public:
		ReadMapBamProcessor(size_t read_prefix_length, const std::string & read_param_names,
		                    const std::string & gtf_path);
	};
}