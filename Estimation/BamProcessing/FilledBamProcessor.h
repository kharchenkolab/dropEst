#pragma once
#include "BamProcessor.h"

namespace Estimation
{
namespace BamProcessing
{
	class FilledBamProcessor : public BamProcessor
	{
	private:

	protected:
		virtual bool get_read_params(const BamTools::BamAlignment &alignment, Tools::ReadParameters &read_params) const override;

	public:
		FilledBamProcessor(size_t read_prefix_length, const std::string & gtf_path);
	};
}
}
