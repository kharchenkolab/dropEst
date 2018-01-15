#include "ReadParametersEfficient.h"

namespace Estimation
{
	ReadParametersEfficient::ReadParametersEfficient(const Tools::ReadParameters &parameters,
	                                                 StringIndexer &barcode_indexer, StringIndexer &umi_indexer,
	                                                 StringIndexer &umi_quality_indexer)
		: _is_empty(parameters.is_empty())
		, _cell_barcode_id(barcode_indexer.add(parameters.cell_barcode()))
		, _umi_id(umi_indexer.add(parameters.umi()))
		, _umi_quality_id(umi_quality_indexer.add(parameters.umi_quality()))
		, _pass_quality_threshold(parameters.pass_quality_threshold())
	{}

	Tools::ReadParameters ReadParametersEfficient::parameters(StringIndexer &barcode_indexer, StringIndexer &umi_indexer,
	                                                          StringIndexer &umi_quality_indexer) const
	{
		return Tools::ReadParameters(barcode_indexer.get_value(this->_cell_barcode_id),
		                             umi_indexer.get_value(this->_umi_id), "",
		                             umi_quality_indexer.get_value(this->_umi_quality_id),
		                             this->_pass_quality_threshold);
	}
}
