#pragma once

#include "Tools/Logs.h"
#ifdef R_LIBS
#include <RcppArmadillo.h>
#endif

namespace Estimation
{
	namespace Results
	{
		class IRTableProvider
		{
		protected:
		#ifdef R_LIBS
			virtual void save_r_table(const std::string &filename) const = 0;
		#endif
		public:
			void save_rds(const std::string &filename) const
			{
			#ifdef R_LIBS
				save_r_table(filename);
			#else
				L_ERR << "Can't print rds without RCpp";
			#endif
			}

		};
	}
}