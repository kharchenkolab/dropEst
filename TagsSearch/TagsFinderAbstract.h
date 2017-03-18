#pragma once

namespace TagsSearch
{
	class TagsFinderAbstract
	{
	public:
		virtual void run(bool save_read_names) = 0;
	};
}