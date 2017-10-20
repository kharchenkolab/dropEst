#pragma once

#include <string>
#include "Interval.h"

namespace TestTools
{
	struct testGtf;
	struct testGeneMerge;
}

namespace Tools
{
	namespace GeneAnnotation
	{
		class GtfRecord : public Interval
		{
			friend struct TestTools::testGtf;
			friend struct TestTools::testGeneMerge;

		public:
			enum RecordType
			{
				NONE,
				INTRON,
				EXON
			};

		private:
			std::string _chr_name = "";
			std::string _gene_id = "";
			std::string _gene_name = "";
			std::string _transcript_id = "";
			RecordType _type;

		public:
			GtfRecord();

			GtfRecord(const std::string &chr_name, const std::string &gene_id, const std::string &gene_name,
			          coord_t start_pos, coord_t end_pos, RecordType type, const std::string &transcript_id = "");

			const std::string &chr_name() const;

			const std::string &gene_id() const;

			const std::string &gene_name() const;

			const std::string &transcript_id() const;

			RecordType type() const;

			bool is_valid() const;

			bool operator<(const GtfRecord &other) const;
		};
	}
}