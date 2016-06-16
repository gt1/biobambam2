/*
    libmaus2
    Copyright (C) 2016 German Tischler

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <libmaus2/bambam/BamHeader.hpp>
#include <libmaus2/util/LineBuffer.hpp>
#include <libmaus2/bambam/SamInfo.hpp>
#include <libmaus2/util/ArgParser.hpp>

int main(int argc, char *argv[])
{
	try
	{
		libmaus2::util::ArgParser const arg(argc,argv);

		libmaus2::util::LineBuffer LB(std::cin);
		char const * a;
		char const * e;

		std::ostringstream headerostr;
		bool haveheader = false;
		libmaus2::bambam::BamHeader::unique_ptr_type Pheader;
		libmaus2::bambam::SamInfo::unique_ptr_type Pinfo;

		while ( LB.getline(&a,&e) )
		{
			if ( e-a && a[0] == '@' )
			{
				if ( ! haveheader )
				{
					headerostr << std::string(a,e) << "\n";
				}
				else
				{
					libmaus2::exception::LibMausException lme;
					lme.getStream() << "[E] header line in data part" << std::endl;
					lme.finish();
					throw lme;
				}
			}
			else if ( e-a )
			{
				if ( ! Pinfo )
				{
					libmaus2::bambam::BamHeader::unique_ptr_type Theader(new libmaus2::bambam::BamHeader(headerostr.str()));
					Pheader = UNIQUE_PTR_MOVE(Theader);
					libmaus2::bambam::SamInfo::unique_ptr_type Tinfo(new libmaus2::bambam::SamInfo(*Pheader));
					Pinfo = UNIQUE_PTR_MOVE(Tinfo);

					std::cout << headerostr.str();
				}

				try
				{
					Pinfo->parseSamLine(a,e);

					if ( Pinfo->algn.isMapped() )
					{
						if (
							static_cast<int64_t>(Pinfo->algn.getPos() + Pinfo->algn.getReferenceLength()) >
							static_cast<int64_t>(Pheader->getRefIDLength(Pinfo->algn.getRefID()))
						)
						{
							libmaus2::exception::LibMausException lme;
							lme.getStream() << "[E] alignment extends past end of reference sequence" << std::endl;
							lme.finish();
							throw lme;
						}
					}

					std::cout << std::string(a,e) << "\n";
				}
				catch(std::exception const & ex)
				{
					std::cerr << ex.what() << std::endl;
				}
			}
		}
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
	}
}
