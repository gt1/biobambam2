/**
    biobambam
    Copyright (C) 2009-2013 German Tischler
    Copyright (C) 2011-2013 Genome Research Limited

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
**/
#include <biobambam/KmerPoisson.hpp>
#include <cmath>

double kmerPoisson(
	uint64_t const L, // length of reference sequence
	uint64_t const KA, // number of A bases in query sequence
	uint64_t const KC, // number of C bases in query sequence
	uint64_t const KG, // number of G bases in query sequence
	uint64_t const KT, // number of T bases in query sequence
	uint64_t const n, // number of occurences of query sequence
	double const pA,
	double const pC,
	double const pG,
	double const pT
)
{
	uint64_t const K = KA+KC+KG+KT;
	
	if ( K > L )
		return 0.0;

	double const lambda = (L+1-K) * ::std::pow(pA,static_cast<double>(KA)) * ::std::pow(pC,static_cast<double>(KC)) * ::std::pow(pG,static_cast<double>(KG)) * ::std::pow(pT,static_cast<double>(KT));
	
	double a = 1.0;
	double b = 1.0;
	for ( uint64_t z = 1; z <= n ; ++z )
	{
		a *= lambda;
		b *= z;
	}
	
	double const c = a/b;
	double const d = 1.0 / ::std::exp(lambda);
	
	return c*d;
}
