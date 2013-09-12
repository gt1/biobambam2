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
#if ! defined(BIOBAMBAM_KMERPROB_HPP)
#define BIOBAMBAM_KMERPROB_HPP

#include <libmaus/types/types.hpp>

/**
 * this function computes an approximation of the probability
 * of a query sequence of KA A bases, KC C bases, KG G bases and
 * KT T bases appearing exactly n times in a reference sequence 
 * of length L. This is based on a simplified model using a Poisson
 * distribution (see e.g. Zhou, Xie: Exact Distribution of the Occurence
 * Number for K-tuples Over an Alphabet of Non-Equal Probability Letters)
 * The distribution is based on the single base relative frequences
 * of pA, pC, pG and pT of the respective bases in the reference sequence
 **/
double kmerPoisson(
	uint64_t const L, // length of reference sequence
	uint64_t const KA, // number of A bases in query sequence
	uint64_t const KC, // number of C bases in query sequence
	uint64_t const KG, // number of G bases in query sequence
	uint64_t const KT, // number of T bases in query sequence
	uint64_t const n, // number of occurences of query sequence
	double const pA = 0.25,
	double const pC = 0.25,
	double const pG = 0.25,
	double const pT = 0.25 
);
#endif
