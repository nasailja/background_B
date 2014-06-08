/*
Tests numerical integration using cubature.

Copyright 2014 Ilja Honkonen
All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include "array"
#include "cmath"
#include "cstdlib"
#include "iostream"
#include "limits"

#include "cubature.h"

int main()
{
	// same as odeint.cpp but with cubature

	const auto
		starts{0.0, -3.0, -10.0},
		ends{1.0, 2.0, 30.0};

	for (const auto start: starts) {
	for (const auto end: ends) {

		double
			integral = std::numeric_limits<double>::quiet_NaN(),
			error = std::numeric_limits<double>::quiet_NaN(),
			// test giving external data to pcubature
			value = 3;

		if (
			pcubature(
				1,
				[](
					unsigned,
					const double*,
					void* external_data,
					unsigned,
					double* fval
				){
					*fval = *(static_cast<double*>(external_data));
					return 0;
				},
				static_cast<void*>(&value),
				1,
				&start,
				&end, 
				0,
				0,
				1e-3, 
				ERROR_LINF,
				&integral,
				&error
			) != 0
		) {
			abort();
		}

		const double
			i = integral,
			// exact solution
			e = 3 * (end - start);

		if (std::fabs(i - e) > 1e-9 * std::max(std::fabs(i), std::fabs(e))) {
			std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
				<< "Too large relative error when integrating from " << start
				<< " to " << end << ", exact: " << e << ", numerical: " << i
				<< std::endl;
			abort();
		}
	}}


	for (const auto start: starts) {
	for (const auto end: ends) {

		double
			integral = std::numeric_limits<double>::quiet_NaN(),
			error = std::numeric_limits<double>::quiet_NaN();

		if (
			pcubature(
				1,
				[](
					unsigned,
					const double* x,
					void*,
					unsigned,
					double* fval
				){
					*fval = std::pow(*x, 3);
					return 0;
				},
				NULL,
				1,
				&start,
				&end, 
				0,
				0,
				1e-3, 
				ERROR_LINF,
				&integral,
				&error
			) != 0
		) {
			abort();
		}

		const double
			i = integral,
			e = (std::pow(end, 4) - std::pow(start, 4)) / 4;

		if (std::fabs(i - e) > 1e-9 * std::max(std::fabs(i), std::fabs(e))) {
			std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
				<< "Too large relative error when integrating from " << start
				<< " to " << end << ", exact: " << e << ", numerical: " << i
				<< std::endl;
			abort();
		}
	}}


	for (const auto start: starts) {
	for (const auto end: ends) {

		double
			integral = std::numeric_limits<double>::quiet_NaN(),
			error = std::numeric_limits<double>::quiet_NaN();

		if (
			pcubature(
				1,
				[](
					unsigned,
					const double* x,
					void*,
					unsigned,
					double* fval
				){
					*fval = std::sin(*x);
					return 0;
				},
				NULL,
				1,
				&start,
				&end, 
				0,
				0,
				1e-3, 
				ERROR_LINF,
				&integral,
				&error
			) != 0
		) {
			abort();
		}

		const double
			i = integral,
			e = -cos(end) + cos(start);

		if (std::fabs(i - e) > 1e-6 * std::max(std::fabs(i), std::fabs(e))) {
			std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
				<< "Too large relative error when integrating from " << start
				<< " to " << end << ", exact: " << e << ", numerical: " << i
				<< std::endl;
			abort();
		}
	}}


	/* from https://en.wikipedia.org/wiki/Integral */
	double
		integral = std::numeric_limits<double>::quiet_NaN(),
		error = std::numeric_limits<double>::quiet_NaN(),
		start = -2,
		end = 2;

	if (
		pcubature(
			1,
			[](
				unsigned,
				const double* x,
				void*,
				unsigned,
				double* fval
			){
				*fval = (
					(322 + 3 * (*x) * (98 + (*x) * (37 + (*x)))) / 100
					- 24 * (*x) / (1 + (*x)*(*x))
				) / 5;
				return 0;
			},
			NULL,
			1,
			&start,
			&end, 
			0,
			0,
			1e-3, 
			ERROR_LINF,
			&integral,
			&error
		) != 0
	) {
		abort();
	}

	double
		i = integral,
		e = 94.0 / 25.0;
	if (std::fabs(i - e) > 1e-9 * std::max(std::fabs(i), std::fabs(e))) {
		std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
			<< "Too large relative error when integrating from -2 to 2, exact: "
			<< e << ", numerical: " << i
			<< std::endl;
		abort();
	}


	// integrate f = x*y
	integral = std::numeric_limits<double>::quiet_NaN();
	error = std::numeric_limits<double>::quiet_NaN();
	if (
		pcubature(
			1,
			[](
				unsigned,
				const double* x,
				void*,
				unsigned,
				double* fval
			){
				*fval = (*x) * (*(x + 1));
				return 0;
			},
			NULL,
			2,
			std::array<double, 2>{-1, -1}.data(),
			std::array<double, 2>{2, 2}.data(),
			0,
			0,
			1e-3, 
			ERROR_LINF,
			&integral,
			&error
		) != 0
	) {
		abort();
	}

	i = integral,
	e = 2.25;
	if (std::fabs(i - e) > 1e-9 * std::max(std::fabs(i), std::fabs(e))) {
		std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
			<< "Too large relative error when integrating from -1,-1 to 2,2; exact: " << e << ", numerical: " << i
			<< std::endl;
		abort();
	}

	return EXIT_SUCCESS;
}
