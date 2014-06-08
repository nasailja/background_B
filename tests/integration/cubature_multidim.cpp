/*
Tests multidimensional numerical integration using cubature.

Copyright 2014 Ilja Honkonen

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

bool good_solution(
	const double integral,
	const double exact,
	const double max_relative_error
) {
	if (
		std::fabs(integral - exact)
		> max_relative_error * std::max(std::fabs(integral), std::fabs(exact))
	) {
		return false;
	}

	return true;
}


int main()
{
	constexpr double max_acceptable_error = 1e-9;

	// integration range
	const std::array<double, 3> start{1, 1, 1}, end{2, 2, 2};

	double integral, error, exact;

	/*
	Returns x*y^2*z^3.

	In dimension(s) over which pcubature is integrating, x, y or z
	are given by r, the rest are given in extra_data. Format of
	extra_data: std::pair<std::array<size_t, 3>, std::array<double, 3>>.
	size_t tells the index of the dimension over which integration over
	Nth dimension in r actually goes. The first r_dims of
	size_t's contain valid data. For the rest the corresponding
	double contains the coordinate at which to integrate over the
	rest of dimensions of r.

	For example when integrating x and z from 0 to 1 with y == -2
	integrand could be called by cubature as (999's mark unused data):
	f = integrand(
		2,
		std::array<double, 2>{0, 0}.data(),
		std::pair<
			std::array<size_t, 3>,
			std::array<double, 3>
		>{{0, 2, 999}, {999, -2, 999}}.data(),
		1,
		ret_val.data()
	);
	and when integrating y from -2 to -1 with x == 0.5 and z == 2/3:
	f = integrand(
		1,
		std::array<double, 1>{-2}.data(),
		std::pair<
			std::array<size_t, 3>,
			std::array<double, 3>
		>{{1, 999, 999}, {0.5, 999, 2/3}}.data(),
		1,
		ret_val.data()
	);
	*/
	auto integrand
		= [](
			unsigned r_dims,
			const double* r,
			void* extra_data,
			unsigned f_dims,
			double* f
		) {
			if (r_dims == 0) {
				return -1;
			}

			if (f_dims != 1) {
				return -2;
			}

			if (r == NULL) {
				std::cerr <<  __FILE__ << "(" << __LINE__ << ")" << std::endl;
				abort();
			}

			if (f == NULL) {
				std::cerr <<  __FILE__ << "(" << __LINE__ << ")" << std::endl;
				abort();
			}

			if (extra_data == NULL) {
				std::cerr <<  __FILE__ << "(" << __LINE__ << ")" << std::endl;
				abort();
			}

			const auto integration_info
				= *static_cast<
					std::pair<
						std::array<size_t, 3>,
						std::array<double, 3>
					>*
				>(extra_data);

			std::array<double, 3> real_r = integration_info.second;

			for (size_t i = 0; i < r_dims; i++) {
				if (integration_info.first[i] > 3) {
					std::cerr <<  __FILE__ << "(" << __LINE__ << ")" << std::endl;
					abort();
				}

				real_r[integration_info.first[i]] = *(r + i);
			}

			*f
				= std::pow(real_r[0], 1)
				* std::pow(real_r[1], 2)
				* std::pow(real_r[2], 3);

			return 0;
		};

	std::pair<std::array<size_t, 3>, std::array<double, 3>> integration_parameters;

	// 1d over x at y == 3/2 and z == 3/2
	exact = std::pow(3./2., 6);
	integration_parameters.first[0] = 0;
	integration_parameters.first[1] = 999;
	integration_parameters.first[2] = 999;
	integration_parameters.second[0] = 999;
	integration_parameters.second[1] = 3.0 / 2.0;
	integration_parameters.second[2] = 3.0 / 2.0;
	if (
		pcubature(
			1,
			integrand,
			&integration_parameters,
			1,
			start.data(),
			end.data(),
			0,
			0,
			1e-3,
			ERROR_LINF,
			&integral,
			&error
		) != 0
	) {
		std::cerr <<  __FILE__ << "(" << __LINE__ << ")" << std::endl;
		abort();
	}

	if (not good_solution(integral, exact, max_acceptable_error)) {
		std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
			<< "Too large relative error when integrating from " << start[0]
			<< " to " << end[0]
			<< ", exact: " << exact
			<< ", numerical: " << integral
			<< std::endl;
		abort();
	}

	// 2d over y and z at x == 3/2
	exact = 105.0 / 8.0;
	integration_parameters.first[0] = 1;
	integration_parameters.first[1] = 2;
	integration_parameters.first[2] = 999;
	integration_parameters.second[0] = 3.0 / 2.0;
	integration_parameters.second[1] = 999;
	integration_parameters.second[2] = 999;
	if (
		pcubature(
			1,
			integrand,
			&integration_parameters,
			2,
			start.data(),
			end.data(),
			0,
			0,
			1e-3,
			ERROR_LINF,
			&integral,
			&error
		) != 0
	) {
		std::cerr <<  __FILE__ << "(" << __LINE__ << ")" << std::endl;
		abort();
	}

	if (not good_solution(integral, exact, max_acceptable_error)) {
		std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
			<< "Too large relative error when integrating from "
			<< start[1] << "," << start[2]
			<< " to " << end[1] << "," << end[2]
			<< ", exact: " << exact
			<< ", numerical: " << integral
			<< std::endl;
		abort();
	}

	// integrate over the edges of the unit cube from start to end
	std::array<
		std::pair<
			// exact value for below integration parameters
			double,
			// same as in integration_parameters
			std::pair<
				std::array<size_t, 3>,
				std::array<double, 3>
			>
		>,
		12
	> integration_params1d{{
		// over x at min y, min z
		{3.0 / 2.0, {{0, 999, 999}, {999, 1, 1}}},
		// x at min y, max z
		{12, {{0, 999, 999}, {999, 1, 2}}},
		// x at max y, min z
		{6, {{0, 999, 999}, {999, 2, 1}}},
		// x at max y, max z
		{48, {{0, 999, 999}, {999, 2, 2}}},
		// over y at min x, min z
		{7.0 / 3.0, {{1, 999, 999}, {1, 999, 1}}},
		// y at min x, max z
		{56.0 / 3.0, {{1, 999, 999}, {1, 999, 2}}},
		// y at max y, min z
		{14.0 / 3.0, {{1, 999, 999}, {2, 999, 1}}},
		// y at max x, max z
		{112.0 / 3.0, {{1, 999, 999}, {2, 999, 2}}},
		// over z at min x, min y
		{15.0 / 4.0, {{2, 999, 999}, {1, 1, 999}}},
		// z at min x, max y
		{15, {{2, 999, 999}, {1, 2, 999}}},
		// z at max x, min y
		{15.0 / 2.0, {{2, 999, 999}, {2, 1, 999}}},
		// z at max x, max y
		{30, {{2, 999, 999}, {2, 2, 999}}}
	}};

	for (auto& info: integration_params1d) {

		const double exact = info.first;

		if (
			pcubature(
				1,
				integrand,
				&info.second,
				1,
				start.data(),
				end.data(),
				0,
				0,
				1e-3,
				ERROR_LINF,
				&integral,
				&error
			) != 0
		) {
			std::cerr <<  __FILE__ << "(" << __LINE__ << ")" << std::endl;
			abort();
		}

		if (not good_solution(integral, exact, max_acceptable_error)) {
			std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
				<< "Too large relative error, exact: " << exact
				<< ", numerical: " << integral
				<< std::endl;
			abort();
		}
	}

	// integrate over the faces of the unit cube from start to end
	std::array<
		std::pair<
			double,
			std::pair<
				std::array<size_t, 3>,
				std::array<double, 3>
			>
		>,
		6
	> integration_params2d{{
		// over x and y at min z
		{7.0 / 2.0, {{0, 1, 999}, {999, 999, 1}}},
		// x, y at max z
		{28, {{0, 1, 999}, {999, 999, 2}}},
		// x, z at min y
		{45.0 / 8.0, {{0, 2, 999}, {999, 1, 999}}},
		// x, z at max y
		{45.0 / 2.0, {{0, 2, 999}, {999, 2, 999}}},
		// y, z at min x
		{35.0 / 4.0, {{1, 2, 999}, {1, 999, 999}}},
		// y, z at max x
		{35.0 / 2.0, {{1, 2, 999}, {2, 999, 999}}}
	}};

	for (auto& info: integration_params2d) {

		const double exact = info.first;

		if (
			pcubature(
				1,
				integrand,
				&info.second,
				2,
				start.data(),
				end.data(),
				0,
				0,
				1e-3,
				ERROR_LINF,
				&integral,
				&error
			) != 0
		) {
			std::cerr <<  __FILE__ << "(" << __LINE__ << ")" << std::endl;
			abort();
		}

		if (not good_solution(integral, exact, max_acceptable_error)) {
			std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
				<< "Too large relative error, exact: " << exact
				<< ", numerical: " << integral
				<< std::endl;
			abort();
		}
	}

	return EXIT_SUCCESS;
}
