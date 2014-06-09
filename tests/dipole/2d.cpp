/*
Tests 2d dipole integration using cubature.

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
#include "cstdlib"
#include "Eigen/Core"
#include "iostream"
#include "vector"

#include "dipole.hpp"

using namespace std;
using namespace background_B;

int main()
{
	// radius of earth
	constexpr double Re = 6.3712e6;

	const std::vector<Eigen::Vector3d> centers{
		{5 * Re, 0 * Re, 0 * Re},
		{-5 * Re, 5 * Re, 0},
		{0, -5 * Re, 5 * Re}
	};
	const std::vector<std::array<size_t, 3>> dimensions{
		{0, 1, 9},
		{1, 2, 9},
		{0, 2, 9}
	};
	const std::vector<Eigen::Vector3d> lengths{
		{8 * Re, 8 * Re, 0},
		{6 * Re, 6 * Re, 0},
		{4 * Re, 4 * Re, 0},
		{2 * Re, 2 * Re, 0},
		{1 * Re, 1 * Re, 0},
		{0.5 * Re, 0.5 * Re, 0},
		{0.1 * Re, 0.1 * Re, 0}
	};

	std::vector<
		std::pair<
			Eigen::Vector3d,
			Eigen::Vector3d
		>
	> dipole_moments_positions{{
		get_earth_dipole_moment<Eigen::Vector3d>(),
		{0, 0, 0}
	}};

	for (const auto& center: centers) {
	for (const auto& dimension: dimensions) {
	for (const auto& length: lengths) {

		Eigen::Vector3d start(center);
		start[dimension[0]] -= length[0] / 2.0;
		start[dimension[1]] -= length[1] / 2.0;

		// reference solution
		Eigen::Vector3d reference(0, 0, 0);
		constexpr size_t ref_samples = 1000;
		for (size_t i = 0; i <= ref_samples; i++) {
		for (size_t j = 0; j <= ref_samples; j++) {
			auto ref_pos(start);
			ref_pos[dimension[0]] += double(i) / ref_samples * length[0];
			ref_pos[dimension[1]] += double(j) / ref_samples * length[1];
			reference += get_dipole_field_0d(
				dipole_moments_positions,
				ref_pos
			);
		}}
		reference /= (ref_samples + 1) * (ref_samples + 1);


		Eigen::Vector3d integrated, error;

		std::tie(
			integrated,
			error
		) = get_dipole_field(
			dipole_moments_positions,
			start,
			dimension,
			length
		);

		const double relative_error
			= (reference - integrated).norm()
			/ std::max(reference.norm(), integrated.norm());

		if (relative_error > 1e-2) {
			std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
				<< "Relative error too large, integral:\n" << integrated
				<< "\nreference:\n" << reference
				<< "\nat:\n" << center / Re
				<< "\nfor lengths: " << length[0] / Re << " " << length[1] / Re
				<< " in dimensions: " << dimension[0] << " " << dimension[1]
				<< std::endl;
			abort();
		}
	}}}

	return EXIT_SUCCESS;
}
