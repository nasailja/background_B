/*
Dipolar field function(s).

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

#ifndef BACKGROUND_B_DIPOLE_HPP
#define BACKGROUND_B_DIPOLE_HPP

#include "limits"

#include "cubature.h"

/*!
\mainpage Function(s) for calculating averages of magnetic dipole(s).

\see
background_B::get_dipole_field_0d()
background_B::get_dipole_field()
*/


/*!
Namespace where all functions are defined.
*/
namespace background_B {


/*!
Returns magnetic field from a point dipole.

\param [dipole_moments_positions] List of dipole moments and
their corresponding positions
\param [field_position] Position of the returned total magnetic field

Assumes Vector_T is a 3d Eigen vector or similar.

Example usage:
\code
const Eigen::Vector3d field
	= background_B::get_dipole_field_0d(
		background_B::get_earth_dipole_moment<Eigen::Vector3d>(),
		{0, 0, 0},
		{6.3712e7, 0, 0} // at 10 Earth radii in +x
	);
\endcode

\see background_B::get_dipole_field()
*/
template <class Vector_T> Vector_T get_dipole_field_0d(
	const std::vector<std::pair<Vector_T, Vector_T>>& dipole_moments_positions,
	const Vector_T& field_position
) {
	static_assert(
		Vector_T::SizeAtCompileTime == 3,
		"ERROR: Only 3 component vectors supported"
	);

	Vector_T ret_val;
	ret_val.setZero();

	constexpr double vacuum_permeability = 1.256637061e-6;

	for (const auto dip_mom_pos: dipole_moments_positions) {

		const Vector_T r = field_position - dip_mom_pos.second;

		bool far_enough = false;
		for (size_t i = 0; i < r.size(); i++) {
			if (fabs(r[i]) > std::numeric_limits<typename Vector_T::Scalar>::epsilon()) {
				far_enough = true;
				break;
			}
		}
		if (not far_enough) {
			ret_val += 2.0 / 3.0 * vacuum_permeability * dip_mom_pos.first;
			continue;
		}

		const double r3 = r.norm() * r.squaredNorm();
		const Vector_T
			r_unit = r / r.norm(),
			projected_dip_mom = dip_mom_pos.first.dot(r_unit) * r_unit;

		ret_val += 1e-7 * (3 * projected_dip_mom - dip_mom_pos.first) / r3;

	}

	return ret_val;
}


/*!
Returns average magnetic field over an area from point dipole(s).

\param [dipole_moments_positions] is a list of dipole moments and the
corresponding positions
\param [start] is the starting point of averaging
\param [dimensions] is a list of dimensions over which to integrate as
values in the range [0, 2], with rest of the dimensions as
values > 2. Dimensions over which to integrate must be listed
first and must be listed in ascending order.
\param [lengths] is a list of lengths for which to integrate from start
in dimensions at the corresponding index in dimensions
\param [maxEval] See http://ab-initio.mit.edu/wiki/index.php/Cubature
\param [reqAbsError] See http://ab-initio.mit.edu/wiki/index.php/Cubature
\param [reqRelError] See http://ab-initio.mit.edu/wiki/index.php/Cubature

Assumes Vector_T is a 3d Eigen vector or similar.

In case of error returns negative values for the accuracy of integration.

Example calculating a 1d integral in x dimension:
\code
Eigen::Vector3d integrated, error;

std::tie(
	integrated,
	error
) = get_dipole_field<Eigen::Vector3d>(
	// earth dipole centered at origin
	{{{0, 0, -7.94e22}, {0, 0, 0}}},
	// integrate from (1, 0, 0) earth radii
	{6.3712e6, 0, 0},
	// integrate in x dimension
	{0, 9, 9},
	// integrate for 1 earth radii in +x
	{6.3712e6, 0, 0}
);
\endcode

Example calculating a 2d integral in y and z dimensions:
\code
Eigen::Vector3d integrated, error;

std::tie(
	integrated,
	error
) = get_dipole_field<Eigen::Vector3d>(
	{{{0, 0, -7.94e22}, {0, 0, 0}}},
	{6.3712e6, 0, 0},
	// integrate in y and z dimensions
	{1, 2, 9},
	// integrate for 1 earth radii in +y and +z
	{0, 6.3712e6, 6.3712e6}
);
\endcode
*/
template <class Vector_T> std::pair<Vector_T, Vector_T> get_dipole_field(
	const std::vector<std::pair<Vector_T, Vector_T>>& dipole_moments_positions,
	const Vector_T& start,
	const std::array<size_t, 3>& dimensions,
	const Vector_T& lengths,
	const size_t maxEval = 0,
	const double reqAbsError = 0,
	const double reqRelError = 1e-6
) {
	static_assert(
		Vector_T::SizeAtCompileTime == 3,
		"ERROR: Only 3 component vectors supported"
	);

	if (dimensions[0] > 3) {
		return std::make_pair(
			get_dipole_field_0d(
				dipole_moments_positions,
				start
			),
			Vector_T(0, 0, 0)
		);
	}

	// check given dimensions
	int largest = -1;
	for (size_t i = 0; i < dimensions.size(); i++) {
		if (dimensions[i] > 3) {
			break;
		}

		if (largest < int(dimensions[i])) {
			largest = int(dimensions[i]);
		} else  {
			return std::make_pair(
				Vector_T(0, 0, 0),
				Vector_T(-1, -1, -1)
			);
		}
	}

	std::tuple<
		std::array<size_t, 3>,
		std::array<double, 3>,
		std::vector<std::pair<Vector_T, Vector_T>>
	> integration_parameters;

	std::get<0>(integration_parameters) = dimensions;
	std::get<1>(integration_parameters)[0] = start[0];
	std::get<1>(integration_parameters)[1] = start[1];
	std::get<1>(integration_parameters)[2] = start[2];
	std::get<2>(integration_parameters) = dipole_moments_positions;

	/*
	Get r_dims, normalization and integration
	range in format used by integrand
	*/
	unsigned r_dims = 0;
	double normalization = 1;
	Vector_T cubature_start(start), cubature_end(start);
	for (size_t i = 0; i < dimensions.size(); i++) {
		if (dimensions[i] > 2) {
			break;
		}

		r_dims++;

		normalization *= lengths[dimensions[i]];

		cubature_start[i] = start[dimensions[i]];
		cubature_end[i] = cubature_start[i] + lengths[dimensions[i]];
	}

	Vector_T integral, error;


	/*
	Returns the dipole field(s) at given point for cubature.

	In dimension(s) over which pcubature is integrating, x, y or z
	are given by r, the rest are given in extra_data. The location
	and dipole moments of dipoles are also given in extra_data.
	Format of extra_data:
	std::tuple<
		std::array<size_t, 3>,
		std::array<double, 3>,
		std::vector<std::pair<Vector_T, Vector_T>>
	>
	where:
	size_t array has the dimensions over which cubature will
	integrate. Only the first r_dims values are used, e.g. if
	integrating over y and z, r_dims == 2 and the first two values
	of the array must be 1 and 2.
	double array tells the coordinates over which cubature is not
	integrating, only the values not set by cubature are used e.g.
	if integrating over x and z at y == -1 then the second value
	must be -1 and the others can be left unspecified.
	The vector has all dipole moments and their positions that
	should be included in the integral.

	For example when integrating x and z from 0 to 1 with y == -2
	integrand could be called by cubature as (999's mark unused data):
	f = integrand(
		2,
		std::array<double, 2>{0, 0}.data(),
		std::tuple<
			std::array<size_t, 3>,
			std::array<double, 3>,
			std::vector<std::pair<Vector_T, Vector_T>>
		>{{0, 2, 999}, {999, -2, 999}, {{0, 0, -7.94e22}, {0, 0, 0}}}.data(),
		3,
		ret_val.data()
	);
	and when integrating y from -2 to -1 with x == 0.5 and z == 2/3:
	f = integrand(
		1,
		std::array<double, 1>{-2}.data(),
		std::pair<
			std::array<size_t, 3>,
			std::array<double, 3>
		>{{1, 999, 999}, {0.5, 999, 2/3}, {{0, 0, -7.94e22}, {0, 0, 0}}.data(),
		3,
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

			if (r_dims > 3) {
				return -2;
			}

			if (f_dims != 3) {
				return -3;
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

			const auto integration_parameters
				= *static_cast<
					std::tuple<
						std::array<size_t, 3>,
						std::array<double, 3>,
						std::vector<std::pair<Vector_T, Vector_T>>
					>*
				>(extra_data);

			Vector_T real_r(std::get<1>(integration_parameters).data());

			for (size_t i = 0; i < r_dims; i++) {
				if (std::get<0>(integration_parameters)[i] > 3) {
					std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
						<< " index " << i
						<< ", value " << std::get<0>(integration_parameters)[i]
						<< std::endl;
					abort();
				}

				real_r[std::get<0>(integration_parameters)[i]] = *(r + i);
			}

			Eigen::Map<Vector_T> ret_val(f);

			ret_val = get_dipole_field_0d(
				std::get<2>(integration_parameters),
				real_r
			);

			return 0;
		};


	const int result
		= pcubature(
			3,
			integrand,
			static_cast<void*>(&integration_parameters),
			r_dims,
			cubature_start.data(),
			cubature_end.data(),
			maxEval,
			reqAbsError,
			reqRelError,
			ERROR_LINF,
			integral.data(),
			error.data()
		);

	if (result != 0) {
		std::cerr <<  __FILE__ << "(" << __LINE__ << "): "
			<< "Could not calculate line integral of a dipole, return value: "
			<< result
			<< std::endl;
		abort();
	}

	return std::make_pair(integral / normalization, error);
}


/*!
Returns the Earth's dipole moment pointing in -Z direction.

Magnitude of the moment is 7.94e22 (A m^2), from
seismo.berkeley.edu/~rallen/eps122/lectures/L05.pdf
*/
template <class Vector_T> Vector_T get_earth_dipole_moment()
{
	return {0, 0, -7.94e22};	
}

/*!
Returns the Jupiter's dipole moment pointing in -Z direction.

Magnitude of the moment is 1.43e27 (A m^2), from
http://dx.doi.org/10.1016/j.epsl.2006.08.008
*/
template <class Vector_T> Vector_T get_jupiter_dipole_moment()
{
	return {0, 0, -1.43e27};	
}

/*!
Returns the Saturn's dipole moment pointing in -Z direction.

Magnitude of the moment is 4.61e25 (A m^2), from
http://dx.doi.org/10.1016/j.epsl.2006.08.008
*/
template <class Vector_T> Vector_T get_saturn_dipole_moment()
{
	return {0, 0, -4.61e25};	
}

/*!
Returns the Uranus' dipole moment pointing in -Z direction.

Magnitude of the moment is 3.97e24 (A m^2), from
http://dx.doi.org/10.1016/j.epsl.2006.08.008
*/
template <class Vector_T> Vector_T get_uranus_dipole_moment()
{
	return {0, 0, -3.97e24};	
}

/*!
Returns the Neptune's dipole moment pointing in -Z direction.

Magnitude of the moment is 1.91e24 (A m^2), from
http://dx.doi.org/10.1016/j.epsl.2006.08.008
*/
template <class Vector_T> Vector_T get_neptune_dipole_moment()
{
	return {0, 0, -1.91e24};	
}

} // namespace

#endif
